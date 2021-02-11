# coding: utf-8
#
#    Project: BioCAT SAXS pipeline
#             https://github.com/biocatiit/saxs-pipeline
#
#
#    Principal author:       Jesse Hopkins
#
#    This is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This software is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this software.  If not, see <http://www.gnu.org/licenses/>.

import threading
import time
import multiprocessing as mp
import queue
import os
import collections
import copy
import logging
import traceback

if __name__ != '__main__':
    logger = logging.getLogger(__name__)

import bioxtasraw.RAWAPI as raw

pipeline_path = os.path.abspath(os.path.join('.', __file__, '..', '..'))
if pipeline_path not in os.sys.path:
    os.sys.path.append(pipeline_path)

from pipeline.analysis import analyze
from pipeline.reduction import reduce_data


class pipeline_thread(threading.Thread):

    def __init__(self, cmd_q, ret_q, abort_event, pipeline_settings, ret_lock=threading.RLock()):
        threading.Thread.__init__(self)
        self.daemon = True
        self.name = 'pipeline_ctrl_thread'

        logger.info("Initializing pipeline control thread")

        self._cmd_q = cmd_q
        self._ret_q = ret_q
        self._ret_lock = ret_lock
        self._abort_event = abort_event
        self._stop_event = threading.Event()

        self._ret_lock.acquire()
        self.data_dir = None
        self.output_dir = None
        self.profiles_dir = None
        self.analysis_dir = None

        self.fprefix = ''

        self.exp_total = 0
        self.exp_processed  = 0
        self.exp_being_processed = 0

        self.experiments = {}
        self.current_experiment = ''

        self.active = False

        self._commands = {'set_data_dir'    : self._set_data_dir,
            'set_fprefix'       : self._set_fprefix,
            'set_data_dir_and_fprefix'  : self._set_data_dir_and_fprefix,
            'set_output_dir'    : self._set_output_dir,
            'set_profiles_dir'  : self._set_profiles_dir,
            'set_analysis_dir'  : self._set_analysis_dir,
            'add_reduction_process' : self._start_reduction_process,
            'add_analysis_process'  : self._start_analysis_process,
            'load_raw_settings' : self._load_raw_settings,
            'start_experiment'  : self._start_experiment,
            'stop_experiment'   : self._stop_experiment,
            'update_pipeline_settings': self._update_pipeline_settings,
            'update_analysis_args'  : self._set_analysis_args,
            }

        #Initialize reduction and analysis threads
        self.reduction_processes = []
        self.analysis_processes = []

        self.num_reduction_processes = 0
        self.num_analysis_processes = 0

        self.manager = mp.Manager()

        #Reduction queues, locks, etc
        self.r_cmd_q = self.manager.Queue()
        self.r_ret_q = self.manager.Queue()
        self.r_cmd_lock = self.manager.Lock()
        self.r_ret_lock = self.manager.Lock()
        self.r_abort_event = self.manager.Event()

        #Processing queues, locks, etc
        self.a_cmd_q = self.manager.Queue()
        self.a_ret_q = self.manager.Queue()
        self.a_cmd_lock = self.manager.Lock()
        self.a_ret_lock = self.manager.Lock()
        self.a_abort_event = self.manager.Event()

        self.mp_log_lock = self.manager.Lock()
        self.mp_log_queue = self.manager.Queue()

        self.mp_log_thread = mp_log_thread(self.mp_log_lock, self.mp_log_queue)
        self.mp_log_thread.start()

        self.num_loaded = self.manager.Value('i', 0)
        self.num_averaged = self.manager.Value('i', 0)
        self.pl_settings = self.manager.dict()
        self._update_pipeline_settings(**{kw : pipeline_settings[kw] for kw in pipeline_settings})

        self.raw_settings_file = self.pl_settings['raw_settings_file']
        self.raw_settings = raw.load_settings(self.raw_settings_file)

        # For the moment, current architechture doesn't allow more than 1 reduction process
        # Can revist if necessary, will leave everything in place.
        # for i in range(self.pl_settings['r_procs']):
        #     self._start_reduction_process()

        for i in range(1):
            self._start_reduction_process()

        for i in range(self.pl_settings['a_procs']):
            self._start_analysis_process()

        #Initialize save thread
        self.s_cmd_q = collections.deque()
        self.s_ret_q = collections.deque()
        self.s_abort_event = threading.Event()

        self.save_thread = reduce_data.save_thread(self.s_cmd_q, self.s_ret_q,
            self.s_abort_event, self.raw_settings)

        self.save_thread.start()

        self._ret_lock.release()



    def run(self):
        while True:
            try:
                self.active = False

                if self._stop_event.is_set():
                    break

                if self._abort_event.is_set():
                    self._abort()

                try:
                    cmd, args, kwargs = self._cmd_q.popleft()
                except IndexError:
                    cmd = None

                if cmd is not None:
                    logger.info("Processing cmd '%s' with args: %s and kwargs: %s ",
                        cmd, ', '.join(['{}'.format(a) for a in args]),
                        ', '.join(['{}: {}'.format(kw, item) for kw, item in kwargs.items()]))

                    self._commands[cmd](*args, **kwargs)

                try:
                    with self.r_ret_lock:
                        profile_data = self.r_ret_q.get_nowait()
                except queue.Empty:
                    profile_data = None

                if profile_data is not None:
                    self.active = True

                    self._save_profiles(profile_data)
                    self._add_profiles_to_experiment(profile_data)

                self._check_exp_status()

                try:
                    with self.a_ret_lock:
                        results = self.a_ret_q.get_nowait()
                except queue.Empty:
                    results = None

                if results is not None:
                    self.active = True

                    self._add_analysis_to_experiment(results)

                self._check_analysis_status()

                if not self.active:
                    time.sleep(0.1)

            except Exception:
                logger.error('Error in pipeline thread:\n{}'.format(traceback.format_exc()))
                time.sleep(0.1)

        logger.info("Quitting pipeline control thread")

    def _start_reduction_process(self):
        self._ret_lock.acquire()
        logger.debug('Starting reduction process %i', len(self.reduction_processes)+1)

        managed_internal_q = self.manager.Queue()

        proc = reduce_data.raver_process(self.r_cmd_q, self.r_ret_q,
            self.r_cmd_lock, self.r_ret_lock, self.r_abort_event,
            self.raw_settings_file, self.pl_settings, self.mp_log_lock,
            self.mp_log_queue, managed_internal_q, self.num_loaded,
            self.num_averaged)

        proc.start()

        self.reduction_processes.append(proc)

        self.num_reduction_processes += 1

        self._ret_lock.release()

    def _start_analysis_process(self):
        self._ret_lock.acquire()
        logger.debug('Starting analysis process %i', len(self.analysis_processes)+1)

        proc = analyze.analysis_process(self.a_cmd_q, self.a_ret_q,
            self.a_cmd_lock, self.a_ret_lock, self.a_abort_event,
            self.raw_settings_file, self.mp_log_lock, self.mp_log_queue)

        proc.start()

        self.analysis_processes.append(proc)

        self.num_analysis_processes += 1

        self._ret_lock.release()

    def _set_data_dir(self, data_dir):
        logger.debug('Setting data directory: %s', data_dir)

        self._ret_lock.acquire()
        self.data_dir = os.path.abspath(os.path.expanduser(data_dir))
        self._ret_lock.release()

        if self.pl_settings['use_default_output_dir']:
            output_dir = os.path.join(data_dir,
                self.pl_settings['default_output_dir'])

            self._set_output_dir(output_dir)

        self.r_cmd_q.put_nowait(['set_data_dir', [copy.deepcopy(self.data_dir)], {}])

    def _set_fprefix(self, fprefix):
        logger.debug('Setting file prefix: %s' %fprefix)

        self._ret_lock.acquire()
        self.fprefix = fprefix
        self._ret_lock.release()

        self.r_cmd_q.put_nowait(['set_fprefix', [self.fprefix], {}])

    def _set_data_dir_and_fprefix(self, data_dir, fprefix):
        logger.debug('Setting data directory and file prefix: %s, %s', data_dir,
            fprefix)

        self._ret_lock.acquire()
        self.data_dir = os.path.abspath(os.path.expanduser(data_dir))
        self.fprefix = fprefix
        self._ret_lock.release()

        if self.pl_settings['use_default_output_dir']:
            output_dir = os.path.join(data_dir,
                self.pl_settings['default_output_dir'])

            self._set_output_dir(output_dir)

        self.r_cmd_q.put_nowait(['set_data_dir_and_fprefix', [copy.deepcopy(self.data_dir),
            copy.deepcopy(self.fprefix)], {}])

    def _set_output_dir(self, output_dir):
        logger.debug('Setting output directory: %s', output_dir)

        self._ret_lock.acquire()
        self.output_dir = os.path.abspath(os.path.expanduser(output_dir))
        self._ret_lock.release()

        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)

        if self.current_experiment in self.experiments:
            self.experiments[self.current_experiment].output_dir = output_dir

        if self.pl_settings['use_default_profiles_dir']:
            profiles_dir = os.path.join(output_dir,
                self.pl_settings['default_profiles_dir'])

            self._set_profiles_dir(profiles_dir)

        if self.pl_settings['use_default_analysis_dir']:
            analysis_dir = os.path.join(output_dir,
                self.pl_settings['default_analysis_dir'])

            self._set_analysis_dir(analysis_dir)

    def _set_profiles_dir(self, profiles_dir):
        logger.debug('Setting profiles output directory: %s', profiles_dir)

        self._ret_lock.acquire()
        self.profiles_dir = os.path.abspath(os.path.expanduser(profiles_dir))
        self._ret_lock.release()

        if self.current_experiment in self.experiments:
            self.experiments[self.current_experiment].profiles_dir = profiles_dir

        if not os.path.exists(self.profiles_dir):
            os.mkdir(self.profiles_dir)

    def _set_analysis_dir(self, analysis_dir):
        logger.debug('Setting analysis output directory: %s', analysis_dir)

        self._ret_lock.acquire()
        self.analysis_dir = os.path.abspath(os.path.expanduser(analysis_dir))
        self._ret_lock.release()

        if self.current_experiment in self.experiments:
            self.experiments[self.current_experiment].analysis_dir = analysis_dir

        if not os.path.exists(self.analysis_dir):
            os.mkdir(self.analysis_dir)

    def _load_raw_settings(self, raw_settings_file):
        logger.debug('Loading RAW settings from %s', raw_settings_file)

        self.raw_settings = raw.load_settings(raw_settings_file)

        self.s_cmd_q.append(['set_raw_settings', [self.raw_settings]])

        self._ret_lock.acquire()
        for proc in self.reduction_processes:
            proc.load_settings(raw_settings_file)

        for proc in self.analysis_processes:
            proc.load_settings(raw_settings_file)

        self._ret_lock.release()

    def _start_experiment(self, *args, **kwargs):
        logger.debug('Starting experiment with args: %s and kwargs: %s',
            ', '.join(['{}'.format(a) for a in args]),
            ', '.join(['{}: {}'.format(kw, item) for kw, item in kwargs.items()]))

        exp_name = args[0]
        exp_type = args[1]
        data_dir = args[2]
        fprefix = args[3]

        if exp_name not in self.experiments:
            new_exp = Experiment(exp_name, exp_type, data_dir, fprefix,
                **kwargs)

            self.experiments[exp_name] = new_exp

            self._ret_lock.acquire()
            self.exp_total += 1
            self._ret_lock.release()

        else:
            self.experiments[exp_name].current_fprefix = fprefix

        self.current_experiment = exp_name

        self._ret_lock.acquire()
        self.data_dir = os.path.abspath(os.path.expanduser(data_dir))
        self.fprefix = fprefix
        self._ret_lock.release()

        if self.pl_settings['use_default_output_dir']:
            output_dir = os.path.join(data_dir,
                self.pl_settings['default_output_dir'])

            self._set_output_dir(output_dir)

        self.r_cmd_q.put_nowait(['set_experiment', [copy.deepcopy(self.data_dir),
            copy.deepcopy(self.fprefix), self.current_experiment], {}])

    def _stop_experiment(self, exp_name):
        logger.debug('Stopping experiment %s', exp_name)

        if exp_name in self.experiments:
            self.experiments[exp_name].stop_experiment()

    def _save_profiles(self, profile_data):
        logger.debug('Saving profiles')
        logger.debug('Profile data: {}'.format(profile_data))

        exp_id, profiles = profile_data

        if self.pl_settings['save_raver_profiles']:
            if exp_id in self.experiments:
                profiles_dir = self.experiments[exp_id].profiles_dir
            else:
                profiles_dir = self.profiles_dir

            if profiles_dir is not None:
                self.s_cmd_q.append(['save_profiles', [copy.deepcopy(profiles_dir), profiles]])

    def _add_profiles_to_experiment(self, profile_data):
        logger.debug('Adding profiles to experiment')

        exp_id, profiles = profile_data

        if exp_id in self.experiments:
            self.experiments[exp_id].add_profiles(profiles)

    def _check_exp_status(self):
        logger.debug('Checking experiment status')

        save_proc_data = self.pl_settings['save_processed_data']
        save_report = self.pl_settings['save_report']
        report_type = self.pl_settings['report_type']

        for exp in self.experiments.values():
            if not exp.collection_finished:
                exp.check_exp_timeout(self.pl_settings)

            if exp.collection_finished and exp.analysis_last_modified == -1:
                self.active = True

                logger.info('Experiment %s data collection finished', exp.exp_name)

                self._ret_lock.acquire()
                self.exp_being_processed += 1
                self._ret_lock.release()

                with self.a_cmd_lock:
                    if exp.exp_type == 'SEC':
                        self.a_cmd_q.put_nowait(['make_and_analyze_series',
                            exp.exp_name, [exp.analysis_dir, exp.profiles,
                            save_proc_data, save_report, report_type,
                            exp.output_dir], self._analysis_args])
                        exp.analysis_last_modified = time.time()

                    elif exp.exp_type == 'Batch':
                        self.a_cmd_q.put_nowait(['subtract_and_analyze_batch',
                            exp.exp_name, [exp.analysis_dir, exp.profiles,
                            exp.buffer_profiles, save_proc_data, save_report,
                            report_type, exp.output_dir], self._analysis_args])
                        exp.analysis_last_modified = time.time()

    def _add_analysis_to_experiment(self, results):
        logger.debug('Adding analysis to experiment')

        res_type = results[0]
        exp_id = results[1]

        if exp_id in self.experiments:
            exp = self.experiments[exp_id]

            if res_type == 'sub_series':
                series = results[2]

                exp.series = series

            elif res_type == 'analysis_results':
                res_dict = results[2]

                if res_dict is not None:
                    exp.sub_profile = res_dict['profile']
                    exp.ift = res_dict['ift']
                    exp.dammif_data = res_dict['dammif_data']
                    exp.denss_data = res_dict['denss_data']

                exp.analysis_last_modified = time.time()
                exp.analysis_finished = True

    def _check_analysis_status(self):
        logger.debug('Checking analysis status')

        exp_keys = list(self.experiments.keys())

        for exp_key in exp_keys:
            exp = self.experiments[exp_key]

            if not exp.analysis_finished:
                exp.check_analysis_timeout(self.pl_settings)

            if exp.analysis_finished:
                logger.info('Experiment %s analysis finished', exp.exp_name)

                self._ret_lock.acquire()
                self.exp_being_processed += -1
                self.exp_processed += 1
                self.exp_total += -1
                self._ret_lock.release()

                self._ret_q.append(exp)
                del self.experiments[exp.exp_name]

    def _set_analysis_args(self):
        self._analysis_args = {
            'dammif'        : self.pl_settings['dammif'],
            'denss'         : self.pl_settings['denss'],
            'dam_runs'      : self.pl_settings['dam_runs'],
            'dam_mode'      : self.pl_settings['dam_mode'],
            'damaver'       : self.pl_settings['damaver'],
            'damclust'      : self.pl_settings['damclust'],
            'dam_refine'    : self.pl_settings['dam_refine'],
            'denss_runs'    : self.pl_settings['denss_runs'],
            'denss_mode'    : self.pl_settings['denss_mode'],
            'denss_aver'    : self.pl_settings['denss_aver'],
            'denss_refine'  : self.pl_settings['denss_refine'],
            'use_atsas'     : self.pl_settings['use_atsas'],
            }

        logger.debug('Analysis arguments: {}'.format(self._analysis_args))

    def _update_pipeline_settings(self, *args, **kwargs):

        logger.debug('Updating pipeline settings: %s',
            ', '.join(['{}: {}'.format(kw, item) for kw, item in kwargs.items()]))

        self.pl_settings.update(kwargs)

        self._set_analysis_args()

    def get_num_loaded(self):
        return copy.deepcopy(self.num_loaded.get())

    def get_num_averaged(self):
        return copy.deepcopy(self.num_averaged.get())

    def _abort(self):
        logger.debug('Aborting pipeline')
        self._abort_event.set()

        self._cmd_q.clear()
        self._ret_q.clear()

        self.s_abort_event.set()

        self.r_abort_event.set()
        self.a_abort_event.set()

        r_aborts = 0

        self._ret_lock.acquire()

        while r_aborts < len(self.reduction_processes):
            try:
                with self.r_ret_lock:
                    ret = self.r_ret_q.get_nowait()

                if ret == 'aborted':
                    r_aborts += 1

            except queue.Empty:
                time.sleep(0.1)

        a_aborts = 0

        while a_aborts < len(self.analysis_processes):
            try:
                with self.a_ret_lock:
                    ret = self.a_ret_q.get_nowait()

                if ret == 'aborted':
                    a_aborts += 1

                elif ret is not None:
                    self._add_analysis_to_experiment(ret)

            except queue.Empty:
                time.sleep(0.1)

        self._ret_lock.release()

        self.r_abort_event.clear()
        self.a_abort_event.clear()

        self._abort_event.clear()

    def stop(self):
        """Stops the thread cleanly."""
        self._abort_event.set()

        while self._abort_event.is_set():
            time.sleep(0.1)

        self.save_thread.stop()
        self.save_thread.join()

        self._ret_lock.acquire()

        for proc in self.analysis_processes:
            proc.stop()
            proc.join()

        for proc in self.reduction_processes:
            proc.stop()
            proc.join()

        self._ret_lock.release()

        self.mp_log_thread.stop()
        self.mp_log_thread.join()

        self._stop_event.set()

class mp_log_thread(threading.Thread):

    def __init__(self, log_lock, log_q):
        threading.Thread.__init__(self)
        self.daemon = True
        self.name = ''

        logger.info("Initializing pipeline mp log thread")

        self._log_lock = log_lock
        self._log_queue = log_q
        self._stop_event = threading.Event()

    def run(self):
        while True:
            try:
                if self._stop_event.is_set():
                    break

                try:
                    with self._log_lock:
                        level, log = self._log_queue.get_nowait()
                    func = getattr(logger, level)
                    func(log)
                except queue.Empty:
                    log = None

                if log is None:
                    time.sleep(0.1)
            except Exception:
                logger.error('Error in log thread:\n{}'.format(traceback.format_exc()))
                time.sleep(0.1)

        logger.info("Quitting pipeline mp log thread")

    def stop(self):
        """Stops the thread cleanly."""
        self._stop_event.set()

class Experiment(object):

    def __init__(self, exp_name, exp_type, data_dir, fprefix, **kwargs):
        self.exp_name = exp_name #Should be unique identifier string for the experiment
        self.exp_type = exp_type #e.g. 'SEC' or 'Batch'
        self.data_dir = data_dir

        self.current_fprefix = fprefix

        self.output_dir = None
        self.profiles_dir = None
        self.analysis_dir = None

        self.profiles = []

        self.collection_finished = False
        self.analysis_finished = False
        self.exp_last_modified = -1
        self.analysis_last_modified = -1

        self.sub_profile = None
        self.ift = None
        self.dammif_data = None
        self.denss_data = None
        self.series = None

        if self.exp_type == 'SEC':
            self.num_exps = kwargs['num_exps']

        elif self.exp_type == 'Batch':
            self.buffer_profiles = []
            self.num_sample_exps = kwargs['num_sample_exps']
            self.num_buffer_exps = kwargs['num_buffer_exps']
            self.sample_prefix = kwargs['sample_prefix']
            self.buffer_prefix = kwargs['buffer_prefix']

    def add_profiles(self, new_profiles):
        if not self.collection_finished:
            if self.exp_type == 'SEC':
                self.profiles.extend(new_profiles)

                if len(self.profiles) == self.num_exps:
                    self.collection_finished = True

                self.exp_last_modified = time.time()

            elif self.exp_type == 'Batch':
                for prof in new_profiles:
                    if prof.getParameter('filename').startswith(self.sample_prefix):
                        self.profiles.append(prof)
                    elif prof.getParameter('filename').startswith(self.buffer_prefix):
                        self.buffer_profiles.append(prof)

                if (len(self.profiles) == self.num_sample_exps
                    and len(self.buffer_profiles) == self.num_buffer_exps):
                    self.collection_finished = True

                self.exp_last_modified = time.time()

    def check_exp_timeout(self, pl_settings):
        if self.exp_type == 'SEC':
            timeout = pl_settings['sec_exp_timeout']
        elif self.exp_type == 'Batch':
            timeout = pl_settings['batch_exp_timeout']

        if timeout > 0 and self.exp_last_modified > 0:
            if time.time() - self.exp_last_modified > timeout:
                logger.error('Experiment %s data collection timed out', self.exp_name)
                self.collection_finished = True

    def check_analysis_timeout(self, pl_settings):
        if self.exp_type == 'SEC':
            timeout = pl_settings['sec_analysis_timeout']
        elif self.exp_type == 'Batch':
            timeout = pl_settings['batch_analysis_timeout']

        if timeout > 0 and self.analysis_last_modified > 0:
            if time.time() - self.analysis_last_modified > timeout:
                logger.error('Experiment %s analysis timed out', self.exp_name)
                self.analysis_finished = True

    def stop_experiment(self):
        if self.exp_last_modified > 0:
            self.collection_finished = True

        else:
            self.exp_last_modified = time.time()

if __name__ == '__main__':
    mp.set_start_method('spawn')
