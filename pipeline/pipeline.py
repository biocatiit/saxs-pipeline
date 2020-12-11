import threading
import time
import multiprocessing as mp
import queue
import os
import collections
import copy

import bioxtasraw.RAWAPI as raw

from .analysis import analyze
from .reports import pdf
from .reduction import reduce_data


class pipeline_thread(threading.Thread):

    def __init__(self, cmd_q, ret_q, abort_event, pipeline_settings,
        raw_settings_file):
        threading.Thread.__init__(self)
        self.daemon = True

        self._cmd_q = cmd_q
        self._ret_q = ret_q
        self._abort_event = abort_event
        self._stop_event = threading.Event()

        self.pl_settings = pipeline_settings
        self.raw_settings_file = raw_settings_file
        self.raw_settings = raw.load_settings(self.raw_settings_file)

        self.data_dir = None
        self.output_dir = None
        self.profiles_dir = None
        self.analysis_dir = None

        self.fprefix = ''

        self.experiments = {}
        self.current_experiment = ''

        self.active = False

        self._set_analysis_args()

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
            }

        #Initialize monitoring thread
        self.m_cmd_q = collections.deque()
        self.m_ret_q = collections.deque()
        self.m_abort_event = threading.Event()

        self.monitor_thread = reduce_data.monitor_and_load(self.m_cmd_q,
            self.m_ret_q, self.m_abort_event, self.raw_settings,
            self.pl_settings)

        self.monitor_thread.start()

        #Initialize save thread
        self.s_cmd_q = collections.deque()
        self.s_ret_q = collections.deque()
        self.s_abort_event = threading.Event()

        self.save_thread = reduce_data.save_thread(self.s_cmd_q, self.s_ret_q,
            self.s_abort_event, self.raw_settings)

        self.save_thread.start()

        #Initialize reduction and analysis threads
        self.reduction_processes = []
        self.analysis_processes = []

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

        for i in range(self.pl_settings['r_procs']):
            self._start_reduction_process()

        for i in range(self.pl_settings['a_procs']):
            self._start_analysis_process()

    def run(self):
        while True:
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
                print(cmd)
                self._commands[cmd](*args, **kwargs)

            # Need to wrap this and the ret_q in a while, and get all things in the q?
            try:
                img_data = self.m_ret_q.popleft()
            except IndexError:
                img_data = None

            if img_data is not None:
                exp_id = img_data[-2]
                data_dir = img_data[-1]
                img_data = img_data[:-2]

                if exp_id is None:
                    exp_id = self.current_experiment

                if data_dir is None:
                    data_dir = self.data_dir

                self.active = True
                with self.r_cmd_lock:
                    self.r_cmd_q.put_nowait(['raver_images', exp_id, img_data,
                        {'load_path': data_dir}])

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

    def _start_reduction_process(self):
        proc = reduce_data.raver_process(self.r_cmd_q, self.r_ret_q,
            self.r_cmd_lock, self.r_ret_lock, self.r_abort_event,
            self.raw_settings_file)

        proc.start()

        self.reduction_processes.append(proc)

    def _start_analysis_process(self):
        proc = analyze.analysis_process(self.a_cmd_q, self.a_ret_q,
            self.a_cmd_lock, self.a_ret_lock, self.a_abort_event,
            self.raw_settings_file)

        proc.start()

        self.analysis_processes.append(proc)

    def _set_data_dir(self, data_dir):
        self.data_dir = os.path.abspath(os.path.expanduser(data_dir))

        if self.pl_settings['use_default_output_dir']:
            output_dir = os.path.join(data_dir,
                self.pl_settings['default_output_dir'])

            self._set_output_dir(output_dir)

        self.m_cmd_q.append(['set_data_dir', [self.data_dir]])

    def _set_fprefix(self, fprefix):
        self.fprefix = fprefix

        self.m_cmd_q.append(['set_fprefix', [self.fprefix]])

    def _set_data_dir_and_fprefix(self, data_dir, fprefix):
        self.data_dir = os.path.abspath(os.path.expanduser(data_dir))
        self.fprefix = fprefix

        if self.pl_settings['use_default_output_dir']:
            output_dir = os.path.join(data_dir,
                self.pl_settings['default_output_dir'])

            self._set_output_dir(output_dir)

        self.m_cmd_q.append(['set_data_dir_and_fprefix', [self.data_dir,
            self.fprefix]])

    def _set_output_dir(self, output_dir):
        self.output_dir = os.path.abspath(os.path.expanduser(output_dir))

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
        self.profiles_dir = os.path.abspath(os.path.expanduser(profiles_dir))

        if self.current_experiment in self.experiments:
            self.experiments[self.current_experiment].profiles_dir = profiles_dir

        if not os.path.exists(self.profiles_dir):
            os.mkdir(self.profiles_dir)

    def _set_analysis_dir(self, analysis_dir):
        self.analysis_dir = os.path.abspath(os.path.expanduser(analysis_dir))

        if self.current_experiment in self.experiments:
            self.experiments[self.current_experiment].analysis_dir = analysis_dir

        if not os.path.exists(self.analysis_dir):
            os.mkdir(self.analysis_dir)

    def _load_raw_settings(self, raw_settings_file):
        self.raw_settings = raw.load_settings(raw_settings_file)

        self.m_cmd_q.append(['set_raw_settings', [self.raw_settings]])
        self.s_cmd_q.append(['set_raw_settings', [self.raw_settings]])

        for proc in self.reduction_processes:
            proc.load_settings(raw_settings_file)

        for proc in self.analysis_processes:
            proc.load_settings(raw_settings_file)

    def _start_experiment(self, *args, **kwargs):
        print('starting experiment')
        exp_name = args[0]
        exp_type = args[1]
        data_dir = args[2]
        fprefix = args[3]

        if exp_name not in self.experiments:
            new_exp = Experiment(exp_name, exp_type, data_dir, fprefix,
                **kwargs)

            self.experiments[exp_name] = new_exp

        else:
            self.expeirments[exp_name].current_fprefix = fprefix

        self.current_experiment = exp_name

        self.data_dir = os.path.abspath(os.path.expanduser(data_dir))
        self.fprefix = fprefix

        if self.pl_settings['use_default_output_dir']:
            output_dir = os.path.join(data_dir,
                self.pl_settings['default_output_dir'])

            self._set_output_dir(output_dir)

        self.m_cmd_q.append(['set_experiment', [self.data_dir,
            self.fprefix, self.current_experiment]])

    def _stop_experiment(self, exp_name):
        if exp_name in self.expeirments:
            self.experiments[exp_name].collection_finished = True

    def _save_profiles(self, profile_data):
        print('saving profiles')
        exp_id, profiles = profile_data
        print(exp_id)

        if self.pl_settings['save_raver_profiles']:
            if exp_id in self.experiments:
                profiles_dir = self.experiments[exp_id].profiles_dir
            else:
                profiles_dir = self.profiles_dir

            for profile in profiles:
                print(profile.getParameter('filename'))
            print(profiles_dir)
            if profiles_dir is not None:
                self.s_cmd_q.append(['save_profiles', [profiles_dir, profiles]])

    def _add_profiles_to_experiment(self, profile_data):
        print('adding profiles to experiment')
        exp_id, profiles = profile_data

        if exp_id in self.experiments:
            self.experiments[exp_id].add_profiles(profiles)

    def _check_exp_status(self):

        save_proc_data = self.pl_settings['save_processed_data']
        save_report = self.pl_settings['save_report']
        report_type = self.pl_settings['report_type']

        for exp in self.experiments.values():
            if not exp.collection_finished:
                exp.check_exp_timeout(self.pl_settings)

            if exp.collection_finished and exp.analysis_last_modified == -1:
                self.active = True
                print('experiment finished')
                with self.a_cmd_lock:
                    if exp.exp_type == 'SEC':
                        self.a_cmd_q.put_nowait(['make_and_analyze_series',
                            exp.exp_name, [exp.analysis_dir, exp.profiles,
                            save_proc_data, save_report, report_type,
                            exp.output_dir], self._analysis_args])
                        exp.analysis_last_modified = time.time()

                    elif exp.exp_type == 'Batch':
                        pass
                        #batch processing goes here, not developed yet

    def _add_analysis_to_experiment(self, results):
        print('adding analysis to experiment')
        res_type = results[0]
        exp_id = results[1]

        if exp_id in self.experiments:
            exp = self.experiments[exp_id]

            if res_type == 'sub_series':
                series = results[2]

                exp.series = series

            elif res_type == 'analysis_results':
                res_dict = results[2]

                exp.sub_profile = res_dict['profile']
                exp.ift = res_dict['ift']
                exp.dammif_data = res_dict['dammif_data']
                exp.denss_data = res_dict['denss_data']

                exp.analysis_last_modified = time.time()
                exp.analysis_finished = True

    def _check_analysis_status(self):
        exp_keys = list(self.experiments.keys())

        for exp_key in exp_keys:
            exp = self.experiments[exp_key]

            if not exp.analysis_finished:
                exp.check_analysis_timeout(self.pl_settings)

            if exp.analysis_finished:
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

    def _update_pipeline_settings(self, *args, **kwargs):
        for key in kwargs:
            if key in self.pl_settings:
                self.pl_settings[key] = kwargs[key]

        self._set_analysis_args()

    def _abort(self):
        self._cmd_q.clear()
        self._ret_q.clear()

        self.r_abort_event.set()
        self.a_abort_event.set()

        self._abort_event.clear()

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
                self.collection_finished = True

    def check_analysis_timeout(self, pl_settings):
        if self.exp_type == 'SEC':
            timeout = pl_settings['sec_analysis_timeout']
        elif self.exp_type == 'Batch':
            timeout = pl_settings['batch_analysis_timeout']

        if timeout > 0 and self.analysis_last_modified > 0:
            if time.time() - self.analysis_last_modified > timeout:
                self.analysis_finished = True

if __name__ == '__main__':
    mp.set_start_method('spawn')


"""
Todos:
3) Need to think about whether to wrap the calls to all the return queues to check for multiple items, so it moves everything along at once
6) Batch mode processing!
"""
