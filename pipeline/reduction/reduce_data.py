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

import os.path
import multiprocessing as mp
import queue
import threading
import collections
import time
import logging
import traceback
import requests
import json

if __name__ != '__main__':
    logger = logging.getLogger(__name__)

# import watchdog.utils.dirsnapshot as dirsnapshot

import bioxtasraw.RAWAPI as raw
import bioxtasraw.SASExceptions as SASExceptions

pipeline_path = os.path.abspath(os.path.join('.', __file__, '..', '..', '..'))
if pipeline_path not in os.sys.path:
    os.sys.path.append(pipeline_path)

from pipeline.reduction import eiger_stream_client

def load_images(filenames, settings):
    img_list, imghdr_list = raw.load_images(filenames, settings)

    return img_list, imghdr_list, filenames

def load_counter_values(filenames, settings, new_filenames=[]):
    """
    Note: when filenames are provided, this expects the pseudo-filenames that
    RAW generates for files with multiple images. So if you loaded 10 images
    from my_data.hdf5, you'd have my_data_00000.hdf5, my_data_00001.hdf5, etc.
    The base image would need to be repeated 10 times in the filenames list.
    Also, if you provide new_filenames, every filename in the filenames list
    should have a corresponding entry in the list ('' if no new_filename, such
    as for a .tiff file).
    e.g. filenames = ['my_img.tiff', 'my_img2.hdf5', 'my_img2.hdf5']
    new_filenames = ['', 'my_img2_00000.hdf5', my_img2_00001.hdf5']
    """
    counter_list = raw.load_counter_values(filenames, settings, new_filenames)

    return counter_list

def load_images_and_counters(filenames, settings):

    img_list = []
    img_hdr_list = []
    counter_list = []
    filename_list = []
    failed = []

    for fname in filenames:
        imgs, img_hdrs = raw.load_images([fname,], settings)

        if all(map(lambda img: img is not None, imgs)):
            img_list.extend(imgs)
            img_hdr_list.extend(img_hdrs)

            new_fnames = []

            for i in range(len(imgs)):
                if len(imgs) > 1:
                    temp_fname = os.path.split(fname)[1].split('.')
                    if len(temp_fname) > 1:
                        temp_fname[-2] = temp_fname[-2] + '_%05i' %(i+1)
                    else:
                        temp_fname[0] = temp_fname[0] + '_%05i' %(i+1)

                    new_fname = '.'.join(temp_fname)
                else:
                    new_fname = os.path.split(fname)[1]


                new_fnames.append(new_fname)

            counters = raw.load_counter_values([fname for n in new_fnames],
                settings, new_fnames)

            counter_list.extend(counters)

            filename_list.extend(new_fnames)

        else:
            failed.append(fname)

    return img_list, img_hdr_list, counter_list, filename_list, failed

class monitor_and_load(threading.Thread):

    def __init__(self, cmd_q, ret_q, abort_event, raw_settings,
        pipeline_settings, log_lock=threading.Lock(), log_queue=queue.Queue()):
        threading.Thread.__init__(self)
        self.daemon = True
        self.name = 'monitor_and_load_thread'

        if mp.current_process().name == 'MainProcess':
            self.in_main = True
        else:
            self.in_main = False

        self._log_lock = log_lock
        self._log_queue = log_queue

        self._log('debug', "Initializing pipeline monitor and load thread")

        self._cmd_q = cmd_q
        self._ret_q = ret_q
        self._abort_event = abort_event
        self._stop_event = threading.Event()

        self.raw_settings = raw_settings
        self.pl_settings = pipeline_settings

        self._monitor_cmd_q = collections.deque()
        self._monitor_ret_q = collections.deque()
        self._monitor_abort = threading.Event()

        self._commands = {'set_data_dir'    : self._set_data_dir,
            'set_fprefix': self._set_fprefix,
            'set_data_dir_and_fprefix': self._set_data_dir_fprefix,
            'set_ret_size'  : self._set_ret_size,
            'set_raw_settings'  : self._set_raw_settings,
            'set_experiment'    : self._set_experiment,
            }

        self.waiting_for_header = []
        self.header_timeout = 10
        self.ret_every = 10

        self._exp_id = None
        self.data_dir = None
        self._fprefix = None

    def run(self):
        if self.pl_settings['data_source'] == 'Files':
            self._log('debug', 'Starting file monitor')
            self._monitor_thread = monitor_thread(self._monitor_cmd_q,
                self._monitor_ret_q, self._monitor_abort, self.pl_settings,
                self._log_lock, self._log_queue)

        else:
            self._log('debug', "Starting Eiger Stream monitor")
            requests.put('http://{}/stream/api/1.8.0/config/mode'.format(
                self.pl_settings['eiger_stream_ip']), data='{"value": "enabled"}')

            self._monitor_thread = eiger_stream_client.EigerStreamClient(
                self.pl_settings['eiger_stream_ip'],
                self.pl_settings['eiger_stream_port'],
                self._monitor_ret_q)

            self._eiger_parser = eiger_stream_client.EigerStreamParser()

        self._monitor_thread.start()

        while True:
            try:
                if self._stop_event.is_set():
                    break

                try:
                    cmd, args = self._cmd_q.popleft()
                except IndexError:
                    cmd = None

                if cmd is not None:
                    self._log('info', "Processing cmd '{}' with args: {}".format(cmd,
                        ', '.join(['{}'.format(a) for a in args])))
                    self._commands[cmd](*args)

                if self._abort_event.is_set():
                    self._abort()

                new_failed = []
                new_succeded = []
                got_new_images = False

                for item in self.waiting_for_header:
                    if self._abort_event.is_set():
                        break

                    img, initial_t, exp_id, data_dir = item

                    if exp_id is None:
                        exp_id = self._exp_id

                    if data_dir is None:
                        data_dir = self.data_dir

                    if time.time() - initial_t > self.header_timeout:
                        new_failed.append(item)
                        self._log('info', 'Failed to load image')
                    else:
                        imgs, fnames, img_hdrs, counters = self._load_images(img)

                        got_new_images = True

                        if len(imgs) > 0:
                            new_succeded.append(item)

                            self._ret_q.append([imgs, fnames, img_hdrs, counters,
                                exp_id, data_dir])

                if self._abort_event.is_set():
                    self._abort()
                    new_failed = []
                    new_succeded = []

                for item in new_failed:
                    img, initial_t, exp_id, data_dir = item

                    for i in range(len(self.waiting_for_header)):
                        item2 = self.waiting_for_header[i]
                        img2, initial_t2, exp_id2, data_dir2 = item2

                        if initial_t == initial_t2 and exp_id == exp_id2:
                            self.waiting_for_header.pop(i)
                            break

                    # if item in self.waiting_for_header:
                    #     self.waiting_for_header.remove(item)

                for item in new_succeded:
                    img, initial_t, exp_id, data_dir = item

                    for i in range(len(self.waiting_for_header)):
                        item2 = self.waiting_for_header[i]
                        img2, initial_t2, exp_id2, data_dir2 = item2

                        if initial_t == initial_t2 and exp_id == exp_id2:
                            self.waiting_for_header.pop(i)
                            break

                try:
                    new_images = []
                    while len(self._monitor_ret_q) > 0:
                        img_data = self._monitor_ret_q.popleft()
                        new_images.append(img_data)
                except IndexError:
                    new_images = []

                for img_data in new_images:
                    new_imgs = []
                    new_img_hdrs = []
                    new_counters = []
                    new_fnames = []

                    if self.pl_settings['data_source'] == 'Files':
                        exp_id = img_data[0]
                        data_dir = img_data[1]
                        new_img_data = img_data[2]

                        if exp_id is None:
                            exp_id = self._exp_id

                        if data_dir is None:
                            data_dir = self.data_dir

                    else:
                        exp_id = self._exp_id
                        data_dir = self.data_dir
                        nimg, header = self._eiger_parser.decodeFrames(img_data)

                        if nimg is not None:
                            try:
                                header2 = json.loads(img_data[3].bytes)
                                header.update(header2)
                            except Exception:
                                logger.error('Error parsing Eiger header:\n{}'.format(traceback.format_exc()))

                            # Check md5 hash of image to verify no corruption?
                            md5_hash = header['hash']

                            new_img_data = [[nimg, header]]
                        else:
                            new_img_data = []


                    for img in new_img_data:
                        if self._abort_event.is_set():
                            break

                        imgs, fnames, img_hdrs, counters = self._load_images(img, new_image=True)

                        if len(imgs) > 0:
                            new_imgs.extend(imgs)
                            new_img_hdrs.extend(img_hdrs)
                            new_counters.extend(counters)
                            new_fnames.extend(fnames)

                            got_new_images = True

                        else:
                            self.waiting_for_header.append((img, time.time(), exp_id, data_dir))

                        if len(new_imgs) >= self.ret_every:
                            self._ret_q.append([new_imgs, new_fnames, new_img_hdrs,
                                new_counters, exp_id, data_dir])
                            new_imgs = []
                            new_img_hdrs = []
                            new_counters = []
                            new_fnames = []

                    if new_imgs:
                        self._ret_q.append([new_imgs, new_fnames, new_img_hdrs,
                            new_counters, exp_id, data_dir])

                if self._abort_event.is_set():
                    self._abort()
                    new_failed = []
                    new_succeded = []

                if not got_new_images:
                    time.sleep(0.1)

            except Exception:
                self._log('error', "Error in monitor and load thread:\n{}".format(traceback.format_exc()))

        self._log('debug', "Quitting pipeline monitor_and_load thread")

    def _set_data_dir(self, data_dir):
        self._log('debug', 'Setting data directory: {}'.format(data_dir))

        if self.pl_settings['data_source'] == 'Files':
            self._monitor_cmd_q.append(['set_data_dir', [data_dir,]])

        self.data_dir = data_dir

    def _set_fprefix(self, fprefix):
        if self.pl_settings['data_source'] == 'Files':
            self._monitor_cmd_q.append(['set_fprefix', [fprefix,]])

        self._fprefix = fprefix

    def _set_data_dir_fprefix(self, data_dir, fprefix):
        self._log('debug', 'Setting data directory: {}'.format(data_dir))

        if self.pl_settings['data_source'] == 'Files':
            self._monitor_cmd_q.append(['set_data_dir_and_fprefix', [data_dir, fprefix]])

        self.data_dir = data_dir
        self._fprefix = fprefix

    def _set_experiment(self, data_dir, fprefix, exp_id):
        self._log('debug', 'Setting experiment to {}, data directory to {}'.format(exp_id, data_dir))
        self._exp_id = exp_id
        self.data_dir = data_dir
        self._fprefix = fprefix

        if self.pl_settings['data_source'] == 'Files':
            self._monitor_cmd_q.append(['set_experiment', [data_dir, fprefix, exp_id]])

    def _set_ret_size(self, size):
        self._log('debug', 'Setting return chunk size to {}'.format(size))
        self.ret_every = int(size)

    def _load_image_and_counter(self, img_name):

        file_size = -1
        while file_size != os.path.getsize(img_name):
            file_size = os.path.getsize(img_name)

        try:
            imgs, img_hdrs, counters, fnames, failed = load_images_and_counters([img_name,],
                self.raw_settings)

            self._log('debug', 'Loaded image {}'.format(img_name))

        except SASExceptions.HeaderLoadError:
            imgs = img_hdrs = counters = fnames = []

        return imgs, fnames, img_hdrs, counters

    def _load_counter(self, filenames, new_filenames):
        try:
            counters = load_counter_values(filenames, self.raw_settings, new_filenames)
        except SASExceptions.HeaderLoadError:
            counters = [{}]

        return counters

    def _load_images(self, img, new_image=False):
        if self.pl_settings['data_source'] == 'Files':
            imgs, fnames, img_hdrs, counters  = self._load_image_and_counter(img)

        else:
            imgs = [img[0]]
            img_hdrs = [img[1]]

            if img_hdrs[0] is not None:
                frame = int(img_hdrs[0]['frame']) + 1

            else:
                frame = 1


            fname = '{}_{:06}.h5'.format(self._fprefix, frame)
            fnames = [fname]
            ctr_base_fnames = [os.path.join(self.data_dir,
                '{}_data_000001.h5'.format(self._fprefix))]

            if new_image:
                self._log('debug', 'Received image {} from Eiger'.format(fname))

            counters = self._load_counter(ctr_base_fnames, fnames)

        if len(counters) == 0 or len(counters[0]) == 0:
            imgs = []
            fnames = []
            img_hdrs = []

        return imgs, fnames, img_hdrs, counters

    def _set_raw_settings(self, settings):
        self._log('debug', 'Setting RAW settings')

        self.raw_settings = settings

    def _abort(self):
        self._monitor_abort.set()

        self._cmd_q.clear()
        self._ret_q.clear()

        self._monitor_cmd_q.clear()
        self._monitor_ret_q.clear()

        self.waiting_for_header = []

        self._abort_event.clear()

    def _log(self, level, msg):
        if self.in_main:
            func = getattr(logger, level)
            func(msg)

        else:
            with self._log_lock:
                self._log_queue.put_nowait((level, '{} - {}'.format(self.name, msg)))

    def stop(self):
        """Stops the process cleanly."""
        self._monitor_thread.stop()
        self._monitor_thread.join()

        self._stop_event.set()

class monitor_thread(threading.Thread):

    def __init__(self, cmd_q, ret_q, abort_event, pipeline_settings,
        log_lock, log_queue):
        threading.Thread.__init__(self)
        self.daemon = True
        self.name = 'monitor_thread'

        if mp.current_process().name == 'MainProcess':
            self.in_main = True
        else:
            self.in_main = False

        self._log_lock = log_lock
        self._log_queue = log_queue

        self._log('debug', "Initializing pipeline monitor thread")

        self._cmd_q = cmd_q
        self._ret_q = ret_q
        self._abort_event = abort_event
        self._stop_event = threading.Event()

        self.pl_settings = pipeline_settings
        self.fprefix = ''

        self._exp_id = None
        self.data_dir = None

        self._commands = {'set_data_dir': self._set_data_dir,
            'set_fprefix': self._set_fprefix,
            'set_data_dir_and_fprefix': self._set_data_dir_fprefix,
            'set_experiment': self._set_experiment,
            }

    def run(self):
        self.dir_snapshot = []

        while True:
            try:
                if self._stop_event.is_set():
                    break

                if self._abort_event.is_set():
                    self._abort()

                try:
                    cmd, args = self._cmd_q.popleft()
                except IndexError:
                    cmd = None

                if cmd is not None:
                    self._log('info', "Processing cmd '{}' with args: {}".format(cmd,
                        ', '.join(['{}'.format(a) for a in args])))

                    self._commands[cmd](*args)

                if self.data_dir is not None and os.path.exists(self.data_dir):
                    new_images = self._check_for_new_files()

                    if self._abort_event.is_set():
                        new_images = []
                        self._abort()

                    if new_images:
                        self._ret_q.append([self._exp_id, self.data_dir, new_images])
                    else:
                        time.sleep(0.1)
                else:
                    time.sleep(0.1)

            except Exception:
                self._log('error', "Error in monitor thread:\n{}".format(traceback.format_exc()))

        self._log('debug', "Quitting pipeline monitor thread")

    def _set_data_dir(self, data_dir):
        self._log('debug', 'Setting data directory: {}'.format(data_dir))

        self.data_dir = data_dir

        self.dir_snapshot = []

    def _set_fprefix(self, fprefix):
        self._log('debug', 'Setting file prefix: {}'.format(fprefix))

        self.fprefix = fprefix

    def _set_data_dir_fprefix(self, data_dir, fprefix):
        self._set_fprefix(fprefix)
        self._set_data_dir(data_dir)

    def _set_experiment(self, data_dir, fprefix, exp_id):
        self._log('debug', 'Setting experiment: {}'.format(exp_id))

        self._exp_id = exp_id
        self._set_data_dir_fprefix(data_dir, fprefix)

    def _check_for_new_files(self):

        new_files = []

        with os.scandir(self.data_dir) as files:
            for f in files:
                if f.is_file():
                    f_path = f.path
                    if f_path not in self.dir_snapshot:
                        new_files.append(f_path)
                        self.dir_snapshot.append(f_path)

        new_images = []

        for f in new_files:
            if self._abort_event.is_set():
                new_images = []
                self._abort()

            if (os.path.splitext(f)[1] in self.pl_settings['image_exts']
                and os.path.split(f)[1].startswith(self.fprefix)):
                new_images.append(os.path.abspath(os.path.expanduser(f)))

        if new_images:
            self._log('info', 'New images found: {}'.format(', '.join(new_images)))

        return new_images

    def _abort(self):
        self._cmd_q.clear()
        self._ret_q.clear()

        if self.data_dir is not None and os.path.exists(self.data_dir):
            self.dir_snapshot = [f.path for f in os.scandir(self.data_dir) if f.is_file()]
        else:
            self.dir_snapshot = []

        self._abort_event.clear()

    def _log(self, level, msg):
        if self.in_main:
            func = getattr(logger, level)
            func(msg)

        else:
            with self._log_lock:
                self._log_queue.put_nowait((level, '{} - {}'.format(self.name, msg)))

    def stop(self):
        """Stops the thread cleanly."""
        self._stop_event.set()


class raver_process(mp.Process):

    def __init__(self, cmd_q, ret_q, cmd_lock, ret_lock, abort_event,
        raw_settings_file, pipeline_settings, log_lock, log_queue,
        managed_internal_q, num_loaded, num_averaged):
        mp.Process.__init__(self)
        self.daemon = True

        self._cmd_q = cmd_q

        self._internal_cmd_q = managed_internal_q
        self._ret_q = ret_q
        self._cmd_lock = cmd_lock
        self._ret_lock = ret_lock
        self._abort_event = abort_event
        self._stop_event = mp.Event()

        self._log_lock = log_lock
        self._log_queue = log_queue
        self.raw_settings_file = raw_settings_file

        self.pl_settings = pipeline_settings

        self._commands = {'raver_images': self._raver_images,
            'set_data_dir'    : self._set_data_dir,
            'set_fprefix'       : self._set_fprefix,
            'set_data_dir_and_fprefix'  : self._set_data_dir_and_fprefix,
            'set_experiment'    : self._set_experiment,
            'load_settings'     : self._load_settings,
            }

        self.ret_every = 10

        self.data_dir = None

        self._exp_id = ''
        self.fprefix = ''

        self.num_loaded = num_loaded
        self.num_averaged = num_averaged

        self.active = False

    def run(self):
        self._log('debug', "Starting pipeline radial averaging process")

        self.raw_settings = raw.load_settings(self.raw_settings_file)

        #Initialize monitoring thread
        self.m_cmd_q = collections.deque()
        self.m_ret_q = collections.deque()
        self.m_abort_event = threading.Event()

        self.monitor_thread = monitor_and_load(self.m_cmd_q,
            self.m_ret_q, self.m_abort_event, self.raw_settings,
            self.pl_settings, self._log_lock, self._log_queue)

        self.monitor_thread.start()

        while True:
            self.active = False

            try:
                if self._stop_event.is_set():
                    break

                if self._abort_event.is_set():
                    self._abort()

                try:
                    with self._cmd_lock:
                        cmd, args, kwargs = self._cmd_q.get_nowait()
                except queue.Empty:
                    try:
                        cmd, args, kwargs = self._internal_cmd_q.get_nowait()
                    except queue.Empty:
                        cmd = None

                if cmd is not None:
                    self._log('debug', "Processing cmd '%s'" %(cmd))
                    self.active = True

                    self._commands[cmd](*args, **kwargs)

                try:
                    img_data = self.m_ret_q.popleft()
                except IndexError:
                    img_data = None

                if img_data is not None:
                    exp_id = img_data[-2]
                    data_dir = img_data[-1]
                    img_data = img_data[:-2]

                    if exp_id is None:
                        exp_id = self._exp_id

                    if data_dir is None:
                        data_dir = self.data_dir

                    self.num_loaded.set(self.num_loaded.value + len(img_data[0]))

                    self.active = True

                    self._raver_images(exp_id, img_data, {'load_path': data_dir})

                if not self.active:
                    time.sleep(0.1)

            except Exception:
                self._log('error', "Error in raver process:\n{}".format(traceback.format_exc()))

            except KeyboardInterrupt:
                self._abort()
                break

        self.monitor_thread.stop()
        self.monitor_thread.join()

        self._log('debug', "Quiting reduction process")

    def _set_data_dir(self, data_dir):
        self._log('debug', 'Setting data directory: %s' %(data_dir))

        self.data_dir = os.path.abspath(os.path.expanduser(data_dir))

        self.m_cmd_q.append(['set_data_dir', [self.data_dir]])

    def _set_fprefix(self, fprefix):
        self._log('debug', 'Setting file prefix: %s' %fprefix)

        self.fprefix = fprefix

        self.m_cmd_q.append(['set_fprefix', [self.fprefix]])

    def _set_data_dir_and_fprefix(self, data_dir, fprefix):
        self._log('debug', 'Setting data directory and file prefix: %s, %s' %(data_dir,
            fprefix))

        self.data_dir = os.path.abspath(os.path.expanduser(data_dir))
        self.fprefix = fprefix

        self.m_cmd_q.append(['set_data_dir_and_fprefix', [self.data_dir,
            self.fprefix]])

    def _set_experiment(self, data_dir, fprefix, exp_id):
        self._log('debug', 'Setting experiment: %s' %(exp_id))

        self._exp_id = exp_id

        self.data_dir = os.path.abspath(os.path.expanduser(data_dir))
        self.fprefix = fprefix

        self.m_cmd_q.append(['set_experiment', [self.data_dir,
            self.fprefix, self._exp_id]])

    def _raver_images(self, *args, **kwargs):
        self._log('debug', 'Radially averaging images')
        img_hdrs = []
        all_counters = []
        load_path = ''

        exp_id = args[0]
        img_data = args[1]
        imgs = img_data[0]
        if len(img_data) >1:
            names = img_data[1]
        if len(img_data) > 2:
            img_hdrs = img_data[2]
        if len(img_data) > 3:
            all_counters = img_data[3]

        if 'img_hdrs' in kwargs:
            img_hdrs = kwargs['img_hdrs']
        if 'all_counters' in kwargs:
            all_counters = kwargs['all_counters']
        if 'names' in kwargs:
            names = kwargs['names']
        if 'load_path' in kwargs:
            load_path = kwargs['load_path']

        profiles = []

        for i, image in enumerate(imgs):
            if self._abort_event.is_set():
                self._abort()
                break

            try:
                img_hdr = img_hdrs[i]
            except Exception:
                img_hdr = {}
            try:
                counters = all_counters[i]
            except Exception:
                counters = {}

            name = names[i]

            profile = raw.integrate_image(image, self.raw_settings, name,
                img_hdr=img_hdr, counters=counters, load_path=load_path)

            profiles.append(profile)

            if len(profiles) >= self.ret_every:
                self.num_averaged.set(self.num_averaged.value + len(profiles))

                with self._ret_lock:
                    self._ret_q.put_nowait([exp_id, profiles])
                profiles = []

        if len(profiles) > 0:
            self.num_averaged.set(self.num_averaged.value + len(profiles))

            with self._ret_lock:
                self._ret_q.put_nowait([exp_id, profiles])

    def load_settings(self, settings_file):
        self._internal_cmd_q.put_nowait(['load_settings', [settings_file], {}])

    def _load_settings(self, settings_file):
        self._log('debug', 'Loading RAW settings from %s' %(settings_file))
        self.raw_settings_file = settings_file
        self.raw_settings = raw.load_settings(settings_file)

        self.m_cmd_q.append(['set_raw_settings', [self.raw_settings]])

    def _log(self, level, msg):
        with self._log_lock:
            self._log_queue.put_nowait((level, '{} - {}'.format(self.name, msg)))

    def _abort(self):
        while True:
            try:
                with self._cmd_lock:
                    cmd, args, kwargs = self._cmd_q.get_nowait()
            except queue.Empty:
                break

        with self._ret_lock:
            self._ret_q.put_nowait('aborted')

        while self._abort_event.is_set():
            time.sleep(0.1)

        self.m_abort_event.set()

    def stop(self):
        """Stops the thread cleanly."""

        self._stop_event.set()


class save_thread(threading.Thread):
    #It's possible this would be more efficient as it's own process, but I'll
    #leave it here as a thread for now, and do some testing if this is a limitation

    def __init__(self, cmd_q, ret_q, abort_event, raw_settings):
        threading.Thread.__init__(self)
        self.daemon = True
        self.name = 'save_thread'

        logger.debug("Initializing pipeline save profiles thread")

        self._cmd_q = cmd_q
        self._ret_q = ret_q
        self._abort_event = abort_event
        self._stop_event = threading.Event()

        self.raw_settings = raw_settings

        self._commands = {'save_profiles': self._save_profiles,
            'set_raw_settings'  : self._set_raw_settings,
            }

    def run(self):
        while True:
            if self._stop_event.is_set():
                break

            if self._abort_event.is_set():
                self._abort()

            try:
                cmd, args = self._cmd_q.popleft()
            except IndexError:
                cmd = None

            if cmd is not None:
                logger.debug("Processing cmd '%s'", cmd)
                self._commands[cmd](*args)

            else:
                time.sleep(0.1)

        logger.debug("Quitting save thread")

    def _save_profiles(self, output_dir, profiles):
        logger.debug('Saving profiles')

        for profile in profiles:
            raw.save_profile(profile, datadir=output_dir,
                settings=self.raw_settings)

        logger.debug('Finished saving profiles')

    def _set_raw_settings(self, settings):
        logger.debug('Setting RAW settings')
        self.raw_settings = settings

    def _abort(self):
        self._cmd_q.clear()
        self._ret_q.clear()

        self._abort_event.clear()

    def stop(self):
        """Stops the thread cleanly."""
        self._stop_event.set()

