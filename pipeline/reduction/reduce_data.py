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
import multiprocessing
import queue
import threading
import collections
import time
import logging

if __name__ != '__main__':
    logger = logging.getLogger(__name__)

import watchdog.utils.dirsnapshot as dirsnapshot

import bioxtasraw.RAWAPI as raw
import bioxtasraw.SASExceptions as SASExceptions

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
        imgs, img_hdrs = raw.load_images(filenames, settings)

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
        pipeline_settings):
        threading.Thread.__init__(self)
        self.daemon = True
        self.name = 'monitor_and_load_thread'

        logger.debug("Initializing pipeline monitor and load thread")

        self._cmd_q = cmd_q
        self._ret_q = ret_q
        self._abort_event = abort_event
        self._stop_event = threading.Event()

        self.raw_settings = raw_settings
        self.pipeline_settings = pipeline_settings

        self._monitor_cmd_q = collections.deque()
        self._monitor_ret_q = collections.deque()
        self._monitor_abort = threading.Event()
        self._monitor_thread = monitor_thread(self._monitor_cmd_q,
            self._monitor_ret_q, self._monitor_abort, self.pipeline_settings)
        self._monitor_thread.start()

        self._commands = {'set_data_dir'    : self._set_data_dir,
            'set_fprefix': self._set_fprefix,
            'set_data_dir_and_fprefix': self._set_data_dir_fprefix,
            'set_ret_size'  : self._set_ret_size,
            'set_raw_settings'  : self._set_raw_settings,
            'set_experiment'    : self._set_experiment,
            }

        self.waiting_for_header = []
        self.failed = []
        self.header_timeout = 10
        self.ret_every = 10

        self._exp_id = None
        self.data_dir = None

    def run(self):
        while True:
            if self._stop_event.is_set():
                break

            try:
                cmd, args = self._cmd_q.popleft()
            except IndexError:
                cmd = None

            if cmd is not None:
                logger.info("Processing cmd '%s' with args: %s", cmd,
                    ', '.join(['{}'.format(a) for a in args]))
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
                else:
                    imgs, fnames, img_hdrs, counters  = self._load_image_and_counter(img)

                    got_new_images = True

                    if len(imgs) > 0:
                        new_succeded.append(item)

                        self._ret_q.append([imgs, img_hdrs, counters, fnames,
                            exp_id, data_dir])



            if self._abort_event.is_set():
                self._abort()
                new_failed = []
                new_succeded = []

            for item in new_failed:
                self.failed.append(item)
                self.waiting_for_header.remove(item)

            for item in new_succeded:
                self.waiting_for_header.remove(item)

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

                exp_id = img_data[0]
                data_dir = img_data[1]
                imgs = img_data[2]

                if exp_id is None:
                    exp_id = self._exp_id

                if data_dir is None:
                    data_dir = self.data_dir

                for img in imgs:
                    if self._abort_event.is_set():
                        break

                    imgs, fnames, img_hdrs, counters  = self._load_image_and_counter(img)

                    if len(imgs) > 0:
                        new_imgs.extend(imgs)
                        new_img_hdrs.extend(img_hdrs)
                        new_counters.extend(counters)
                        new_fnames.extend(fnames)

                        got_new_images = True

                    else:
                        self.waiting_for_header.append((img, time.time(), exp_id, data_dir))

                    if len(imgs) >= self.ret_every:
                        self._ret_q.append([new_imgs, new_img_hdrs,
                            new_counters, new_fnames, exp_id, data_dir])
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

        logger.debug("Quitting pipeline monitor_and_load thread")

    def _set_data_dir(self, data_dir):
        logger.debug('Setting data directory: %s', data_dir)

        self._monitor_cmd_q.append(['set_data_dir', [data_dir,]])

        self.data_dir = data_dir

    def _set_fprefix(self, fprefix):
        self._monitor_cmd_q.append(['set_fprefix', [fprefix,]])

    def _set_data_dir_fprefix(self, data_dir, fprefix):
        logger.debug('Setting data directory: %s', data_dir)

        self._monitor_cmd_q.append(['set_data_dir_and_fprefix', [data_dir, fprefix]])

        self.data_dir = data_dir

    def _set_experiment(self, data_dir, fprefix, exp_id):
        logger.debug('Setting experiment to %s, data directory to %s', exp_id, data_dir)
        self._exp_id = exp_id
        self.data_dir = data_dir

        self._monitor_cmd_q.append(['set_experiment', [data_dir, fprefix, exp_id]])

    def _set_ret_size(self, size):
        logger.debug('Setting return chunk size to %i', size)
        self.ret_every = int(size)

    def _load_image_and_counter(self, img_name):
        logger.debug('Loading image %s', img_name)

        file_size = -1
        while file_size != os.path.getsize(img_name):
            file_size = os.path.getsize(img_name)

        try:
            imgs, img_hdrs, counters, fnames, failed = load_images_and_counters([img_name,],
                self.raw_settings)
        except SASExceptions.HeaderLoadError:
            imgs = img_hdrs = counters = fnames = failed = []

        if len(failed) > 0:
            logger.error('Failed to load image %s', img_name)

        return imgs, fnames, img_hdrs, counters

    def _set_raw_settings(self, settings):
        logger.debug('Setting RAW settings')

        self.raw_settings = settings

    def _abort(self):
        self._monitor_abort.set()

        self._cmd_q.clear()
        self._ret_q.clear()

        self._monitor_cmd_q.clear()
        self._monitor_ret_q.clear()

        self.waiting_for_header = []
        self.failed = []

        self._abort_event.clear()


    def stop(self):
        """Stops the process cleanly."""
        self._monitor_thread.stop()
        self._monitor_thread.join()

        self._stop_event.set()

class monitor_thread(threading.Thread):

    def __init__(self, cmd_q, ret_q, abort_event, pipeline_settings):
        threading.Thread.__init__(self)
        self.daemon = True
        self.name = 'monitor_thread'

        logger.debug("Initializing pipeline monitor thread")

        self._cmd_q = cmd_q
        self._ret_q = ret_q
        self._abort_event = abort_event
        self._stop_event = threading.Event()

        self.pipeline_settings = pipeline_settings
        self.fprefix = ''

        self._exp_id = None
        self.data_dir = None
        self.dir_snapshot = dirsnapshot.EmptyDirectorySnapshot()

        self._commands = {'set_data_dir': self._set_data_dir,
            'set_fprefix': self._set_fprefix,
            'set_data_dir_and_fprefix': self._set_data_dir_fprefix,
            'set_experiment': self._set_experiment,
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
                logger.info("Processing cmd '%s' with args: %s", cmd,
                    ', '.join(['{}'.format(a) for a in args]))

                self._commands[cmd](*args)

            if self.data_dir is not None:
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

        logger.debug("Quitting pipeline monitor thread")

    def _set_data_dir(self, data_dir):
        logger.debug('Setting data directory: %s', data_dir)

        self.data_dir = data_dir

        self.dir_snapshot = dirsnapshot.EmptyDirectorySnapshot()

    def _set_fprefix(self, fprefix):
        logger.debug('Setting file prefix: %s', fprefix)

        self.fprefix = fprefix

    def _set_data_dir_fprefix(self, data_dir, fprefix):
        self._set_fprefix(fprefix)
        self._set_data_dir(data_dir)

    def _set_experiment(self, data_dir, fprefix, exp_id):
        logger.debug('Setting experiment: %s', exp_id)

        self._exp_id = exp_id
        self._set_data_dir_fprefix(data_dir, fprefix)

    def _check_for_new_files(self):

        new_snapshot = dirsnapshot.DirectorySnapshot(self.data_dir, False)

        diff = dirsnapshot.DirectorySnapshotDiff(self.dir_snapshot, new_snapshot)

        self.dir_snapshot = new_snapshot

        new_images = []

        for f in diff.files_created:
            if self._abort_event.is_set():
                new_images = []
                self._abort()

            if (os.path.splitext(f)[1] in self.pipeline_settings['image_exts']
                and os.path.split(f)[1].startswith(self.fprefix)):
                new_images.append(os.path.abspath(os.path.expanduser(f)))

        if new_images:
            logger.debug('New images found: %s', ', '.join(new_images))

        return new_images

    def _abort(self):
        self._cmd_q.clear()
        self._ret_q.clear()

        self.dir_snapshot = dirsnapshot.DirectorySnapshot(self.data_dir, False)

        self._abort_event.clear()

    def stop(self):
        """Stops the thread cleanly."""
        self._stop_event.set()


class raver_process(multiprocessing.Process):

    def __init__(self, cmd_q, ret_q, cmd_lock, ret_lock, abort_event,
        raw_settings_file, log_lock, log_queue):
        multiprocessing.Process.__init__(self)
        self.daemon = True

        self._cmd_q = cmd_q
        self._ret_q = ret_q
        self._cmd_lock = cmd_lock
        self._ret_lock = ret_lock
        self._abort_event = abort_event
        self._stop_event = multiprocessing.Event()

        self._log_lock = log_lock
        self._log_queue = log_queue

        self._log('debug', "Initializing pipeline radial averaging process")

        self.raw_settings = raw.load_settings(raw_settings_file)

        self._commands = {'raver_images': self._raver_images,
            }

        self.ret_every = 10

    def run(self):
        while True:
            if self._stop_event.is_set():
                break

            if self._abort_event.is_set():
                self._abort()

            try:
                with self._cmd_lock:
                    cmd, exp_id, args, kwargs = self._cmd_q.get_nowait()
            except queue.Empty:
                cmd = None

            if cmd is not None:
                self._log('debug', "Processing cmd '%s'" %(cmd))

                self._commands[cmd](exp_id, *args, **kwargs)

            else:
                time.sleep(0.1)

        self._log('debug', "Quiting reduction process")

    def _raver_images(self, exp_id, *args, **kwargs):
        self._log('debug', 'Radially averaging images')
        img_hdrs = []
        all_counters = []
        load_path = ''

        imgs = args[0]
        if len(args) >1:
            names = args[1]
        if len(args) > 2:
            img_hdrs = args[2]
        if len(args) > 3:
            all_counters = args[3]


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
                with self._ret_lock:
                    self._ret_q.put_nowait([exp_id, profiles])
                profiles = []

        if len(profiles) > 0:
            with self._ret_lock:
                self._ret_q.put_nowait([exp_id, profiles])

    def load_settings(self, settings_file):
        self._log('debug', 'Loading RAW settings from %s' %(settings_file))

        self.raw_settings = raw.load_settings(settings_file)

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

"""
To do:

3) Consider using more advanced functions from the watchdog package:
https://blog.magrathealabs.com/filesystem-events-monitoring-with-python-9f5329b651c3
http://thepythoncorner.com/dev/how-to-create-a-watchdog-in-python-to-look-for-filesystem-changes/
https://www.geeksforgeeks.org/create-a-watchdog-in-python-to-look-for-filesystem-changes/

"""
