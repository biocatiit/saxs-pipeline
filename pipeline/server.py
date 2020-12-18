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

from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import object, range, map
from io import open

import threading
import logging
import logging.handlers as handlers
import collections
import traceback
import time
import sys
import os
import multiprocessing as mp

if __name__ != '__main__':
    logger = logging.getLogger(__name__)

import zmq
import wx

from . import pipeline
from . import settings


class ControlServer(threading.Thread):
    """

    """

    def __init__(self, ip, port, pipeline_cmd_q, pipeline_ret_q, name='ControlServer'):
        """
        Initializes the custom thread. Important parameters here are the
        list of known commands ``_commands`` and known pumps ``known_pumps``.

        :param collections.deque command_queue: The queue used to pass commands to
            the thread.

        :param threading.Event stop_event: An event that is set when the thread
            needs to abort, and otherwise is not set.
        """
        threading.Thread.__init__(self, name=name)
        self.daemon = True

        logger.info("Initializing control server: %s", self.name)

        self.ip = ip
        self.port = port
        self.pipeline_cmd_q = pipeline_cmd_q
        self.pipeline_ret_q = pipeline_ret_q

        self._stop_event = threading.Event()

        self.context = zmq.Context()
        self.socket = self.context.socket(zmq.PAIR)
        self.socket.bind("tcp://{}:{}".format(self.ip, self.port))


    def run(self):
        """
        Custom run method for the thread.
        """
        while True:
            try:
                if self.socket.poll(10) > 0:
                    logger.debug("Getting new command")
                    command = self.socket.recv_json()
                else:
                    command = None
            except Exception:
                command = None

            if self._stop_event.is_set():
                logger.debug("Stop event detected")
                break

            if command is not None:
                cmd = command['command']
                get_response = command['response']
                logger.info("Processing cmd '%s' with args: %s and kwargs: %s ",
                    cmd[0], ', '.join(['{}'.format(a) for a in cmd[1]]),
                    ', '.join(['{}: {}'.format(kw, item) for kw, item in cmd[2].items()]))

                try:
                    if cmd[0] == 'ping':
                        answer = 'ping received'

                    else:
                        self.pipeline_cmd_q.append(cmd)

                        if get_response:
                            start_time = time.time()
                            while len(self.pipeline_ret_q) == 0 and time.time()-start_time < 5:
                                time.sleep(0.01)

                            if len(self.pipeline_ret_q) == 0:
                                answer = ''
                            else:
                                answer = self.pipeline_ret_q.popleft()
                        else:
                            answer = 'cmd sent'

                    if answer == '':
                        logger.exception('No response received from device')
                    else:
                        logger.debug('Sending command response: %s', answer)
                        self.socket.send_json(answer)

                except Exception:
                    cmd = command['command']
                    msg = ("Failed to run command '%s' with args: %s and "
                        "kwargs: %s. Exception follows:" %(cmd[0],
                        ', '.join(['{}'.format(a) for a in cmd[1]]),
                        ', '.join(['{}:{}'.format(kw, item) for kw, item in cmd[2].items()])))
                    logger.exception(msg)
                    logger.exception(traceback.print_exc())

            else:
                time.sleep(0.01)

        if self._stop_event.is_set():
            self._stop_event.clear()

        logger.info("Quitting control server: %s", self.name)

    def stop(self):
        """Stops the thread cleanly."""
        # logger.info("Starting to clean up and shut down pump control thread: %s", self.name)
        self.socket.unbind("tcp://{}:{}".format(self.ip, self.port))
        self.socket.close()
        self.context.destroy()

        self._stop_event.set()

if __name__ == '__main__':
    #Example usage
    mp.set_start_method('spawn')

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    h1 = logging.StreamHandler(sys.stdout)
    h1.setLevel(logging.DEBUG)
    # h1.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(threadName)s - %(levelname)s - %(message)s')
    h1.setFormatter(formatter)
    logger.addHandler(h1)

    cmd_q = collections.deque()
    ret_q = collections.deque()
    abort_event = threading.Event()

    pl_settings = settings.Settings()

    raw_settings_file = 'SAXS.cfg'

    pl_thread = pipeline.pipeline_thread(cmd_q, ret_q, abort_event,
        pl_settings, raw_settings_file)

    pl_thread.start()

    port = '5556'
    ip = '192.168.1.23'

    pipeline_server = ControlServer(ip, port, cmd_q, ret_q, 'SAXSPipeServer')
    pipeline_server.start()

    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        pipeline_server.stop()
        pipeline_server.join()

    logger.info("Quitting server")
