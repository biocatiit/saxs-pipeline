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
import collections
import traceback
import time
import sys

if __name__ != '__main__':
    logger = logging.getLogger(__name__)

import zmq


class ControlClient(threading.Thread):
    """

    """

    def __init__(self, ip, port, command_queue, answer_queue, abort_event,
        timeout_event, name='ControlClient'):
        """
        """
        threading.Thread.__init__(self, name=name)
        self.daemon = True

        logger.info("Starting control client: %s", self.name)

        self.ip = ip
        self.port = port
        self.command_queue = command_queue
        self.answer_queue = answer_queue
        self._abort_event = abort_event
        self._stop_event = threading.Event()
        self.timeout_event = timeout_event

        self.connect_error = 0

        self.heartbeat = 60
        self.last_ping = 0

        self.resend_missed_commands_on_reconnect = True
        self.missed_cmds = collections.deque()

    def run(self):
        """
        Custom run method for the thread.
        """

        logger.info("Connecting to %s on port %s", self.ip, self.port)
        self.context = zmq.Context()
        self.socket = self.context.socket(zmq.PAIR)
        self.socket.set(zmq.LINGER, 0)
        self.socket.connect("tcp://{}:{}".format(self.ip, self.port))

        while True:
            try:
                if not self.socket.closed:
                    if time.time() - self.last_ping > self.heartbeat:
                        self.last_ping = time.time()
                        self._ping()
                else:
                    if self.socket.closed:
                        self._ping()

                if len(self.command_queue) > 0:
                    # logger.debug("Getting new command")
                    command = self.command_queue.popleft()
                else:
                    command = None

                if self._abort_event.is_set():
                    logger.debug("Abort event detected")
                    self._abort()
                    command = None

                if self._stop_event.is_set():
                    logger.debug("Stop event detected")
                    break

                if command is not None:
                    # logger.debug("For device %s, processing cmd '%s' with args: %s and kwargs: %s ", device, cmd[0], ', '.join(['{}'.format(a) for a in cmd[1]]), ', '.join(['{}:{}'.format(kw, item) for kw, item in cmd[2].items()]))

                    if not self.socket.closed:
                        self._send_cmd(command)

                    elif self.resend_missed_commands_on_reconnect:
                        self.missed_cmds.append(command)

                else:
                    time.sleep(0.01)

            except Exception:
                logger.error('Error in client thread:\n{}'.format(traceback.format_exc()))

        if self._stop_event.is_set():
            self._stop_event.clear()
        else:
            self._abort()

        if not self.socket.closed:
            self.socket.disconnect("tcp://{}:{}".format(self.ip, self.port))
            self.socket.close(0)

        self.context.destroy(0)
        logger.info("Quitting remote client thread: %s", self.name)

    def _send_cmd(self, command):
        cmd = command['command']
        get_response = command['response']

        try:
            self.socket.send_json(command)

            start_time = time.time()
            while self.socket.poll(10) == 0 and time.time()-start_time < 5:
                pass

            if self.socket.poll(10) > 0:
                answer = self.socket.recv_json()
            else:
                answer = ''

            if answer == '':
                raise zmq.ZMQError(msg="Could not get a response from the server")
            else:
                self.connect_error = 0

            # logger.debug('Command response: %s' %(answer))

            if get_response:
                self.answer_queue.append(answer)

        except zmq.ZMQError:
            cmd = command['command']
            msg = ("Pipeline failed to run command '%s' with args: %s and "
                "kwargs: %s. Timeout or other ZMQ error." %(cmd[0],
                ', '.join(['{}'.format(a) for a in cmd[1]]),
                ', '.join(['{}:{}'.format(kw, item) for kw, item in cmd[2].items()])))
            logger.error(msg)
            self.connect_error += 1
            self._ping()
            if not self.timeout_event.set():
                self.answer_queue.append(None)

            self.missed_cmds.append(command)

        except Exception:
            cmd = command['command']
            msg = ("Pipeline failed to run command '%s' with args: %s "
                "and kwargs: %s. Exception follows:" %(cmd[0],
                ', '.join(['{}'.format(a) for a in cmd[1]]),
                ', '.join(['{}:{}'.format(kw, item) for kw, item in cmd[2].items()])))
            logger.error(msg)
            logger.error(traceback.print_exc())
            self.connect_error += 1

            self.missed_cmds.append(command)

        if self.connect_error > 5 and not self.timeout_event.is_set():
            msg = ('5 consecutive failures to run a pipeline command.')
            logger.error(msg)
            logger.error("Connection timed out")
            self.timeout_event.set()

    def _ping(self):
        # logger.debug("Checking if server is active")
        cmd = {'device': 'server', 'command': ('ping', (), {}), 'response': False}

        connect_tries = 0

        if not self.socket.closed:
            while connect_tries < 5:
                self.socket.send_json(cmd)

                start_time = time.time()
                while self.socket.poll(10) == 0 and time.time()-start_time < 1:
                    pass

                if self.socket.poll(10) > 0:
                    answer = self.socket.recv_json()
                else:
                    answer = ''

                if answer == 'ping received':
                    logger.debug("Connection to server verified")
                    connect_tries = 5
                else:
                    logger.error("Could not get a response from the server")
                    connect_tries = connect_tries+1

                    if connect_tries == 5:
                        logger.error("Connection timed out")
                        self.timeout_event.set()
                        self.connect_error = 6
                        self.socket.disconnect("tcp://{}:{}".format(self.ip, self.port))
                        self.socket.close(0)

        else:
            self.socket = self.context.socket(zmq.PAIR)
            self.socket.set(zmq.LINGER, 0)
            self.socket.connect("tcp://{}:{}".format(self.ip, self.port))

            while connect_tries < 5:
                self.socket.send_json(cmd)

                start_time = time.time()
                while self.socket.poll(10) == 0 and time.time()-start_time < 0.1:
                    pass

                if self.socket.poll(10) > 0:
                    answer = self.socket.recv_json()
                else:
                    answer = ''

                if answer == 'ping received':
                    logger.debug("Connection to server verified")
                    connect_tries = 5
                    self.timeout_event.clear()
                    self.connect_error = 0

                    if self.resend_missed_commands_on_reconnect:
                        while len(self.missed_cmds) > 0 and not self.timeout_event.is_set():
                            cmd = self.missed_cmds.popleft()
                            self._send_cmd(cmd)

                else:
                    connect_tries = connect_tries+1

            if self.timeout_event.is_set():
                self.socket.disconnect("tcp://{}:{}".format(self.ip, self.port))
                self.socket.close(0)



    def _abort(self):
        logger.info("Aborting remote client thread %s current and future commands", self.name)
        self.command_queue.clear()
        self._abort_event.clear()
        logger.debug("Remote client thread %s aborted", self.name)

    def stop(self):
        """Stops the thread cleanly."""
        logger.info("Starting to clean up and shut down remote client thread: %s", self.name)

        self._stop_event.set()

if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    h1 = logging.StreamHandler(sys.stdout)
    # h1.setLevel(logging.DEBUG)
    h1.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(threadName)s - %(levelname)s - %(message)s')
    h1.setFormatter(formatter)
    logger.addHandler(h1)

    port = '5556'
    ip = '164.54.204.37'

    cmd_q = collections.deque()
    return_q = collections.deque()
    abort_event = threading.Event()
    timeout_event = threading.Event()

    control_client = ControlClient(ip, port, cmd_q, return_q, abort_event,
        timeout_event, name='PipelineCtrlClient')
    control_client.start()

    data_dir = '.'
    fprefix = 'test'
    num_exps = 1000
    exp_cmd = ('start_experiment', ['sec_saxs_2', 'SEC', data_dir, fprefix],
        {'num_exps': num_exps})

    exp_client_cmd = {'command': exp_cmd, 'response': False}
    cmd_q.append(exp_client_cmd)

    control_client.stop()
