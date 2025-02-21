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
import json
import requests

if __name__ != '__main__':
    logger = logging.getLogger(__name__)

import numpy as np
import zmq
import bitshuffle
import lz4


class EigerStreamClient(threading.Thread):
    """

    """

    def __init__(self, ip, port, data_queue, name='EigerStreamClient'):
        """
        """
        threading.Thread.__init__(self, name=name)
        self.daemon = True

        logger.info("Starting Eiger stream client: %s", self.name)

        self.ip = ip
        self.port = port
        self.data_queue = data_queue
        self._stop_event = threading.Event()

        #Should I try to set the areadetector/eiger variables so that stream is setup properly?

    def run(self):
        """
        Custom run method for the thread.
        """

        logger.info("Connecting to %s on port %s", self.ip, self.port)
        self.context = zmq.Context()
        self.socket = self.context.socket(zmq.PULL)
        self.socket.set(zmq.LINGER, 0)
        self.socket.connect("tcp://{}:{}".format(self.ip, self.port))

        while True:
            try:
                try:
                    frames = self.socket.recv_multipart(zmq.NOBLOCK, copy=False)
                    self.data_queue.append(frames)
                except zmq.ZMQError:
                    pass

                if self._stop_event.is_set():
                    logger.debug("Stop event detected")
                    break

            except Exception:
                logger.error('Error in client thread:\n{}'.format(traceback.format_exc()))

        if self._stop_event.is_set():
            self._stop_event.clear()

        if not self.socket.closed:
            self.socket.disconnect("tcp://{}:{}".format(self.ip, self.port))
            self.socket.close(0)

        self.context.destroy(0)
        logger.info("Quitting remote client thread: %s", self.name)


    def stop(self):
        """Stops the thread cleanly."""
        logger.info("Starting to clean up and shut down remote client thread: %s", self.name)

        self._stop_event.set()

class EigerStreamParser():

    def decodeFrames(self, frames, return_header=True):
        """
        decode and proces EIGER ZMQ stream frames
        """
        try:
            header = json.loads(frames[0].bytes)
        except Exception:
            logger.error('Error parsing Eiger header:\n{}'.format(traceback.format_exc()))

        logger.debug(header)

        if header["htype"].startswith("dimage-"):
            data =self._decodeImage(frames)
        elif header["htype"].startswith("dheader-"):
            #I don't think we need to do this, and it will just waste time
            # self._decodeHeader(frames, header)
            data = None
        elif header["htype"].startswith("dseries_end"):
            #I don't think we need to do this, and it will just waste time
            # self._decodeEndOfSeries(frames)
            data = None
        else:
            logger.error("Unexpected message from Eiger stream")
            data = None

        if return_header:
            ret = [data, header]
        else:
            ret = [data, None]

        return ret

    def _decodeImage(self, frames):
        """
        decode ZMQ image frames
        """

        info = json.loads(frames[1].bytes) # info dict

        logger.debug(info)

        if info["encoding"] == "lz4<": #TODO: soft code flag
            data = self.readLZ4(frames[2], info["shape"], info["type"])
        elif info["encoding"] == "bs32-lz4<": #TODO: soft code flag
            data = self.readBSLZ4(frames[2], info["shape"], info["type"])
        elif info["encoding"] =="bs16-lz4<":
            data = self.readBSLZ4(frames[2], info["shape"], info["type"])
        else:
            raise IOError("[ERR] encoding %s is not implemented" %info["encoding"])

        return data

    def _decodeEndOfSeries(self, frames):
        logger.debug('End of Eiger series %s', json.loads(frames[0].bytes))
        return True

    def _decodeHeader(self, frames, header):
        """
        decode and process ZMQ header frames
        """
        if header["header_detail"]:
            logger.debug(header['header_detail'])

        if header["header_detail"] != "none":
            for key, value in json.loads(frames[1].bytes).iteritems():
                logger.debug(key, value)

        if header["header_detail"] == "all":
            if json.loads(frames[2].bytes)["htype"].startswith("dflatfield"):
                pass
            if json.loads(frames[4].bytes)["htype"].startswith("dpixelmask"):
                pass
            if json.loads(frames[6].bytes)["htype"].startswith("dcountrate"):
                pass

    def readBSLZ4(self, frame, shape, dtype):
        """
        unpack bitshuffle-lz4 compressed frame and return np array image data
        frame: zmq frame
        shape: image shape
        dtype: image data type
        """

        data = frame.bytes
        blob = np.fromstring(data[12:], dtype=np.uint8)
        # blocksize is big endian uint32 starting at byte 8, divided by element size
        dt = np.dtype(dtype)
        blocksize = np.ndarray(shape=(), dtype=">u4", buffer=data[8:12])/dt.itemsize
        imgData = bitshuffle.decompress_lz4(blob, shape[::-1], dt, blocksize)

        return imgData

    def readLZ4(self, frame, shape, dtype):
        """
        unpack lz4 compressed frame and return np array image data
        frame: zmq frame
        dataSize: data size of single value
        """

        dtype = np.dtype(dtype)
        dataSize = dtype.itemsize*shape[0]*shape[1] # bytes * image size

        imgData = lz4.block.decompress(struct.pack('<I', dataSize) + frame.bytes)

        return np.reshape(np.fromstring(imgData, dtype=dtype), shape[::-1])

if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    h1 = logging.StreamHandler(sys.stdout)
    h1.setLevel(logging.DEBUG)
    # h1.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(threadName)s - %(levelname)s - %(message)s')
    h1.setFormatter(formatter)
    logger.addHandler(h1)

    port = '9999'
    ip = '164.54.204.140'

    r = requests.put('http://164.54.204.140/stream/api/1.8.0/config/mode',
    data='{"value": "enabled"}')
    print(r)
    print(r.text)


    return_q = collections.deque()

    control_client = EigerStreamClient(ip, port, return_q)
    control_client.start()

    parser = EigerStreamParser()

    num_frames = 0

    try:
        while True:
            if len(return_q) > 0:
                frames = return_q.popleft()
                data = parser.decodeFrames(frames)
                if data is not None:
                    num_frames += 1
                    logger.info('Got frame %i', num_frames)
    except KeyboardInterrupt:
        pass

    control_client.stop()
    control_client.join()
