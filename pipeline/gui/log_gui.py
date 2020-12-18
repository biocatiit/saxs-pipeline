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

import logging
import traceback

if __name__ != '__main__':
    logger = logging.getLogger(__name__)

import wx

class LogPanel(wx.Panel):

    def __init__(self, parent, settings, *args, **kwargs):
        wx.Panel.__init__(self, parent, *args, **kwargs)

        self.settings = settings
        self.main_frame = self.FindWindowByName('PipelineFrame')

        self._create_layout()

        self.txt_handler = CustomConsoleHandler(self.log, 100000)

        self.txt_handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(threadName)s - %(levelname)s - %(message)s')
        # formatter = logging.Formatter('%(asctime)s - %(threadName)s - %(levelname)s - %(message)s')
        self.txt_handler.setFormatter(formatter)

        my_logger = logging.getLogger()
        my_logger.addHandler(self.txt_handler)

    def _FromDIP(self, size):
        # This is a hack to provide easy back compatibility with wxpython < 4.1
        try:
            return self.FromDIP(size)
        except Exception:
            return size

    def _create_layout(self):

        ctrl_sizer = wx.StaticBoxSizer(wx.VERTICAL, self, 'Log Controls')
        ctrl_box = ctrl_sizer.GetStaticBox()

        log_level_label = wx.StaticText(ctrl_box, label='Log level:')
        log_level = wx.Choice(ctrl_box, choices=['Debug', 'Info', 'Warning',
            'Error', 'Critical'])
        log_level.SetStringSelection('Debug')
        log_level.Bind(wx.EVT_CHOICE, self._on_log_level)

        log_size_label = wx.StaticText(ctrl_box, label='Log length:')
        self.log_size = wx.TextCtrl(ctrl_box, value='100000')
        log_size_set = wx.Button(ctrl_box, label='Set')
        log_size_set.Bind(wx.EVT_BUTTON, self._on_log_size)

        ctrl_grid_sizer = wx.FlexGridSizer(cols=3, vgap=self._FromDIP(5),
            hgap=self._FromDIP(5))
        ctrl_grid_sizer.Add(log_level_label, flag=wx.ALIGN_CENTER_VERTICAL)
        ctrl_grid_sizer.Add(log_level, flag=wx.ALIGN_CENTER_VERTICAL)
        ctrl_grid_sizer.Add((1,1))
        ctrl_grid_sizer.Add(log_size_label, flag=wx.ALIGN_CENTER_VERTICAL)
        ctrl_grid_sizer.Add(self.log_size, flag=wx.ALIGN_CENTER_VERTICAL)
        ctrl_grid_sizer.Add(log_size_set, flag=wx.ALIGN_CENTER_VERTICAL)

        ctrl_sizer.Add(ctrl_grid_sizer, flag=wx.ALL|wx.EXPAND,
            border=self._FromDIP(5))


        log_sizer = wx.StaticBoxSizer(wx.VERTICAL, self, 'Pipeline Log')
        log_box = log_sizer.GetStaticBox()

        self.log = wx.TextCtrl(log_box, style=wx.TE_MULTILINE|wx.TE_READONLY)

        log_sizer.Add(self.log, proportion=1, flag=wx.ALL|wx.EXPAND,
            border=self._FromDIP(5))

        top_sizer = wx.BoxSizer(wx.VERTICAL)
        top_sizer.Add(ctrl_sizer, flag=wx.EXPAND|wx.ALL, border=self._FromDIP(5))
        top_sizer.Add(log_sizer, proportion=1,
            flag=wx.LEFT|wx.RIGHT|wx.BOTTOM|wx.EXPAND, border=self._FromDIP(5))

        self.SetSizer(top_sizer)

    def _on_log_level(self, evt):
        log_level = eval('logging.{}'.format(evt.GetString().upper()))
        self.txt_handler.setLevel(log_level)

    def _on_log_size(self, evt):
        try:
            length = int(self.log_size.GetValue())
        except Exception:
            length = None

        if length is not None:
            self.txt_handler.set_max_length(length)

    def on_exit(self):
        pass


class CustomConsoleHandler(logging.StreamHandler):
    """
    Logs python logging statements to a wxpython textctrl. Based on code from:
    https://www.blog.pythonlibrary.org/2013/08/09/wxpython-how-to-redirect-pythons-logging-module-to-a-textctrl/
    """

    def __init__(self, textctrl, length=-1):
        """"""
        logging.StreamHandler.__init__(self)
        self.textctrl = textctrl
        self.length = length

    def emit(self, record):
        """Constructor"""
        msg = self.format(record)
        txt_length = self.textctrl.GetLastPosition()

        if self.length != -1 and  txt_length > self.length:
            to_remove = len(msg) + (txt_length-self.length)
            final_length =  min(to_remove, txt_length)
            wx.CallAfter(self.textctrl.Remove, 0, final_length)
            wx.CallAfter(self.textctrl.SetInsertionPointEnd)

        wx.CallAfter(self.textctrl.WriteText, msg + "\n")
        self.flush()

    def set_max_length(self, length):
        self.length = length
        txt_length = self.textctrl.GetLastPosition()

        if self.length != -1 and  txt_length > self.length:
            to_remove = txt_length - self.length
            final_length =  min(to_remove, txt_length)
            wx.CallAfter(self.textctrl.Remove, 0, final_length)
            wx.CallAfter(self.textctrl.SetInsertionPointEnd)
