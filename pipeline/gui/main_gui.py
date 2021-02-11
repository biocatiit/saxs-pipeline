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
import logging.handlers as handlers
import signal
import sys
import os
import multiprocessing as mp
import traceback
import collections
import threading
import platform

import wx

pipeline_path = os.path.abspath(os.path.join('.', __file__, '..', '..', '..'))
if pipeline_path not in os.sys.path:
    os.sys.path.append(pipeline_path)

from pipeline.gui import gui_utils
from pipeline import settings
from pipeline.gui import settings_gui
from pipeline.gui import log_gui
from pipeline.gui import status_gui
import pipeline

if __name__ != '__main__':
    logger = logging.getLogger(__name__)

class PipelineFrame(wx.Frame):
    """
    """
    def __init__(self, server_ip, server_port, *args, **kwargs):
        """
        """
        super(PipelineFrame, self).__init__(None, *args, **kwargs)
        logger.debug('Setting up the BioFrame')

        self.component_sizers = {}
        self.component_panels = {}

        self.Bind(wx.EVT_CLOSE, self._on_exit)

        self.settings = settings.Settings()

        self.standard_paths = wx.StandardPaths.Get()

        self.local_data_dir = self.standard_paths.GetUserLocalDataDir().replace('main_gui', 'pipeline')

        if not os.path.exists(self.local_data_dir):
            os.mkdir(self.local_data_dir)

        self._component_panels = []

        self.pipeline_cmd_q = None
        self.pipeline_ret_q = None
        self.pipeline_abort_event = None
        self.pipeline_thread = None
        self.pipeline_server = None

        self.server_ip = server_ip
        self.server_port = server_port

        self._create_layout()

        self.sleep_inhibit = gui_utils.SleepInhibit()
        self.sleep_inhibit.on()
        self._bind_signals()



        client_display = wx.GetClientDisplayRect()
        minsize = (min(400, client_display.Width), min(600, client_display.Height))
        self.SetMinSize(self._FromDIP(minsize))

        self.Fit()

        size = (min(450, client_display.Width), min(600, client_display.Height))
        self.SetSize(self._FromDIP(size))

        self.Raise()
        self.Show()

        self._on_startup()

        if os.path.exists(self.settings['raw_settings_file']):
            wx.CallAfter(self._start_pipeline)


    def _FromDIP(self, size):
        # This is a hack to provide easy back compatibility with wxpython < 4.1
        try:
            return self.FromDIP(size)
        except Exception:
            return size

    def _create_layout(self):
        top_panel = wx.Panel(self)

        notebook = wx.Notebook(top_panel)

        self.log_panel = log_gui.LogPanel(notebook, self.settings, name='LogPanel')
        self.settings_panel = settings_gui.SettingsPanel(notebook, self.settings,
            name='SettingsPanel')
        self.status_panel = status_gui.StatusPanel(notebook, self.settings,
            name='StatusPanel')

        self._component_panels.append(self.settings_panel)
        self._component_panels.append(self.log_panel)
        self._component_panels.append(self.status_panel)

        notebook.AddPage(self.status_panel, 'Status')
        notebook.AddPage(self.settings_panel, 'Settings')
        notebook.AddPage(self.log_panel, 'Log')

        top_panel_sizer = wx.BoxSizer()
        top_panel_sizer.Add(notebook, proportion=1, flag=wx.EXPAND)
        top_panel.SetSizer(top_panel_sizer)

        top_sizer = wx.BoxSizer()
        top_sizer.Add(top_panel, proportion=1, flag=wx.EXPAND)

        self.SetSizer(top_sizer)

    def _on_startup(self):

        settings_file = os.path.join(self.local_data_dir, 'backup.pcfg')

        if os.path.exists(settings_file):
            dlg = wx.MessageDialog(parent=self,
                message='Load last used configuration?',
                caption='Restore configuration',
                style=wx.YES_NO|wx.ICON_QUESTION)
            answer = dlg.ShowModal()
            dlg.Destroy()

            if answer == wx.ID_YES:
                settings.load_settings(self.settings, settings_file)
                self.settings_panel.set_values()

    def _bind_signals(self):
        """Binds sigint and sigterm so the pipeline can gracefully exit"""

        self.orig_sigint_handler = signal.getsignal(signal.SIGINT)
        self.orig_sigterm_handler = signal.getsignal(signal.SIGTERM)
        signal.signal(signal.SIGINT, gui_utils.signal_handler)
        signal.signal(signal.SIGTERM, gui_utils.signal_handler)

        #This is a bit silly, but seems necessary on MacOS if the pipeline is a background process
        self.heartbeat = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self._on_heartbeat, self.heartbeat)
        self.heartbeat.Start(5000)

    def _on_heartbeat(self, evt):
        pass

    def create_file_dialog(self, mode, name='Config files', ext='*.pcfg', dir_msg=''):
        f = None

        if mode == wx.FD_OPEN:
            filters = name + ' (' + ext + ')|' + ext + '|All files (*.*)|*.*'
            dialog = wx.FileDialog( None, style=mode, wildcard=filters)

        elif mode == wx.FD_SAVE:
            filters = name + ' ('+ext+')|'+ext
            dialog = wx.FileDialog( None, style=mode|wx.FD_OVERWRITE_PROMPT,
                wildcard=filters)
        else:
            dialog = wx.DirDialog(self, message=dir_msg)

        # Show the dialog and get user input
        if dialog.ShowModal() == wx.ID_OK:
            f = dialog.GetPath()

        # Destroy the dialog
        dialog.Destroy()

        return f

    def _start_pipeline(self):
        logger.debug('Starting the pipeline')
        self.pipeline_cmd_q = collections.deque()
        self.pipeline_ret_q = collections.deque()
        self.pipeline_abort_event = threading.Event()
        self.pipeline_ret_lock = threading.RLock()

        self.pipeline_thread = pipeline.control.pipeline_thread(self.pipeline_cmd_q,
            self.pipeline_ret_q, self.pipeline_abort_event, self.settings,
            self.pipeline_ret_lock)

        self.pipeline_thread.start()

        self.status_panel.start_btn.Disable()
        self.status_panel.abort_btn.Enable()
        self.status_panel.pipeline_status.SetLabel('Running')

        self.pipeline_server = pipeline.server.ControlServer(self.server_ip,
            self.server_port, self.pipeline_cmd_q, self.pipeline_ret_q,
            'SAXSPipeServer')

        self.pipeline_server.start()

    def _on_exit(self, evt):
        self._cleanup_and_quit()

    def cleanup_and_quit_forced(self):
        signal.signal(signal.SIGINT, self.orig_sigint_handler)
        signal.signal(signal.SIGTERM, self.orig_sigterm_handler)

        self._cleanup_and_quit()

    def _cleanup_and_quit(self):
        logger.debug('Closing the PipelineFrame')

        logger.removeHandler(self.log_panel.txt_handler)

        if self.pipeline_server is not None:
            self.pipeline_server.stop()
            self.pipeline_server.join()

        if self.pipeline_thread is not None:
            self.pipeline_thread.stop()
            self.pipeline_thread.join()

        sys.excepthook = sys.__excepthook__

        try:
            dialog_list = wx.GetApp().modal_dialog_hook.all_modal_dialogs

            for dialog in dialog_list:
                dialog.EndModal(wx.ID_CANCEL)
                dialog.Destroy()

            self.heartbeat.Stop()

            self._closing = True
            self.sleep_inhibit.force_off()

        except Exception:
            pass

        finally:
            for panel in self.component_panels:
                panel.on_exit()

            for w in wx.GetTopLevelWindows():
                if w != self:
                    w.Destroy()

            self.Destroy()


class MyApp(wx.App):
    """The top level wx.App that we subclass to add an exceptionhook."""

    def OnInit(self):
        """Initializes the app. Calls the :class:`MainFrame`"""

        standard_paths = wx.StandardPaths.Get() #Can't do this until you start the wx app
        info_dir = standard_paths.GetUserLocalDataDir().replace('main_gui', 'pipeline')

        if not os.path.exists(info_dir):
            os.mkdir(info_dir)

        h2 = handlers.RotatingFileHandler(os.path.join(info_dir, 'pipeline.log'),
            maxBytes=10e6, backupCount=100, delay=True)
        # h2.setLevel(logging.INFO)
        h2.setLevel(logging.DEBUG)
        formatter2 = logging.Formatter('%(asctime)s - %(threadName)s - %(levelname)s - %(message)s')
        h2.setFormatter(formatter2)

        logger.addHandler(h2)

        # sys.excepthook = self.ExceptionHook

        self.modal_dialog_hook = gui_utils.MyModalDialogHook()
        self.modal_dialog_hook.Register()

        title = 'SAXS Pipeline'

        server_ip = '164.54.204.82'
        # server_ip = '192.168.1.14'
        server_port = '5556'

        if len(sys.argv) == 2:
            server_ip = sys.argv[1]
        elif len(sys.argv) == 3:
            server_ip = sys.argv[1]
            server_port = sys.argv[2]

        frame = PipelineFrame(server_ip, server_port, title=title,
            name='PipelineFrame')
        self.SetTopWindow(frame)

        return True

    def BringWindowToFront(self):
        """
        Overwrites this default method to deal with the possibility that it
        is called when the frame is closed.
        """
        try: # it's possible for this event to come when the frame is closed
            self.GetTopWindow().Raise()
        except:
            pass

    def ExceptionHook(self, errType, value, trace):
        """
        Creates an exception hook that catches all uncaught exceptions and informs
        users of them, instead of letting the program crash. From
        http://apprize.info/python/wxpython/10.html
        """
        err = traceback.format_exception(errType, value, trace)
        errTxt = "\n".join(err)
        msg = ("An unexpected error has occurred, please report it to the "
                "developers. You may need to restart RAW to continue working"
                "\n\nError:\n%s" %(errTxt))

        logger.error('Error in GUI:\n{}'.format(traceback.format_exc()))

        if self and self.IsMainLoopRunning():
            if not self.HandleError(value):
                wx.CallAfter(wx.lib.dialogs.scrolledMessageDialog, None, msg,
                    "Unexpected Error")
        else:
            sys.stderr.write(msg)

    def HandleError(self, error):
        """
        Override in subclass to handle errors

        :return: True to allow program to continue running without showing error.
            False to show the error.
        """
        return False

def main():
    gui_utils.setup_thread_excepthook()

    logger.debug('Setting up wx app')

    app = MyApp(0)
    app.MainLoop()

if __name__ == '__main__':
    #if platform.system() == 'Darwin':
    mp.set_start_method('spawn')

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    h1 = logging.StreamHandler(sys.stdout)
    # h1.setLevel(logging.DEBUG)
    h1.setLevel(logging.INFO)
    # formatter = logging.Formatter('%(asctime)s - %(name)s - %(threadName)s - %(levelname)s - %(message)s')
    formatter = logging.Formatter('%(asctime)s - %(threadName)s - %(levelname)s - %(message)s')
    h1.setFormatter(formatter)
    logger.addHandler(h1)

    try:
        for module_logger in logging.root.manager.loggerDict.keys():
            if not module_logger.startswith('pipeline'):
                logging.getLogger(module_logger).setLevel(60)

    except Exception:
        traceback.print_exc()

    main()
