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
import os

if __name__ != '__main__':
    logger = logging.getLogger(__name__)

import wx

class StatusPanel(wx.Panel):

    def __init__(self, parent, settings, *args, **kwargs):
        wx.Panel.__init__(self, parent, *args, **kwargs)

        self.settings = settings
        self.main_frame = self.FindWindowByName('PipelineFrame')

        self._create_layout()

        self.data_dir = None
        self.output_dir = None
        self.profiles_dir = None
        self.analysis_dir = None
        self.fprefix = ''
        self.num_rproc = ''
        self.num_aproc = ''
        self.images_loaded = ''
        self.images_averaged = ''
        self.current_exps = ''
        self.exp_being_processed = ''
        self.processed_exps = ''

        self._status_timer = wx.Timer()
        self._status_timer.Bind(wx.EVT_TIMER, self._on_status_timer)
        self._status_timer.Start(1000)

    def _FromDIP(self, size):
        # This is a hack to provide easy back compatibility with wxpython < 4.1
        try:
            return self.FromDIP(size)
        except Exception:
            return size

    def _create_layout(self):
        ctrl_sizer = wx.StaticBoxSizer(wx.VERTICAL, self, 'Pipeline Control')
        ctrl_box = ctrl_sizer.GetStaticBox()

        status_lbl = wx.StaticText(ctrl_box, label='Pipeline Status:')
        self.pipeline_status = wx.StaticText(ctrl_box, label='Not Started')

        self.start_btn = wx.Button(ctrl_box, label='Start')
        self.abort_btn = wx.Button(ctrl_box, label='Abort')

        self.start_btn.Bind(wx.EVT_BUTTON, self._on_start)
        self.abort_btn.Bind(wx.EVT_BUTTON, self._on_abort)
        self.abort_btn.Disable()

        btn_sizer = wx.BoxSizer(wx.HORIZONTAL)
        btn_sizer.Add(self.start_btn, flag=wx.RIGHT, border=self._FromDIP(5))
        btn_sizer.Add(self.abort_btn)

        ctrl_grid_sizer = wx.GridBagSizer(vgap=self._FromDIP(5), hgap=self._FromDIP(5))
        ctrl_grid_sizer.Add(status_lbl, (0,0), flag=wx.ALIGN_CENTER_VERTICAL)
        ctrl_grid_sizer.Add(self.pipeline_status, (0,1), span=(0,2),
            flag=wx.ALIGN_CENTER_VERTICAL)
        ctrl_grid_sizer.Add(btn_sizer, (1,0), span=(0,3),
            flag=wx.ALIGN_CENTER_HORIZONTAL)


        # raver_label = wx.StaticText(ctrl_box, label='Radial average processes:')
        # self.raver_ctrl = wx.StaticText(ctrl_box, size=self._FromDIP((40, -1)))
        # raver_add = wx.Button(ctrl_box, label='Add')

        analysis_label = wx.StaticText(ctrl_box, label='Analysis processes:')
        self.analysis_ctrl = wx.StaticText(ctrl_box, size=self._FromDIP((40, -1)))
        analysis_add = wx.Button(ctrl_box, label='Add')

        # raver_add.Bind(wx.EVT_BUTTON, self._add_raver_proc)
        analysis_add.Bind(wx.EVT_BUTTON, self._add_analysis_proc)

        ctrl_grid_sizer_2 = wx.FlexGridSizer(cols=3, vgap=self._FromDIP(5),
            hgap=self._FromDIP(5))
        # ctrl_grid_sizer_2.Add(raver_label, flag=wx.ALIGN_CENTER_VERTICAL)
        # ctrl_grid_sizer_2.Add(self.raver_ctrl, flag=wx.ALIGN_CENTER_VERTICAL)
        # ctrl_grid_sizer_2.Add(raver_add, flag=wx.ALIGN_CENTER_VERTICAL)
        ctrl_grid_sizer_2.Add(analysis_label, flag=wx.ALIGN_CENTER_VERTICAL)
        ctrl_grid_sizer_2.Add(self.analysis_ctrl, flag=wx.ALIGN_CENTER_VERTICAL)
        ctrl_grid_sizer_2.Add(analysis_add, flag=wx.ALIGN_CENTER_VERTICAL)


        fprefix_label = wx.StaticText(ctrl_box, label='File prefix:')
        self.fprefix_ctrl = wx.TextCtrl(ctrl_box, style=wx.TE_READONLY)
        fprefix_btn = wx.Button(ctrl_box, label='Set')
        fprefix_btn.Bind(wx.EVT_BUTTON, self._on_fprefix_set)

        dd_label = wx.StaticText(ctrl_box, label='Data directory:')
        self.data_dir_ctrl = wx.TextCtrl(ctrl_box, style=wx.TE_READONLY)
        dd_btn = wx.Button(ctrl_box, label='Select')
        dd_btn.Bind(wx.EVT_BUTTON, self._on_dd_select)

        od_label = wx.StaticText(ctrl_box, label='Output directory:')
        self.output_dir_ctrl = wx.TextCtrl(ctrl_box, style=wx.TE_READONLY)
        od_btn = wx.Button(ctrl_box, label='Select')
        od_btn.Bind(wx.EVT_BUTTON, self._on_od_select)

        pd_label = wx.StaticText(ctrl_box, label='Profiles directory:')
        self.profiles_dir_ctrl = wx.TextCtrl(ctrl_box, style=wx.TE_READONLY)
        pd_btn = wx.Button(ctrl_box, label='Select')
        pd_btn.Bind(wx.EVT_BUTTON, self._on_pd_select)

        ad_label = wx.StaticText(ctrl_box, label='Analysis directory:')
        self.analysis_dir_ctrl = wx.TextCtrl(ctrl_box, style=wx.TE_READONLY)
        ad_btn = wx.Button(ctrl_box, label='Select')
        ad_btn.Bind(wx.EVT_BUTTON, self._on_ad_select)

        dir_grid_sizer = wx.FlexGridSizer(cols=3, vgap=self._FromDIP(5),
            hgap=self._FromDIP(5))
        dir_grid_sizer.Add(fprefix_label, flag=wx.ALIGN_CENTER_VERTICAL)
        dir_grid_sizer.Add(self.fprefix_ctrl, flag=wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
        dir_grid_sizer.Add(fprefix_btn, flag=wx.ALIGN_CENTER_VERTICAL)
        dir_grid_sizer.Add(dd_label, flag=wx.ALIGN_CENTER_VERTICAL)
        dir_grid_sizer.Add(self.data_dir_ctrl, flag=wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
        dir_grid_sizer.Add(dd_btn, flag=wx.ALIGN_CENTER_VERTICAL)
        dir_grid_sizer.Add(od_label, flag=wx.ALIGN_CENTER_VERTICAL)
        dir_grid_sizer.Add(self.output_dir_ctrl, flag=wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
        dir_grid_sizer.Add(od_btn, flag=wx.ALIGN_CENTER_VERTICAL)
        dir_grid_sizer.Add(pd_label, flag=wx.ALIGN_CENTER_VERTICAL)
        dir_grid_sizer.Add(self.profiles_dir_ctrl, flag=wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
        dir_grid_sizer.Add(pd_btn, flag=wx.ALIGN_CENTER_VERTICAL)
        dir_grid_sizer.Add(ad_label, flag=wx.ALIGN_CENTER_VERTICAL)
        dir_grid_sizer.Add(self.analysis_dir_ctrl, flag=wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
        dir_grid_sizer.Add(ad_btn, flag=wx.ALIGN_CENTER_VERTICAL)
        dir_grid_sizer.AddGrowableCol(1)

        ctrl_sizer.Add(ctrl_grid_sizer, flag=wx.EXPAND|wx.ALL,
            border=self._FromDIP(5))
        ctrl_sizer.Add(ctrl_grid_sizer_2, flag=wx.EXPAND|wx.LEFT|wx.RIGHT|wx.BOTTOM,
            border=self._FromDIP(5))
        ctrl_sizer.Add(dir_grid_sizer, flag=wx.EXPAND|wx.LEFT|wx.RIGHT|wx.BOTTOM,
            border=self._FromDIP(5))


        status_sizer = wx.StaticBoxSizer(wx.VERTICAL, self, 'Pipeline Status')
        status_box = status_sizer.GetStaticBox()

        loaded_label = wx.StaticText(status_box, label='Images loaded:')
        self.loaded_ctrl = wx.StaticText(status_box, size=self._FromDIP((40, -1)))

        averaged_label = wx.StaticText(status_box, label='Images averaged:')
        self.averaged_ctrl = wx.StaticText(status_box, size=self._FromDIP((40, -1)))

        exp_current_label = wx.StaticText(status_box, label='Ongoing experiments:')
        self.exp_current_ctrl = wx.StaticText(status_box, size=self._FromDIP((40, -1)))

        exp_in_proc_label = wx.StaticText(status_box, label='Experments processing:')
        self.exp_in_proc_ctrl = wx.StaticText(status_box, size=self._FromDIP((40, -1)))

        exp_proc_label = wx.StaticText(status_box, label='Finished experiments:')
        self.exp_proc_ctrl = wx.StaticText(status_box, size=self._FromDIP((40, -1)))

        status_grid_sizer = wx.FlexGridSizer(cols=2, vgap=self._FromDIP(5),
            hgap=self._FromDIP(5))
        status_grid_sizer.Add(loaded_label, flag=wx.ALIGN_CENTER_VERTICAL)
        status_grid_sizer.Add(self.loaded_ctrl, flag=wx.ALIGN_CENTER_VERTICAL)
        status_grid_sizer.Add(averaged_label, flag=wx.ALIGN_CENTER_VERTICAL)
        status_grid_sizer.Add(self.averaged_ctrl, flag=wx.ALIGN_CENTER_VERTICAL)
        status_grid_sizer.Add(exp_current_label, flag=wx.ALIGN_CENTER_VERTICAL)
        status_grid_sizer.Add(self.exp_current_ctrl, flag=wx.ALIGN_CENTER_VERTICAL)
        status_grid_sizer.Add(exp_in_proc_label, flag=wx.ALIGN_CENTER_VERTICAL)
        status_grid_sizer.Add(self.exp_in_proc_ctrl, flag=wx.ALIGN_CENTER_VERTICAL)
        status_grid_sizer.Add(exp_proc_label, flag=wx.ALIGN_CENTER_VERTICAL)
        status_grid_sizer.Add(self.exp_proc_ctrl, flag=wx.ALIGN_CENTER_VERTICAL)

        status_sizer.Add(status_grid_sizer, flag=wx.EXPAND|wx.ALL)

        top_sizer = wx.BoxSizer(wx.VERTICAL)
        top_sizer.Add(ctrl_sizer, flag=wx.EXPAND)
        top_sizer.Add(status_sizer, flag=wx.EXPAND)

        self.SetSizer(top_sizer)

    def _on_start(self, evt):
        if os.path.exists(self.settings['raw_settings_file']):
            wx.CallAfter(self.main_frame._start_pipeline)
        else:
            msg = ("Cannot find RAW settings file:\n{}".format(
                self.settings['raw_settings_file']))
            dlg = wx.MessageDialog(parent=self,
                message=msg,
                caption='Cannot start pipeline without a RAW config file',
                style=wx.OK|wx.ICON_QUESTION)
            dlg.ShowModal()
            dlg.Destroy()

    def _on_abort(self, evt):
        if self.main_frame.pipeline_abort_event is not None:
            self.main_frame.pipeline_abort_event.set()

    def _add_raver_proc(self, evt):
        if self.main_frame.pipeline_cmd_q is not None:
            self.main_frame.pipeline_cmd_q.append(['add_reduction_process', [], {}])

    def _add_analysis_proc(self, evt):
        if self.main_frame.pipeline_cmd_q is not None:
            self.main_frame.pipeline_cmd_q.append(['add_analysis_process', [], {}])

    def update_pipeline_settings(self):
        if self.main_frame.pipeline_cmd_q is not None:
            self.main_frame.pipeline_cmd_q.append(['update_pipeline_settings',
                [], {kw : self.settings[kw] for kw in self.settings}])

    def _on_dd_select(self, evt):
        self._on_select_dir('data')

    def _on_od_select(self, evt):
        self._on_select_dir('output')

    def _on_pd_select(self, evt):
        self._on_select_dir('profiles')

    def _on_ad_select(self, evt):
        self._on_select_dir('analysis')

    def _on_select_dir(self, dir_type):
        path = self.main_frame.create_file_dialog('dir',
            dir_msg='Select {} directory'.format(dir_type))

        if path is not None:
            if dir_type == 'data':
                self.data_dir_ctrl.SetValue(path)
                if self.main_frame.pipeline_cmd_q is not None:
                    self.main_frame.pipeline_cmd_q.append(['set_data_dir',
                        [path], {}])
            elif dir_type == 'output':
                self.output_dir_ctrl.SetValue(path)
                if self.main_frame.pipeline_cmd_q is not None:
                    self.main_frame.pipeline_cmd_q.append(['set_output_dir',
                        [path], {}])
            elif dir_type == 'profiles':
                self.profiles_dir_ctrl.SetValue(path)
                if self.main_frame.pipeline_cmd_q is not None:
                    self.main_frame.pipeline_cmd_q.append(['set_profiles_dir',
                        [path], {}])
            elif dir_type == 'analysis':
                self.analysis_dir_ctrl.SetValue(path)
                if self.main_frame.pipeline_cmd_q is not None:
                    self.main_frame.pipeline_cmd_q.append(['set_analysis_dir',
                        [path], {}])

    def _on_fprefix_set(self, evt):
        fprefix = self.fprefix_ctrl.GetValue()

        if self.main_frame.pipeline_cmd_q is not None:
            self.main_frame.pipeline_cmd_q.append(['set_fprefix', [fprefix], {}])

    def _on_status_timer(self, evt):
        self._get_status()

    def _get_status(self):
        pt = self.main_frame.pipeline_thread

        if pt is not None:
            pl = self.main_frame.pipeline_ret_lock
            pl.acquire()
            data_dir = pt.data_dir
            output_dir = pt.output_dir
            profiles_dir = pt.profiles_dir
            analysis_dir = pt.analysis_dir
            fprefix = pt.fprefix
            num_rproc = str(len(pt.reduction_processes))
            num_aproc = str(len(pt.analysis_processes))
            images_loaded = str(pt.get_num_loaded())
            images_averaged = str(pt.get_num_averaged())
            current_exps = str(pt.exp_total)
            exp_being_processed = str(pt.exp_being_processed)
            processed_exps = str(pt.exp_processed)
            pl.release()

            if data_dir != self.data_dir:
                self.data_dir_ctrl.SetValue(data_dir)
                self.data_dir = data_dir

            if output_dir != self.output_dir:
                self.output_dir_ctrl.SetValue(output_dir)
                self.output_dir = output_dir

            if profiles_dir != self.profiles_dir:
                self.profiles_dir_ctrl.SetValue(profiles_dir)
                self.profiles_dir = profiles_dir

            if analysis_dir != self.analysis_dir:
                self.analysis_dir_ctrl.SetValue(analysis_dir)
                self.analysis_dir = analysis_dir

            if fprefix != self.fprefix:
                self.fprefix_ctrl.SetValue(fprefix)
                self.fprefix = fprefix

            # if num_rproc != self.num_rproc:
            #     self.raver_ctrl.SetLabel(num_rproc)
            #     self.num_rproc = num_rproc

            if num_aproc != self.num_aproc:
                self.analysis_ctrl.SetLabel(num_aproc)
                self.num_aproc = num_aproc

            if images_loaded != self.images_loaded:
                self.loaded_ctrl.SetLabel(images_loaded)
                self.images_loaded = images_loaded

            if images_averaged != self.images_averaged:
                self.averaged_ctrl.SetLabel(images_averaged)
                self.images_averaged = images_averaged

            if current_exps != self.current_exps:
                self.exp_current_ctrl.SetLabel(current_exps)
                self.current_exps = current_exps

            if exp_being_processed != self.exp_being_processed:
                self.exp_in_proc_ctrl.SetLabel(exp_being_processed)
                self.exp_being_processed = exp_being_processed

            if processed_exps != self.processed_exps:
                self.exp_proc_ctrl.SetLabel(processed_exps)
                self.processed_exps = processed_exps

    def on_exit(self):
        self._status_timer.Stop()
