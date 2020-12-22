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
import wx.lib.scrolledpanel as scrolled

pipeline_path = os.path.abspath(os.path.join('.', __file__, '..', '..', '..'))
if pipeline_path not in os.sys.path:
    os.sys.path.append(pipeline_path)

from pipeline import settings

class SettingsPanel(wx.Panel):

    def __init__(self, parent, settings, *args, **kwargs):
        wx.Panel.__init__(self, parent, *args, **kwargs)

        self.settings = settings
        self.main_frame = self.FindWindowByName('PipelineFrame')

        self._create_layout()

    def _FromDIP(self, size):
        # This is a hack to provide easy back compatibility with wxpython < 4.1
        try:
            return self.FromDIP(size)
        except Exception:
            return size

    def _create_layout(self):

        self.settings_sub_panel = SettingsSubPanel(self, self.settings)

        button_sizer = wx.BoxSizer(wx.HORIZONTAL)

        load_button = wx.Button(self, label='Load Settings')
        save_button = wx.Button(self, label='Save Settings')
        apply_button = wx.Button(self, label='Apply to Pipeline')

        load_button.Bind(wx.EVT_BUTTON, self._on_load)
        save_button.Bind(wx.EVT_BUTTON, self._on_save)
        apply_button.Bind(wx.EVT_BUTTON, self._on_apply)

        button_sizer.Add(load_button, flag=wx.LEFT|wx.TOP|wx.BOTTOM,
            border=self._FromDIP(5))
        button_sizer.Add(save_button, flag=wx.LEFT|wx.TOP|wx.BOTTOM,
            border=self._FromDIP(5))
        button_sizer.Add(apply_button, flag=wx.LEFT|wx.TOP|wx.BOTTOM,
            border=self._FromDIP(5))

        top_sizer = wx.BoxSizer(wx.VERTICAL)
        top_sizer.Add(self.settings_sub_panel, proportion=1, flag=wx.EXPAND)
        top_sizer.Add(button_sizer, flag=wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_HORIZONTAL,
            border=self._FromDIP(5))

        self.SetSizer(top_sizer)

    def set_values(self):
        self.settings_sub_panel.set_values()

    def get_values(self):
        self.settings_sub_panel.get_values()

    def _on_load(self, evt):
        fname = self.main_frame.create_file_dialog(wx.FD_OPEN)

        if fname is not None:
            success, msg = settings.load_settings(self.settings, fname)

        else:
            success = False
            msg = ''

        if not success and msg != '':
            dlg = wx.MessageDialog(self, msg, "Error loading settings",
                style=wx.ICON_ERROR|wx.OK)
            dlg.ShowModal()
            dlg.Destroy()

        self.set_values()
        self._update()

    def _on_save(self, evt):
        self._apply()

        fname = self.main_frame.create_file_dialog(wx.FD_SAVE)

        if fname is not None:
            success = settings.save_settings(self.settings, fname,
                backup_path=self.main_frame.standard_paths.GetUserLocalDataDir())

        else:
            success = True

        if not success:
            msg = "Failed to save pipeline settings."
            dlg = wx.MessageDialog(self, msg, "Error saving settings",
                style=wx.ICON_ERROR|wx.OK)
            dlg.ShowModal()
            dlg.Destroy()

    def _on_apply(self, evt):
        self._apply()

    def _apply(self):
        self.get_values()
        self._update()

    def _update(self):
        self.main_frame.status_panel.update_pipeline_settings()

    def on_exit(self):
        pass

class SettingsSubPanel(scrolled.ScrolledPanel):

    def __init__(self, parent, settings, *args, **kwargs):

        scrolled.ScrolledPanel.__init__(self, parent, *args, **kwargs)
        self.SetScrollRate(20,20)

        self.settings = settings
        self.settings_panel = self.FindWindowByName('SettingsPanel')

        self._create_layout()

        self.set_values()

    def _FromDIP(self, size):
        # This is a hack to provide easy back compatibility with wxpython < 4.1
        try:
            return self.FromDIP(size)
        except Exception:
            return size

    def _create_layout(self):
        sections = []
        for key in self.settings:
            if self.settings.get_full(key)[2] not in sections:
                sections.append(self.settings.get_full(key)[2])

        if 'Startup' in sections:
            sections.remove('Startup')
            sections.insert(0, 'Startup')

        top_sizer = wx.BoxSizer(wx.VERTICAL)

        for section in sections:
            sizer = wx.StaticBoxSizer(wx.VERTICAL, self, section)
            box = sizer.GetStaticBox()

            grid_sizer = wx.GridBagSizer(hgap=self._FromDIP(5),
                vgap=self._FromDIP(5))
            row=0

            for key in self.settings:
                values = self.settings.get_full(key)

                if values[2] == section:

                    if values[1].lower() == 'bool':
                        widget = wx.CheckBox(box, label=values[3], name=key)

                        grid_sizer.Add(widget, (row, 0), span=(1,2))

                    elif values[1].lower() == 'choice':
                        label = wx.StaticText(box, label=values[3][0])
                        ctrl = wx.Choice(box, choices=values[3][1], name=key)

                        grid_sizer.Add(label, (row, 0),
                            flag=wx.ALIGN_CENTER_VERTICAL)
                        grid_sizer.Add(ctrl, (row, 1),
                            flag=wx.ALIGN_CENTER_VERTICAL)

                    else:
                        label = wx.StaticText(box, label=values[3])
                        ctrl = wx.TextCtrl(box, name=key)

                        grid_sizer.Add(label, (row, 0),
                            flag=wx.ALIGN_CENTER_VERTICAL)
                        grid_sizer.Add(ctrl, (row, 1),
                            flag=wx.ALIGN_CENTER_VERTICAL)

                        if key == 'raw_settings_file':
                            select_btn = wx.Button(box, label='Select')
                            select_btn.Bind(wx.EVT_BUTTON, self._on_choose_raw)

                            grid_sizer.Add(select_btn, (row, 2),
                                flag=wx.ALIGN_CENTER_VERTICAL)

                    row += 1

            sizer.Add(grid_sizer, flag=wx.EXPAND)
            top_sizer.Add(sizer, border=self._FromDIP(5), flag=wx.ALL|wx.EXPAND)

        self.SetSizer(top_sizer)

    def set_values(self):
        for key in self.settings:
            ctrl = self.FindWindowByName(key)

            values = self.settings.get_full(key)

            if values[1].lower() == 'bool':
                val = values[0]
            else:
                val = str(values[0])

            if values[1].lower() == 'choice':
                ctrl.SetStringSelection(val)
            else:
                ctrl.SetValue(val)

    def get_values(self):
        for key in self.settings:
            ctrl = self.FindWindowByName(key)

            values = self.settings.get_full(key)

            if values[1].lower() == 'choice':
                val = ctrl.GetStringSelection()
            else:
                val = ctrl.GetValue()

            if values[1].lower() == 'list':
                val = eval(val)
            elif values[1].lower() == 'int':
                val = int(val)
            elif values[1].lower() == 'float':
                val = float(val)

            self.settings[key] = val

    def _on_choose_raw(self, evt):
        fname = self.settings_panel.main_frame.create_file_dialog(wx.FD_OPEN,
            name='RAW Config Files', ext='*.cfg')

        if fname is not None:
            ctrl = self.FindWindowByName('raw_settings_file')
            ctrl.SetValue(fname)
