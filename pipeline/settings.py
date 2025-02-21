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

import copy
import json
import traceback
import os
import logging

if __name__ != '__main__':
    logger = logging.getLogger(__name__)

import numpy as np

class Settings(object):

    def __init__(self):

        # Settings entries are self describing, items have the following values:
        # [Actual useful value, type for GUI, section for GUI, label for GUI]
        # If the type is a choice, the GUI label is actually (label, choices)

        self._params = {
            'image_exts'                : [['.tiff', '.tif'], 'list', 'Data', 'Image extensions:'], #used to filter files in the monitor and load thread
            'r_procs'                   : [1, 'int', 'Startup', 'Radial averaging processes:'], #Number of reduction processes to start by default ***Note: for now limited to 1
            'a_procs'                   : [1, 'int', 'Startup', 'Analysis processes:'], #Number of analysis processes to start by default
            'raw_settings_file'         : ['', 'str', 'Startup', 'RAW settings file:'], #RAW settings file to give to the pipeline on startup
            'default_output_dir'        : ['..', 'str', 'Data', 'Default output directory:'], #The default output directory location, relative to the data directory
            'use_default_output_dir'    : [True, 'bool', 'Data', 'Use default output directory'], #If the default output directory location should be used, or if it will be set manually
            'default_profiles_dir'      : ['profiles', 'str', 'Data', 'Default profiles directory:'], #The default profiles directory to save radially averaged .dat files to, relative to the output directory
            'use_default_profiles_dir'  : [True, 'bool', 'Data', 'Use default profiles directory'], #If the default profiles directory location should be used, or if it will be set manually
            'default_analysis_dir'      : ['analysis', 'str', 'Data', 'Default analysis directory:'], #The default analysis directory to save radially averaged .dat files to, relative to the output directory
            'use_default_analysis_dir'  : [True, 'bool', 'Data', 'Use default analysis directory'], #If the default analysis directory location should be used, or if it will be set manually
            'save_raver_profiles'       : [True, 'bool', 'Data', 'Save radially averaged profiles'], #Save profiles after radial averaging as individual .dat files
            'data_source'               : ['Eiger Stream', 'choice', 'Data', ('Data source:', ['Files', 'Eiger Stream'])], #If data source is files on disk, or an Eiger stream client
            'eiger_stream_ip'           : ['164.54.204.140', 'str', 'Data', 'Eiger IP:'], #Eiger DCU IP
            'eiger_stream_port'         : ['9999', 'str', 'Data', 'Eiger port:'], #Eiger DCU Stream port
            'sec_exp_timeout'           : [600, 'float', 'Analysis', 'SEC data collection timeout [s]:'], #Number of seconds without new data that causes the pipeline to assume the experiment has finished but the pipeline was not informed. -1 is infinite
            'batch_exp_timeout'         : [600, 'float', 'Analysis', 'Batch data collection timeout [s]:'], #Number of seconds without new data that causes the pipeline to assume the experiment has finished but the pipeline was not informed. -1 is infinite
            'other_exp_timeout'         : [600, 'float', 'Analysis', 'Other data collection timeout [s]:'], #Number of seconds without new data that causes the pipeline to assume the experiment has finished but the pipeline was not informed. -1 is infinite
            'sec_analysis_timeout'      : [600, 'float', 'Analysis', 'SEC analysis timeout [s]:'], #Number of seconds to wait for analysis to finish before assuming it failed. -1 is infinite
            'batch_analysis_timeout'    : [600, 'float', 'Analysis', 'Batch analysis timeout [s]:'], #Number of seconds to wait for analysis to finish before assuming it failed. -1 is infinite
            'other_analysis_timeout'    : [600, 'float', 'Analysis', 'Other analysis timeout [s]:'], #Number of seconds to wait for analysis to finish before assuming it failed. -1 is infinite
            'save_processed_data'       : [True, 'bool', 'Data', 'Save processed data'], #Save processed data in the analysis directory
            'save_report'               : [True, 'bool', 'Analysis', 'Save report'], #Save a report on the analysis
            'report_type'               : ['pdf', 'choice', 'Analysis', ('Report type:', ['pdf'])], #Report type, right now only pdf is allowed
            'denss'                     : [True, 'bool', 'Analysis', 'Run DENSS'], #Do DENSS reconstructions
            'denss_runs'                : [4, 'int', 'Analysis', 'Number of DENSS models:'], #Number of DENSS models
            'denss_mode'                : ['Fast', 'choice', 'Analysis', ('DENSS mode:', ['Fast', 'Slow', 'Membrane'])], #DENSS mode
            'denss_aver'                : [True, 'bool', 'Analysis', 'Run DENSS average'], #Average denss models
            'denss_refine'              : [False, 'bool', 'Analysis', 'Run DENSS refinement'], #Refine denss models
            'use_atsas'                 : [True, 'bool', 'Analysis', 'Use ATSAS'], #Use ATSAS
            'dammif'                    : [True, 'bool', 'Analysis', 'Run DAMMIF'], #Do dammif reconstructions
            'dam_runs'                  : [3, 'int', 'Analysis', 'Number of DAMMIF models:'], #Number of dammif models
            'dam_mode'                  : ['Fast', 'choice', 'Analysis', ('DAMMIF mode:', ['Fast', 'Slow'])], #Dammif mode
            'damaver'                   : [True, 'bool', 'Analysis', 'Run DAMAVER'], #Run damaver
            'damclust'                  : [False, 'bool', 'Analysis', 'Run DAMCLUST'], #Run damclust
            'dam_refine'                : [False, 'bool', 'Analysis', 'Run DAMMIN refinement'], #Run a dammin refinement
            }

    def __getitem__(self, key):
        return self._params[key][0]

    def __setitem__(self, key, value):
        self._params[key][0] = value

    def __iter__(self):
        return iter(self._params)

    def get(self, key):
        return self._params[key][0]

    def get_all(self):
        return self._params

    def get_full(self, key):
        return self._params[key]


def load_settings(settings, path):
    loaded_param = read_settings(path)

    if loaded_param is None:
        success = False

    else:
        for each_key in loaded_param:
            if each_key in settings:
                settings[each_key] = copy.copy(loaded_param[each_key])

        default_settings = Settings()

        for key in default_settings:
            if key not in loaded_param:
                settings[key] = default_settings[key]

        success = True

    msg = ''

    return success, msg

def save_settings(settings, savepath, save_backup=True, backup_path=''):

    exclude_keys = []

    save_dict = {}

    for each_key in settings:
        if each_key not in exclude_keys:
            save_dict[each_key] = settings[each_key]

    save_dict = copy.deepcopy(save_dict)

    success = write_settings(savepath, save_dict)

    if success:
        dummy_settings = Settings()

        test_settings = load_settings(dummy_settings, savepath)

        if test_settings is None:
            os.remove(savepath)
            success = False

    if success and save_backup:
        backup_file = os.path.join(backup_path, 'backup.pcfg')

        success = write_settings(backup_file, save_dict)

        if success:
            dummy_settings = Settings()

            test_settings = load_settings(dummy_settings, backup_file)

            if test_settings is None:
                os.remove(backup_file)
                success = False

    return success

def write_settings(filename, settings):
    try:
        with open(filename, 'w', encoding='utf-8') as f:
            settings_str = json.dumps(settings, indent=4, sort_keys=True,
                cls=MyEncoder, ensure_ascii=False)

            f.write(settings_str)
        success = True

    except Exception:
        logger.exception('Failed to save settings to %s, error follows', filename)
        logger.exception(traceback.print_exc())
        success = False

    return success

def read_settings(filename):

    try:
        with open(filename, 'r') as f:
            settings = f.read()
        settings = dict(json.loads(settings))
    except Exception:
        logger.exception('Failed to load settings from %s, error follows', filename)
        logger.exception(traceback.print_exc())
        settings = None

    return settings

class MyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(MyEncoder, self).default(obj)
