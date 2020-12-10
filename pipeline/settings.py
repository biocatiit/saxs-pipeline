class Settings(object):

    def __init__(self):

        self._params = {
            'image_exts'                : [['.tiff', '.tif'], 'list'], #used to filter files in the monitor and load thread
            'r_procs'                   : [1, 'int'], #Number of reduction processes to start by default
            'a_procs'                   : [1, 'int'], #Number of analysis processes to start by default
            'default_output_dir'        : ['..', 'str'], #The default output directory location, relative to the data directory
            'use_default_output_dir'    : [True, 'bool'], #If the default output directory location should be used, or if it will be set manually
            'default_profiles_dir'      : ['./profiles', 'str'], #The default profiles directory to save radially averaged .dat files to, relative to the output directory
            'use_default_profiles_dir'  : [True, 'bool'], #If the default profiles directory location should be used, or if it will be set manually
            'default_analysis_dir'      : ['./analysis', 'str'], #The default analysis directory to save radially averaged .dat files to, relative to the output directory
            'use_default_analysis_dir'  : [True, 'bool'], #If the default analysis directory location should be used, or if it will be set manually
            'save_raver_profiles'       : [True, 'bool'], #Save profiles after radial averaging as individual .dat files
            'sec_exp_timeout'           : [120, 'float'], #Number of seconds without new data that causes the pipeline to assume the experiment has finished but the pipeline was not informed. -1 is infinite
            'batch_exp_timeout'         : [600, 'float'], #Number of seconds without new data that causes the pipeline to assume the experiment has finished but the pipeline was not informed. -1 is infinite
            'sec_analysis_timeout'      : [600, 'float'], #Number of seconds to wait for analysis to finish before assuming it failed. -1 is infinite
            'batch_analysis_timeout'    : [600, 'float'], #Number of seconds to wait for analysis to finish before assuming it failed. -1 is infinite
            'save_processed_data'       : [True, 'bool'], #Save processed data in the analysis directory
            'save_report'               : [True, 'bool'], #Save a report on the analysis
            'report_type'               : ['pdf', 'choice'], #Report type, right now only pdf is allowed
            'dammif'                    : [True, 'bool'], #Do dammif reconstructions
            'denss'                     : [True, 'bool'], #Do DENSS reconstructions
            'dam_runs'                  : [3, 'int'], #Number of dammif models
            'dam_mode'                  : ['Fast', 'choice'], #Dammif mode
            'damaver'                   : [True, 'bool'], #Run damaver
            'damclust'                  : [False, 'bool'], #Run damclust
            'dam_refine'                : [False, 'bool'], #Run a dammin refinement
            'denss_runs'                : [3, 'int'], #Number of DENSS models
            'denss_mode'                : ['Fast', 'choice'], #DENSS mode
            'denss_aver'                : [True, 'bool'], #Average denss models
            'denss_refine'              : [False, 'bool'], #Refine denss models
            'use_atsas'                 : [True, 'bool'], #Use ATSAS
            }

    def __getitem__(self, key):
        return self._params[key][0]

    def __setitem__(self, key, value):
        self._params[key][0] = value

    def get(self, key):
        return self._params[key][0]
