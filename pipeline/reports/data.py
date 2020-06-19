
import collections

GuinierData = collections.namedtuple('Guinier', ['Rg', 'I0', 'Rg_err',
    'I0_err', 'n_min', 'n_max', 'q_min', 'q_max', 'qRg_min', 'qRg_max', 'r_sq'],
    defaults=[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1])

AbsMWData = collections.namedtuple('Abs_MW', ['MW', 'Buffer_density',
    'Protein_density', 'Partial_specific_volume'], defaults=[-1, -1, -1, -1])
I0MWData = collections.namedtuple('I0_MW', ['MW'], defaults=[-1])
VpMWData = collections.namedtuple('PV_MW', ['MW', 'Density', 'q_max',
    'Porod_volume_corrected', 'Porod_volume', 'cutoff'], defaults=[-1, -1, -1,
    -1, -1, ''])
VcMWData = collections.namedtuple('Vc_MW', ['MW', 'Type', 'q_max',
    'Volume_of_correlation', 'cutoff'], defaults=[-1, '', -1, -1, ''])
SSMWData = collections.namedtuple('Shape_and_size_MW', ['MW', 'Dmax', 'Shape'],
    defaults=[-1, -1, ''])
BayesMWData = collections.namedtuple('Bayes_MW', ['MW', 'Probability',
    'Confidence_interval_lower', 'Confidence_interval_upper',
    'Confidence_interval_probability'], defaults=[-1, -1, -1, -1, -1])

BIFTData = collections.namedtuple('BIFT', ['Dmax', 'Rg', 'I0', 'Dmax_err',
    'Rg_err', 'I0_err', 'Chi_sq', 'q_min', 'q_max', 'Evidence', 'log_alpha',
    'Evidence_err', 'log_alpha_err'], defaults=[-1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1])

GNOMData = collections.namedtuple('GNOM', ['Dmax', 'Rg', 'I0', 'Rg_err',
    'I0_err', 'Chi_sq', 'Total_estimate', 'Quality', 'q_min', 'q_max'],
    defaults=[-1, -1, -1, -1, -1, -1, -1, '', -1, -1])

Metadata = collections.namedtuple('Metadata', ['Sample_to_detector_distance',
    'Wavelength', 'Exposure_time', 'Exposure_period', 'Flow_rate',
    'Detector', 'Instrument', 'Absolute_scale'], defaults=[-1, -1, -1, -1, -1,
    '', '', False])

class SECData(object):
    """
    The goal of this class is to contain all of the information about a SEC
    experiment in a single python object, with easily accessible info.
    """

    def __init__(self, secm):
        """
        Takes as input a RAW SECM object, and extracts parameters from there.
        """

        self.filename = secm.getParameter('filename')

        self.buffer_range = secm.buffer_range
        self.sample_range = secm.sample_range

        self.rg = secm.rg_list
        self.rg_err = secm.rger_list
        self.i0 = secm.i0_list
        self.i0_err = secm.i0er_list
        self.vpmw = secm.vpmw_list
        self.vcmw = secm.vcmw_list
        self.vcmw_err = secm.vcmwer_list

        self.has_calc_data = secm.calc_has_data

        self.time = secm.time
        self.frames = secm.plot_frame_list

        if secm.baseline_subtracted_sasm_list:
            self.total_i = secm.total_i_bcsub
            self.mean_i = secm.mean_i_bcsub

            self.baseline_corrected = True
            self.subtracted = True
        elif secm.subtracted_sasm_list:
            self.total_i = secm.total_i_sub
            self.mean_i = secm.mean_i_sub

            self.baseline_corrected = False
            self.subtracted = True
        else:
            self.total_i = secm.total_i
            self.mean_i = secm.mean_i

            self.baseline_corrected = False
            self.subtracted = False

class SAXSData(object):
    """
    The goal of this class is to contain all the information about a single
    SAXS scattering profile in a single python object, with easily accessible
    info.
    """

    # Left side is key in RAW sasm analysis dict, right side is key in data namedtuple
    _guinier_trans = {
        'Rg'        : 'Rg',
        'I0'        : 'I0',
        'Rg_err'    : 'Rg_err',
        'I0_err'    : 'I0_err',
        'nStart'    : 'n_min',
        'nEnd'      : 'n_max',
        'qStart'    : 'q_min',
        'qEnd'      : 'q_max',
        'qRg_min'   : 'qRg_min',
        'qRg_max'   : 'qRg_max',
        'rsq'       : 'r_sq',
        }

    _absmw_trans = {
        'MW'                        : 'MW',
        'Density_buffer'            : 'Buffer_density',
        'Density_dry_protein'       : 'Protein_density',
        'Partial_specific_volume'   : 'Partial_specific_volume',
        }

    _i0mw_trans = {
        'MW'    : 'MW',
        }

    _vpmw_trans = {
        'MW'                : 'MW',
        'Density'           : 'Density',
        'Q_max'             : 'q_max',
        'VPorod_corrected'  : 'Porod_volume_corrected',
        'VPorod'            : 'Porod_volume',
        'Cutoff'            : 'cutoff'
        }

    _vcmw_trans = {
        'MW'        : 'MW',
        'Type'      : 'Type',
        'Q_max'     : 'q_max',
        'Vcor'      : 'Volume_of_correlation',
        'Cutoff'    : 'cutoff',
        }

    _ssmw_trans = {
        'MW'    : 'MW',
        'Dmax'  : 'Dmax',
        'Shape' : 'Shape',
        }

    _bayesmw_trans = {
        'MW'                            : 'MW',
        'MWProbability'                 : 'Probability',
        'ConfidenceIntervalLower'       : 'Confidence_interval_lower',
        'ConfidenceIntervalUpper'       : 'Confidence_interval_upper',
        'ConfidenceIntervalProbability' : 'Confidence_interval_probability',
        }

    _mw_trans = collections.OrderedDict([
        ('Absolute',            'Absolute'),
        ('I(0)Concentration',   'Reference'),
        ('PorodVolume',         'Porod_volume'),
        ('VolumeOfCorrelation', 'Volume_of_correlation'),
        ('ShapeAndSize',        'Shape_and_size'),
        ('DatmwBayes',          'Bayesian'),
        ])

    _mw_methods = {
        'Absolute'              : (_absmw_trans, AbsMWData),
        'Reference'             : (_i0mw_trans, I0MWData),
        'Porod_volume'          : (_vpmw_trans, VpMWData),
        'Volume_of_correlation' : (_vcmw_trans, VcMWData),
        'Shape_and_size'        : (_ssmw_trans, SSMWData),
        'Bayesian'              : (_bayesmw_trans, BayesMWData),
        }

    _bift_trans = {
        'Dmax'              : 'Dmax',
        'Real_Space_Rg'     : 'Rg',
        'Real_Space_I0'     : 'I0',
        'Dmax_Err'          : 'Dmax_err',
        'Real_Space_Rg_Err' : 'Rg_err',
        'Real_Space_I0_Err' : 'I0_err',
        'ChiSquared'        : 'Chi_sq',
        'qStart'            : 'q_min',
        'qEnd'              : 'q_max',
        'Evidence'          : 'Evidence',
        'LogAlpha'          : 'log_alpha',
        'Evidence_Err'      : 'Evidence_err',
        'LogAlpha_Err'      : 'log_alpha_err',
        }

    _gnom_trans = {
        'Dmax'                      : 'Dmax',
        'Real_Space_Rg'             : 'Rg',
        'Real_Space_I0'             : 'I0',
        'Real_Space_Rg_Err'         : 'Rg_err',
        'Real_Space_I0_Err'         : 'I0_err',
        'GNOM_ChiSquared'           : 'Chi_sq',
        'Total_Estimate'            : 'Total_estimate',
        'GNOM_Quality_Assessment'   : 'Quality',
        'qStart'                    : 'q_min',
        'qEnd'                      : 'q_max',
        }

    _calib_trans = {
        'Sample-to-detector distance (mm)'  : 'Sample_to_detector_distance',
        'Wavelength (A)'                    : 'Wavelength',
        }

    _counters_trans = {
        'Flow rate (ml/min)'        : 'Flow_rate',
        'Exposure time/frame (s)'   : 'Exposure_time',
        'Exposure_period/frame_s'   : 'Exposure_period',
        'Instrument'                : 'Instrument',
        }

    _metadata_trans = {
        'Detector'  : 'Detector',
        }


    def __init__(self, sasm):
        """
        Takes as input a RAW SASM object, and extracts parameters from there.
        """
        self.filename = sasm.getParameter('filename')

        self.q = sasm.getQ()
        self.i = sasm.getI()
        self.err = sasm.getErr()

        all_params = sasm.getAllParameters()

        if 'analysis' in all_params:
            self._analysis_data = sasm.getParameter('analysis')
        else:
            self._analysis_data = {}

        if 'history' in all_params:
            self._history_data = sasm.getParameter('history')
        else:
            self._history_data = {}

        if 'counters' in all_params:
            self._counters_data = sasm.getParameter('counters')
        else:
            self._counters_data = {}

        if 'calibration_params' in all_params:
            self._calibration_data = sasm.getParameter('calibration_params')
        else:
            self._calibration_data = {}

        if 'metadata' in all_params:
            self._metadata = sasm.getParameter('metadata')
        else:
            self._metadata = {}

        if 'normalizations' in all_params:
            self._normalization_data = sasm.getParameter('normalizations')
        else:
            self._normalization_data = {}

        self._extract_analysis_data()
        self._extract_metadata()

    def _extract_analysis_data(self):
        """
        Extracts data from the sasm analysis dictionary into sets of named tuples
        defined at the top of this method.

        If you want to add more modify data types, add or modify the appropriate
        named tuple, and then add or modify the translation method in the class
        definition, such as _guinier_trans. This should transparently handle
        missing keys by setting the value to -1 or ''.
        """

        # Grab Guinier data
        if 'guinier' in self._analysis_data:
            guinier_analysis = self._analysis_data['guinier']

            data_dict = {}

            for key, value in self._guinier_trans.items():
                if key in guinier_analysis:
                    val = guinier_analysis[key]

                    if key != 'nStart' and key != 'nEnd':
                        val = float(val)
                    else:
                        val = int(val)

                    data_dict[value] = val

        self.guinier_data = GuinierData(**data_dict)

        #Grab MW data
        self.mw_data = collections.OrderedDict()

        for key, value in self._mw_trans.items():
            trans_dict, data_method = self._mw_methods[value]

            data_dict = {}

            if 'molecularWeight' in self._analysis_data:
                mw_analysis = self._analysis_data['molecularWeight']

                if key in mw_analysis:
                    for data_key, data_item in trans_dict.items():
                        if data_key in mw_analysis[key]:
                            data_dict[data_item] = mw_analysis[key][data_key]

            data_tuple = data_method(**data_dict)

            self.mw_data[value] = data_tuple

        # Grab BIFT data
        if 'BIFT' in self._analysis_data:
            bift_analysis = self._analysis_data['BIFT']

            data_dict = {}

            for key, value in self._bift_trans.items():
                if key in bift_analysis:
                    data_dict[value] = float(bift_analysis[key])

        self.bift_data = BIFTData(**data_dict)

        # Grab GNOM data
        if 'GNOM' in self._analysis_data:
            gnom_analysis = self._analysis_data['GNOM']

            data_dict = {}

            for key, value in self._gnom_trans.items():
                if key in gnom_analysis:
                    val = gnom_analysis[key]

                    if key != 'GNOM_Quality_Assessment':
                        val = float(val)

                    data_dict[value] = val

        self.gnom_data = GNOMData(**data_dict)


    def _extract_metadata(self):
        """
        Extracts metadata from the sasm header, calibration_params, metadata,
        and normalizations dictionaries into a named tuple defined at the
        top of the method.

        See _extract_analysis_data for how to modify what's read in.
        """
        metadata_dict = {}

        for key, value in self._calib_trans.items():
            if key in self._calibration_data:
                metadata_dict[value] = self._calibration_data[key]

        for key, value in self._counters_trans.items():
            if key in self._counters_data:
                metadata_dict[value] = self._counters_data[key]

        for key, value in self._metadata_trans.items():
            if key in self._metadata:
                metadata_dict[value] = self._metadata[key]

        if 'Absolute_scale' in self._normalization_data:
            metadata_dict['Absolute_scale'] = True

        self.metadata = Metadata(**metadata_dict)

class IFTData(object):
    """
    The goal of this class is to contain all the information about a IFT
    in a single python ojbect, with easily accessible info.
    """

    def __init__(self, iftm):
        """
        Takes as input a RAW IFTM object, and extracts parameters from there.
        """

        self.filename = iftm.getParameter('filename')

        self.r = iftm.r
        self._p_orig = iftm.p
        self._p_err_orig = iftm.err

        self.q = iftm.q_orig
        self.i = iftm.i_orig
        self.i_err = iftm.err_orig

        self.i_fit = iftm.i_fit

        self.q_extrap = iftm.q_extrap
        self.i_extrap = iftm.i_extrap

        self.dmax = iftm.getParameter('dmax')
        self.rg = iftm.getParameter('rg')
        self.i0 = iftm.getParameter('i0')
        self.rg_err = iftm.getParameter('rger')
        self.i0_err = iftm.getParameter('i0er')
        self.chi_sq = iftm.getParameter('chisq')

        self.type = iftm.getParameter('algorithm')

        if self.type == 'BIFT':
            self.dmax_err = iftm.getParameter('dmaxer')
        elif self.type == 'GNOM':
            self.total_estimate = iftm.getParameter('TE')
            self.quality = iftm.getParameter('quality')

        self.p = self._p_orig/self.i0
        self.p_err = self._p_err_orig/self.i0
