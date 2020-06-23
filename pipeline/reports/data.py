import collections

import numpy as np

from . import utils

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
    'Wavelength', 'Exposure_time', 'Exposure_period', 'Flow_rate', 'Detector',
    'Instrument', 'Absolute_scale', 'File_prefix', 'Date', 'RAW_version',
    'q_range'], defaults=[-1, -1, -1, -1, -1, '', '', False, '', '', '', ''])

class SECData(object):
    """
    The goal of this class is to contain all of the information about a SEC
    experiment in a single python object, with easily accessible info.
    """

    _calib_trans = {
        'Sample-to-detector distance (mm)'  : 'Sample_to_detector_distance',
        'Sample_Detector_Distance'          : 'Sample_to_detector_distance',
        'Wavelength (A)'                    : 'Wavelength',
        'Wavelength'                        : 'Wavelength',
        }

    _counters_trans = {
        'Flow rate (ml/min)'        : 'Flow_rate',
        'LC_flow_rate_mL/min'       : 'Flow_rate',
        'Exposure time/frame (s)'   : 'Exposure_time',
        'Exposure_time/frame_s'     : 'Exposure_time',
        'Exposure_period/frame_s'   : 'Exposure_period',
        'Instrument'                : 'Instrument',
        'File_prefix'               : 'File_prefix',
        'Date'                      : 'Date',
        }

    _metadata_trans = {
        'Detector'  : 'Detector',
        }

    def __init__(self, secm):
        """
        Takes as input a RAW SECM object, and extracts parameters from there.
        """

        self.filename = secm.getParameter('filename')

        self.buffer_range = secm.buffer_range
        self.sample_range = secm.sample_range
        self.baseline_start_range = secm.baseline_start_range
        self.baseline_end_range = secm.baseline_end_range
        self.baseline_type = secm.baseline_type

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

        self.get_metadata(secm)
        self.get_efa_data(secm)

    def get_metadata(self, secm):
        metadata_dict = {}

        first_prof = secm.getSASM(0)

        all_params = first_prof.getAllParameters()

        metadata_dict = {}

        if 'calibration_params' in all_params:
            calibration_params = first_prof.getParameter('calibration_params')

            for key, value in self._calib_trans.items():
                if key in calibration_params:
                    metadata_dict[value] = calibration_params[key]

        if 'counters' in all_params:
            counters = first_prof.getParameter('counters')

            for key, value in self._counters_trans.items():
                if key in counters:
                    metadata_dict[value] = counters[key]

        if 'metadata' in all_params:
            metadata = first_prof.getParameter('metadata')

            for key, value in self._metadata_trans.items():
                if key in metadata:
                    metadata_dict[value] = metadata[key]

        if 'normalizations' in all_params:
            normalizations = first_prof.getParameter('normalizations')

            if 'Absolute_scale' in normalizations:
                metadata_dict['Absolute_scale'] = True

        if 'raw_version' in all_params:
            metadata_dict['RAW_version'] = all_params['raw_version']

        q_i = first_prof.getQ()[0]
        q_f = first_prof.getQ()[-1]
        metadata_dict['q_range'] = '{} to {}'.format(utils.text_round(q_i, 4),
            utils.text_round(q_f, 2))

        self.metadata = Metadata(**metadata_dict)

    def get_efa_data(self, secm):
        analysis_dict = secm.getParameter('analysis')

        if 'efa' in analysis_dict:
            efa_dict = analysis_dict['efa']

            self.efa_done = True
            self.efa_ranges = efa_dict['ranges']
            self.efa_start = efa_dict['fstart']
            self.efa_end = efa_dict['fend']
            self.efa_nsvs = efa_dict['nsvs']

        else:
            self.efa_done = False
            self.efa_ranges = []
            self.efa_start = ''
            self.efa_end = ''
            self.efa_nsvs = ''



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
        'VPorod_Corrected'  : 'Porod_volume_corrected',
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
        'Sample_Detector_Distance'          : 'Sample_to_detector_distance',
        'Wavelength (A)'                    : 'Wavelength',
        'Wavelength'                        : 'Wavelength',
        }

    _counters_trans = {
        'Flow rate (ml/min)'        : 'Flow_rate',
        'LC_flow_rate_mL/min'       : 'Flow_rate',
        'Exposure time/frame (s)'   : 'Exposure_time',
        'Exposure_time/frame_s'     : 'Exposure_time',
        'Exposure_period/frame_s'   : 'Exposure_period',
        'Instrument'                : 'Instrument',
        'File_prefix'               : 'File_prefix',
        'Date'                      : 'Date',
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
        self._extract_metadata(all_params)

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
        data_dict = {}

        if 'guinier' in self._analysis_data:
            guinier_analysis = self._analysis_data['guinier']

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
        data_dict = {}

        if 'BIFT' in self._analysis_data:
            bift_analysis = self._analysis_data['BIFT']

            for key, value in self._bift_trans.items():
                if key in bift_analysis:
                    data_dict[value] = float(bift_analysis[key])

        self.bift_data = BIFTData(**data_dict)

        # Grab GNOM data
        data_dict = {}

        if 'GNOM' in self._analysis_data:
            gnom_analysis = self._analysis_data['GNOM']

            for key, value in self._gnom_trans.items():
                if key in gnom_analysis:
                    val = gnom_analysis[key]

                    if key != 'GNOM_Quality_Assessment':
                        val = float(val)

                    data_dict[value] = val

        self.gnom_data = GNOMData(**data_dict)


    def _extract_metadata(self, all_params):
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

        if 'raw_version' in all_params:
            metadata_dict['RAW_version'] = all_params['raw_version']

        q_i = self.q[0]
        q_f = self.q[-1]
        metadata_dict['q_range'] = '{} to {}'.format(utils.text_round(q_i, 4),
            utils.text_round(q_f, 2))

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

        self.a_score = -1
        self.a_cats = -1
        self.a_interp = ''

        self.metadata = Metadata()

class EFAData(object):
    """
    Contains information about EFA that's not contained within the series analysis
    dictionary.
    """

    def __init__(self, frames, conc, chi, forward_efa, backward_efa):
        self.frames = frames
        self.conc = conc
        self.chi = chi
        self.forward_efa = forward_efa
        self.backward_efa = backward_efa
        self.n_svals = self.conc.shape[1]

class DammifData(object):
    """
    Contains information about DAMMIF run from the .csv file RAW saves.
    """

    def __init__(self, prefix, program, mode, sym, aniso, num, damaver,
        damclust, refined, nsd, nsd_std, included, res, res_err, clusters,
         rep_model):

        self.prefix = prefix
        self.program = program
        self.mode = mode
        self.sym = sym
        self.aniso = aniso
        self.num = num
        self.damaver = damaver
        self.damclust = damclust
        self.refined = refined
        self.nsd = nsd
        self.nsd_std = nsd_std
        self.included = included
        self.res = res
        self.res_err = res_err
        self.clusters = clusters
        self.rep_model = rep_model


def parse_efa_file(filename):
    with open(filename, 'r') as f:
        data = f.readlines()

    conc_idx = -1
    chi_idx = -1
    fwd_idx = -1
    bck_idx = -1
    svd_idx = -1

    for j, line in enumerate(data):
        if 'Concentration Matrix' in line:
            conc_idx = j

        elif 'Rotation Chi^2' in line:
            chi_idx = j

        elif 'Forward EFA Results' in line:
            fwd_idx = j

        elif 'Backward EFA Results' in line:
            bck_idx = j

        elif 'Singular Value Results' in line:
            svd_idx = j

    if conc_idx >= 0 and chi_idx >=0:
        conc_data = data[conc_idx+2:chi_idx-1]
    else:
        conc_data = []

    if chi_idx >=0 and fwd_idx >=0:
        chi_data =data[chi_idx+2:fwd_idx-1]
    else:
        chi_data = []

    if fwd_idx >=0 and bck_idx >=0:
        fwd_data =data[fwd_idx+2:bck_idx-1]
    else:
        fwd_data = []

    if bck_idx >=0 and svd_idx >=0:
        bck_data =data[bck_idx+2:svd_idx-1]
    else:
        bck_data = []

    frames = []
    conc = []
    chi = []
    fwd = []
    bck = []

    for line in conc_data:
        data = line.split(',')
        frame = int(float(data[0]))
        temp_data= list(map(float, data[1:]))

        frames.append(frame)
        conc.append(temp_data)

    frames = np.array(frames)
    conc = np.array(conc)

    for line in chi_data:
        data = line.split(',')
        temp_data= float(data[1])

        chi.append(temp_data)

    chi = np.array(chi)

    for line in fwd_data:
        data = line.split(',')
        temp_data= list(map(float, data[1:]))
        fwd.append(temp_data)

    fwd = np.array(fwd)

    for line in bck_data:
        data = line.split(',')
        temp_data= list(map(float, data[1:]))

        bck.append(temp_data)

    bck = np.array(bck)

    efa_data = EFAData(frames, conc, chi, fwd, bck)

    return efa_data


def parse_dammif_file(filename):
    with open(filename, 'r') as f:
        data = f.readlines()

    prefix = ''
    program = ''
    mode = ''
    sym = ''
    aniso = ''
    num = -1
    damaver = False
    damclust = False
    refined = False
    nsd = -1
    nsd_std = -1
    included = -1
    res = -1
    res_err = -1
    clusters = -1
    rep_model = -1

    for line in data:
        if 'Program used' in line:
            program = line.split(':')[-1].strip()
        elif 'Mode:' in line:
            mode = line.split(':')[-1].strip()
        elif 'Symmetry' in line:
            sym = line.split(':')[-1].strip()
        elif 'Anisometry' in line:
            aniso = line.split(':')[-1].strip()
        elif 'Total number' in line:
            num = int(line.split(':')[-1].strip())
        elif 'Used DAMAVER' in line:
            damaver = bool(line.split(':')[-1].strip())
        elif 'Reinfed with DAMMIN' in line:
            refined = bool(line.split(':')[-1].strip())
        elif 'Used DAMCLUST' in line:
            damclust = bool(line.split(':')[-1].strip())
        elif 'Mean NSD' in line:
            nsd = float(line.split(':')[-1].strip())
        elif 'Stdev. NSD' in line:
            nsd_std = float(line.split(':')[-1].strip())
        elif 'DAMAVER Included' in line:
            included = int(line.split(':')[-1].strip().split(' ')[0].strip())
        elif 'Representative mode' in line:
            rep_model = int(line.split(':')[-1].strip())
        elif 'Ensemble resolution' in line:
            res_data = line.split(':')[-1].strip()
            res = float(res_data.split('+')[0].strip())
            res_err = float(res_data.split('-')[1].strip().split(' ')[0].strip())
        elif 'Number of clusters' in line:
            clusters = int(line.split(':')[-1].strip())
        elif 'Output prefix' in line:
            prefix = line.split(':')[-1].strip()

    dammif_data = DammifData(prefix, program, mode, sym, aniso, num, damaver,
        damclust, refined, nsd, nsd_std, included, res, res_err, clusters,
        rep_model)

    return dammif_data
