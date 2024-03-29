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

from collections import OrderedDict, defaultdict
import tempfile
import os

from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Table, Image,
    XPreformatted, KeepTogether, TableStyle)
import matplotlib.pyplot as plt

import bioxtasraw.RAWAPI as raw

from .plots import overview_plot as overview_plot
from .plots import efa_plot as efa_plot
from .utils import text_round as text_round
from . import data
from .data import SECData as SECData

temp_files = []

def generate_report(fname, datadir, profiles, ifts, series, extra_data=None):
    """
    Inputs a list of profile data, a list of ift data, and a list of series
    data to be included in the report. Makes a PDF report.
    """

    elements = []

    overview = generate_overview(profiles, ifts, series)
    elements.extend(overview)

    exp = generate_exp_params(profiles, ifts, series)
    elements.extend(exp)

    if len(series) > 0:
        s_elements = generate_series_params(profiles, ifts, series, extra_data)
        elements.extend(s_elements)

    if len(profiles) > 0:
        guinier = generate_guinier_params(profiles, ifts, series)
        elements.extend(guinier)

        mw = generate_mw_params(profiles, ifts, series)
        elements.extend(mw)

    if len(ifts) > 0:
        if any(ift.type == 'GNOM' for ift in ifts):
            gnom = generate_gnom_params(profiles, ifts, series)
            elements.extend(gnom)

        if any(ift.type == 'BIFT' for ift in ifts):
            bift = generate_bift_params(profiles, ifts, series)
            elements.extend(bift)

    if (extra_data is not None and 'dammif' in extra_data and
        len(extra_data['dammif']) > 0 and
        any(dam_data is not None for dam_data in extra_data['dammif'])):
        dammif = generate_dammif_params(extra_data['dammif'])
        elements.extend(dammif)

    if (extra_data is not None and 'denss' in extra_data and
        len(extra_data['denss']) > 0 and
        any(dam_data is not None for dam_data in extra_data['denss'])):
        denss = generate_denss_params(extra_data['denss'])
        elements.extend(denss)

    datadir = os.path.abspath(os.path.expanduser(datadir))
    fname = '{}.pdf'.format(os.path.splitext(fname)[0])

    doc = SimpleDocTemplate(os.path.join(datadir, fname), pagesize=letter,
        leftMargin=1*inch, rightMargin=1*inch, topMargin=1*inch,
        bottomMargin=1*inch)
    doc.build(elements)

    global temp_files

    for fname in temp_files:
        if os.path.isfile(fname):
            os.remove(fname)

    temp_files = []

def generate_overview(profiles, ifts, series):
    """
    Generates the overview portion of the report. Returns a list of flowables.
    """
    styles = getSampleStyleSheet()

    elements = []

    name_list = []
    date_list = []

    if len(series) > 0:
        data_list = series
    elif len(profiles) > 0:
        data_list = profiles
    elif len(ifts) > 0:
        data_list = ifts
    else:
        data_list = []

    for s in data_list:
        if ('File_prefix' in s.metadata._fields
            and getattr(s.metadata, 'File_prefix') != ''):
            name_list.append(getattr(s.metadata, 'File_prefix'))
        else:
            name_list.append(s.filename)

        if 'Date' in s.metadata._fields and getattr(s.metadata, 'Date') != '':
            date_list.append(':'.join(getattr(s.metadata, 'Date').split(':')[:-1]))
        else:
            date_list.append('N/A')

    name_str = ', '.join(name_list)
    date_str = ', '.join(date_list)

    if 'N/A' in name_str:
        title_text = 'SAXS data overview'.format(name_str)
    else:
        title_text = '{} SAXS data overview'.format(name_str)
    ov_title = Paragraph(title_text, styles['Heading1'])

    ov_text = Paragraph('Summary:', styles['Heading2'])

    if any([date != 'N/A' for date in date_list]):
        summary_text = ('Data name(s): {}\n'
            'Collection date(s): {}'.format(name_str, date_str))
    else:
         summary_text = ('Data name(s): {}\n'.format(name_str))

    ov_summary = XPreformatted(summary_text, styles['Normal'])

    elements.append(ov_title)
    elements.append(ov_text)
    elements.append(ov_summary)


    # Make overview figure
    if len(profiles) > 0:
        has_profiles = True
        has_rg = any([prof.guinier_data.Rg > 0 and prof.guinier_data.I0 > 0 for prof in profiles])
    else:
        has_profiles = False
        has_rg = False

    if len(ifts) > 0:
        has_ifts = True
    else:
        has_ifts = False

    if len(series) > 0:
        has_series = True
    else:
        has_series = False

    if 'N/A' in name_str:
        caption = ('SAXS data summary figure.')
    else:
        caption = ('SAXS data summary figure for {}. '.format(name_str))

    if len(series) == 1:
        series_label = ('Series intensity (blue, left axis) vs. frame, and, if '
            'available, Rg vs. frame (red, right axis). Green shaded regions '
            'are buffer regions, purple shaded regions are sample regions.')
    else:
        series_label = ('Series intensity (left axis) vs. frame, and, if '
            'available, Rg vs. frame (right axis). Green shaded regions are '
            'buffer regions, purple shaded regions are sample regions.')

    profile_label = ('Scattering profile(s) on a log-lin scale.')
    guinier_label = ('Guinier fit(s) (top) and fit residuals (bottom).')
    kratky_label = ('Normalized Kratky plot. Dashed lines show where a '
        'globular system would peak.')
    ift_label = ('P(r) function(s), normalized by I(0).')


    img_width = 6

    if has_profiles and has_rg and has_ifts and has_series:
        img_height = 6

        caption = ('{} a) {} b) {} c) {} d) {} e) {}'.format(caption,
            series_label, profile_label, guinier_label, kratky_label,
            ift_label))

    elif has_profiles and has_rg and has_ifts and not has_series:
        img_height = 4

        caption = ('{} a) {} b) {} c) {} d) {}'.format(caption, profile_label,
            guinier_label, kratky_label, ift_label))

    elif has_profiles and has_rg and not has_ifts and has_series:
        img_height = 6

        caption = ('{} a) {} b) {} c) {} d) {}'.format(caption, series_label,
            profile_label, guinier_label, kratky_label))

    elif has_profiles and has_rg and not has_ifts and not has_series:
        img_height = 4

        caption = ('{} a) {} b) {} c) {}'.format(caption, profile_label,
            guinier_label, kratky_label))

    elif has_profiles and not has_rg and has_ifts and has_series:
        img_height = 4

        caption = ('{} a) {} b) {} c) {}'.format(caption,
            series_label, profile_label, ift_label))

    elif has_profiles and not has_rg and has_ifts and not has_series:
        img_height = 2

        caption = ('{} a) {} b) {}'.format(caption, profile_label,
            ift_label))

    elif has_profiles and not has_rg and not has_ifts and has_series:
        img_height = 4

        caption = ('{} a) {} b) {}'.format(caption, series_label,
            profile_label))

    elif has_profiles and not has_rg and not has_ifts and not has_series:
        img_height = 2

        caption = ('{} {}'.format(caption, profile_label))

    elif not has_profiles and has_ifts and has_series:
        img_height = 4

        caption = ('{} a) {} b) {}'.format(caption, series_label, ift_label))

    elif not has_profiles and not has_ifts and has_series:
        img_height = 2

        caption = ('{} {}'.format(caption, series_label))

    elif not has_profiles and has_ifts and not has_series:
        img_height = 2

        caption = ('{} {}'.format(caption, ift_label))

    else:
        img_height = 0

    if img_height > 0:
        ov_plot = overview_plot(profiles, ifts, series,
            img_width=img_width, img_height=img_height)

        ov_figure = make_figure(ov_plot.figure, caption, img_width, img_height,
            styles)

        elements.append(ov_figure)


    # Make overview table
    table_pairs = [
        ('', 'name'),
        ('Guinier Rg', 'gu_rg'),
        ('Guinier I(0)', 'gu_i0'),
        ('M.W. (Vp)', 'mw_vp'),
        ('M.W. (Vc)', 'mw_vc'),
        ('M.W. (S&S)', 'mw_ss'),
        ('M.W. (Bayes)', 'mw_bayes'),
        ('M.W. (Abs.)', 'mw_abs'),
        ('M.W. (Ref.)', 'mw_ref'),
        ('GNOM Dmax', 'gn_dmax'),
        ('GNOM Rg', 'gn_rg'),
        ('GNOM I(0)', 'gn_i0'),
        ('BIFT Dmax', 'b_dmax'),
        ('BIFT Rg', 'b_rg'),
        ('BIFT I(0)', 'b_i0'),
        ]

    table_dict = OrderedDict()

    table_data = []

    required_data = ['']


    for profile in profiles:

        filename = profile.filename

        if profile.guinier_data.Rg != -1:
            rg = profile.guinier_data.Rg
            rg_err = profile.guinier_data.Rg_err
            i0 = profile.guinier_data.I0
            i0_err = profile.guinier_data.I0_err

            guinier_rg = '{} +/- {}'.format(text_round(rg, 2),
                text_round(rg_err, 2))
            guinier_i0 = '{} +/- {}'.format(text_round(i0, 2),
                text_round(i0_err, 2))
        else:
            guinier_rg = ''
            guinier_i0 = ''


        mw_data = defaultdict(str)

        for mw_key in profile.mw_data:
            mw_val = profile.mw_data[mw_key]

            if mw_val.MW != -1 and mw_val.MW != '':
                val = text_round(mw_val.MW, 1)
            else:
                val = ''

            if mw_key == 'Volume_of_correlation':
                table_key = 'mw_vc'
            elif mw_key == 'Porod_volume':
                table_key = 'mw_vp'
            elif mw_key == 'Shape_and_size':
                table_key = 'mw_ss'
            elif mw_key == 'Bayesian':
                table_key = 'mw_bayes'
            elif mw_key == 'Absolute':
                table_key = 'mw_abs'
            elif mw_key == 'Reference':
                table_key = 'mw_ref'

            mw_data[table_key] = val


        if profile.gnom_data.Dmax != -1:
            dmax = profile.gnom_data.Dmax
            rg = profile.gnom_data.Rg
            rg_err = profile.gnom_data.Rg_err
            i0 = profile.gnom_data.I0
            i0_err = profile.gnom_data.I0_err

            gnom_dmax = '{}'.format(text_round(dmax, 0))

            gnom_rg = '{} +/- {}'.format(text_round(rg, 2),
                text_round(rg_err, 2))
            gnom_i0 = '{} +/- {}'.format(text_round(i0, 2),
                text_round(i0_err, 2))

        else:
            gnom_dmax = ''
            gnom_rg = ''
            gnom_i0 = ''

        if profile.bift_data.Dmax != -1:
            dmax = profile.bift_data.Dmax
            dmax_err = profile.bift_data.Dmax_err
            rg = profile.bift_data.Rg
            rg_err = profile.bift_data.Rg_err
            i0 = profile.bift_data.I0
            i0_err = profile.bift_data.I0_err

            bift_dmax = '{} +/- {}'.format(text_round(dmax, 0),
                text_round(dmax_err, 0))

            bift_rg = '{} +/- {}'.format(text_round(rg, 2),
                text_round(rg_err, 2))
            bift_i0 = '{} +/- {}'.format(text_round(i0, 2),
                text_round(i0_err, 2))

        else:
            bift_dmax = ''
            bift_rg = ''
            bift_i0 = ''

        for header, key in table_pairs:
            if key == 'name':
                value = filename
            elif key == 'gu_rg':
                value = guinier_rg
            elif key == 'gu_i0':
                value = guinier_i0
            elif key.startswith('mw'):
                value =mw_data[key]
            elif key == 'gn_dmax':
                value = gnom_dmax
            elif key == 'gn_rg':
                value = gnom_rg
            elif key == 'gn_i0':
                value = gnom_i0
            elif key == 'b_dmax':
                value = bift_dmax
            elif key == 'b_rg':
                value = bift_rg
            elif key == 'b_i0':
                value = bift_i0

            if header in table_dict:
                table_dict[header].append(value)
            else:
                table_dict[header] = [value]

    for header, values in table_dict.items():
        if header in required_data:
            table_entry = [header]
            table_entry.extend(values)
            table_data.append(table_entry)
        else:
            if any(val != '' for val in values):
                table_entry = [header]
                table_entry.extend(values)
                table_data.append(table_entry)

    if len(table_data) > 1:
        ov_table = Table(table_data, spaceBefore=0.25*inch, spaceAfter=0.1*inch)

        table_style = TableStyle(
            [('LINEBELOW', (0, 0), (-1, 0), 1, colors.black),
            ('LINEAFTER', (0, 0), (0,-1), 1, colors.black),
            ])

        ov_table.setStyle(table_style)
        ov_table.hAlign = 'LEFT'

        if 'N/A' in name_str:
            table_caption = ('SAXS data summary table. ')
        else:
            table_caption = ('SAXS data summary table for {}. '.format(name_str))

        table_text = Paragraph(table_caption, styles['Normal'])

        ov_table = KeepTogether([ov_table, table_text])

        elements.append(ov_table)

    return elements


def generate_exp_params(profiles, ifts, series):
    styles = getSampleStyleSheet()

    name_list = []

    if len(series) > 0:
        data_list = series
    elif len(profiles) > 0:
        data_list = profiles
    else:
        data_list = []

    for s in data_list:
        if isinstance(s, SECData):
            if ('File_prefix' in s.metadata._fields
            and getattr(s.metadata, 'File_prefix') != ''):
                name_list.append(getattr(s.metadata, 'File_prefix'))
            else:
                name_list.append(s.filename)
        else:
            name_list.append(s.filename)

    exp_text = Paragraph('Experimental parameters:', styles['Heading2'])

    table_pairs = [
        ('', 'name'),
        ('Date', 'Date'),
        ('Instrument', 'Instrument'),
        ('Experiment Type', 'Experiment_type'),
        ('Column', 'Column'),
        ('Mixer', 'Mixer'),
        ('Sample', 'Sample'),
        ('Buffer', 'Buffer'),
        ('Temperature [C]', 'Temperature'),
        ('Loaded volume [uL]', 'Loaded_volume'),
        ('Concentration [mg/ml]', 'Concentration'),
        ('Detector', 'Detector'),
        ('Wavelength (A)', 'Wavelength'),
        ('Camera length (m)', 'Sample_to_detector_distance'),
        ('q-measurement range (1/A)', 'q_range'),
        ('Exposure time (s)', 'Exposure_time'),
        ('Exposure period (s)', 'Exposure_period'),
        ('Flow rate (ml/min)', 'Flow_rate'),
        ('Attenuation', 'Transmission'),
        ('RAW version', 'RAW_version'),
        ('Notes', 'Notes'),
        ]

    table_dict = OrderedDict()

    table_data = []

    required_data = ['']

    for j, s in enumerate(data_list):
        for header, key in table_pairs:
            if key in s.metadata._fields:
                value = getattr(s.metadata, key)
            else:
                value = ''

            if key == 'Wavelength':
                if value != '' and value != -1:
                    value = text_round(value, 3)
                else:
                    value = ''

            elif key == 'Sample_to_detector_distance':
                if value != '' and value != -1:
                    value = float(value)/1000.
                    value = text_round(value, 3)
                else:
                    value = ''

            elif key == 'name':
                value = name_list[j]

            elif key == 'Date':
                value = ':'.join(value.split(':')[:-1])

            elif key == 'Transmission':
                if value != -1:
                    if float(value) == 1:
                        value = 'None'
                    else:
                        value = str(round(1./float(value),4))
                else:
                    value = ''

            else:
                if value != -1:
                    value = str(value)
                else:
                    value = ''

            if header in table_dict:
                table_dict[header].append(value)
            else:
                table_dict[header] = [value]


    for header, values in table_dict.items():
        if header in required_data:
            table_entry = [header]
            table_entry.extend(values)
            table_data.append(table_entry)
        else:
            if any(val != '' for val in values):
                table_entry = [header]
                table_entry.extend(values)
                table_data.append(table_entry)

    if len(table_data) > 1:
        exp_table = Table(table_data)

        table_style = TableStyle(
            [('LINEBELOW', (0, 0), (-1, 0), 1, colors.black),
            ('LINEAFTER', (0, 0), (0,-1), 1, colors.black),
            ])

        exp_table.setStyle(table_style)
        exp_table.hAlign = 'LEFT'

        elements = [exp_text, exp_table]

    else:
        elements = []

    return elements


def generate_series_params(profiles, ifts, series, extra_data):
    styles = getSampleStyleSheet()

    has_series_data = False

    name_list = []

    for s in series:
        if ('File_prefix' in s.metadata._fields
            and getattr(s.metadata, 'File_prefix') != ''):
            name_list.append(getattr(s.metadata, 'File_prefix'))
        elif s.filename != '':
            name_list.append(s.filename)
        else:
            name_list.append('N/A')

    series_text = Paragraph('Series:', styles['Heading2'])


    table_pairs = [
        ('', 'name'),
        ('Buffer range', 'buf_range'),
        ('Sample range', 'sam_range'),
        ('Baseline correction', 'baseline_type'),
        ('Baseline start range', 'baseline_start'),
        ('Baseline end range', 'baseline_end'),
        ]

    table_dict = OrderedDict()

    table_data = []

    required_data = []


    for j, s in enumerate(series):
        buffer_range = ', '.join(['{} to {}'.format(*br) for br in s.buffer_range])
        sample_range = ', '.join(['{} to {}'.format(*sr) for sr in s.sample_range])

        if s.baseline_corrected:
            baseline_type = s.baseline_type
            baseline_start = '{} to {}'.format(*s.baseline_start_range)
            baseline_end = '{} to {}'.format(*s.baseline_end_range)
        else:
            baseline_type = ''
            baseline_start = ''
            baseline_end = ''

        for header, key in table_pairs:
            if key == 'name':
                value = name_list[j]
            elif key == 'buf_range':
                value = buffer_range
            elif key == 'sam_range':
                value = sample_range
            elif key == 'baseline_type':
                value = baseline_type
            elif key == 'baseline_start':
                value = baseline_start
            elif key == 'baseline_end':
                value = baseline_end

            if header in table_dict:
                table_dict[header].append(value)
            else:
                table_dict[header] = [value]

    for header, values in table_dict.items():
        if header in required_data:
            table_entry = [header]
            table_entry.extend(values)
            table_data.append(table_entry)
        else:
            if any(val != '' for val in values):
                table_entry = [header]
                table_entry.extend(values)
                table_data.append(table_entry)

    if len(table_data) > 1:
        series_table = Table(table_data)

        table_style = TableStyle(
            [('LINEBELOW', (0, 0), (-1, 0), 1, colors.black),
            ('LINEAFTER', (0, 0), (0,-1), 1, colors.black),
            ])

        series_table.setStyle(table_style)
        series_table.hAlign = 'LEFT'
        series_table = KeepTogether([series_text, series_table])

        elements = [series_table]

        has_series_data = True

    else:
        elements = [series_text]

    #EFA
    efa_elements = []

    efa_table_pairs = [
        ('', 'name'),
        ('EFA data range', 'efa_range'),
        ('Number of components', 'nsvs'),
        ]

    efa_table_dict = OrderedDict()

    efa_table_data = []

    efa_required_data = ['', 'EFA data range', 'Number of components']

    for j, s in enumerate(series):
        if s.efa_done:
            if len(series) > 1:
                efa_title = Paragraph('{} EFA results:'.format(name_list[j]),
                    styles['Heading3'])
            else:
                efa_title = Paragraph('EFA results:', styles['Heading3'])

            # Make EFA table
            for header, key in efa_table_pairs:
                if key == 'name':
                    value = name_list[j]
                elif key == 'efa_range':
                    value = '{} to {}'.format(s.efa_start, s.efa_end)
                elif key == 'nsvs':
                    value = '{}'.format(s.efa_nsvs)
                else:
                    value = ''

                if header in efa_table_dict:
                    efa_table_dict[header].append(value)
                else:
                    efa_table_dict[header] = [value]

            for k, efa_range in enumerate(s.efa_ranges):
                header = 'Component {}'.format(k)
                value = '{} to {}'.format(*efa_range)

                if header in efa_table_dict:
                    efa_table_dict[header].append(value)
                else:
                    efa_table_dict[header] = [value]


            for header, values in efa_table_dict.items():
                if header in efa_required_data:
                    efa_table_entry = [header]
                    efa_table_entry.extend(values)
                    efa_table_data.append(efa_table_entry)
                else:
                    if any(val != '' for val in values):
                        efa_table_entry = [header]
                        efa_table_entry.extend(values)
                        efa_table_data.append(efa_table_entry)

            efa_table = Table(efa_table_data)

            table_style = TableStyle(
                [('LINEBELOW', (0, 0), (-1, 0), 1, colors.black),
                ('LINEAFTER', (0, 0), (0,-1), 1, colors.black),
                ])

            efa_table.setStyle(table_style)
            efa_table.hAlign = 'LEFT'
            # efa_table = KeepTogether([efa_table])

            # Make EFA plot

            efa_plot_panel = efa_plot(s)

            img_width = 6
            img_height = 2

            efa_caption = ('EFA deconvolution results. a) The full series '
                'intensity (blue), the selected intensity range for EFA '
                '(black), and (if available) Rg values (red). b) The selected '
                'intensity range for EFA (black), and the individual component '
                'ranges for deconvolution, with component range 0 starting at '
                'the top left, and component number increasing in descending '
                'order to the right.')

            if s.efa_extra_data:
                efa_caption = efa_caption + (' c) Mean chi^2 values between the '
                    'fit of the EFA deconvolution and the original data. d) '
                    'Area normalized concentration profiles for each component. '
                    'Colors match the component range colors in b.')

                efa_caption = efa_caption + (' e) Deconvolved scattering '
                    'profiles. Colors match the component range colors in '
                    'b and the concentration range colors in d.')

                img_height = 6

            efa_figure = make_figure(efa_plot_panel.figure, efa_caption, img_width,
                img_height, styles)

            efa_elements.append(efa_title)
            efa_elements.append(efa_table)
            efa_elements.append(efa_figure)

            has_series_data = True

    elements.extend(efa_elements)

    #REGALS
    regals_elements = []

    regals_table_pairs = [
        ('', 'name'),
        ('REGALS data range', 'regals_range'),
        ('Number of components', 'nsvs'),
        ('Convergence criteria', 'conv_type'),
        ('Max iterations', 'max_iter'),
        ('Min iterations', 'min_iter'),
        ('Convergence Tolerance', 'tol'),
        ]

    regals_table_dict = OrderedDict()

    regals_table_data = []

    regals_required_data = ['', 'REGALS data range', 'Number of components']

    regals_comp_table_dict = OrderedDict()
    regals_comp_table_data = []
    regals_comp_required_data = []

    regals_profile_component_table_pairs = [
        ('Component', 'name'),
        ('Profile type', 'type'),
        ('Prof. Lambda', 'lambda'),
        ('Prof. Control points', 'Nw'),
        ('Dmax', 'dmax'),
        ('Zero at R=0', 'is_zero_at_r0'),
        ('Zero at Dmax', 'is_zero_at_dmax'),
        ]

    regals_conc_component_table_pairs = [
        ('Concentration type', 'type'),
        ('Conc. Lambda', 'lambda'),
        ('Start', 'xmin'),
        ('End', 'xmax'),
        ('Conc. Control points', 'Nw'),
        ('Zero at start', 'is_zero_at_xmin'),
        ('Zero at end', 'is_zero_at_xmax'),
        ]

    for j, s in enumerate(series):
        if s.regals_done:
            if len(series) > 1:
                regals_title = Paragraph('{} REGALS results:'.format(name_list[j]),
                    styles['Heading3'])
            else:
                regals_title = Paragraph('REGALS results:', styles['Heading3'])

            # Make REGALS table
            for header, key in regals_table_pairs:
                if key == 'name':
                    value = name_list[j]
                elif key == 'regals_range':
                    value = '{} to {}'.format(s.regals_start, s.regals_end)
                elif key == 'nsvs':
                    value = '{}'.format(s.regals_nsvs)
                elif key == 'conv_type':
                    value = '{}'.format(s.regals_run_settings['conv_type'])
                elif key == 'max_iter':
                    value = '{}'.format(s.regals_run_settings['max_iter'])
                elif key == 'min_iter':
                    if s.regals_run_settings['conv_type'] != 'Iterations':
                        value = '{}'.format(s.regals_run_settings['min_iter'])
                    else:
                        value = ''
                elif key == 'tol':
                    if s.regals_run_settings['conv_type'] != 'Iterations':
                        value = '{}'.format(s.regals_run_settings['tol'])
                    else:
                        value = ''
                else:
                    value = ''

                if header in regals_table_dict:
                    regals_table_dict[header].append(value)
                else:
                    regals_table_dict[header] = [value]

            for header, values in regals_table_dict.items():
                if header in regals_comp_required_data:
                    regals_table_entry = [header]
                    regals_table_entry.extend(values)
                    regals_table_data.append(regals_table_entry)
                else:
                    if any(val != '' for val in values):
                        regals_table_entry = [header]
                        regals_table_entry.extend(values)
                        regals_table_data.append(regals_table_entry)

            regals_table = Table(regals_table_data)

            table_style = TableStyle(
                [('LINEBELOW', (0, 0), (-1, 0), 1, colors.black),
                ('LINEAFTER', (0, 0), (0,-1), 1, colors.black),
                ])

            regals_table.setStyle(table_style)
            regals_table.hAlign = 'LEFT'
            # regals_table = KeepTogether([regals_table])


            #Make REGALS compoment table
            for k, comp_settings in enumerate(s.regals_component_settings):
                prof_comp = comp_settings[0]
                conc_comp = comp_settings[1]

                for header, key in regals_profile_component_table_pairs:
                    if key == 'name':
                        value = '{}'.format(k)
                    elif key == 'lambda':
                        value = '{}'.format(text_round(prof_comp[key], 2))
                    elif key in prof_comp:
                        value = '{}'.format(prof_comp[key])
                    elif key in prof_comp['kwargs']:
                        value = '{}'.format(prof_comp['kwargs'][key])
                    else:
                        value = ''

                    if header in regals_comp_table_dict:
                        regals_comp_table_dict[header].append(value)
                    else:
                        regals_comp_table_dict[header] = [value]

                for header, key in regals_conc_component_table_pairs:
                    if key == 'lambda':
                        value = '{}'.format(text_round(conc_comp[key], 2))
                    elif key in conc_comp:
                        value = '{}'.format(conc_comp[key])
                    elif key in conc_comp['kwargs']:
                        value = '{}'.format(conc_comp['kwargs'][key])
                    else:
                        value = ''

                    if header in regals_comp_table_dict:
                        regals_comp_table_dict[header].append(value)
                    else:
                        regals_comp_table_dict[header] = [value]

            for header, values in regals_comp_table_dict.items():
                if header in regals_required_data:
                    regals_table_entry = [header]
                    regals_table_entry.extend(values)
                    regals_comp_table_data.append(regals_table_entry)
                else:
                    if any(val != '' for val in values):
                        regals_table_entry = [header]
                        regals_table_entry.extend(values)
                        regals_comp_table_data.append(regals_table_entry)

            regals_comp_table = Table(regals_comp_table_data)

            table_style = TableStyle(
                [('LINEBELOW', (0, 0), (-1, 0), 1, colors.black),
                ('LINEAFTER', (0, 0), (0,-1), 1, colors.black),
                ])

            regals_comp_table.setStyle(table_style)
            regals_comp_table.hAlign = 'LEFT'

            # Make REGALS plot

            regals_plot_panel = efa_plot(s, is_regals=True)

            img_width = 6
            img_height = 6

            regals_caption = ('REGALS deconvolution results. a) The full series '
                'intensity (blue), the selected intensity range for regals '
                '(black), and (if available) Rg values (red). b) The selected '
                'intensity range for regals (black), and the individual component '
                'ranges for deconvolution, with component range 0 starting at '
                'the top left, and component number increasing in descending '
                'order to the right. c) Mean chi^2 values between the '
                'fit of the regals deconvolution and the original data. d) '
                'Area normalized concentration profiles for each component. '
                'Colors match the component range colors in b. e) Deconvolved '
                'scattering profiles. Colors match the component range colors in '
                'b.')

            if len(s.regals_ifts) > 0:
                regals_caption = regals_caption + (' f) P(r) functions. Colors '
                    'match the component range colors in b.')

            regals_figure = make_figure(regals_plot_panel.figure, regals_caption,
                img_width, img_height, styles)

            regals_spacer = Paragraph('<br/><br/>', styles['Normal'])
            regals_elements.append(regals_title)
            regals_elements.append(regals_table)
            regals_elements.append(regals_spacer)
            regals_elements.append(regals_comp_table)
            regals_elements.append(regals_figure)

            has_series_data = True

    elements.extend(regals_elements)

    if not has_series_data:
        elements = []

    return elements

def generate_guinier_params(profiles, ifts, series):
    styles = getSampleStyleSheet()

    guinier_text = Paragraph('Guinier:', styles['Heading2'])

    absolute = [profile.metadata.Absolute_scale for profile in profiles]
    if all(absolute):
        i0_label = 'I(0) [1/cm]'
    else:
        i0_label = 'I(0) [Arb.]'

    table_pairs = [
        ('', 'name'),
        ('Rg [A]', 'rg'),
        (i0_label, 'i0'),
        ('q-range [1/A]', 'q_range'),
        ('qmin*Rg', 'qRg_min'),
        ('qmax*Rg', 'qRg_max'),
        ('r^2', 'rsq'),
        ]

    table_dict = OrderedDict()

    table_data = []

    required_data = ['']


    for profile in profiles:
        filename = profile.filename

        if profile.guinier_data.Rg != -1:
            rg = profile.guinier_data.Rg
            rg_err = profile.guinier_data.Rg_err
            i0 = profile.guinier_data.I0
            i0_err = profile.guinier_data.I0_err
            qmin = profile.guinier_data.q_min
            qmax = profile.guinier_data.q_max
            qRg_min = profile.guinier_data.qRg_min
            qRg_max = profile.guinier_data.qRg_max
            rsq = profile.guinier_data.r_sq

            rg = '{} +/- {}'.format(text_round(rg, 2),
                text_round(rg_err, 2))
            i0 = '{} +/- {}'.format(text_round(i0, 2),
                text_round(i0_err, 2))

            q_range = '{} to {}'.format(text_round(qmin, 4), text_round(qmax, 4))
            qmin_Rg = '{}'.format(text_round(qRg_min, 3))
            qmax_Rg = '{}'.format(text_round(qRg_max, 3))
            rsq = '{}'.format(text_round(rsq, 3))
        else:
            rg = ''
            i0 = ''
            q_range = ''
            qmin_Rg = ''
            qmax_Rg = ''
            rsq = ''


        for header, key in table_pairs:
            if key == 'name':
                value = filename
            elif key == 'rg':
                value = rg
            elif key == 'i0':
                value = i0
            elif key == 'q_range':
                value = q_range
            elif key == 'qRg_min':
                value = qmin_Rg
            elif key == 'qRg_max':
                value = qmax_Rg
            elif key == 'rsq':
                value = rsq

            if header in table_dict:
                table_dict[header].append(value)
            else:
                table_dict[header] = [value]

    for header, values in table_dict.items():
        if header in required_data:
            table_entry = [header]
            table_entry.extend(values)
            table_data.append(table_entry)
        else:
            if any(val != '' for val in values):
                table_entry = [header]
                table_entry.extend(values)
                table_data.append(table_entry)

    if len(table_data) > 1:
        guinier_table = Table(table_data)

        table_style = TableStyle(
            [('LINEBELOW', (0, 0), (-1, 0), 1, colors.black),
            ('LINEAFTER', (0, 0), (0,-1), 1, colors.black),
            ])

        guinier_table.setStyle(table_style)
        guinier_table.hAlign = 'LEFT'
        guinier_table = KeepTogether([guinier_text, guinier_table])

        elements = [guinier_table]

    else:
        elements = []

    return elements

def generate_mw_params(profiles, ifts, series):
    styles = getSampleStyleSheet()

    mw_text = Paragraph('Molecular weight:', styles['Heading2'])

    table_pairs = [
        ('', 'name'),
        ('M.W. (Vp) [kDa]', 'mw_vp'),
        ('Porod Volume [A^3]', 'vp'),
        ('M.W. (Vc) [kDa]', 'mw_vc'),
        ('M.W. (S&S) [kDa]', 'mw_ss'),
        ('Shape (S&S)', 'shape'),
        ('Dmax (S&S)', 'dmax'),
        ('M.W. (Bayes) [kDa]', 'mw_bayes'),
        ('Bayes Probability', 'prob'),
        ('Bayes Confidence\nInterval [kDa]', 'ci'),
        ('Bayes C.I. Prob.', 'ci_prob'),
        ('M.W. (Abs.) [kDa]', 'mw_abs'),
        ('M.W. (Ref.) [kDa]', 'mw_ref'),
        ]

    table_dict = OrderedDict()

    table_data = []

    required_data = ['']


    for profile in profiles:

        filename = profile.filename

        mw_data = defaultdict(str)

        for mw_key in profile.mw_data:
            mw_val = profile.mw_data[mw_key]

            if mw_val.MW != -1 and mw_val.MW != '':
                val = text_round(mw_val.MW, 1)
            else:
                val = ''

            if mw_key == 'Volume_of_correlation':
                table_key = 'mw_vc'
            elif mw_key == 'Porod_volume':
                table_key = 'mw_vp'
            elif mw_key == 'Shape_and_size':
                table_key = 'mw_ss'
            elif mw_key == 'Bayesian':
                table_key = 'mw_bayes'
            elif mw_key == 'Absolute':
                table_key = 'mw_abs'
            elif mw_key == 'Reference':
                table_key = 'mw_ref'

            mw_data[table_key] = val

        p_vol = ''

        if 'Porod_volume' in profile.mw_data:
            mw_val = profile.mw_data['Porod_volume']

            if mw_val.MW != -1 and mw_val.MW != '':
                p_vol = '{}'.format(text_round(mw_val.Porod_volume_corrected, 2))

        prob = ''
        ci = ''
        ci_prob = ''

        if 'Bayesian' in profile.mw_data:
            mw_val = profile.mw_data['Bayesian']

            if mw_val.MW != -1 and mw_val.MW != '':
                prob = '{}'.format(text_round(mw_val.Probability, 1))

                ci_lower = mw_val.Confidence_interval_lower
                ci_upper = mw_val.Confidence_interval_upper

                ci = ('{} to {}'.format(text_round(ci_lower, 1),
                    text_round(ci_upper, 1)))

                ci_prob = '{}'.format(text_round(mw_val.Confidence_interval_probability, 1))

        shape = ''
        dmax = ''

        if 'Shape_and_size' in profile.mw_data:
            mw_val = profile.mw_data['Shape_and_size']

            if mw_val.MW != -1 and mw_val.MW != '':
                shape = mw_val.Shape
                dmax = '{}'.format(text_round(mw_val.Dmax, 1))


        for header, key in table_pairs:
            if key == 'name':
                value = filename
            elif key.startswith('mw'):
                value =mw_data[key]
            elif key == 'vp':
                value = p_vol
            elif key == 'prob':
                value = prob
            elif key == 'ci':
                value = ci
            elif key == 'ci_prob':
                value = ci_prob
            elif key == 'shape':
                value = shape
            elif key == 'dmax':
                value = dmax

            if header in table_dict:
                table_dict[header].append(value)
            else:
                table_dict[header] = [value]

    for header, values in table_dict.items():
        if header in required_data:
            table_entry = [header]
            table_entry.extend(values)
            table_data.append(table_entry)
        else:
            if any(val != '' for val in values):
                table_entry = [header]
                table_entry.extend(values)
                table_data.append(table_entry)

    if len(table_data) > 1:
        mw_table = Table(table_data)

        table_style = TableStyle(
            [('LINEBELOW', (0, 0), (-1, 0), 1, colors.black),
            ('LINEAFTER', (0, 0), (0,-1), 1, colors.black),
            ])

        mw_table.setStyle(table_style)
        mw_table.hAlign = 'LEFT'
        mw_table = KeepTogether([mw_text, mw_table])

        elements = [mw_table]

    else:
        elements = []

    return elements

def generate_gnom_params(profiles, ifts, series):
    styles = getSampleStyleSheet()

    gnom_text = Paragraph('GNOM IFT:', styles['Heading2'])

    table_pairs = [
        ('', 'name'),
        ('Dmax [A]', 'dmax'),
        ('Rg [A]', 'rg'),
        ('I(0)', 'i0'),
        ('Chi^2', 'chi_sq'),
        ('Total Estimate', 'te'),
        ('Quality', 'quality'),
        ('q-range [1/A]', 'q_range'),
        ('Ambiguity score', 'a_score'),
        ('Ambiguity cats.', 'a_cats'),
        ('Ambiguity', 'a_interp'),
        ]

    table_dict = OrderedDict()

    table_data = []

    required_data = ['']


    for ift in ifts:

        if ift.type == 'GNOM':
            filename = ift.filename

            dmax = ift.dmax
            rg = ift.rg
            rg_err = ift.rg_err
            i0 = ift.i0
            i0_err = ift.i0_err
            chi_sq = ift.chi_sq
            te = ift.total_estimate
            quality = ift.quality

            dmax = '{}'.format(text_round(dmax, 0))

            rg = '{} +/- {}'.format(text_round(rg, 2),
                text_round(rg_err, 2))
            i0 = '{} +/- {}'.format(text_round(i0, 2),
                text_round(i0_err, 2))

            chi_sq = '{}'.format(text_round(chi_sq, 3))
            te = '{}'.format(text_round(te, 3))

            q_range = '{} to {}'.format(text_round(ift.q[0], 4),
                text_round(ift.q[-1], 4))

            if ift.a_score != -1:
                a_score = '{}'.format(text_round(ift.a_score, 2))
                a_cats = '{}'.format(ift.a_cats, 0)
                a_interp = ift.a_interp
            else:
                a_score = ''
                a_cats = ''
                a_interp = ''

            for header, key in table_pairs:
                if key == 'name':
                    value = filename
                elif key == 'dmax':
                    value = dmax
                elif key == 'rg':
                    value = rg
                elif key == 'i0':
                    value = i0
                elif key == 'chi_sq':
                    value = chi_sq
                elif key == 'te':
                    value = te
                elif key == 'quality':
                    value = quality
                elif key == 'q_range':
                    value = q_range
                elif key == 'a_score':
                    value = a_score
                elif key == 'a_cats':
                    value = a_cats
                elif key == 'a_interp':
                    value = a_interp

                if header in table_dict:
                    table_dict[header].append(value)
                else:
                    table_dict[header] = [value]

    for header, values in table_dict.items():
        if header in required_data:
            table_entry = [header]
            table_entry.extend(values)
            table_data.append(table_entry)
        else:
            if any(val != '' for val in values):
                table_entry = [header]
                table_entry.extend(values)
                table_data.append(table_entry)

    if len(table_data) > 1:
        gnom_table = Table(table_data)

        table_style = TableStyle(
            [('LINEBELOW', (0, 0), (-1, 0), 1, colors.black),
            ('LINEAFTER', (0, 0), (0,-1), 1, colors.black),
            ])

        gnom_table.setStyle(table_style)
        gnom_table.hAlign = 'LEFT'
        gnom_table = KeepTogether([gnom_text, gnom_table])

        elements = [gnom_table]

    else:
        elements = []

    return elements

def generate_bift_params(profiles, ifts, series):
    styles = getSampleStyleSheet()

    bift_text = Paragraph('BIFT:', styles['Heading2'])

    table_pairs = [
        ('', 'name'),
        ('Dmax [A]', 'dmax'),
        ('Rg [A]', 'rg'),
        ('I(0)', 'i0'),
        ('Chi^2', 'chi_sq'),
        ('q-range [1/A]', 'q_range'),
        ]

    table_dict = OrderedDict()

    table_data = []

    required_data = ['']


    for ift in ifts:

        if ift.type == 'BIFT':
            filename = ift.filename

            dmax = ift.dmax
            dmax_err = ift.dmax_err
            rg = ift.rg
            rg_err = ift.rg_err
            i0 = ift.i0
            i0_err = ift.i0_err
            chi_sq = ift.chi_sq

            dmax = '{} +/- {}'.format(text_round(dmax, 1),
                text_round(dmax_err, 1))

            rg = '{} +/- {}'.format(text_round(rg, 2),
                text_round(rg_err, 2))
            i0 = '{} +/- {}'.format(text_round(i0, 2),
                text_round(i0_err, 2))

            chi_sq = '{}'.format(text_round(chi_sq, 3))

            q_range = '{} to {}'.format(text_round(ift.q[0], 4),
                text_round(ift.q[-1], 4))

            for header, key in table_pairs:
                if key == 'name':
                    value = filename
                elif key == 'dmax':
                    value = dmax
                elif key == 'rg':
                    value = rg
                elif key == 'i0':
                    value = i0
                elif key == 'chi_sq':
                    value = chi_sq
                elif key == 'q_range':
                    value = q_range

                if header in table_dict:
                    table_dict[header].append(value)
                else:
                    table_dict[header] = [value]

    for header, values in table_dict.items():
        if header in required_data:
            table_entry = [header]
            table_entry.extend(values)
            table_data.append(table_entry)
        else:
            if any(val != '' for val in values):
                table_entry = [header]
                table_entry.extend(values)
                table_data.append(table_entry)

    if len(table_data) > 1:
        bift_table = Table(table_data)

        table_style = TableStyle(
            [('LINEBELOW', (0, 0), (-1, 0), 1, colors.black),
            ('LINEAFTER', (0, 0), (0,-1), 1, colors.black),
            ])

        bift_table.setStyle(table_style)
        bift_table.hAlign = 'LEFT'
        bift_table = KeepTogether([bift_text, bift_table])

        elements = [bift_table]

    else:
        elements = []

    return elements

def generate_dammif_params(dammif_data):
    styles = getSampleStyleSheet()

    dammif_text = Paragraph('Bead model reconstructions:', styles['Heading2'])

    table_pairs = [
        ('', 'prefix'),
        ('Program', 'program'),
        ('Mode', 'mode'),
        ('Symmetry', 'sym'),
        ('Anisometry', 'aniso'),
        ('Number of reconstructions', 'num'),
        ('Ran DAMAVER', 'damaver'),
        ('Ran DAMCLUST', 'damclust'),
        ('Refined with DAMMIN', 'refined'),
        ('Mean NSD', 'nsd'),
        ('Included models', 'included'),
        ('Resolution (SASRES)', 'res'),
        ('Representative model', 'rep_model'),
        ('Number of clusters', 'clusters'),
        ]

    table_dict = OrderedDict()

    table_data = []

    required_data = ['']


    for info in dammif_data:
        if info is not None:
            for header, key in table_pairs:
                value = getattr(info, key)

                if value == -1:
                    value = ''

                else:
                    if key == 'nsd':
                        value = '{} +/- {}'.format(text_round(value, 3),
                            text_round(info.nsd_std, 3))
                    elif key == 'res':
                        value = '{} +/- {}'.format(text_round(value, 0),
                            text_round(info.res_err, 0))
                    elif key == 'included':
                        value = '{} of {}'.format(value, info.num)
                    elif not isinstance(value, str):
                        value = '{}'.format(value)

                if header in table_dict:
                    table_dict[header].append(value)
                else:
                    table_dict[header] = [value]

    for header, values in table_dict.items():
        if header in required_data:
            table_entry = [header]
            table_entry.extend(values)
            table_data.append(table_entry)
        else:
            if any(val != '' for val in values):
                table_entry = [header]
                table_entry.extend(values)
                table_data.append(table_entry)

    if len(table_data) > 1:
        dammif_table = Table(table_data)

        table_style = TableStyle(
            [('LINEBELOW', (0, 0), (-1, 0), 1, colors.black),
            ('LINEAFTER', (0, 0), (0,-1), 1, colors.black),
            ])

        dammif_table.setStyle(table_style)
        dammif_table.hAlign = 'LEFT'
        dammif_table = KeepTogether([dammif_text, dammif_table])

        elements = [dammif_table]

    else:
        elements = []

    return elements

def generate_denss_params(denss_data):
    styles = getSampleStyleSheet()

    denss_text = Paragraph('Electron density reconstructions:', styles['Heading2'])

    table_pairs = [
        ('', 'prefix'),
        ('Mode', 'mode'),
        ('Symmetry', 'sym'),
        ('Symmetry axis', 'sym_axis'),
        ('Symmetry factpr', 'sym_factor'),
        ('Number of reconstructions', 'num'),
        ('Averaged', 'average'),
        ('Refined', 'refined'),
        ('Mean RSC', 'rsc'),
        ('Included models', 'included'),
        ('Resolution', 'res'),
        ]

    table_dict = OrderedDict()

    table_data = []

    required_data = ['']


    for info in denss_data:
        if info is not None:
            for header, key in table_pairs:
                value = getattr(info, key)

                if value == -1:
                    value = ''

                else:
                    if key == 'rsc':
                        value = '{} +/- {}'.format(text_round(value, 3),
                            text_round(info.rsc_std, 3))
                    elif key == 'res':
                        value = '{}'.format(text_round(value, 1))
                    elif key == 'included':
                        value = '{} of {}'.format(value, info.num)
                    elif not isinstance(value, str):
                        value = '{}'.format(value)

                if header in table_dict:
                    table_dict[header].append(value)
                else:
                    table_dict[header] = [value]

    for header, values in table_dict.items():
        if header in required_data:
            table_entry = [header]
            table_entry.extend(values)
            table_data.append(table_entry)
        else:
            if any(val != '' for val in values):
                table_entry = [header]
                table_entry.extend(values)
                table_data.append(table_entry)

    if len(table_data) > 1:
        denss_table = Table(table_data)

        table_style = TableStyle(
            [('LINEBELOW', (0, 0), (-1, 0), 1, colors.black),
            ('LINEAFTER', (0, 0), (0,-1), 1, colors.black),
            ])

        denss_table.setStyle(table_style)
        denss_table.hAlign = 'LEFT'
        denss_table = KeepTogether([denss_text, denss_table])

        elements = [denss_table]

    else:
        elements = []

    return elements

def make_figure(figure, caption, img_width, img_height, styles):
    """
    Takes as input matplotlib figure, a string, and image width and height in
    inches and returns a flowable with image and caption thnat will stay together.
    """
    datadir = os.path.abspath(os.path.expanduser(tempfile.gettempdir()))
    filename = tempfile.NamedTemporaryFile(dir=datadir).name
    filename = os.path.split(filename)[-1] + '.png'
    filename = os.path.join(datadir, filename)

    global temp_files
    temp_files.append(filename) #Note defined at a module level

    figure.savefig(filename, dpi=300)
    plt.close(figure)
    image = Image(filename, img_width*inch, img_height*inch, lazy=2)
    image.hAlign = 'CENTER'

    text = Paragraph(caption, styles['Normal'])

    return_fig = KeepTogether([image, text])

    return return_fig

def make_report_from_files(name, out_dir, data_dir, profiles, ifts, series,
    dammif_data=None, denss_data=None):

    data_dir = os.path.abspath(os.path.expanduser(data_dir))

    for j in range(len(profiles)):
        profiles[j] = os.path.abspath(os.path.expanduser(os.path.join(data_dir,
            profiles[j])))

    for j in range(len(ifts)):
        ifts[j] = os.path.abspath(os.path.expanduser(os.path.join(data_dir,
            ifts[j])))

    for j in range(len(series)):
        series[j] = os.path.abspath(os.path.expanduser(os.path.join(data_dir,
            series[j])))

    if dammif_data is not None:
        for j in range(len(dammif_data)):
            if dammif_data[j] is not None:
                dammif_data[j] = os.path.abspath(os.path.expanduser(
                    os.path.join(data_dir, dammif_data[j])))


    profiles = raw.load_profiles(profiles)
    profile_data = [data.SAXSData(profile) for profile in profiles]

    ifts = raw.load_ifts(ifts)
    ift_data = [data.IFTData(ift) for ift in ifts]

    for j, ift in enumerate(ift_data):
        if ift.type == 'GNOM':
            try:
                a_score, a_cats, a_interp = raw.ambimeter(ifts[j])

                ift.a_score = a_score
                ift.a_cats = a_cats
                ift.a_interp = a_interp
            except Exception:
                pass

    series = raw.load_series(series)
    series_data = [data.SECData(s) for s in series]

    dammif_results = []

    if dammif_data is not None:
        for dammif_file in dammif_data:
            if dammif_file is not None:
                results = data.parse_dammif_file(dammif_file)
            else:
                results = None

            dammif_results.append(results)

    denss_results = []

    if denss_data is not None:
        for denss_file in denss_data:
            if denss_file is not None:
                results = data.parse_denss_file(denss_file)
            else:
                results = None

            denss_results.append(results)

    extra_data = {'dammif': dammif_results, 'denss': denss_results}

    generate_report(name, out_dir, profile_data, ift_data, series_data,
        extra_data)


def make_report_from_data(name, out_dir, profiles, ifts, series,
    dammif_data=None, denss_data=None):

    profile_data = [data.SAXSData(profile) for profile in profiles]
    ift_data = [data.IFTData(ift) for ift in ifts]
    series_data = [data.SECData(s) for s in series]

    for j, ift in enumerate(ift_data):
        if ift.type == 'GNOM':
            try:
                a_score, a_cats, a_interp = raw.ambimeter(ifts[j])

                ift.a_score = a_score
                ift.a_cats = a_cats
                ift.a_interp = a_interp
            except Exception:
                pass

    if dammif_data is not None:
        extra_data = {'dammif': dammif_data}
    else:
        extra_data = {'dammif': []}

    if denss_data is not None:
        extra_data['denss'] = denss_data
    else:
        extra_data['denss'] = []

    generate_report(name, out_dir, profile_data, ift_data, series_data,
        extra_data)
