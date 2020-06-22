from collections import OrderedDict, defaultdict

from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, Image, XPreformatted, KeepTogether, TableStyle, CondPageBreak

from . import plots

from .utils import text_round as text_round
from . import utils

def generate_report(profiles, ifts, series, extra_data=None):
    """
    Inputs a list of profile data, a list of ift data, and a list of series
    data to be included in the report. Makes a PDF report.
    """

    elements = []

    overview = generate_overview(profiles, ifts, series)
    elements.extend(overview)
    elements.append(CondPageBreak(2*inch))

    exp = generate_exp_params(profiles, ifts, series)
    elements.extend(exp)
    elements.append(CondPageBreak(2*inch))

    s_elements = generate_series_params(profiles, ifts, series, extra_data)
    elements.extend(s_elements)
    elements.append(CondPageBreak(2*inch))

    # guinier = generate_guinier_params(profiles, ifts, series)
    # elements.extend(guinier)
    # elements.append(CondPageBreak(2*inch))

    # mw = generate_mw_params(profiles, ifts, series)
    # elements.extend(mw)
    # elements.append(CondPageBreak(2*inch))

    # gnom = generate_gnom_params(profiles, ifts, series)
    # elements.extend(gnom)
    # elements.append(CondPageBreak(2*inch))

    # bift = generate_bift_params(profiles, ifts, series)
    # elements.extend(bift)
    # elements.append(CondPageBreak(2*inch))

    # dammif = generate_dammif_params(profiles, ifts, series)
    # elements.extend(dammif)
    # elements.append(CondPageBreak(2*inch))

    doc = SimpleDocTemplate('test.pdf', pagesize=letter, leftMargin=1*inch,
        rightMargin=1*inch, topMargin=1*inch, bottomMargin=1*inch)
    doc.build(elements)

def generate_overview(profiles, ifts, series):
    """
    Generates the overview portion of the report. Returns a list of flowables.
    """
    styles = getSampleStyleSheet()

    name_list = []
    date_list = []

    for s in series:
        if 'prefix' in s.metadata:
            name_list.append(s.metadata['prefix'])
        else:
            name_list.append(s.filename)

        if 'date' in s.metadata:
            date_list.append(':'.join(s.metadata['date'].split(':')[:-1]))
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

    summary_text = ('Data name(s): {}\n'
        'Collection date(s): {}'.format(name_str, date_str))
    ov_summary = XPreformatted(summary_text, styles['Normal'])


    # table_data = []

    # table_header = ['', 'Guinier Rg', 'Guinier I(0)']

    # mw_keys = []

    # for profile in profiles:
    #     for key, value in profile.mw_data.items():
    #         if value.MW != -1 and value.MW != '':
    #             if key == 'Volume_of_correlation':
    #                 header = 'M.W. (Vc)'
    #             elif key == 'Porod_volume':
    #                 header = 'M.W. (Vp)'
    #             elif key == 'Shape_and_size':
    #                 header = 'M.W. (S&S)'
    #             elif key == 'Bayesian':
    #                 header = 'M.W. (Bayes)'
    #             else:
    #                 header = key

    #             if key not in mw_keys:
    #                 mw_keys.append(key)

    #             if header not in table_header:
    #                 table_header.append(header)

    # table_header.extend(['GNOM Dmax', 'GNOM Rg', 'GNOM I(0)'])

    table_pairs = [
        ('', 'name'),
        ('Guinier Rg', 'gu_rg'),
        ('Guinier I(0)', 'gu_i0'),
        ('M.W. (Vp)', 'mw_vp'),
        ('M.W. (Vc)', 'mw_vc'),
        ('M.W. (S&S)', 'mw_ss'),
        ('M.W. (Bayes)', 'mw_bayes'),
        ('M.W. Abs.', 'mw_abs'),
        ('M.W. Ref.', 'mw_ref'),
        ('GNOM Dmax', 'gn_dmax'),
        ('GNOM Rg', 'gn_rg'),
        ('GNOM I(0)', 'gn_i0'),
        ('BIFT Dmax', 'b_dmax'),
        ('BIFT Rg', 'b_rg'),
        ('BIFT I(0)', 'b_i0'),
        ]

    table_dict = OrderedDict()

    table_data = []

    required_data = ['', 'Guinier Rg', 'Guinier I(0)', 'M.W. (Vp)', 'M.W. (Vc)']


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

    ov_table = Table(table_data, spaceBefore=0.25*inch, spaceAfter=0.1*inch)

    table_style = TableStyle(
        [('INNERGRID', (0, 0), (-1, -1), 1, colors.black),
        ('ALIGN', (1,0), (-1, -1), 'CENTER')
        ])

    ov_table.setStyle(table_style)
    ov_table.hAlign = 'LEFT'

    if 'N/A' in name_str:
        table_caption = ('SAXS data summary table. ')
    else:
        table_caption = ('SAXS data summary table for {}. '.format(name_str))

    table_text = Paragraph(table_caption, styles['Normal'])

    ov_table = KeepTogether([ov_table, table_text])

    img_width = 6
    img_height = 6

    ov_plot = plots.series_overview_plot(profiles, ifts, series,
        img_width=img_width, img_height=img_height)

    if 'N/A' in name_str:
        caption = ('SAXS data summary figure. ')
    else:
        caption = ('SAXS data summary figure for {}. '.format(name_str))

    caption = caption + ('a) Series intensity (left axis) vs. '
        'frame, and if available, Rg vs. frame (right axis). Green '
        'shaded regions are buffer regions, purple shaded regions are sample '
        'regions. b) Scattering profile(s) on a log-lin scale. c) Guinier fit(s) '
        '(top) and fit residuals (bottom). d) Normalized Kratky plot. Dashed '
        'lines show where a globuar system would peak. e) P(r) function(s), '
        'normalized by I(0).')

    ov_figure = make_figure(ov_plot.figure, caption, img_width, img_height,
        styles)

    return [ov_title, ov_text, ov_summary, ov_figure, ov_table]


def generate_exp_params(profiles, ifts, series):
    styles = getSampleStyleSheet()

    name_list = []
    date_list = []

    for s in series:
        if 'prefix' in s.metadata:
            name_list.append(s.metadata['prefix'])
        else:
            name_list.append(s.filename)

        if 'date' in s.metadata:
            date_list.append(':'.join(s.metadata['date'].split(':')[:-1]))
        else:
            date_list.append('')

    exp_text = Paragraph('Experimental parameters:', styles['Heading2'])


    table_pairs = [
        ('', 'name'),
        ('Date', 'date'),
        ('Instrument', 'instrument'),
        ('Detector', 'detector'),
        ('Wavelength (A)', 'wavelength'),
        ('Camera length (m)', 'distance'),
        ('q-measurement range (1/A)', 'q_range'),
        ('Exposure time (s)', 'exp_time'),
        ('Exposure period (s)', 'exp_period'),
        ('Flow rate (ml/min)', 'flow_rate'),
        ('RAW version', 'version'),
        ]

    table_dict = OrderedDict()

    table_data = []

    required_data = ['', 'Date', 'Wavelength (A)', 'Camera length (m)',
        'Exposure time (s)']

    for j, s in enumerate(series):

        for header, key in table_pairs:
            if key in s.metadata:
                value = s.metadata[key]
            else:
                value = ''

            if key == 'wavelength':
                if value != '' and value != -1:
                    value = text_round(value, 3)
                else:
                    value = ''

            elif key == 'distance':
                if value != '' and value != -1:
                    value = float(value)/1000.
                    value = text_round(value, 3)
                else:
                    value = ''

            elif key == 'name':
                value = name_list[j]

            elif key == 'date':
                value = date_list[j]

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

    exp_table = Table(table_data, spaceBefore=0.25*inch, spaceAfter=0.1*inch)

    table_style = TableStyle(
        [('LINEBELOW', (0, 0), (-1, 0), 1, colors.black),
        ('LINEAFTER', (0, 0), (0,-1), 1, colors.black),
        ])

    exp_table.setStyle(table_style)
    exp_table.hAlign = 'LEFT'

    return [exp_text, exp_table]


def generate_series_params(profiles, ifts, series, extra_data):
    styles = getSampleStyleSheet()

    name_list = []

    for s in series:
        if 'prefix' in s.metadata:
            name_list.append(s.metadata['prefix'])
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

    required_data = ['', 'Buffer range', 'Sample range', 'Baseline correction',]


    for j, s in enumerate(series):
        buffer_range = ', '.join(['{} to {}'.format(*br) for br in s.buffer_range])
        sample_range = ', '.join(['{} to {}'.format(*sr) for sr in s.sample_range])

        if s.baseline_corrected:
            baseline_type = s.baseline_type
            baseline_start = '{} to {}'.format(*s.baseline_start_range)
            baseline_end = '{} to {}'.format(*s.baseline_end_range)
        else:
            baseline_type = 'None'
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

    series_table = Table(table_data, spaceBefore=0.25*inch, spaceAfter=0.1*inch)

    table_style = TableStyle(
        [('LINEBELOW', (0, 0), (-1, 0), 1, colors.black),
        ('LINEAFTER', (0, 0), (0,-1), 1, colors.black),
        ])

    series_table.setStyle(table_style)
    series_table.hAlign = 'LEFT'
    series_table = KeepTogether([series_text, series_table])

    efa_elements = []

    for j, s in enumerate(series):
        if s.efa_done:
            if len(series) > 1:
                efa_title = Paragraph('{} EFA results:'.format(name_list[j]),
                    styles['Heading3'])
            else:
                efa_title = Paragraph('EFA results:', styles['Heading3'])

            if extra_data is not None:
                extra_efa_data = extra_data['efa'][j]
            else:
                extra_efa_data = None

            efa_plot = plots.efa_plot(s, extra_efa_data)

    return [series_table]

def make_figure(figure, caption, img_width, img_height, styles):
    """
    Takes as input matplotlib figure, a string, and image width and height in
    inches and returns a flowable with image and caption thnat will stay together.
    """

    figure.savefig('fig.png', dpi=300)
    image = Image('fig.png', img_width*inch, img_height*inch)
    image.hAlign = 'CENTER'

    text = Paragraph(caption, styles['Normal'])

    return_fig = KeepTogether([image, text])

    return return_fig


