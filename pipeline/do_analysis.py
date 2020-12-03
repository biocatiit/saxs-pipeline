import bioxtasraw.RAWAPI as raw

from .analysis import analyze
from .reports import pdf


def analyze_and_report_from_processed_files(name, out_dir, data_dir, profiles, ifts, raw_settings=None):

    if raw_settings is None:
        raw_settings = raw.__default_settings

    print('Doing analysis')
    profiles, ifts = analyze.analyze_files(out_dir, data_dir, profiles, ifts, raw_settings)

    series = []

    print('Generating report')
    pdf.make_report_from_data(name, out_dir, profiles, ifts, series, dammif_data=None)

    for profile in profiles:
        raw.save_profile(profile, datadir=out_dir, settings=raw_settings)

    for ift in ifts:
        raw.save_ift(ift, datadir=out_dir)
