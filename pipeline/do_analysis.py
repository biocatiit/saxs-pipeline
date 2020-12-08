import bioxtasraw.RAWAPI as raw

from .analysis import analyze
from .reports import pdf


def analyze_and_report_from_processed_files(name, out_dir, data_dir, profiles,
    ifts, raw_settings=None):

    if raw_settings is None:
        raw_settings = raw.__default_settings

    profiles, ifts, all_dammif_data, all_denss_data = analyze.analyze_files(out_dir,
        data_dir, profiles, ifts, raw_settings)

    series = []

    pdf.make_report_from_data(name, out_dir, profiles, ifts, series,
        dammif_data=all_dammif_data, denss_data=all_denss_data)

    for profile in profiles:
        raw.save_profile(profile, datadir=out_dir, settings=raw_settings)

    for ift in ifts:
        raw.save_ift(ift, datadir=out_dir)

def analyze_and_report_from_data(name, out_dir, profiles, ifts, raw_settings=None):
    if raw_settings is None:
        raw_settings = raw.__default_settings

    profiles, ifts, all_dammif_data, all_denss_data = analyze.analyze_data(out_dir,
        profiles, ifts, raw_settings)

    series = []

    pdf.make_report_from_data(name, out_dir, profiles, ifts, series,
        dammif_data=all_dammif_data, denss_data=all_denss_data)

    for profile in profiles:
        raw.save_profile(profile, datadir=out_dir, settings=raw_settings)

    for ift in ifts:
        raw.save_ift(ift, datadir=out_dir)



"""
Example use:
import bioxtasraw.RAWAPI as raw

pipeline_path = os.path.abspath(os.path.expanduser('~/Documents/software_dev/saxs-pipeline'))
if pipeline_path not in os.sys.path:
    os.sys.path.append(pipeline_path)

from pipeline import do_analysis
import pipeline.analysis

#From files:
do_analysis.analyze_and_report_from_processed_files('test_analysis.pdf', '.',
    '../data', ['GI_sub_clean.dat'], [])


#From data in memory:
profiles = raw.load_profiles(['../data/GI_sub_clean.dat'])

do_analysis.analyze_and_report_from_data('test_analysis_mem.pdf', '.', profiles,
    [])

"""
