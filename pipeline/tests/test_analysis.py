import os
import multiprocessing as mp

pipeline_path = os.path.abspath(os.path.expanduser('~/Documents/software_dev/saxs-pipeline/'))
if pipeline_path not in os.sys.path:
    os.sys.path.append(pipeline_path)

from pipeline.analysis import analyze

import bioxtasraw.RAWAPI as raw


if __name__ == '__main__':
    analysis_args = {
        'dammif'        : True,
        'denss'         : True,
        'dam_runs'      : 3,
        'dam_mode'      : 'Fast',
        'damaver'       : True,
        'damclust'      : False,
        'dam_refine'    : False,
        'denss_runs'    : 4,
        'denss_mode'    : 'Fast',
        'denss_aver'    : True,
        'denss_refine'  : False,
        'use_atsas'     : True,
        }

    cmd_q = mp.Queue()
    ret_q = mp.Queue()
    cmd_lock = mp.Lock()
    ret_lock = mp.Lock()
    abort_event = mp.Event()
    log_lock = mp.Lock()
    log_q = mp.Queue()

    raw_settings_file = './20220316_SAXS.cfg'

    print('Starting analysis process')

    my_proc = analyze.analysis_process(cmd_q, ret_q, cmd_lock, ret_lock, abort_event,
            raw_settings_file, log_lock, log_q)

    my_proc.start()

    series_name = ['./TM01_011.hdf5']

    print('Loading series')

    series = raw.load_series(series_name)[0]
    profiles = series.getAllSASMs()

    print('Sending command')

    cmd_q.put_nowait(['make_and_analyze_series', 'test', ['.', profiles, True,
        True, 'pdf', '.'], analysis_args])

    while True:
        try:
            ret = ret_q.get_nowait()
            print(ret)
        except Exception:
            pass

        try:
            with log_lock:
                ret = log_q.get_nowait()
                print(ret)
        except Exception:
            pass
