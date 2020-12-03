import os
import copy

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import bioxtasraw.RAWAPI as raw
import bioxtasraw.SASFileIO as SASFileIO

from ..reports import data as report_data

def find_dmax(profile, settings):
    analysis_dict = profile.getParameter('analysis')
    try:
        rg = float(analysis_dict['guinier']['Rg'])
    except Exception:
        rg = -1

    if rg != -1:

        # Get Dmax from DATCLASS
        try:
            dc_dmax = float(analysis_dict['molecularWeight']['ShapeAndSize']['Dmax'])
        except Exception:
            dc_dmax = -1

        if dc_dmax == -1:
            try:
                dc_mw, dc_shape, dc_dmax = raw.mw_datclass(profile)
            except Exception:
                dc_dmax = -1

        if dc_dmax != -1:
            dmax = int(round(dc_dmax))

        else:
            #Calculate the IFT using BIFT
            try:
                bift_dmax = float(analysis_dict['BIFT']['Dmax'])
            except Exception:
                bift_dmax = -1

            if bift_dmax == -1:
                try:
                    (bift, bift_dmax, bift_rg, bift_i0, bift_dmax_err,
                    bift_rg_err, bift_i0_err, bift_chi_sq, bift_log_alpha,
                    bift_log_alpha_err, bift_evidence,
                    bift_evidence_err) = raw.bift(profile, settings=settings)
                except Exception:
                    bift_dmax = -1

            #Calculate the IFT using DATGNOM
            try:
                (datgnom_ift, datgnom_dmax, datgnom_rg, datgnom_i0,
                datgnom_rg_err, datgnom_i0_err, datgnom_total_est,
                datgnom_chi_sq, datgnom_alpha, datgnom_quality) = raw.datgnom(profile)
            except Exception:
                datgnom_dmax = -1

            if bift_dmax != -1 and datgnom_dmax != -1:
                dmax = np.mean([bift_dmax, datgnom_dmax])

            elif bift_dmax != -1:
                dmax = bift_dmax*0.79

            elif datgnom_dmax != -1:
                dmax = datgnom_dmax*1.2

            else:
                dmax = -1

            if dmax != -1:
                dmax = round(dmax)

        if dmax != -1:
            # Refine if Dmax is too long
            ift_results = raw.gnom(profile, dmax)
            ift = ift_results[0]
            dmax_start = dmax

            while dmax > dmax_start*0.5 and np.any(ift.p[-20:] < 0):
                dmax = dmax -1
                ift_results = raw.gnom(profile, dmax)
                ift = ift_results[0]

            if dmax_start == dmax:
                #Refine if Dmax is too short
                ift_results = raw.gnom(profile, dmax, dmax_zero=False)
                ift_unforced = ift_results[0]
                dmax_start = dmax

                while dmax < dmax_start*1.5 and ift_unforced.p[-1]>0.01*ift_unforced.p.max():
                    dmax = dmax +1
                    ift_results = raw.gnom(profile, dmax, dmax_zero=False)
                    ift_unforced = ift_results[0]

    else:
        dmax = -1

    return dmax

def calc_mw(profile, settings):
    if 'Conc' in profile.getAllParameters():
        try:
            conc = float(profile.getParameter('Conc'))
        except Exception:
            conc = -1

    else:
        conc = -1

    if conc != -1:
        raw.mw_ref(profile, conc, settings=settings)
        raw.mw_abs(profile, conc, settings=settings)

    raw.mw_vp(profile, settings=settings)
    raw.mw_vc(profile, settings=settings)
    raw.mw_bayes(profile)
    raw.mw_datclass(profile)

def run_dammif(ift, settings, out_dir, nruns=5, mode='Fast', average=True,
    cluster=False, refine=False):
    #Create individual bead model reconstructions
    filename = ift.getParameter('filename').replace(' ', '_')
    prefix = os.path.splitext(filename)[0][:30]

    temp_ift = copy.deepcopy(ift)
    temp_ift.setParameter('filename', prefix+'.out')

    dammif_dir = os.path.join(out_dir, os.path.splitext(filename)[0]+'_dammif')
    if not os.path.exists(dammif_dir):
        os.mkdir(dammif_dir)
    else:
        files = os.listdir(dammif_dir)
        for f in files:
            os.remove(os.path.join(dammif_dir, f))

    chi_sq_vals = []
    rg_vals = []
    dmax_vals = []
    mw_vals = []
    ev_vals = []
    model_data = []

    write_ift = False
    ift_name = '{}.out'.format(prefix)

    raw.save_ift(temp_ift, datadir=dammif_dir)

    for i in range(nruns):
        chi_sq, rg, dmax, mw, ev = raw.dammif(temp_ift,
            '{}_{:02d}'.format(prefix, i+1), dammif_dir, mode=mode,
            write_ift=write_ift, ift_name=ift_name)

        chi_sq_vals.append(chi_sq)
        rg_vals.append(rg)
        dmax_vals.append(dmax)
        mw_vals.append(mw)
        ev_vals.append(ev)


    #Average the bead model reconstructions
    damaver_files = ['{}_{:02d}-1.pdb'.format(prefix, i+1) for i in range(nruns)]
    if average and nruns>1:
        (mean_nsd, stdev_nsd, rep_model, result_dict, res, res_err,
            res_unit) = raw.damaver(damaver_files, prefix, dammif_dir)

        nsd_inc = 0
        ex_items = []

        for key in result_dict:
            if result_dict[key][0].lower() == 'include':
                nsd_inc += 1
            else:
                ex_items.append(key)

        r_num = int(os.path.splitext(rep_model)[0].split('_')[-1].split('-')[0])
        nsd_data = [('Mean NSD:', mean_nsd),
                ('Stdev. NSD:', stdev_nsd),
                ('DAMAVER Included:', nsd_inc, 'of', nruns),
                ('Representative model:', r_num),
                ]

        if ex_items:
            nsd_data.append(('Excluded Models:', ' ,'.join(ex_items)))

        res_data = [('Ensemble resolution:', res, '+/-', res_err, res_unit)]

    else:
        nsd_data = []
        res_data = []

    for i, name in enumerate(damaver_files):
        if average and nruns>1:
            nsd = result_dict[name][1]
        else:
            nsd = ''
        model_data.append([i+1, chi_sq_vals[i], rg_vals[i], dmax_vals[i],
            ev_vals[i], mw_vals[i], nsd])

    if average and nruns>1:
        damaver_name = os.path.join(dammif_dir, prefix+'_damaver.pdb')
        damfilt_name = os.path.join(dammif_dir, prefix+'_damfilt.pdb')

        atoms, header, model_info= SASFileIO.loadPDBFile(damaver_name)
        model_data.append(['damaver', model_info['chisq'], model_info['rg'],
            model_info['dmax'], model_info['excluded_volume'], model_info['mw'],
            ''])

        atoms, header, model_info = SASFileIO.loadPDBFile(damfilt_name)
        model_data.append(['damfilt', model_info['chisq'], model_info['rg'],
            model_info['dmax'], model_info['excluded_volume'], model_info['mw'],
            ''])

    #Cluster the bead model reconstructions
    if cluster and nruns > 1:
        cluster_list, distance_list = raw.damclust(damaver_files, prefix,
            dammif_dir)

        clust_num = ('Number of clusters:', len(cluster_list))

        clist_data = []

        for cl in cluster_list:
            if cl.dev == -1:
                isolated = 'Y'
                dev = ''
            else:
                isolated = 'N'
                dev = str(cl.dev)

            clist_data.append([cl.num, isolated, cl.rep_model, dev])

        dlist_data = [list(map(str, d_data)) for d_data in distance_list]

    else:
        clust_num = []
        clist_data = []
        dlist_data = []

    #Refine the bead model
    do_refine = average and refine and nruns > 1
    if do_refine:
        chi_sq, rg, dmax, mw, ev = raw.dammin(temp_ift, 'refine_{}'.format(prefix),
            dammif_dir, 'Refine', initial_dam='{}_damstart.pdb'.format(prefix),
            write_ift=False, ift_name=ift_name)

        chi_sq_vals.append(chi_sq)
        rg_vals.append(rg)
        dmax_vals.append(dmax)
        mw_vals.append(mw)
        ev_vals.append(ev)

        model_data.append(['refine', chi_sq, rg, dmax, ev, mw, ''])

    ambi_score, ambi_cats, ambi_eval = raw.ambimeter(temp_ift, datadir=dammif_dir,
        write_ift=False, filename=ift_name)

    ambi_data = [('Compatible shape categories:', ambi_cats),
        ('Ambiguity score:', ambi_score), ('AMBIMETER says:', ambi_eval)]

    setup_data = [('Input file:', filename), ('Output prefix:', prefix),
        ('Output directory:', os.path.abspath(os.path.expanduser(dammif_dir))),
        ('Program used:', 'DAMMIF'), ('Mode:', mode), ('Symmetry:', 'P1'),
        ('Anisometry:', 'Unknown'), ('Total number of reconstructions:', nruns),
        ('Used DAMAVER:', average and nruns>1), ('Refined with DAMMIN:', do_refine),
        ('Used DAMCLUST:', cluster and nruns>1),
        ]

    model_plots = make_dammif_figures(dammif_dir, prefix, nruns, do_refine)

    save_path = os.path.join(dammif_dir, '{}_dammif_results.csv'.format(prefix))

    save_data = SASFileIO.saveDammixData(save_path, ambi_data, nsd_data,
        res_data, clust_num, clist_data, dlist_data, model_data, setup_data,
        model_plots)

    dammif_data = report_data.parse_dammif_file(None, save_data.split('\n'))

    return dammif_data

def make_dammif_figures(path, prefix, nruns, refine):
    figures = []
    for num in range(1, nruns+1):
        fprefix = '%s_%s' %(prefix, str(num).zfill(2))
        fir_name = os.path.join(path, fprefix+'.fir')

        sasm, fit_sasm = SASFileIO.loadFitFile(fir_name)

        fig = plt.Figure((5,4), 75)

        gridspec = mpl.gridspec.GridSpec(2, 1, height_ratios=[1, 0.3])

        q = sasm.getQ()
        i = sasm.getI()
        i_fit = fit_sasm.getI()
        err = sasm.getErr()

        residual = i - i_fit
        residual = residual/err

        ax0 = fig.add_subplot(gridspec[0])
        ax0.semilogy(q, i, 'bo')
        ax0.semilogy(q, i_fit, 'r')
        ax0.set_xlabel('q')
        ax0.set_ylabel('I(q)')

        ax1 = fig.add_subplot(gridspec[1])
        ax1.axhline(0, color='k', linewidth=1.0)
        ax1.plot(q, residual, 'bo')
        ax1.set_xlabel('q')
        ax1.set_ylabel('$\Delta I(q)/\sigma (q)$')

        fig.subplots_adjust(left = 0.1, bottom = 0.12, right = 0.95, top = 0.95,
            hspace=0.25)
        fig.set_facecolor('white')

        figures.append(['{}'.format(num+1), [fig]])

    if refine:
        fir_name = os.path.join(path, 'refine_'+prefix+'.fir')

        sasm, fit_sasm = SASFileIO.loadFitFile(fir_name)

        fig = plt.Figure((5,4), 75)

        gridspec = mpl.gridspec.GridSpec(2, 1, height_ratios=[1, 0.3])

        q = sasm.getQ()
        i = sasm.getI()
        i_fit = fit_sasm.getI()
        err = sasm.getErr()

        residual = i - i_fit
        residual = residual/err

        ax0 = fig.add_subplot(gridspec[0])
        ax0.semilogy(q, i, 'bo')
        ax0.semilogy(q, i_fit, 'r')
        ax0.set_xlabel('q')
        ax0.set_ylabel('I(q)')

        ax1 = fig.add_subplot(gridspec[1])
        ax1.axhline(0, color='k', linewidth=1.0)
        ax1.plot(q, residual, 'bo')
        ax1.set_xlabel('q')
        ax1.set_ylabel('$\Delta I(q)/\sigma (q)$')

        fig.subplots_adjust(left = 0.1, bottom = 0.12, right = 0.95, top = 0.95,
            hspace=0.25)
        fig.set_facecolor('white')

        figures.append(['refined', [fig]])

    return figures

def run_denss(ift, settings, out_dir, nruns=5, mode='Fast', average=True,
    refine=False):
    #Create individual bead model reconstructions
    filename = ift.getParameter('filename').replace(' ', '_')
    prefix = os.path.splitext(filename)[0]

    denss_dir = os.path.join(out_dir, os.path.splitext(filename)[0]+'_denss')
    if not os.path.exists(denss_dir):
        os.mkdir(denss_dir)
    else:
        files = os.listdir(denss_dir)
        for f in files:
            os.remove(os.path.join(denss_dir, f))

    rhos = []
    chi_vals = []
    rg_vals = []
    support_vol_vals = []
    sides = []
    fit_data = []

    denss_results = []

    for i in range(nruns):
        (rho, chi_sq, rg, support_vol, side, q_fit, I_fit, I_extrap,
            err_extrap, all_chi_sq, all_rg, all_support_vol) = raw.denss(ift,
            '{}_{:02d}'.format(prefix, i+1), denss_dir, mode=mode)

        rhos.append(rho)
        chi_vals.append(chi_sq)
        rg_vals.append(rg)
        support_vol_vals.append(support_vol)
        sides.append(side)
        fit_data.append([q_fit, I_fit, I_extrap, err_extrap])

        denss_results.append([q_fit, I_extrap, err_extrap, q_fit, I_fit,
            all_chi_sq, all_rg, all_support_vol])


    #Average the electron reconstructions
    rsc_data = []

    if average and nruns > 1:
        (average_rho, mean_cor, std_cor, threshold, res, scores,
            fsc) = raw.denss_average(np.array(rhos), side,
            '{}_average'.format(prefix), denss_dir)

        rsc_data.append(('Mean RSC:', round(mean_cor, 4)))
        rsc_data.append(('Stdev. RSC:', round(std_cor, 4)))
        rsc_data.append(('Number of models included:', np.sum(scores>threshold)))
        rsc_data.append(('Total number of models:', nruns))

        if np.sum(scores<=threshold) > 0:
            rsc_data.append(('Excluded Models:', ' ,'.join(np.argwhere(scores<=threshold))))

        res_data = [('Fourier Shell Correlation Resolution (Angstrom):', res)]

        average_results = {'fsc': fsc}

    else:
        rsc_data = []
        res_data = []

        average_results = None

    model_data = []

    if average and nruns > 1:
        for i, score in enumerate(scores):
            model_data.append([i+1, round(chi_vals[i], 5), round(rg_vals[i], 2),
                round(support_vol_vals[i], 2), round(score, 4)])
    else:
        for i in range(nruns):
            model_data.append([i+1, round(chi_vals[i], 5), round(rg_vals[i], 2),
                round(support_vol_vals[i], 2), ''])

    #Refine the electron density
    do_refine = refine and average and nruns > 1

    if do_refine:
        (refined_rho, refined_chi_sq, refined_rg, refined_support_vol,
            refined_side, refined_q_fit, refined_I_fit, refined_I_extrap,
            refined_err_extrap, refined_all_chi_sq, refined_all_rg,
            refined_all_support_vol) = raw.denss(ift, '{}_refine'.format(prefix),
            denss_dir, mode=mode, initial_model=average_rho)

        model_data.append(['Refine', round(refined_chi_sq, 5), round(refined_rg, 2),
            round(refined_support_vol, 2), ''])

        refined_results = [refined_q_fit, refined_I_extrap,
            refined_err_extrap, refined_q_fit, refined_I_fit,
            refined_all_chi_sq, refined_all_rg, refined_all_support_vol]

    else:
        refined_results = None

    if ift.getParameter('algorithm') == 'GNOM':
        ambi_score, ambi_cats, ambi_eval = raw.ambimeter(ift)

        ambi_data = [('Compatible shape categories:', ambi_cats),
            ('Ambiguity score:', ambi_score), ('AMBIMETER says:', ambi_eval)]

    else:
        ambi_data = []

    setup_data = [('Input file:', filename),
        ('Output prefix:', prefix),
        ('Output directory:', os.path.abspath(os.path.expanduser(denss_dir))),
        ('Mode:', mode),
        ('Total number of reconstructions:', nruns),
        ('Number of electrons in molecule:', ''),
        ('Averaged:', average and nruns > 1),
        ('Refined:', do_refine),
        ('Symmetry applied:', 'False')
        ]

    model_plots = make_denss_figures(ift, denss_results, average_results,
        refined_results)

    filename = prefix + '_denss_results.csv'
    save_path = os.path.join(denss_dir, filename)

    save_data = SASFileIO.saveDenssData(save_path, ambi_data, res_data,
        model_plots, setup_data, rsc_data, model_data)

    denss_data = report_data.parse_denss_file(None, save_data.split('\n'))

    return denss_data

def make_denss_figures(iftm, denss_results_list, average_results=None,
    refine_results=None):

    figures = []

    for num, denss_results in enumerate(denss_results_list):
        fig, fig2 = denss_plot(iftm, denss_results)

        figures.append(['{}'.format(num+1), [fig, fig2]])

    if average_results is not None:
        fig = plt.Figure((3.25,2.5))

        res = average_results['fsc'][:, 0]
        fsc = average_results['fsc'][:, 1]

        ax0 = fig.add_subplot(111)
        ax0.plot(res, fsc, 'bo-')
        ax0.set_xlabel('Resolution ($\\AA^{-1}$)', fontsize='small')
        ax0.set_ylabel('Fourier Shell Correlation', fontsize='small')
        ax0.tick_params(labelsize='x-small')

        fig.subplots_adjust(left = 0.1, bottom = 0.12, right = 0.95, top = 0.95)
        fig.set_facecolor('white')

        figures.append(['Average', [fig]])

    if refine_results is not None:
        fig, fig2 = denss_plot(iftm, denss_results)

        figures.append(['Refine', [fig, fig2]])


    return figures

def denss_plot(iftm, denss_results):
    fig = plt.Figure((3.25,2.5))

    if iftm.getParameter('algorithm') == 'GNOM':
        q = iftm.q_extrap
        I = iftm.i_extrap

        ext_pts = len(I)-len(iftm.i_orig)
        sigq = np.empty_like(I)
        sigq[:ext_pts] = I[:ext_pts]*np.mean((iftm.err_orig[:10]/iftm.i_orig[:10]))
        sigq[ext_pts:] = I[ext_pts:]*(iftm.err_orig/iftm.i_orig)
    else:
        q = iftm.q_orig
        I = iftm.i_fit
        sigq = I*(iftm.err_orig/iftm.i_orig)
    #handle sigq values whose error bounds would go negative and be missing on the log scale
    sigq2 = np.copy(sigq)
    sigq2[sigq>I] = I[sigq>I]*.999

    qdata = denss_results[0]
    Idata = denss_results[1]
    qbinsc = denss_results[3]
    Imean = denss_results[4]

    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[3,1])

    ax0 = fig.add_subplot(gs[0])
    ax0.plot(q[q<=qdata[-1]], I[q<=qdata[-1]], 'k-', label='Smoothed Exp. Data')
    ax0.plot(qdata, Idata, 'bo',alpha=0.5,label='Interpolated Data')
    ax0.plot(qbinsc, Imean,'r.',label='Scattering from Density')
    handles,labels = ax0.get_legend_handles_labels()
    handles = [handles[2], handles[0], handles[1]]
    labels = [labels[2], labels[0], labels[1]]
    ymin = np.min(np.hstack((I,Idata,Imean)))
    ymax = np.max(np.hstack((I,Idata,Imean)))
    ax0.set_ylim([0.5*ymin,1.5*ymax])
    ax0.legend(handles,labels, fontsize='small')
    ax0.semilogy()
    ax0.set_ylabel('I(q)', fontsize='small')
    ax0.tick_params(labelbottom=False, labelsize='x-small')

    ax1 = fig.add_subplot(gs[1])
    ax1.axhline(0, color='k', linewidth=1.0)
    ax1.plot(qdata, np.log10(Imean[qbinsc==qdata])-np.log10(Idata), 'ro-')
    ylim = ax1.get_ylim()
    ymax = np.max(np.abs(ylim))
    ax1.set_ylim([-ymax,ymax])
    ax1.yaxis.major.locator.set_params(nbins=5)
    xlim = ax0.get_xlim()
    ax1.set_xlim(xlim)
    ax1.set_ylabel('Residuals', fontsize='small')
    ax1.set_xlabel(r'q ($\mathrm{\AA^{-1}}$)', fontsize='small')
    ax1.tick_params(labelsize='x-small')

    fig.subplots_adjust(left = 0.2, bottom = 0.15, right = 0.95, top = 0.95)
    fig.set_facecolor('white')


    fig2 = plt.Figure((3.25,2.5))

    chi = denss_results[5]
    rg = denss_results[6]
    vol = denss_results[7]

    ax0 = fig2.add_subplot(311)
    ax0.plot(chi[chi>0])
    ax0.set_ylabel('$\chi^2$', fontsize='small')
    ax0.semilogy()
    ax0.tick_params(labelbottom=False, labelsize='x-small')

    ax1 = fig2.add_subplot(312)
    ax1.plot(rg[rg>0])
    ax1.set_ylabel('Rg', fontsize='small')
    ax1.tick_params(labelbottom=False, labelsize='x-small')

    ax2 = fig2.add_subplot(313)
    ax2.plot(vol[vol>0])
    ax2.set_xlabel('Step', fontsize='small')
    ax2.set_ylabel('Support Volume ($\mathrm{\AA^{3}}$)', fontsize='small')
    ax2.semilogy()
    ax2.tick_params(labelsize='x-small')

    fig2.subplots_adjust(left = 0.2, bottom = 0.15, right = 0.95, top = 0.95)
    fig2.set_facecolor('white')

    return fig, fig2

def model_free_analysis(profile, settings):
    raw_results = raw.auto_guinier(profile, settings=settings)

    raw_rg = raw_results[0]

    if raw_rg != -1:
        calc_mw(profile, settings)

        dmax = find_dmax(profile, settings)

        if dmax != -1:
            ift_results = raw.gnom(profile, dmax)
            ift = ift_results[0]

        else:
            ift = None

    else:
        ift = None

    return profile, ift

def model_based_analysis(ift, settings, out_dir, dammif=True, denss=True,
    dam_runs=3, dam_mode='Fast', dam_aver=True, dam_clust=False,
    dam_refine=False, denss_runs=3, denss_mode='Fast', denss_aver=True,
    denss_refine=False):
    if dammif:
        dammif_data = run_dammif(ift, settings, out_dir, nruns=dam_runs,
            mode=dam_mode, average=dam_aver, cluster=dam_clust,
            refine=dam_refine)

    else:
        dammif_data = None

    if denss:
        denss_data = run_denss(ift, settings, out_dir, nruns=denss_runs,
            mode=denss_mode, average=denss_aver, refine=denss_refine)

    else:
        denss_data = None

    return dammif_data, denss_data


def analyze_files(out_dir, data_dir, profiles, ifts, raw_settings, dammif=True,
    denss=True, dam_runs=3, dam_mode='Fast', dam_aver=True, dam_clust=False,
    dam_refine=False, denss_runs=3, denss_mode='Fast', denss_aver=True,
    denss_refine=False):
    data_dir = os.path.abspath(os.path.expanduser(data_dir))

    for j in range(len(profiles)):
        profiles[j] = os.path.abspath(os.path.expanduser(os.path.join(data_dir,
            profiles[j])))

    for j in range(len(ifts)):
        ifts[j] = os.path.abspath(os.path.expanduser(os.path.join(data_dir,
            ifts[j])))

    profiles = raw.load_profiles(profiles)


    new_ifts = []

    all_dammif_data = []
    all_denss_data = []

    for i in range(len(profiles)):
        profiles[i], ift = model_free_analysis(profiles[i], raw_settings)

        if ift is not None:
            new_ifts.append(ift)

            if dammif:
                ift_results = raw.gnom(profiles[i], ift.getParameter('dmax'),
                    cut_dam=True)
                cut_ift = ift_results[0]

                dammif_data, _ = model_based_analysis(cut_ift, raw_settings,
                    out_dir, dammif=dammif, denss=False, dam_runs=dam_runs,
                    dam_aver=dam_aver, dam_clust=dam_clust,
                    dam_refine=dam_refine, denss_runs=denss_runs,
                    denss_aver=denss_aver, denss_refine=denss_refine)

                if dammif_data is not None:
                    all_dammif_data.append(dammif_data)

            if denss:
                _, denss_data = model_based_analysis(ift, raw_settings, out_dir,
                    dammif=False, denss=denss, dam_runs=dam_runs,
                    dam_aver=dam_aver, dam_clust=dam_clust,
                    dam_refine=dam_refine, denss_runs=denss_runs,
                    denss_aver=denss_aver, denss_refine=denss_refine)

                if denss_data is not None:
                    all_denss_data.append(denss_data)


    ifts = raw.load_ifts(ifts)

    for ift in ifts:
        dammif_data, denss_data = model_based_analysis(ift, raw_settings,
            out_dir, dammif=dammif, denss=denss, dam_runs=dam_runs,
            dam_aver=dam_aver, dam_clust=dam_clust, dam_refine=dam_refine,
            denss_runs=denss_runs, denss_aver=denss_aver,
            denss_refine=denss_refine)

        if dammif_data is not None:
            all_dammif_data.append(dammif_data)

        if denss_data is not None:
            all_denss_data.append(denss_data)

    all_ifts = new_ifts + ifts

    return profiles, all_ifts, all_dammif_data, all_denss_data


def analyze_data(out_dir, profiles, ifts, raw_settings, dammif=True,
    denss=False, dam_runs=3, dam_mode='Fast', dam_aver=True, dam_clust=False,
    dam_refine=False, denss_runs=3, denss_mode='Fast', denss_aver=True,
    denss_refine=False):

    new_ifts = []

    all_dammif_data = []
    all_denss_data = []

    for i in range(len(profiles)):
        profiles[i], ift = model_free_analysis(profiles[i], raw_settings)

        if ift is not None:
            new_ifts.append(ift)

            if dammif:
                ift_results = raw.gnom(profiles[i], ift.getParameter('dmax'),
                    cut_dam=True)
                cut_ift = ift_results[0]

                dammif_data, _ = model_based_analysis(cut_ift, raw_settings,
                    out_dir, dammif=dammif, denss=False, dam_runs=dam_runs,
                    dam_aver=dam_aver, dam_clust=dam_clust,
                    dam_refine=dam_refine, denss_runs=denss_runs,
                    denss_aver=denss_aver, denss_refine=denss_refine)

                if dammif_data is not None:
                    all_dammif_data.append(dammif_data)

            if denss:
                _, denss_data = model_based_analysis(ift, raw_settings, out_dir,
                    dammif=False, denss=denss, dam_runs=dam_runs,
                    dam_aver=dam_aver, dam_clust=dam_clust,
                    dam_refine=dam_refine, denss_runs=denss_runs,
                    denss_aver=denss_aver, denss_refine=denss_refine)

                if denss_data is not None:
                    all_denss_data.append(denss_data)

    for ift in ifts:
        dammif_data, denss_data = model_based_analysis(ift, raw_settings,
            out_dir, dammif=dammif, denss=denss, dam_runs=dam_runs,
            dam_aver=dam_aver, dam_clust=dam_clust, dam_refine=dam_refine,
            denss_runs=denss_runs, denss_aver=denss_aver,
            denss_refine=denss_refine)

        if dammif_data is not None:
            all_dammif_data.append(dammif_data)

        if denss_data is not None:
            all_denss_data.append(denss_data)

    all_ifts = new_ifts + ifts

    return profiles, all_ifts, all_dammif_data, all_denss_data
