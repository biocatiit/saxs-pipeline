import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

mpl.rc('font', size = 8.0, family='Arial')
mpl.rc('legend', frameon=False, fontsize='medium')
mpl.rc('axes', labelsize='medium', linewidth=1, facecolor='white',
    axisbelow=False, labelpad=2.5, xmargin = 0.015, ymargin = 0.02)
mpl.rc('xtick', labelsize='medium', top=True, direction='in')
mpl.rc('ytick', labelsize='medium', right=True, direction='in')
mpl.rc('lines', linewidth=1)
mpl.rc('mathtext', default='regular')

def make_patch_spines_invisible(ax):
    # from https://matplotlib.org/examples/pylab_examples/multiple_yaxis_with_spines.html
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

class series_overview_plot(object):
    """
    Makes an overview plot with 4 panels. a) Series intensity and rg.
    b) Log-lin profiles. c) Normalized Kratky profiles. d) P(r)
    """

    def __init__(self, profiles, ifts, series, int_type='Total',
        series_data='Rg', img_width=6, img_height=6):

        self.profiles = profiles
        self.ifts = ifts
        self.series = series

        self.int_type = int_type
        self.series_data = series_data

        self.figure = plt.figure(figsize=(img_width, img_height))
        self.gs = self.figure.add_gridspec(3, 2)

        self._make_series_plot()
        self._make_profile_plot()
        self._make_guinier_plot()
        self._make_kratky_plot()
        self._make_ift_plot()

        self.figure.subplots_adjust(left=0.1, right=0.90, wspace=0.3,
            bottom=0.07, top=0.98, hspace=0.3)

    def _make_series_plot(self, row=0):
        ax = self.figure.add_subplot(self.gs[row, :])
        ax2 = ax.twinx()

        ax.spines["right"].set_visible(False)
        make_patch_spines_invisible(ax2)
        ax2.spines["right"].set_visible(True)

        lines = []

        for series in self.series:
            x_data = series.frames

            if self.int_type == 'Total':
                y_data = series.total_i
            elif self.int_type == 'Mean':
                y_data = series.mean_i

            if series.has_calc_data:
                if self.series_data == 'Rg':
                    y2_data = series.rg
                elif self.series_data == 'I0':
                    y2_data = series.i0
                elif self.series_data == 'MW_Vc':
                    y2_data = series.vcmw
                elif self.series_data == 'MW_Vp':
                    y2_data = series.vpmw

            line, = ax.plot(x_data, y_data, '-', label=series.filename)

            if series.has_calc_data:
                line2, = ax2.plot(x_data[y2_data>0], y2_data[y2_data>0], 'o',
                    label='{} {}'.format(series.filename, self.series_data),
                    markersize=1)

            lines.append(line)
            lines.append(line2)

        if len(self.series) == 1:
            series = self.series[0]

            int_line = lines[0]
            ax.yaxis.label.set_color(int_line.get_color())
            ax.tick_params(axis='y', colors=int_line.get_color())
            ax.spines['left'].set_color(int_line.get_color())

            calc_color = '#D7191C'

            if series.has_calc_data:
                calc_line = lines[-1]
                calc_line.set_markerfacecolor(calc_color)
                calc_line.set_markeredgecolor(calc_color)

                ax2.yaxis.label.set_color(calc_color)
                ax2.tick_params(axis='y', colors=calc_color)
                ax2.spines['right'].set_color(calc_color)


            if series.buffer_range:
                for buf in series.buffer_range:
                    ax.axvspan(buf[0], buf[1], color='#2ca02c', alpha=0.5)

            if series.sample_range:
                for sam in series.sample_range:
                    ax.axvspan(sam[0], sam[1], color='#B879CB', alpha=0.5)


        labels = [l.get_label() for l in lines]

        if len(self.series) > 1:
            ax.legend(lines, labels, fontsize='small')

        ax.set_xlabel('Frames')

        ax.set_ylabel('{} Intensity [Arb.]'.format(self.int_type))

        if self.series_data == 'Rg':
            ax2.set_ylabel(r'Rg [$\AA$]')
        elif self.series_data == 'I0':
            ax2.set_ylabel('I(0)')
        elif self.series_data == 'MW_Vc':
            ax2.set_ylabel('MW (Vc) [kDa]')
        elif self.series_data == 'MW_Vp':
            ax2.set_ylabel('MW (Vp) [kDa]')

        ax.text(-0.05, 1.0, 'a', transform = ax.transAxes, fontweight='bold',
            size='large')

    def _make_profile_plot(self, row=1, column=0):
        ax = self.figure.add_subplot(self.gs[row, column])

        ax.set_yscale('log')

        for profile in self.profiles:
            ax.plot(profile.q, profile.i, markersize=1, label=profile.filename)

        if len(self.profiles) > 1:
            ax.legend(fontsize='small')

        absolute = [profile.metadata.Absolute_scale for profile in self.profiles]

        ax.set_xlabel(r'q [$\AA^{-1}$]')

        if all(absolute):
            ax.set_ylabel(r'Intensity [$cm^{-1}$]')
        else:
            ax.set_ylabel('Intensity [Arb.]')

        ax.text(-.15,1.0, 'b', transform = ax.transAxes, fontweight='bold',
            size='large')

    def _make_guinier_plot(self, row=1, column=1):
        gssub = mpl.gridspec.GridSpecFromSubplotSpec(2, 1, self.gs[row, column],
            height_ratios=[1, 0.3], hspace=0.03)

        ax = self.figure.add_subplot(gssub[0])
        res_ax = self.figure.add_subplot(gssub[1], sharex=ax)

        plt.setp(ax.get_xticklabels(), visible=False)


        ax.set_yscale('log')
        # ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(tick_formatter))
        ax.yaxis.set_minor_formatter(plt.NullFormatter())



        for profile in self.profiles:
            fit = guinier_fit(profile.q, profile.guinier_data.Rg,
                profile.guinier_data.I0)
            n_min = profile.guinier_data.n_min
            n_max = profile.guinier_data.n_max

            ax.plot(profile.q[:n_max]**2, profile.i[:n_max], 'o', markersize=3,
                label=profile.filename)
            ax.plot(profile.q[n_min:n_max]**2, fit[n_min:n_max], '-', color='k')
            if n_min > 0:
                ax.plot(profile.q[:n_min]**2, fit[:n_min], '--', color='0.6')

            res_y = (profile.i[n_min:n_max]-fit[n_min:n_max])/profile.err[n_min:n_max]
            res_ax.plot(profile.q[n_min:n_max]**2, res_y, 'o', markersize=3)

        res_ax.axhline(0, color='0')

        if len(self.profiles) > 1:
            ax.legend(fontsize='small')

        absolute = [profile.metadata.Absolute_scale for profile in self.profiles]

        if all(absolute):
            ax.set_ylabel(r'Intensity [$cm^{-1}$]')
        else:
            ax.set_ylabel('Intensity [Arb.]')

        res_ax.set_xlabel(r'q$^2$ [$\AA^{-2}$]')
        res_ax.set_ylabel(r'$\Delta$I/$\sigma$')

        ax.text(-.15,1.0, 'c', transform = ax.transAxes, fontweight='bold',
            size='large')

    def _make_kratky_plot(self, row=2, column=0):
        ax = self.figure.add_subplot(self.gs[row, column])

        ax.axvline(np.sqrt(3), 0, 1, linestyle = 'dashed', color='0.6')
        ax.axhline(3/np.e, 0, 1, linestyle = 'dashed', color='0.6')

        for profile in self.profiles:
            i0 = profile.guinier_data.I0
            rg = profile.guinier_data.Rg

            qRg = profile.q*rg

            ax.plot(qRg, qRg**2*profile.i/i0, markersize=1,
                label=profile.filename)

        ax.axhline(0, color='k')

        if len(self.profiles) > 1:
            ax.legend(fontsize='small')

        ax.set_xlabel(r'q$R_g$')
        ax.set_ylabel(r'(q$R_g$)$^2$I(q)/I(0)')

        ymin, ymax = ax.get_ylim()
        ax.set_ylim(max(-0.1, ymin), ymax)

        ax.text(-.15,1.0, 'd', transform = ax.transAxes, fontweight='bold',
            size='large')

    def _make_ift_plot(self, row=2, column=1):
        ax = self.figure.add_subplot(self.gs[row, column])
        ax.axhline(0, color='k')
        # plt.setp(ax.get_yticklabels(), visible=False)

        for ift in self.ifts:
            ax.plot(ift.r, ift.p, markersize=1, label=ift.filename)

        if len(self.profiles) > 1:
            ax.legend(fontsize='small')

        ax.set_xlabel(r'r [$\AA$]')
        ax.set_ylabel('P(r)/I(0)')

        ax.text(-.15,1.0, 'e', transform = ax.transAxes, fontweight='bold',
            size='large')

class efa_plot(object):
    """
    Makes an overview plot with 4 panels. a) Series intensity and rg.
    b) Log-lin profiles. c) Normalized Kratky profiles. d) P(r)
    """

    def __init__(self, series, efa_data, int_type='Total',
        series_data='Rg', img_width=6, img_height=6):

        self.series = series
        self.efa_data = efa_data

        self.int_type = int_type
        self.series_data = series_data


        if efa_data is None:
            if img_width == 6 and img_height == 6:
                self.figure = plt.figure(figsize=(6, 2))
            else:
                self.figure = plt.figure(figsize=(img_width, img_height))

            self.gs = self.figure.add_gridspec(1, 2)

        else:
            self.figure = plt.figure(figsize=(img_width, img_height))

            self.gs = self.figure.add_gridspec(3, 2)

        self._make_series_plot()
        self._make_efa_range_plot()

        self.figure.subplots_adjust(left=0.1, right=0.90, wspace=0.3,
            bottom=0.07, top=0.98, hspace=0.3)

    def _make_series_plot(self, row=0, column=0):
        ax = self.figure.add_subplot(self.gs[row, column])
        ax2 = ax.twinx()

        ax.spines["right"].set_visible(False)
        make_patch_spines_invisible(ax2)
        ax2.spines["right"].set_visible(True)

        x_data = self.series.frames

        if self.int_type == 'Total':
            y_data = self.series.total_i
        elif self.int_type == 'Mean':
            y_data = self.series.mean_i

        if self.series.has_calc_data:
            if self.series_data == 'Rg':
                y2_data = self.series.rg
            elif self.series_data == 'I0':
                y2_data = self.series.i0
            elif self.series_data == 'MW_Vc':
                y2_data = self.series.vcmw
            elif self.series_data == 'MW_Vp':
                y2_data = self.series.vpmw

        int_line, = ax.plot(x_data, y_data, '-', label=self.series.filename,
            color='k')

        if self.series.has_calc_data:
            calc_line, = ax2.plot(x_data[y2_data>0], y2_data[y2_data>0], 'o',
                label='{} {}'.format(self.series.filename, self.series_data),
                markersize=1)


        ax.yaxis.label.set_color(int_line.get_color())
        ax.tick_params(axis='y', colors=int_line.get_color())
        ax.spines['left'].set_color(int_line.get_color())

        calc_color = '#D7191C'

        if self.series.has_calc_data:
            calc_line.set_markerfacecolor(calc_color)
            calc_line.set_markeredgecolor(calc_color)

        ax2.yaxis.label.set_color(calc_color)
        ax2.tick_params(axis='y', colors=calc_color)
        ax2.spines['right'].set_color(calc_color)

        ax.set_xlabel('Frames')

        ax.set_ylabel('{} Intensity [Arb.]'.format(self.int_type))

        if self.series_data == 'Rg':
            ax2.set_ylabel(r'Rg [$\AA$]')
        elif self.series_data == 'I0':
            ax2.set_ylabel('I(0)')
        elif self.series_data == 'MW_Vc':
            ax2.set_ylabel('MW (Vc) [kDa]')
        elif self.series_data == 'MW_Vp':
            ax2.set_ylabel('MW (Vp) [kDa]')

        ax.text(-0.15, 1.0, 'a', transform = ax.transAxes, fontweight='bold',
            size='large')

    def _make_efa_range_plot(self, row=0, column=1):
        start = int(self.series.efa_start)
        end = int(self.series.efa_end)
        ranges = self.series.efa_ranges

        frame_data = self.series.frames

        if self.int_type == 'Total':
            int_data = self.series.total_i
        elif self.int_type == 'Mean':
            int_data = self.series.mean_i

        ax = self.figure.add_subplot(self.gs[row, column])
        plt.setp(ax.get_yticklabels(), visible=False)

        ax.plot(frame_data[start:end+1], int_data[start:end+1], '-', color='k')
        ax.set_prop_cycle(None)

        for i in range(len(ranges)):
            color = next(ax._get_lines.prop_cycler)['color']

            ax.annotate('', xy=(ranges[i][0], 0.975-0.05*(i)),
                xytext=(ranges[i][1], 0.975-0.05*(i)),
                xycoords=('data', 'axes fraction'),
                arrowprops = dict(arrowstyle = '<->', color = color))

            ax.axvline(ranges[i][0], 0, 0.975-0.05*(i), linestyle='dashed',
                color=color)
            ax.axvline(ranges[i][1], 0, 0.975-0.05*(i), linestyle='dashed',
                color=color)


def guinier_fit(q, rg, i0):
    return i0*np.exp(-rg**2*q**2/3)

def tick_formatter(x, pos):
    return "{}".format(x)
