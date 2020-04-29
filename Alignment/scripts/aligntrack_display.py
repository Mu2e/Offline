
import sys
import uproot


import numpy as np
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec
from matplotlib import rc
import matplotlib.backends.backend_pdf

mpl.use('Agg')
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage[cm]{sfmath}']
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'cm'


def gauss(x, a, mean, sigma):
    return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))

def hist_meanresid(ax, hist, label=''):
    y, bin_edges = hist

    bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])

    line = ax.errorbar(
        bin_centers,
        y,
        yerr = y*0,#**0.5,
        marker = 'x',
        fmt='x-',
        linewidth=0.5,
        drawstyle = 'steps-mid',
        label=label
    )
    ax.set_xlim(bin_edges[0], bin_edges[-1]);
    ax.set_ylim(-50, 50)
    ax.hlines(0, 0,36,linestyles='dashed')
    ax.set_title('Mean residuals per plane')

    ax.set_xlabel('Plane ID', fontsize=9)
    ax.set_ylabel('Mean residuals over plane ($\\mu m$)', fontsize=9)
    ax.xaxis.set_ticks(np.arange(0, 36.2, 1))
    ax.tick_params(which='major', direction='in')
    ax.xaxis.grid(which='major', linestyle='--')

    return line[0]

def hist_poiserr(ax, hist, label='', fitgauss=False):
    y, bin_edges = hist

    bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])
    lab = 'no fit'

    if fitgauss:
        try:
            (A,mu,sigma), pcov = curve_fit(gauss, bin_centers, y, )
            lab = '$\\mu = %.2f, \\sigma = %.2f$' % (mu,sigma)

            ax.plot(bin_centers, gauss(bin_centers, A,mu,sigma),  'r-.', lw=0.5)
        except:
            pass


    line = ax.errorbar(
        bin_centers,
        y,
        yerr = y**0.5,
        marker = '.',
        fmt='.-',
        linewidth=0.5,
        drawstyle = 'steps-mid',
        label=lab
    )

    ax.hlines(0, np.min(bin_edges), np.max(bin_edges), linestyles='dashed')

    ax.set_xlim(bin_edges[0], bin_edges[-1]);
    ax.set_ylabel('Entries / %.2f' % (abs(bin_edges[0] - bin_edges[-1])/(len(bin_edges)-1)))

    ymin, ymax = ax.get_ylim()
    this_ymin, this_ymax = 0, y.max() * 1.25

    if this_ymin < ymin:
        ymin = this_ymin
    if this_ymax > ymax:
        ymax = this_ymax

    ax.set_ylim(ymin, ymax)
    ax.xaxis.set_label_coords(0.95, -0.125)
    ax.minorticks_on()
    ax.tick_params(which='major', direction='in')
    ax.tick_params(which='minor', direction='in')
    ax.grid(which='major')

    ax.legend(fontsize="xx-small")

    return line[0]

def main():
    filename = sys.argv[1]

    fig = plt.figure(figsize=(20.841,12.195), dpi=100)
    gs = gridspec.GridSpec(4, 5, hspace=0.4)

    evt_counts = {}
    lines={}
    for filename in sys.argv[1:]:
        events = uproot.open(filename)["AlignTrackCollector/tracks"]

        evt_counts[filename] = len(events)
        filename += ' (%d selected events)' % evt_counts[filename]

        variables_eventlevel = {
            "chisq": (50, 0, 5, False),
            "chisq_doca": (50, 0, 5, False),

            "pvalue": (20, 0, 1, False),

            "nHits": (50, 0, 50, False),

            "planes_trav": (10, 0, 10, False),
            "panels_trav": (20, 0, 20, False),
        }

        variables_hitlevel = {

            "doca_resid": (50, -1.5, 1.5, True),
            "time_resid": (50, -15,15, True),
            "doca_resid_err": (50, 0, 0.3, False),
            "doca": (50, 0, 3, False),
            "time": (50, 0, 2000, False),

            "drift_res": (50, 0,20, False),

            "plane": (36, 0, 36, False),
            "panel": (216, 0, 216, False),

            "pull_doca": (50, -7.5, 7.5, True),
            "pull_hittime": (50, -7.5, 7.5, True)

        }

        allvars = {**variables_eventlevel, **variables_hitlevel}

        cuts = {
           #"doca_cut": lambda df_ev,df_hit: (df_ev, df_hit[df_hit['doca'] > 0.15])
           #"pval_cut": lambda df_ev,df_hit: (df_ev[df_ev['pvalue'] < 0.1], df_hit[df_hit['pvalue'] < 0.1])
        }

        plots = {}
        planeresid = None
        df_evt_only = events.pandas.df("*", namedecode="utf-8", flatten=None)
        df_hitlevel = events.pandas.df("*", namedecode="utf-8", flatten=True)

        df_evt_only['chisq_doca'] = df_evt_only['chisq_doca'] / df_evt_only['ndof']

        for _, cut_fn in cuts.items():
            df_evt_only, df_hitlevel = cut_fn(df_evt_only, df_hitlevel)

        # N.B. residuals here are converted to microns
        counts, edges = np.histogram(df_hitlevel["plane"], weights=(df_hitlevel["doca_resid"] * 1000.0), bins=36, range=(-0.5,35.5))
        counts = np.divide(counts, np.histogram(df_hitlevel["plane"], bins=36, range=(-0.5,35.5))[0])


        if planeresid is None:
            planeresid = counts,edges
        else:
            planeresid = planeresid[0] + counts, edges


        for name, (bins, xmin, xmax,fitgauss) in allvars.items():
            if name in variables_hitlevel:
                v = df_hitlevel[name]
            else:
                v = df_evt_only[name]

            counts, edges = np.histogram(v, bins=bins, range=(xmin,xmax))

            if name not in plots:
                plots[name] = counts, edges
            else:
                plots[name] = plots[name][0] + counts, edges


        ax = plt.subplot(gs[:4])
        line = hist_meanresid(ax,planeresid)

        for i, (name, (bins, xmin, xmax,fitgauss)) in enumerate(allvars.items()):
            ax = plt.subplot(gs[i+4])
            line = hist_poiserr(ax, plots[name], fitgauss=fitgauss)
            qual = 'tracks'
            if name in variables_hitlevel:
                qual = 'hits'

            ax.set_title('%s (%s)' % (name.replace('_', ' '), qual))
            ax.set_xlabel(name.replace('_', ' '), fontsize=9)

            if filename not in lines:
                lines[filename] = line


    fig.suptitle('Collected Tracks (comparing %d files)' % len(sys.argv[1:]), fontsize=20)
    fig.legend(lines.values(), lines.keys(), "upper right")

    plt.tight_layout()

    #plt.show()

    pdf = matplotlib.backends.backend_pdf.PdfPages("aligntrack_display.pdf")
    for fig_ in range(1, fig.number+1):
        pdf.savefig( fig_ )
    pdf.close()

if __name__=='__main__':
    main()
