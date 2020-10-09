

import matplotlib.pyplot as plt
import io,sys,os
import numpy as np
from pathlib import Path

import mp2prod
from matplotlib.backends.backend_pdf import PdfPages
def multipage(filename, figs=None, dpi=200):
    pp = PdfPages(filename)
    if figs is None:
        figs = [plt.figure(n) for n in plt.get_fignums()]
    for fig in figs:
        fig.savefig(pp, format='pdf')
    pp.close()

def get_alignconsts(filename, tablename):
    iobuf = io.StringIO()
    read_lines = False
    with open(filename, 'r') as f:
        for line in f.readlines():
            if 'TABLE' in line:
                if line.strip().split()[-1] == tablename:
                    read_lines = True
                    continue
                else:
                    read_lines = False
                    continue

            if read_lines:
                iobuf.write(line)
    iobuf.seek(0)

    return np.genfromtxt(iobuf,delimiter=',')

def get_alignconsts_errors(filename, tablename):
    if '_in.txt' in filename: return None
    print ("reading corresponding millepede.res file for errors")
    # we need to look in the same folder for a millepede.res file
    path = Path(filename)

    mpres_path = os.path.join(path.parent, 'millepede.res')

    alignconsts = mp2prod.AlignmentConstants()
    alignconsts.read_mp_file(mpres_path)

    x_shift_errors = []
    y_shift_errors = []
    z_shift_errors = []
    xrot = []
    yrot = []
    zrot = []
    for plane in range(0,36):
        x_shift_errors.append(float(alignconsts.tables['TrkAlignPlane'].errors[plane][0]))
        y_shift_errors.append(float(alignconsts.tables['TrkAlignPlane'].errors[plane][1]))
        z_shift_errors.append(float(alignconsts.tables['TrkAlignPlane'].errors[plane][2]))

        xrot.append(float(alignconsts.tables['TrkAlignPlane'].errors[plane][3]))
        yrot.append(float(alignconsts.tables['TrkAlignPlane'].errors[plane][4]))
        zrot.append(float(alignconsts.tables['TrkAlignPlane'].errors[plane][5]))

    return x_shift_errors, y_shift_errors, z_shift_errors, xrot, yrot, zrot


def pairwise(t):
    it = iter(t)
    return zip(it,it)


for param_col in [0,1,2,3,4,5]:
#    plt.figure(param_col)
    to_plot = []
    col_str = ['x-shift', 'y-shift','z-shift', 'x-rot','y-rot','z-rot'][param_col]
    for label, filename in pairwise(sys.argv[1:]):
        shifts_errors = get_alignconsts_errors(filename, 'TrkAlignPlane')
        shifts = get_alignconsts(filename, 'TrkAlignPlane')
        shifts = shifts[:,1+param_col]*1000 # convert to micrometers
        if shifts_errors is None:
            shifts_errors = np.zeros((36,))
        else:
            shifts_errors = np.array(shifts_errors[param_col])*1000

        print (shifts_errors)
        to_plot += [(label, shifts, shifts_errors)]


    plane_id = np.arange(0,35.5,1)

    # pulls for x-shifts, y-shifts, z-shifts
    with plt.style.context(['science', 'grid','no-latex']):
        f_=plt.figure(param_col+6)
        f_.set_size_inches(10,5)
        for filename, shifts, shifts_errors in to_plot:
            plt.hist(shifts/shifts_errors, 
                bins=np.arange(-20,20,4), 
                range=(-20,20), 
                align='left', 
                histtype='step', 
                linewidth=1.2, 
                label=filename) 

        plt.xlim(-20,20)
        plt.ylim(0, 40)
        plt.yticks(np.arange(0,40,4))
        plt.xlabel('%s / error (all planes)' % col_str)
        plt.ylabel('count / 4.0')

        plt.title('%s pull histogram (all planes) after alignment fit' % col_str)
        plt.legend(fontsize='small')


    # shifts per plane and pulls underneath
        f_, axs = plt.subplots(2, num=param_col, gridspec_kw={'height_ratios': [10, 5]})
        f_.set_size_inches(10,5)
        maxmag = 0
        for filename, shifts, shifts_errors in to_plot:
            #plt.scatter(plane_id, shifts, label=filename)
            axs[0].errorbar(plane_id, shifts, yerr=shifts_errors,label=filename,
                        marker = '.',
                        fmt='x-',
                        linewidth=0.5,
                        drawstyle = 'steps-mid')

            maxmagt = np.max(np.abs(shifts))*1.15

            if maxmagt > maxmag:
                maxmag = maxmagt
        axs[0].set_ylim(-maxmag,maxmag)
        axs[0].axhline(0,0,36, color='k',label='Nominal')
        axs[0].set_ylabel('%s (um)' % col_str)
        if 'rot' in col_str:
            axs[0].set_ylabel('%s (mrad)' % col_str)
        axs[0].set_title('%s after alignment fit' % col_str)
        axs[0].legend(fontsize='small')
        maxmag = 10
        for filename, shifts, shifts_errors in to_plot:
            #plt.scatter(plane_id, shifts, label=filename)
            pulls = shifts/shifts_errors
            axs[1].bar(np.arange(0,35.5,1), pulls,label=filename)
            # axs[1].errorbar(np.arange(0,35.5,1), pulls, yerr=np.arange(0,35.5,1)*0,
            #             label=filename,
            #             marker = '.',
            #             fmt='x-',
            #             linewidth=0.5,
            #             drawstyle = 'steps-mid')
            try:
                pullsnonan = pulls[np.isfinite(pulls)]

                tmaxmag = np.max(np.abs(pullsnonan))*1.05
                if tmaxmag > maxmag:
                    maxmag = maxmag
                axs[1].set_ylim(-maxmag,maxmag)
            except:
                axs[1].set_ylim(-10,10)

        axs[1].set_xlabel('Plane ID')
        axs[1].set_ylabel('%s / error' % col_str)


        axs[1].set_xticks(np.arange(0,36,2))
        axs[0].set_xticks(np.arange(0,36,2))


multipage('alignment_multipage.pdf')

plt.show()

