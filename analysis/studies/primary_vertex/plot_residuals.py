# Plot primary vertex coordinates


import os
import sys
import uproot
import argparse
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

# global pyplot settings
plt.rc("text", usetex=True)
plt.rc("font", family="serif")


if __name__=='__main__':

    # get input files
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputfiles', required=True, nargs='+')
    args = parser.parse_args()

    # settings
    treename = 'events'
    outputdir = 'output_plots'
    branches_to_read = [
      'genEventType',
      'PV_x',
      'PV_y',
      'PV_z',
      'GenPV_x',
      'GenPV_y',
      'GenPV_z',
    ]

    # read input files
    batches = []
    for idx, inputfile in enumerate(args.inputfiles):
        print(f'Reading file {idx+1} / {len(args.inputfiles)}', end='\r')
        readstr = ':'.join([inputfile, treename])
        try:
            with uproot.open(readstr) as f:
                batches.append(f.arrays(branches_to_read))
        except:
            msg = f'WARNING: could not read tree {readstr}, skipping...'
            print(msg)
    events = ak.concatenate(batches)

    # define categories
    categories = {
        'ud': ((events['genEventType'] >= 1) & (events['genEventType'] <= 2)),
        's': (events['genEventType'] == 3),
        'c': (events['genEventType'] == 4),
        'b': (events['genEventType'] == 5)
    }
    category_names = list(categories.keys())

    labeldict = {
        'ud': r'$u\overline{u}, d\overline{d}$',
        's': 's',
        'c': 'c',
        'b': 'b',
    }

    colordict = {
        'ud': 'mediumblue',
        's': 'darkorchid',
        'c': 'red',
        'b': 'dodgerblue',
    }

    # make output directory
    if not os.path.exists(outputdir): os.makedirs(outputdir)

    # loop over variables to plot
    fields = events.fields
    for variable in fields:
        if 'Gen' in variable: continue
        matching_variable = 'Gen'+variable
        if matching_variable not in fields:
            continue
            #msg = f'Expected variable {matching_variable} not found.'
            #raise Exception(msg)
        coordinate = variable.split('_')[-1]

        # loop over categories and aggregate data
        data = {}
        for category_name, mask in categories.items():

            # get data
            xdata = events[matching_variable].to_numpy()[mask]
            ydata = events[variable].to_numpy()[mask]
            data[category_name] = ydata - xdata

        # determine binning
        alldata = np.concatenate(list(data.values()))
        #xmin = np.quantile(alldata, 0.01)
        #xmax = np.quantile(alldata, 0.99)
        xmin = -0.02
        xmax = 0.02
        bins = np.linspace(xmin, xmax, num=51)

        # clip outliers
        #for category_name, values in data.items():
        #    data[category_name] = np.clip(values, a_min=xmin, a_max=xmax)

        # make plot
        fig, ax = plt.subplots()
        for category_name, values in data.items():
            # make histogram
            hist = np.histogram(values, bins=bins)[0]
            errors = np.sqrt(hist)
            binwidths = bins[1:] - bins[:-1]
            integral = np.sum(np.multiply(hist, binwidths))
            hist = hist.astype(float) / integral
            errors = errors.astype(float) / integral
            # plot histogram
            ax.stairs(hist+errors, baseline=hist-errors, edges=bins,
                        fill=True, color=colordict[category_name], alpha=0.3)
            ax.stairs(hist, edges=bins, label=labeldict[category_name], 
                        color=colordict[category_name], alpha=1, linewidth=2)
        ax.set_ylabel('Events (normalized)', fontsize=17)
        ax.set_xlabel(f'{coordinate} residual [cm]', fontsize=17)
        ax.legend(fontsize=17)
        # print widths
        txt = f''
        for category_name, values in data.items():
            # remove outliers
            values = values[((values > xmin) & (values < xmax))]
            txt += category_name + r': RMS = ' + '{:.1f} $\mu m$\n'.format(np.std(values)*1e4)
        txt = txt.strip(' \t\n')
        text = ax.text(0.05, 0.95, txt, va='top', transform=ax.transAxes)
        text.set_bbox({'facecolor': 'white', 'alpha': 0.5})
        figname = os.path.join(outputdir, f'{variable}_residual.png')
        fig.savefig(figname)
