# Do model evaluation for jet flavour tagging

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

thisdir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(thisdir)

from tools import read_file
from plot_roc_multi import make_roc_curves

# global pyplot settings
plt.rc("text", usetex=True)
plt.rc("font", family="serif")


if __name__=='__main__':

    # read command line args
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputfiles', required=True, nargs='+')
    args = parser.parse_args()

    # hard coded settings
    treename = 'Events'
    signal_categories = {
        's': {
            'label_branch': 'recojet_isS',
            'score_branch': 'score_recojet_isS',
            'color': 'blue',
            'label': r'$s$-jets'
        }
    }
    background_categories = {
        'ud': {
            'label_branch': 'recojet_isUDG',
            'score_branch': 'score_recojet_isUDG',
            'color': 'green',
            'label': r'$ud$-jets'
        }
    }
    all_categories = {**signal_categories, **background_categories}

    # find all branches to read
    branches_to_read = (
        [cat['label_branch'] for cat in all_categories.values()]
        + [cat['score_branch'] for cat in all_categories.values()]
    )

    # loop over input files
    roc_curves = []
    aucs = []
    for inputfile in args.inputfiles:
        print(f'Running on input file {inputfile}...')

        # load events
        events = read_file(
                   inputfile,
                   treename = treename,
                   branches = branches_to_read)
        keys = list(events.keys())
        nevents = len(events[list(events.keys())[0]])
        print('Read events file with following properties:')
        print(f'  - Keys: {keys}.')
        print(f'  - Number of events: {nevents}.')

        # define output directory
        outputdir = inputfile.replace('.root', '_plots')

        # make roc curve
        roc_curve, auc = make_roc_curves(
                            events,
                            signal_categories,
                            background_categories)
        roc_curves.append(roc_curve)
        aucs.append(auc)

    # make a plot
    labels = ['$s$ vs $ud$, all input features', '$s$ vs $ud$, no dE/dx']
    stylelist = [None, ':']
    fig, ax = plt.subplots()
    for idx in range(len(args.inputfiles)):

        # get roc curve 
        key = ('s', 'ud')
        signal_category_settings = signal_categories['s']
        background_category_settings = background_categories['ud']
        efficiency_sig, efficiency_bkg = roc_curves[idx][key]
        auc = aucs[idx][key]

        # make a plot of the ROC curve
        label = labels[idx]
        label += ' (AUC: {:.2f})'.format(auc)
        color = 'forestgreen'
        ax.plot(efficiency_bkg, efficiency_sig,
        color=color, linewidth=3, label=label, linestyle=stylelist[idx])

    # add random guessing line
    dummy_efficiency = np.linspace(0, 1, num=101)
    ax.plot(dummy_efficiency, dummy_efficiency,
      color='darkblue', linewidth=3, linestyle='--', label='Random guessing')

    # add aleph logo
    docms = True
    extracmstext = 'Archived Sim.'
    if docms:
        cmstext = r'$\bf{ALEPH}$'
        if extracmstext is not None:
            for part in extracmstext.split(' '): cmstext += r' $\it{' + f' {part}' + r'}$'
        cmstext_in_box = True # maybe later add as argument
        if cmstext_in_box:
            text = ax.text(0.02, 0.98, cmstext,
                    ha='left', va='top', fontsize=20, transform=ax.transAxes)
            text.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='white'))
            # modify the axis range to accommodate the CMS text
            #yscale = ax.get_ylim()[1] - ax.get_ylim()[0]
            #ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] + yscale*0.2)
        else:
            ax.text(0., 1., cmstext,
                    ha='left', va='bottom', fontsize=20, transform=ax.transAxes)

    # other plot settings
    ax.set_xlabel('Background pass-through', fontsize=22)
    ax.set_ylabel('Signal efficiency', fontsize=22)
    ax.tick_params(axis='both', labelsize=17)
    ax.grid(which='both')
    ax.set_ylim((-0.05, 1.1))
    legend_in_box = False # maybe later add as argument
    if legend_in_box: leg = ax.legend(fontsize=17)
    else: 
        leg = ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=17)
    # save the figure
    figname = os.path.join(outputdir, 'roc.png')
    fig.savefig(figname, bbox_extra_artists=(leg,), bbox_inches='tight')
    fig.savefig(figname.replace('.png', '.pdf'), bbox_extra_artists=(leg,), bbox_inches='tight')
    print(f'Saved figure {figname}.')
