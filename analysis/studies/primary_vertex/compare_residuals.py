# Plot primary vertex coordinates
# and correlation between gen and reco primary vertex.


import os
import sys
import uproot
import argparse
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy


if __name__=='__main__':

    # get input files
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--baseline', required=True, nargs='+')
    parser.add_argument('-t', '--test', required=True, nargs='+')
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

    # plot aesthetic settings
    labeldict = {}
    for varname in branches_to_read:
        coord = varname.split('_')[-1]
        descr = varname
        if varname.startswith('PV_'): descr = f'Reconstructed PV {coord}-coordinate [cm]'
        elif varname.startswith('GenPV_'): descr = f'Generated PV {coord}-coordinate [cm]'
        labeldict[varname] = descr

    # read input files
    eventdict = {}
    inputdict = {'baseline': args.baseline, 'test': args.test}
    set_names = list(inputdict.keys())
    for set_name, inputfiles in inputdict.items():
        print(f'Now reading {set_name}...')
        batches = []
        for idx, inputfile in enumerate(inputfiles):
            print(f'Reading file {idx+1} / {len(inputfiles)}', end='\r')
            readstr = ':'.join([inputfile, treename])
            try:
                with uproot.open(readstr) as f:
                    batches.append(f.arrays(branches_to_read))
            except:
                msg = f'WARNING: could not read tree {readstr}, skipping...'
                print(msg)
        eventdict[set_name] = ak.concatenate(batches)

    # define categories
    categories = {}
    category_names = []
    for set_name, events in eventdict.items():
        this_categories = {
            'all': np.ones(len(events)).astype(bool),
            'uds': ((events['genEventType'] >= 1) & (events['genEventType'] <= 3)),
            'c': (events['genEventType'] == 4),
            'b': (events['genEventType'] == 5)
        }
        category_names = list(this_categories.keys())
        categories[set_name] = this_categories

    # make output directory
    if not os.path.exists(outputdir): os.makedirs(outputdir)

    # loop over variables to plot
    fields = eventdict['baseline'].fields
    for variable in fields:
        if 'Gen' in variable: continue
        matching_variable = 'Gen'+variable
        if matching_variable not in fields:
            continue
            #msg = f'Expected variable {matching_variable} not found.'
            #raise Exception(msg)

        # loop over sets and categories and aggregate data
        data = {}
        for set_name, events in eventdict.items():
            data[set_name] = {}
            for category_name, mask in categories[set_name].items():

                # get data
                xdata = events[matching_variable].to_numpy()[mask]
                ydata = events[variable].to_numpy()[mask]
                data[set_name][category_name] = ydata - xdata

        # make plots per category
        for category_name in category_names:

            # determine binning
            thisdata = {set_name: data[set_name][category_name] for set_name in set_names}
            alldata = np.concatenate(list(thisdata.values()))
            xmin = np.quantile(alldata, 0.01)
            xmax = np.quantile(alldata, 0.99)
            bins = np.linspace(xmin, xmax, num=51)

            # clip outliers
            for set_name, values in thisdata.items():
                thisdata[set_name] = np.clip(values, a_min=xmin, a_max=xmax)

            # make plot
            fig, ax = plt.subplots()
            for set_name, values in thisdata.items():
                ax.hist(values, bins=bins, density=True, histtype='step', linewidth=2, label=set_name)
            ax.set_ylabel('Events (normalized)', fontsize=12)
            ax.set_xlabel(labeldict[variable] + ' residual', fontsize=12)
            ax.legend()
            txt = f'Flavour category: {category_name}\n'
            for set_name, values in thisdata.items():
                txt += set_name + r': $\sigma$ = ' + '{:.3e}\n'.format(np.std(values))
            txt = txt.strip(' \t\n')
            text = ax.text(0.05, 0.95, txt, va='top', transform=ax.transAxes)
            text.set_bbox({'facecolor': 'white', 'alpha': 0.5})
            figname = os.path.join(outputdir, f'{variable}_residual_{category_name}.png')
            fig.savefig(figname)

        # make plots per category
        for set_name in set_names:

            # determine binning
            thisdata = deepcopy(data[set_name])
            alldata = np.concatenate(list(thisdata.values()))
            xmin = np.quantile(alldata, 0.01)
            xmax = np.quantile(alldata, 0.99)
            bins = np.linspace(xmin, xmax, num=51)

            # clip outliers
            for category_name, values in thisdata.items():
                thisdata[category_name] = np.clip(values, a_min=xmin, a_max=xmax)

            # make plot
            fig, ax = plt.subplots()
            for category_name, values in thisdata.items():
                ax.hist(values, bins=bins, density=True, histtype='step', linewidth=2, label=category_name)
            ax.set_ylabel('Events (normalized)', fontsize=12)
            ax.set_xlabel(labeldict[variable] + ' residual', fontsize=12)
            ax.legend()
            txt = f'Set: {set_name}\n'
            for category_name, values in thisdata.items():
                txt += category_name + r': $\sigma$ = ' + '{:.3e}\n'.format(np.std(values))
            txt = txt.strip(' \t\n')
            text = ax.text(0.05, 0.95, txt, va='top', transform=ax.transAxes)
            text.set_bbox({'facecolor': 'white', 'alpha': 0.5})
            figname = os.path.join(outputdir, f'{variable}_residual_{set_name}.png')
            fig.savefig(figname)

