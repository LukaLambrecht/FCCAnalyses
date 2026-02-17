# Plot primary vertex coordinates
# and correlation between gen and reco primary vertex.


import os
import sys
import uproot
import argparse
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt


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
    namedict = {}
    for varname in branches_to_read:
        coord = varname.split('_')[-1]
        descr = varname
        if varname.startswith('PV_'): descr = f'Reconstructed PV {coord}-coordinate [cm]'
        elif varname.startswith('GenPV_'): descr = f'Generated PV {coord}-coordinate [cm]'
        namedict[varname] = descr

    # read input files
    eventdict = {}
    inputdict = {'baseline': args.baseline, 'test': args.test}
    for name, inputfiles in inputdict.items():
        print(f'Now reading {name}...')
        batches = []
        for idx, inputfile in enumerate(inputfiles):
            print(f'Reading file {idx+1} / {len(inputfiles)}', end='\r')
            readstr = ':'.join([inputfile, treename])
            with uproot.open(readstr) as f:
                batches.append(f.arrays(branches_to_read))
        eventdict[name] = ak.concatenate(batches)

    # define categories
    categories = {}
    category_names = []
    for name, events in eventdict.items():
        this_categories = {
            'all': np.ones(len(events)).astype(bool),
            'b': (events['genEventType'] == 5)
        }
        category_names = list(this_categories.keys())
        categories[name] = this_categories

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

        # loop over categories
        for category_name in category_names:

            # get data
            data = {}
            for name, events in eventdict.items():
                mask = categories[name][category_name]
                xdata = events[matching_variable].to_numpy()[mask]
                ydata = events[variable].to_numpy()[mask]
                data[name] = ydata - xdata
            alldata = np.concatenate(list(data.values()))

            # determine binning
            xmin = np.quantile(alldata, 0.01)
            xmax = np.quantile(alldata, 0.99)
            bins = np.linspace(xmin, xmax, num=51)

            # make plot
            fig, ax = plt.subplots()
            for name, values in data.items():
                ax.hist(values, bins=bins, density=True, histtype='step', linewidth=2, label=name)
            ax.set_ylabel('Events (normalized)', fontsize=12)
            ax.set_xlabel(namedict[variable] + ' residual', fontsize=12)
            ax.legend()
            txt = ''
            for name, values in data.items():
                txt += name + r': $\sigma$ = ' + '{:.3e}\n'.format(np.std(values))
            txt = txt.strip(' \t\n')
            text = ax.text(0.05, 0.95, txt, va='top', transform=ax.transAxes)
            text.set_bbox({'facecolor': 'white', 'alpha': 0.5})
            figname = os.path.join(outputdir, f'{variable}_residual_{category_name}.png')
            fig.savefig(figname)
