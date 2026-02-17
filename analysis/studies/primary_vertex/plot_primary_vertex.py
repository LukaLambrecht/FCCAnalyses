# Plot primary vertex coordinates
# and correlation between gen and reco primary vertex.


import os
import sys
import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt


if __name__=='__main__':

    # settings
    inputfiles = sys.argv[1:]
    treename = 'events'
    outputdir = 'output_plots'
    branches_to_read = [
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
    batches = []
    for idx, inputfile in enumerate(inputfiles):
        print(f'Reading file {idx+1} / {len(inputfiles)}', end='\r')
        readstr = ':'.join([inputfile, treename])
        with uproot.open(readstr) as f:
            batches.append(f.arrays(branches_to_read))
    events = ak.concatenate(batches)

    # make output directory
    if not os.path.exists(outputdir): os.makedirs(outputdir)

    # loop over individual variables to plot
    for variable in events.fields:
        print(f'Running on variable {variable}...')

        # get data
        data = events[variable].to_numpy()
        mean = np.mean(data)
        std = np.std(data)

        # determine suitable binning
        minv = np.quantile(data, 0.01)
        maxv = np.quantile(data, 0.99)
        bins = np.linspace(minv, maxv, num=51)

        # make figure
        fig, ax = plt.subplots()
        ax.hist(data, bins=bins, density=True, histtype='step', linewidth=2, label=variable)

        # plot aesthetics
        ax.set_ylabel('Events (normalized)', fontsize=12)
        ax.set_xlabel(namedict[variable], fontsize=12)
        infotxt = 'Avg: {:.3e}'.format(mean) + '\n' + 'Std: {:.3e}'.format(std)
        ax.text(0.05, 0.95, infotxt, fontsize=12, va='top', transform=ax.transAxes)
        #ax.legend()

        # save figure
        fig.tight_layout()
        figname = os.path.join(variable+'.png')
        fig.savefig(figname)

        # same with log scale
        ax.set_yscale('log')
        fig.tight_layout()
        figname = os.path.join(variable+'.png')
        fig.savefig(figname)

        # close figures to save memory
        plt.close()

    # loop over correlations to plot
    for variable in events.fields:
        if 'Gen' in variable: continue
        matching_variable = 'Gen'+variable
        if matching_variable not in events.fields:
            msg = f'Expected variable {matching_variable} not found.'
            raise Exception(msg)

        # make scatter plot
        xdata = events[matching_variable].to_numpy()
        ydata = events[variable].to_numpy()
        xmin = np.quantile(xdata, 0.01)
        xmax = np.quantile(xdata, 0.99)
        ymin = np.quantile(ydata, 0.01)
        ymax = np.quantile(ydata, 0.99)
        fig, ax = plt.subplots()
        ax.scatter(xdata, ydata, c='b', s=5, label='MC events')
        xpred = np.linspace(xmin, xmax, num=2)
        ypred = xpred
        ax.plot(xpred, ypred, color='r', linestyle='dashed', label='Average expectation (reco = gen)')
        ax.set_xlabel(matching_variable)
        ax.set_ylabel(variable)
        ax.legend(framealpha=0.9)
        ax.set_xlim((xmin, xmax))
        ax.set_ylim((ymin, ymax))
        fig.tight_layout()
        figname = os.path.join(outputdir, variable+'_vs_'+matching_variable+'_scatter.png')
        fig.savefig(figname)

        # make density histogram
        xbins = np.linspace(xmin, xmax, num=51)
        ybins = np.linspace(ymin, ymax, num=51)
        fig, ax = plt.subplots()
        ax.hist2d(xdata, ydata, bins=(xbins, ybins), density=True)
        ax.plot(xpred, ypred, color='r', linestyle='dashed', label='Average expectation (reco = gen)')
        ax.legend(framealpha=0.9, fontsize=12)
        ax.set_xlabel(namedict[matching_variable], fontsize=12)
        ax.set_ylabel(namedict[variable], fontsize=12)
        fig.tight_layout()
        figname = os.path.join(outputdir, variable+'_vs_'+matching_variable+'_density.png')
        fig.savefig(figname)

        # make residual histogram
        data = ydata - xdata
        xmin = np.quantile(data, 0.01)
        xmax = np.quantile(data, 0.99)
        bins = np.linspace(xmin, xmax, num=51)
        fig, ax = plt.subplots()
        ax.hist(data, bins=bins, density=True)
        ax.set_ylabel('Events (normalized)', fontsize=12)
        ax.set_xlabel(namedict[variable] + ' residual', fontsize=12)
        text = ax.text(0.05, 0.95, r'$\sigma$ = ' + '{:.3e}'.format(np.std(data)), va='top', transform=ax.transAxes)
        text.set_bbox({'facecolor': 'white', 'alpha': 0.5})
        figname = os.path.join(outputdir, variable+'_vs_'+matching_variable+'_residual.png')
        fig.savefig(figname)
