import os
import sys
import uproot
import awkward as ak
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


if __name__=='__main__':

    # settings
    inputfiles = sys.argv[1:]
    treename = 'events'
    outputdir = 'output_plots'

    variables = [
      'JetsConstituents_pt',
      'JetsConstituents_pz',
      'JetsConstituents_charge',
      'JetsConstituents_pdgId',

      'JetsConstituents_dEdx_pads_type',
      'JetsConstituents_dEdx_pads_value',
      'JetsConstituents_dEdx_pads_error',
      'JetsConstituents_dEdx_wires_type',
      'JetsConstituents_dEdx_wires_value',
      'JetsConstituents_dEdx_wires_error',      
    ]

    # make output dir if needed
    if not os.path.exists(outputdir): os.makedirs(outputdir)

    # read input files
    batches = []
    branches_to_read = variables
    for idx, inputfile in enumerate(inputfiles):
        print(f'Reading file {idx+1} / {len(inputfiles)}', end='\r')
        readstr = ':'.join([inputfile, treename])
        with uproot.open(readstr) as f:
            batches.append(f.arrays(branches_to_read))
    events = ak.concatenate(batches)
    print(f'Read {len(events)} entries.')

    # make kinematic mask (optional)
    kinematic_mask = np.ones(len(events)).astype(bool)
    print(f'Made kinematic mask with {np.sum(kinematic_mask)} / {len(kinematic_mask)} entries passing.')

    # make charge mask (optional)
    charge_mask = (events['JetsConstituents_charge'] != 0)

    # make categories
    categories = {
        'pion': np.abs(events['JetsConstituents_pdgId'])==211,
        'kaon': np.abs(events['JetsConstituents_pdgId'])==321,
        'proton': np.abs(events['JetsConstituents_pdgId'])==2212,
        'muon': np.abs(events['JetsConstituents_pdgId'])==13,
        'electron': np.abs(events['JetsConstituents_pdgId'])==11
    }

    labeldict = {
        'pion': r'$\pi$',
        'kaon': r'$K$',
        'proton': r'$p$',
        'muon': r'$\mu$',
        'electron': r'$e$'
    }

    # loop over pads and wires
    for system in ['pads', 'wires']:
        print(f'Running on system {system}...')

        # make good measurement mask
        measurement_mask = (events[f'JetsConstituents_dEdx_{system}_type'] == 0)
        total_mask = ((kinematic_mask) & (charge_mask) & (measurement_mask))

        # loop over categories
        category_data = {}
        for category_label, category_mask in categories.items():
        
            # get data
            values = events[f'JetsConstituents_dEdx_{system}_value'][((total_mask) & (category_mask))]
            pt = events['JetsConstituents_pt'][((total_mask) & (category_mask))]
            pz = events['JetsConstituents_pz'][((total_mask) & (category_mask))]
            p = np.sqrt(np.square(pt) + np.square(pz))
            
            # strategies for flattening per-constituent data
            strategy = 'flatten' # choose from "flatten" or "leading"

            # approach 1: take all constituents
            if strategy=='flatten':
                values = ak.flatten(values, axis=None)
                p = ak.flatten(p, axis=None)

            else: raise Exception(f'Strategy {strategy} not recognized.')

            # parsing
            values = values.to_numpy()
            p = p.to_numpy()
            if np.isnan(values).any():
                msg = 'WARNING: replacing NaN by 0...'
                print(msg)
                np.nan_to_num(values, copy=False, nan=0)

            # optional: ignore dummy values
            mask = (values > 0.1).astype(bool)
            values = values[mask]
            p = p[mask]

            # add to struct
            category_data[category_label] = (p, values)

        # make figure
        fig, ax = plt.subplots()
        for category_label, data in category_data.items():
            ax.scatter(data[0], data[1], s=1, label=labeldict[category_label], alpha=0.1)

        # plot aesthetics
        ax.set_ylabel('dE/dx', fontsize=15)
        ax.set_xlabel('Momentum [Gev]', fontsize=15)
        ax.grid(which='both', axis='both')
        leg = ax.legend(fontsize=15, loc='upper right')
        for lh in leg.legend_handles:
            lh.set_alpha(1)
            lh.set_sizes([5])
        ax.set_xscale('log')
        ax.set_ylim((0, 10))
        ax.set_xlim((0.1, 10))
        ax.tick_params(labelsize=12)

        # save figure
        fig.tight_layout()
        outputfile = os.path.join(outputdir, f'dedx_{system}_scatter.png')
        fig.savefig(outputfile)

        # make 1D histograms in slices of momentum
        pbins = [0, 0.5, 1, 1.5, 2, 3, 4, 5, 10, 20, 30]
        for pidx in range(len(pbins)-1):
            plow = pbins[pidx]
            phigh = pbins[pidx+1]

            # set binning
            vbins = np.linspace(0.5, 2.5, num=51)
            if plow == 0: vbins = np.linspace(0.5, 10, num=51)
            if plow == 0.5: vbins = np.linspace(0.5, 3, num=51)

            # make the figure
            fig, ax = plt.subplots()
            for category_label, data in category_data.items():
                # get values
                p = data[0]
                mask = ((p > plow) & (p < phigh))
                values = data[1][mask]
                # make histogram
                hist = np.histogram(values, bins=vbins)[0]
                errors = np.sqrt(hist)
                binwidths = vbins[1:] - vbins[:-1]
                integral = np.sum(np.multiply(hist, binwidths))
                hist = hist.astype(float) / integral
                errors = errors.astype(float) / integral
                # set color
                temp = ax.stairs(hist+errors, baseline=hist-errors, edges=vbins, fill=True, alpha=0.3)
                color = temp.get_facecolor()
                ax.stairs(hist, edges=vbins, label=category_label, color=color, alpha=1, linewidth=2)

            # plot aesthetics
            ax.set_ylabel('Number of particles (normalized)', fontsize=12)
            ax.set_xlabel('$dE/dx$', fontsize=12)
            ax.grid(which='both', axis='both')
            ax.legend()
            ax.text(0.05, 0.95, f'{plow} < p < {phigh} [Gev]', va='top', transform=ax.transAxes)

            # save figure
            fig.tight_layout()
            outputfile = os.path.join(outputdir, f'dedx_{system}_pslice_{plow}_{phigh}.png')
            fig.savefig(outputfile)
            plt.close()

