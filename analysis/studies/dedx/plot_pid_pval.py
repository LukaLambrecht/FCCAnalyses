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
      'JetsConstituents_dEdx_wires_type',

      'JetsConstituents_PID_pval_pads_ele',
      'JetsConstituents_PID_pval_pads_mu',
      'JetsConstituents_PID_pval_pads_pi',
      'JetsConstituents_PID_pval_pads_kaon',
      'JetsConstituents_PID_pval_pads_proton',

      'JetsConstituents_PID_pval_wires_ele',
      'JetsConstituents_PID_pval_wires_mu',
      'JetsConstituents_PID_pval_wires_pi',
      'JetsConstituents_PID_pval_wires_kaon',
      'JetsConstituents_PID_pval_wires_proton',
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

    # loop over pads and wires
    for system in ['pads', 'wires']:
        print(f'Running on system {system}...')

        # make good measurement mask
        measurement_mask = (events[f'JetsConstituents_dEdx_{system}_type'] == 0)
        total_mask = ((kinematic_mask) & (charge_mask) & (measurement_mask))

        # loop over variables
        for pidtype in ['ele', 'mu', 'pi', 'kaon', 'proton']:
            varname = f'JetsConstituents_PID_pval_{system}_{pidtype}'

            # loop over categories
            category_data = {}
            for category_label, category_mask in categories.items():
        
                # get data
                values = events[varname][((total_mask) & (category_mask))]
            
                # strategies for flattening per-constituent data
                strategy = 'flatten' # choose from "flatten" or "leading"

                # approach 1: take all constituents
                if strategy=='flatten':
                    values = ak.flatten(values, axis=None)

                else: raise Exception(f'Strategy {strategy} not recognized.')

                # parsing
                values = values.to_numpy()
                if np.isnan(values).any():
                    msg = 'WARNING: replacing NaN by 0...'
                    print(msg)
                    np.nan_to_num(values, copy=False, nan=0)

                # add to struct
                category_data[category_label] = values

            # make figure
            fig, ax = plt.subplots()
            bins = np.linspace(-1, 1, num=101)
            for category_label, data in category_data.items():
                hist = np.histogram(data, bins=bins)[0]
                errors = np.sqrt(hist)
                binwidths = bins[1:] - bins[:-1]
                integral = np.sum(np.multiply(hist, binwidths))
                if integral > 0:
                    hist = hist / integral
                    errors = errors / integral
                temp = ax.stairs(hist+errors, baseline=hist-errors, edges=bins, fill=True, alpha=0.2)
                color = temp.get_facecolor()
                ax.stairs(hist, edges=bins, linewidth=3, label=category_label, color=color, alpha=1)

            # plot aesthetics
            text = ax.text(0.05, 0.9, f'{system}\n{pidtype} PID p-value', fontsize=12, transform=ax.transAxes)
            text.set_bbox(dict(facecolor='white', alpha=0.5))
            ax.set_ylabel('Jet constituents (normalized)', fontsize=12)
            ax.set_xlabel('Signed p-value', fontsize=12)
            ax.grid(which='both', axis='both')
            leg = ax.legend(fontsize=12)
            ax.set_yscale('log')

            # save figure
            fig.tight_layout()
            outputfile = os.path.join(outputdir, f'pid_pval_{system}_{pidtype}.png')
            fig.savefig(outputfile)
            plt.close()

            # make another figure with absolute value and log scale
            fig, ax = plt.subplots()
            bins = np.logspace(-5, 0, num=101)
            for category_label, data in category_data.items():
                hist = np.histogram(np.abs(data), bins=bins)[0]
                errors = np.sqrt(hist)
                binwidths = bins[1:] - bins[:-1]
                integral = np.sum(np.multiply(hist, binwidths))
                if integral > 0:
                    hist = hist / integral
                    errors = errors / integral
                temp = ax.stairs(hist+errors, baseline=hist-errors, edges=bins, fill=True, alpha=0.2)
                color = temp.get_facecolor()
                ax.stairs(hist, edges=bins, linewidth=3, label=category_label, color=color, alpha=1)

            # plot aesthetics
            text = ax.text(0.05, 0.9, f'{system}\n{pidtype} PID p-value', fontsize=12, transform=ax.transAxes)
            text.set_bbox(dict(facecolor='white', alpha=0.5))
            ax.set_ylabel('Jet constituents (normalized)', fontsize=12)
            ax.set_xlabel('Absolute p-value', fontsize=12)
            ax.grid(which='both', axis='both')
            leg = ax.legend(fontsize=12)
            ax.set_yscale('log')
            ax.set_xscale('log')

            # save figure
            fig.tight_layout()
            outputfile = os.path.join(outputdir, f'pid_pval_{system}_{pidtype}_abs.png')
            fig.savefig(outputfile)
            plt.close()
