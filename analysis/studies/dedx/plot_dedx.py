import os
import sys
import uproot
import awkward as ak
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# global pyplot settings
plt.rc("text", usetex=True)
plt.rc("font", family="serif")

# local imports
from parametrization import get_parametrized_curve


if __name__=='__main__':

    # settings
    inputfiles = sys.argv[1:]
    treename = 'events'
    outputdir = 'output_plots'
    isdata = ('output_data' in inputfiles[0]) # to make more robust

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

    # temp: limit number of events
    #events = events[:3000]

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
    if isdata:
        categories = {
            'data': np.ones(len(events)).astype(bool)
        }

    labeldict = {
        'pion': r'$\pi$',
        'kaon': r'$K$',
        'proton': r'$p$',
        'muon': r'$\mu$',
        'electron': r'$e$',
        'data': 'Data'
    }

    colordict = {
        'pion': 'mediumblue',
        'kaon': 'darkorchid',
        'proton': 'red',
        'muon': 'dodgerblue',
        'electron': 'darkturquoise',
        'data': 'mediumblue'
    }

    # loop over pads and wires
    for system in ['pads', 'wires']:
        print(f'Running on system {system}...')

        # get data
        allvalues = events[f'JetsConstituents_dEdx_{system}_value']
        allerrors = events[f'JetsConstituents_dEdx_{system}_error']
        alltypes = events[f'JetsConstituents_dEdx_{system}_type']

        # make good measurement mask
        measurement_mask = (alltypes == 0)
        
        # optional: add statistical significance to measurement mask
        #significance_mask = (np.divide(allerrors, allvalues) < 1/np.sqrt(50))
        #measurement_mask = ((measurement_mask) & (significance_mask))

        # make total mask to apply
        total_mask = ((kinematic_mask) & (charge_mask) & (measurement_mask))

        # loop over categories
        category_data = {}
        for category_label, category_mask in categories.items():
        
            # get data
            values = allvalues[((total_mask) & (category_mask))]
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

        # print final number of points that will be plotted
        print('Number of measurement points:')
        for category_name, data in category_data.items():
            print(f'  - {category_name}: {len(data[0])}')
        print(f'  -> total: {sum([len(v[0]) for v in category_data.values()])}')

        # make figure
        fig, ax = plt.subplots()
        for category_label, data in category_data.items():
            ax.scatter(data[0], data[1], s=1, c=colordict[category_label], label=labeldict[category_label], alpha=0.1)

        # optional: add parametrizations
        paxis = np.logspace(-1, 2, num=100)
        if isdata:
            for category_label in labeldict.keys():
                if category_label=='data': continue
                y = get_parametrized_curve(paxis, species=category_label, subsystem=system)
                #color = colordict[category_label]
                color = 'black'
                ax.plot(paxis, y, color=color)
        else:
            for category_label in category_data.keys():
                y = get_parametrized_curve(paxis, species=category_label, subsystem=system)
                #color = colordict[category_label]
                color = 'black'
                ax.plot(paxis, y, color=color)

        # plot aesthetics
        ax.set_ylabel('dE/dx (normalized to MIPs)', fontsize=17)
        ax.set_xlabel('Particle momentum [Gev]', fontsize=17)
        ax.grid(which='both', axis='both')
        leg = ax.legend(fontsize=17, loc='upper right')
        for lh in leg.legend_handles:
            lh.set_alpha(1)
            lh.set_sizes([25])
        ax.set_xscale('log')
        ax.set_ylim((0, 5))
        ax.set_xlim((0.3, 60))
        ax.tick_params(labelsize=17)

        docms = True
        if docms:
            cmstext = r'$\bf{ALEPH}$'
            extracmstext = 'Archived Sim.'
            if extracmstext is not None:
                for part in extracmstext.split(' '): cmstext += r' $\it{' + f' {part}' + r'}$'
            cmstext_in_box = False # maybe later add as argument
            if cmstext_in_box:
                ax.text(0.02, 0.98, cmstext,
                    ha='left', va='top', fontsize=20, transform=ax.transAxes)
                # modify the axis range to accommodate the CMS text
                if logscale:
                    yscale = ax.get_ylim()[1]/ax.get_ylim()[0]
                    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*yscale**(0.2))
                else:
                    yscale = ax.get_ylim()[1] - ax.get_ylim()[0]
                    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] + yscale*0.2)
            else:
                ax.text(0., 1., cmstext,
                        ha='left', va='bottom', fontsize=20, transform=ax.transAxes)

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
                # plot histogram
                ax.stairs(hist+errors, baseline=hist-errors, edges=vbins, fill=True, color=colordict[category_label], alpha=0.3)
                ax.stairs(hist, edges=vbins, label=labeldict[category_label], color=colordict[category_label], alpha=1, linewidth=2)

            # plot aesthetics
            ax.set_ylabel('Number of particles (normalized)', fontsize=17)
            ax.set_xlabel('dE/dx (normalized to MIPs)', fontsize=17)
            ax.grid(which='both', axis='both')
            ax.legend(fontsize=17)
            text = ax.text(0.05, 0.95, f'{plow} $<$ p $<$ {phigh} [Gev]', va='top',
                            transform=ax.transAxes, fontsize=17,
                            bbox={'facecolor': 'white', 'alpha': 0.5})
            ax.tick_params(labelsize=17)
            ymin, ymax = ax.get_ylim()
            ax.set_ylim((ymin, ymax*1.3))

            # save figure
            fig.tight_layout()
            outputfile = os.path.join(outputdir, f'dedx_{system}_pslice_{plow}_{phigh}.png')
            fig.savefig(outputfile)
            plt.close()

