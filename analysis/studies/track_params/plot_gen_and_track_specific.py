import os
import sys
import json
import uproot
import argparse
import numpy as np
import awkward as ak
from fnmatch import fnmatch
import matplotlib.pyplot as plt

from plot_gen_and_track_generic import get_track_curve
from plot_gen_and_track_generic import get_track_data_from_event
from plot_gen_and_track_generic import get_genpart_data_from_event
    

if __name__=='__main__':

    # settings
    inputfile = sys.argv[1]
    eventidx = int(sys.argv[2])
    sposrange = 5
    snegrange = 0
    xlim = None
    ylim = None
    zlim = None
    xlim = (-3, 3)
    ylim = (-3, 3)
    zlim = (-3, 3)

    # set branches to read
    branches_to_read = [
        'GenParticle_px',
        'GenParticle_py',
        'GenParticle_pz',
        'GenParticle_x',
        'GenParticle_y',
        'GenParticle_z',
        'GenParticle_pdgId',

        'Tracks_d0',
        'Tracks_z0',
        'Tracks_phi0',
        'Tracks_omega',
        'Tracks_tanlambda',

        'TrackToMCMap'
    ]

    # read input files
    treename = 'events'
    readstr = ':'.join([inputfile, treename])
    with uproot.open(readstr) as f:
        events = f.arrays(branches_to_read, entry_start=eventidx, entry_stop=eventidx+1)
    print(f'Read {len(events)} entries.')

    track_data = get_track_data_from_event(events, sposrange=sposrange, snegrange=snegrange)
    genpart_data = get_genpart_data_from_event(events, sposrange=sposrange)
    tracktomc = np.squeeze(events['TrackToMCMap'].to_numpy())

    # find tracks and gen particles according to specific selections
    trackids = []
    genpartids = []
    for track_idx, genpart_idx in enumerate(tracktomc):
        if genpart_idx >= len(genpart_data): continue
        if np.abs(genpart_data[genpart_idx]['genpart_pdgid']) != 211: continue
        print(track_data[track_idx]['refpoint_coords'])
        print(genpart_data[genpart_idx]['refpoint_coords'])
        print('---')
        trackids.append(track_idx)
        genpartids.append(genpart_idx)

    # make a plot
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.axis('equal')

    # add the tracks
    label = 'Track helix'
    for track_idx, track_params in enumerate(track_data):

        if track_idx not in trackids: continue

        color = 'blue'
        x, y, z = track_params['track_coords']
        ax.plot(x, y, z, color=color, label=label)
        x, y, z = track_params['ext_coords']
        ax.plot(x, y, z, color=color)
        label = None

    # add the gen particles
    label = 'Gen particle'
    for genpart_idx, genpart_params in enumerate(genpart_data):
    
        if genpart_idx not in genpartids: continue

        x, y, z = genpart_params['genpart_coords']
        ax.plot(x, y, z, color='red', linestyle='--', label=label)
        label = None

    # add the z-axis / beamline
    zlims = ax.get_zlim()
    ax.plot([0, 0], [0, 0], [zlims[0], zlims[1]], color='black', linestyle='dashed')

    # plot aesthetics
    if xlim is not None: ax.set_xlim(xlim)
    if ylim is not None: ax.set_ylim(xlim)
    if zlim is not None: ax.set_zlim(zlim)
    ax.set_xlabel('x [cm]', fontsize=12)
    ax.set_ylabel('y [cm]', fontsize=12)
    ax.set_zlabel('z [cm]', fontsize=12)
    ax.legend(fontsize=12)
  
    fig.tight_layout()
    fig.savefig('test.png', dpi=300)
