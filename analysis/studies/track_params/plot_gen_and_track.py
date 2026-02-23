import os
import sys
import json
import uproot
import argparse
import numpy as np
import awkward as ak
from fnmatch import fnmatch
import matplotlib.pyplot as plt


def get_track_curve(d0, z0, phi0, omega, tanlambda, s):
    
    # parameter parsing
    x0 = -d0 * np.sin(phi0)
    y0 = d0 * np.cos(phi0)

    # make parametrized coordinates
    x = x0 + 1/omega * (np.sin(phi0 + omega*s) - np.sin(phi0))
    y = y0 - 1/omega * (np.cos(phi0 + omega*s) - np.cos(phi0))
    z = z0 + s * tanlambda

    # also get direction vector at the reference point
    dxds = np.cos(phi0)
    dyds = np.sin(phi0)
    dzds = tanlambda
    direction = np.array([dxds, dyds, dzds])
    direction = direction / np.sqrt(np.sum(np.square(direction)))

    return ( (x, y, z), (x0, y0, z0), direction )


if __name__=='__main__':

    # settings
    inputfile = sys.argv[1]
    eventidx = int(sys.argv[2])
    trackids = [int(el.strip(' ')) for el in sys.argv[3].split(',')]
    sposrange = 5
    snegrange = 0
    #xlim = (-0.07, 0.07)
    #ylim = (-0.07, 0.07)
    xlim = (-3, 3)
    ylim = (-3, 3)
    zlim = None
    zlim = (-3, 3)

    # set branches to read
    branches_to_read = [
        'GenParticle_px',
        'GenParticle_py',
        'GenParticle_pz',
        'GenParticle_x',
        'GenParticle_y',
        'GenParticle_z',

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

    # get track parameters
    d0 = np.squeeze(events['Tracks_d0'].to_numpy())
    z0 = np.squeeze(events['Tracks_z0'].to_numpy())
    phi0 = np.squeeze(events['Tracks_phi0'].to_numpy())
    omega = np.squeeze(events['Tracks_omega'].to_numpy())
    tanlambda = np.squeeze(events['Tracks_tanlambda'].to_numpy())

    # make parametric curves for tracks
    track_data = []
    spos = np.linspace(0, sposrange, num=100)
    sneg = np.linspace(-snegrange, 0, num=100)
    for track_idx in range(len(d0)):

        # get parameters for this track
        this_d0 = d0[track_idx]
        this_z0 = z0[track_idx]
        this_phi0 = phi0[track_idx]
        this_omega = omega[track_idx]
        this_tanlambda = tanlambda[track_idx]

        # skip neutral constituents
        if np.abs(this_d0 + 9) < 1e-12: continue

        # get track curve
        res = get_track_curve(this_d0, this_z0, this_phi0, this_omega, this_tanlambda, spos)
        (track_coords, refpoint_coords, refpoint_direction) = res
        ext_coords = get_track_curve(this_d0, this_z0, this_phi0, this_omega, this_tanlambda, sneg)[0]

        # add results to list
        track_data.append({
            'track_idx': track_idx,
            'track_coords': track_coords,
            'ext_coords': ext_coords,
            'refpoint_coords': refpoint_coords,
            'refpoint_direction': refpoint_direction
        })

    # get gen particle parameters
    genpart_x = np.squeeze(events['GenParticle_x'].to_numpy())
    genpart_y = np.squeeze(events['GenParticle_y'].to_numpy())
    genpart_z = np.squeeze(events['GenParticle_z'].to_numpy())
    genpart_px = np.squeeze(events['GenParticle_px'].to_numpy())
    genpart_py = np.squeeze(events['GenParticle_py'].to_numpy())
    genpart_pz = np.squeeze(events['GenParticle_pz'].to_numpy())

    # make parametric curves for gen particles
    genpart_data = []
    for genpart_idx in range(len(genpart_px)):

            # get parameters for this genpart
            this_x = genpart_x[genpart_idx]
            this_y = genpart_y[genpart_idx]
            this_z = genpart_z[genpart_idx]
            this_px = genpart_px[genpart_idx]
            this_py = genpart_py[genpart_idx]
            this_pz = genpart_pz[genpart_idx]

            # get genpart line
            genpart_coords = (
                this_x + spos * this_px,
                this_y + spos * this_py,
                this_z + spos * this_pz
            )
        
            # add results to list
            genpart_data.append({
                'genpart_idx': genpart_idx,
                'genpart_coords': genpart_coords,
            })

    # get links
    tracktomc = np.squeeze(events['TrackToMCMap'].to_numpy())
    genpartids = [tracktomc[trackidx] for trackidx in trackids]

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
