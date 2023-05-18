import matplotlib
import matplotlib.pyplot as plt
import h5py
import argparse
import numpy as np
import twoBytwo_defs
import threshold_backgrounds
import auxiliary
import glob
import json
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)
# threshold backgrounds: ==> only addressing charged pions for now...!!!
# (1) any charged pions in event are sufficiently short such that they are undetectable --> less than 3 cm (6-7 pixels)
# (2) one pi0 gamma escapes FV
# (3) pi0 gammas inbalance: gamma1>>gamma2
# (4) pi0 gammas are colinear

def files_processed(processed_files, total_files=1023, \
                    production_pot=1e19, target_pot=2.5e19):
    return target_pot/((processed_files*production_pot)/total_files)
    


def charged_pion_threshold(d, threshold, scale_factor):
    cp_length={} # longest charged pion per neutrino interaction
    for key in d.keys():
        if d[key]['pdg']==abs(211):
            contained_length=d[key]['contained_length']
            if contained_length>0.:
                spill_id = key.split('-')[0]
                vertex_id = key.split('-')[1]
                temp=(spill_id, vertex_id)
                if temp not in cp_length: cp_length[temp]=0.
                if contained_length>cp_length[temp]:
                    cp_length[temp]=contained_length

    count_bgd=0
    for key in cp_length.keys():
        if cp_length[key]<=threshold: count_bgd+=1
    print(count_bgd,' irreducible CC backgrounds with ',threshold,' cm tracking threshold')

    lbins=np.linspace(0,200,41); sbins=np.linspace(0,10,11)
    fig, ax = plt.subplots(figsize=(6,6))

    meson_length=[cp_length[key] for key in cp_length.keys()]
    weight=[scale_factor]*len(meson_length)

    ax.hist(meson_length, bins=lbins, weights=weight, histtype='step', color='b')
    ax.axvline(x=threshold, linestyle='dashed', color='orange', label='3 cm tracking threshold')
    ax.set_xlim(0,200)
    ax.legend(loc='lower right')
    ax.set_xlabel(r'Maximum 2x2-contained $\pi^{\pm}$ Track Length [cm]')
    ax.set_ylabel(r'$\nu$ Interactions / 5 cm')
    
    axtwin = ax.twinx()
    axtwin.hist([cp_length[key] for key in cp_length.keys()], bins=lbins, histtype='step',\
                cumulative=True, density=True, linestyle='solid', color='r')
    axtwin.set_ylabel('Cumulative Probability')
    
    ax1 = plt.axes([0,0,1,1])
    ip =InsetPosition(ax, [0.4, 0.3, 0.5, 0.5])
    ax1.set_axes_locator(ip)
    mark_inset(ax, ax1, loc1=2, loc2=4, fc='none', ec='0.5')
    
    ax1.hist(meson_length, bins=sbins, weights=weight, histtype='step', color='b')
    ax1.set_xlim(0,10)
    ax1.axvline(x=threshold, linestyle='dashed', color='orange', label='tracking threshold')
#    ax1.set_xlabel(r'$\pi^{\pm}$ Length [cm]')
#    ax1.set_ylabel(r'$\nu$ Interactions')

    ax.text(115,2420,r'NuMI ME RHC 2.5$\times$10$^{19}$ POT')
    ax.text(0,2420,r'$\nu$CC')
    plt.savefig('irreducible_cc_beam_bgd.png')


def main(cc_json_file, tracking_threshold, n_files_processed):
    f = open(cc_json_file)
    cc_pion_dict=json.load(f)
    scale_factor = files_processed(n_files_processed)
    charged_pion_threshold(cc_pion_dict, tracking_threshold, scale_factor)
    


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-cc', '--cc_json_file', default=None, required=False, type=str, \
                        help='''string corresponding to the path of the CC pion backgrounds JSON file''')
    parser.add_argument('-t', '--tracking_threshold', default=3., required=False, type=float, \
                        help='''Tracking threshold in track length [cm]''')
    parser.add_argument('-n', '--n_files_processed', default=1, required=True, type=int, \
                        help='''File count of number of files processed in production sample''')
    args = parser.parse_args()
    main(**vars(args))
