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


def files_processed(processed_files, total_files=1023, \
                    production_pot=1e19, target_pot=2.5e19):
    return target_pot/((processed_files*production_pot)/total_files)



def charged_pion_threshold(d, threshold, scale_factor):
    cp_count={} # charged pion count above tracking threshold per neutrino interaction
    for key in d.keys():
        if d[key]['pdg']==abs(211):
            contained_length=d[key]['contained_length']
            if contained_length>threshold:
                spill_id = key.split('-')[0]
                vertex_id = key.split('-')[1]
                temp=(spill_id, vertex_id)
                if temp not in cp_count: cp_count[temp]=0.
                cp_count[temp]+=1

    cp_length={} # single track MIP length per neutrino interaction
    for key in d.keys():
        spill_id = key.split('-')[0]
        vertex_id = key.split('-')[1]
        temp=(spill_id, vertex_id)
        if temp in cp_count.keys():
            if cp_count[temp]!=1: continue
            contained_length=d[key]['contained_length']
            cp_length[temp]=contained_length            

    count_bgd=0
    for key in cp_length.keys():
        if cp_length[key]>threshold: count_bgd+=1
    print(count_bgd,' irreducible NC backgrounds with ',threshold,' cm tracking threshold')

    bins=np.linspace(0,200,41)
    fig, ax = plt.subplots(figsize=(6,6))

    meson_length=[cp_length[key] for key in cp_length.keys()]
    weight=[scale_factor]*len(meson_length)
    ax.hist(meson_length, bins=bins, weights=weight, histtype='step', color='b')
    ax.axvline(x=threshold, linestyle='dashed', color='orange', label='3 cm tracking threshold')
    ax.set_xlim(0,200)
    ax.legend(loc='lower right')
    ax.set_xlabel(r'Maximum 2x2-contained $\pi^{\pm}$ Track Length [cm]')
    ax.set_ylabel(r'$\nu$ Interactions / 5 cm')
    
    axtwin = ax.twinx()
    axtwin.hist(meson_length, bins=bins, weights=weight, histtype='step',\
                cumulative=True, density=True, linestyle='solid', color='r')
    axtwin.set_ylabel('Cumulative Probability')

    ax.text(115,2250,r'NuMI ME RHC 2.5$\times$10$^{19}$ POT')
    ax.text(0,2250,r'$\nu$NC')

    plt.savefig('irreducible_nc_beam_bgd.png')

    

def main(nc_json_file, tracking_threshold, n_files_processed):
    f = open(nc_json_file)
    nc_pion_dict=json.load(f)
    scale_factor = files_processed(n_files_processed)
    charged_pion_threshold(nc_pion_dict, tracking_threshold, scale_factor)
    


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-nc', '--nc_json_file', default=None, required=False, type=str, \
                        help='''string corresponding to the path of the NC pion backgrounds JSON file''')
    parser.add_argument('-t', '--tracking_threshold', default=3., required=False, type=float, \
                        help='''Tracking threshold in track length [cm]''')
    parser.add_argument('-n', '--n_files_processed', default=1, required=True, type=int, \
                        help='''File count of number of files processed in production sample''')
    args = parser.parse_args()
    main(**vars(args))
