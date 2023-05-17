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


def charged_pion_threshold(d, threshold):
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

    lbins=np.linspace(0,200,51); sbins=np.linspace(0,10,51)
    fig, ax = plt.subplots(figsize=(8,8))

    ax.hist([cp_length[key] for key in cp_length.keys()], bins=lbins, histtype='step', color='b')
    ax.axvline(x=threshold, linestyle='dashed', color='orange', label='tracking threshold')
    ax.set_xlim(0,200)
    axtwin = ax.twinx()
    axtwin.hist([cp_length[key] for key in cp_length.keys()], bins=lbins, histtype='step',\
                cumulative=True, density=True, linestyle='solid', color='r')

    axtwin.set_ylabel('Cumulative Probability')
    ax.legend(loc='upper right')
    ax.set_xlabel('Length [cm]')
    ax.set_ylabel(r'$\nu$ Interactions')
    ax.set_title('Maximum Charged Pion Contained Length'+'\n'+\
                 r'Per $\nu$ NC Interaction')
    plt.savefig('irreducible_nc_beam_bgd.png')

    

def main(nc_json_file, tracking_threshold):
    f = open(nc_json_file)
    nc_pion_dict=json.load(f)
    charged_pion_threshold(nc_pion_dict, tracking_threshold)
    


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-nc', '--nc_json_file', default=None, required=False, type=str, \
                        help='''string corresponding to the path of the NC pion backgrounds JSON file''')
    parser.add_argument('-t', '--tracking_threshold', default=3., required=False, type=float, \
                        help='''Tracking threshold in track length [cm]''')
    args = parser.parse_args()
    main(**vars(args))
