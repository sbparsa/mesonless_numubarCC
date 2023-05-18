import numpy as np
import matplotlib.pyplot as plt
import json
import argparse


def files_processed(processed_files, total_files=1023, \
                    production_pot=1e19, target_pot=2.5e19):
    return target_pot/((processed_files*production_pot)/total_files)


def plot_stacked_histo(signal, signal_factor, wrong_sign, wrong_sign_factor, \
                       metric, bins, xlabel, ylabel, figname, leg_location):
    fig, ax = plt.subplots(figsize=(6,6))

    s = [signal[key][metric] for key in signal.keys()]
    s_weight = [signal_factor]*len(s)
    w = [wrong_sign[key][metric] for key in wrong_sign.keys()]
    w_weight = [wrong_sign_factor]*len(w)

    ax.hist(s, bins=bins, histtype='bar', stacked=True, \
            label=r'mesonless $\bar{\nu}_\mu$ CC')
    ax.hist(w, bins=bins, histtype='bar', stacked=True, \
            label=r'mesonless $\nu_\mu$ CC')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(loc=leg_location)
    ax.grid(True)
    plt.savefig(figname+'.png')
    


def main(signal, n_signal, wrong_sign, n_wrong_sign):
    f = open(signal)
    signal_dict=json.load(f)
    signal_sf = files_processed(50)

    f = open(wrong_sign)
    wrong_sign_dict=json.load(f)
    wrong_sign_sf = files_processed(50)                     

    plot_stacked_histo(signal_dict, signal_sf, \
                       wrong_sign_dict, wrong_sign_sf, \
                       'nu_energy', np.linspace(0,1e4,51), \
                       r'$\nu$ Energy [MeV]', r'$\nu$ Interactions / 200 MeV',\
                       'stacked_nu_energy', 'upper right')

#    plot_stacked_histo(signal_dict, signal_sf, \
#                       wrong_sign_dict, wrong_sign_sf, \
#                       'q2', np.linspace(0,1e4,51), \
#                       r'$Q^2$ [MeV]', r'$\nu$ Interactions / 200 MeV',\
#                       'stacked_q2', 'upper right')

    plot_stacked_histo(signal_dict, signal_sf, \
                       wrong_sign_dict, wrong_sign_sf, \
                       'mom', np.linspace(0,1e4,51), \
                       r'Muon Candidate Momentum [MeV/c]', r'$\nu$ Interactions / 200 MeV/c',\
                       'stacked_mu_momentum', 'upper right')

    plot_stacked_histo(signal_dict, signal_sf, \
                       wrong_sign_dict, wrong_sign_sf, \
                       'ang', np.linspace(0,np.pi,51), \
                       r'$\theta_\mu$ [radians]', r'$\nu$ Interactions',\
                       'stacked_mu_angle', 'upper right')

    plot_stacked_histo(signal_dict, signal_sf, \
                       wrong_sign_dict, wrong_sign_sf, \
                       'vtx_x', np.linspace(-600,600,51), \
                       r'$\nu$ Vertex X Position [cm]', r'$\nu$ Interactions',\
                       'stacked_vertex_x', 'upper right')

    plot_stacked_histo(signal_dict, signal_sf, \
                       wrong_sign_dict, wrong_sign_sf, \
                       'vtx_y', np.linspace(-500,500,51), \
                       r'$\nu$ Vertex Y Position [cm]', r'$\nu$ Interactions',\
                       'stacked_vertex_y', 'upper right')

    plot_stacked_histo(signal_dict, signal_sf, \
                       wrong_sign_dict, wrong_sign_sf, \
                       'vtx_z', np.linspace(-2000,1000,51), \
                       r'$\nu$ Vertex Z Position [cm]', r'$\nu$ Interactions',\
                       'stacked_vertex_z', 'upper right')
    
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s','--signal', default='signal_dict.json', \
                        type=str, help='''signal JSON''')
    parser.add_argument('-ns','--n_signal', default=50, type=int, \
                        help='''number of files processed for signal JSON''')
    parser.add_argument('-w','--wrong_sign', default='wrong_sign_bkg_dict.json', \
                        type=str, help='''wrong sign background JSON''')
    parser.add_argument('-nw','--n_wrong_sign', default=50, type=int, \
                        help='''number of files processed for wrong sign background JSON''')
    args = parser.parse_args()
    main(**vars(args))
