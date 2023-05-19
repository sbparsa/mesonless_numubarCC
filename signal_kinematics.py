import matplotlib
import matplotlib.pyplot as plt
import h5py
import glob
import json
import argparse
import numpy as np
import twoBytwo_defs
import auxiliary
import signal_characterization as sig_char
from plot_signal_muons import plot_muons
from plot_signal_hadrons import plot_hadrons

def main(sim_dir, input_type, file_limit):

    test_count = 0

    ### NOTE: Current POT scaling is based on MiniRun3 larnd file situation
    if int(file_limit) < 1000.: 
        scale_factor = (1./(int(file_limit)/1000.))*2.5
    else:
        scale_factor = 2.5

    muon_dict = dict() # Initialize muon dictionary
    hadron_dict = dict() # Initialize hadron dictionary

    ws_muon_dict = dict() # Initialize wrong sign muon dictionary
    ws_hadron_dict = dict() # Initialize wrong sign hadron dictionary
    
    # Dictionaries for combining with other background explorations
    signal_dict = dict() # Initialize dictionary for signal muons for full comparison
    wrong_sign_bkg_dict = dict() # Initialize dictionary for wrong sign bkg muons for full comparison
    
    file_ext = '' ## Changes based on input type

    if input_type == 'larnd': 
        file_ext = '.LARNDSIM.h5'
    elif input_type == 'edep':
        file_ext = '.EDEPSIM.h5'

    for sim_file in glob.glob(sim_dir+'/*'+file_ext):

        if test_count ==int(file_limit) : break
        test_count+=1

        if (test_count % 5 == 0):
            print("Processing file: ", str(test_count))

        sim_h5 = h5py.File(sim_file,'r')

        ### partition file by spill
        unique_spill = np.unique(sim_h5['trajectories']['eventID'])
        for spill_id in unique_spill:

            ghdr, gstack, traj, vert, seg = auxiliary.get_spill_data(sim_h5, spill_id, input_type)

            ### partition by vertex ID within beam spill
            for v_i in range(len(vert['vertexID'])):

                vert_pos= [vert['x_vert'][v_i], vert['y_vert'][v_i], vert['z_vert'][v_i]]
                vert_in_active_LAr = twoBytwo_defs.fiducialized_vertex( vert_pos )

                ##### REQUIRE: neutrino vertex in LAr active volume #####
                if vert_in_active_LAr==False: continue

                #print('Vertex in active LAr.')

                vert_id = vert['vertexID'][v_i]

                nu_mu_bar = auxiliary.signal_nu_pdg(ghdr, vert_id)
                nu_mu = auxiliary.wrong_sign_nu_pdg(ghdr, vert_id)
                is_cc = auxiliary.signal_cc(ghdr, vert_id)
                mesonless = auxiliary.signal_meson_status(gstack, vert_id)
                fv_particle_origin=twoBytwo_defs.fiducialized_particle_origin(traj, vert_id)

                ##### REQUIRE: (A) nu_mu_bar, (B) CC, (C) NO pions present, (D) final state particle start point in FV
                if nu_mu_bar==True and is_cc==True and mesonless==True and fv_particle_origin==True:
                    sig_char.muon_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, muon_dict, wrong_sign=False)
                    sig_char.hadron_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, hadron_dict)
                    sig_char.get_truth_dict(spill_id, vert_id, ghdr, gstack, traj, vert, seg, signal_dict)
                elif nu_mu==True and is_cc==True and mesonless==True and fv_particle_origin==True:
                    sig_char.muon_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, ws_muon_dict, wrong_sign=True)
                    sig_char.hadron_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, ws_hadron_dict)
                    sig_char.get_truth_dict(spill_id, vert_id, ghdr, gstack, traj, vert, seg, wrong_sign_bkg_dict)

    # PLOT: Signal Event Info      
    plot_muons(muon_dict, scale_factor, sig_bkg = 0)
    plot_hadrons(hadron_dict, scale_factor, sig_bkg = 0)

    # PLOT: Wrong Sign Background Event Info   
    plot_muons(ws_muon_dict, scale_factor, sig_bkg = 3)
    plot_hadrons(ws_hadron_dict, scale_factor, sig_bkg = 3)

    auxiliary.save_dict_to_json(signal_dict, "signal_dict", True)
    auxiliary.save_dict_to_json(wrong_sign_bkg_dict, "wrong_sign_bkg_dict", True)
    auxiliary.save_dict_to_json(muon_dict, "muon_dict", True)
    #auxiliary.save_dict_to_json(hadron_dict, "hadron_dict", True) # TODO: Make Hadron Dict saveable

    signal_count = len(signal_dict)*scale_factor
    wrong_sign_bkg_count = len(wrong_sign_bkg_dict)*scale_factor
    outfile = open('signal_wrong_sign_bkg_event_counts.txt', "w")
    outfile.writelines(["Signal Events (scaled to 2.5e19 POT): "+str(signal_count)+"\n", \
                        "Wrong Sign Background Events (scaled to 2.5e19 POT): "+str(wrong_sign_bkg_count)+"\n", \
                        "Number of files used to get count: "+str(file_limit)+"\n"])
    outfile.close()


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--sim_dir', default=None, required=True, type=str, help='''string corresponding to the path of the directory containing edep-sim or larnd ouput simulation file to be considered''')
    parser.add_argument('-l', '--file_limit', default=None, required=True, type=str, help='''int corresponding to the maximum number of edep-sim or larnd ouput simulation files to be considered''')
    #parser.add_argument('-f', '--sim_file', default=None, required=False, type=str, help='''string corresponding to the path of the edep-sim ouput simulation file to be considered''')
    parser.add_argument('-t', '--input_type', default='larnd', choices=['edep', 'larnd'], type=str, help='''string corresponding to the output file type: edep or larnd''')
    args = parser.parse_args()
    main(**vars(args))
