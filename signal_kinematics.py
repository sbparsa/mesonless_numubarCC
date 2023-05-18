import matplotlib
import matplotlib.pyplot as plt
import h5py
import glob
import json
import argparse
import numpy as np
import twoBytwo_defs
import signal_characterization_and_plotting as sig_char_plot

nu_signal_pdg=-14
meson_pdg={111,211,-211,130,310,311,321,-321,221,331}


'''TO DO:
(1) Dump muon dictionary into json file
'''

def print_keys_attributes(sim_h5):
    print(sim_h5.keys())
    print('GENIE HDR: ',sim_h5['genie_hdr'].dtype)
    print('GENIE STACK: ',sim_h5['genie_stack'].dtype)
    print('TRACKS: ', sim_h5['tracks'].dtype)
    print('TRAJECTORIES', sim_h5['trajectories'].dtype)
    print('VERTICES', sim_h5['vertices'].dtype)


def get_spill_data(sim_h5, spill_id):
    ### mask data if not spill under consideration
    ghdr_spill_mask = sim_h5['genie_hdr'][:]['eventID']==spill_id
    gstack_spill_mask = sim_h5['genie_stack'][:]['eventID']==spill_id
    traj_spill_mask = sim_h5['trajectories'][:]['eventID']==spill_id
    vert_spill_mask = sim_h5['vertices'][:]['eventID']==spill_id
    seg_spill_mask = sim_h5['tracks'][:]['eventID']==spill_id

    ### apply spill mask
    ghdr = sim_h5['genie_hdr'][ghdr_spill_mask]
    gstack = sim_h5['genie_stack'][gstack_spill_mask]
    traj = sim_h5['trajectories'][traj_spill_mask]
    vert = sim_h5['vertices'][vert_spill_mask]
    seg = sim_h5['tracks'][seg_spill_mask]
    
    return ghdr, gstack, traj, vert, seg

    
def signal_nu_pdg(ghdr, vert_id):
    ghdr_vert_mask = ghdr['vertexID']==vert_id
    ghdr_nu_interaction = ghdr[ghdr_vert_mask]['nu_pdg']
    if ghdr_nu_interaction[0]==nu_signal_pdg: return True
    else: return False

def wrong_sign_nu_pdg(ghdr, vert_id):
    ghdr_vert_mask = ghdr['vertexID']==vert_id
    ghdr_nu_interaction = ghdr[ghdr_vert_mask]['nu_pdg']
    if ghdr_nu_interaction[0]== -1*nu_signal_pdg: return True
    else: return False
    

def tuple_key_to_string(d):
    out={}
    for key in d.keys():
        string_key=""
        max_length=len(key)
        for i in range(max_length):
            if i<len(key)-1: string_key+=str(key[i])+"-"
            else: string_key+=str(key[i])
        out[string_key]=d[key]
    return out                    
            


def save_dict_to_json(d, name, if_tuple):
    with open(name+".json", "w") as outfile:
        if if_tuple==True:
            updated_d = tuple_key_to_string(d)
            json.dump(updated_d, outfile, indent=4)
        else:
            json.dump(d, outfile, indent=4)


def signal_cc(ghdr, vert_id):
    ghdr_vert_mask = ghdr['vertexID']==vert_id
    return ghdr[ghdr_vert_mask]['isCC'][0]

    

def signal_meson_status(gstack, vert_id):
    gstack_vert_mask = gstack['vertexID']==vert_id
    gstack_pdg_set = set(gstack[gstack_vert_mask]['part_pdg'])
    #print("Event PDG Stack:", gstack_pdg_set)
    if len(meson_pdg.intersection(gstack_pdg_set))==0: return True
    else: return False


def main(sim_dir, input_type, file_limit):

    test_count = 0

    ### NOTE: Current POT scaling is based on MiniRun3 larnd file situation
    if int(file_limit) < 1000.: 
        scale_factor = (1./(int(file_limit)/1000.))*2.5
    else:
        scale_factor = 2.5

    muon_dict = dict() # Initialize muon dictionary
    hadron_dict = dict() # Initialize hadron dictionary
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

            ghdr, gstack, traj, vert, seg = get_spill_data(sim_h5, spill_id)

            ### partition by vertex ID within beam spill
            for v_i in range(len(vert['vertexID'])):

                vert_pos= [vert['x_vert'][v_i], vert['y_vert'][v_i], vert['z_vert'][v_i]]
                vert_in_active_LAr = twoBytwo_defs.fiducialized_vertex( vert_pos )

                ##### REQUIRE: neutrino vertex in LAr active volume #####
                if vert_in_active_LAr==False: continue

                #print('Vertex in active LAr.')

                vert_id = vert['vertexID'][v_i]

                nu_mu_bar = signal_nu_pdg(ghdr, vert_id)
                nu_mu = wrong_sign_nu_pdg(ghdr, vert_id)
                is_cc = signal_cc(ghdr, vert_id)
                mesonless = signal_meson_status(gstack, vert_id)
                fv_particle_origin=twoBytwo_defs.fiducialized_particle_origin(traj, vert_id)

                ##### REQUIRE: (A) nu_mu_bar, (B) CC, (C) NO pions present, (D) final state particle start point in FV
                if nu_mu_bar==True and is_cc==True and mesonless==True and fv_particle_origin==True:
                    sig_char_plot.muon_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, muon_dict)
                    sig_char_plot.hadron_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, hadron_dict)
                    sig_char_plot.get_truth_dict(spill_id, vert_id, ghdr, gstack, traj, vert, seg, signal_dict)
                elif nu_mu==True and is_cc==True and mesonless==True and fv_particle_origin==True:
                    sig_char_plot.get_truth_dict(spill_id, vert_id, ghdr, gstack, traj, vert, seg, wrong_sign_bkg_dict)

                
    sig_char_plot.plot_muons(muon_dict, scale_factor)
    sig_char_plot.plot_hadrons(hadron_dict, scale_factor)

    save_dict_to_json(signal_dict, "signal_dict", True)
    save_dict_to_json(wrong_sign_bkg_dict, "wrong_sign_bkg_dict", True)



if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--sim_dir', default=None, required=True, type=str, help='''string corresponding to the path of the directory containing edep-sim or larnd ouput simulation file to be considered''')
    parser.add_argument('-l', '--file_limit', default=None, required=True, type=str, help='''int corresponding to the maximum number of edep-sim or larnd ouput simulation files to be considered''')
    #parser.add_argument('-f', '--sim_file', default=None, required=False, type=str, help='''string corresponding to the path of the edep-sim ouput simulation file to be considered''')
    parser.add_argument('-t', '--input_type', default='larnd', choices=['edep', 'larnd'], type=str, help='''string corresponding to the output file type: edep or larnd''')
    args = parser.parse_args()
    main(**vars(args))
