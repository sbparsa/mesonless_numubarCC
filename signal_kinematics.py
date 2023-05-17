import matplotlib
import matplotlib.pyplot as plt
import h5py
import glob
import json
import argparse
import numpy as np
import twoBytwo_defs
import auxiliary
import signal_characterization_and_plotting as sig_char_plot

nu_signal_pdg=-14
pion_pdg={111,211,-211}


'''TO DO:
(1) Dump muon dictionary into json file
'''


def signal_nu_pdg(ghdr, vert_id):
    ghdr_vert_mask = ghdr['vertexID']==vert_id
    ghdr_nu_interaction = ghdr[ghdr_vert_mask]['nu_pdg']
    if ghdr_nu_interaction[0]==nu_signal_pdg: return True
    else: return False
    


def signal_cc(ghdr, vert_id):
    ghdr_vert_mask = ghdr['vertexID']==vert_id
    return ghdr[ghdr_vert_mask]['isCC'][0]

    

def signal_pion_status(gstack, vert_id):
    gstack_vert_mask = gstack['vertexID']==vert_id
    gstack_pdg_set = set(gstack[gstack_vert_mask]['part_pdg'])
    #print("Event PDG Stack:", gstack_pdg_set)
    if len(pion_pdg.intersection(gstack_pdg_set))==0: return True
    else: return False



def main(sim_dir, input_type):

    test_count = 0

    muon_dict = dict() # Initialize muon dictionary
    hadron_dict = dict() # Initialize hadron dictionary
    
    file_ext = '' ## Changes based on input type

    if input_type == 'larnd': 
        file_ext = '.LARNDSIM.h5'
    elif input_type == 'edep':
        file_ext = '.EDEPSIM.h5'

    for sim_file in glob.glob(sim_dir+'/*'+file_ext):

        if test_count ==20 : break
        test_count+=1

        if (test_count % 5 == 0):
            print("Processing file: ", str(test_count))

        sim_h5 = h5py.File(sim_file,'r')

        ### partition file by spill
        unique_spill = np.unique(sim_h5['trajectories']['eventID'])
        for spill_id in unique_spill:

            ghdr, gstack, traj, vert, seg = auxiliary.get_spill_data(sim_h5, spill_id)

            ### partition by vertex ID within beam spill
            for v_i in range(len(vert['vertexID'])):

                vert_pos= [vert['x_vert'][v_i], vert['y_vert'][v_i], vert['z_vert'][v_i]]
                vert_in_active_LAr = twoBytwo_defs.fiducialized_vertex( vert_pos )

                ##### REQUIRE: neutrino vertex in LAr active volume #####
                if vert_in_active_LAr==False: continue

                #print('Vertex in active LAr.')

                vert_id = vert['vertexID'][v_i]

                nu_mu_bar = signal_nu_pdg(ghdr, vert_id)
                is_cc = signal_cc(ghdr, vert_id)
                pionless = signal_pion_status(gstack, vert_id)
                fv_particle_origin=twoBytwo_defs.fiducialized_particle_origin(traj, vert_id)

                ##### REQUIRE: (A) nu_mu_bar, (B) CC, (C) NO pions present, (D) final state particle start point in FV
                if nu_mu_bar==True and is_cc==True and pionless==True and fv_particle_origin==True:
                    sig_char_plot.muon_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, muon_dict)
                    sig_char_plot.hadron_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, hadron_dict)
                    ### TO DO: Hadron dict

                
    sig_char_plot.plot_muons(muon_dict)
    sig_char_plot.plot_hadrons(hadron_dict)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--sim_dir', default=None, required=True, type=str, help='''string corresponding to the path of the directory containing edep-sim or larnd ouput simulation file to be considered''')
    #parser.add_argument('-f', '--sim_file', default=None, required=False, type=str, help='''string corresponding to the path of the edep-sim ouput simulation file to be considered''')
    parser.add_argument('-t', '--input_type', default='larnd', choices=['edep', 'larnd'], type=str, help='''string corresponding to the output file type: edep or larnd''')
    args = parser.parse_args()
    main(**vars(args))
