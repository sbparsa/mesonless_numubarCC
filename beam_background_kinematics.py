import matplotlib
import matplotlib.pyplot as plt
import h5py
import argparse
import numpy as np
import twoBytwo_defs
import threshold_backgrounds

nu_signal_pdg=-14
pion_pdg={111,211,-211}

'''TO DO:
(1) Normalize to 2.5x10^19 POT
(2) Check pi0 visible energy containment calculation
(3) Optional multi-file consuming with glob
(4) Plot: pion multiplicity per spill (total, by PDG)
(5) Plot: parent PDG pie chart: add latex labels
'''


def print_keys_attributes(sim_h5):
    print(sim_h5.keys(),'\n')
    print('GENIE HDR: ',sim_h5['genie_hdr'].dtype,'\n')
    print('GENIE STACK: ',sim_h5['genie_stack'].dtype,'\n')
    print('SEGMENTS: ', sim_h5['segments'].dtype,'\n')
    print('TRAJECTORIES', sim_h5['trajectories'].dtype,'\n')
    print('VERTICES', sim_h5['vertices'].dtype)

    

def get_spill_data(sim_h5, spill_id):
    ### mask data if not spill under consideration
    ghdr_spill_mask = sim_h5['genie_hdr'][:]['spillID']==spill_id
    gstack_spill_mask = sim_h5['genie_stack'][:]['spillID']==spill_id
    traj_spill_mask = sim_h5['trajectories'][:]['spillID']==spill_id
    vert_spill_mask = sim_h5['vertices'][:]['spillID']==spill_id
    seg_spill_mask = sim_h5['segments'][:]['spillID']==spill_id

    ### apply spill mask
    ghdr = sim_h5['genie_hdr'][ghdr_spill_mask]
    gstack = sim_h5['genie_stack'][gstack_spill_mask]
    traj = sim_h5['trajectories'][traj_spill_mask]
    vert = sim_h5['vertices'][vert_spill_mask]
    seg = sim_h5['segments'][seg_spill_mask]
    
    return ghdr, gstack, traj, vert, seg



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
    if len(pion_pdg.intersection(gstack_pdg_set))==0: return True
    else: return False


    
def main(sim_file, input_type):

    sim_h5 = h5py.File(sim_file,'r')
    pion_dict = dict()

#    print_keys_attributes(sim_h5)
    
    ### partition file by spill
    unique_spill = np.unique(sim_h5['trajectories']['spillID'])
    for spill_id in unique_spill:

        ghdr, gstack, traj, vert, seg = get_spill_data(sim_h5, spill_id)

        ### partition by vertex ID within beam spill
        for v_i in range(len(vert['vertexID'])):
            vert_pos= [vert['x_vert'][v_i], vert['y_vert'][v_i], vert['z_vert'][v_i]]
            vert_in_active_LAr = twoBytwo_defs.fiducialized_vertex( vert_pos )

            ##### REQUIRE neutrino vertex in LAr active volume #####
            if vert_in_active_LAr==False: continue

            vert_id = vert['vertexID'][v_i]
            
            nu_mu_bar = signal_nu_pdg(ghdr, vert_id)
            is_cc = signal_cc(ghdr, vert_id)
            pionless = signal_pion_status(gstack, vert_id)
            fv_particle_origin=twoBytwo_defs.fiducialized_particle_origin(traj, vert_id)
                        
            ##### THRESHOLD BACKGROUNDS #####
            ##### REQUIRE: (A) nu_mu_bar, (B) CC, (C) pions present, (D) final state particle start point in FV
            if nu_mu_bar==True and is_cc==True and pionless==False and fv_particle_origin==True:
                threshold_backgrounds.pion_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, pion_dict)
                
    threshold_backgrounds.plot_threshold_backgrounds(pion_dict)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--sim_file', default=None, required=True, type=str, help='''string corresponding to the path of the edep-sim ouput simulation file to be considered''')
    parser.add_argument('-t', '--input_type', default='edep', choices=['edep', 'larnd'], type=str, help='''string corresponding to the output file type: edep or larnd''')
    args = parser.parse_args()
    main(**vars(args))
