import matplotlib
import matplotlib.pyplot as plt
import h5py
import argparse
import numpy as np
import twoBytwo_defs

nu_signal_pdg=-14
pion_pdg={111,211,-211}



def print_keys_attributes(sim_h5):
    print(sim_h5.keys())
    print('GENIE HDR: ',sim_h5['genie_hdr'].dtype)
    print('GENIE STACK: ',sim_h5['genie_stack'].dtype)
    print('SEGMENTS: ', sim_h5['segments'].dtype)
    print('TRAJECTORIES', sim_h5['trajectories'].dtype)
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



def fiducialized_vertex(vert_pos):
    flag=False; x_drift_flag=False; y_vertical_flag=False; z_beam_flag=False
    for i in range(3):
        for i_bounds, bounds in enumerate(twoBytwo_defs.tpc_bounds(i)):
            if vert_pos[i]>bounds[0] and vert_pos[i]<bounds[1]:
                if i==0: x_drift_flag=True; break
                if i==1: y_vertical_flag=True
                if i==2: z_beam_flag=True
    if x_drift_flag==True and y_vertical_flag==True and z_beam_flag==True: flag=True
    return flag



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
    print(gstack_pdg_set)
    if len(pion_pdg.intersection(gstack_pdg_set))==0: return True
    else: return False



def main(sim_file, input_type):

    sim_h5 = h5py.File(sim_file,'r')
#    print_keys_attributes(sim_h5)
#    return

    ### partition file by spill
    unique_spill = np.unique(sim_h5['trajectories']['spillID'])
    for spill_id in unique_spill:

        ghdr, gstack, traj, vert, seg = get_spill_data(sim_h5, spill_id)

        ### partition by vertex ID within beam spill
        for v_i in range(len(vert['vertexID'])):
            vert_pos= [vert['x_vert'][v_i], vert['y_vert'][v_i], vert['z_vert'][v_i]]
            vert_in_active_LAr = fiducialized_vertex( vert_pos )

            nu_mu_bar = signal_nu_pdg(ghdr, vert['vertexID'][v_i])
            is_cc = signal_cc(ghdr, vert['vertexID'][v_i])
            pionless = signal_pion_status(gstack, vert['vertexID'][v_i])
            


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--sim_file', default=None, required=True, type=str, help='''string corresponding to the path of the edep-sim ouput simulation file to be considered''')
    parser.add_argument('-t', '--input_type', default='edep', choices=['edep', 'larnd'], type=str, help='''string corresponding to the output file type: edep or larnd''')
    args = parser.parse_args()
    main(**vars(args))
