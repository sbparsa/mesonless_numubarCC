import twoBytwo_defs
import numpy as np
import json

##### HDF5 FILE PARSING-------------------------------------                    


def print_keys_attributes(sim_h5):
    print(sim_h5.keys(),'\n')
    print('GENIE HDR: ',sim_h5['genie_hdr'].dtype,'\n')
    print('GENIE STACK: ',sim_h5['genie_stack'].dtype,'\n')
    print('SEGMENTS: ', sim_h5['segments'].dtype,'\n')
    print('TRAJECTORIES', sim_h5['trajectories'].dtype,'\n')
    print('VERTICES', sim_h5['vertices'].dtype)



def get_spill_data(sim_h5, spill_id):
    ghdr_spill_mask = sim_h5['genie_hdr'][:]['spillID']==spill_id
    gstack_spill_mask = sim_h5['genie_stack'][:]['spillID']==spill_id
    traj_spill_mask = sim_h5['trajectories'][:]['spillID']==spill_id
    vert_spill_mask = sim_h5['vertices'][:]['spillID']==spill_id
    seg_spill_mask = sim_h5['segments'][:]['spillID']==spill_id

    ghdr = sim_h5['genie_hdr'][ghdr_spill_mask]
    gstack = sim_h5['genie_stack'][gstack_spill_mask]
    traj = sim_h5['trajectories'][traj_spill_mask]
    vert = sim_h5['vertices'][vert_spill_mask]
    seg = sim_h5['segments'][seg_spill_mask]

    return ghdr, gstack, traj, vert, seg
