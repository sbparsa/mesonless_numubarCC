import twoBytwo_defs
import numpy as np
import json

neutral_pdg=[111] #, 22] #, 2112] # add K0, rho0, eta0?
meson_pdg={111,211,-211,130,310,311,321,-321,221,331}
nu_signal_pdg=-14

##### HDF5 FILE PARSING-------------------------------------


def print_keys_attributes(sim_h5, input_type):
    print(sim_h5.keys(),'\n')
    if input_type=='edep':
        print('GENIE HDR: ',sim_h5['genie_hdr'].dtype,'\n')
        print('GENIE STACK: ',sim_h5['genie_stack'].dtype,'\n')
        print('SEGMENTS: ', sim_h5['segments'].dtype,'\n')
        print('TRAJECTORIES', sim_h5['trajectories'].dtype,'\n')
        print('VERTICES', sim_h5['vertices'].dtype)
    elif input_type=='larnd':
        print('GENIE HDR: ',sim_h5['genie_hdr'].dtype,'\n')
        print('GENIE STACK: ',sim_h5['genie_stack'].dtype,'\n')
        print('LIGHT DAT: ',sim_h5['light_dat'].dtype,'\n')
        print('LIGHT TRIGGER: ',sim_h5['light_trig'].dtype,'\n')
        print('LIGHT WAVEFORM: ',sim_h5['light_wvfm'].dtype,'\n')
        print('MC PACKETS ASSN: ',sim_h5['mc_packets_assn'].dtype,'\n')
        print('MESSAGES: ',sim_h5['messages'].dtype,'\n')
        print('PACKETS: ',sim_h5['packets'].dtype,'\n')
        print('TRACKS: ',sim_h5['tracks'].dtype,'\n')
        print('TRAJECTORIES', sim_h5['trajectories'].dtype,'\n')
        print('VERTICES', sim_h5['vertices'].dtype)
        
        
    
def get_spill_data(sim_h5, spill_id, input_type):
    parse_var=''; seg_var=''
    if input_type=='edep': parse_var='spillID'; seg_var='segments'
    elif input_type=='larnd': parse_var='eventID'; seg_var='tracks'        

    ghdr_spill_mask = sim_h5['genie_hdr'][:][parse_var]==spill_id
    gstack_spill_mask = sim_h5['genie_stack'][:][parse_var]==spill_id
    traj_spill_mask = sim_h5['trajectories'][:][parse_var]==spill_id
    vert_spill_mask = sim_h5['vertices'][:][parse_var]==spill_id
    seg_spill_mask = sim_h5[seg_var][:][parse_var]==spill_id

    ghdr = sim_h5['genie_hdr'][ghdr_spill_mask]
    gstack = sim_h5['genie_stack'][gstack_spill_mask]
    traj = sim_h5['trajectories'][traj_spill_mask]
    vert = sim_h5['vertices'][vert_spill_mask]
    seg = sim_h5[seg_var][seg_spill_mask]

    return ghdr, gstack, traj, vert, seg
    


##### OUTPUT DICTIONARY TO JSON------------------------------


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


##### SIGNAL BASIC CHARACTERISTICS--------------------------
def signal_nu_pdg(ghdr, vert_id):
    ghdr_vert_mask = ghdr['vertexID']==vert_id
    ghdr_nu_interaction = ghdr[ghdr_vert_mask]['nu_pdg']
    if ghdr_nu_interaction[0]==nu_signal_pdg: return True
    else: return False


def signal_cc(ghdr, vert_id):
    ghdr_vert_mask = ghdr['vertexID']==vert_id
    return ghdr[ghdr_vert_mask]['isCC'][0]



def signal_meson_status(gstack, vert_id):
    gstack_vert_mask = gstack['vertexID']==vert_id
    gstack_pdg_set = set(gstack[gstack_vert_mask]['part_pdg'])
    if len(meson_pdg.intersection(gstack_pdg_set))==0: return True
    else: return False

def wrong_sign_nu_pdg(ghdr, vert_id):
    ghdr_vert_mask = ghdr['vertexID']==vert_id
    ghdr_nu_interaction = ghdr[ghdr_vert_mask]['nu_pdg']
    if ghdr_nu_interaction[0]== -1*nu_signal_pdg: return True
    else: return False


    
##### FIND PARENT PDG --------------------------------------
def find_parent_pdg(parent_id, vertex_id, traj, ghdr):
    if parent_id==-1:
        ghdr_mask = ghdr['vertexID']==vertex_id
        parent_pdg=ghdr[ghdr_mask]['nu_pdg']
    else:
        parent_mask = traj['trackID']==parent_id
        parent_pdg = traj[parent_mask]['pdgId']
    if parent_pdg==[]: parent_pdg=[0]
    return parent_pdg

            

##### SAME PDG CONNECTED TRAJECTORIES -----------------------


def same_pdg_connected_trajectories(track_pdg, track_id, vertex_assoc_traj,\
                                    traj, ghdr):
    trackid_set = {track_id}
    
    # walk up the family tree
    this_pdg = track_pdg
    this_track_id = track_id

    while this_pdg==track_pdg:
        particle_mask = vertex_assoc_traj['trackID'] == this_track_id
        this_pdg = find_parent_pdg(vertex_assoc_traj[particle_mask]['parentID'],
                                   vertex_assoc_traj[particle_mask]['vertexID'],
                                   traj, ghdr)
        if this_pdg==track_pdg:
            this_track_id = vertex_assoc_traj[particle_mask]['parentID'].tolist()[0]
            trackid_set.add(this_track_id)

    # walk down the family tree
    this_pdg = track_pdg
    this_track_id = track_id
    while this_pdg==track_pdg:
        particle_mask = vertex_assoc_traj['parentID'] == this_track_id
        child_particle=vertex_assoc_traj[particle_mask]
        if len(child_particle)==1:
            this_pdg = child_particle['pdgId']
            if child_particle['pdgId']==track_pdg:
                this_track_id = child_particle['trackID'].tolist()[0]
                trackid_set.add(this_track_id)
        else: break

    return trackid_set
    

            
##### FIDUCIAL VOLUME/ TOTAL VOLUME ENERGY DEPOSITION -------


def tpc_edep_charged_e(trackID, traj, seg):
    seg_id_mask=seg['trackID']==trackID
    contained_e={}
    for i in range(8): contained_e[i]=0.
    for sg in seg[seg_id_mask]:
        tpc_fv=twoBytwo_defs.tpc_vertex([(sg['x_start']+sg['x_end'])/2.,
                                         (sg['y_start']+sg['y_end'])/2.,
                                         (sg['z_start']+sg['z_end'])/2.])
        for key in tpc_fv.keys():
            if tpc_fv[key]==True:
                contained_e[key]+=sg['dE']
    return contained_e


            
def fv_edep_charged_e(trackID, traj, seg):
    seg_id_mask=seg['trackID']==trackID
    contained_e=0.
    for sg in seg[seg_id_mask]:
        if twoBytwo_defs.fiducialized_vertex([(sg['x_start']+sg['x_end'])/2.,
                                              (sg['y_start']+sg['y_end'])/2.,
                                              (sg['z_start']+sg['z_end'])/2.]):
            contained_e+=sg['dE']
    return contained_e



def total_edep_charged_e(trackID, traj, seg):
    seg_id_mask=seg['trackID']==trackID
    total_e=0.
    for sg in seg[seg_id_mask]: total_e+=sg['dE']
    return total_e



def fv_edep_neutral_e(trackID, traj, seg):
    flag=True; daughter_track_ids=set()
    temp_track_id=[trackID]
    while flag==True:
        second_temp_track_id=[]
        for ttid in temp_track_id:
            parent_mask = traj['parentID']==trackID
            daughter_id=traj[parent_mask]['trackID']
            if len(daughter_id)!=0:
                for did in daughter_id:
                    daughter_track_ids.add(did)
                    second_temp_track_id.append(did)
            else: flag=False
        temp_track_id=second_temp_track_id
    contained_e=0
    print('[CONTAINED EDEP] NEUTRAL daughter track IDs: ', daughter_track_ids)
    for dti in daughter_track_ids:
        contained_e+=fv_edep_charged_e(dti, traj, seg)
    return contained_e



def total_edep_neutral_e(trackID, traj, seg):
    flag=True; daughter_track_ids=set()
    temp_track_id=[trackID]
    while flag==True:
        second_temp_track_id=[]
        for ttid in temp_track_id:
            parent_mask = traj['parentID']==trackID
            daughter_id=traj[parent_mask]['trackID']
            if len(daughter_id)!=0:
                for did in daughter_id:
                    daughter_track_ids.add(did)
                    second_temp_track_id.append(did)
            else: flag=False
        temp_track_id=second_temp_track_id
    total_e=0
    print('[TOTAL EDEP] NEUTRAL daughter track IDs: ', daughter_track_ids)
    for dti in daughter_track_ids:
        total_e+=total_edep_charged_e(dti, traj, seg)
    return total_e


def tpc_contained_energy(pdg, trackID, traj, seg):
    if pdg in neutral_pdg:
        temp={}
        for i in range(8): temp[i]=0.
        return temp
    else: return tpc_edep_charged_e(trackID, traj, seg)

            

def fv_contained_energy(pdg, trackID, traj, seg):
    if pdg in neutral_pdg: return 0.#fv_edep_neutral_e(trackID, traj, seg)
    else: return fv_edep_charged_e(trackID, traj, seg)


    
def total_energy(pdg, trackID, traj, seg):
    if pdg in neutral_pdg: return 0.#total_edep_neutral_e(trackID, traj, seg)
    else: return total_edep_charged_e(trackID, traj, seg)



##### FIDUCIAL VOLUME/ TOTAL VOLUME LENGTH ------------------------
    


def fv_edep_charged_length(trackID, traj, seg):
    seg_id_mask=seg['trackID']==trackID
    contained_length=0.
    for sg in seg[seg_id_mask]:
        if twoBytwo_defs.fiducialized_vertex([(sg['x_start']+sg['x_end'])/2.,
                                              (sg['y_start']+sg['y_end'])/2.,
                                              (sg['z_start']+sg['z_end'])/2.]):
            contained_length+=np.sqrt( (sg['x_start']-sg['x_end'])**2+
                                       (sg['y_start']-sg['y_end'])**2.+
                                       (sg['z_start']-sg['z_end'])**2. )
    return contained_length



def total_edep_charged_length(trackID, traj, seg):
    seg_id_mask=seg['trackID']==trackID
    contained_length=0.
    for sg in seg[seg_id_mask]:
            contained_length+=np.sqrt( (sg['x_start']-sg['x_end'])**2+
                                       (sg['y_start']-sg['y_end'])**2.+
                                       (sg['z_start']-sg['z_end'])**2. )
    return contained_length



def fv_contained_length(pdg, trackID, traj, seg):
    if pdg in neutral_pdg: return 0.
    else: return fv_edep_charged_length(trackID, traj, seg)

    

def total_length(pdg, trackID, traj, seg):
    if pdg in neutral_pdg: return 0.
    else: return total_edep_charged_length(trackID, traj, seg)




##### MISCELLANEOUS ----------------------------------------------------------------------------

def np_array_of_array_to_flat_list(a):
    b = list(a)
    return [list(c)[0] for c in b]
