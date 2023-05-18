import matplotlib
import matplotlib.pyplot as plt
import twoBytwo_defs
import numpy as np
import auxiliary


'''TO DO:
(2) Add text labelling for pdg ids  (NEW METHOD)
'''

def get_truth_dict(spill_id, vert_id, ghdr, gstack, traj, vert, seg, signal_dict):

    traj_vert_mask = traj['vertexID']==vert_id
    final_states = traj[traj_vert_mask]

    ghdr_vert_mask = ghdr['vertexID']==vert_id
    truth_level_summ = ghdr[ghdr_vert_mask]

    mom = truth_level_summ['lep_mom'] # Save truth-level outgoing muon momentum
    ang = truth_level_summ['lep_ang'] *np.pi / 180. # Save truth-level muon angle with beam
    nu_energy = truth_level_summ['Enu'] # Save truth-level neutrino energy
    q2 = truth_level_summ['Q2'] # Save truth-level interaction 4-momentum squared
    vtx = truth_level_summ['vertex'] # Save truth-level vertex information
    vtx_x = vtx[0][0]
    vtx_y = vtx[0][1]
    vtx_z = vtx[0][2]
                                                                                                                                          
    signal_dict[(spill_id,vert_id)]=dict(
        nu_energy=float(nu_energy),
        q2 = float(q2),
        mom=float(mom), 
        ang=float(ang),
        vtx_x = float(vtx_x), 
        vtx_y = float(vtx_y), 
        vtx_z = float(vtx_z))
    return
        
def muon_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, muon_dict):

    traj_vert_mask = traj['vertexID']==vert_id
    final_states = traj[traj_vert_mask]

    ghdr_vert_mask = ghdr['vertexID']==vert_id
    truth_level_summ = ghdr[ghdr_vert_mask]

    mom = truth_level_summ['lep_mom'] # Save truth-level outgoing muon momentum
    ang = truth_level_summ['lep_ang'] *np.pi / 180. # Save truth-level muon angle with beam
    nu_energy = truth_level_summ['Enu'] # Save truth-level neutrino energy
    q2 = truth_level_summ['Q2'] # Save truth-level interaction 4-momentum squared

    total_edep=0.; contained_edep=0.; total_length=0.; contained_length=0.
    #print("PDG IDs of F.S. Particles:", final_states['pdgId'])
    gstack_vert_mask = gstack['vertexID']==vert_id
    gstack_pdg_set = set(gstack[gstack_vert_mask]['part_pdg'])
    #print("Event PDG Stack:", gstack_pdg_set)
    for fs in final_states:

        if fs['pdgId'] != -13: continue

        #print("Particle is muon.")

        pdg = fs['pdgId'] # *** pdg ***     
        parent_pdg = auxiliary.find_parent_pdg(fs['parentID'],vert_id, traj, ghdr)# *** parent pdg ***
        
        track_id = fs['trackID']
        seg_id_mask = seg['trackID']==track_id
        total_edep += sum(np.array(seg[seg_id_mask]['dE'])) # *** total visible energy ***

        for sg in seg[seg_id_mask]:
            total_length+=np.sqrt((sg['x_start']-sg['x_end'])**2+
                                  (sg['y_start']-sg['y_end'])**2+
                                  (sg['z_start']-sg['z_end'])**2) # *** total length ***                                                                                                              
        
        # Save contained energy and length of muon track
        for sg in seg[seg_id_mask]:
            if twoBytwo_defs.fiducialized_vertex([(sg['x_start']+sg['x_end'])/2.,
                                                  (sg['y_start']+sg['y_end'])/2.,
                                                  (sg['z_start']+sg['z_end'])/2.] ):
                contained_edep+=sg['dE'] # *** contained visible energy ***                                                                                                                           
                contained_length+=np.sqrt((sg['x_start']-sg['x_end'])**2+
                                              (sg['y_start']-sg['y_end'])**2+
                                              (sg['z_start']-sg['z_end'])**2) # *** contained length ***


     #   # *** Characterize Muon Endpoint/Containment ***
        end_pt = fs['xyz_end']
        start_pt = fs['xyz_start']
        end_pt_loc = twoBytwo_defs.particle_end_loc(start_pt, end_pt)
                                                                                                                                  
    muon_dict[(spill_id,vert_id)]=dict(
        pdg=int(pdg),
        parent_pdg=int(parent_pdg),
        total_edep=float(total_edep),
        contained_edep=float(contained_edep),
        total_length=float(total_length),
        contained_length=float(contained_length),#end_pt_loc = end_pt_loc,
        mom=float(mom), 
        ang=float(ang),
        nu_energy=float(nu_energy),
        q2 = float(q2),
        end_pt_loc = str(end_pt_loc))
    return


def hadron_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, hadron_dict):
        
    traj_vert_mask = traj['vertexID']==vert_id
    final_states = traj[traj_vert_mask]

    leptons_abs_pdg = [11, 12, 13, 14, 15, 16]
    
    ### GET: Hadron multiplicity for this vertex
    #print("PDG IDs of F.S. Particles:", final_states['pdgId'])
    gstack_vert_mask = gstack['vertexID']==vert_id
    gstack_vert = gstack[gstack_vert_mask] # Get particles from interaction
    gstack_pdg_set = set(gstack_vert['part_pdg'])

    gstack_vert_fs_mask = gstack_vert['part_status']==1 # Exclude initial state particles

    gstack_vert_fs = gstack_vert[gstack_vert_fs_mask]['part_pdg'] # Get final state particle PDG IDs

    gstack_vert_fs_hadrons = [fsp for fsp in gstack_vert_fs if abs(fsp) not in leptons_abs_pdg and fsp != 22] # Exclude f.s. leptons and photons 
    gstack_vert_fs_pdg_set = set(gstack_vert_fs_hadrons) # Get f.s. particle PDG IDs in set
    #gstack_vert_fs_hadrons = gstack_vert_fs[gstack_lep_phot_mask]

    #print("Event PDG Stack:", gstack_vert['part_pdg'])
    #print("Final State Hadrons:", gstack_vert_fs_hadrons)
    #print("Final State PDG Stack:", gstack_vert_fs)
    #print("Event PDG Stack Set:", gstack_pdg_set)

    hadron_mult = len(gstack_vert_fs_hadrons)
    n_mult = 0
    p_mult = 0
    other_had_mult = 0

    for had in range(hadron_mult):

        if gstack_vert_fs_hadrons[had] == 2112:
            n_mult+=1
        elif gstack_vert_fs_hadrons[had] == 2212:
            p_mult+=1
        else:
            other_had_mult+=1

    ### GET: Total and contained energy deposits of hadrons
    total_edep=0.; contained_edep=0.; total_length=0.; contained_length=0.
    for fs in final_states:

        if abs(fs['pdgId']) in leptons_abs_pdg: continue # No leptons
        if fs['pdgId'] > 1000000000: continue # No nuclei
        if fs['pdgId'] == 22: continue # No photons
        
        track_id = fs['trackID']
        seg_id_mask = seg['trackID']==track_id
        total_edep += sum(np.array(seg[seg_id_mask]['dE'])) # *** total visible energy ***
        #print('List of visible hadron energy components:', np.array(seg[seg_id_mask]['dE']))

        for sg in seg[seg_id_mask]:
            total_length+=np.sqrt((sg['x_start']-sg['x_end'])**2+
                                  (sg['y_start']-sg['y_end'])**2+
                                  (sg['z_start']-sg['z_end'])**2) # *** total length ***                                                                                                              
        
        # Save contained energy and length of hadrons
        #print('Initial Contained Edep:', contained_edep)
        for sg in seg[seg_id_mask]:
            if twoBytwo_defs.fiducialized_vertex([(sg['x_start']+sg['x_end'])/2.,
                                                  (sg['y_start']+sg['y_end'])/2.,
                                                  (sg['z_start']+sg['z_end'])/2.] ):
                contained_edep+=sg['dE'] # *** contained visible energy ***  
                #print("Contained Hadron Energy component:", sg['dE'])                                                                                                                         
                contained_length+=np.sqrt((sg['x_start']-sg['x_end'])**2+
                                              (sg['y_start']-sg['y_end'])**2+
                                              (sg['z_start']-sg['z_end'])**2) # *** contained length ***
    
    #print('Total Hadron Edep calculated:', total_edep)          
    #print('Contained Hadron Edep calculated:', contained_edep)
            
    hadron_dict[(spill_id,vert_id, track_id)]=dict(
        hadron_mult = hadron_mult,
        neutron_mult = n_mult,
        proton_mult = p_mult,
        other_had_mult = other_had_mult,
        hadron_pdg = gstack_vert_fs_hadrons,
        hadron_pdg_set = gstack_vert_fs_pdg_set,
        total_edep=total_edep,
        contained_edep=contained_edep,
        total_length=total_length,
        contained_length=contained_length)
    return
                                                 


    