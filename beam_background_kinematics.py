import matplotlib
import matplotlib.pyplot as plt
import h5py
import argparse
import numpy as np
import twoBytwo_defs

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



def fiducialized_particle_origin(traj, vert_id):
    traj_vert_mask = traj['vertexID']==vert_id
    final_states = traj[traj_vert_mask]
    for fs in final_states:
        if fiducialized_vertex(fs['xyz_start'])==True:
            return True
    return False


    
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


    
def pion_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, pion_dict):
    traj_vert_mask = traj['vertexID']==vert_id
    final_states = traj[traj_vert_mask]

    for fs in final_states:
        if fs['pdgId'] not in [111,211,-211]: continue

        pdg = fs['pdgId'] # *** pdg ***

        parent_id = fs['parentID']
        parent = final_states['trackID']==parent_id
        parent_pdg = final_states['pdgId'][0] # *** parent pdg ***
        
        track_id = fs['trackID']
        total_edep=0.; contained_edep=0.; total_length=0.; contained_length=0.
        
        if abs(pdg)==211:
            seg_id_mask = seg['trackID']==track_id
            total_edep = sum(seg[seg_id_mask]['dE']) # *** total visible energy ***            
            
            for sg in seg[seg_id_mask]:
                total_length+=np.sqrt((sg['x_start']-sg['x_end'])**2+
                                      (sg['y_start']-sg['y_end'])**2+
                                      (sg['z_start']-sg['z_end'])**2) # *** total length ***

            for sg in seg[seg_id_mask]:
                if fiducialized_vertex( [(sg['x_start']+sg['x_end'])/2.,
                                         (sg['y_start']+sg['y_end'])/2.,
                                         (sg['z_start']+sg['z_end'])/2.] ):
                    contained_edep+=sg['dE'] # *** contained visible energy ***
                    contained_length+=np.sqrt((sg['x_start']-sg['x_end'])**2+
                                              (sg['y_start']-sg['y_end'])**2+
                                              (sg['z_start']-sg['z_end'])**2) # *** contained length ***

        ### save all energy depositions hailing from pi0 track ID
        if pdg==111:
            flag=True; daughter_track_ids=set()
            temp_track_id=[track_id]
            while flag==True:
                #print('DAUGHTER TRACK ID SET: ',daughter_track_ids)
                second_temp_track_id=[]
                for ttid in temp_track_id:
                    parent_mask = final_states['parentID']==ttid
                    daughter_id=final_states[parent_mask]['trackID']
                    if len(daughter_id)!=0:
                        for did in daughter_id:
                            daughter_track_ids.add(did)
                            second_temp_track_id.append(did)
                    else:
                        flag=False
                temp_track_id=second_temp_track_id

            for dti in daughter_track_ids:
                seg_id_mask = seg['trackID']==dti
                total_edep+=sum(seg[seg_id_mask]['dE']) # *** total visible energy ***
        
                for sg in seg[seg_id_mask]:
                    if fiducialized_vertex( [(sg['x_start']+sg['x_end'])/2.,
                                             (sg['y_start']+sg['y_end'])/2.,
                                             (sg['z_start']+sg['z_end'])/2.] ):
                        contained_edep+=sg['dE'] # *** contained visible energy ***

#        print(pdg,'\t',parent_pdg,'\t',total_edep,' MeV\t',contained_edep,' MeV\t',
#              total_length,' cm\t',contained_length,' cm')
        pion_dict[(spill_id,vert_id, track_id)]=dict(
            pdg=pdg,
            parent_pdg=parent_pdg,
            total_edep=total_edep,
            contained_edep=contained_edep,
            total_length=total_length,
            contained_length=contained_length)
    return



def plot_threshold_backgrounds(d):
    fig0, ax0 = plt.subplots(1,3,figsize=(12,4))
    bins=np.linspace(0,1000,50)
    ax0[0].hist([d[key]['total_edep'] for key in d.keys() if d[key]['pdg']==211],
                bins=bins, label='total', histtype='step')
    ax0[1].hist([d[key]['total_edep'] for key in d.keys() if d[key]['pdg']==-211],
                bins=bins, label='total', histtype='step')
    ax0[2].hist([d[key]['total_edep'] for key in d.keys() if d[key]['pdg']==111],
                bins=bins, label='total', histtype='step')
    ax0[0].hist([d[key]['contained_edep'] for key in d.keys() if d[key]['pdg']==211],
                bins=bins, label='contained',histtype='step')
    ax0[1].hist([d[key]['contained_edep'] for key in d.keys() if d[key]['pdg']==-211],
                bins=bins, label='contained', histtype='step')
    ax0[2].hist([d[key]['contained_edep'] for key in d.keys() if d[key]['pdg']==111],
                bins=bins, label='contained', histtype='step')
    for i in range(3):
        ax0[i].set_xlabel('Visible Energy [MeV]')
        ax0[i].set_ylabel('Count / 20 MeV')
        ax0[i].grid(True)
        if i==0: ax0[i].set_title(r'$\pi^+$'); ax0[i].legend()
        if i==1: ax0[i].set_title(r'$\pi^-$')
        if i==2: ax0[i].set_title(r'$\pi^0$')
    plt.show()

    fig1, ax1 = plt.subplots(figsize=(6,4))
    bins=np.linspace(0,1,20)
    ax1.hist([d[key]['contained_edep']/d[key]['total_edep'] for key in d.keys() if d[key]['pdg']==211 and d[key]['total_edep']!=0],
             bins=bins, label=r'$\pi^+$', histtype='step')
    ax1.hist([d[key]['contained_edep']/d[key]['total_edep'] for key in d.keys() if d[key]['pdg']==-211 and d[key]['total_edep']!=0],
             bins=bins, label=r'$\pi^-$', histtype='step')
    ax1.hist([d[key]['contained_edep']/d[key]['total_edep'] for key in d.keys() if d[key]['pdg']==111 and d[key]['total_edep']!=0],
             bins=bins, label=r'$\pi^0$', histtype='step')
    ax1.set_xlabel('Visible Energy Containment Fraction')
    ax1.set_ylabel('Count')
    ax1.legend()
    ax1.grid(True)
    plt.show()

    fig2, ax2 = plt.subplots(1,3,figsize=(12,4))
    piplus_parent_pdg_list=[d[key]['parent_pdg'] for key in d.keys() if d[key]['pdg']==211]
    piplus_parent_pdg_set=set(piplus_parent_pdg_list)
    piplus_particle_count=[(pdg, piplus_parent_pdg_list.count(pdg)) for pdg in piplus_parent_pdg_set]
    piplus_fraction=[100*(i[1]/len(piplus_parent_pdg_list)) for i in piplus_particle_count]
    print('pi+: ',piplus_particle_count)
    ax2[0].pie(piplus_fraction)

    piminus_parent_pdg_list=[d[key]['parent_pdg'] for key in d.keys() if d[key]['pdg']==-211]
    piminus_parent_pdg_set=set(piminus_parent_pdg_list)
    piminus_particle_count=[(pdg, piminus_parent_pdg_list.count(pdg)) for pdg in piminus_parent_pdg_set]
    piminus_fraction=[100*(i[1]/len(piminus_parent_pdg_list)) for i in piminus_particle_count]
    print('pi-: ',piminus_particle_count)
    ax2[1].pie(piminus_fraction)

    pi0_parent_pdg_list=[d[key]['parent_pdg'] for key in d.keys() if d[key]['pdg']==111]
    pi0_parent_pdg_set=set(pi0_parent_pdg_list)
    pi0_particle_count=[(pdg, pi0_parent_pdg_list.count(pdg)) for pdg in pi0_parent_pdg_set]
    pi0_fraction=[100*(i[1]/len(pi0_parent_pdg_list)) for i in pi0_particle_count]
    print('pi0: ',pi0_particle_count)
    ax2[2].pie(pi0_fraction)
    ax2[0].set_title(r'$\pi^+$')
    ax2[1].set_title(r'$\pi^-$')
    ax2[2].set_title(r'$\pi^0$')
    plt.show()

    fig3, ax3 = plt.subplots(1,2,figsize=(8,4))
    bins=np.linspace(0,400,40)
    ax3[0].hist([d[key]['total_length'] for key in d.keys() if d[key]['pdg']==211],
                bins=bins, label='total', histtype='step')
    ax3[1].hist([d[key]['total_length'] for key in d.keys() if d[key]['pdg']==-211],
                bins=bins, label='total', histtype='step')
    ax3[0].hist([d[key]['contained_length'] for key in d.keys() if d[key]['pdg']==211],
                bins=bins, label='contained',histtype='step')
    ax3[1].hist([d[key]['contained_length'] for key in d.keys() if d[key]['pdg']==-211],
                bins=bins, label='contained', histtype='step')
    for i in range(2):
        ax3[i].set_xlabel('Track Length [cm]')
        ax3[i].set_ylabel('Count / 10 cm')
        ax3[i].grid(True)
        if i==0: ax3[i].set_title(r'$\pi^+$'); ax3[i].legend()
        if i==1: ax3[i].set_title(r'$\pi^-$')
    plt.show()
    

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
            vert_in_active_LAr = fiducialized_vertex( vert_pos )

            ##### REQUIRE neutrino vertex in LAr active volume #####
            if vert_in_active_LAr==False: continue

            vert_id = vert['vertexID'][v_i]
            
            nu_mu_bar = signal_nu_pdg(ghdr, vert_id)
            is_cc = signal_cc(ghdr, vert_id)
            pionless = signal_pion_status(gstack, vert_id)
            fv_particle_origin=fiducialized_particle_origin(traj, vert_id)
                        
            ##### THRESHOLD BACKGROUNDS #####
            ##### REQUIRE: (A) nu_mu_bar, (B) CC, (C) pions present, (D) final state particle start point in FV
            if nu_mu_bar==True and is_cc==True and pionless==False and fv_particle_origin==True:
                pion_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, pion_dict)
                
    plot_threshold_backgrounds(pion_dict)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--sim_file', default=None, required=True, type=str, help='''string corresponding to the path of the edep-sim ouput simulation file to be considered''')
    parser.add_argument('-t', '--input_type', default='edep', choices=['edep', 'larnd'], type=str, help='''string corresponding to the output file type: edep or larnd''')
    args = parser.parse_args()
    main(**vars(args))
