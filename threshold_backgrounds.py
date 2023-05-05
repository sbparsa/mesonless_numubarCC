import twoBytwo_defs


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
                if twoBytwo_defs.fiducialized_vertex( [(sg['x_start']+sg['x_end'])/2.,
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
                    if twoBytwo_defs.fiducialized_vertex( [(sg['x_start']+sg['x_end'])/2.,
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
