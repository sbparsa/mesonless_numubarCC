import matplotlib
import matplotlib.pyplot as plt
import twoBytwo_defs
import numpy as np

def dirt_muon_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, dirt_dict):

    traj_vert_mask = traj['vertexID']==vert_id
    final_states = traj[traj_vert_mask]

    ghdr_vert_mask = ghdr['vertexID']==vert_id
    truth_level_summ = ghdr[ghdr_vert_mask]

    ## thg=is works if event is cc
    # Save truth-level outgoing muon momentum   
    true_mom = truth_level_summ['lep_mom']
    
    # Save truth-level muon angle with beam                                                                                 
    true_angle = truth_level_summ['lep_ang']

    # Save truth-level outgoing muon energy                                                                                
    true_energy = truth_level_summ['Elep']
    
    # Save truth-level neutrino energy                                                                                   
    nu_energy = truth_level_summ['Enu']

    # Save truth-level interaction 4-momentum squared                                                                     
    q_sq = truth_level_summ['Q2']

    total_edep=0.; contained_edep=0.; total_length=0.; contained_length=0.
    

    #print("PDG IDs of F.S. Particles:", final_states['pdgId'])
    gstack_vert_mask = gstack['vertexID']==vert_id
    gstack_pdg_set = set(gstack[gstack_vert_mask]['part_pdg'])
    #print("Event PDG Stack:", gstack_pdg_set)
    
    for fs in final_states:
        if fs['pdgId'] not in [13, -13]: continue # [111,211,-211]: continue

        ##place a cut for backgrounds
        #print('muon start point', fs['xyz_start'])
        #print('muon end point', fs['xyz_end'])
        if( twoBytwo_defs.fiducialized_vertex(fs['xyz_start'])):
            print('This is a dirt background')
            print('Neutrino vertex: ', truth_level_summ['vertex'])
            print('muon start point: ', fs['xyz_start'])
            print('muon end point: ', fs['xyz_end'])
        else:
            #print('started out of FV')
            continue

        pdg = fs['pdgId'] # *** pdg ***
        
        parent_id = fs['parentID']
        parent = final_states['trackID']==parent_id

        parent_pdg = final_states[parent]['pdgId'] # *** parent pdg ***

        if parent_id ==-1:
            ghdr_mask = ghdr['vertexID']==fs['vertexID']
            parent_pdg = ghdr[ghdr_mask]['nu_pdg']

        parent_pdg = parent_pdg.tolist()[0]
        #print('check parent pdg: ', parent_pdg)
        
        track_id = fs['trackID']
        total_edep=0.; contained_edep=0.; total_length=0.; contained_length=0.

        #### checking dirt muons
        if abs(pdg)==13:
            seg_id_mask = seg['trackID']==track_id
            total_edep = sum(seg[seg_id_mask]['dE']) # *** total visible energy ***
            contained_edep = 0
            contained_length = 0
            
            for sg in seg[seg_id_mask]:
                total_length+=np.sqrt((sg['x_start']-sg['x_end'])**2+
                                      (sg['y_start']-sg['y_end'])**2+
                                      (sg['z_start']-sg['z_end'])**2) # *** total length *** This is not correct

            for sg in seg[seg_id_mask]:
                if twoBytwo_defs.fiducialized_vertex( [(sg['x_start']+sg['x_end'])/2.,
                                         (sg['y_start']+sg['y_end'])/2.,
                                         (sg['z_start']+sg['z_end'])/2.] ):
                    contained_edep+=sg['dE'] # *** contained visible energy ***

                    contained_length+=np.sqrt((sg['x_start']-sg['x_end'])**2+
                                              (sg['y_start']-sg['y_end'])**2+
                                              (sg['z_start']-sg['z_end'])**2) # *** contained length ***

        if contained_edep>5:
            print(pdg,'\t',parent_pdg,'\t',total_edep,' MeV\t',contained_edep,' MeV\t', total_length,' cm\t',contained_length,' cm')


            dirt_dict[(spill_id,vert_id, track_id)]=dict(
                pdg=int(pdg),
                parent_pdg=parent_pdg,
                total_edep=total_edep,
                contained_edep=contained_edep,
                total_length=total_length,
                contained_length=contained_length,
                true_mom=true_mom,
                true_angle=true_angle,
                true_energy=true_energy,
                nu_energy=nu_energy,
                q_sq = q_sq)
            
            ## end if
    return


def plot_dirt_backgrounds(d , sig_bkg =0):

    # DEFINE: Plotting muon kinematics for signal or background events
    
    sample_type = ''
    if sig_bkg == 0:
        sample_type = 'signal'
    elif sig_bkg == 1:
        sample_type = 'dirt_bkg'
    elif sig_bkg == 2:
        sample_type = 'beam_bkg'
    else:
        return "Error: plot_muons function given undefined signal/background definition"
    
    #PLOT: total visible energy + contained visible energy
    fig0, ax0 = plt.subplots(1,2,figsize=(8,4))
    bins=np.linspace(0,1000,50)
    ax0[0].hist([d[key]['total_edep'] for key in d.keys() if d[key]['pdg']==-13],
                bins=bins, label='total', histtype='step')
    ax0[1].hist([d[key]['total_edep'] for key in d.keys() if d[key]['pdg']==13],
                bins=bins, label='total', histtype='step')

    ax0[0].hist([d[key]['contained_edep'] for key in d.keys() if d[key]['pdg']==-13],
                bins=bins, label='contained',histtype='step')
    ax0[1].hist([d[key]['contained_edep'] for key in d.keys() if d[key]['pdg']==13],
                bins=bins, label='contained', histtype='step')

    for i in range(2):
        ax0[i].set_xlabel('Visible Energy [MeV]')
        ax0[i].set_ylabel('Count / 20 MeV')
        ax0[i].grid(True)
        if i==0: ax0[i].set_title(r'$\mu^+$'); ax0[i].legend()
        if i==1: ax0[i].set_title(r'$\mu^-$')    
    plt.savefig(sample_type + '_muon_vis_edep.png')

    #PLOT:c containment fraction
    fig1, ax1 = plt.subplots(figsize=(6,4))
    bins=np.linspace(0,1,20)
    ax1.hist([d[key]['contained_edep']/d[key]['total_edep'] for key in d.keys() if d[key]['pdg']==13 and d[key]['total_edep']!=0],
             bins=bins, label=r'$\mu^+$', histtype='step')
    ax1.hist([d[key]['contained_edep']/d[key]['total_edep'] for key in d.keys() if d[key]['pdg']==-13 and d[key]['total_edep']!=0],
             bins=bins, label=r'$\mu^-$', histtype='step')

    ax1.set_xlabel('Visible Energy Containment Fraction')
    ax1.set_ylabel('Count')
    ax1.legend()
    ax1.grid(True)
    #plt.show()
    plt.savefig(sample_type+'_muon_edep_containment_fraction,png')

    #PLOT: parent pdg pi chart
    fig2, ax2 = plt.subplots(1,2,figsize=(8,4))
    muplus_parent_pdg_list=[d[key]['parent_pdg'] for key in d.keys() if d[key]['pdg']==-13]
    print('check parent pdg lisr: ', muplus_parent_pdg_list)

    muplus_parent_pdg_set=set(muplus_parent_pdg_list)
    muplus_particle_count=[(pdg, muplus_parent_pdg_list.count(pdg)) for pdg in muplus_parent_pdg_set]
    muplus_fraction=[100*(i[1]/len(muplus_parent_pdg_list)) for i in muplus_particle_count]
    muplus_pdg_label=[str(i[0]) for i in muplus_particle_count]
    print('$\mu^+$: ',muplus_particle_count)
    ax2[0].pie(muplus_fraction, labels=muplus_pdg_label, autopct='%1.1f%%')
    ax2[0].set_title(r'$\mu^+$')
    
    muminus_parent_pdg_list=[d[key]['parent_pdg'] for key in d.keys() if d[key]['pdg']==13]
    muminus_parent_pdg_set=set(muminus_parent_pdg_list)
    muminus_particle_count=[(pdg, muminus_parent_pdg_list.count(pdg)) for pdg in muminus_parent_pdg_set]
    muminus_fraction=[100*(i[1]/len(muminus_parent_pdg_list)) for i in muminus_particle_count]
    muminus_pdg_label=[str(i[0]) for i in muminus_particle_count]
    print('$\mu^-$: ',muminus_particle_count)
    ax2[1].pie(muminus_fraction, labels=muminus_pdg_label, autopct='%1.1f%%')
    ax2[1].set_title(r'$\mu^-$')
    #plt.show()
    plt.savefig(sample_type+'_muon_parent_pdg.png')

    #PLOT: Visible length
    fig3, ax3 = plt.subplots(1,2,figsize=(8,4))
    bins=np.linspace(0,400,40)
    ax3[0].hist([d[key]['total_length'] for key in d.keys() if d[key]['pdg']==-13],
                bins=bins, label='total', histtype='step')
    ax3[1].hist([d[key]['total_length'] for key in d.keys() if d[key]['pdg']==13],
                bins=bins, label='total', histtype='step')
    ax3[0].hist([d[key]['contained_length'] for key in d.keys() if d[key]['pdg']==-13],
                bins=bins, label='contained',histtype='step')
    ax3[1].hist([d[key]['contained_length'] for key in d.keys() if d[key]['pdg']==13],
                bins=bins, label='contained', histtype='step')
    for i in range(2):
        ax3[i].set_xlabel('Track Length [cm]')
        ax3[i].set_ylabel('Count / 10 cm')
        ax3[i].grid(True)
        if i==0: ax3[i].set_title(r'$\mu^+$'); ax3[i].legend()
        if i==1: ax3[i].set_title(r'$\mu^-$')

    #plt.show()
    plt.savefig(sample_type+'_muon_track_length')
    
    # PLOT: truth-level outgoing muon (lepton) angle                                                                
    fig2, ax2 = plt.subplots(figsize=(6,4))
    data2 = np.cos((np.pi / 180.)*np.array([d[key]['true_angle'] for key in d.keys()]))
    bins2 = np.linspace(0.75,1,51)
    ax2.hist(data2, bins=bins2, histtype='step')
    ax2.set_xlabel(r"Cosine of Outgoing Muon Angle with Beam Direction")
    ax2.set_ylabel("Count")
    plt.savefig(sample_type+"_events_outgoing_muon_angle_truth.png")

    # PLOT: truth-level outgoing muon (lepton) momentum                                                                
    fig3, ax3 = plt.subplots(figsize=(6,4))
    data3 = np.array([d[key]['true_mom'] for key in d.keys()])
    bins3 = np.linspace(0,30000,21)
    ax3.hist(data3, bins=bins3, histtype='step')
    ax3.set_xlabel(r"Outgoing Muon Momentum [MeV]")
    ax3.set_ylabel("Count")
    plt.savefig(sample_type+"_events_outgoing_muon_momentum_truth.png")

    # PLOT: truth-level 4-momentum squared of interaction                                                                
    fig4, ax4 = plt.subplots(figsize=(6,4))
    data4 = np.array([d[key]['q_sq'] for key in d.keys()])
    bins4 = np.linspace(0,3000000,21)
    ax4.hist(data4, bins=bins4, histtype='step')
    ax4.set_xlabel(r"4-Momentum Transfer Squared [MeV$^2$]")
    ax4.set_ylabel("Count")
    plt.savefig(sample_type+"_events_qsq_truth.png")

    # PLOT: truth-level neutrino energy of interaction                                                                  
    fig5, ax5 = plt.subplots(figsize=(6,4))
    data5 = np.array([d[key]['nu_energy'] for key in d.keys()])
    bins5 = np.linspace(0,20000,21)
    ax5.hist(data5, bins=bins5, histtype='step')
    ax5.set_xlabel(r"Incident Neutrino Energy [MeV]")
    ax5.set_ylabel("Count")
    plt.savefig(sample_type+"_events_nu_energy_truth.png")
