import matplotlib
import matplotlib.pyplot as plt
import twoBytwo_defs
import numpy as np


'''TO DO:
(2) Add text labelling for pdg ids  (NEW METHOD)
'''
        
def muon_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, muon_dict):

    traj_vert_mask = traj['vertexID']==vert_id
    final_states = traj[traj_vert_mask]

    ghdr_vert_mask = ghdr['vertexID']==vert_id
    truth_level_summ = ghdr[ghdr_vert_mask]

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

        if fs['pdgId'] != -13: continue

        #print("Particle is muon.")

        pdg = fs['pdgId'] # *** pdg ***     
        if fs['parentID']==-1:
            ghdr_vertexID_mask=ghdr['vertexID']==fs['vertexID']
            parent_pdg=ghdr[ghdr_vertexID_mask]['nu_pdg']
        else:
            parent_trackID_mask = traj['trackID']==fs['parentID']
            parent_pdg=traj[parent_trackID_mask]['pdgId']
        parent_pdg = parent_pdg.tolist()[0] # *** parent pdg ***
        
        track_id = fs['trackID']
        seg_id_mask = seg['trackID']==track_id
        total_edep = sum(seg[seg_id_mask]['dE']) # *** total visible energy ***

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
                                                                                                                                          
    muon_dict[(spill_id,vert_id)]=dict(
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
    return


    ''' TO DO: (1) Define visible/contained energy for hadrons from one vertex
               (2) Define visible/contained length for hadrons from one vertex
               (3) Define variables related to hadron ids/multiplicities from interaction'''
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

    ### GET: Total and contained energy deposits of hadrons
    total_edep=0.; contained_edep=0.; total_length=0.; contained_length=0.
    for fs in final_states:

        if abs(fs['pdgId']) in leptons_abs_pdg: continue # No leptons
        if fs['pdgId'] > 1000000000: continue # No nuclei
        if fs['pdgId'] == 22: continue # No photons
        
        track_id = fs['trackID']

        seg_id_mask = seg['trackID']==track_id
        total_edep = sum(seg[seg_id_mask]['dE']) # *** total visible energy ***

        for sg in seg[seg_id_mask]:
            total_length+=np.sqrt((sg['x_start']-sg['x_end'])**2+
                                  (sg['y_start']-sg['y_end'])**2+
                                  (sg['z_start']-sg['z_end'])**2) # *** total length ***                                                                                                              
        
        # Save contained energy and length of hadrons
        for sg in seg[seg_id_mask]:
            if twoBytwo_defs.fiducialized_vertex([(sg['x_start']+sg['x_end'])/2.,
                                                  (sg['y_start']+sg['y_end'])/2.,
                                                  (sg['z_start']+sg['z_end'])/2.] ):
                contained_edep+=sg['dE'] # *** contained visible energy ***                                                                                                                           
                contained_length+=np.sqrt((sg['x_start']-sg['x_end'])**2+
                                              (sg['y_start']-sg['y_end'])**2+
                                              (sg['z_start']-sg['z_end'])**2) # *** contained length ***
            
    hadron_dict[(spill_id,vert_id, track_id)]=dict(
        hadron_mult = hadron_mult,
        hadron_pdg = gstack_vert_fs_hadrons,
        hadron_pdg_set = gstack_vert_fs_pdg_set,
        total_edep=total_edep,
        contained_edep=contained_edep,
        total_length=total_length,
        contained_length=contained_length)
    return
                                      
# PLOT: Muon kinematics
#       sig_bkg is an int such that 0 == signal, 1 == 'dirt' backgrounds, 2 == 'beam' backgrounds
def plot_muons(d, sig_bkg = 0):
    
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
        
                                                  
    # PLOT: total visible energy + contained visible energy
    fig0, ax0 = plt.subplots(figsize=(12,4))
    data0tot = np.array([d[key]['total_edep'] for key in d.keys()])
    data0cont = np.array([d[key]['contained_edep'] for key in d.keys()])
    bins0=np.linspace(0,1000,50)
    ax0.hist(data0tot, bins=bins0, label='total', histtype='step')
    ax0.hist(data0cont, bins=bins0, label='contained',histtype='step')
    ax0.set_xlabel('Visible Energy [MeV]')
    ax0.set_ylabel('Count')
    ax0.set_title(r'Muon Energy')
    ax0.legend()
    ax0.grid(True)
    plt.savefig(sample_type+"_events_muon_visible_energy.png")

    # PLOT: muon energy containment fraction
    fig1, ax1 = plt.subplots(figsize=(6,4))
    data1 = np.array([d[key]['contained_edep']/d[key]['total_edep'] for key in d.keys() if d[key]['total_edep']!=0])
    bins1=np.linspace(0,1,20)
    ax1.hist(data1, bins=bins1, histtype='step')
    ax1.set_xlabel('Visible Muon Energy Containment Fraction')
    ax1.set_ylabel('Count')
    ax1.grid(True)       
    plt.savefig(sample_type+"_events_muon_containment_fraction.png")    
    
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

# PLOT: Hadron kinematics
#       sig_bkg is an int such that 0 == signal, 1 == 'dirt' backgrounds, 2 == 'beam' backgrounds
def plot_hadrons(d, sig_bkg = 0):
    
    # DEFINE: Plotting muon kinematics for signal or background events
    sample_type = ''
    if sig_bkg == 0: 
        sample_type = 'signal'
    elif sig_bkg == 1:
        sample_type = 'dirt_bkg'
    elif sig_bkg == 2:
        sample_type = 'beam_bkg'
    else: 
        return "Error: plot_hadrons function given undefined signal/background definition"
    
    # PLOT: total visible energy + contained visible energy
    fig0, ax0 = plt.subplots(figsize=(12,4))
    data0tot = np.array([d[key]['total_edep'] for key in d.keys()])
    data0cont = np.array([d[key]['contained_edep'] for key in d.keys()])
    bins0=np.linspace(0,1000,50)
    ax0.hist(data0tot, bins=bins0, label='total', histtype='step')
    ax0.hist(data0cont, bins=bins0, label='contained',histtype='step')
    ax0.set_xlabel('Hadron Visible Energy [MeV]')
    ax0.set_ylabel('Count')
    ax0.set_title(r'Muon Energy')
    ax0.legend()
    ax0.grid(True)
    plt.savefig(sample_type+"_events_hadron_visible_energy.png")

    # PLOT: hadron energy containment fraction
    fig1, ax1 = plt.subplots(figsize=(6,4))
    data1 = np.array([d[key]['contained_edep']/d[key]['total_edep'] for key in d.keys() if d[key]['total_edep']!=0])
    bins1=np.linspace(0,1,20)
    ax1.hist(data1, bins=bins1, histtype='step')
    ax1.set_xlabel('Visible Hadron Energy Containment Fraction')
    ax1.set_ylabel('Count')
    ax1.grid(True)       
    plt.savefig(sample_type+"_events_hadron_energy_containment_fraction.png")

    # PLOT: truth-level 4-momentum squared of interaction
    fig2, ax2 = plt.subplots(figsize=(6,4))
    data2 = np.array([d[key]['hadron_mult'] for key in d.keys()])
    bins2 = np.linspace(-0.5,20.5,22)
    ax2.hist(data2, bins=bins2, histtype='step')
    ax2.set_xlabel(r"Hadron Multiplicity")
    ax2.set_ylabel("Count") 
    plt.savefig(sample_type+"_events_hadron_multiplicity_truth.png")    

    # PLOT: truth-level 4-momentum squared of interaction
    fig3, ax3 = plt.subplots(figsize=(6,4))
    hadron_fs_pdg_list=[sorted(d[key]['hadron_pdg_set']) for key in d.keys()]
    hadron_fs_pdg_set=set(tuple(pdg) for pdg in hadron_fs_pdg_list)
    #print("Hadron PDG List:", hadron_fs_pdg_list)
    #print("Hadron PDG Set:", hadron_fs_pdg_set)
    hadron_fs_pdg_count=[(pdg_set, hadron_fs_pdg_list.count(list(pdg_set))) for pdg_set in hadron_fs_pdg_set]
    hadron_fs_pdg_fraction=[100*(i[1]/len(data2)) for i in hadron_fs_pdg_count]
    hadron_fs_pdg_labels=[''.join(str(j)+"," for j in i[0]) for i in hadron_fs_pdg_count]
    #print("Number of Events:", len(hadron_fs_pdg_list))
    #print("Hadron FS PDG Count:", hadron_fs_pdg_count)
    #print("Hadron FS PDG Fractions:", hadron_fs_pdg_fraction)
    #print("Hadron FS PDG Labels:", hadron_fs_pdg_labels)
    ax3.pie(hadron_fs_pdg_fraction, labels=hadron_fs_pdg_labels, autopct='%1.1f%%')
    ax3.set_title(r"Final State Hadron PDG IDs")
    plt.savefig(sample_type+"_events_hadron_pdg_ids_truth.png")    
    