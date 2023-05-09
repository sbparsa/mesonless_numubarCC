import matplotlib
import matplotlib.pyplot as plt
import twoBytwo_defs
import numpy as np


'''TO DO:
(1) Add method to save dict with hadron info (per interaction) (NEW METHOD)
(2) Add plotting for hadrons (NEW METHOD)
'''
        
def muon_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, muon_dict):

    traj_vert_mask = traj['vertexID']==vert_id
    final_states = traj[traj_vert_mask]

    ghdr_vert_mask = ghdr['vertexID']==vert_id
    truth_level_summ = ghdr[ghdr_vert_mask]

    for fs in final_states:
        if fs['pdgId'] != 13: continue

        pdg = fs['pdgId'] # *** pdg ***                                                                                                                                                                   
        parent_id = fs['parentID']
        parent = final_states['trackID']==parent_id
        parent_pdg = final_states[parent]['pdgId'] # *** parent pdg ***
        
        track_id = fs['trackID']
        total_edep=0.; contained_edep=0.; total_length=0.; contained_length=0.

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
             
        # Save truth-level outgoing muon momentum
        true_mom = truth_level_summ['lep_mom']

        # Save truth-level muon angle with beam
        true_angle = truth_level_summ['lep_ang']

        # Save truth-level outgoing muon energy
        true_energy = truth_level_summ['Elep']
                                                                                                                                          
        muon_dict[(spill_id,vert_id, track_id)]=dict(
            pdg=pdg,
            parent_pdg=parent_pdg,
            total_edep=total_edep,
            contained_edep=contained_edep,
            total_length=total_length,
            contained_length=contained_length,
            true_mom=true_mom, 
            true_angle=true_angle,
            true_energy=true_energy)
    return


    ''' TO DO: (1) Define visible/contained energy for hadrons from one vertex
               (2) Define visible/contained length for hadrons from one vertex
               (3) Define variables related to hadron ids/multiplicities from interaction'''
#def signal_hadron_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, hadron_dict):
#        
#    traj_vert_mask = traj['vertexID']==vert_id
#
#    final_states = traj[traj_vert_mask]
#
#    leptons_and_pions_abs_pdg = [11, 12, 13, 14, 15, 16, 111, 211]
#            
#    for fs in final_states:
#        if abs(fs['pdgId']) in leptons_and_pions_abs_pdg: continue
#        hadron_dict[(spill_id,vert_id, track_id)]=dict(
#            pdg=pdg,
#            parent_pdg=parent_pdg,
#            total_edep=total_edep,
#            contained_edep=contained_edep,
#            total_length=total_length,
#            contained_length=contained_length)
#    return
                                      
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
    bins0=np.linspace(0,1000,50)
    ax0.hist([d[key]['total_edep'] for key in d.keys()],
                bins=bins0, label='total', histtype='step')
    ax0.hist([d[key]['contained_edep'] for key in d.keys()],
                bins=bins0, label='contained',histtype='step')
    ax0.set_xlabel('Visible Energy [MeV]')
    ax0.set_ylabel('Count / 20 MeV')
    ax0.set_title(r'Muon Energy')
    ax0.legend()
    ax0.grid(True)
    plt.show()
    plt.savefig(sample_type+"_events_muon_visible_energy.png")

    # PLOT: muon energy containment fraction
    fig1, ax1 = plt.subplots(figsize=(6,4))
    bins1=np.linspace(0,1,20)
    ax1.hist([d[key]['contained_edep']/d[key]['total_edep'] for key in d.keys() if d[key]['pdg']==13 and d[key]['total_edep']!=0],
             bins=bins1, histtype='step')
    ax1.set_xlabel('Visible Muon Energy Containment Fraction')
    ax1.set_ylabel('Count')
    ax1.legend()
    ax1.grid(True)
    plt.show()      
    plt.savefig(sample_type+"_events_muon_containment_fraction.png")    
    
    # PLOT: truth-level outgoing muon (lepton) angle 
    fig2, ax2 = plt.subplots(figsize=(6,4))
    bins2 = np.linspace(0,66,34)
    ax2.hist([d[key]['true_angle'] for key in d.keys()], bins=bins2)
    ax2.set_xlabel(r"Outgoing Muon Angle with Beam Direction [$^\circ$]")
    ax2.set_ylabel("Count")
    plt.show()
    plt.savefig(sample_type+"_events_outgoing_muon_angle_truth.png")   

    # PLOT: truth-level outgoing muon (lepton) momentum 
    fig3, ax3 = plt.subplots(figsize=(6,4))
    bins3 = np.linspace(0,50,51)
    ax3.hist([d[key]['true_mom']/1000 for key in d.keys()], bins=bins3)
    ax3.set_xlabel(r"Outgoing Muon Momentum [GeV]")
    ax3.set_ylabel("Count")
    plt.show()  
    plt.savefig(sample_type+"_events_outgoing_muon_momentum_truth.png")  

    # PLOT: truth-level outgoing muon (lepton) energy 
    # I think this is redundant with muon momentum plot, but adding just to see
    fig4, ax4 = plt.subplots(figsize=(6,4))
    bins4 = np.linspace(0,50,51)
    ax4.hist([d[key]['true_energy']/1000 for key in d.keys()], bins=bins3)
    ax4.set_xlabel(r"Outgoing Muon Energy [GeV]")
    ax4.set_ylabel("Count")
    plt.show()  
    plt.savefig(sample_type+"_events_outgoing_muon_energy_truth.png")                      