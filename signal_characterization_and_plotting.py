import matplotlib
import matplotlib.pyplot as plt
import twoBytwo_defs
import numpy as np


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
        if fs['parentID']==-1:
            ghdr_vertexID_mask=ghdr['vertexID']==fs['vertexID']
            parent_pdg=ghdr[ghdr_vertexID_mask]['nu_pdg']
        else:
            parent_trackID_mask = traj['trackID']==fs['parentID']
            parent_pdg=traj[parent_trackID_mask]['pdgId']
        parent_pdg = parent_pdg.tolist()[0] # *** parent pdg ***
        
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
     #   end_pt = fs['xyz_end']
     #   end_pt_loc = ''
#
     #   if twoBytwo_defs.fiducialized_vertex(end_pt):
     #       end_pt_loc = 'Contained in 2x2'
     #   elif twoBytwo_defs.MINERvA_vertex(end_pt):
     #       end_pt_loc = 'Contained in MINERvA'
     #   else:
#
     #       x_MINERvA = twoBytwo_defs.MINERvA_bounds(0)
     #       y_MINERvA = twoBytwo_defs.MINERvA_bounds(1)
     #       z_MINERvA = twoBytwo_defs.MINERvA_bounds(2)
#
     #       x_tpc = twoBytwo_defs.tpc_bounds(0)
     #       y_tpc = twoBytwo_defs.tpc_bounds(1)
     #       z_tpc = twoBytwo_defs.tpc_bounds(2)
#
     #       # Check whether endpoint Z is between 2x2 and MINERvA
     #       if end_pt[2]<z_MINERvA[0] and end_pt[2]>z_tpc[1]:
     #           if end_pt[0]>x_MINERvA[0] and end_pt[0]<x_MINERvA[1] and end_pt[1]>y_MINERvA[0] and end_pt[1]<y_MINERvA[1]:
     #               end_pt_loc = 'Stops between 2x2 and MINERvA'
     #           else:
     #               end_pt_loc = 'Exits out the side of 2x2'
     #       # Check leaving back or side of MINERvA
     #       else:
     #           start_pt = fs['xyz_start']
     #           traj_vector = end_pt - start_pt
     #           if end_pt[2]>z_MINERvA[0] and end_pt[2]<z_MINERvA[1]:
#
     #               # if track goes through MINERvA vs. not through
     #           #elif track ends beyond MINERVA z bounds
     #           else:
     #               end_pt_loc = 'Undefined track endpoint location'
#
     #                                                                                                                                     
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
        q2 = float(q2))
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
        hadron_pdg = gstack_vert_fs_hadrons,
        hadron_pdg_set = gstack_vert_fs_pdg_set,
        total_edep=total_edep,
        contained_edep=contained_edep,
        total_length=total_length,
        contained_length=contained_length)
    return
                                      
# PLOT: Muon kinematics
#       sig_bkg is an int such that 0 == signal, 1 == 'dirt' backgrounds, 2 == 'beam' backgrounds
def plot_muons(d, scale_factor, sig_bkg = 0):
    
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
    fig0, ax0 = plt.subplots(figsize=(8,4))
    data0tot = np.array([d[key]['total_edep'] for key in d.keys()])
    data0cont = np.array([d[key]['contained_edep'] for key in d.keys()])
    counts0tot, bins0tot = np.histogram(data0tot, bins=np.linspace(0,400,20))
    counts0cont, bins0cont = np.histogram(data0cont, bins=np.linspace(0,400,20))
    ax0.hist(bins0tot[:-1], bins=bins0tot, weights = counts0tot*scale_factor, label='total', histtype='step')
    ax0.hist(bins0cont[:-1], bins=bins0cont, weights = counts0cont*scale_factor, label='contained',histtype='step')
    ax0.set_xlabel('Visible Muon Energy [MeV]')
    ax0.set_ylabel('Count / 20 MeV')
    #ax0.set_title(r'Muon Energy')
    ax0.legend()
    ax0.grid(True)
    plt.savefig(sample_type+"_events_muon_visible_energy.png")

    # PLOT: muon energy containment fraction
    fig1, ax1 = plt.subplots(figsize=(6,4))
    data1 = np.array([d[key]['contained_edep']/d[key]['total_edep'] for key in d.keys() if d[key]['total_edep']!=0])
    counts1, bins1 = np.histogram(data1, bins=np.linspace(0,1,20))
    ax1.hist(bins1[:-1], bins=bins1, weights = counts1*scale_factor, histtype='step')
    ax1.set_xlabel('Visible Muon Energy Containment Fraction')
    ax1.set_ylabel('Count / 0.05')
    ax1.grid(True)       
    plt.savefig(sample_type+"_events_muon_containment_fraction.png")    
    
    # PLOT: truth-level outgoing muon (lepton) angle 
    fig2, ax2 = plt.subplots(figsize=(6,4))
    data2 = np.cos(np.array([d[key]['ang'] for key in d.keys()]))
    counts2, bins2 =np.histogram(data2, bins=np.linspace(0.75,1,51))
    ax2.hist(bins2[:-1], bins=bins2, weights = counts2*scale_factor, histtype='step')
    ax2.set_xlabel(r"Cosine of Outgoing Muon Angle with Beam Direction")
    ax2.set_ylabel("Count / 0.005")
    plt.savefig(sample_type+"_events_outgoing_muon_angle_truth.png")   

    # PLOT: truth-level outgoing muon (lepton) momentum 
    fig3, ax3 = plt.subplots(figsize=(6,4))
    data3 = np.array([d[key]['mom'] for key in d.keys()])/1000.
    counts3, bins3 =np.histogram(data3, bins=np.linspace(0,30,31))
    ax3.hist(bins3[:-1], bins=bins3, weights = counts3*scale_factor, histtype='step')
    ax3.set_xlabel(r"Outgoing Muon Momentum [GeV]")
    ax3.set_ylabel("Count / GeV")
    plt.savefig(sample_type+"_events_outgoing_muon_momentum_truth.png")  

    # PLOT: truth-level 4-momentum squared of interaction
    fig4, ax4 = plt.subplots(figsize=(6,4))
    data4 = np.array([d[key]['q2'] for key in d.keys()]) / 1000000.
    counts4, bins4 =np.histogram(data4, bins=np.linspace(0,4,41))
    ax4.hist(bins4[:-1], bins=bins4, weights = counts4*scale_factor, histtype='step')
    ax4.set_xlabel(r"4-Momentum Transfer Squared [GeV$^2$]")
    ax4.set_ylabel(r"Count / 0.1 GeV$^2$") 
    plt.savefig(sample_type+"_events_qsq_truth.png")  

    # PLOT: truth-level neutrino energy of interaction
    fig5, ax5 = plt.subplots(figsize=(6,4))
    data5 = np.array([d[key]['nu_energy'] for key in d.keys()]) / 1000.
    counts5, bins5 =np.histogram(data5, bins=np.linspace(0,30,31))
    ax5.hist(bins5[:-1], bins=bins5, weights = counts5*scale_factor, histtype='step')
    ax5.set_xlabel(r"Incident Neutrino Energy [GeV]")
    ax5.set_ylabel("Count / GeV") 
    plt.savefig(sample_type+"_events_nu_energy_truth.png")                       

# PLOT: Hadron kinematics
#       sig_bkg is an int such that 0 == signal, 1 == 'dirt' backgrounds, 2 == 'beam' backgrounds
def plot_hadrons(d, scale_factor, sig_bkg = 0):
    
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
    fig0, ax0 = plt.subplots(figsize=(8,4))
    data0tot = np.array([d[key]['total_edep'] for key in d.keys()])
    data0cont = np.array([d[key]['contained_edep'] for key in d.keys()])
    counts0tot, bins0tot = np.histogram(data0tot, bins=np.linspace(0,800,40))
    counts0cont, bins0cont = np.histogram(data0cont, bins=np.linspace(0,800,40))
    ax0.hist(bins0tot[:-1], bins=bins0tot, weights = counts0tot*scale_factor, label='total', histtype='step')
    ax0.hist(bins0cont[:-1], bins=bins0cont, weights = counts0cont*scale_factor, label='contained',histtype='step')
    ax0.set_xlabel('Hadron Visible Energy [MeV]')
    ax0.set_ylabel('Count / 20 MeV')
    ax0.legend()
    ax0.grid(True)
    plt.savefig(sample_type+"_events_hadron_visible_energy.png")

    # PLOT: hadron energy containment fraction
    fig1, ax1 = plt.subplots(figsize=(6,4))
    data1 = np.array([d[key]['contained_edep']/d[key]['total_edep'] for key in d.keys() if d[key]['total_edep']!=0])
    counts1, bins1 = np.histogram(data1, bins=np.linspace(0,1,20))
    ax1.hist(bins1[:-1], bins=bins1, weights = counts1*scale_factor, histtype='step')
    ax1.set_xlabel('Visible Hadron Energy Containment Fraction')
    ax1.set_ylabel('Count / 0.05')
    ax1.grid(True)       
    plt.savefig(sample_type+"_events_hadron_energy_containment_fraction.png")

    # PLOT: truth-level 4-momentum squared of interaction
    fig2, ax2 = plt.subplots(figsize=(6,4))
    data2 = np.array([d[key]['hadron_mult'] for key in d.keys()])
    counts2, bins2 =np.histogram(data2, bins=np.linspace(0,25,26))
    ax2.hist(bins2[:-1], bins=bins2, weights = counts2*scale_factor, histtype='step')
    ax2.set_xlabel(r"Hadron Multiplicity")
    ax2.set_ylabel("Count / Hadron") 
    plt.savefig(sample_type+"_events_hadron_multiplicity_truth.png")    

    # PLOT: truth-level 4-momentum squared of interaction
    #       ** no scale factor applied because we're looking at fractions anyways ** 
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
    