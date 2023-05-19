import matplotlib
import matplotlib.pyplot as plt
import h5py
import glob
import json
import argparse
import numpy as np
import twoBytwo_defs
import auxiliary
import signal_characterization as sig_char

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
    elif sig_bkg == 3:
        sample_type = 'wrong_sign_bkg'
    else: 
        return "Error: plot_muons function given undefined signal/background definition"
    
        
                                                  
    # PLOT: total visible energy + contained visible energy
    fig0, ax0 = plt.subplots(figsize=(8,4))
    data0tot = np.array([d[key]['total_edep'] for key in d.keys()])
    data0cont = np.array([d[key]['contained_edep'] for key in d.keys()])
    counts0tot, bins0tot = np.histogram(data0tot, bins=np.linspace(0,400,20))
    counts0cont, bins0cont = np.histogram(data0cont, bins=np.linspace(0,400,20))
    ax0.hist(bins0tot[:-1], bins=bins0tot, weights = counts0tot*scale_factor, label='total', histtype='step')
    ax0.hist(bins0cont[:-1], bins=bins0cont, weights = counts0cont*scale_factor, label='contained',histtype='step', linestyle='--')
    ax0.set_xlabel('Visible Muon Energy [MeV]')
    ax0.set_ylabel('Count / 20 MeV')
    #ax0.set_title(r'Muon Energy')
    ax0.legend()
    ax0.set_yscale('log')
    ax0.grid(True)
    plt.savefig(sample_type+"_events_muon_visible_energy.png")
    plt.close(fig0)

    # PLOT: muon energy containment fraction
    fig1, ax1 = plt.subplots(figsize=(6,4))
    data1 = np.array([d[key]['contained_edep']/d[key]['total_edep'] for key in d.keys() if d[key]['total_edep']!=0])
    counts1, bins1 = np.histogram(data1, bins=np.linspace(0,1,20))
    ax1.hist(bins1[:-1], bins=bins1, weights = counts1*scale_factor, histtype='step')
    ax1.set_xlabel('Visible Muon Energy Containment Fraction')
    ax1.set_ylabel('Count / 0.05')
    ax1.grid(True)       
    plt.savefig(sample_type+"_events_muon_containment_fraction.png")
    plt.close(fig1)    
    
    # PLOT: truth-level outgoing muon (lepton) angle 
    fig2, ax2 = plt.subplots(figsize=(6,4))
    data2 = np.array([d[key]['ang'] for key in d.keys()])
    counts2, bins2 =np.histogram(data2, bins=np.linspace(-0.45,0.75,61))
    ax2.hist(bins2[:-1], bins=bins2, weights = counts2*scale_factor, histtype='step')
    ax2.set_xlabel(r"Outgoing Muon Angle with Beam Direction")
    ax2.set_ylabel("Count / 0.02 Rad")
    plt.savefig(sample_type+"_events_outgoing_muon_angle_truth.png")
    plt.close(fig2)   

    # PLOT: truth-level outgoing muon (lepton) momentum 
    fig3, ax3 = plt.subplots(figsize=(6,4))
    data3 = np.array([d[key]['mom'] for key in d.keys()])/1000.
    counts3, bins3 =np.histogram(data3, bins=np.linspace(0,15,16))
    ax3.hist(bins3[:-1], bins=bins3, weights = counts3*scale_factor, histtype='step')
    ax3.set_xlabel(r"Outgoing Muon Momentum [GeV]")
    ax3.set_ylabel("Count / GeV")
    plt.savefig(sample_type+"_events_outgoing_muon_momentum_truth.png")
    plt.close(fig3)  

    # PLOT: truth-level 4-momentum squared of interaction
    fig4, ax4 = plt.subplots(figsize=(6,4))
    data4 = np.array([d[key]['q2'] for key in d.keys()]) / 1000000.
    counts4, bins4 =np.histogram(data4, bins=np.linspace(0,5,51))
    ax4.hist(bins4[:-1], bins=bins4, weights = counts4*scale_factor, histtype='step')
    ax4.set_xlabel(r"4-Momentum Transfer Squared [GeV$^2$]")
    ax4.set_ylabel(r"Count / 0.1 GeV$^2$") 
    plt.savefig(sample_type+"_events_qsq_truth.png")
    plt.close(fig4)  

    # PLOT: truth-level neutrino energy of interaction
    fig5, ax5 = plt.subplots(figsize=(6,4))
    data5 = np.array([d[key]['nu_energy'] for key in d.keys()]) / 1000.
    counts5, bins5 =np.histogram(data5, bins=np.linspace(0,15,16))
    ax5.hist(bins5[:-1], bins=bins5, weights = counts5*scale_factor, histtype='step')
    ax5.set_xlabel(r"Incident Neutrino Energy [GeV]")
    ax5.set_ylabel("Count / GeV") 
    plt.savefig(sample_type+"_events_nu_energy_truth.png")
    plt.close(fig5)      

    # PLOT: truth-level neutrino energy of interaction STACKED HIST BY END_PT_LOC
    loc_labels = [twoBytwo_defs.particle_end_loc_dict[k] for k in twoBytwo_defs.particle_end_loc_dict.keys()]
    data6f = []; data6u = []; data6d = []; data6b = []; data6s = []; data6p = []
    data7f = []; data7u = []; data7d = []; data7b = []; data7s = []; data7p = []
    data8f = []; data8u = []; data8d = []; data8b = []; data8s = []; data8p = []
    data9f = []; data9u = []; data9d = []; data9b = []; data9s = []; data9p = []
    for key in d.keys():
        if d[key]['end_pt_loc'] == 'f':
            data6f.append(d[key]['nu_energy'] / 1000.)
            data7f.append(d[key]['q2'] / 1000000.)
            data8f.append(d[key]['ang'])
            data9f.append(d[key]['mom'] / 1000.)
        elif d[key]['end_pt_loc'] == 'd':
            data6d.append(d[key]['nu_energy'] / 1000.)
            data7d.append(d[key]['q2'] / 1000000.)
            data8d.append(d[key]['ang'])
            data9d.append(d[key]['mom'] / 1000.)
        elif d[key]['end_pt_loc'] == 'b':
            data6b.append(d[key]['nu_energy'] / 1000.)
            data7b.append(d[key]['q2'] / 1000000.)
            data8b.append(d[key]['ang'])
            data9b.append(d[key]['mom'] / 1000.)
        elif d[key]['end_pt_loc'] == 's':
            data6s.append(d[key]['nu_energy'] / 1000.)
            data7s.append(d[key]['q2'] / 1000000.)   
            data8s.append(d[key]['ang'])
            data9s.append(d[key]['mom'] / 1000.)          
        elif d[key]['end_pt_loc'] == 'p':
            data6p.append(d[key]['nu_energy'] / 1000.)
            data7p.append(d[key]['q2'] / 1000000.)
            data8p.append(d[key]['ang'])
            data9p.append(d[key]['mom'] / 1000.)
        elif d[key]['end_pt_loc'] == 'u':
            data6u.append(d[key]['nu_energy'] / 1000.)
            data7u.append(d[key]['q2'] / 1000000.)
            data8u.append(d[key]['ang'])
            data9u.append(d[key]['mom'] / 1000.)

    fig6, ax6 = plt.subplots(figsize=(9,6))
    bins6 = np.linspace(0,15,16)
    counts6f, bins6f =np.histogram(np.array(data6f), bins=bins6)
    counts6d, bins6d =np.histogram(np.array(data6d), bins=bins6)
    counts6b, bins6b =np.histogram(np.array(data6b), bins=bins6)
    counts6s, bins6s =np.histogram(np.array(data6s), bins=bins6)
    counts6p, bins6p =np.histogram(np.array(data6p), bins=bins6)
    counts6u, bins6u =np.histogram(np.array(data6u), bins=bins6)
    ax6.hist((bins6f[:-1],bins6d[:-1],bins6b[:-1],bins6s[:-1],bins6p[:-1],bins6u[:-1]), bins=bins6, \
        weights = (counts6f*scale_factor,counts6d*scale_factor,counts6b*scale_factor,counts6s*scale_factor,counts6p*scale_factor,counts6u*scale_factor), histtype='bar', label=loc_labels, stacked='True')
    ax6.set_xlabel(r"Incident Neutrino Energy [GeV]")
    ax6.set_ylabel("Count / GeV") 
    ax6.legend(loc='upper right')
    ax6.set_title('Signal Event Neutrino Energy Spectrum by Muon Track End Behavior')
    plt.savefig(sample_type+"_events_nu_energy_truth_stacked_by_muon_end_loc.png")
    plt.close(fig6)   

    fig7, ax7 = plt.subplots(figsize=(9,6))
    bins7 = np.linspace(0,5,51)
    counts7f, bins7f =np.histogram(np.array(data7f), bins=bins7)
    counts7d, bins7d =np.histogram(np.array(data7d), bins=bins7)
    counts7b, bins7b =np.histogram(np.array(data7b), bins=bins7)
    counts7s, bins7s =np.histogram(np.array(data7s), bins=bins7)
    counts7p, bins7p =np.histogram(np.array(data7p), bins=bins7)
    counts7u, bins7u =np.histogram(np.array(data7u), bins=bins7)
    ax7.hist((bins7f[:-1],bins7d[:-1],bins7b[:-1],bins7s[:-1],bins7p[:-1],bins7u[:-1]), bins=bins7, \
        weights = (counts7f*scale_factor,counts7d*scale_factor,counts7b*scale_factor,counts7s*scale_factor,counts7p*scale_factor,counts7u*scale_factor), histtype='bar', label=loc_labels, stacked='True')
    ax7.legend(loc='upper right')
    ax7.set_ylabel(r"Count / 0.1 GeV$^2$") 
    ax7.set_xlabel(r"4-Momentum Transfer Squared [GeV$^2$]")
    ax7.set_title(r'Signal Event Q$^2$ by Muon Track End Behavior')
    plt.savefig(sample_type+"_events_qsq_truth_stacked_by_muon_end_loc.png") 
    plt.close(fig7) 

    fig8, ax8 = plt.subplots(figsize=(9,6))
    bins8 = np.linspace(-0.45,0.75,61)
    counts8f, bins8f =np.histogram(np.array(data8f), bins=bins8)
    counts8d, bins8d =np.histogram(np.array(data8d), bins=bins8)
    counts8b, bins8b =np.histogram(np.array(data8b), bins=bins8)
    counts8s, bins8s =np.histogram(np.array(data8s), bins=bins8)
    counts8p, bins8p =np.histogram(np.array(data8p), bins=bins8)
    counts8u, bins8u =np.histogram(np.array(data8u), bins=bins8)
    ax8.hist((bins8f[:-1],bins8d[:-1],bins8b[:-1],bins8s[:-1],bins8p[:-1],bins8u[:-1]), bins=bins8, \
        weights = (counts8f*scale_factor,counts8d*scale_factor,counts8b*scale_factor,counts8s*scale_factor,counts8p*scale_factor,counts8u*scale_factor), histtype='bar', label=loc_labels, stacked='True')
    ax8.legend(loc='upper right')
    ax8.set_title('Signal Event Outgoing Muon Angle by Muon Track End Behavior')
    ax8.set_xlabel(r"Outgoing Muon Angle with Beam Direction")
    ax8.set_ylabel("Count / 0.02 Rad")
    plt.savefig(sample_type+"_events_muon_angle_truth_stacked_by_muon_end_loc.png") 
    plt.close(fig8) 

    fig9, ax9 = plt.subplots(figsize=(9,6))
    bins9 = np.linspace(0,15,16)
    counts9f, bins9f =np.histogram(np.array(data9f), bins=bins9)
    counts9d, bins9d =np.histogram(np.array(data9d), bins=bins9)
    counts9b, bins9b =np.histogram(np.array(data9b), bins=bins9)
    counts9s, bins9s =np.histogram(np.array(data9s), bins=bins9)
    counts9p, bins9p =np.histogram(np.array(data9p), bins=bins9)
    counts9u, bins9u =np.histogram(np.array(data9u), bins=bins9)
    ax9.hist((bins9f[:-1],bins9d[:-1],bins9b[:-1],bins9s[:-1],bins9p[:-1],bins9u[:-1]), bins=bins9, \
        weights = (counts9f*scale_factor,counts9d*scale_factor,counts9b*scale_factor,counts9s*scale_factor,counts9p*scale_factor,counts9u*scale_factor), histtype='bar', label=loc_labels, stacked='True')
    ax9.legend(loc='upper right')
    ax9.set_title('Signal Event Outgoing Muon Momentum by Muon Track End Behavior')
    ax9.set_xlabel(r"Outgoing Muon Momentum [GeV]")
    ax9.set_ylabel("Count / GeV")
    plt.savefig(sample_type+"_events_outgoing_muon_momentum_truth_stacked_by_muon_end_loc.png")
    plt.close(fig9)  