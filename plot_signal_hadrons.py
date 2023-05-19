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
    elif sig_bkg == 3:
        sample_type = 'wrong_sign_bkg'
    else: 
        return "Error: plot_hadrons function given undefined signal/background definition"
    
    # PLOT: total visible energy + contained visible energy
    fig0, ax0 = plt.subplots(figsize=(8,4))
    data0tot = np.array([d[key]['total_edep'] for key in d.keys()])
    data0cont = np.array([d[key]['contained_edep'] for key in d.keys()])
    counts0tot, bins0tot = np.histogram(data0tot, bins=np.linspace(0,800,40))
    counts0cont, bins0cont = np.histogram(data0cont, bins=np.linspace(0,800,40))
    ax0.hist(bins0tot[:-1], bins=bins0tot, weights = counts0tot*scale_factor, label='total', histtype='step')
    ax0.hist(bins0cont[:-1], bins=bins0cont, weights = counts0cont*scale_factor, label='contained',histtype='step', linestyle='--')
    ax0.set_xlabel('Hadron Visible Energy [MeV]')
    ax0.set_ylabel('Count / 20 MeV')
    ax0.set_yscale('log')
    ax0.legend()
    ax0.grid(True)
    plt.savefig(sample_type+"_events_hadron_visible_energy.png")
    plt.close(fig0)

    # PLOT: hadron energy containment fraction
    fig1, ax1 = plt.subplots(figsize=(6,4))
    data1 = np.array([d[key]['contained_edep']/d[key]['total_edep'] for key in d.keys() if d[key]['total_edep']!=0])
    counts1, bins1 = np.histogram(data1, bins=np.linspace(0,1,20))
    ax1.hist(bins1[:-1], bins=bins1, weights = counts1*scale_factor, histtype='step')
    ax1.set_xlabel('Visible Hadron Energy Containment Fraction')
    ax1.set_ylabel('Count / 0.05')
    ax1.grid(True)       
    plt.savefig(sample_type+"_events_hadron_energy_containment_fraction.png")
    plt.close(fig1)

    # PLOT: hadron multiplicity
    fig2, ax2 = plt.subplots(figsize=(6,4))
    data2 = np.array([d[key]['hadron_mult'] for key in d.keys()])
    counts2, bins2 =np.histogram(data2, bins=np.linspace(0,25,26))
    ax2.hist(bins2[:-1], bins=bins2, weights = counts2*scale_factor, histtype='step')
    ax2.set_xlabel(r"Hadron Multiplicity")
    ax2.set_ylabel("Count / Hadron") 
    plt.savefig(sample_type+"_events_hadron_multiplicity_truth.png")
    plt.close(fig2)    

    # PLOT: truth-level 4-momentum squared of interaction
    #       ** no scale factor applied because we're looking at fractions anyways ** 
    fig3, ax3 = plt.subplots(figsize=(6,4))
    hadron_fs_pdg_list=[sorted(d[key]['hadron_pdg_set']) for key in d.keys()]
    hadron_fs_pdg_set=set(tuple(pdg) for pdg in hadron_fs_pdg_list)
    #print("Hadron PDG List:", hadron_fs_pdg_list)
    #print("Hadron PDG Set:", hadron_fs_pdg_set)
    hadron_fs_pdg_count=[(pdg_set, hadron_fs_pdg_list.count(list(pdg_set))) for pdg_set in hadron_fs_pdg_set]
    hadron_fs_pdg_fraction=[100*(i[1]/len(data2)) for i in hadron_fs_pdg_count]
    hadron_fs_pdg_labels=['+'.join(str(auxiliary.hadron_pdg_dict[j]) for j in i[0]) for i in hadron_fs_pdg_count]
    #print("Number of Events:", len(hadron_fs_pdg_list))
    #print("Hadron FS PDG Count:", hadron_fs_pdg_count)
    #print("Hadron FS PDG Fractions:", hadron_fs_pdg_fraction)
    #print("Hadron FS PDG Labels:", hadron_fs_pdg_labels)
    ax3.pie(hadron_fs_pdg_fraction, labels=hadron_fs_pdg_labels, autopct='%1.1f%%')
    ax3.set_title(r"Final State Hadrons in Signal Events")
    plt.savefig(sample_type+"_events_hadron_pdg_ids_truth.png")
    plt.close(fig3)    

    # PLOT: other hadron multiplicity
    fig4, ax4 = plt.subplots(figsize=(6,4))
    data4 = np.array([d[key]['other_had_mult'] for key in d.keys()])
    counts4, bins4 =np.histogram(data4, bins=np.linspace(0,10,11))
    ax4.hist(bins4[:-1], bins=bins4, weights = counts4*scale_factor, histtype='step')
    ax4.set_xlabel(r"Other Hadron Multiplicity")
    ax4.set_ylabel("Count / Other Hadron") 
    plt.savefig(sample_type+"_events_other_hadron_multiplicity_truth.png")
    plt.close(fig4)    

    # PLOT: neutron multiplicity
    fig5, ax5 = plt.subplots(figsize=(6,4))
    data5 = np.array([d[key]['neutron_mult'] for key in d.keys()])
    counts5, bins5 =np.histogram(data5, bins=np.linspace(0,20,21))
    ax5.hist(bins5[:-1], bins=bins5, weights = counts5*scale_factor, histtype='step')
    ax5.set_xlabel(r"Neutron Multiplicity")
    ax5.set_ylabel("Count / Neutron") 
    plt.savefig(sample_type+"_events_neutron_multiplicity_truth.png")
    plt.close(fig5)    

    # PLOT: proton multiplicity
    fig6, ax6 = plt.subplots(figsize=(6,4))
    data6 = np.array([d[key]['proton_mult'] for key in d.keys()])
    counts6, bins6 =np.histogram(data6, bins=np.linspace(0,20,21))
    ax6.hist(bins6[:-1], bins=bins6, weights = counts6*scale_factor, histtype='step')
    ax6.set_xlabel(r"Proton Multiplicity")
    ax6.set_ylabel("Count / Proton") 
    plt.savefig(sample_type+"_events_proton_multiplicity_truth.png")
    plt.close(fig6)   

    
    # PLOT: Fractions of Events with diff numbers of protons
    #       ** no scale factor applied because we're looking at fractions anyways ** 
    fig7, ax7 = plt.subplots(figsize=(6,4))
    p_mult_list = []
    total_np_events = 0
    for key in d.keys():
        if d[key]['other_had_mult']==0 and d[key]['proton_mult']>0 and d[key]['neutron_mult']>0:
            p_mult_list.append(d[key]['proton_mult'])
            total_np_events+=1
    if total_np_events >0:
        p_mult_count=[(mult_count, p_mult_list.count(mult_count)) for mult_count in np.arange(np.max(np.array(p_mult_list)))]
        #print("P mult count:", p_mult_count)
        p_mult_fraction=[100*(i[1]/total_np_events) for i in p_mult_count if i[1]>0]
        p_mult_labels=[str(i[0]) for i in p_mult_count if i[1]>0]
        #print("P mult labels:", p_mult_labels)
        ax7.pie(p_mult_fraction, labels=p_mult_labels, autopct='%1.1f%%')
        ax7.set_title(r"Proton Multiplicity in Signal Events with Neutrons and Protons")
        plt.savefig(sample_type+"_events_proton_mult_in_pn_events_truth.png")
    plt.close(fig7)    

    # PLOT: Max proton length for n p events
    fig8, ax8 = plt.subplots(figsize=(8,4))
    p_tot_lens = []
    p_cont_lens = []
    events_w_protons = 0
    for key in d.keys():
        if d[key]['proton_mult']>0:
            p_tot_lens.append(d[key]['max_p_total_length'])
            p_cont_lens.append(d[key]['max_p_contained_length'])
            events_w_protons+=1
    if events_w_protons >0:  
        data8tot = np.array(p_tot_lens)
        data8cont = np.array(p_cont_lens)
        counts8tot, bins8tot = np.histogram(data8tot, bins=np.linspace(0,60,60))
        counts8cont, bins8cont = np.histogram(data8cont, bins=np.linspace(0,60,60))
        ax8.hist(bins8tot[:-1], bins=bins8tot, weights = counts8tot*scale_factor, label='total', histtype='step')
        ax8.hist(bins8cont[:-1], bins=bins8cont, weights = counts8cont*scale_factor, label='contained',histtype='step', linestyle='--')
        ax8.set_xlabel(r"Length [cm]")
        ax8.set_title("Length of Longest Proton Track in Signal Events with Neutrons and Protons")
        ax8.set_ylabel("Count / cm") 
        ax8.legend()
        plt.savefig(sample_type+"_events_max_proton_length_in_pn_events_truth.png")   
    plt.close(fig8)