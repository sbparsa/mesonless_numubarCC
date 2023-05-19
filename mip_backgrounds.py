import matplotlib
import matplotlib.pyplot as plt
import twoBytwo_defs
import auxiliary
import numpy as np


def primaries(spill_id, vert_id, ghdr, gstack, traj, vert, seg, primary_dict):
    traj_mask=traj['vertexID']==vert_id
    final_states=traj[traj_mask]

    track_ids_seen=set()
    for fs in final_states:
        if fs['parentID']!=-1: continue
        
        pdg = fs['pdgId']
        if pdg in [111, 211, -211]: continue

        track_id = fs['trackID']
        if track_id in track_ids_seen: continue

        track_id_set = auxiliary.same_pdg_connected_trajectories(pdg, track_id, \
                                                                 final_states, \
                                                                 traj, ghdr)
        track_ids_seen.update(track_id_set)

        total_edep, contained_edep, total_length, contained_length=[0. for i in range(4)]
        for tis in track_id_set:
            total_edep += auxiliary.total_energy(pdg, tis, traj, seg)
            contained_edep += auxiliary.fv_contained_energy(pdg, tis, traj, seg)
            total_length += auxiliary.total_length(pdg, tis, traj, seg)
            contained_length += auxiliary.fv_contained_length(pdg, tis, traj, seg)        
        
        primary_dict[(spill_id,vert_id, track_id)]=dict(
            pdg=int(pdg),
            total_edep=float(total_edep),
            contained_edep=float(contained_edep),
            total_length=float(total_length),
            contained_length=float(contained_length))
        
    return primary_dict



def pion_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, pion_dict):
    traj_vert_mask = traj['vertexID']==vert_id
    final_states = traj[traj_vert_mask]

    ghdr_mask=ghdr['vertexID']==vert_id
    nu_energy=ghdr[ghdr_mask]['Enu']
    q2=ghdr[ghdr_mask]['Q2']
    mom=ghdr[ghdr_mask]['lep_mom']
    ang=ghdr[ghdr_mask]['lep_ang']
    
    vtx_mask = vert['vertexID']==vert_id
    vtx_x=vert[vtx_mask]['x_vert']
    vtx_y=vert[vtx_mask]['y_vert']
    vtx_z=vert[vtx_mask]['z_vert']

    track_ids_seen=set()
    for fs in final_states:
        pdg = fs['pdgId']
        if pdg not in [111,211,-211]: continue
        
        track_id = fs['trackID']        
        if track_id in track_ids_seen: continue
        
        track_id_set = auxiliary.same_pdg_connected_trajectories(pdg, track_id, \
                                                                 final_states, \
                                                                 traj, ghdr)
        track_ids_seen.update(track_id_set)

        parent_pdg = auxiliary.find_parent_pdg(fs['parentID'], fs['vertexID'], traj, ghdr)

        total_edep, contained_edep, total_length, contained_length=[0. for i in range(4)]
        for tis in track_id_set:
            total_edep += auxiliary.total_energy(pdg, tis, traj, seg)
            contained_edep += auxiliary.fv_contained_energy(pdg, tis, traj, seg)
            total_length += auxiliary.total_length(pdg, tis, traj, seg)
            contained_length += auxiliary.fv_contained_length(pdg, tis, traj, seg)
            
        pion_dict[(spill_id,vert_id, track_id)]=dict(
            pdg=int(pdg),
            parent_pdg=int(parent_pdg),
            total_edep=float(total_edep),
            contained_edep=float(contained_edep),
            total_length=float(total_length),
            contained_length=float(contained_length),
            nu_energy=float(nu_energy),
            q2=float(q2),
            mom=float(mom),
            ang=float(ang),
            vtx_x=float(vtx_x),
            vtx_y=float(vtx_y),
            vtx_z=float(vtx_z))
        
    return



def plot_threshold_backgrounds(d, edep=True, contained_energy_ratio=True, \
                               parent_pdg=True, length=True, \
                               vertex_mult=True, spill_mult=True, save=True):
    if edep: pion_edep(d, save)
    if contained_energy_ratio: pion_contained_energy_ratio(d, save)    
    if parent_pdg: pion_parent_pdg(d, save)
    if length: pion_length(d, save)
    if vertex_mult: vertex_multiplicity(d, save)
    if spill_mult: spill_multiplicity(d, save)

    
    
def pion_edep(d, save):
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
    if save==True: plt.savefig("pion_edep.png")
    else: plt.show()



def pion_contained_energy_ratio(d, save):    
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
    if save==True: plt.savefig("pion_contained_edep_ratio.png")
    else: plt.show()



def pion_parent_pdg(d, save):    
    fig2, ax2 = plt.subplots(1,3,figsize=(12,4))
    piplus_parent_pdg_list=[d[key]['parent_pdg'] for key in d.keys() if d[key]['pdg']==211]
    piplus_parent_pdg_set=set(piplus_parent_pdg_list)
    piplus_particle_count=[(pdg, piplus_parent_pdg_list.count(pdg)) for pdg in piplus_parent_pdg_set]
    piplus_fraction=[100*(i[1]/len(piplus_parent_pdg_list)) for i in piplus_particle_count]
    piplus_pdg_label=[str(i[0]) for i in piplus_particle_count]
    print('pi+: ',piplus_particle_count)
    ax2[0].pie(piplus_fraction, labels=piplus_pdg_label, autopct='%1.1f%%')
    ax2[0].set_title(r'$\pi^+$')
    
    piminus_parent_pdg_list=[d[key]['parent_pdg'] for key in d.keys() if d[key]['pdg']==-211]
    piminus_parent_pdg_set=set(piminus_parent_pdg_list)
    piminus_particle_count=[(pdg, piminus_parent_pdg_list.count(pdg)) for pdg in piminus_parent_pdg_set]
    piminus_fraction=[100*(i[1]/len(piminus_parent_pdg_list)) for i in piminus_particle_count]
    piminus_pdg_label=[str(i[0]) for i in piminus_particle_count]
    print('pi-: ',piminus_particle_count)
    ax2[1].pie(piminus_fraction, labels=piminus_pdg_label, autopct='%1.1f%%')
    ax2[1].set_title(r'$\pi^-$')
    
    pi0_parent_pdg_list=[d[key]['parent_pdg'] for key in d.keys() if d[key]['pdg']==111]
    pi0_parent_pdg_set=set(pi0_parent_pdg_list)
    pi0_particle_count=[(pdg, pi0_parent_pdg_list.count(pdg)) for pdg in pi0_parent_pdg_set]
    pi0_fraction=[100*(i[1]/len(pi0_parent_pdg_list)) for i in pi0_particle_count]
    pi0_pdg_label=[str(i[0]) for i in pi0_particle_count]
    print('pi0: ',pi0_particle_count)
    ax2[2].pie(pi0_fraction, labels=pi0_pdg_label, autopct='%1.1f%%')
    ax2[2].set_title(r'$\pi^0$')
    if save==True: plt.savefig("pion_parent_pdg.png")
    else: plt.show()



def pion_length(d, save):
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
    if save==True: plt.savefig("pion_length.png")
    else: plt.show()



def vertex_multiplicity(d, save):
    fig4, ax4 = plt.subplots(figsize=(6,4))
    vertex_pion_mult, vertex_piplus_mult, vertex_piminus_mult, vertex_pi0_mult=[[] for i in range(4)]
    current_vertex_id=-1;
    pion_ctr=0; piplus_ctr=0; piminus_ctr=0; pi0_ctr=0
    for key in d.keys():
        if key[1]!=current_vertex_id:
            if current_vertex_id!=-1:
                vertex_pion_mult.append(pion_ctr); pion_ctr=0
                vertex_piplus_mult.append(piplus_ctr); piplus_ctr=0
                vertex_piminus_mult.append(piminus_ctr); piminus_ctr=0
                vertex_pi0_mult.append(pi0_ctr); pi0_ctr=0
            current_vertex_id=key[1]
        else:
            pion_ctr+=1
            if d[key]['pdg']==211: piplus_ctr+=1
            if d[key]['pdg']==-211: piminus_ctr+=1
            if d[key]['pdg']==111: pi0_ctr+=1
    bins=np.linspace(0,9,10)
    ax4.hist(vertex_pion_mult, bins=bins, label=r'All $\pi$', alpha=0.33)
    ax4.hist(vertex_piplus_mult, bins=bins, label=r'$\pi^+$', histtype='step', linewidth=2)
    ax4.hist(vertex_piminus_mult, bins=bins, label=r'$\pi^-$', histtype='step', linewidth=2)
    ax4.hist(vertex_pi0_mult, bins=bins, label=r'$\pi^0$', histtype='step', linewidth=2)
    ax4.grid(True)
    ax4.set_xlabel('Multiplicity')
    ax4.set_ylabel(r'Count per $\nu$ Interaction')
    ax4.legend()
    if save==True: plt.savefig("pion_vertex_multiplicity.png")
    else: plt.show()



def spill_multiplicity(d, save):
    fig5, ax5 = plt.subplots(figsize=(6,4))
    spill_pion_mult, spill_piplus_mult, spill_piminus_mult, spill_pi0_mult=[[] for i in range(4)]
    current_spill_id=-1;
    pion_ctr=0; piplus_ctr=0; piminus_ctr=0; pi0_ctr=0
    for key in d.keys():
        if key[0]!=current_spill_id:
            if current_spill_id!=-1:
                spill_pion_mult.append(pion_ctr); pion_ctr=0
                spill_piplus_mult.append(piplus_ctr); piplus_ctr=0
                spill_piminus_mult.append(piminus_ctr); piminus_ctr=0
                spill_pi0_mult.append(pi0_ctr); pi0_ctr=0
            current_spill_id=key[0]
        else:
            pion_ctr+=1
            if d[key]['pdg']==211: piplus_ctr+=1
            if d[key]['pdg']==-211: piminus_ctr+=1
            if d[key]['pdg']==111: pi0_ctr+=1
    bins=np.linspace(0,9,10)
    ax5.hist(spill_pion_mult, bins=bins, label=r'All $\pi$', alpha=0.33)
    ax5.hist(spill_piplus_mult, bins=bins, label=r'$\pi^+$', histtype='step', linewidth=2)
    ax5.hist(spill_piminus_mult, bins=bins, label=r'$\pi^-$', histtype='step', linewidth=2)
    ax5.hist(spill_pi0_mult, bins=bins, label=r'$\pi^0$', histtype='step', linewidth=2)
    ax5.grid(True)
    ax5.set_xlabel('Multiplicity')
    ax5.set_ylabel(r'Count per Beam Spill')
    ax5.legend()
    if save==True: plt.savefig("pion_spill_multiplicity.png")
    else: plt.show()
    
