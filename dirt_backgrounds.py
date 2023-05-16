import matplotlib
import matplotlib.pyplot as plt
import twoBytwo_defs
import numpy as np

def dirt_muon_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, dirt_dict):
    traj_vert_mask = traj['vertexID']==vert_id
    final_states = traj[traj_vert_mask]

    for fs in final_states:
        if fs['pdgId'] not in [13, -13]: continue # [111,211,-211]: continue

        pdg = fs['pdgId'] # *** pdg ***
        
        parent_id = fs['parentID']
        parent = final_states['trackID']==parent_id

        parent_pdg = final_states[parent]['pdgId'] # *** parent pdg ***

        if parent_id ==-1:
            ghdr_mask = ghdr['vertexID']==fs['vertexID']
            parent_pdg = ghdr[ghdr_mask]['nu_pdg']

        parent_pdg = parent_pdg.tolist()[0]
        print('check parent pdg: ', parent_pdg)
        
        track_id = fs['trackID']
        total_edep=0.; contained_edep=0.; total_length=0.; contained_length=0.

        #### checking dirt muons
        if abs(pdg)==13:
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

        if contained_edep>5:
            print(pdg,'\t',parent_pdg,'\t',total_edep,' MeV\t',contained_edep,' MeV\t', total_length,' cm\t',contained_length,' cm')


            dirt_dict[(spill_id,vert_id, track_id)]=dict(
                pdg=pdg,
                parent_pdg=parent_pdg,
                total_edep=total_edep,
                contained_edep=contained_edep,
                total_length=total_length,
                contained_length=contained_length)
        ## end if
    return


def plot_dirt_backgrounds(d):
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
    #plt.show()
    plt.savefig('vis_edep.png')

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
    plt.savefig('vis_edep_containment,png')

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
    plt.savefig('parent_pdg.png')

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
    plt.savefig('track_length')
