import matplotlib
import matplotlib.pyplot as plt
import h5py
import argparse
import numpy as np
import twoBytwo_defs
import mip_backgrounds
import auxiliary
import glob

nu_signal_pdg=-14
pion_pdg={111,211,-211}


def main(sim_dir, input_type):
    cc_dict, nc_dict, cc_primaries_dict, nc_primaries_dict = [dict() for i in range(4)]
    file_ctr=0

    file_ext = '' ### modified by commandline argument
    if input_type=='larnd': file_ext='.LARNDSIM.h5'
    elif input_type=='edep': file_ext='.EDEPSIM.h5'
    
    for sim_file in glob.glob(sim_dir+'/*'+file_ext):
        file_ctr+=1
        if file_ctr>100: break
#        if file_ctr%10==0: print('FILE #: ',file_ctr)
        print('FILE #: ',file_ctr)
        sim_h5 = h5py.File(sim_file,'r')
    
        ### partition file by spill
        unique_spill = np.unique(sim_h5['trajectories']['spillID'])
        for spill_id in unique_spill:
            ghdr, gstack, traj, vert, seg = auxiliary.get_spill_data(sim_h5, spill_id, input_type)

            ### partition by vertex ID within beam spill
            for v_i in range(len(vert['vertexID'])):
                vert_pos= [vert['x_vert'][v_i], vert['y_vert'][v_i], vert['z_vert'][v_i]]
                vert_in_active_LAr = twoBytwo_defs.fiducialized_vertex( vert_pos )

                ##### REQUIRE neutrino vertex in LAr active volume #####
                if vert_in_active_LAr==False: continue

                vert_id = vert['vertexID'][v_i]
            
                nu_mu_bar = auxiliary.signal_nu_pdg(ghdr, vert_id)
                is_cc = auxiliary.signal_cc(ghdr, vert_id)
                pionless = auxiliary.signal_meson_status(gstack, vert_id)
                fv_particle_origin=twoBytwo_defs.fiducialized_particle_origin(traj, vert_id)
                        
                ##### THRESHOLD BACKGROUNDS #####
                if is_cc==True and pionless==False and fv_particle_origin==True:
                    mip_backgrounds.pion_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, cc_dict)
#                    mip_backgrounds.primaries(spill_id, vert_id, ghdr, gstack, traj, vert, seg, cc_primaries_dict)

                ##### PID BACKGROUNDS #####
#                if is_cc==False and pionless==False and fv_particle_origin==True:
#                    mip_backgrounds.pion_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, nc_dict)
#                    mip_backgrouns.primaries(spill_id, vert_id, ghdr, gstack, traj, vert, seg, nc_primaries_dict)

    auxiliary.save_dict_to_json(cc_dict, 'cc_pion_backgrounds', True)
    auxiliary.save_dict_to_json(cc_primaries_dict, 'cc_primaries', True)
#    auxiliary.save_dict_to_json(nc_dict, 'nc_pion_backgrounds', True)
#    auxiliary.save_dict_to_json(nc_primaries_dict, 'nc_primaries', True)
    


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--sim_dir', default=None, required=True, type=str, \
                        help='''string corresponding to the path of the directory containing edep-sim or larnd ouput simulation file(s)''')
    parser.add_argument('-t', '--input_type', default='edep', choices=['edep', 'larnd'], type=str, \
                        help='''string corresponding to the output file type: edep or larnd''')
    args = parser.parse_args()
    main(**vars(args))
