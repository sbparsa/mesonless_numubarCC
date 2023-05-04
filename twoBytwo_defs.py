import numpy as np

#### geometry definitions copied straight from https://github.com/DUNE/2x2_sim/blob/main/validation/edepsim_validation.py

NDHallwidths = [1000.,550.,2000.] # cm

def tpc_bounds(i):
            """A sad little function that returns the bounds of each 2x2 tpc in one dimension.
            The dimension is chosen by i: 0, 1, 2 -> x, y, z.
            Values are taken from 2x2_sim/run-edep-sim/geometry/Merged2x2MINERvA_v2"""
            
            active_tpc_widths = [30.6, 130., 64.] # cm
            
            # The positions in cm of the center of each tpc relative to a module center.
            # There are two tpcs for each module.            
            tpcs_relative_to_module = [[-15.7,0.,0.], [15.7, 0., 0.]]

            # The positions in cm of each of the four modules, relative to the 2x2 center position.
            modules_relative_to_2x2= [[-33.5,0.,-33.5],
                                      [33.5,0.,-33.5],
                                      [-33.5,0.,33.5],
                                      [33.5,0.,33.5]]
            
            # The position of the 2x2 center, relative to the center of the ND hall
            detector_center = [0.,52.25,0.]
        
            # Get the tpc bounds relative to the tpc center in the ith coordinates
            tpc_bounds = np.array([-active_tpc_widths[i]/2., active_tpc_widths[i]/2.])
            
            tpc_bounds_relative_to_2x2 = []
            for tpc in tpcs_relative_to_module:
                tpc_bound_relative_to_module = tpc_bounds + tpc[i]
                for module in modules_relative_to_2x2:
                    bound = tpc_bound_relative_to_module + module[i]
                    tpc_bounds_relative_to_2x2.append(bound)
                    
            bounds_relative_to_NDhall = np.array(tpc_bounds_relative_to_2x2) + detector_center[i]
            
            return np.unique(bounds_relative_to_NDhall, axis = 0)


def MINERvA_bounds(i):
            """A sadder littler function that returns the bounds of the MINERvA detector for a given
            dimension i: 0, 1, 2 -> x, y, z.
            For now, I take the detector to just simply be two monolithic hexagonal prisms,
            downstream and upstream of the 2x2 modules. 
            """
            
            # Taken from the gdml file.
            MINERvA_center = [0., 43., -654.865]
                        
            # From GDML file, the length of one side of the outer detector in cm
            side_length = 199.439876988864
            
            # Properties of a hexagon, given a side length
            long_diameter = 2*side_length
            short_diameter = np.sqrt(3)/2.*long_diameter 
            
            width = short_diameter 
            height = long_diameter
            
            # The detector bounds, in each dimention. The bounds of the z dimension for each prism
            # were obtained by looking at the central positions of the first and last HCal frame 
            # position for each slab. 
            detector_bounds = np.array([ 
                [[-width/2., width/2.], [-height/2., height/2.], [403.27185, 451.18895]],
                [[-width/2., width/2.], [-height/2., height/2.], [858.54105, 997.2314]]
            ])
            
            bounds_relative_to_NDhall = [] 
            for bound in detector_bounds:
                
                bounds_relative_to_NDhall.append(bound[i] + MINERvA_center[i])
            
            return np.unique(np.array(bounds_relative_to_NDhall), axis = 0)
