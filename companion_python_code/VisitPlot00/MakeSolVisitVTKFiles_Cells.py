###################################################################################################
## 06/10/2021  Alex Alekseenko
##
## This code will read saved solutions produced by Fortran FFB code and will save each solution
##  as a file plottable in VisIt or Paraview. This subroutine will make a .vtk rectangular grid file
##
#####################################################################################################

### TRIMMING PARAMETERS
MM = 41
Mtrim = 0
###
err_step = 0.005
###

import numpy as np
import my_readwrite
import my_vtk_tools

#########################
## First we will set path to a directory where a CLEAN solution is saved
######################################
path = 'F:\\temp\\cleaned_100_SVZero_Runs\\runs06152021\\damporthogcomplimnt\\432\\k65'
#######################################
# open a file that contains information about the nodal points.
#filename = 'E:\\SimulationsLearningCollision\\take3_M41\\good\\080\\CollTrnDta180_1su1sv1sw41MuUU41MvVU41MwWU_nodes.dat'
#nodes_u, nodes_v, nodes_w, nodes_gwts = my_readwrite.read_nodes(filename)
#########################################################################

#######################################
# open a file that contains information about the nodal points.
filename = 'F:\\SimulationsLearningCollision\\take3_M41\\good\\080\\CollTrnDta180_1su1sv1sw41MuUU41MvVU41MwWU_grids.dat'
grids_cap_u, grids_cap_v, grids_cap_w, grids_u_tmp, grids_v_tmp, grids_w_tmp = my_readwrite.read_grids(filename)
#########################################################################
# Next we need to extract the data for the first grid.
#########################################################################
mu=grids_cap_u[0,0]
mv=grids_cap_v[0,0]
mw=grids_cap_w[0,0]
## next we truncate the grids arrays to only include the first grid
grids_u = grids_u_tmp[:,0:mu]
grids_v = grids_v_tmp[:,0:mu]
grids_w = grids_w_tmp[:,0:mu]
## now we can use these values in the subroutine that writes VTK files
######################
print("checking directory for solution save files:", path)
names, timearry = my_readwrite.my_get_soltn_file_names_after_time_wtime(path, -0.1)
tt1 = names.sort()
tt2 = timearry.sort()
print("found {0} solution saves files".format(len(timearry)))
######################
for iii in range(0, len(names) - 1):
    solarry_0, collarry_0, solsize = my_readwrite.my_read_sol_coll_trim(names[iii], MM, Mtrim)
    ####
    my_vtk_tools.writeVTKsolcollcellvalsRLG(solarry_0, collarry_0, grids_u, grids_v, grids_w,timearry[iii],path)
############# ALL DONE
