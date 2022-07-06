###################################################################################################
## 08/08/2020  Alex Alekseenko
##
## This code will read saved solutions produced by Fortran FFB code and will save each solution
##  as a file plottable in VisIt or Paraview
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

#########################
## First we will set path to a directory where a CLEAN solution is saved
######################################
path = 'E:\\SimulationsLearningCollision\\take3_M41\\good\\164'
#######################################
# open a file that contains information about the nodal points.
filename = 'E:\\SimulationsLearningCollision\\take3_M41\\good\\080\\CollTrnDta180_1su1sv1sw41MuUU41MvVU41MwWU_nodes.dat'
nodes_u, nodes_v, nodes_w, nodes_gwts = my_readwrite.read_nodes(filename)
#########################################################################

######################
print("checking directory for solution save files:", path)
names, timearry = my_readwrite.my_get_soltn_file_names_after_time_wtime(path, -0.1)
tt1 = names.sort()
tt2 = timearry.sort()
print("found {0} solution saves files".format(len(timearry)))
######################
n1 = nodes_u.shape[1]  # total number of nodal points
for iii in range(0, len(names) - 1):
    solarry_0, collarry_0, solsize = my_readwrite.my_read_sol_coll_trim(names[iii], MM, Mtrim)
    #### Open a .vtk to write ASCII data
    #### First, we need to get a name of the file.
    vtk_file_path = path + '\\plotsol_{0:8.8f}.vtp'.format(timearry[iii])
    save_file = open(vtk_file_path, 'w')
    ## Now the file is opened, we will write stuff into it.
    save_file.write('<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\"> \n')
    save_file.write('  <PolyData>\n')  # ints ' + 'I{0:02d}J{1:02d}K{2:02d}\n'.format(ii1,ii2,ii4))
    save_file.write('    <Piece NumberOfPoints=\"{0}\" NumberOfVerts=\"0\" '.format(n1) +
                    'NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n')
    save_file.write('      <Points>\n')
    save_file.write('        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n')
    data_line = ''
    for ii3 in range(n1):
        data_line = data_line + '{0:8.5f} '.format(nodes_u[0, ii3])
        data_line = data_line + '{0:8.5f} '.format(nodes_v[0, ii3])
        data_line = data_line + '{0:8.5f} '.format(nodes_w[0, ii3])
        data_line = data_line + '  '
        if (ii3 + 1) % 10 == 0 and ii3 > 0:
            save_file.write(data_line + '\n')
            data_line = ''
    save_file.write(data_line + '\n')
    save_file.write('        </DataArray>\n')
    save_file.write('      </Points>\n')
    save_file.write('      <PointData Scalars=\"solution\">\n')
    print('Solution is about to be recorded time={0}. max val={1}, min val={2}'.format(timearry[iii],
                                                                                       np.amax(solarry_0[0, :]),
                                                                                       np.amin(solarry_0[0, :])))
    save_file.write('        <DataArray type=\"Float32\" Name=\"solution\" format=\"ascii\">\n')
    data_line = ''
    for ii4 in range(0, n1):
        data_line = data_line + '{0:8.5f} '.format(
            solarry_0[0, ii4])  # this may will ran outside of array bounds if array is inconsistent
        if (ii4 + 1) % 20 == 0 and ii4 > 0:
            save_file.write(data_line + '\n')
            data_line = ''
    save_file.write(data_line + '\n')
    save_file.write('        </DataArray>\n')
    print('Solution is about to be recorded time={0}. max val={1} min val={2}'.format(timearry[iii],
                                                                                      np.amax(
                                                                                          np.absolute(solarry_0[0, :])),
                                                                                      np.amin(np.absolute(
                                                                                          solarry_0[0, :]))))
    save_file.write('        <DataArray type=\"Float32\" Name=\"solABS\" format=\"ascii\">\n')
    data_line = ''
    for ii4 in range(0, n1):
        if np.absolute(solarry_0[0, ii4]) > 0.0001:
            data_line = data_line + '{0:8.5f} '.format(np.absolute(solarry_0[0, ii4]))
        else:
            data_line = data_line + '{0:8.5f} '.format(0.0001)
        if (ii4 + 1) % 20 == 0 and ii4 > 0:
            save_file.write(data_line + '\n')
            data_line = ''
    save_file.write(data_line + '\n')
    save_file.write('        </DataArray>\n')
    save_file.write('      </PointData>\n')
    #############################
    save_file.write('    </Piece>\n')
    save_file.write('  </PolyData>\n')
    save_file.write('</VTKFile>\n')
    #### all done, the points are written. Close the file
    save_file.close()
    #####################################################################################################
    ###
    ### NEXT WE WRITE COLLISION SAVE IN A SEPARATE FILE.
    ###
    #####################################################################################################
    #### First, we need to get a name of the file.
    vtk_file_path = path + '\\plotcoll{0:8.8f}.vtp'.format(timearry[iii])
    save_file = open(vtk_file_path, 'w')
    ## Now the file is opened, we will write stuff into it.
    save_file.write('<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\"> \n')
    save_file.write('  <PolyData>\n')  # ints ' + 'I{0:02d}J{1:02d}K{2:02d}\n'.format(ii1,ii2,ii4))
    save_file.write('    <Piece NumberOfPoints=\"{0}\" NumberOfVerts=\"0\" '.format(n1) +
                    'NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n')
    save_file.write('      <Points>\n')
    save_file.write('        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n')
    data_line = ''
    for ii3 in range(n1):
        data_line = data_line + '{0:8.5f} '.format(nodes_u[0, ii3])
        data_line = data_line + '{0:8.5f} '.format(nodes_v[0, ii3])
        data_line = data_line + '{0:8.5f} '.format(nodes_w[0, ii3])
        data_line = data_line + '  '
        if (ii3 + 1) % 10 == 0 and ii3 > 0:
            save_file.write(data_line + '\n')
            data_line = ''
    save_file.write(data_line + '\n')
    save_file.write('        </DataArray>\n')
    save_file.write('      </Points>\n')
    save_file.write('      <PointData Scalars=\"solution\">\n')
    print('Solution is about to be recorded time={0}. max val={1}, min val={2}'.format(timearry[iii],
                                                                                       np.amax(collarry_0[0, :]),
                                                                                       np.amin(collarry_0[0, :])))
    save_file.write('        <DataArray type=\"Float32\" Name=\"colloperat\" format=\"ascii\">\n')
    data_line = ''
    for ii4 in range(0, n1):
        data_line = data_line + '{0:8.5f} '.format(
            collarry_0[0, ii4])  # this may will ran outside of array bounds if array is inconsistent
        if (ii4 + 1) % 20 == 0 and ii4 > 0:
            save_file.write(data_line + '\n')
            data_line = ''
    save_file.write(data_line + '\n')
    save_file.write('        </DataArray>\n')
    print('Solution is about to be recorded time={0}. max val={1} min val={2}'.format(timearry[iii],
                                                                                      np.amax(np.absolute(
                                                                                          collarry_0[0, :])), np.amin(
            np.absolute(collarry_0[0, :]))))
    save_file.write('        <DataArray type=\"Float32\" Name=\"collABS\" format=\"ascii\">\n')
    data_line = ''
    for ii4 in range(0, n1):
        if np.absolute(collarry_0[0, ii4]) > 0.00001:
            data_line = data_line + '{0:8.5f} '.format(np.absolute(collarry_0[0, ii4]))
        else:
            data_line = data_line + '{0:8.5f} '.format(0.0001)
        if (ii4 + 1) % 20 == 0 and ii4 > 0:
            save_file.write(data_line + '\n')
            data_line = ''
    save_file.write(data_line + '\n')
    save_file.write('        </DataArray>\n')
    save_file.write('      </PointData>\n')
    #############################
    save_file.write('    </Piece>\n')
    save_file.write('  </PolyData>\n')
    save_file.write('</VTKFile>\n')
    #### all done, the points are written. Close the file
    save_file.close()
############# ALL DONE
