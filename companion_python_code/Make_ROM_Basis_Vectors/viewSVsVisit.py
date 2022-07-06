##############################
###  07/27/2020 A. Alekseenko
###  This is main program for making .vtp files to view with visit and powerview
##############################

### TRIMMING PARAMETERS
MM = 41
Mtrim = 0

#########################
### Read pickled SVD

import numpy as np
import pickle
import my_readwrite
import my_vtk_tools
path = "F:\\temp\\cleaned_SVs_100_PlusMaxwell"
### There is a file with all singular vectors and with just the first 100
### Please pay attention that the trimming is correct
open_full= False
if open_full:
    save_file = open(path+'M41_Trim7_Sol_SVD.dat', 'rb')
    Vh = pickle.load(save_file)
    s = pickle.load(save_file)
    U = pickle.load(save_file)
    save_file.close()
    save_file = open(path+'M41_Trim7_Sol_SVD_100.dat', 'wb')
    pickle.dump(Vh[0:100, :], save_file)
    pickle.dump(s[0:100], save_file)
    pickle.dump(U[:,0:100], save_file)
    save_file.close()

########################
else:
    save_file = open(path + '\\SpHom100Trim0SVs_NoMxwellPickle.dat', 'rb')
#    Vh = pickle.load(save_file)   # we have files that only contain singular vectors/bases
#    s = pickle.load(save_file)
    svect = pickle.load(save_file)
    save_file.close()
#    for ii in range(100):
#        print(s[ii])

###########################################################
## In this piece we will create a .vtp file that contains 3D mesh and 3D singular vectors.
##
###########################################################
make_svplot_file = True
sv_start = 40
sv_end = 60
if make_svplot_file:
    n1=MM**3
    ##### read the nodal points #############
    filename = 'F:\\SimulationsLearningCollision\\take3_M41\\good\\080\\CollTrnDta180_1su1sv1sw41MuUU41MvVU41MwWU_nodes.dat'
    nodes_u, nodes_v, nodes_w, nodes_gwts = my_readwrite.read_nodes(filename)
    #### Open a .vtk to write ASCII data
    #### First, we need to get a name of the file.
    vtk_file_path = path + '\\SpHom100NoMxwell_sv{0:02d}_{1:02d}_plotting.vtp'.format(sv_start, sv_end)
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
        data_line = data_line + '{0:8.5f} '.format(nodes_u[0,ii3])
        data_line = data_line + '{0:8.5f} '.format(nodes_v[0,ii3])
        data_line = data_line + '{0:8.5f} '.format(nodes_w[0,ii3])
        data_line = data_line + '  '
        if (ii3 + 1) % 10 == 0 and ii3 > 0:
            save_file.write(data_line + '\n')
            data_line = ''
    save_file.write(data_line + '\n')
    save_file.write('        </DataArray>\n')
    save_file.write('      </Points>\n')
    save_file.write('      <PointData Scalars=\"SingVect001\">\n')
    for ii3 in range(sv_start,sv_end):
        print('num of SV={0} max val={1} min val={2}'.format(ii3,np.amax(svect[:,ii3]),np.amin(svect[:,ii3])))
        save_file.write('        <DataArray type=\"Float32\" Name=\"SingVect{0:03d}\" format=\"ascii\">\n'.format(ii3+1))
        data_line = ''
        for ii4 in range(0, n1):
            data_line = data_line + '{0:8.5f} '.format(svect[ii4,ii3])
            if (ii4 + 1) % 20 == 0 and ii4 > 0:
               save_file.write(data_line + '\n')
               data_line = ''
        save_file.write(data_line + '\n')
        save_file.write('        </DataArray>\n')
    for ii3 in range(sv_start, sv_end):
        print('num of SV={0} max val={1} min val={2}'.format(ii3, np.amax(np.absolute(svect[:, ii3])), np.amin(np.absolute(svect[:, ii3]))))
        save_file.write(
            '        <DataArray type=\"Float32\" Name=\"SingVectABS{0:03d}\" format=\"ascii\">\n'.format(ii3 + 1))
        data_line = ''
        for ii4 in range(0, n1):
            if np.absolute(svect[ii4, ii3])>0.0001:
                data_line = data_line + '{0:8.5f} '.format(np.absolute(svect[ii4, ii3]))
            else:
                data_line = data_line + '{0:8.5f} '.format(0.0001)
            if (ii4 + 1) % 20 == 0 and ii4 > 0:
                save_file.write(data_line + '\n')
                data_line = ''
        save_file.write(data_line + '\n')
        save_file.write('        </DataArray>\n')
    save_file.write('      </PointData>\n')
    ####### Next we will record trajectories -- lines
    #############################
    #############################
    save_file.write('    </Piece>\n')
    save_file.write('  </PolyData>\n')  # ints ' + 'I{0:02d}J{1:02d}K{2:02d}\n'.format(ii1,ii2,ii4))
    save_file.write('</VTKFile>\n')
    #### all done, the points are written. Close the file
    save_file.close()