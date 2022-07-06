##############################
###  05/29/2020 A. Alekseenko
###  This is main program for studying SVD of solutions
### This requires Python 2.7
##############################
import numpy as np

##############################
# Import visit's Python interface
#import sys
#sys.path.append("C:\\Program Files\\LLNL\\VisIt 3.1.2\\lib\\site-packages")

#sys.path.append("C:\\Program Files\\LLNL\\VisIt 3.1.2\\lib")
#import visit_writer

#import visit
#visit.Launch()

#
###############################

### TRIMMING PARAMETERS
MM = 41
Mtrim = 7

#########################
### Read pickled SVD

import pickle
import my_readwrite
import my_vtk_tools

### There is a file with all singular vectors and with just the first 100
### Please pay attention that the trimming is correct
open_full= False
if open_full:
    save_file = open('E:\\SimulationsLearningCollision\\StSol0\\M41_Trim7_Sol_SVD.dat', 'rb')
    Vh = pickle.load(save_file)
    s = pickle.load(save_file)
    U = pickle.load(save_file)
    save_file.close()
    save_file = open('E:\\SimulationsLearningCollision\\StSol0\\M41_Trim7_Sol_SVD_100.dat', 'wb')
    pickle.dump(Vh[0:100, :], save_file)
    pickle.dump(s[0:100], save_file)
    pickle.dump(U[:,0:100], save_file)
    save_file.close()

########################
else:
    save_file = open('E:\\SimulationsLearningCollision\\StSol0\\M41_Trim7_Sol_SVD_100.dat', 'rb')
    Vh = pickle.load(save_file)
    s = pickle.load(save_file)
    svect = pickle.load(save_file)
    save_file.close()
    print(s[0:100])

###########################################################
## In this piece we will create a directory in which runs of solution from beginning to the steady state
## will be summarized. Each saved soliution will be converted to truncated SVD basis pepresentation and an
## array will be created with the solution. We will refer to this array as the trajectory of the solution in
## the SVD basis. Additional iformation in the array will be the time and entropy of the VDF for each saved solution
##
###########################################################
create_pickled_SVD_solutions_arrays = False
if create_pickled_SVD_solutions_arrays:
      #####################################
      ## First we will prepare list of directories that contain single runs of solutions
      ######################################
      # path = 'G:\\SimulationsLearningCollision\\take2_accurate\\'
      path = 'E:\\SimulationsLearningCollision\\take3_M41'

      dirs = my_readwrite.my_get_dirs_names(path)

      #####################################
      ## select a subset of solutions randomly
      #####################################
      import random

      num_runs_avail = len(dirs)  # total available runs
      use_shorter_sample = False
      if use_shorter_sample:
        num_runs = 10  # number of runs to select
        runs_dir_names = []  # prepare an empty array that will keep names of the solutions to study
        #### GIVE IT A RUN
        for x in range(1000):
           ii = random.randint(0, num_runs_avail - 1)
        ######################

        for x in range(num_runs):
           ii = random.randint(0, num_runs_avail - 1)
           runs_dir_names.append(dirs[ii])
        ######################
      else:
        runs_dir_names = dirs
        num_runs = num_runs_avail  # number of runs to select
      ######################
      # We need some prep work. We will need to get nodes and weights of the quadratures associated with the solution loaded.
      # these will be used to evaluate the solution entropy.
      # these points can be obtained from the DG-Boltzmann solution saves of mesh parameters:
      filename = 'E:\\SimulationsLearningCollision\\take3_M41\\080\\CollTrnDta180_1su1sv1sw41MuUU41MvVU41MwWU_nodes.dat'
      nodes_u, nodes_v, nodes_w, nodes_gwts = my_readwrite.read_nodes(filename)
      #########################################################################

      #########################################################################
      #  Prepare rad decay filter factor array
      nodes_u_trim, new_size = my_readwrite.solution_trim(nodes_u, MM, Mtrim)
      nodes_v_trim, new_size = my_readwrite.solution_trim(nodes_v, MM, Mtrim)
      nodes_w_trim, new_size = my_readwrite.solution_trim(nodes_w, MM, Mtrim)
      nodes_gwts_trim, new_size = my_readwrite.solution_trim(nodes_gwts, MM, Mtrim)
      #########################################################################
      Mat_Moments = np.zeros((new_size, 5))
      Mat_Moments[:, 0] = nodes_gwts_trim[0, :]
      Mat_Moments[:, 1] = nodes_gwts_trim[0, :] * nodes_u_trim[0, :]
      Mat_Moments[:, 2] = nodes_gwts_trim[0, :] * nodes_v_trim[0, :]
      Mat_Moments[:, 3] = nodes_gwts_trim[0, :] * nodes_w_trim[0, :]
      scrp_array = (nodes_u_trim ** 2 + nodes_v_trim ** 2 + nodes_w_trim ** 2)
      Mat_Moments[:, 4] = nodes_gwts_trim[0, :] * scrp_array[0, :]
      ############################################################################
      ZZ =np.matmul(svect.T, Mat_Moments)
      u,s,vh = np.linalg.svd(ZZ,full_matrices=False)
      Q_m = u[:, 0:4]
      ############################################################################
      atemp = np.ones((1, new_size)) * .0000001  # technical array to replace small and zero values to evaluate log
      ######################
      # Next we are going to create an array where solution will be encoded into SVD coordinates. The time stamp will be added to the solution
      # Also, solution's entropy will be computed and added as additional entry
      # FORMAT of the array First index: goes along the time saves. Second intex: goes over
      # singular vectors and entropy and time, in that sequence
      #######################
      n2 = svect.shape[1]
      import os
      for ii in range(num_runs):
          ####################################
          ### From all available testing solutions get an entire spatially homogeneous relaxation run
          ### now we get file names into an array
          names, timearry = my_readwrite.my_get_soltn_file_names_after_time_wtime(runs_dir_names[ii], -0.1)
          tt1 = names.sort()
          tt2 = timearry.sort()
          svectT = svect.T
          ### Next we load the solutions from these files into a solution array.
          if len(names)>0:
             bigSVDarray = np.zeros((len(names), n2 + 2))  # created an empty array
             CalcGrad = np.zeros((len(names), n2))
             VarGrad = np.zeros((len(names), n2))
             bigColSVDarry = np.zeros((len(names), n2))
             for ii1 in range(len(names)):
                 solarry, collarry, solsize = my_readwrite.my_read_sol_coll_trim(names[ii1], MM, Mtrim)  # this is the first file
                 bigSVDarray[ii1,0:n2] = np.matmul(solarry, svect)  # created an empty array
                 bigColSVDarry[ii1, 0:n2] = np.matmul(collarry, svect)  # created an empty array
                 bigSVDarray[ii1,n2] = np.sum(solarry[0,:]*np.log(np.fmax(solarry[0,:],atemp))*nodes_gwts_trim[0,:])
                 bigSVDarray[ii1,n2+1] = timearry[ii1]
                 for ii2 in range(n2):
                    CalcGrad[ii1,ii2] = np.sum(svectT[ii2,:]*(1+np.log(np.fmax(solarry[0,:],atemp)))*nodes_gwts_trim[0,:])
                 VarGrad[ii1,:] = CalcGrad[ii1,:] - np.matmul(Q_m,np.matmul(Q_m.T,CalcGrad.T[:,ii1]))
                 ##############################################################
                 ## Get values of moments for the initial data
                 ##############################################################
                 if ii1==0:
                     moment_vector_0 = np.matmul(solarry, Mat_Moments)
                     #print(np.matmul(Q_m.T,VarGrad.T),np.matmul(Q_m.T,Q_m))   ## DEBUG PRINTOUT
             #######################################
             ## the data array was created. Now we need to record it to a pickle file.
             ###### SORTING BAD AND GOOD DATA based on conservation properties of the solution
             ###### some solutions violate the conservation laws a lot because fo the
             ###### numerical errors in the code. THese solutions will be sorted in a BAD folder.
             ###### other solutions will be sorted into a good folder.
             ###### sollarry must have the last solution in the bunch, the steady state. This depends whether the sorting trick above worked.
             moment_vector = np.matmul(solarry, Mat_Moments)   # compute the concerved moments of the solution. Next check concervation
             err0mom = np.abs(moment_vector[0,0]-moment_vector_0[0,0])
             err1mom = np.abs(moment_vector[0,1]-moment_vector_0[0,1])
             err2mom = np.abs(moment_vector[0,2]-moment_vector_0[0,2])
             err3mom = np.abs(moment_vector[0,3]-moment_vector_0[0,3])
             err4mom = np.abs(moment_vector[0,4]-moment_vector_0[0,4])
             ###########################################################################################
             ###### ATTENTION: HARD CODED CONSTANTS, this may change on a different class of solutions
             ############################################################################################
             if err0mom<0.005 and err1mom<0.005 and err2mom<0.005 and err3mom<0.005 and err4mom<0.01:
                 newpath = runs_dir_names[ii].replace('take3_M41', 'take3_M41_SVDsols\\good', 1)
             else:
                 newpath = runs_dir_names[ii].replace('take3_M41', 'take3_M41_SVDsols\\bad', 1)
             ##### END SORTING solutions
             os.mkdir(newpath)
             scrlen = len(runs_dir_names[ii])
             svd_file_prefix = names[0][scrlen:scrlen+14]
             svd_file_path = newpath+svd_file_prefix+'svd_SolET_array.dat'
             save_file = open(svd_file_path, 'wb')
             pickle.dump(bigSVDarray, save_file)
             pickle.dump(bigColSVDarry, save_file)
             pickle.dump(CalcGrad, save_file)
             pickle.dump(VarGrad, save_file)
             save_file.close()

####################################################################
# In this piece we will create a .vtk file where different 2D data will
# be stored ready for plotting. The FIELD VTK seemingly is not supported by VisIt
#####################################################################
create_vtk_files_trajectories_FIELD=False
if create_vtk_files_trajectories_FIELD:
    #first we will get the list of all directories where there is pickled dat
    # that can be used to create such vtk files.
    path = 'F:\\SimulationsLearningCollision\\take3_M41_SVDsols'
    dirs = my_readwrite.my_get_dirs_names(path)
    #Now we will call a subroutine that will scroll through all pickled
    # SVD encoded trajectories and will create a bunch of files and place them in path/VisItData/ to
    # be imported by visit.
    for ii in range(len(dirs)):
        files = my_vtk_tools.my_get_svdsols_file_names(dirs[ii])
        save_file = open(files[0], 'rb')
        BigArray = pickle.load(save_file)
        save_file.close()
        #BigArray contains SVD-encoded trajectories. We will now take a 2D slice/projection of a trajectory and record it into
        # a .vtk file and place it inside the same directory.
        n1=BigArray.shape[0]
        n2=BigArray.shape[1]
        ############################################
        iin = 5  # number of dimensions that will be graphed (all pairs will be produced)
        ############################################
        #### Open a .vtk to write ASCII data
        #### First, we need to get a name of the file.
        scrlen = len(dirs[ii])
        svd_file_prefix = files[0][scrlen+1:scrlen + 14]
        vtk_file_path = dirs[ii] + '\\'+ svd_file_prefix + '.vtk'
        save_file = open(vtk_file_path, 'w')
        ## Now the file is opened, we will write stuff into it.
        save_file.write('# vtk DataFile Version 2.0 \n')
        save_file.write('SVD Coded Traj ' + svd_file_prefix + '\n')
        save_file.write('ASCII\n')
        save_file.write('DATASET FIELD ' + svd_file_prefix + 'SVD_SolET_Columns {0}\n'.format(n2))
        for ii1 in range(n2-2):
            save_file.write(svd_file_prefix +'_COORD_{0} 1 {1} float\n'.format(ii1, n1))
            data_line = ''
            for ii2 in range(n1):
                data_line = data_line +'{0:8.5f} '.format(BigArray[ii2, ii1])
                if (ii2+1) % 10 == 0 and ii2 > 0:
                    save_file.write(data_line+'\n')
                    data_line = ''
            save_file.write(data_line+'\n')
            save_file.write('\n')
        ################ Now we save the entropy
        save_file.write(svd_file_prefix + '_Entropy 1 {0} float\n'.format(n1))
        data_line = ''
        for ii2 in range(n1):
            data_line = data_line + '{0:8.5f} '.format(BigArray[ii2, n2 - 2])
            if (ii2+1) % 10 == 0 and ii2 > 0:
                save_file.write(data_line+'\n')
                data_line = ''
        save_file.write(data_line+'\n')
        save_file.write('\n')
        ################# Now the time
        save_file.write(svd_file_prefix + '_Time 1 {0} float\n'.format(n1))
        data_line = ''
        for ii2 in range(n1):
            data_line = data_line + '{0:8.5f} '.format(BigArray[ii2, n2 - 2])
            if (ii2+1) % 10 == 0 and ii2 > 0:
                save_file.write(data_line+'\n')
                data_line = ''
        save_file.write(data_line+'\n')
        ################## All done. Close the file
        save_file.close()
####################################################################
# In this piece we will create a .vtk file where different 2D data will
# be stored ready for plotting. one part of data is 2D scalar plots. We will select a number of
# dimensions to plot and will be considering all possible pairs of variables. For each pair, we will
# prepare a curve plot of a curve with a value (most likely color-coded).
#####################################################################
create_vtk_files_2Dtrajectories=False
if create_vtk_files_2Dtrajectories:
    #first we will get the list of all directories where there is pickled dat
    # that can be used to create such vtk files.
    path = 'F:\\SimulationsLearningCollision\\take3_M41_SVDsols'
    dirs = my_readwrite.my_get_dirs_names(path)
    #Now we will call a subroutine that will scroll through all pickled
    # SVD encoded trajectories and will create a bunch of files and place them in path/VisItData/ to
    # be imported by visit.
    for ii in range(len(dirs)):
        files = my_vtk_tools.my_get_svdsols_file_names(dirs[ii])
        save_file = open(files[0], 'rb')
        BigArray = pickle.load(save_file)
        save_file.close()
        #BigArray contains SVD-encoded trajectories. We will now take a 2D slice/projection of a trajectory and record it into
        # a .vtk file and place it inside the same directory.
        n1=BigArray.shape[0]
        n2=BigArray.shape[1]
        ############################################
        iin = 2  # number of dimensions that will be graphed (all pairs will be produced)
        ############################################
        #### Open a .vtk to write ASCII data
        #### First, we need to get a name of the file.
        scrlen = len(dirs[ii])
        svd_file_prefix = files[0][scrlen+1:scrlen + 14]
        for ii1 in range(iin):
            for ii2 in range(ii1+1,iin):
                vtk_file_path = dirs[ii] + '\\'+ svd_file_prefix + '2D_I{0:02d}J{1:02d}.vtk'.format(ii1,ii2)
                save_file = open(vtk_file_path, 'w')
                ## Now the file is opened, we will write stuff into it.
                save_file.write('# vtk DataFile Version 2.0 \n')
                save_file.write('SVD Coded Traj ' + svd_file_prefix + 'I{0:00}J{0:00}\n'.format(ii1,ii2))
                save_file.write('ASCII\n')
                save_file.write('DATASET UNSTRUCTURED_GRID\n')
                save_file.write('POINTS {0} float\n'.format(n1))
                data_line = ''
                for ii3 in range(n1):
                    data_line = data_line + '{0:8.5f} '.format(BigArray[ii3, ii1])
                    data_line = data_line + '{0:8.5f} '.format(BigArray[ii3, ii2])
                    data_line = data_line + '  '
                    if (ii3+1) % 10 == 0 and ii3 > 0:
                       save_file.write(data_line+'\n')
                       data_line = ''
                save_file.write(data_line+'\n')
                save_file.write('\n')
                #############################
                ####### 2D points are written
                ####### Next we write the trajectory --- cell
                ####### there will be only one cell in this file
                #############################
                save_file.write('CELLS 1 {0}\n'.format(n1+1))
                data_line = '{0} '.format(n1)
                for ii3 in range(1,n1+1):
                    data_line = data_line + '{0} '.format(ii3-1)
                    if (ii3 + 1) % 20 == 0 and ii3 > 0:
                        save_file.write(data_line + '\n')
                        data_line = ''
                save_file.write(data_line + '\n')
                save_file.write('\n')
                #################################
                ### The cell is written. Now we write the type of the cell
                #################################
                save_file.write('CELL_TYPES 1\n'.format(n1 + 1))
                save_file.write('4\n')
                save_file.write('\n')
                #################################
                ### next we write the nodal values of the entropy
                ### that is scalar function that is
                ### defined on each point/node
                #################################
                save_file.write('POINT_DATA {0}\n'.format(n1))
                save_file.write('SCALARS Entropy float 1\n')
                save_file.write('LOOKUP_TABLE default\n')
                data_line = ''
                for ii3 in range(0,n1):
                    data_line = data_line + '{0:8.5f} '.format(BigArray[ii3, n2-2])
                    if (ii3 + 1) % 20 == 0 and ii3 > 0:
                        save_file.write(data_line + '\n')
                        data_line = ''
                save_file.write(data_line + '\n')
                save_file.write('\n')
                #### all done, the trajectory is written. Close the file
                save_file.close()
        ########################## Do it for all cobinations of indices i,j =1,..., inn
    ########################## Doe it for all folders.
####################################################################
# In this piece we will create a .vtk file where different 2D data will
# be stored ready for plotting. one part of data is 3D scalar plots. We will select a number of
# dimensions to plot and will be considering all possible pairs of variables. For each pair, we will
# prepare a curve plot of a curve with a value (most likely color-coded).
#####################################################################
create_vtk_files_3Dtrajectories=False
if create_vtk_files_3Dtrajectories:
    #first we will get the list of all directories where there is pickled dat
    # that can be used to create such vtk files.
    path = 'F:\\SimulationsLearningCollision\\take3_M41_SVDsols'
    dirs = my_readwrite.my_get_dirs_names(path)
    ###################################################
    ## Pull out a subset of trajectories
    ####################################################
    use_shorter_sample = True
    if use_shorter_sample:
        num_runs = 5  # number of runs to select
        runs_dir_names = []  # prepare an empty array that will keep names of the solutions to study
        #### GIVE IT A RUN
        for x in range(1000):
            ii = random.randint(0, num_runs_avail - 1)
        ######################

        for x in range(num_runs):
            ii = random.randint(0, num_runs_avail - 1)
            runs_dir_names.append(dirs[ii])
        #####################
        dirs = runs_dir_names
    ######################

    ####################################################
    #Now we will call a subroutine that will scroll through all pickled
    # SVD encoded trajectories and will create a bunch of files and place them in path/VisItData/ to
    # be imported by visit.
    for ii in range(len(dirs)):
        files = my_vtk_tools.my_get_svdsols_file_names(dirs[ii])
        save_file = open(files[0], 'rb')
        BigArray = pickle.load(save_file)
        save_file.close()
        #BigArray contains SVD-encoded trajectories. We will now take a 2D slice/projection of a trajectory and record it into
        # a .vtk file and place it inside the same directory.
        n1=BigArray.shape[0]
        n2=BigArray.shape[1]
        ############################################
        iin = 3  # number of dimensions that will be graphed (all pairs will be produced)
        ############################################
        #### Open a .vtk to write ASCII data
        #### First, we need to get a name of the file.
        scrlen = len(dirs[ii])
        svd_file_prefix = files[0][scrlen+1:scrlen + 14]
        for ii1 in range(iin):
            for ii2 in range(ii1+1,iin):
                for ii4 in range(ii2+1,iin):
                 vtk_file_path = dirs[ii] + '\\'+ svd_file_prefix + '3D_I{0:02d}J{1:02d}K{2:02d}.vtk'.format(ii1,ii2,ii4)
                 save_file = open(vtk_file_path, 'w')
                 ## Now the file is opened, we will write stuff into it.
                 save_file.write('# vtk DataFile Version 2.0 \n')
                 save_file.write('SVD Coded Traj ' + svd_file_prefix + 'I{0:02d}J{1:02d}K{2:02d}\n'.format(ii1,ii2,ii4))
                 save_file.write('ASCII\n')
                 save_file.write('DATASET UNSTRUCTURED_GRID\n')
                 save_file.write('POINTS {0} float\n'.format(n1))
                 data_line = ''
                 for ii3 in range(n1):
                    data_line = data_line + '{0:8.5f} '.format(BigArray[ii3, ii1])
                    data_line = data_line + '{0:8.5f} '.format(BigArray[ii3, ii2])
                    data_line = data_line + '{0:8.5f} '.format(BigArray[ii3, ii4])
                    data_line = data_line + '  '
                    if (ii3+1) % 10 == 0 and ii3 > 0:
                       save_file.write(data_line+'\n')
                       data_line = ''
                 save_file.write(data_line+'\n')
                 save_file.write('\n')
                 #############################
                 ####### 2D points are written
                 ####### Next we write the trajectory --- cell
                 ####### there will be only one cell in this file
                 #############################
                 save_file.write('CELLS 1 {0}\n'.format(n1+1))
                 data_line = '{0} '.format(n1)
                 for ii3 in range(1,n1+1):
                    data_line = data_line + '{0} '.format(ii3-1)
                    if (ii3 + 1) % 20 == 0 and ii3 > 0:
                        save_file.write(data_line + '\n')
                        data_line = ''
                 save_file.write(data_line + '\n')
                 save_file.write('\n')
                 #################################
                 ### The cell is written. Now we write the type of the cell
                 #################################
                 save_file.write('CELL_TYPES 1\n')
                 save_file.write('4\n')
                 save_file.write('\n')
                 #################################
                 ### next we write the nodal values of the entropy
                 ### that is scalar function that is
                 ### defined on each point/node
                 #################################
                 save_file.write('POINT_DATA {0}\n'.format(n1))
                 save_file.write('SCALARS Entropy float 1\n')
                 save_file.write('LOOKUP_TABLE default\n')
                 data_line = ''
                 for ii3 in range(0,n1):
                    data_line = data_line + '{0:8.5f} '.format(BigArray[ii3, n2-2])
                    if (ii3 + 1) % 20 == 0 and ii3 > 0:
                        save_file.write(data_line + '\n')
                        data_line = ''
                 save_file.write(data_line + '\n')
                 save_file.write('\n')
                 #### all done, the trajectory is written. Close the file
                 save_file.close()
        ########################## Do it for all cobinations of indices i,j,k =1,..., inn
    ########################## Doe it for all folders.
####################################################################
# In this piece we will create a .vtk file where different 3D data will
# be stored ready for plotting. All SVD coded trajectories will be bunched in one array and
# the points will be recorded as a point mesh .
#####################################################################
create_vtk_files_3Dsolpoints=True
if create_vtk_files_3Dsolpoints:
    #first we will get the list of all directories where there is pickled dat
    # that can be used to create such vtk files.
    path = 'E:\\SimulationsLearningCollision\\take3_M41_SVDsols\\good'
    dirs = my_readwrite.my_get_dirs_names(path)
    ###################################################
    ## Pull out a subset of trajectories
    ####################################################
    use_shorter_sample = True
    if use_shorter_sample:
        import random
        num_runs = 5  # number of runs to select
        runs_dir_names = []  # prepare an empty array that will keep names of the solutions to study
        #### GIVE IT A RUN
        for x in range(1000):
            ii = random.randint(0, len(dirs) - 1)
        ######################

        for x in range(num_runs):
            ii = random.randint(0, len(dirs) - 1)
            runs_dir_names.append(dirs[ii])
        #####################
        dirs = runs_dir_names
    ######################

    #Now we will call a subroutine that will scroll through all pickled
    # SVD encoded trajectories and bunch then up
    ############################################
    # first one:
    files = my_vtk_tools.my_get_svdsols_file_names(dirs[0])
    save_file = open(files[0], 'rb')
    BigArray = pickle.load(save_file)
    recs = [BigArray.shape[0]]
    BigCArry = pickle.load(save_file)
    CalcGrad = pickle.load(save_file)
    VarGrad = pickle.load(save_file)
    save_file.close()
    # and now the rest
    for ii in range(1,len(dirs)):
        files = my_vtk_tools.my_get_svdsols_file_names(dirs[ii])
        save_file = open(files[0], 'rb')
        BA_temp = pickle.load(save_file)
        BC_temp = pickle.load(save_file)
        CG_temp = pickle.load(save_file)
        VG_temp = pickle.load(save_file)
        save_file.close()
        BigArray = np.concatenate((BigArray, BA_temp), axis=0)
        BigCArry = np.concatenate((BigCArry, BC_temp), axis=0)
        CalcGrad = np.concatenate((CalcGrad, CG_temp), axis=0)
        VarGrad = np.concatenate((VarGrad, VG_temp), axis=0)
        recs.append(BA_temp.shape[0])
    # now BigArray contains SVD-encoded trajectories of all foind pickled SVD files.
    #We will now take a 3D slices/projections of a point dataset and record it into
    # a .vtk file and place it inside the same directory.

    n1=BigArray.shape[0]
    n2=BigArray.shape[1]
    ############################################
    iin = 5  # number of dimensions that will be graphed (all pairs will be produced)
    ############################################
    for ii1 in range(iin):
        for ii2 in range(ii1+1,iin):
             for ii4 in range(ii2+1,iin):
                 #### Open a .vtk to write ASCII data
                 #### First, we need to get a name of the file.
                 vtk_file_path = path + '\\SVDcodedSolsPoints_3D_I{0:02d}J{1:02d}K{2:02d}.vtp'.format(ii1,ii2,ii4)
                 save_file = open(vtk_file_path, 'w')
                 ## Now the file is opened, we will write stuff into it.
                 save_file.write('<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\"> \n')
                 save_file.write('  <PolyData>\n')                                       #ints ' + 'I{0:02d}J{1:02d}K{2:02d}\n'.format(ii1,ii2,ii4))
                 save_file.write('    <Piece NumberOfPoints=\"{0}\" NumberOfVerts=\"0\" '.format(n1) +
                                 'NumberOfLines=\"{0}\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n'.format(len(recs)))
                 save_file.write('      <Points>\n')
                 save_file.write('        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n')
                 data_line = ''
                 for ii3 in range(n1):
                    data_line = data_line + '{0:8.5f} '.format(BigArray[ii3, ii1])
                    data_line = data_line + '{0:8.5f} '.format(BigArray[ii3, ii2])
                    data_line = data_line + '{0:8.5f} '.format(BigArray[ii3, ii4])
                    data_line = data_line + '  '
                    if (ii3+1) % 10 == 0 and ii3 > 0:
                       save_file.write(data_line+'\n')
                       data_line = ''
                 save_file.write(data_line+'\n')
                 save_file.write('        </DataArray>\n')
                 save_file.write('      </Points>\n')
                 save_file.write('      <PointData Scalars=\"entropy\" Vectors=\"CalcGrad\">\n')
                 save_file.write('        <DataArray type=\"Float32\" Name=\"entropy\" format=\"ascii\">\n')
                 data_line = ''
                 for ii3 in range(0, n1):
                     data_line = data_line + '{0:8.5f} '.format(BigArray[ii3, n2 - 2])
                     if (ii3 + 1) % 20 == 0 and ii3 > 0:
                         save_file.write(data_line + '\n')
                         data_line = ''
                 save_file.write(data_line + '\n')
                 save_file.write('        </DataArray>\n')
                 save_file.write('        <DataArray type=\"Int32\" Name=\"NumTrajPoints\" format=\"ascii\">\n')
                 data_line = ''
                 nnn = 0
                 for ii3 in range(len(recs)):
                     for ii5 in range(recs[ii3]):
                      nnn = nnn + 1
                      data_line = data_line + '{0} '.format(nnn)
                      if (nnn) % 20 == 0:
                         save_file.write(data_line + '\n')
                         data_line = ''
                 save_file.write(data_line + '\n')
                 save_file.write('        </DataArray>\n')
                 #############################
                 ## we now add two vector fields.
                 ## The first vecotr filed is a calculus gradient of entropy. The second is a projection of the
                 ## caclulus gradient onto a mass-conservative tangent space.
                 #############################
                 save_file.write('        <DataArray type=\"Float32\" Name=\"CalcGrad\" NumberOfComponents=\"3\" format=\"ascii\">\n')
                 data_line = ''
                 for ii3 in range(0, n1):
                     data_line = data_line + '{0:8.5f} '.format(CalcGrad[ii3, ii1])
                     data_line = data_line + '{0:8.5f} '.format(CalcGrad[ii3, ii2])
                     data_line = data_line + '{0:8.5f} '.format(CalcGrad[ii3, ii4])
                     if (ii3 + 1) % 10 == 0 and ii3 > 0:
                         save_file.write(data_line + '\n')
                         data_line = ''
                 save_file.write(data_line + '\n')
                 save_file.write('        </DataArray>\n')
                 save_file.write('        <DataArray type=\"Float32\" Name=\"VarGrad\"  NumberOfComponents=\"3\" format=\"ascii\">\n')
                 data_line = ''
                 for ii3 in range(0, n1):
                     data_line = data_line + '{0:8.5f} '.format(VarGrad[ii3, ii1])
                     data_line = data_line + '{0:8.5f} '.format(VarGrad[ii3, ii2])
                     data_line = data_line + '{0:8.5f} '.format(VarGrad[ii3, ii4])
                     if (ii3 + 1) % 10 == 0 and ii3 > 0:
                         save_file.write(data_line + '\n')
                         data_line = ''
                 save_file.write(data_line + '\n')
                 save_file.write('        </DataArray>\n')
                 save_file.write('        <DataArray type=\"Float32\" Name=\"CollOper\"  NumberOfComponents=\"3\" format=\"ascii\">\n')
                 data_line = ''
                 for ii3 in range(0, n1):
                     data_line = data_line + '{0:8.5f} '.format(BigCArry[ii3, ii1])
                     data_line = data_line + '{0:8.5f} '.format(BigCArry[ii3, ii2])
                     data_line = data_line + '{0:8.5f} '.format(BigCArry[ii3, ii4])
                     if (ii3 + 1) % 10 == 0 and ii3 > 0:
                         save_file.write(data_line + '\n')
                         data_line = ''
                 save_file.write(data_line + '\n')
                 save_file.write('        </DataArray>\n')
                 save_file.write('      </PointData>\n')
                 ####### Next we will record trajectories -- lines
                 #############################
                 save_file.write('      <Lines>\n')
                 save_file.write('        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n')
                 data_line = ''
                 for ii3 in range(0, n1):
                     data_line = data_line + '{0} '.format(ii3)
                     if (ii3 + 1) % 20 == 0 and ii3 > 0:
                         save_file.write(data_line + '\n')
                         data_line = ''
                 save_file.write(data_line + '\n')
                 save_file.write('        </DataArray>\n')
                 save_file.write('        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n')
                 nnn = 0
                 data_line=''
                 for ii3 in range(len(recs)):
                     nnn = nnn + recs[ii3]
                     data_line = data_line + '{0} '.format(nnn)
                     if (ii3 + 1) % 20 == 0 and ii3 > 0:
                         save_file.write(data_line + '\n')
                         data_line = ''
                 save_file.write(data_line + '\n')
                 save_file.write('        </DataArray>\n')
                 save_file.write('      </Lines>\n')
                 save_file.write('      <CellData Scalars=\"TrajNum\">\n')
                 save_file.write('        <DataArray type=\"Int32\" Name=\"TrajNum\" format=\"ascii\">\n')
                 data_line = ''
                 for ii3 in range(len(recs)):
                     data_line = data_line + '{0} '.format(ii3)
                     if (ii3 + 1) % 20 == 0 and ii3 > 0:
                         save_file.write(data_line + '\n')
                         data_line = ''
                 save_file.write(data_line + '\n')
                 save_file.write('        </DataArray>\n')
                 save_file.write('      </CellData>\n')
                 ####### trajectories ar written
                 #############################
                 save_file.write('    </Piece>\n')
                 save_file.write('  </PolyData>\n')  # ints ' + 'I{0:02d}J{1:02d}K{2:02d}\n'.format(ii1,ii2,ii4))
                 save_file.write('</VTKFile>\n')
                 #### all done, the points are written. Close the file
                 save_file.close()
        ########################## Do it for all cobinations of indices i,j,k =1,..., inn
    ########################## Doe it for all folders.


###############################
print(range(5))
###########
