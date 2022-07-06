##################################
###  07/28/2020 A. Alekseenko
###  This code will run over the folders of computed solutions to
### spatially homogeneous relaxation and will select solutions
### that satisfy conservation law within given accuracy.
### These solutions will be then copied to a new folder, that contains such solutions
### Another action that will be performed is that some timesaves will not be copied
### The time saves will be selected by the value of the deviation from the steady state -- L1 Error, referredas the Error.
### A scale in Error will be designed.
### Error will be evaluated on the first save. The First save is copied into the clean adata folder regardless.
### The new trashhold for the error is estblished by subtracting a aset step value, err_step from the value of error for the
### first save. The saved solution that fall below the treshhold is copied to the folder and the new treshhold is established.
### this process will continue until we reach the last save that is copied regardless it falling below the treshhold.
###
###
##################################

### TRIMMING PARAMETERS
MM = 41
Mtrim = 0
###
err_step=0.005
###

import numpy as np
import my_readwrite

#########################
## First we will prepare list of directories that contain single runs of solutions
######################################
path = 'F:\\temp\\sphomruns\\bad'
allsolsubfolder = 'bad'
cleanedsolsubfolder = 'bad_cl'
dirs = my_readwrite.my_get_dirs_names(path)
num_runs = len(dirs)  # total available runs
######################
# We need some prep work. We will need to get nodes and weights of the quadratures associated with the solution loaded.
# these will be used to evaluate the solution entropy.
# these points can be obtained from the DG-Boltzmann solution saves of mesh parameters:
filename = 'F:\\SimulationsLearningCollision\\take3_M41\\good\\080\\CollTrnDta180_1su1sv1sw41MuUU41MvVU41MwWU_nodes.dat'
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
########################################################################
import my_distributions
import os
############################################################################
# now we start to sort through the folders
for ii in range(num_runs):
    ####################################
    ### From all available testing solutions get an entire spatially homogeneous relaxation run
    ### now we get file names into an array
    print("checking directory", dirs[ii])
    names, timearry = my_readwrite.my_get_soltn_file_names_after_time_wtime(dirs[ii], -0.1)
    tt1 = names.sort()
    tt2 = timearry.sort()
    #### Got names of all solution saved files from folder dirs[ii]
    #### Now will check conservation by looking at the first and last saved solution:
    solarry_0, collarry_0, solsize = my_readwrite.my_read_sol_coll_trim(names[0], MM, Mtrim)
    moment_vector_0 = np.matmul(solarry_0, Mat_Moments)
    solarry_N, collarry_N, solsize = my_readwrite.my_read_sol_coll_trim(names[-1], MM, Mtrim)
    moment_vector_N = np.matmul(solarry_N, Mat_Moments)
    #### compare the moments
    err0mom = np.abs(moment_vector_N[0, 0] - moment_vector_0[0, 0])
    err1mom = np.abs(moment_vector_N[0, 1] - moment_vector_0[0, 1])
    err2mom = np.abs(moment_vector_N[0, 2] - moment_vector_0[0, 2])
    err3mom = np.abs(moment_vector_N[0, 3] - moment_vector_0[0, 3])
    err4mom = np.abs(moment_vector_N[0, 4] - moment_vector_0[0, 4])
    ###########################################################################################
    ###### ATTENTION: HARD CODED CONSTANTS, this may change on a different class of solutions
    ############################################################################################
    ####
    if err0mom < 0.001 and err1mom < 0.001 and err2mom < 0.001 and err3mom < 0.001 and err4mom < 0.005:
        newpath = dirs[ii].replace(allsolsubfolder, cleanedsolsubfolder, 1)   ##########  Will need to be changed for different folders
        os.mkdir(newpath)
        ############################
        print('Found good solution: '+ dirs[ii]+'\n')
        #### solution dirs[0] is already in solarray_0, so we will do nto have to read it.
        ####  ready to move the files
        sol_maxwell = my_distributions.maxwellian(moment_vector_0,nodes_u_trim,nodes_v_trim,nodes_w_trim)  # evaluate truncated maxwellian with the same macroparameters
        sol_maxwell =  np.array(sol_maxwell).reshape((1,(MM-2*Mtrim)**3))
        L1_Err_sol = np.sum(np.absolute(solarry_0-sol_maxwell)*nodes_gwts_trim)/moment_vector_0[0,0]
        L1_Err_trshhold = L1_Err_sol - err_step # next treshhold for L1_error (L1 norm of deviation from maxwellian)
        #### copy the file into a clean folder
        copyfilename = names[0].replace(allsolsubfolder, cleanedsolsubfolder, 1)
        os.system('copy '+names[0] + ' ' + copyfilename)
        #######################################
        print("Copying first save:"+names[0]+'\n')
        print("L1_err_sol=",L1_Err_sol)
        #######################################
        for iii in range (1,len(names)-1):
            solarry_0, collarry_0, solsize = my_readwrite.my_read_sol_coll_trim(names[iii], MM, Mtrim)
            moment_vector_0 = np.matmul(solarry_0, Mat_Moments)
            sol_maxwell = my_distributions.maxwellian(moment_vector_0, nodes_u_trim, nodes_v_trim,nodes_w_trim)  # evaluate truncated maxwellian with the same macroparameters
            sol_maxwell = np.array(sol_maxwell).reshape((1,(MM - 2 * Mtrim) ** 3))
            L1_Err_sol = np.sum(np.absolute(solarry_0-sol_maxwell)*nodes_gwts_trim)/moment_vector_0[0,0]  # (L1 norm of deviation from maxwellian)
            if L1_Err_sol<L1_Err_trshhold:
                copyfilename = names[iii].replace(allsolsubfolder, cleanedsolsubfolder, 1)
                os.system('copy ' + names[iii] + ' ' + copyfilename)
                L1_Err_trshhold = L1_Err_sol - err_step
                ##########################################
                print("Copying first save:" + names[iii] + '\n')
                print("L1_err_sol=", L1_Err_sol)
                ##########################################
        #### the last save solution is copied anyway
        copyfilename = names[-1].replace(allsolsubfolder, cleanedsolsubfolder, 1)
        os.system('copy ' + names[-1] + ' ' + copyfilename)
        #### all done here.
        print("Copying last save:")
    #### all done here
##########
print(range(5))