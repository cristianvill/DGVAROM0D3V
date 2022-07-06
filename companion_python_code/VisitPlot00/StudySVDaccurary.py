##################################
###  08/28/2020 A. Alekseenko
###  This code will read a soluton and project and
### read an array of singular vectors of the solution data matrix.
### Then the code will project a recover the solution and evaluate a bunch of
### norms of the difference.
###
###
####################################

### TRIMMING PARAMETERS
MM = 41
Mtrim = 0

#########################
### Read pickled SVD

import numpy as np
import pickle
import my_readwrite
import my_vtk_tools
path = "E:\\temp\\cleaned_SVs_100"
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
    save_file = open(path + '\\sol_SVD_100.dat', 'rb')
    Vh = pickle.load(save_file)
    s = pickle.load(save_file)
    svect = pickle.load(save_file)
    save_file.close()
    for ii in range(100):
        print(s[ii])

##########################
DI100 = np.matmul(svect.T,svect) - np.eye(100)
maxrec=np.amax(DI100)
minrec=np.amin(DI100)
###########################


# Next we will read the solution
#sname="E:\\SimulationsLearningCollision\\take3_M41\\good\\432\\CollTrnDta432_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0000000000_SltnColl.dat"
sname="E:\\temp\\sphomruns\\good\\432\\CollTrnDta432_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0000000000_SltnColl.dat"
solarry, collarry, solsize = my_readwrite.my_read_sol_coll_trim(sname, MM, Mtrim)
#####################
#results will be written into a text file
out_file = "accuracry432.txt"
save_file = open(out_file, 'w')
#####################
# Set upthe list of value for k
k_list =[16,27,35,42,53]
for i in k_list:
    U = svect[:,0:i-1]
    solprojrec = np.matmul(np.matmul(solarry,U),U.T)
    ## Now we will evaluate the errors
    Linf_err_abs = np.amax(np.absolute(solarry-solprojrec))
    Linf_err_rel = np.amax(np.absolute(solarry - solprojrec))/np.amax(np.absolute(solarry))

    L1_err_abs = (np.sum(np.absolute(solarry - solprojrec)))
    L1_err_rel = (np.sum(np.absolute(solarry - solprojrec)))/np.sum(np.absolute(solarry))

    L2_err_abs = np.sqrt(np.sum((solarry - solprojrec)**2))
    L2_err_rel = np.sqrt(np.sum((solarry - solprojrec)**2))/np.sqrt(np.sum((solarry)**2))
    save_file.write("\n")
    save_file.write("k={0},\n".format(i))
    save_file.write("Linf_abs={0:14.12f}, Linf_rel={1:14.12f}, \n ".format(Linf_err_abs,Linf_err_rel) )
    save_file.write("L1_abs={0:14.12f}, L1_rel={1:14.12f}, \n ".format(L1_err_abs,L1_err_rel) )
    save_file.write("L2_abs={0:14.12f}, L2_rel={1:14.12f}, \n ".format(L2_err_abs,L2_err_rel) )
    save_file.write("\n")
save_file.close()