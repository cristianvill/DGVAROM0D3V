#####################################
## This program collects data on spatially homogeneous relaxations and
## saves it in a pickled file. It also makes SDV on data on solutions
## and saves it as a pickled file and also as a file readble by fortran
#####################################

### TRIMMING PARAMETERS
MM = 41
Mtrim = 0
### Directory that keeps solutions (in subfolders)
path='/pylon5/ms5phtp/alekseen/cleaned_data'
#### Directory to save the pickled files ########
pathpickled='/pylon5/ms5phtp/alekseen/cleaned_data_work'


##############
import numpy as np
import my_readwrite

### get names of all files in this folder, and its subfolders as one giant list
names = my_readwrite.my_get_soltn_file_names_time (path, 0.31)
zz = len(names)
print("Found {0:03d} files".format(zz))

#######################################################################
# let us load training data, which will be instances of solutions.
# we read the files, get solutions and values of collision operator out of them and
# copy solution into a list of arrays. Each array contains one instance of the solution
### TRIMMED read:
sol_data_train, col_data_train, solsize = my_readwrite.my_read_sol_coll_trim(names[0],MM, Mtrim)  # this is the first file
# let us add more data points/read a few more files
for i in range(zz):
    solarry, colarry, solsize1 = my_readwrite.my_read_sol_coll_trim(names[i], MM, Mtrim)
    assert solsize1==solsize  # a check that all arrays have the same size -- should be true, but still..
    sol_data_train = np.concatenate((sol_data_train, solarry), axis = 0)
    col_data_train = np.concatenate((col_data_train, colarry), axis = 0)
# done creating a list of training arrays
#########################################
## Now we will save the arrays in pickled files
import pickle
save_file = open(pathpickled+'sol_data_train.dat', 'wb')
pickle.dump(sol_data_train, save_file)
save_file.close()
save_file = open(pathpickled+'col_data_train.dat', 'wb')
pickle.dump(col_data_train, save_file)
save_file.close()
#####################################################################################
## NEXT We will try to do some linear algebra on these arrays #######################
from scipy import linalg
U, s, Vh =  linalg.svd(sol_data_train.T,full_matrices=False)
print("the first 100 singular vectors: \n", s[0:100])
save_file = open(pathpickled + 'sol_SVD_100.dat', 'wb')
pickle.dump(Vh[0:100, :], save_file)
pickle.dump(s[0:100], save_file)
pickle.dump(U[:, 0:100], save_file)
save_file.close()
#######################################################################################
#### Check accuracy of the truncated SDV:
Vhtr = Vh[0:100,:]
Utr = U[:, 0:100]
str = s[0:100]
err = np.amax(np.absolute( sol_data_train  -  np.matmul(Utr,np.matmul(np.diag(str), Vhtr)) ))
#######################################################################################
print("saved the 100 singular vectors: \n", s[0:100])

exit()
