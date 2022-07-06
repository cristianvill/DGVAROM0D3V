#####################################
###  09/10/2020 A. Alekseenko
###  This code will investigate the asymptotic stability of the ROM model and will
###  Explore conservtive projection.
###
#####################################

### TRIMMING PARAMETERS
MM = 41
Mtrim = 0
###
k=12
###

import numpy as np
import my_readwrite


# first we load the three index array of the ROM kernel. The third index is the non-symmetric index -- the output index
CKern_all = my_readwrite.readromkrnl("E:\\temp\\cleaned100_MkKrnl\\M41Trim0\\SpHomT2D1M41Trm0SVKrnl_k53.dat")
######################
##  Our first exercise  is to compute the linearization term matrix B^{j.}_{.k}in the error equation
##
##  \partial_{t} e_{k} = B^{j}_{k}f e_{j} + A_{..k}^{ij.} e_{i}e_{j}
##   Formula for the B_{jk} is A^{i..}_{.jk} \hat{u}_{i} , where \hat{u}_{i} is the projection of the steady state maxwellian
## Thus we will need to load projection matrix, the mech, evaluate maxwellian on the mesh, project the maxwellian and
## compute the array convolution.
######################

######################
# We need some prep work. We will need to get nodes and weights of the quadratures associated with the solution loaded.
# these will be used to evaluate the solution entropy.
# these points can be obtained from the DG-Boltzmann solution saves of mesh parameters:
filename = 'E:\\SimulationsLearningCollision\\take3_M41\\good\\080\\CollTrnDta180_1su1sv1sw41MuUU41MvVU41MwWU_nodes.dat'
nodes_u, nodes_v, nodes_w, nodes_gwts = my_readwrite.read_nodes(filename)
#########################################################################
#########################################################################
Mat_Moments = np.zeros((nodes_u.shape[1], 5))
Mat_Moments[:, 0] = nodes_gwts[0,:]
Mat_Moments[:, 1] = nodes_gwts[0,:] * nodes_u[0, :]
Mat_Moments[:, 2] = nodes_gwts[0,:] * nodes_v[0,:]
Mat_Moments[:, 3] = nodes_gwts[0,:] * nodes_w[0,:]
scrp_array = (nodes_u ** 2 + nodes_v ** 2 + nodes_w ** 2)
Mat_Moments[:, 4] = nodes_gwts[0,:] * scrp_array[0,:]
########################################################################
##########################################################################
import my_distributions
moment_vector_0 = np.array([[1.0,0.0,0.0,0.0,0.3]]) # ATTENTION: Last entry is energy, not temperature
sol_maxwell = my_distributions.maxwellian(moment_vector_0,nodes_u,nodes_v,nodes_w)  # evaluate truncated maxwellian with the same macroparameters

#########################################################################
### Load the projection vectors
#########################################################################
path = "E:\\temp\\cleaned_SVs_100\\sol_SVD_100.dat"
### There is a file with all singular vectors and with just the first 100
### Please pay attention that the trimming is correct
import pickle
save_file = open(path, 'rb')
Vh = pickle.load(save_file)
s = pickle.load(save_file)
svect = pickle.load(save_file)
save_file.close()
for ii in range(100):
   print(s[ii])
#########################################################################
## we will need only a subset of the projection vectors
#########################################################################
# We create an array that collects all the eigenvalues
all_eigs =np.zeros((54,51))
k_all = CKern_all.shape[0]

############################################################
#
make_csv_eigenvalues=False
if make_csv_eigenvalues:
 for k in range(3,54):
  if k<k_all:
      CKern=CKern_all[0:k,0:k,0:k]
  else:
      CKern=CKern_all
      k=k_all
  Proj = svect[:,0:k]
  ###
  uhat = np.matmul(sol_maxwell,Proj)
  Bmat = np.tensordot(CKern,uhat,axes=([0],[0]))
  ###
  ## Now we will study Bmat
  evals, evects = np.linalg.eig(Bmat)
  all_eigs[0,k-3] = k
  all_eigs[1:k+1,k-3]=evals[:]
 ### All done, Now export Bmat to .csv file
 my_readwrite.my_writeMatrix_cvs(all_eigs)

####################################################
#
study_steasy_state=True
if study_steasy_state:
    k = 53
    if k < k_all:
        CKern = CKern_all[0:k, 0:k, 0:k]
    else:
        CKern = CKern_all
        k = k_all
    Proj = svect[:, 0:k]
    ###
    uhat = np.matmul(sol_maxwell, Proj)
    Bmat = np.tensordot(CKern, uhat, axes=([0], [0]))
    ###
    # check conservaive properties of projection-recovery
    moms_recover = np.matmul(np.matmul(Proj, uhat),Mat_Moments)
    ###
    ## Now we will study Bmat
    evals, evects = np.linalg.eig(Bmat)
    evects_real = np.real(evects)
    all_eigs[0, k - 3] = k
    all_eigs[1:k + 1, k - 3] = evals[:]
    ############### NOW WE WILL TAKE A LOOK at is the steady state for Boltzmann is really a steady state for the system.
    st_state = np.tensordot(Bmat, uhat, axes=([0], [0]))
    project_evects = np.matmul(uhat, np.real(evects))
    #########################
    ## next we will try to look at the steady state that the system is converging to ....
    path = 'E:\\temp\\cleamed_100_ROMruns\\results\\413_record\\CollTrnDta413_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.2715000000_SltnColl.dat'
    solarry, collarry, solsize = my_readwrite.my_read_sol_coll_trim(path, MM, Mtrim)
    uhat1 = np.matmul(solarry, Proj)
    stst_project_evects = np.matmul(uhat1, np.real(evects))
    duhat=np.absolute(uhat1-uhat)/np.absolute(uhat)
    proj_err=np.absolute(stst_project_evects-project_evects)/np.absolute(project_evects)


    for i in range(k):
        print(i,evals[i])
range(5)

