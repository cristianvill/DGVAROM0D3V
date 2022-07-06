##############################
###  12/30/2020 A. Alekseenko
###  This is a collection of code pieces to
###  create a ROM basis and to save it for Fortran to read.
###
##############################

##############
import numpy as np
import my_readwrite   ## import the entire module

################################
##
## First we open the pickled file with 100 Singular vectors from the "100 cleaned solutions" set.
##
################################
MM=41
Mtrim=0
##
import pickle
path = "F:\\temp\\cleaned_SVs_100"
save_file = open(path + '\\sol_SVD_100.dat', 'rb')
Vh = pickle.load(save_file)
s = pickle.load(save_file)
svect = pickle.load(save_file)
save_file.close()
for ii in range(100):
    print(s[ii])
##################################
##
## now we have the 100 singular vectors from 100 cleaned solutions data set loaded in numpy array svect
##
##################################
## Our next step is to replace the first singular vector with the steady state maxwellian
##################################
# to add a maxwellian, we need the mesh from our dicscrete solution.
# specifically, we need arrays of velocity points and gauss weights. these
# points can be obtained from the DG-Boltzmann solution saves:
filename='F:\\SimulationsLearningCollision\\take3_M41\\good\\080\\CollTrnDta180_1su1sv1sw41MuUU41MvVU41MwWU_nodes.dat'
nodes_u, nodes_v, nodes_w, nodes_gwts = my_readwrite.read_nodes(filename)
##################################
# now we evaluate the Maxwellian
import distributions
Mxwl_vdf = distributions.dimless_mxwln(nodes_u, nodes_v, nodes_w,1,0.0,0.0,0.0,0.2)
####################################
# let us compare the first singular vector and the maxwellian
# we will compute the orthogonal compliment of maxwellian to the zeros vector and check its norm
# and than take the l2 norm of the difference.
####################################
dot_prod = np.sum(svect.T[0,:]*Mxwl_vdf)
u0_norm2 = np.sum(svect.T[0,:]**2)
tempr = np.sum(((nodes_u)**2 + (nodes_v)**2 + (nodes_w)**2)*Mxwl_vdf*nodes_gwts)*2/3
u0_perp = Mxwl_vdf - dot_prod/u0_norm2*svect.T[0,:]
u0_perp_norm = np.sqrt(np.sum(u0_perp**2))/np.sqrt(np.sum(Mxwl_vdf**2))
####################################
## if this norm is small, we can try to replace the first singluar vector. If the norm is large, say more than .2
## then we will keep both vectors in the bunch.
##
## Well...  it turns out that the projection has about 13% of the maxwellian. It is not negligible.
## Still two possible approaches: replacve the first sv with maxwellian and add the Maxwellian upfront
## I guess we will have to try both cases.
####################################
## both arrays need to be be run through QR so that the columns are orthonormal
Z001 = np.concatenate((Mxwl_vdf.T,svect[:,1:]), axis=1)
arryRplFirst2Mxwell, R_m = np.linalg.qr(Z001)
##################### Save Singular Vectors into a binary array to be read in Fortran
# UNcomment if need to save on disk
#my_readwrite.save_SVs_fortr("SpHom100Trim0SVs_RplFirstMxwellReadFort.dat",arryRplFirst2Mxwell.T,MM,Mtrim)
### Save in pickled form
save_file = open('SpHom100Trim0SVs_RplFirstMxwellPickle.dat', 'wb')
pickle.dump(arryRplFirst2Mxwell, save_file)
save_file.close()
#####################################
Z001 = np.concatenate((Mxwl_vdf.T,svect[:,:]), axis=1)
arryMxwelladdSVs, R_m = np.linalg.qr(Z001)
##################### Save Singular Vectors into a binary array to be read in Fortran
# UNcomment if need to save on disk
#my_readwrite.save_SVs_fortr("SpHom100Trim0SVs_addMxwellReadFort.dat",arryMxwelladdSVs.T,MM,Mtrim)
### Save in pickled form
save_file = open('SpHom100Trim0SVs_addMxwellPickle.dat', 'wb')
pickle.dump(arryMxwelladdSVs, save_file)
save_file.close()
###################### We save a copy that does not include a maxwellian ############
# Uncomment if need to save on disk
my_readwrite.save_SVs_fortr("SpHom100Trim0SVs_NoMxwellReadFort.dat",arryMxwelladdSVs[:,1:].T,MM,Mtrim)
### Save in pickled form
save_file = open('SpHom100Trim0SVs_NoMxwellPickle.dat', 'wb')
pickle.dump(arryMxwelladdSVs[:,1:], save_file)
save_file.close()
#####################################


#####################################
## Next we will create a momnet free basis it will consist of the Maxwellian (l2-normalized) and
# the Singular vectors, that were perturbed a bit so that they satisfy the zero conservative moments
# the perturbations are inserted at points where the value of the singular vectors is
# above a treshhold. After that, the conservative vectors need to be re-orthogonalized.
# While we perturb the vectors we will generally lose orthogonality to the first vecotor of the bunch,
# the maxwellian, we can make re-orthogonalized. It would be great if we could keep it orthogonal to the
# Maxwellian.
#####################################
# we start from creating the moment-kernels
Mat_Moments = np.zeros((arryMxwelladdSVs.shape[0],5))
Mat_Moments[:,0] = nodes_gwts[0,:]
Mat_Moments[:,1] = nodes_gwts[0,:]*nodes_u[0,:]
Mat_Moments[:,2] = nodes_gwts[0,:]*nodes_v[0,:]
Mat_Moments[:,3] = nodes_gwts[0,:]*nodes_w[0,:]
scrp_array = (nodes_u**2+nodes_v**2+nodes_w**2)
Mat_Moments[:,4] = nodes_gwts[0,:]*scrp_array[0,:]
### Now matrix Mat_Moments[:,2] contains the vectors that are kernels of the first five moments -- which are the coserved moments.
## Next we will orthogonalize these
Q_m, R_m = np.linalg.qr(Mat_Moments)
## Next we will develop a basis in which subtract projections of vectors in arryMxwelladdSVs[1:,:]
## matrix of violation of constraints
Bmat = (-1)*np.matmul(Mat_Moments.T,arryMxwelladdSVs[:,1:])
### Next we will compute a perturbation to the singular vectors that will make them to satisfy the constraint -- it
### will offset the quantities accumulated in Bmat. The pertubations will be performed where the singular vector is
### larger than a threshold value. This way, we will not perturb local supportness of the singular vectors.
### The constructed perturbations will make the vectors conserved moments -free (will have zero mass, momentum and energy)
### it will also be orthrogonal to the first vector.
### Here is the basic math assuming that SVs are perturbed in all components. Generalization to masked perturbation is easy
### M - matrix of moment kernels
### S - matrix of singular vectors to be processed.
### b = -M^T S -- the amount of violation of conserved moments
### f - first vector --- Maxwellian
### Need to find D such that
### M^T (S+D)=0
### f^T D = 0
##############
### Solution: Perturbation of minimal l2 error:
##############
### Let W be the marix made out of, W=[M f]
### Consider SVD of W = UGV^T
### Seek least l2-norm solution to W^T d_i = [b_i,0]^T column-wise
### where d_{i} is a columnt of D and b_{i} is the corresponding column of b
### Assume d_i=Uz_{i} (linear combination of columns of U)
### then z_{i}=G^{-1} V^T [b_i,0]
######### Now we will implement this, except we also add a mask based on the absolute value of the SV.
Wmat = np.concatenate((Mat_Moments,arryMxwelladdSVs[:,0:1]), axis=1)
#############
Dmat = np.zeros((arryMxwelladdSVs.shape[0],arryMxwelladdSVs.shape[1]-1))  # create a space.
for i in range(Dmat.shape[1]):
    # we will process one vector at a time.
    e_i = np.concatenate((Bmat[:,i],np.zeros((1))), axis=0)
    ##### Now we will implement the mask
    trshld = 0.01   ### threshhold for the mask -- zero mask is applied to values less trshld by abs. val.
    max_value = np.amax(np.absolute(arryMxwelladdSVs[:,1+i]))
    index_map = np.zeros((arryMxwelladdSVs.shape[0]))
    ###
    Wmat_big = np.zeros((0,6))
    for j in range(arryMxwelladdSVs.shape[0]):
        if abs(arryMxwelladdSVs[j, 1+i]) >= trshld * max_value:
            Wmat_big = np.concatenate((Wmat_big , Wmat[j:j+1,:]), axis=0)
            index_map[j] = Wmat_big.shape[0]
        else:
            index_map[j] = 0
    # Now, Wmat_big is a copy of Wmat with some entries masked.
    U, g, Vh = np.linalg.svd(Wmat_big, full_matrices=False)
    z_i = np.matmul(Vh, e_i) / g
    d_i_big = np.matmul(U, z_i)
    # now we produce the full masked perturbation:
    for j in range(arryMxwelladdSVs.shape[0]):
        if index_map[j]>0:
          Dmat[j, i] = d_i_big[int(index_map[j])-1]
    # all done witht he perturbation
## Diagnostic
SVZero, R_m = np.linalg.qr(arryMxwelladdSVs[:,1:]+Dmat)
Z002=np.matmul(Wmat.T,SVZero)
###

##################### Save Singular Vectors into a binary array to be read in Fortran
# UNcomment if need to save on disk
my_readwrite.save_SVs_fortr("SpHom100Trim0SVs_SVZeroBasisReadFort.dat",SVZero.T,MM,Mtrim)
save_file = open('SpHom100Trim0SVs_SVZeroBasisPickle.dat', 'wb')
pickle.dump(SVZero, save_file)
save_file.close()
#####################################
print(range(5))