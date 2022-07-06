#####################################
###  09/10/2020 A. Alekseenko
###  This code will investigate the linear stability of the SV Zero ROM model and will explore
###  correcting the linear operator
###
#####################################


import numpy as np
import my_readwrite

####################################################
## First we need to import the linear kernel #######
BKern_all, nn = my_readwrite.my_read_linkrnl("F:\\temp\\cleaned_SVs_100_PlusMaxwell\\SVZeroCln100M41Tr0KrnlP5_Bkrnl.dat")
##
#########################################################################
## we will need only a subset of the projection vectors
#########################################################################
# We create an array that collects all the eigenvalues
all_eigs =np.zeros((nn+1,nn))

############################################################
#
from scipy.sparse.linalg import eigs
make_csv_eigenvalues=True
if make_csv_eigenvalues:
 for k in range(0,nn):
  if k < nn-1:
      BKern=BKern_all[0:k+1,0:k+1]
  else:
      BKern=BKern_all
      k = nn-1
  ###
  #evals, evects = np.linalg.eig(BKern)

  evals, evects = eigs(BKern,k+1)
  all_eigs[0,k] = k+1
  all_eigs[1:k+2,k]=np.array(evals,ndmin=1).reshape([k+1])
 ### All done, Now export Bmat to .csv file
 my_readwrite.my_writeMatrix_cvs(all_eigs)
