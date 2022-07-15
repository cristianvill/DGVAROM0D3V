!
!  DGV_dta_drvn_bol_macro.f90
!
! This module contains routines that are involved in the evaluation of the spatial operator using data driven ROM models 
! or used to develope pre-computed arrays that are used in ROM models
!
! he main purpoose of this module is to avoid circular dependencies with high level DGV modules. 
!
!
! This module contains subroutines that are dependend on other high level DGV modules. There is an additional module, 
! DGV_data_driven_boltz.f90 that contains subrounites that are not dependent on other high level DGV modules.
!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module DGV_dta_drvn_bol_macro
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
   implicit none


                                    

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Mk_Coll_Ker_SVD_Basis_MACRO_Serial
!
! This is an adaptation of the original MPI-based subroutine that used O(n^8) methods to 
! single node based FFB based model. The FFB uses Fourier transform and it is not clear if
! it will be easy to parallelize the code using OpenMP. The loop will probably have to be re-structured. We will 
! start with just making is serial (although we do not know how the LAPACK's FFT will behave -- it may will be OpenMP
!
!
!
! This subroutine is a cut and paste from the code that implements 
! developmen of collision kernel for global data driven basis for a class of solutions
! The data driven basis is obtained by performing an SVD on collection of solutions to a class 
! of problems. 
!
! The subroutin will read a file with singular vectors. These singular vectors will constutie 
! the optimal basis to represent this particular class of solutions. 
!
! The subroutine will compute a three dimensional array Asvd(i,j,k) where i,j,k run over 
! indices of the SVD singular vectors. 
! 
!
! Asvd(i,j,k) = (mech related const)*\int_{R^3} \int_{R^3}\int_{R^3} v_{i} (u-\xi) v_{j} (w-\xi) A(u, w) \v_{k}(\xi) du dw d\xi 
!
! To increase scalability of the process, the Asvd will be stored in the memory and on the hard drive as a single index array.
! 
! we note that in the case of symmetric formulation of the Glerkin projectiuon of the collision operator, Asvd(i,j,k) = Asvd(j,i,k) 
! so, we will only keep the j>=i entries
!!
! The enumeration will correcpond to the following indexing look of the i,j,k values. Suppose we need to compute an Asvd array for m basis vectors. 
!
!for l = 1:m
!
!  j = l
!  for k = 1,l-1
!    for i = 1:j 
!      ADD ENTRY TO Asvd 
!    end for 
!  end for
!
!  k = l
!  for j=1,l
!    for i=1,j
!      ADD ENTRY TO Asvd
!    end for
!  end for 
!  
!end for         
! 
! when using Asvd, we will know m --- the number of basis vectors used. This will determine 
! the total number of entries we will need: (m*m*(m+1)/2)
! 
! The entries of Asvd may will be written in chunks. Each chunk will have a two single index values saved at the beginning of the file: 
! si_start of the first entrie and si_end of the last entrie. From this index we can figure the starting combinations of the (i,j,k) to 
! compute. Following the above single index enumeration algorithm, 
! the index i should be changed the fastest, then index j and then index k. The algorithm to figure out the starting position the loop 
! is below. 
!
! for k_start = 1:M
!    if (k_start+1)^2(k_start +2)/2 > si_start
!     exit
! end for     
! si = k_start^2(k_start+1)/2
! for k_end = k_start : M
!    if (k_end)^2(k_end + 1)/2 > si_end
!     exit
! end for     
!
!
!
!for l = k_start:k_end
!  
!  for k = k_start,l-1
!    j = k
!    for i = 1:j 
!      si=si+1
!      if s1>=si_start and si<= si_end then
!          ADD ENTRY TO Asvd 
!      end if 
!    end for 
!  end for
!
!  k = l   
!  for j=1,l
!    for i=1,j
!      si=si+1
!      if s1>=si_start and si<= si_end then
!          ADD ENTRY TO Asvd 
!      end if 
!    end for
!  end for 
!
!end for         

!
! When requesting the chunk to be computed, 




!
! The total storage is O(m^3).
! The indexing of the one dimensional array will be such as to allow for 
!
! 
! 
! evaluation of Boltzmann collision operator for a number of distribution function
! The distribution functions are read from the file with the filename prefix set up in the 
! DtaDrvnBparameters.dat. Then evaluation of the colision operator is called using MPI parallel 
! subrouines. This code is to be exectured on the master node (irank=0). The matchind subrouine 
! FlyTrap0DMk_Coll_Dta_DGV_MPI have to be called on the slave nodes to perform the evaluation of
! the collision operator. The subroutines on both the master and the slave nodes are dependent 
! on the set up subroutines that set up supplimentary arrays necessary for computations. Those are 
! called in InitDGV1D.
!
! sol_array is a global variable
!
! In the code there is a nbnumber of commented MPI barriers to debug the communication. These should match the barriers in FlyTrap0DMk_Coll_Dta_DGV_MPI
! 
! Also look for debug barriers in the code for master and slave nodes.
! 
! vriables svdname,chunkname,si_start,si_end are public variables of the module
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 

subroutine Mk_Coll_Ker_SVD_Basis_MACRO_Serial(pad)

use DGV_data_driven_boltz, only:si_start,si_end,SVfile_name,SVKrnl_name,SetDGVParams_DataDrivenBLZM,WriteCKrnl_DtaDrvnB,&
                                     ReadSingVect_DtaDrvnB
use DGV_collision_mod

real (DP), dimension (:,:), pointer :: SVects => null() ! array to keep singular vectors.
real (DP), dimension (:), pointer :: CKrnl => null()  ! array to keep the collision kernel components -- using single index and the conventions to enumerate as above. 
real (DP), dimension (:), allocatable :: coll ! 

!!!!!!!!!!!!!!!!!
integer ::  pad ! size of zero padding should be consistent with pre-copmputed FF transformed kernel 
!
integer (I4B) :: gg, zz, loc_alloc_stat, time_step, err_count ! to keep them allocation status
!
integer (I4B) :: nsvs, mm ! scrap variable to keep the number of available singular values annd number od components in the singular values
integer (I4B) :: k_start, k_end, si, i,j,k,m ! scrap indices
real :: pr_time_1, pr_time_2 ! variable to calculate the processor time 
!  
!!!!!!!!!!!!!!!!!
real (DP), dimension (:), pointer :: sv1 => null(),sv2 => null() ! scrap arrays to hold singular vectors

!!!!!!!!!!!!!!!!!
call SetDGVParams_DataDrivenBLZM("DtaDrvnBparameters.dat",0) ! this one only needs to be called at the master node.

!!!!!!!!!!!!!!!!!!
call ReadSingVect_DtaDrvnB(SVfile_name,SVects) ! This will allocate the array of singular vectors and use the space to read the singular vectors from a file on the hard drive. 
                             ! Name of the file is in variable svdname 
                             ! Values of singular vecotors are in array pointed to by SVects 

!!!!!!!!!!!!!!!!!!!!! CHECK. Not an essential piece !!!!!!!!!!!!!!!!!!!!!!!
!!!!! Check that we read correctly from the file. Can remove later 
mm=size(SVects,1)   ! number of velocity point in each VDF
nsvs=size(SVects,2) ! number of components in the singular vectors. Should be equal to Mu*su*Mv*sv*Mw*sw 

!!!! Allocate the array to keep values of the collision kernel
allocate (CKrnl(si_start:si_end), coll(1:mm), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "Make_Collision_Kernel_SVD_Basis_MACRO: Allocation error for variables (CKrnl,coll)"
 end if 
 CKrnl = 0.0
 !!!!
 ! ready to start feeling in the collision kernel array
 !!!!!!!!!!!!!!!!

! try to measure time to run 
call cpu_time (pr_time_1) ! End the account of time
 
!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!! Main loop 
!!! The loop will actually be several loops. 
!!! We need to find the triples of indices i,j,k that correspond the 
!!! Specific value of the single index and then use it to populate the array CKrnl
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do k_start = 1,nsvs
   if ((k_start+1)*(((k_start+1)*(k_start+2))/2) > si_start) then 
       exit
   end if 
end do    
si = (k_start-1)*((k_start*(k_start-1))/2)
do k_end = k_start,nsvs
   if (k_end*((k_end*(k_end + 1))/2) >= si_end) then 
      exit
   end if 
end do     
!!! we now know that index k will be between k_start and k_end which will 
!!! give us an idea on what components to compute
!!!!!!!!!!!!!!!!!!!!!
! DEBUG
print *,"values k_start,k_end,si_start,si_end, si", k_start, k_end, si_start, si_end, si
! end DEBUG
!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!! DEBUG !!!!!!!!!!
!print *, "processor irank=", irank, "mark 1 "
!call flush(6)
!call MPI_Barrier(MPI_COMM_WORLD, ierr) 
!!!!!!!!!!!! END DEBUG !!!!!!!!!!!!

do m = k_start,k_end ! we only focus on these values of k
   ! for each new value of k we add "back slab" and then the "right side slab"
   ! Back slab goes first: 
   j = m
   do k = 1,m-1
    do i = 1,m 
      si=si+1
      if ((si>=si_start) .and. (si<= si_end)) then
         !!!! Evaluate the entry of the collision kernel 
         ! may will need to apply decomposition for the first vector of some other combinations of vectors????
         sv1 => SVects(:,i)
         sv2 => SVects(:,j)
         ! call EvalFourierColTwoF_OMP(sv1,sv2,coll) !! CHANGE TO AN APPROPRIATE CALL
         call EvalFourierColTwoZeroPadF_OMP(sv1,sv2,coll,pad) !!  CHANGE TO AN APPROPRIATE CALL
         sv1 => SVects(:,k)
         CKrnl(si) = sum(coll*sv1)
         !CKrnl(si) = sum(coll*SVects(:,k))
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         call flush(6)
         print *,"Make_Collision_Kernel_SVD_Basis_MACRO: (B) make a collision kernel entry i,j,k,si,CKrnl(si):",i,j,k,si,CKrnl(si)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end if 
    end do 
   end do 
   ! right side slab is next  
   k = m   
   do j=1,k
    do i=1,j
      si=si+1
      if ((si >= si_start) .and. (si <= si_end)) then
         !!!! Evaluate the entry of the collision kernel 
         ! may will need to apply decomposition for the first vector of some other combinations of vectors????
         sv1 => SVects(:,i)
         sv2 => SVects(:,j)
         ! call EvalFourierColTwoF_OMP(sv1,sv2,coll) !! CHANGE TO AN APPROPRIATE CALL
         call EvalFourierColTwoZeroPadF_OMP(sv1,sv2,coll,pad) !!  CHANGE TO AN APPROPRIATE CALL
         sv1 => SVects(:,k)
         CKrnl(si) = sum(coll*sv1)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         call flush(6)
         print *,"Make_Collision_Kernel_SVD_Basis_MACRO: (A) make a collision kernel entry i,j,k,si,CKrnl(si):", i,j,k, si,CKrnl(si) 
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end if 
    end do 
   end do  
  ! all done for this value of m
end do         
!
! Write the values of the collision kernel on the hard drive 
call WriteCKrnl_DtaDrvnB(SVKrnl_name,CKrnl,si_start,si_end) 
call flush(6)
print *,"Make_Collision_Kernel_SVD_Basis_MACRO: Wrote the CKrnl on the hard drive. si_start, CKrnl(si_start,si_start+5):", & 
                                si_start, CKrnl(si_start:si_start+5) 
nullify(sv1,sv2)
deallocate(CKrnl,coll,SVects)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call cpu_time (pr_time_2) ! End the account of time
print *, "Make_Collision_Data_Array_Solutions_MACRO: Processor time lapsed in seconds for main run:", pr_time_2 - pr_time_1 ! Report the account of time
! we now test evaluation of the basis functions:
end subroutine Mk_Coll_Ker_SVD_Basis_MACRO_Serial


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Mk_Coll_LinKer_SVZrBasis_MACRO_Serial(pad)
!
! This subroutine will compute the linear component of the collision operator for use with 
! SVD Zero Basis. Recall that SVD Zero Basis is obtained by first performing SVD on solution data. 
! Then we add exact maxwellian as the first vector of the basis , then we orthogonalize, then we 
! perturb the basis vectors from #1 to the lastso that they have zero mass, momentum and energy. 
! This basis is convenienc to defive the main variable to be the deviation from Maxwellian
! 
! The collision operator can be expressed as 
! 
!\partial_{t} e_{k} = \sum_{k''} e_{k''} \hat{B}_{kk''} + \sum_{k',k''} e_{k'}e_{k''}\hat{A}_{k',k'',k}
!
! Where  \hat{A}_{j,j,k} is stored in the array Asvd(i,j,k) where i,j,k run over 
! indices of the SVD singular vectors. 
! 
! Asvd(i,j,k) = SVD PROJ((mesh related const)*\int_{R^3} \int_{R^3}\int_{R^3} v_{i} (u-\xi) v_{j} (w-\xi) A(u, w) \v_{k}(\xi) du dw d\xi )
!
! Here is a better formula for Asvd(k',k'',k)
!
! \hat{A}_{k',k'',k} = \sum_{j} H_{kj} ( \sum_{j',j''} H_{k',j'-j} H_{k'',j''-j} A_{j',j''})
!
! Similarly, the linear portion uses matrix \hat{B}_{k,k''} defiend as follows:
!
!\hat{B}_{k,k''} = 2 \sum_{j} H_{kj} ( \sum_{j',j''} f_{M,j'-j} H_{k'',j''-j} )
! where f_{M,j'-j} is the discrete exact maxwellian
!
! \hat{B}_{k,k''} will be stored in the array BKrnl(k,k'')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 

subroutine Mk_Coll_LinKer_SVZrBasis_MACRO_Serial(pad)

use DGV_commvar, only: nodes_u, nodes_v, nodes_w
use DGV_data_driven_boltz, only:SVfile_name,SVLinKrnl_name,WriteBKrnl_DtaDrvnB,ReadSingVect_DtaDrvnB,SetDGVParams_DataDrivenBLZM
use DGV_collision_mod
use DGV_distributions_mod

real (DP), dimension (:,:), pointer :: SVects => null() ! array to keep singular vectors.
real (DP), dimension (:,:), pointer :: BKrnl => null()  ! array to keep the collision kernel components -- using single index and the conventions to enumerate as above. 
real (DP), dimension (:), allocatable :: coll, fm ! 

!!!!!!!!!!!!!!!!!
integer ::  pad ! size of zero padding should be consistent with pre-copmputed FF transformed kernel 
!
integer (I4B) :: gg, zz, loc_alloc_stat, time_step, err_count ! to keep them allocation status and such
!
integer (I4B) :: nsvs, mm ! scrap variable to keep the number of available singular values annd number od components in the singular values
integer (I4B) :: i,j,k,m ! scrap indices
real :: pr_time_1, pr_time_2 ! variable to calculate the processor time 
!  
!!!!!!!!!!!!!!!!!
real (DP), dimension (:), pointer :: sv1 => null() ! scrap arrays to hold singular vectors

!!!!!!!!!!!!!!!!!
call SetDGVParams_DataDrivenBLZM("DtaDrvnBparameters.dat",0) ! this one only needs to be called at the master node.

!!!!!!!!!!!!!!!!!!
call ReadSingVect_DtaDrvnB(SVfile_name,SVects) ! This will allocate the array of singular vectors and use the space to read the singular vectors from a file on the hard drive. 
                             ! Name of the file is in variable svdname 
                             ! Values of singular vecotors are in array pointed to by SVects 

!!!!!!!!!!!!!!!!!!!!! CHECK. Not an essential piece !!!!!!!!!!!!!!!!!!!!!!!
!!!!! Check that we read correctly from the file. Can remove later 
mm=size(SVects,1)   ! number of velocity point in each VDF
nsvs=size(SVects,2) ! number of components in the singular vectors. Should be equal to Mu*su*Mv*sv*Mw*sw 

!!!! Allocate the array to keep values of the collision kernel
allocate (BKrnl(nsvs,nsvs), coll(1:mm), fm(1:mm), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "Mk_Coll_LinKer_SVZrBasis_MACRO_Serial: Allocation error for variables (BKrnl,coll,fm)"
 end if 
 BKrnl = 0.0_DP
 !!!!
 ! ready to start feeling in the collision kernel array
 !!!!!!!!!!!!!!!!

! try to measure time to run 
call cpu_time (pr_time_1) ! End the account of time
 
!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!! Main loop 
!!! The loop will actually be several loops. 
!!! We need to find the triples of indices i,j,k that correspond the 
!!! Specific value of the single index and then use it to populate the array CKrnl
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!! First we compute the appropriate Maxwellian.
!!!!!!! ATTENTION: hard coded macroparameters !!!!!!!
fm = maxwelveldist(0.2_DP,0.0_DP,0.0_DP,0.0_DP,2.0_DP,nodes_u, nodes_v, nodes_w) ! Note that there is a factor of 2 that needs to be introduced. Hence . 
                                                                                 ! density is doubled. 
!!!!!!! 

do i = 1,nsvs
 sv1 => SVects(:,i)
 ! call EvalFourierColTwoF_OMP(sv1,sv2,coll) !! CHANGE TO AN APPROPRIATE CALL
 call EvalFourierColTwoZeroPadF_OMP(fm,sv1,coll,pad) !!  CHANGE TO AN APPROPRIATE CALL
 do j = 1,nsvs
  sv1 => SVects(:,j)
  BKrnl(j,i) = sum(coll*sv1)
 end do    
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 print *,"Mk_Coll_LinKer_SVZrBasis_MACRO_Serial: make a linear collision kernel entry BKrnl(:,i), i, max, min :", &
                                            i, maxval(BKrnl(:,i)), minval(BKrnl(:,i))
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end do 
!!!!!!!!!!!!!!!
! Write the values of the linear collision kernel on the hard drive 
call WriteBKrnl_DtaDrvnB(SVLinKrnl_name,BKrnl) 
print *,"Mk_Coll_LinKer_SVZrBasis_MACRO_Serial: Wrote the BKrnl on the hard drive." 
nullify(sv1)
deallocate(BKrnl,coll,fm,SVects)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call cpu_time (pr_time_2) ! End the account of time
print *, "Mk_Coll_LinKer_SVZrBasis_MACRO_Serial: Processor time lapsed in seconds for main run:", pr_time_2 - pr_time_1 ! Report the account of time
! we now test evaluation of the basis functions:
end subroutine Mk_Coll_LinKer_SVZrBasis_MACRO_Serial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Mk_Coll_Ker_SVD_Basis_MACRO_Optml
!!
!! This subroutine computes components of the binary ROM collision kernel
!!
!! This subroutin is a modification of Mk_Coll_Ker_SVD_Basis_MACRO_Serial. The principal change is the sequence in which 
!! the component are computed. 
!!
!! The subroutine will compute a three dimensional array Asvd(i,j,k) where i,j,k run over 
!! indices of the SVD singular vectors. 
!! 
!!
!! Asvd(i,j,k) = (mech related const)*\int_{R^3} \int_{R^3}\int_{R^3} v_{i} (u-\xi) v_{j} (w-\xi) A(u, w) \v_{k}(\xi) du dw d\xi 
!!
!! The components are computed using the following loop: 
!!
!! k_e -- the target size of the ROM basis
!!
!! For i=i_b,i_e (to compute in chunks) 
!!  For j=1,i 
!!    Compute the collision operator C_{ij,\xi}=\int_{R^3}\int_{R^3} v_{i} (u-\xi) v_{j} (w-\xi) A(u, w) \v_{k}(\xi)
!!    on the discrete grid it will be C_{ij,\alpha}
!!    Then start pojecting. The final array is idexed by a single index with a very specific enumeration (See below) to match 
!!    that enumeration, the following two looks in k are used: 
!!    For k=1,i-1 
!!      Compute projection of C_{ij,\alpha} on v_{k},i.e, Asvd(i,j,k) = \sum_{\alpha} C_{i,j,\alpha} v_{k,\alpha}
!!      then compute the value of the single index: 
!!      s= (i-1)^2 i /2 +(k-1)i+j for that projection 
!!      record the component in the single index array
!!    end for k
!!    For k = i,k_e
!!      Compute projection of C_{ij,\alpha} on v_{k},i.e, Asvd(i,j,k) = \sum_{\alpha} C_{i,j,\alpha} v_{k,\alpha}
!!      then compute the value of the single index: 
!!      s = (k-1)^2 k /2 +(k-1)k + (i-1)i/2 + j = (k+1)k(k-1)/2 + (i-1)i/2 + j
!!      for that projection 
!!      record the component in the single index array
!!    end for k
!!  end for j 
!! end for i
!!
!!      
!! To increase scalability of the process, the Asvd will be stored in the memory and on the hard drive as a single index array.
!! 
!! we note that in the case of symmetric formulation of the Glerkin projectiuon of the collision operator, Asvd(i,j,k) = Asvd(j,i,k) 
!! so, we will only keep the j>=i entries
!!
!! The enumeration will correcpond to the following indexing look of the i,j,k values. Suppose we need to compute an Asvd array for m basis vectors. 
!!
!!for l = 1:m
!!
!!  j = l
!!  for k = 1,l-1
!!    for i = 1:j 
!!      ADD ENTRY TO Asvd 
!!    end for 
!!  end for
!!
!!  k = l
!!  for j=1,l
!!    for i=1,j
!!      ADD ENTRY TO Asvd
!!    end for
!!  end for 
!!  
!!end for         
!! 
!! when using Asvd, we will know m --- the number of basis vectors used. This will determine 
!! the total number of entries we will need: (m*m*(m+1)/2)
!! 
!! In this subroutine, only a portion of Asvd may be computed. However, the computed values will not 
!! correspond to consequtive values of single intex 
!! 
!! To enable coputation in portions, we will create the full size Asvd array using the target 
!! value k_end and fill it with zeros. Then, the computed values will overwrtie the zeros in the 
!! array in the entries corresponding to single index. If there is another portion of the array computed, 
!! e.g., for diffeent range of i, that other portion will have non-zero entries in different locations. 
!! to combine the two part we need to ADD the two arrays. (THis will be implemented in Python
!! 
!! 
!! The total storage is O(m^3).
!! The indexing of the one dimensional array will be such as to allow for 
!!
!! 
!! 
!! vriables svdname,chunkname,k_end, i_start, i_end are public variables of the module
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine Mk_Coll_Ker_SVD_Basis_MACRO_Optml(pad)

use DGV_data_driven_boltz, only: k_tgt,i_start,i_end,SVfile_name,SVKrnl_name,SetDGVParams_DataDrivenBLZM,&
                                 ReadSingVect_DtaDrvnB,WriteCKrnl_DtaDrvnB
use DGV_collision_mod

real (DP), dimension (:,:), pointer :: SVects => null() ! array to keep singular vectors.
real (DP), dimension (:), pointer :: CKrnl => null()  ! array to keep the collision kernel components -- using single index and the conventions to enumerate as above. 
real (DP), dimension (:), allocatable :: coll ! 

!!!!!!!!!!!!!!!!!
integer ::  pad ! size of zero padding should be consistent with pre-copmputed FF transformed kernel 
!
integer (I4B) :: gg, zz, loc_alloc_stat, time_step, err_count ! to keep them allocation status
!
integer (I4B) :: nsvs, mm ! scrap variable to keep the number of available singular values annd number od components in the singular values
integer (I4B) :: si, i,j,k ! scrap indices
real :: pr_time_1, pr_time_2 ! variable to calculate the processor time 
!  
!!!!!!!!!!!!!!!!!
real (DP), dimension (:), pointer :: sv1 => null(),sv2 => null() ! scrap arrays to hold singular vectors

!!!!!!!!!!!!!!!!!
call SetDGVParams_DataDrivenBLZM("DtaDrvnBparameters.dat",0) ! this one only needs to be called at the master node.

!!!!!!!!!!!!!!!!!!
call ReadSingVect_DtaDrvnB(SVfile_name,SVects) ! This will allocate the array of singular vectors and use the space to read the singular vectors from a file on the hard drive. 
                             ! Name of the file is in variable svdname 
                             ! Values of singular vecotors are in array pointed to by SVects 

!!!!!!!!!!!!!!!!!!!!! CHECK. Not an essential piece !!!!!!!!!!!!!!!!!!!!!!!
!!!!! Check that we read correctly from the file. Can remove later 
mm=size(SVects,1)   ! number of velocity point in each VDF
nsvs=size(SVects,2) ! number of components in the singular vectors. Should be equal to Mu*su*Mv*sv*Mw*sw 

!!!! Allocate the array to keep values of the collision kernel
allocate (CKrnl(1:k_tgt*((k_tgt*(k_tgt+1))/2)), coll(1:mm), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "Mk_Coll_Ker_SVD_Basis_MACRO_Optml: Allocation error for variables (CKrnl,coll)"
 end if 
 CKrnl = 0.0
 !!!!
 ! ready to start feeling in the collision kernel array
 !!!!!!!!!!!!!!!!

! try to measure time to run 
call cpu_time (pr_time_1) ! End the account of time
 
do i = i_start,i_end
 do j = 1,i
  !!!! Evaluate the entry of the collision kernel 
  ! may will need to apply decomposition for the first vector of some other combinations of vectors????
  sv1 => SVects(:,i)
  sv2 => SVects(:,j)
  ! call EvalFourierColTwoF_OMP(sv1,sv2,coll) !! CHANGE TO AN APPROPRIATE CALL
  call EvalFourierColTwoZeroPadF_OMP(sv1,sv2,coll,pad) !!  CHANGE TO AN APPROPRIATE CALL
  ! Now we need to project on all possible ROM basis functions. 
  !!!!!!!!!!!!!!!!!!!!!!!!!!! 
  !!! Main loop in k
  !!! The loop will actually be two loops. 
  !!! We need to find the single index that correspond to the triples of indices i,j,k and then 
  !!! the single index to populate the array CKrnl
  !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do k = 1,i-1
   ! compute the single index
   si = (i-1)*((i-1)*i)/2  + (k-1)*i + j
   ! do the projection 
   sv1 => SVects(:,k)
   CKrnl(si) = sum(coll*sv1)
!!!! Debug !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         call flush(6)
!         print *,"Mk_Coll_Ker_SVD_Basis_MACRO_Optml: (A) make a collision kernel entry i,j,k,si,CKrnl(si):",i,j,k,si,CKrnl(si)
!!!!! end debug !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end do 
  !!! done with this portion of projections
  do k = i,k_tgt
   ! compute the single index
   si = (k-1)*((k*(k+1))/2) + ((i-1)*i)/2  + j
   ! do the projection 
   sv1 => SVects(:,k)
   CKrnl(si) = sum(coll*sv1)
!!!! Debug !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         call flush(6)
!         print *,"Mk_Coll_Ker_SVD_Basis_MACRO_Optml: (B) make a collision kernel entry i,j,k,si,CKrnl(si):",i,j,k,si,CKrnl(si)
!!!!! end debug !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end do 
  !!!!!!!!!!! Finished with projections for this combination of i and j 
!!!! Debug !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         call flush(6)
         print *,"Mk_Coll_Ker_SVD_Basis_MACRO_Optml: done with combination i,j",i,j
         print *,"Mk_Coll_Ker_SVD_Basis_MACRO_Optml: last record k,si,CKrnl(si):",k,si,CKrnl(si)
!!!!! end debug !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 end do 
end do 
   
! Write the values of the collision kernel on the hard drive 
call WriteCKrnl_DtaDrvnB(SVKrnl_name,CKrnl,1,k_tgt*((k_tgt*(k_tgt+1))/2)) 
!!!!!! Debug !!!!!!!!!!
call flush(6)
print *,"Mk_Coll_Ker_SVD_Basis_MACRO_Optml: Wrote the CKrnl on the hard drive. si, CKrnl(si):", & 
                                si, CKrnl(si) 
!!!!!!!!!!! end debug !!!!!!!!!!!
nullify(sv1,sv2)
deallocate(CKrnl,coll,SVects)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call cpu_time (pr_time_2) ! End the account of time
print *, "Mk_Coll_Ker_SVD_Basis_MACRO_Optml: Processor time lapsed in seconds for main run:", pr_time_2 - pr_time_1 ! Report the account of time
! we now test evaluation of the basis functions:
end subroutine Mk_Coll_Ker_SVD_Basis_MACRO_Optml


end module DGV_dta_drvn_bol_macro