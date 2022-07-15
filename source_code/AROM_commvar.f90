! 7/11/2022 Alex Alekseenko, CSUN, alexander.alekseenko@csun.edu 
!           Cristian Villatoro, cristian.villatoro.193@my.csun.edu
!  AROM_commvar.f90 -- module containing definitions of global variables for the library 
!  implementing ROM and Adaptive ROM models for collision operator. 
!						subroutines and modules of the library will access these variables. 
!						also, subroutines that use the library will access these variables. 
!                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!

module AROM_commvar
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Global ROM Variables 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


real (DP), dimension (:,:), allocatable :: sol_array ! arrays to store values of solution files that we read from the hard drive and the corresponding 
                                      ! array to store the outputs of the collsiion operator
real (DP), dimension (:,:), allocatable :: Projector ! These is the matrix that will do projection and recovery of the ROM model
real (DP), dimension (:,:), allocatable :: LinROMKrnl ! These is the matrix will contain linear kernel for the ROM model. Only is used if SV Zero Basis is used.

real (DP), dimension (:), allocatable :: ROMKrnl ! This is an array that will keep the components of the singular vector collision kernel

!!!!!!!!!!!!!!!!!!!!!!!!!!
character (len=132) :: name_solutions_file_read  ! name of the file that contains solutions/VDFs that will be read from the drive to prepare data 
character (len=132) :: SVKrnl_name ! the variable to store the name of the output collision kernel data file 
character (len=132) :: SVLinKrnl_name ! the variable to store the name of the Linear collision kernel data file -- used with the SV Zero Basis 
character (len=132) :: SVfile_name ! the variable to store the name of the file where svd vectors are stored.
integer (I4B) :: si_start ! the number of the first entry of the collision kernel array to compute -- using single inmdex and the above numbering conventions
integer (I4B) :: si_end   ! the number of the last entry of the collision kernel array to compute 
integer (I4B) :: num_SVKernl_chnks     ! components of signular Vector kernel may be computed in portions and stored in multiple files, call chunks. This parameter give the total number of chunks. 
                                       ! chunk enumeration starts from 0.
                                       
integer (I4B) :: k_tgt, i_start, i_end ! k_tgt - target size of the ROM basis. i_start and i_end --- parameters determining the portiopn of the ROM binary 
                                       ! kernel to be computed using the Mk_Coll_Ker_SVD_Basis_MACRO_Optml subroutine.                                        
logical :: flag_SVZeroBasisInUse    ! this variable will keep a flag indicating that SV  Zero Basis is in use -- this will call for different collision operator, compared to other bases.                                    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Adaptive ROM Variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module AROM_commvar

