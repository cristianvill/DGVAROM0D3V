!
! TnrColl_mod.f90
!
!
! This module contains customized subroutines that are used to run the FFB03V code to produce training data 
! to learn the Boltzmann collision operator.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module TrnColl_Subs

use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
   implicit none

real (DP), parameter, private :: pi25DT = 3.141592653589793238462643d0


contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! WriteInitMxwls_0D3VTrnColl(maxwells)
!! 
!! This subroutine will save content of array Maxwells to the hard drive.
!! The array contains macroparameters of the Maxwells streams that 
!! was used to generate the initial data for the solution
!!!!!!!!!!!!!!!!!!!!!!!!!!''

subroutine WriteInitMxwls_0D3VTrnColl(maxwells) 

use BGKVDCF0D_readwrite

intrinsic Trim

real (DP), dimension (:) :: maxwells ! the array containing macroparamters of Maxwellians. 
                                     ! sequences of 5 macroparameters in order: density, temperature, bulk velocity , following by the 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                    
!
character (len=132) :: file_name ! the variable to store the file name 

! first, we prepare the file name to store the solution
call MakeBaseNameBGK0D2(file_name, "results/")
file_name = trim(Adjustl(file_name))//"_InitMaxwelsArry.dat"  ! this file will keep the array 
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="UNFORMATTED", access="SEQUENTIAL")

write (15) size(maxwells,1) ! record the length of the array
write (15) maxwells ! record the order of the time integrator.
! end save 
close (15)
end subroutine WriteInitMxwls_0D3VTrnColl  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteSolnColl_TrnColl
!
! This subroutine will damp on the disk a value of the solution and the corresponding value of the collision operator along with some 
! parameters.
!  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteSolnColl_TrnColl (suff,f,fcol,t)

use BGKVDCF0D_readwrite

intrinsic Trim
!
character (len=*) :: suff ! the string to hold the suffix
real (DP), dimension (:), intent(in) :: f,fcol ! the first array has the solution at time t, the second array has the collision operator for that solution 
real (DP), intent (in) :: t ! value of time when the solution is saved
!
character (len=132) :: file_name ! the variable to store the file name 

! first, we prepare the file name to store the solution
call MakeBaseNameBGK0D2(file_name, "results/")
file_name = trim(Adjustl(file_name))//"_time"//trim(Adjustl(suff))//"_SltnColl.dat"  ! this file will keep the array 
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="UNFORMATTED", access="SEQUENTIAL")

write (15) t ! record the current time and the used dt
write (15) size(f,1) ! record the length of the arrays.
! Next we save the solution .... 
write (15) f,fcol
! end save 
close (15)
end subroutine WriteSolnColl_TrnColl 

                                     

end module TrnColl_Subs 