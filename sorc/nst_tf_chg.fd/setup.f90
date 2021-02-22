!> @file
!!
!! This module contains variables to run nst_tf_chg.
!!
!! @author Xu Li @date 2017-03-13
module setup

  integer :: nsmth !< the number of 9-point smooth
  integer :: istyp
  
contains

  !> Initialize grided related variables to default values. Set
  !! defaults for observation related variables.
  !!
  !! @author Xu Li @date 2017-03-13
  subroutine init_setup
    implicit none

!   Initialize arrays used in namelist obs_input 
    nsmth = 0
    istyp = 0
  end subroutine init_setup
  
end module setup


