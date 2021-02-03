program gsl_oro_data

!--------------------------------------------------------------------
! This program calls subroutines which calculate the parameters
! required for the subgrid-scale orographic gravity-wave drag (GWDO)
! scheme on the FV3 grid.  These parameters are for the small-scale
! GWD and Beljaars (2004) turbulent orographic form drag (TOFD)
! schemes of the GSL drag suite.  The output fields are:
! var, con, ol{1,2,3,4} and oa{1,2,3,4}
! or in FV3 parlance:
! stddev, convexity, ol{1,2,3,4} and oa{1,2,3,4}
! These variables are output to netCDF.
!
! Note:  This program works for both the global FV3GFS cubed
!        sphere, i.e., for tiles 1 through 6, and for the
!        regional lam, i.e., for tile 7
!
! Author:  Michael Toy -- NOAA/GSL   January 27, 2021
! Based on code by Michael Duda provided by NCAR/MMM
!--------------------------------------------------------------------

use gsl_oro_data_sm_scale, only: calc_gsl_oro_data_sm_scale
use gsl_oro_data_lg_scale, only: calc_gsl_oro_data_lg_scale
use gsl_oro_data_sm_scale_sar, only: calc_gsl_oro_data_sm_scale_sar
use gsl_oro_data_lg_scale_sar, only: calc_gsl_oro_data_lg_scale_sar

implicit none


character(len=2) :: tile_num   ! tile number entered by user
character(len=7) :: res_indx   ! grid-resolution index, e.g., 96, 192, 384, 768,
                               ! etc. entered by user
integer :: tile_num_int

logical :: duplicate_oro_data_file   ! flag for whether oro_data_ls file is a duplicate
                   ! of oro_data_ss due to minimum grid size being less than 7.5km



! Read in FV3GFS grid info
print *
print *, "Enter tile number:"
read (5,*) tile_num
print *
print *, "Enter grid-resolution index:"
read (5,*) res_indx
print *
print *, "Creating tile oro_data for tile number: ", tile_num
print *, "Grid resolution = ", res_indx
print *


read(tile_num,*) tile_num_int  ! integer form of tile_num

if ( tile_num_int.le.6 ) then  ! tile is a global tile

   call calc_gsl_oro_data_sm_scale(tile_num,res_indx,duplicate_oro_data_file)

   print *, "duplicate_oro_data_file =", duplicate_oro_data_file
   print *

   if ( .not.duplicate_oro_data_file ) then
      call calc_gsl_oro_data_lg_scale(tile_num,res_indx)
   end if

elseif ( tile_num_int.eq.7 ) then  ! tile is a regional tile

   ! print *, "this is a regional tile"
   ! print *

   call calc_gsl_oro_data_sm_scale_sar(tile_num,res_indx,            &
                                       duplicate_oro_data_file)

   print *, "duplicate_oro_data_file =", duplicate_oro_data_file

   if ( .not.duplicate_oro_data_file ) then
      call calc_gsl_oro_data_lg_scale_sar(tile_num,res_indx)
   end if

else

   print *, "Error: ", tile_num, " is not a valid tile number"

end if


print *
print *, "End program gsl_oro_data"
print *


end program gsl_oro_data
