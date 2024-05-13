!> @file
!! @brief Create orographic (oro_data) files for use by GSL drag suite
!! @author Michael Toy, NOAA/GSL
!! @date 2021-03-12
!!
!! Program GSL_ORO_DATA
!!
!! This program calls subroutines which calculate the parameters
!! required for the GSL subgrid-scale orographic gravity-wave drag (GWDO)
!! suite on the FV3 grid.  These parameters are for the small-scale
!! GWD (Tsiringakis et al., 2017) and turbulent orographic form drag (TOFD)
!! (Beljaars et al., 2004) schemes of the GSL drag suite.
!! The output fields are:
!! - stddev      standard deviation of subgrid-scale topograpy
!! - convexity   convexity (kurtosis) of subgrid-scale topography
!! - ol{1,2,3,4} orographic effective lengths of subgrid-scale topography
!!   for 4 orientations: 1-westerly, 2-southerly, 3-southwesterly, 4-northwesterly
!! - oa{1,2,3,4} orographic asymmetries of subgrid-scale topography
!!   for 4 orientations: 1-westerly, 2-southerly, 3-southwesterly, 4-northwesterly
!!
!! Note:  This program works for both the global FV3GFS cubed
!!        sphere, i.e., for tiles 1 through 6, (and 7 if nested
!!        grid) (halo.eq.-999 for no halo), and for the stand-alone
!!        regional lam (tile 7 and halo.ne.-999)
!!        If a halo number is given, this is only to specify the
!!        Cxxx_grid.halox data used for input.  The oro_data files
!!        are always "halo0" output.
!!
!! Based on code by Michael Duda provided by NCAR/MMM

!> Brief description of program:  Creates orographic (oro_data) files
!! needed by the GSL drag suite physics parameterization 
!!
!! @author Michaei Toy, NOAA/GSL
!! @return 0 for success, error code otherwise.
program gsl_oro_data

use omp_lib

use gsl_oro_data_sm_scale, only: calc_gsl_oro_data_sm_scale, timef
use gsl_oro_data_lg_scale, only: calc_gsl_oro_data_lg_scale

implicit none

character(len=2) :: tile_num   ! tile number entered by user
character(len=7) :: res_indx   ! grid-resolution index, e.g., 96, 192, 384, 768,
                               ! etc. entered by user
character(len=4) :: halo       ! halo value entered by user (for input grid data)

logical :: duplicate_oro_data_file   ! flag for whether oro_data_ls file is a duplicate
                   ! of oro_data_ss due to minimum grid size being less than 7.5km

integer :: tid, nthreads

real :: tbeg, tend

! Read in FV3GFS grid info
print *
print *, "Enter tile number:"
read (5,*) tile_num
print *
print *, "Enter grid-resolution index:"
read (5,*) res_indx
print *
print *, "Enter halo number (-999 for no halo):"
read (5,*) halo
print *
print *, "Creating tile oro_data for tile number: ", tile_num
print *, "Grid resolution = ", res_indx
print *, "Halo = ", halo
print *

!$OMP PARALLEL PRIVATE(TID)
  tid = omp_get_thread_num()
  if (tid==0) then
    nthreads = omp_get_num_threads()
    print*,'Number of threads = ', nthreads
  endif
!$OMP END PARALLEL

tbeg=timef()
call calc_gsl_oro_data_sm_scale(tile_num,res_indx,halo,duplicate_oro_data_file)
tend=timef()

print*,'timing of calc_gsl_oro_data_sm_scale              ',tend-tbeg

print *, "duplicate_oro_data_file =", duplicate_oro_data_file
print *

if ( .not.duplicate_oro_data_file ) then
   tbeg=timef()
   call calc_gsl_oro_data_lg_scale(tile_num,res_indx,halo)
   tend=timef()
   print*,'timing of calc_gsl_oro_data_lg_scale              ',tend-tbeg
end if


print *
print *, "End program gsl_oro_data"
print *


end program gsl_oro_data
