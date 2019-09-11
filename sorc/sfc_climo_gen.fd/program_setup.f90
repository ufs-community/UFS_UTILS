 module program_setup

!--------------------------------------------------------------------------
! module documentation block
!
! Module: program setup
!   pgrmmr: gayno           org: w/np2           date: 2018
!
! Abstract: Set up program execution
!
! Usage: use program_setup
!
! Public Subroutines:
! -------------------
! read_setup_namelist          Reads configuration namelist
!
! Public variables:
! -----------------
!
! Here 'input' indicates variables associated with the input source data 
! and 'mdl' indicates variables associated with the fv3 model grid.
!
! To not process a surface field, set its 'input_file' variable to NULL.
! However, vegetation type must always be processed as it defines
! landice points.
!
! halo                         Number of row/cols defining the lateral
!                              boundary halo.  Used for regionanl
!                              nests.
! input_facsf_file             File containing input fractional
!                              coverage data for strong zenith angle
!                              dependent albedo.
! input_vegetation_greenness_  File containing input vegetation
! file                         greenness data.
! input_leaf_area_index_file   File containing input leaf area index   
!                              data.
! input_maximum_snow_albedo_   File containing input maximum snow
! file                         albedo data.
! input_snowfree_albedo_file   File containing input snow-free 
!                              albedo data.
! input_soil_type_file         File containing input soil type data.
! input_slope_type_file        File containing input slope type data.
! input_substrate_temperature_ File containing input soil substrate
! file                         temperature data.
! input_vegetation_type_file   File containing input vegetation type data.
! leaf_area_index_method       Interpolation method for leaf area index.
!                              Conservative or bilinear (default).
! maximum snow albedo_method   Interpolation method for max snow albedo.
!                              Conservative or bilinear (default).
! mosaic_file_mdl              Model grid mosaic file
! orog_dir_mdl                 Directory containing the model grid
!                              orography files.
! orog_files_mdl               Model grid orography filenames.
! snowfree_albedo_method       Interpolation method for snowfree albedo.
!                              Conservative or bilinear (default).
! vegetation_greenness_        Interpolation method for vegetation 
! method                       greenness.  Conservative or bilinear.
!                              Default is bilinear.
!--------------------------------------------------------------------------

 implicit none

 private

 character(len=500), public   :: input_leaf_area_index_file = "NULL"
 character(len=500), public   :: input_facsf_file = "NULL"
 character(len=500), public   :: input_substrate_temperature_file = "NULL"
 character(len=500), public   :: input_maximum_snow_albedo_file = "NULL"
 character(len=500), public   :: input_snowfree_albedo_file = "NULL"
 character(len=500), public   :: input_slope_type_file = "NULL"
 character(len=500), public   :: input_soil_type_file = "NULL"
 character(len=500), public   :: input_vegetation_type_file = "NULL"
 character(len=500), public   :: input_vegetation_greenness_file = "NULL"
 character(len=500), public   :: mosaic_file_mdl = "NULL"
 character(len=500), public   :: orog_dir_mdl = "NULL"
 character(len=500), public   :: orog_files_mdl(6) = "NULL"

 character(len=50), public    :: leaf_area_index_method='bilinear'
 character(len=50), public    :: maximum_snow_albedo_method='bilinear'
 character(len=50), public    :: snowfree_albedo_method='bilinear'
 character(len=50), public    :: vegetation_greenness_method='bilinear'

 integer, public              :: halo = 0

 public :: read_setup_namelist

 contains

 subroutine read_setup_namelist(localpet)

!-----------------------------------------------------------------------
!  subroutine documentation block
!
! Subroutine: read setup namelist
!   prgmmr: gayno          org: w/np2           date: 2018
!
! Abstract: Read program setup namelist
!
! Usage:  call read_setup_namelist (localpet)
!
!   input argument list:
!     localpet               mpi task number
!-----------------------------------------------------------------------

 implicit none

 include 'mpif.h'

 integer, intent(in)   :: localpet

 integer               :: ierr

 namelist /config/ input_facsf_file, input_substrate_temperature_file, &
                   input_maximum_snow_albedo_file, input_snowfree_albedo_file, &
                   input_slope_type_file, input_soil_type_file, &
                   input_leaf_area_index_file, input_vegetation_type_file, &
                   input_vegetation_greenness_file, mosaic_file_mdl, &
                   orog_dir_mdl, orog_files_mdl, halo, &
                   vegetation_greenness_method, leaf_area_index_method, &
                   maximum_snow_albedo_method, snowfree_albedo_method

 print*,"- READ SETUP NAMELIST, LOCALPET: ", localpet

 open(41, file="./fort.41", iostat=ierr, err=900)
 read(41, nml=config, iostat=ierr, err=901)
 close (41)

 return

 900 print*,'- FATAL ERROR OPENING CONFIG NAMELIST'
 print*,'- IOSTAT IS: ', ierr
 call mpi_abort(mpi_comm_world, 10, ierr)

 901 print*,'- FATAL ERROR READING CONFIG NAMELIST'
 print*,'- IOSTAT IS: ', ierr
 call mpi_abort(mpi_comm_world, 11, ierr)

 end subroutine read_setup_namelist

 end module program_setup
