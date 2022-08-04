!> @file
!! @brief Program setup.
!! @author George Gayno @date 2018

!> Set up program execution
!!
!! Public variables:
!!
!! Here 'input' indicates variables associated with the input source data 
!! and 'mdl' indicates variables associated with the fv3 model grid.
!!
!! To not process a surface field, set its 'input_file' variable to NULL.
!! However, vegetation type must always be processed as it defines
!! landice points.
!!
!! @author George Gayno @date 2018
 module program_setup

 implicit none

 private

 character(len=500), public   :: input_leaf_area_index_file = "NULL" !< File containing input leaf area index data.
 character(len=500), public   :: input_facsf_file = "NULL" !< File containing input fractional
                                                           !! coverage data for strong zenith angle
                                                           !! dependent albedo.
 character(len=500), public   :: input_substrate_temperature_file = "NULL" !< File containing input soil substrate temperature data.
 character(len=500), public   :: input_maximum_snow_albedo_file = "NULL" !< File containing input maximum snow albedo data.
 character(len=500), public   :: input_snowfree_albedo_file = "NULL" !< File containing input snow-free albedo data.
 character(len=500), public   :: input_slope_type_file = "NULL" !< File containing input slope type data.
 character(len=500), public   :: input_soil_type_file = "NULL" !< File containing input soil type data.
 character(len=500), public   :: input_vegetation_type_file = "NULL" !< File containing input vegetation type data.
 character(len=500), public   :: input_vegetation_greenness_file = "NULL" !< File containing input vegetation greenness data.
 character(len=500), public   :: mosaic_file_mdl = "NULL" !< Model grid mosaic file.
 character(len=500), public   :: orog_dir_mdl = "NULL" !< Directory containing the model grid orography files.
 character(len=500), public   :: orog_files_mdl(6) = "NULL" !< Model grid orography filenames.

 character(len=50), public    :: leaf_area_index_method='bilinear' !< Interpolation method for leaf area index. Conservative or bilinear (default).
 character(len=50), public    :: maximum_snow_albedo_method='bilinear' !< Interpolation method for max snow albedo. Conservative or bilinear (default).
 character(len=50), public    :: snowfree_albedo_method='bilinear' !< Interpolation method for snowfree albedo. Conservative or bilinear (default).
 character(len=50), public    :: vegetation_greenness_method='bilinear' !< Interpolation method for vegetation greenness. Conservative or bilinear (default).

 integer, public              :: halo = 0 !< Number of row/cols defining the lateral
                                          !! boundary halo. Used for regional nests.

 public :: read_setup_namelist

 contains

!> Read program setup namelist
!!
!! @param[in] localpet mpi task number
!! @author George Gayno @date 2018
 subroutine read_setup_namelist(localpet)

 use mpi

 implicit none

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
