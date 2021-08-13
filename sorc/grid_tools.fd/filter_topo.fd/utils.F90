!> @file
!! @brief Utility routines.
!! @author GFDL programmer

!> Module that contains general utility routines.
!!
!! @author GFDL programmer
 module utils

 implicit none

 public

 character(len=512) :: topo_file = "orog" !< Path/name of the topography (or orography) file.
 character(len=128) :: topo_field = "orog_filt" !< NetCDF record name of the filtered
                                                !! topography (or orography).
 character(len=128) :: mask_field = "slmsk" !< NetCDF record name of the land/sea mask.
 character(len=512) :: grid_file = "atmos_mosaic.nc" !< Path/name of the grid mosaic file.

 logical :: zero_ocean = .true.  !< If true, no diffusive flux into water/ocean
                                 !! area (preserve islands).
 logical :: nested = .false.     !< If true, process a global grid with a nest.
 logical :: regional = .false.   !< If true, process a stand-alone regional grid.

 integer :: grid_type = 0     !< Grid type. 0 for a gnomonic grid.
 
 real    :: stretch_fac = 1.0 !< Grid stretching factor.
 real    :: res = 48.         !< The 'CRES' resolution.

 contains

!> Read the program namelist file. Then, write the namelist
!! variables to standard output.
!!
!! @author GFDL Programmer
 subroutine read_namelist

 implicit none

 integer :: stdunit = 6, unit=7, io_status
 logical :: opened

 namelist /filter_topo_nml/ topo_file, topo_field, mask_field, grid_file, zero_ocean, &
        stretch_fac, res, nested, grid_type, regional

 do
   inquire( unit=unit, opened=opened )
   if( .NOT.opened )exit
   unit = unit + 1
   if( unit.EQ.100 )call handle_err(-1, 'Unable to locate unit number.' )
 end do

 open( unit=unit, file='input.nml', iostat=io_status )
 read( unit,filter_topo_nml, iostat=io_status )
 close(unit)

 if (io_status > 0) call handle_err(-1, 'Error reading input.nml')

 write (stdunit, nml=filter_topo_nml)

 end subroutine read_namelist

!> Prints an error message to standard output,
!! then halts program execution with a
!! bad status.
!!
!! @param[in] status Error status code.
!! @param[in] string Error message.
!! @author GFDL Programmer
 subroutine handle_err(status, string)

 implicit none

#include <netcdf.inc>

 integer,          intent(in) :: status
 character(len=*), intent(in) :: string
 character(len=256) :: errmsg

 if (status .ne. nf_noerr) then
   errmsg = nf_strerror(status)
   errmsg = trim(errmsg) // " " // trim(string)
   print *, "FATAL ERROR:"
   print *, trim(errmsg)
   error stop 'Stopped'
 endif

 end subroutine handle_err

 end module utils
