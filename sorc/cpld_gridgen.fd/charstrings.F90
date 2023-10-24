!> @file
!! @brief Define required character string variables
!! @author Denise.Worthen@noaa.gov
!!
!> This module contains the character string variables
!! @author Denise.Worthen@noaa.gov

module charstrings

  use gengrid_kinds, only : CL,CM,CS

  implicit none

  character(len=CL) :: dirsrc                                               !< The source directory containing the fix files for MOM6
  character(len=CL) :: dirout                                               !< The directory where output files will be written
  character(len=CL) :: fv3dir                                               !< The directory containing the FV3 mosaic files
  character(len=CS) :: res                                                  !< The Ocean/Ice resolution, e.g. 500 (5deg), 100 (1deg)
  !! 050 (1/2deg), 025 (1/4deg)
  character(len=CS) :: atmres                                               !< The ATM resolution, e.g. C96, C192, C384
  character(len=CL) :: logmsg                                               !< An informational message

  character(len=CL) :: maskfile  = 'ocean_mask.nc'                          !< The name of the MOM6 mask file
  character(len=CS) :: maskname  = 'mask'                                   !< The variable name of the mask field
  character(len=CL) :: editsfile                                            !< The name of the topo edits file (resolution specific)

  character(len=CL) :: topofile                                             !< The name of the MOM6 bathymetry file
  character(len=CS) :: toponame  = 'depth'                                  !< The name of the bathymetry field

  character(len=CL) :: history                                              !< A documentation string
  character(len=CS) :: cdate                                                !< The date stamp of file creation

  character(len= 2), dimension(4) :: staggerlocs = (/'Ct','Cu','Cv','Bu'/)  !< The named stagger locations of the grid

end module charstrings
