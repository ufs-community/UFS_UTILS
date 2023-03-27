!> @file
!! @brief Define the variables for output
!! @author Denise.Worthen@noaa.gov
!!
!> This module defines the attributes for variables written to the tripole, cice and scrip grid files
!! @author Denise.Worthen@noaa.gov

module vartypedefs

  use charstrings, only : CL, CM, CS

  implicit none

  integer, parameter :: maxvars = 20    !< The maximum number of variables written to a file

  type :: vardefs
     character(len=CM)   ::  var_name    !< A variable name
     character(len=CM)   :: long_name    !< A variable's long name
     character(len=CM)   :: unit_name    !< A variable's unit
     character(len= 2)   ::  var_type    !< A variable's type
     character(len=CM)   ::  vertices    !< A variable's vertices
  end type vardefs

  type(vardefs) ::    fixvars(maxvars)  !< Attribute definitions for the variables written to the main tripole file
  type(vardefs) ::   cicevars(maxvars)  !< Attribute definitions for the variables written to the CICE grid file
  type(vardefs) ::  scripvars(maxvars)  !< Attribute definitions for the variables written to any SCRIP file

contains

  !> Define the variables written to the tripole grid file
  !!
  !! @author Denise.Worthen@noaa.gov

  subroutine fixvars_typedefine

    ! local variables
    integer :: ii = 0

    !default
    fixvars(:)%var_type  = 'r8'
    fixvars(:)%vertices  = ''

    ii = ii + 1
    fixvars(ii)%var_name  = 'lonCt'
    fixvars(ii)%long_name = 'Longitude of center (Ct) points'
    fixvars(ii)%unit_name = 'degrees_east'
    fixvars(ii)%vertices  = 'lonCt_vert'

    ii = ii + 1
    fixvars(ii)%var_name  = 'latCt'
    fixvars(ii)%long_name = 'Latitude of center (Ct) points'
    fixvars(ii)%unit_name = 'degrees_north'
    fixvars(ii)%vertices  = 'latCt_vert'

    ii = ii + 1
    fixvars(ii)%var_name  = 'lonCv'
    fixvars(ii)%long_name = 'Longitude of meridional velocity (Cv) points'
    fixvars(ii)%unit_name = 'degrees_east'
    fixvars(ii)%vertices  = 'lonCv_vert'

    ii = ii + 1
    fixvars(ii)%var_name  = 'latCv'
    fixvars(ii)%long_name = 'Latitude of meridional velocity (Cv) points'
    fixvars(ii)%unit_name = 'degrees_north'
    fixvars(ii)%vertices  = 'latCv_vert'

    ii = ii + 1
    fixvars(ii)%var_name  = 'lonCu'
    fixvars(ii)%long_name = 'Longitude of zonal velocity (Cu) points'
    fixvars(ii)%unit_name = 'degrees_east'
    fixvars(ii)%vertices  = 'lonCu_vert'

    ii = ii + 1
    fixvars(ii)%var_name  = 'latCu'
    fixvars(ii)%long_name = 'Latitude of zonal velocity (Cu) points'
    fixvars(ii)%unit_name = 'degrees_north'
    fixvars(ii)%vertices  = 'latCu_vert'

    ii = ii + 1
    fixvars(ii)%var_name  = 'lonBu'
    fixvars(ii)%long_name = 'Longitude of corner (Bu) points'
    fixvars(ii)%unit_name = 'degrees_east'
    fixvars(ii)%vertices  = 'lonBu_vert'

    ii = ii + 1
    fixvars(ii)%var_name  = 'latBu'
    fixvars(ii)%long_name = 'Latitude of corner (Bu) points'
    fixvars(ii)%unit_name = 'degrees_north'
    fixvars(ii)%vertices  = 'latBu_vert'

    ii = ii + 1
    fixvars(ii)%var_name  = 'lonCt_vert'
    fixvars(ii)%long_name = 'Longitude Vertices of Ct points'
    fixvars(ii)%unit_name = 'degrees_east'

    ii = ii + 1
    fixvars(ii)%var_name  = 'latCt_vert'
    fixvars(ii)%long_name = 'Latitude Vertices of Ct points'
    fixvars(ii)%unit_name = 'degrees_north'

    ii = ii + 1
    fixvars(ii)%var_name  = 'lonCu_vert'
    fixvars(ii)%long_name = 'Longitude Vertices of Cu points'
    fixvars(ii)%unit_name = 'degrees_east'

    ii = ii + 1
    fixvars(ii)%var_name  = 'latCu_vert'
    fixvars(ii)%long_name = 'Latitude Vertices of Cu points'
    fixvars(ii)%unit_name = 'degrees_north'

    ii = ii + 1
    fixvars(ii)%var_name  = 'lonCv_vert'
    fixvars(ii)%long_name = 'Longitude Vertices of Cv points'
    fixvars(ii)%unit_name = 'degrees_east'

    ii = ii + 1
    fixvars(ii)%var_name  = 'latCv_vert'
    fixvars(ii)%long_name = 'Latitude Vertices of Cv points'
    fixvars(ii)%unit_name = 'degrees_north'

    ii = ii + 1
    fixvars(ii)%var_name  = 'lonBu_vert'
    fixvars(ii)%long_name = 'Longitude Vertices of Bu points'
    fixvars(ii)%unit_name = 'degrees_east'

    ii = ii + 1
    fixvars(ii)%var_name  = 'latBu_vert'
    fixvars(ii)%long_name = 'Latitude Vertices of Bu points'
    fixvars(ii)%unit_name = 'degrees_north'

  end subroutine fixvars_typedefine
  !> Define the variables written to the CICE grid file
  !!
  !! @author Denise.Worthen@noaa.gov

  subroutine cicevars_typedefine

    ! local variables
    integer :: ii = 0

    !default
    cicevars(:)%var_type = 'r8'
    cicevars(:)%vertices  = ''

    ii = ii + 1
    cicevars(ii)%var_name  = 'ulon'
    cicevars(ii)%long_name = 'Longitude of corner (Bu) points'
    cicevars(ii)%unit_name = 'radians'

    ii = ii + 1
    cicevars(ii)%var_name  = 'ulat'
    cicevars(ii)%long_name = 'Latitude of corner (Bu) points'
    cicevars(ii)%unit_name = 'radians'

    ii = ii + 1
    cicevars(ii)%var_name  = 'hte'
    cicevars(ii)%long_name = 'Distance between corner (Bu) points, east face'
    cicevars(ii)%unit_name = 'cm'

    ii = ii + 1
    cicevars(ii)%var_name  = 'htn'
    cicevars(ii)%long_name = 'Distance between corner (Bu) points, north face'
    cicevars(ii)%unit_name = 'cm'

    ii = ii + 1
    cicevars(ii)%var_name  = 'angle'
    cicevars(ii)%long_name = 'Angle at corner (Bu) points'
    cicevars(ii)%unit_name = 'radians'

    ii = ii + 1
    cicevars(ii)%var_name  = 'kmt'
    cicevars(ii)%long_name = 'ocean fraction at T-cell centers'
    cicevars(ii)%unit_name = 'none'
    cicevars(ii)%var_type  = 'i4'

  end subroutine cicevars_typedefine
  !> Define the variables written to any SCRIP grid file
  !!
  !! @author Denise.Worthen@noaa.gov

  subroutine scripvars_typedefine

    ! local variables
    integer :: ii = 0

    !default
    scripvars(:)%long_name  = ''
    scripvars(:)%var_type   = 'r8'
    scripvars(:)%vertices   = ''

    ii = ii + 1
    scripvars(ii)%var_name  = 'grid_center_lat'
    scripvars(ii)%unit_name = 'degrees'

    ii = ii + 1
    scripvars(ii)%var_name  = 'grid_center_lon'
    scripvars(ii)%unit_name = 'degrees'

    ii = ii + 1
    scripvars(ii)%var_name  = 'grid_corner_lat'
    scripvars(ii)%unit_name = 'degrees'

    ii = ii + 1
    scripvars(ii)%var_name  = 'grid_corner_lon'
    scripvars(ii)%unit_name = 'degrees'

  end subroutine scripvars_typedefine
end module vartypedefs
