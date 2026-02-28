!> @file
!! @brief Define and allocate required grid variables
!! @author Denise.Worthen@noaa.gov
!!
!> This module contains the grid variables
!! @author Denise.Worthen@noaa.gov

module grdvars

  use gengrid_kinds, only : dbl_kind, real_kind, int_kind

  implicit none

  real(kind=dbl_kind), parameter ::      pi = 3.14159265358979323846_dbl_kind  !< the value of PI
  real(kind=dbl_kind), parameter :: deg2rad = pi/180.0_dbl_kind                !< degree to radian conversion
  real(kind=dbl_kind), parameter ::  rearth = 6371.0_dbl_kind                  !< earth radius (km)

  integer :: ni                                                    !< i-dimension of output grid
  integer :: nj                                                    !< j-dimension of output grid
  integer :: npx                                                   !< i or j-dimension of fv3 tile

  integer :: nx                                                    !< i-dimension of MOM6 supergrid
  integer :: ny                                                    !< j-dimension of MOM6 supergrid

  logical :: editmask                                              !< flag indicating whether the MOM6 land mask
                                                                   !! should be edited. Default is false.
  logical :: debug                                                 !< flag indicating whether grid information
                                                                   !! should be printed for debugging purposes
                                                                   !! Default is false.
  logical :: do_postwgts                                           !< flag indicating whether then ESMF weights to
                                                                   !! regrid from the tripole grid to a rectilinear
                                                                   !! grid should be generated. Default is false.
  logical :: roottask                                              !< flag indicating whether this is the roottask

  integer, parameter :: nv = 4                                     !< the number of vertices for each stagger location
  integer, parameter :: ncoord = 2*4                               !< the number of coord pairs (lat,lon) for each of
                                                                   !! 4 stagger locations
  integer, parameter :: nverts = 2*4                               !< the number of coord pairs (lat,lon) for the
                                                                   !! vertices of each stagger location
  integer, parameter ::  nvars = ncoord + nverts                   !< the total number of cooridinate variables


  real(dbl_kind)     :: sg_maxlat                                  !< the maximum latitute present in the supergrid
                                                                   !! file
  integer(int_kind)  :: ipole(2)                                   !< the i-index for both pole locations
                                                                   !! along the top-most row

  integer, parameter, dimension(nv) :: iVertNE = (/0, -1, -1,  0/) !< The i-offsets defining the vertices of a point on the NE-Arakawa C-grid
  integer, parameter, dimension(nv) :: jVertNE = (/0,  0, -1, -1/) !< The j-offsets defining the vertices of a point on the NE-Arakawa C-grid

  type :: staggers
     real(dbl_kind), allocatable :: lat(:,:)       !< The latitude of the center grid points of a C-grid location
     real(dbl_kind), allocatable :: lon(:,:)       !< The longitudee of the center grid points of a C-grid location
     real(dbl_kind), allocatable :: latvert(:,:,:) !< The latitude of the corners grid points of a C-grid location
     real(dbl_kind), allocatable :: lonvert(:,:,:) !< The longitude of the corners grid points of a C-grid location
     real(dbl_kind), allocatable :: xlat(:)        !< The latitude of a stagger location at either j=0 or j=jmax+1
     real(dbl_kind), allocatable :: xlon(:)        !< The longitude of a stagger location at either j=0 or j=jmax+1
     integer, allocatable        :: iVert(:)       !< The i-index off-set array defining the indices on the stagger grid
                                                   ! which provides the vertices.
     integer, allocatable        :: jVert(:)       !< The j-index off-set array defining the indices on the stagger grid
                                                   ! which provides the vertices
                                                   ! Bu grid->Ct vertices, Ct grid->Bu vertices
                                                   ! Cu grid->Cv vertices, Cv grid->Cu vertices
  end type staggers
  type(staggers) :: Ct, Cu, Cv, Bu

  type :: grid
     type(staggers) :: Ct
     type(staggers) :: Cu
     type(staggers) :: Cv
     type(staggers) :: Bu
     ! MOM6 fields
     real(dbl_kind),  allocatable, dimension(:,:) :: areaCt !< The grid areas of the Ct grid cell in m2
     real(dbl_kind),  allocatable, dimension(:,:) :: anglet !< The rotation angle on Ct points (opposite sense from angle)
     real(dbl_kind),  allocatable, dimension(:,:) :: angle  !< The rotation angle on Bu points
     real(dbl_kind),  allocatable, dimension(:,:) :: angchk !< The rotation angle on Ct points, as calculated by CICE
                                                            !! internally using angle on Bu
     real(dbl_kind),  allocatable, dimension(:) :: xangCt   !< The rotation angle on the Ct grid points on the opposite
                                                            !! side of the tripole seam
     real(real_kind), allocatable, dimension(:,:) :: wet4   !< The ocean mask from a MOM6 mask file, stored as real*4 (nd)
     real(dbl_kind),  allocatable, dimension(:,:) :: wet8   !< The ocean mask from a MOM6 mask file, stored as real*8 (nd)
     real(real_kind), allocatable, dimension(:,:) :: dp4    !< The ocean depth from a MOM6 topog file, stored as real*4 (m)
     real(dbl_kind),  allocatable, dimension(:,:) :: dp8    !< The ocean depth from a MOM6 topog file, stored as real*8 (m)
     ! CICE6 fields
     real(dbl_kind),  allocatable, dimension(:,:) :: ulon   !< The longitude points (on the Bu grid) for CICE6
                                                            !! (radians)
     real(dbl_kind),  allocatable, dimension(:,:) :: ulat   !< The latitude points (on the Bu grid) for CICE6
                                                            !! (radians)
     real(dbl_kind),  allocatable, dimension(:,:) ::  htn   !< The grid cell width in centimeters of the CICE6
                                                            !! grid in the x-direction (i-dimension)
     real(dbl_kind),  allocatable, dimension(:,:) ::  hte   !< The grid cell width in centimeters of the CICE6
                                                            !! grid in the y-direction (j-dimension)
  end type grid

  real(kind=real_kind), parameter :: minimum_depth = 9.5           !< The minimum depth for MOM6
  real(kind=real_kind), parameter :: maximum_depth = 6500.0        !< The maximum depth for MOM6
  real(kind=real_kind), parameter :: masking_depth = 0.0           !< The masking depth for MOM6. Depths shallower than
                                                                   !! minimum_depth but deeper than masking_depth are
                                                                   !! rounded to minimum_depth
  real(kind=real_kind), parameter :: maximum_lat = 88.0            !< The maximum latitude for water points for WW3

  ! ATM resolutions
  integer, parameter :: maxatmres = 10                             !< The maximum number of possible ATM resolutions
  integer, allocatable, dimension(:) :: catm                       !< The ATM resolutions for mapped ocean masks

contains
  !> Allocate grid variables
  !!
  !! @author Denise Worthen

  subroutine allocate_all()

    allocate(grid%Ct%lat(ni,nj), grid%Ct%lon(ni,nj))
    allocate(grid%Cu%lat(ni,nj), grid%Cu%lon(ni,nj))
    allocate(grid%Cv%lat(ni,nj), grid%Cv%lon(ni,nj))
    allocate(grid%Bu%lat(ni,nj), grid%Bu%lon(ni,nj))

    allocate(grid%Ct%iVert(nv), grid%Ct%jVert(nv))
    allocate(grid%Cu%iVert(nv), grid%Cu%jVert(nv))
    allocate(grid%Cv%iVert(nv), grid%Cv%jVert(nv))
    allocate(grid%Bu%iVert(nv), grid%Bu%jVert(nv))
    allocate(grid%areaCt(ni,nj), grid%anglet(ni,nj), grid%angle(ni,nj), grid%angchk(ni,nj))

    allocate(grid%Ct%latvert(ni,nj,nv), grid%Ct%lonvert(ni,nj,nv))
    allocate(grid%Cu%latvert(ni,nj,nv), grid%Cu%lonvert(ni,nj,nv))
    allocate(grid%Cv%latvert(ni,nj,nv), grid%Cv%lonvert(ni,nj,nv))
    allocate(grid%Bu%latvert(ni,nj,nv), grid%Bu%lonvert(ni,nj,nv))

    allocate(grid%Ct%xlon(ni), grid%Ct%xlat(ni))
    allocate(grid%Cu%xlon(ni), grid%Cu%xlat(ni))
    allocate(grid%Cv%xlon(ni), grid%Cv%xlat(ni))
    allocate(grid%Bu%xlon(ni), grid%Bu%xlat(ni))
    allocate(grid%xangCt(ni))

    allocate(grid%wet4(ni,nj))
    allocate(grid%wet8(ni,nj))

    allocate(grid%dp4(ni,nj))
    allocate(grid%dp8(ni,nj))

    allocate(grid%ulon(ni,nj), grid%ulat(ni,nj))
    allocate(grid%htn(ni,nj), grid%hte(ni,nj))

  end subroutine allocate_all

end module grdvars
