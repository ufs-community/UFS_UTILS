!> @file
!! @brief Define and allocate required grid variables
!! @author Denise.Worthen@noaa.gov
!!
!> This module contains the grid variables
!! @author Denise.Worthen@noaa.gov

module grdvars

  use gengrid_kinds, only : dbl_kind, real_kind, int_kind

  implicit none

  integer :: ni                                                    !< i-dimension of output grid
  integer :: nj                                                    !< j-dimension of output grid
  integer :: npx                                                   !< i or j-dimension of fv3 tile

  integer :: nx                                                    !< i-dimension of MOM6 supergrid
  integer :: ny                                                    !< j-dimension of MOM6 supergrid

  logical :: editmask                                              !< flag indicating whether the MOM6 land mask
                                                                   !! should be edited. Default is false
  logical :: debug                                                 !< flag indicating whether grid information
                                                                   !! should be printed for debugging purposes
                                                                   !! Default is false
  logical :: do_postwgts                                           !< flag indicating whether then ESMF weights to
                                                                   !! regrid from the tripole grid to a rectilinear
                                                                   !! grid should be generated. Default is false.
  logical :: roottask                                              !< flag indicating whether this is the roottask

  integer, parameter :: nv = 4.                                    !< the number of vertices for each stagger location
  integer, parameter :: ncoord = 2*4.                              !< the number of coord pairs (lat,lon) for each of
                                                                   !! 4 stagger locations
  integer, parameter :: nverts = 2*4.                              !< the number of coord pairs (lat,lon) for the
                                                                   !! vertices of each stagger location
  integer, parameter ::  nvars = ncoord + nverts                   !< the total number of cooridinate variables


  real(dbl_kind)     :: sg_maxlat                                  !< the maximum latitute present in the supergrid
                                                                   !! file
  integer(int_kind)  :: ipole(2)                                   !< the i-index for both pole locations
                                                                   !! along the top-most row

  integer, parameter, dimension(nv) :: iVertCt = (/0, -1, -1,  0/) !< The i-offsets of the Bu grid at each Ct(i,j)
                                                                   !! which determine the 4 vertices of each Ct grid
                                                                   !! grid point in i
  integer, parameter, dimension(nv) :: jVertCt = (/0,  0, -1, -1/) !< The j-offsets of the Bu grid at each Ct(i,j)
                                                                   !! which determine the 4 vertices of each Ct
                                                                   !! grid point in j
  integer, dimension(nv) :: iVertCv                                !< The i-offsets of the Cu grid at each Cv(i,j)
                                                                   !! which determine the 4 vertices of each Cv
                                                                   !! grid point in i
  integer, dimension(nv) :: jVertCv                                !< The j-offsets of the Cu grid at each Cv(i,j)
                                                                   !! which determine the 4 vertices of each Cv
                                                                   !! grid point in j
  integer, dimension(nv) :: iVertCu                                !< The i-offsets of the Cv grid at each Cu(i,j)
                                                                   !! which determine the 4 vertices of each Cu
                                                                   !! grid point in i
  integer, dimension(nv) :: jVertCu                                !< The j-offsets of the Cv grid at each Cu(i,j)
                                                                   !! which determine the 4 vertices of each Cu
                                                                   !! grid point in j
  integer, dimension(nv) :: iVertBu                                !< The i-offsets of the Ct grid at each Bu(i,j)
                                                                   !! which determine the 4 vertices of each Bu
                                                                   !! grid point in i
  integer, dimension(nv) :: jVertBu                                !< The j-offsets of the Ct grid at each Bu(i,j)
                                                                   !! which determine the 4 vertices of each Bu
                                                                   !! grid point in j
  ! Super-grid source grid variables
  real(dbl_kind), allocatable, dimension(:,:)   :: x               !< The longitudes of the MOM6 supergrid
  real(dbl_kind), allocatable, dimension(:,:)   :: y               !< The latitudes of the MOM6 supergrid
  real(dbl_kind), allocatable, dimension(:,:)   :: dx              !< The grid cell width in meters of the supergrid
                                                                   !! in the x-direction (i-dimension)
  real(dbl_kind), allocatable, dimension(:,:)   :: dy              !< The grid cell width in meters of the supergrid
                                                                   !! in the y-direction (j-dimension)

  ! Output grid variables
  real(dbl_kind), allocatable, dimension(:,:) :: latCt             !< The latitude of the center (tracer) grid points
                                                                   !! on the C-grid
  real(dbl_kind), allocatable, dimension(:,:) :: lonCt             !< The longitude of the center (tracer) grid
                                                                   !! points on the C-grid
  real(dbl_kind), allocatable, dimension(:,:) :: latCv             !< The latitude of the v-velocity grid points on
                                                                   !! the C-grid
  real(dbl_kind), allocatable, dimension(:,:) :: lonCv             !< The longitude of the v-velocity grid points on
                                                                   !! the C-grid
  real(dbl_kind), allocatable, dimension(:,:) :: latCu             !< The latitude of the u-velocity grid points on
                                                                   !! the C-grid
  real(dbl_kind), allocatable, dimension(:,:) :: lonCu             !< The longitude of the u-velocity grid points on
                                                                   !! the C-grid
  real(dbl_kind), allocatable, dimension(:,:) :: latBu             !< The latitude of the corner points on the C-grid.
                                                                   !! These are equivalent to u,v velocity grid
                                                                   !! points on the B-grid
  real(dbl_kind), allocatable, dimension(:,:) :: lonBu             !< The longitude of the corner points on the
                                                                   !! C-grid. These are equivalent to u,v velocity
                                                                   !! grid points on the B-grid
  real(dbl_kind), allocatable, dimension(:,:) :: areaCt            !< The grid areas of the Ct grid cell in m2
  real(dbl_kind), allocatable, dimension(:,:) :: anglet            !< The rotation angle on Ct points (opposite sense
                                                                   !! from angle)
  real(dbl_kind), allocatable, dimension(:,:) :: angle             !< The rotation angle on Bu points
  real(dbl_kind), allocatable, dimension(:,:) :: angchk            !< The rotation angle on Ct points, as calculated by
                                                                   !! CICE internally using angle on Bu

  real(dbl_kind), allocatable, dimension(:,:,:) :: latCt_vert      !< The latitudes of the 4 vertices of each Ct grid
                                                                   !! point
  real(dbl_kind), allocatable, dimension(:,:,:) :: lonCt_vert      !< The longitudes of the 4 vertices of each Ct
                                                                   !! grid point

  real(dbl_kind), allocatable, dimension(:,:,:) :: latCv_vert      !< The latitudes of the 4 vertices of each Cv grid
                                                                   !! point
  real(dbl_kind), allocatable, dimension(:,:,:) :: lonCv_vert      !< The longitudes of the 4 vertices of each Cv
                                                                   !! grid point

  real(dbl_kind), allocatable, dimension(:,:,:) :: latCu_vert      !< The latitudes of the 4 vertices of each Cu grid
                                                                   !! point
  real(dbl_kind), allocatable, dimension(:,:,:) :: lonCu_vert      !< The longitudes of the 4 vertices of each Cu
                                                                   !! grid point

  real(dbl_kind), allocatable, dimension(:,:,:) :: latBu_vert      !< The latitudes of the 4 vertices of each Bu grid
                                                                   !! point
  real(dbl_kind), allocatable, dimension(:,:,:) :: lonBu_vert      !< The longitudes of the 4 vertices of each Bu
                                                                   !! grid point


  real(dbl_kind), allocatable, dimension(:) :: xlonCt              !< The longitude of the Ct grid points on the
                                                                   !! opposite side of the tripole seam
  real(dbl_kind), allocatable, dimension(:) :: xlatCt              !< The latitude of the Ct grid points on the
                                                                   !! opposite side of the tripole seam
  real(dbl_kind), allocatable, dimension(:) :: xangCt              !< The rotation angle on the Ct grid points on the
                                                                   !! opposite side of the tripole seam

  real(dbl_kind), allocatable, dimension(:) :: xlonCu              !< The longitude of the Cu grid points on the
                                                                   !! opposite side of the tripole seam
  real(dbl_kind), allocatable, dimension(:) :: xlatCu              !< The latitude of the Cu grid points on the
                                                                   !! opposite side of the tripole seam
  real(dbl_kind), allocatable, dimension(:) :: dlatBu              !< The latitude spacing between Bu points at the
                                                                   !! grid bottom
  real(dbl_kind), allocatable, dimension(:) :: dlatCv              !< The latitude spacing between Cv points at the
                                                                   !! grid bottom
  ! MOM6 fix fields
  real(real_kind), allocatable, dimension(:,:) :: wet4             !< The ocean mask from a MOM6 mask file, stored as
                                                                   !! real*4 (nd)
  real(dbl_kind),  allocatable, dimension(:,:) :: wet8             !< The ocean mask from a MOM6 mask file, stored as
                                                                   !! real*8 (nd)

  real(real_kind), allocatable, dimension(:,:) :: dp4              !< The ocean depth from a MOM6 topog file, stored
                                                                   !! as real*4 (m)
  real(dbl_kind),  allocatable, dimension(:,:) :: dp8              !< The ocean depth from a MOM6 topog file, stored
                                                                   !! as real*8 (m)

  ! CICE6 fields
  real(dbl_kind), allocatable, dimension(:,:) :: ulon              !< The longitude points (on the Bu grid) for CICE6
                                                                   !! (radians)
  real(dbl_kind), allocatable, dimension(:,:) :: ulat              !< The latitude points (on the Bu grid) for CICE6
                                                                   !! (radians)
  real(dbl_kind), allocatable, dimension(:,:) ::  htn              !< The grid cell width in centimeters of the CICE6
                                                                   !! grid in the x-direction (i-dimension)
  real(dbl_kind), allocatable, dimension(:,:) ::  hte              !< The grid cell width in centimeters of the CICE6
                                                                   !! grid in the y-direction (j-dimension)

  real(kind=real_kind), parameter :: minimum_depth = 9.5           !< The minimum depth for MOM6
  real(kind=real_kind), parameter :: maximum_depth = 6500.0        !< The maximum depth for MOM6
  real(kind=real_kind), parameter :: masking_depth = 0.0           !< The masking depth for MOM6. Depths shallower than
                                                                   !! minimum_depth but deeper than masking_depth are
                                                                   !! rounded to minimum_depth
  real(kind=real_kind), parameter :: maximum_lat = 88.0            !< The maximum latitude for water points for WW3

contains
  !> Allocate grid variables
  !!
  !! @author Denise Worthen

  subroutine allocate_all

    allocate( x(0:nx,0:ny),  y(0:nx,0:ny) )
    allocate(  dx(nx,0:ny), dy(0:nx,ny) )

    allocate( latCt(ni,nj), lonCt(ni,nj) )
    allocate( latCv(ni,nj), lonCv(ni,nj) )
    allocate( latCu(ni,nj), lonCu(ni,nj) )
    allocate( latBu(ni,nj), lonBu(ni,nj) )

    allocate( areaCt(ni,nj), anglet(ni,nj), angle(ni,nj), angchk(ni,nj))

    allocate( latCt_vert(ni,nj,nv), lonCt_vert(ni,nj,nv) )
    allocate( latCv_vert(ni,nj,nv), lonCv_vert(ni,nj,nv) )
    allocate( latCu_vert(ni,nj,nv), lonCu_vert(ni,nj,nv) )
    allocate( latBu_vert(ni,nj,nv), lonBu_vert(ni,nj,nv) )

    allocate( xlonCt(ni), xlatCt(ni), xangCt(ni) )
    allocate( xlonCu(ni), xlatCu(ni) )
    allocate( dlatBu(ni), dlatCv(ni) )

    allocate( wet4(ni,nj) )
    allocate( wet8(ni,nj) )

    allocate(  dp4(ni,nj) )
    allocate(  dp8(ni,nj) )

    allocate( ulon(ni,nj), ulat(ni,nj) )
    allocate(  htn(ni,nj),  hte(ni,nj) )

  end subroutine allocate_all

end module grdvars
