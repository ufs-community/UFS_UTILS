module gengrid_utils

  use gengrid_kinds, only : dbl_kind, int_kind
  use grdvars      , only : deg2rad, rearth, nv

implicit none

private

public calc_dist
public reshape_staggers

contains
  !> Calculate the distance between two lat/lon points
  !!
  !! @param[in]  lat1,lon1,lat2,lon2   !< the lat and lon of 2 points, in degrees
  !! @return     distance              !< the distance between the 2 points
  function calc_dist(lat1, lon1, lat2, lon2) result(distance)

    real(dbl_kind), intent(in) :: lat1, lon1, lat2, lon2
    real(dbl_kind) :: distance
    real(dbl_kind) :: dlat, dlon, a, c, phi1, phi2

    phi1 = lat1 * deg2rad
    phi2 = lat2 * deg2rad
    dlat = (lat2 - lat1) * deg2rad
    dlon = (lon2 - lon1) * deg2rad

    ! Haversine formula
    a = sin(dlat/2.0_dbl_kind)**2 + cos(phi1) * cos(phi2) * sin(dlon/2.0_dbl_kind)**2
    c = 2.0_dbl_kind * atan2(sqrt(a), sqrt(1.0_dbl_kind - a))

    distance = rearth * c
  end function calc_dist
    !> Get center and corner grid points for a given stagger location
  !!
  !! @param[in]  iind                    the start/end index in the i-dimension
  !! @param[in]  jind                    the start/end index in the j-dimension
  !! @param[in]  lon, lat, mask          2D lat,lon centers and mask for the stagger
  !! @param[in]  lonvert, latvert        3D lat,lon vertices  for the stagger
  !! @param[out] cnlons, cnlats, cnmask  1D center lons,lats and mask
  !! @param[out] crlons, crlats          2D corner lons/lats
  !!
  !! @author Denise.Worthen@noaa.gov
  subroutine reshape_staggers(iind, jind, lon, lat, mask, lonvert, latvert, cnlons, cnlats, cnmask, crlons, crlats)
    integer,           intent(in)  :: iind(:), jind(:)
    real(dbl_kind),    intent(in)  :: lon(:,:), lat(:,:)
    integer(int_kind), intent(in)  :: mask(:,:)
    real(dbl_kind),    intent(in)  :: lonvert(:,:,:), latvert(:,:,:)
    real(dbl_kind),    intent(out) :: cnlons(:), cnlats(:)
    integer(int_kind), intent(out) :: cnmask(:)
    real(dbl_kind),    intent(out) :: crlons(:,:), crlats(:,:)

    integer :: idim, jdim, n
    integer :: ib, ie, jb, je

    ib = iind(1); ie = iind(2)
    jb = jind(1); je = jind(2)
    idim = ie - ib + 1
    jdim = je - jb + 1

    cnlons = reshape(    lon(ib:ie, jb:je), (/idim*jdim/))
    cnlats = reshape(    lat(ib:ie, jb:je), (/idim*jdim/))
    cnmask = reshape(   mask(ib:ie, jb:je), (/idim*jdim/))
    crlats = reshape(latvert(ib:ie, jb:je, :), (/nv, idim*jdim/), order=(/2,1/))
    crlons = reshape(lonvert(ib:ie, jb:je, :), (/nv, idim*jdim/), order=(/2,1/))

  end subroutine reshape_staggers
end module gengrid_utils
