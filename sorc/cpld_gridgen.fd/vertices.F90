!> @file
!! @brief Fill the vertices for any stagger location
!! @author Denise.Worthen@noaa.gov
!!
!> This module fills the vertices for any stagger location. The i,j indices of the source lat and lon
!! arrays are modified by the vertex offsets to give the latitudes and longitues of each vertex for the
!! desired stagger lcation. The routine fill_vertices will fill the vertex values of all non-boundary
!! rows for any stagger locations. For the Ct and Cu grids, the fill_bottom routine fills the bottom
!! most vertex values. For the Cv and Bu grids, the routine fill_top fills the topmost vertex values using
!! the values from across the tripole seam.
!! @author Denise.Worthen@noaa.gov

module vertices

  use gengrid_kinds, only : dbl_kind
  use grdvars,       only : ni,nj,nv

  implicit none

  contains
!> Fill the vertices for any stagger location between bounding j-rows
!!
!! @param[in] jbeg   the beginning row
!! @param[in] jend   the ending row
!! @param[in] iVert  the i-offset applied to the i-index of a stagger grid
!! @param[in] jVert  the j-offset applied to the j-index of a stagger grid
!! @param[in] lat    the latitudes of the stagger grid which define each vertex
!! @param[in] lon    the longitudes of the stagger grid which define each vertex
!! @param[out] latvert   the latitudes of each vertex
!! @param[out] lonvert   the longitudes of each vertex
!! @author Denise.Worthen@noaa.gov
  
  subroutine fill_vertices(jbeg,jend,iVert,jVert,lat,lon,latvert,lonvert)

                              integer, intent( in) :: jbeg,jend
                              integer, intent( in) :: iVert(nv), jVert(nv)
  real(dbl_kind), dimension(ni,nj),    intent( in) ::  lat, lon

  real(dbl_kind), dimension(ni,nj,nv), intent(out) :: latvert,lonvert

  ! local variables
  integer :: i,j,n,ii,jj

  do j = jbeg,jend
   do i = 1,ni
    do n = 1,nv
      ii = i + iVert(n); jj = j + jVert(n)
      if(ii .eq.    0)ii = ni
      if(ii .eq. ni+1)ii = 1
      latvert(i,j,n)   = lat(ii,jj)
      lonvert(i,j,n)   = lon(ii,jj)
    enddo
   enddo
  enddo
  end subroutine fill_vertices

!> Fill the vertices for a stagger location along the bottom j-row
!!
!! @param[in] iVert  the i-offset applied to the i-index of a stagger grid
!! @param[in] jVert  the j-offset applied to the j-index of a stagger grid
!! @param[in] lat    the latitudes of the stagger grid which define each vertex
!! @param[in] lon    the longitudes of the stagger grid which define each vertex
!! @param[in] dlat   the approximate latitude along the bottom-most row
!! @param[out] latvert   the latitudes of each vertex
!! @param[out] lonvert   the longitudes of each vertex
!! @author Denise.Worthen@noaa.gov

  subroutine fill_bottom(iVert,jVert,lat,lon,latvert,lonvert,dlat)

                              integer, intent( in) :: iVert(nv), jVert(nv)
  real(dbl_kind), dimension(ni,nj),    intent( in) ::  lat, lon
  real(dbl_kind), dimension(ni),       intent( in) ::  dlat

  real(dbl_kind), dimension(ni,nj,nv), intent(out) :: latvert,lonvert

  ! local variables
  integer :: i,j,n,ii,jj

  ! fill in grid bottom (j=1)
  ! vertices 1,2 are available
  ! vertices 3,4 must be set manually
      j = 1
   do i = 1,ni
    do n = 1,2
      ii = i + iVert(n); jj = j + jVert(n)
      if(ii .eq.    0)ii = ni
      if(ii .eq. ni+1)ii = 1
      latvert(i,j,n)   = lat(ii,jj)
      lonvert(i,j,n)   = lon(ii,jj)
    enddo
    do n = 3,4
      ii = i + iVert(n)
      if(ii .eq.    0)ii = ni
      if(ii .eq. ni+1)ii = 1
      latvert(i,j, n) =  dlat(ii)
    enddo
      lonvert(i,j, 3) = lonvert(i,j,2)
      lonvert(i,j, 4) = lonvert(i,j,1)
   enddo
   end subroutine fill_bottom

!> Fill the vertices for a stagger location along the top j-row
!!
!! @param[in] iVert  the i-offset applied to the i-index of a stagger grid
!! @param[in] jVert  the j-offset applied to the j-index of a stagger grid
!! @param[in] lat    the latitudes of the stagger grid which define each vertex
!! @param[in] lon    the longitudes of the stagger grid which define each vertex
!! @param[in] xlat   the latitude across the tripole seam
!! @param[in] xlon   the longitude across the tripole seam
!! @param[out] latvert   the latitudes of each vertex
!! @param[out] lonvert   the longitudes of each vertex
!! @author Denise.Worthen@noaa.gov

   subroutine fill_top(iVert,jVert,lat,lon,latvert,lonvert,xlat,xlon)

                              integer, intent( in) :: iVert(nv), jVert(nv)
  real(dbl_kind), dimension(ni,nj),    intent( in) ::  lat,  lon
  real(dbl_kind), dimension(ni),       intent( in) :: xlat, xlon
  
  real(dbl_kind), dimension(ni,nj,nv), intent(out) :: latvert,lonvert

  ! local variables
  integer :: i,j,n,ii,jj

  ! fill in grid top (j=nj)
  ! vertices 3,4 are available
  ! vertices 1,2 must be set manually using 'across seam' values
      j = nj
   do i = 1,ni
    do n = 3,4
      ii = i + iVert(n); jj = j + jVert(n)
      if(ii .eq.    0)ii = ni
      if(ii .eq. ni+1)ii = 1
      latvert(i,j,n)  = lat(ii,jj)
      lonvert(i,j,n)  = lon(ii,jj)
    enddo
    do n = 1,2
      ii = i + iVert(n)
      if(ii .eq.    0)ii = ni
      if(ii .eq. ni+1)ii = 1
      latvert(i,j,n)  = xlat(ii)
      lonvert(i,j,n)  = xlon(ii)
   enddo
  enddo
      !latCv_vert(i,j, 1) = latCv_vert(i,j,4)
      !latCv_vert(i,j, 2) = latCv_vert(i,j,3)
      !lonCv_vert(i,j, 1) = lonCv_vert(i,j,4)+240.d0
      !lonCv_vert(i,j, 2) = lonCv_vert(i,j,3)+240.d0
  end subroutine fill_top
end module vertices
