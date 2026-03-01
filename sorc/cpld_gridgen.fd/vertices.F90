!> @file
!! @brief Fill the vertices for any stagger location
!! @author Denise.Worthen@noaa.gov
!!
!> This module fills the vertices for any stagger location. The i,j indices of the source lat and lon
!! arrays are modified by the vertex offsets to give the latitudes and longitues of each vertex for the
!! desired stagger lcation.
!! @author Denise.Worthen@noaa.gov
module vertices

  use gengrid_kinds, only : dbl_kind
  use grdvars,       only : nv

  implicit none

  private

  public :: fill_vertices

  interface fill_vertices
     module procedure fill_vertices_bottom
     module procedure fill_vertices_top
  end interface fill_vertices

  ! module local variables
  integer :: ib,ie,jb,je
  integer :: i,j,n,ii,jj

contains
  !> Fill the vertices for any stagger location, inclusive of bottom-most row
  !!
  !! @param[in] iVert,jVert        the i and j-offset applied to the indices of a stagger grid
  !! @param[in] lat,lon            the lats and lons of the stagger grid which define each vertex
  !! @param[in] latbot,lonbot      the lats and lons outside the bottom edge of the grid
  !! @param[out] latvert,latvert   the lats and lons of each vertex
  !! @author Denise.Worthen@noaa.gov
  subroutine fill_vertices_bottom(iVert,jVert,lat,lon,latbot,lonbot,latvert,lonvert,lb)

    integer,        intent(in)  :: lb
    integer,        intent(in)  :: iVert(nv), jVert(nv)
    real(dbl_kind), intent(in)  :: lat(:,:), lon(:,:)
    real(dbl_kind), intent(in)  :: latbot(:), lonbot(:)
    real(dbl_kind), intent(out) :: latvert(:,:,:), lonvert(:,:,:)

    ib = lbound(lat,1); ie = ubound(lat,1)
    jb = lbound(lat,2); je = ubound(lat,2)

    do j = jb,je
       do i = ib,ie
          do n = 1,nv
             ii = i + iVert(n); jj = j + jVert(n)
             if(ii .eq.    0)ii = ie
             if(ii .eq. ie+1)ii = 1
             if(jj .eq.    0) then
                latvert(i,j,n)   = latbot(ii)
                lonvert(i,j,n)   = lonbot(ii)
             else
                latvert(i,j,n)   = lat(ii,jj)
                lonvert(i,j,n)   = lon(ii,jj)
             end if
          enddo
       enddo
    enddo
  end subroutine fill_vertices_bottom

  !> Fill the vertices for any stagger location, inclusive of top-most row
  !!
  !! @param[in] iVert,jVert        the i and j-offset applied to the indices of a stagger grid
  !! @param[in] lat,lon            the lats and lons of the stagger grid which define each vertex
  !! @param[in] lattop,lontop      the lats and lons outside the top edge of the grid
  !! @param[out] latvert,latvert   the lats and lons of each vertex
  !! @author Denise.Worthen@noaa.gov
  subroutine fill_vertices_top(iVert,jVert,lat,lon,lattop,lontop,latvert,lonvert)

    integer,        intent(in)  :: iVert(nv), jVert(nv)
    real(dbl_kind), intent(in)  :: lat(:,:), lon(:,:)
    real(dbl_kind), intent(in)  :: lattop(:), lontop(:)
    real(dbl_kind), intent(out) :: latvert(:,:,:),lonvert(:,:,:)

    ib = lbound(lat,1); ie = ubound(lat,1)
    jb = lbound(lat,2); je = ubound(lat,2)

    do j = jb,je
       do i = ib,ie
          do n = 1,nv
             ii = i + iVert(n); jj = j + jVert(n)
             if(ii .eq.    0)ii = ie
             if(ii .eq. ie+1)ii = 1
             if(jj .eq. je+1)then
                latvert(i,j,n)   = lattop(ii)
                lonvert(i,j,n)   = lontop(ii)
             else
                latvert(i,j,n)   = lat(ii,jj)
                lonvert(i,j,n)   = lon(ii,jj)
             endif
          enddo
       enddo
    enddo
  end subroutine fill_vertices_top
end module vertices
