!> @file
!! @brief Utility routines.
!! @author George Gayno NOAA/EMC

!> Module containing utility routines.
!! @author George Gayno NOAA/EMC
 module utils

 public remap_coef

 contains
 
!> Generate the weights and index of the grids used in the bilinear
!! interpolation.
!!
!! This routine was taken from the forecast model -
!! ./atmos_cubed_sphere/tools/fv_treat_da_inc.f90.
!!
!! @param[in] is Start index in x-direction of the source array.
!! @param[in] ie End index in x-direction of the source array.
!! @param[in] js Start index in y-direction of the source array.
!! @param[in] je End index in y-direction of the source array.
!! @param[in] im x-dimension of the source array.
!! @param[in] jm y-dimension of the source array.
!! @param[in] lon 1-d array of longitudes (in radians).
!! @param[in] lat 1-d array of latitudes (in radians).
!! @param[in] agrid 2-d array for lon [agrid(:,:,1)] & lat
!! [agrid(:,:,2)] (in radians).
!! @param[out] s2c Bi-linear interpolation weights of the four nearby
!! grids of the source array.
!! @param[out] id1 Index 1 in x-direction of the nearby grids of
!! the source array.
!! @param[out] id2 Index 2 in x-direction of the nearby grids of
!! the source array.
!! @param[out] jdc Index in y-direction of the nearby grid of the
!! source array.
!! @author Xu Li
 SUBROUTINE REMAP_COEF( is, ie, js, je,&
      im, jm, lon, lat, id1, id2, jdc, s2c, agrid )

    implicit none
    integer, intent(in):: is, ie, js, je
    integer, intent(in):: im, jm
    real,    intent(in):: lon(im), lat(jm)
    real,    intent(out):: s2c(is:ie,js:je,4)
    integer, intent(out), dimension(is:ie,js:je):: id1, id2, jdc
    real,    intent(in):: agrid(is:ie,js:je,2)
    ! local:
    real :: rdlon(im)
    real :: rdlat(jm)
    real:: a1, b1
    real, parameter :: pi = 3.1415926
    integer i,j, i1, i2, jc, i0, j0
    do i=1,im-1
      rdlon(i) = 1. / (lon(i+1) - lon(i))
    enddo
    rdlon(im) = 1. / (lon(1) + 2.*pi - lon(im))

    do j=1,jm-1
      rdlat(j) = 1. / (lat(j+1) - lat(j))
    enddo

    ! * Interpolate to cubed sphere cell center
    do 5000 j=js,je

      do i=is,ie

        if ( agrid(i,j,1)>lon(im) ) then
          i1 = im;     i2 = 1
          a1 = (agrid(i,j,1)-lon(im)) * rdlon(im)
        elseif ( agrid(i,j,1)<lon(1) ) then
          i1 = im;     i2 = 1
          a1 = (agrid(i,j,1)+2.*pi-lon(im)) * rdlon(im)
        else
          do i0=1,im-1
            if ( agrid(i,j,1)>=lon(i0) .and. agrid(i,j,1)<=lon(i0+1) ) then
              i1 = i0;  i2 = i0+1
              a1 = (agrid(i,j,1)-lon(i1)) * rdlon(i0)
              go to 111
            endif
          enddo
        endif
111     continue

        if ( agrid(i,j,2)<lat(1) ) then
          jc = 1
          b1 = 0.
        elseif ( agrid(i,j,2)>lat(jm) ) then
          jc = jm-1
          b1 = 1.
        else
          do j0=1,jm-1
            if ( agrid(i,j,2)>=lat(j0) .and. agrid(i,j,2)<=lat(j0+1) ) then
              jc = j0
              b1 = (agrid(i,j,2)-lat(jc)) * rdlat(jc)
              go to 222
            endif
          enddo
        endif
222     continue

        if ( a1<0.0 .or. a1>1.0 .or.  b1<0.0 .or. b1>1.0 ) then
             write(*,*) 'gid=', i,j,a1, b1
        endif

        s2c(i,j,1) = (1.-a1) * (1.-b1)
        s2c(i,j,2) =     a1  * (1.-b1)
        s2c(i,j,3) =     a1  *     b1
        s2c(i,j,4) = (1.-a1) *     b1
        id1(i,j) = i1
        id2(i,j) = i2
        jdc(i,j) = jc
      enddo   !i-loop
5000 continue   ! j-loop

 END SUBROUTINE REMAP_COEF

 end module utils
