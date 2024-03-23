!> @file
!! @brief Determine the rotation angle on center and
!! corner points
!! @author Denise.Worthen@noaa.gov
!!
!> This module finds the rotation angle for at both center and corner points
!! It utilizes the MOM6 function modulo_around_point
!! @author Denise.Worthen@noaa.gov

module angles

  use gengrid_kinds, only : dbl_kind, int_kind
  use grdvars,       only : debug

  implicit none

contains
  !> Find the rotation angle on center (Bu) grid points
  !!
  !! @param[in]  ni           the i-dimension of the grid
  !! @param[in]  nj           the j-dimension of the grid
  !! @param[in]  xangCt       the angle across the tripole seam
  !! @param[in]  anglet       the rotation angle on Ct points
  !! @param[out] angle        the rotation angle on Bu points
  !! @author Denise.Worthen@noaa.gov

  subroutine find_angq(ni,nj,xangCt,anglet,angle)

    integer       , intent(in)  :: ni,nj
    real(dbl_kind), intent(in)  :: xangCt(:)
    real(dbl_kind), intent(in)  :: anglet(:,:)
    real(dbl_kind), intent(out) :: angle(:,:)

    ! local variables
    integer :: i,j,i2

    real(dbl_kind) :: angle_0, angle_w, angle_s, angle_sw
    real(dbl_kind) :: p25 = 0.25

    !---------------------------------------------------------------------
    ! find the angle on corners using the same relationship CICE uses
    ! internally to calculate angles on Ct using angles on Bu
    !
    !           w-----------------0 Ct(i+1,j+1)
    !           |                 |
    !        ----------Bu(i,j)---------- Bu lies on seam at j=nj
    !           |                 |
    !   Ct(i,j) sw----------------s
    !
    !---------------------------------------------------------------------

    angle = 0.0
    do j = 2,nj
       do i = 1,ni-1
          if (j .lt. nj) then
             angle_0  = anglet(i+1,j+1)
             angle_w  = anglet(i,  j+1)
             angle_s  = anglet(i+1,j  )
             angle_sw = anglet(i  ,j  )
          else
             angle_0  = xangCt(i+1  )
             angle_w  = xangCt(i    )
             angle_s  = anglet(i+1,j)
             angle_sw = anglet(i,  j)
          end if
          angle(i,j) = atan2(p25*(sin(angle_0) + sin(angle_w) + sin(angle_s) + sin(angle_sw)), &
                             p25*(cos(angle_0) + cos(angle_w) + cos(angle_s) + cos(angle_sw)))

          if (abs(angle(i,j)) .le. 1.0e-10)angle(i,j) = 0.0
       enddo
    enddo
    angle(ni,:) = -angle(1,:)

  end subroutine find_angq

  !> Verify the rotation angle on center (Ct) grid points using angle on corner
  !! (Bu) grid points
  !!
  !! @param[in]  ni         the i-dimension of the grid
  !! @param[in]  nj         the j-dimension of the grid
  !! @param[in]  angle      the rotation angle on Bu points
  !! @param[out] angchk     the rotation angle on Ct points
  !! @author Denise.Worthen@noaa.gov
  subroutine find_angchk(ni,nj,angle,angchk)

    integer       , intent(in)  :: ni,nj
    real(dbl_kind), intent(in)  :: angle(:,:)
    real(dbl_kind), intent(out) :: angchk(:,:)

    ! local variables
    integer :: i,j
    real(dbl_kind) :: angle_0, angle_w, angle_s, angle_sw
    real(dbl_kind) :: p25 = 0.25

    !---------------------------------------------------------------------
    ! check: calculate anglet from angle on corners as CICE does internally.
    ! since angle changes sign between CICE and MOM6, (-1)*angchk ~ anglet
    !
    !               w-----------------0 Bu(i,j)
    !               |                 |
    !               |     Ct(i,j)     |
    !               |                 |
    !   Bu(i-1,j-1) sw----------------s
    !
    !---------------------------------------------------------------------

    angchk = 0.0
    do j = 2,nj
       do i = 2,ni
          angle_0  = angle(i  ,j  )
          angle_w  = angle(i-1,j  )
          angle_s  = angle(i,  j-1)
          angle_sw = angle(i-1,j-1)
          angchk(i,j) = atan2(p25*(sin(angle_0) + sin(angle_w) + sin(angle_s) + sin(angle_sw)), &
                              p25*(cos(angle_0) + cos(angle_w) + cos(angle_s) + cos(angle_sw)))
       enddo
    enddo
    angchk(1,:) = -angchk(ni,:)

  end subroutine find_angchk

  !> Find the rotation angle on center (Ct) grid points
  !!
  !! @param[in]  ni            the i-dimension of the grid
  !! @param[in]  nj            the j-dimension of the grid
  !! @param[in]  lonBu         the longitudes of the corner grid points
  !! @param[in]  latBu         the latitudes of the corner grid points
  !! @param[in]  lonCt         the longitudes of the center grid points
  !! @param[out] anglet        the rotation angle on Ct points
  !! @author Denise.Worthen@noaa.gov

  subroutine find_ang(ni,nj,lonBu,latBu,lonCt,anglet)

    integer       , intent(in)  :: ni,nj
    real(dbl_kind), intent(in)  :: lonBu(:,:)
    real(dbl_kind), intent(in)  :: latBu(:,:)
    real(dbl_kind), intent(in)  :: lonCt(:,:)
    real(dbl_kind), intent(out) :: anglet(:,:)

    ! local variables
    integer :: i,j,m,n
    integer :: ii,jj

    ! from geolonB fix in MOM6
    real(dbl_kind) :: len_lon ! The periodic range of longitudes, usually 360 degrees.
    real(dbl_kind) :: pi_720deg ! One quarter the conversion factor from degrees to radians.
    real(dbl_kind) :: lonB(2,2)  ! The longitude of a point, shifted to have about the same value.
    real(dbl_kind) :: lon_scale = 0.0

    !---------------------------------------------------------------------
    ! rotation angle for "use_bugs" = false case from MOM6
    ! src/initialization/MOM_shared_initialization.F90 but allow for not
    ! having halo values
    ! note this does not reproduce sin_rot,cos_rot found in MOM6 output
    ! differences are ~O 10-6
    !---------------------------------------------------------------------

    anglet = 0.0
    pi_720deg = atan(1.0) / 180.0
    len_lon = 360.0
    do j=1,nj; do i = 1,ni
       do n=1,2 ; do m=1,2
          jj = J+n-2; ii = I+m-2
          if(jj .eq. 0)jj = 1
          if(ii .eq. 0)ii = ni
          lonB(m,n) = modulo_around_point(LonBu(ii,jj), LonCt(i,j), len_lon)
          !  lonB(m,n) = modulo_around_point(LonBu(I+m-2,J+n-2), LonCt(i,j), len_lon)
       enddo; enddo
       jj = j-1; ii = i-1
       if(jj .eq. 0)jj = 1
       if(ii .eq. 0)ii = ni
       lon_scale = cos(pi_720deg*((LatBu(ii,jj) + LatBu(I,J)) + &
            (LatBu(I,jj) + LatBu(ii,J)) ) )
       anglet(i,j) = atan2(lon_scale*((lonB(1,2) - lonB(2,1)) + (lonB(2,2) - lonB(1,1))), &
            (LatBu(ii,J) - LatBu(I,jj)) + &
            (LatBu(I,J) - LatBu(ii,jj)) )

       !lon_scale = cos(pi_720deg*((LatBu(I-1,J-1) + LatBu(I,J)) + &
       !                           (LatBu(I,J-1) + LatBu(I-1,J)) ) )
       !anglet(i,j) = atan2(lon_scale*((lonB(1,2) - lonB(2,1)) + (lonB(2,2) - lonB(1,1))), &
       !               (LatBu(I-1,J) - LatBu(I,J-1)) + &
       !               (LatBu(I,J) - LatBu(I-1,J-1)) )
    enddo; enddo

  end subroutine find_ang
  ! -----------------------------------------------------------------------------
  !> Return the modulo value of x in an interval [xc-(Lx/2) xc+(Lx/2)]
  !! If Lx<=0, then it returns x without applying modulo arithmetic.
  !!
  !! From src/initialization/MOM_shared_initialization.F90:
  !! @param[in] x   Value to which to apply modulo arithmetic
  !! @param[in] xc  Center of modulo range
  !! @param[in] Lx  Modulo range width
  !! @return x_mod  Value x shifted by an integer multiple of Lx to be close to xc
  function modulo_around_point(x, xc, Lx) result(x_mod)
    use gengrid_kinds, only : dbl_kind

    real(dbl_kind), intent(in) :: x
    real(dbl_kind), intent(in) :: xc
    real(dbl_kind), intent(in) :: Lx
    real(dbl_kind) :: x_mod

    if (Lx > 0.0) then
       x_mod = modulo(x - (xc - 0.5*Lx), Lx) + (xc - 0.5*Lx)
    else
       x_mod = x
    endif
  end function modulo_around_point
end module angles
