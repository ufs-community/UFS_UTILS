!> @file
!! @brief unit test for angles
!! @author Denise.Worthen@noaa.gov
!!
!! Given a 5x5 block of coords, calculate the angle, anglet and angchk
!! and compare against values from MOM6 and CICE6 history output. The
!! 3 5x5 blocks are for the 1deg tripole grid. The first two blocks contain
!! the two polar points. The final block is located in the lower right quad
!! of the tripole region (i=300:304,j=306:310). Because the angle calculations
!! involve points outside of the 5x5 blocks (eg, i+1), not all 5x5 points are
!! checked.
!!
!! @author Denise.Worthen@noaa.gov

  program ftst_find_angq

    use gengrid_kinds, only : dbl_kind, real_kind
    use angles       , only : find_angq, find_ang, find_angchk

    implicit none

    integer       , parameter :: ni = 5, nj = 5, nblocks=3

    ! verification values from model history files, mx100
    ! blocks 1 and 2 center on the left and right tripole
    ! block 3 is in lower right quad of the arctic at (300:304,306:310)
    real(dbl_kind), dimension(ni,nj,nblocks) :: cice6_angle
    real(dbl_kind), dimension(ni,nj,nblocks) :: cice6_anglet
    real(dbl_kind), dimension(ni,nj,nblocks) :: mom6_sinrot

    ! block grid point values
    real(dbl_kind), dimension(ni,nj,nblocks) :: lonct, lonbu, latbu

    ! test values
    real(dbl_kind), dimension(ni,nj,nblocks) :: anglet, angle, angchk
    real(dbl_kind), dimension(ni,nblocks)    :: xangCt

    real(dbl_kind)   :: dmax, diff
    character(len=1) :: ck

    integer :: i,i2,j,k
    integer :: ipole(nblocks)
    logical :: debug = .false.

    !CICE6
    data ((cice6_angle(i,j,1), i=1,ni), j=1,nj) /                                &
         0.55923,      0.29966,      0.00000,      -0.29966,      -0.55923,      &
         0.70064,      0.38859,      0.00000,      -0.38859,      -0.70064,      &
         0.90943,      0.53371,      0.00000,      -0.53371,      -0.90943,      &
         1.20822,      0.75580,      0.00000,      -0.75580,      -1.20822,      &
         1.38117,      0.89328,      0.00000,      -0.89328,      -1.38117/
    data ((cice6_angle(i,j,2), i=1,ni), j=1,nj) /                                &
         0.55923,      0.29966,      0.00000,      -0.29966,      -0.55923,      &
         0.70064,      0.38859,      0.00000,      -0.38859,      -0.70064,      &
         0.90943,      0.53371,      0.00000,      -0.53371,      -0.90943,      &
         1.20822,      0.75580,      0.00000,      -0.75580,      -1.20822,      &
         1.38117,      0.89328,      0.00000,      -0.89328,      -1.38117/
    data ((cice6_angle(i,j,3), i=1,ni), j=1,nj) /                                &
         -1.14949,      -1.16044,      -1.17076,      -1.18049,      -1.18968,   &
         -1.18188,      -1.19217,      -1.20185,      -1.21097,      -1.21957,   &
         -1.21425,      -1.22383,      -1.23284,      -1.24131,      -1.24930,   &
         -1.24649,      -1.25534,      -1.26364,      -1.27145,      -1.27880,   &
         -1.27854,      -1.28661,      -1.29419,      -1.30130,      -1.30799 /

    data ((cice6_anglet(i,j,1), i=1,ni), j=1,nj) /                               &
         0.60483,      0.38996,       0.13521,      -0.13521,      -0.38996,     &
         0.73204,      0.48692,       0.17198,      -0.17198,      -0.48692,     &
         0.90492,      0.63285,       0.23026,      -0.23026,      -0.63285,     &
         1.13597,      0.85126,       0.32134,      -0.32134,      -0.85126,     &
         1.34026,      1.05945,       0.41175,      -0.41175,      -1.05945/
    data ((cice6_anglet(i,j,2), i=1,ni), j=1,nj) /                               &
         0.60483,      0.38996,       0.13521,      -0.13521,      -0.38996,     &
         0.73204,      0.48692,       0.17198,      -0.17198,      -0.48692,     &
         0.90492,      0.63285,       0.23026,      -0.23026,      -0.63285,     &
         1.13597,      0.85126,       0.32134,      -0.32134,      -0.85126,     &
         1.34026,      1.05945,       0.41175,      -0.41175,      -1.05945/
    data ((cice6_anglet(i,j,3), i=1,ni), j=1,nj) /                               &
         -1.12735,      -1.13895,      -1.14989,      -1.16020,      -1.16995,   &
         -1.16004,      -1.17100,      -1.18130,      -1.19102,      -1.20018,   &
         -1.19278,      -1.20303,      -1.21267,      -1.22174,      -1.23029,   &
         -1.22546,      -1.23498,      -1.24391,      -1.25231,      -1.26021,   &
         -1.25800,      -1.26675,      -1.27495,      -1.28264,      -1.28988 /

    ! MOM6 sin_rot
    data ((mom6_sinrot(i,j,1), i=1,ni), j=1,nj) /                                &
         -0.57231,      -0.38355,      -0.13608,       0.13608,       0.38355,   &
         -0.67404,      -0.47460,      -0.17342,       0.17342,       0.47460,   &
         -0.79354,      -0.60678,      -0.23239,       0.23239,       0.60678,   &
         -0.91250,      -0.79512,      -0.32692,       0.32692,       0.79512,   &
         -0.98949,      -0.97272,      -0.43486,       0.43486,       0.97272/
    data ((mom6_sinrot(i,j,2), i=1,ni), j=1,nj) /                                &
         -0.57231,      -0.38355,      -0.13608,       0.13608,       0.38355,   &
         -0.67404,      -0.47460,      -0.17342,       0.17342,       0.47460,   &
         -0.79354,      -0.60678,      -0.23239,       0.23239,       0.60678,   &
         -0.91250,      -0.79512,      -0.32692,       0.32692,       0.79512,   &
         -0.98949,      -0.97272,      -0.43486,       0.43486,       0.97272/
    data ((mom6_sinrot(i,j,3), i=1,ni), j=1,nj) /                                &
         0.90334,       0.90826,       0.91278,       0.91694,       0.92079,    &
         0.91689,       0.92120,       0.92516,       0.92881,       0.93216,    &
         0.92946,       0.93320,       0.93662,       0.93976,       0.94264,    &
         0.94103,       0.94420,       0.94711,       0.94977,       0.95221,    &
         0.95154,       0.95419,       0.95661,       0.95882,       0.96085 /

    ! lon Ct,Bu and latBu
    data ((lonct(i,j,1), i=1,ni), j=1,nj) /                                      &
         -245.26937,    -232.98519,    -218.04609,    -201.95391,    -187.01481, &
         -252.74561,    -239.00064,    -220.46665,    -199.53335,    -180.99936, &
         -262.76769,    -248.28139,    -224.73626,    -195.26374,    -171.71861, &
         -275.86250,    -263.23478,    -234.03925,    -185.96075,    -156.76522, &
         -291.65329,    -286.25383,    -263.72015,    -156.27985,    -133.74617 /
    data ((lonbu(i,j,1), i=1,ni), j=1,nj) /                                      &
         -242.68153,    -227.77943,    -210.00000,    -192.22057,    -177.31847, &
         -251.01234,    -233.49458,    -210.00000,    -186.50542,    -168.98766, &
         -262.99894,    -243.55611,    -210.00000,    -176.44389,    -157.00106, &
         -279.68531,    -263.47519,    -210.00000,    -156.52481,    -140.31469, &
         -300.00000,    -300.00000,    -210.00000,    -120.00000,    -120.00000 /
    data ((latbu(i,j,1), i=1,ni), j=1,nj) /                                      &
         88.99001,      89.10706,      89.14964,      89.10706,      88.99001,   &
         89.16921,      89.31628,      89.37292,      89.31628,      89.16921,   &
         89.31750,      89.50699,      89.58912,      89.50699,      89.31750,   &
         89.41886,      89.66093,      89.79818,      89.66093,      89.41886,   &
         89.45503,      89.72754,      90.00000,      89.72754,      89.45503 /

    data ((lonct(i,j,2), i=1,ni), j=1,nj) /                                      &
         -65.26937,     -52.98519,     -38.04609,     -21.95391,      -7.01481,  &
         -72.74561,     -59.00064,     -40.46665,     -19.53335,      -0.99936,  &
         -82.76769,     -68.28139,     -44.73626,     -15.26374,       8.28139,  &
         -95.86250,     -83.23478,     -54.03925,      -5.96075,      23.23478,  &
        -111.65329,    -106.25383,     -83.72015,      23.72015,      46.25383 /
    data ((lonbu(i,j,2), i=1,ni), j=1,nj) /                                      &
         -62.68153,     -47.77943,     -30.00000,     -12.22057,       2.68153,  &
         -71.01234,     -53.49458,     -30.00000,      -6.50542,      11.01234,  &
         -82.99894,     -63.55611,     -30.00000,       3.55611,      22.99894,  &
         -99.68531,     -83.47519,     -30.00000,      23.47519,      39.68531,  &
        -120.00000,    -120.00000,     -30.00000,      60.00000,      60.00000 /
    data ((latbu(i,j,2), i=1,ni), j=1,nj) /                                      &
         88.99001,      89.10706,      89.14964,      89.10706,      88.99001,   &
         89.16921,      89.31628,      89.37292,      89.31628,      89.16921,   &
         89.31750,      89.50699,      89.58912,      89.50699,      89.31750,   &
         89.41886,      89.66093,      89.79818,      89.66093,      89.41886,   &
         89.45503,      89.72754,      90.00000,      89.72754,      89.45503 /

    data ((lonct(i,j,3), i=1,ni), j=1,nj) /                                      &
         38.08144,      38.86965,      39.62025,      40.33605,      41.01961,   &
         39.67340,      40.41426,      41.11877,      41.78975,      42.42974,   &
         41.27339,      41.96446,      42.62075,      43.24505,      43.83985,   &
         42.87655,      43.51566,      44.12188,      44.69789,      45.24611,   &
         44.47789,      45.06321,      45.61778,      46.14417,      46.64470 /
    data ((lonbu(i,j,3), i=1,ni), j=1,nj) /                                      &
         39.26336,      40.00933,      40.71974,      41.39725,      42.04431,   &
         40.83545,      41.53354,      42.19746,      42.82988,      43.43321,   &
         42.41197,      43.06002,      43.67560,      44.26131,      44.81949,   &
         43.98813,      44.58428,      45.14993,      45.68757,      46.19946,   &
         45.55910,      46.10182,      46.61624,      47.10473,      47.56942 /
    data ((latbu(i,j,3), i=1,ni), j=1,nj) /                                      &
         80.97195,      80.70286,      80.43134,      80.15742,      79.88114,   &
         81.07710,      80.80479,      80.53019,      80.25334,      79.97425,   &
         81.17218,      80.89690,      80.61947,      80.33992,      80.05826,   &
         81.25741,      80.97942,      80.69942,      80.41741,      80.13341,   &
         81.33309,      81.05265,      80.77032,      80.48610,      80.19999 /

    print *,"Starting test of cpld_gridgen routine find_angq"

    angle = 0.0
    anglet = 0.0
    angchk = 0.0
    ipole = 0

    j = nj
    do i = 1,ni
       if(latBu(i,j,1) .eq. 90.00)ipole(1) = i
       if(latBu(i,j,2) .eq. 90.00)ipole(2) = i
    enddo

    !------------------------------------------
    ! find anglet and test against mom6 values
    !------------------------------------------

    call find_ang((/2,5/),(/2,5/),lonBu(:,:,1),latBu(:,:,1),lonCt(:,:,1),anglet(:,:,1))
    call find_ang((/2,5/),(/2,5/),lonBu(:,:,2),latBu(:,:,2),lonCt(:,:,2),anglet(:,:,2))
    call find_ang((/2,5/),(/2,5/),lonBu(:,:,3),latBu(:,:,3),lonCt(:,:,3),anglet(:,:,3))

    do k = 1,nblocks
       write(ck,'(i1.1)')k
       dmax = 1.0e-30
       do j = 1,nj
          do i = 1,ni
             if (abs(anglet(i,j,k)) .gt. 0.0)then
                diff  = abs(anglet(i,j,k) - asin(mom6_sinrot(i,j,k)))
                dmax = max(diff,dmax)
             end if
          end do
       end do
       call passfail(dmax, 1.0e-4, 'MOM6 anglet, block '//trim(ck))
    end do

    if (debug) then
       print *,'anglet '
       print *,'left'
       do j = 1,nj
          print '(5f15.6)',(anglet(i,j,1), i = 1,ni)
       end do
       print *,'right'
       do j = 1,nj
          print '(5f15.6)',(anglet(i,j,2), i = 1,ni)
       end do
       print *,'quad'
       do j = 1,nj
          print '(5f15.6)',(anglet(i,j,3), i = 1,ni)
       end do
    end if

    !--------------------------------------------------------
    !
    !---------------------------------------------------------

    xangCt = 0.0
    do i = 2,ni
       i2 = ipole(2)+(ipole(1)-i)+1
       xangCt(i,1) = -anglet(i2,nj,1)       ! angle changes sign across seam
       xangCt(i,2) = -anglet(i2,nj,2)
       xangCt(i,3) =  anglet(i,nj,3)
    end do

    !-------------------------------------------
    ! find angle and test against cice6 values
    !-------------------------------------------

    call find_angq((/2,5/),(/2,5/),xangCt(:,1),anglet(:,:,1),angle(:,:,1))
    call find_angq((/2,5/),(/2,5/),xangCt(:,2),anglet(:,:,2),angle(:,:,2))
    call find_angq((/2,5/),(/2,4/),xangCt(:,3),anglet(:,:,3),angle(:,:,3))
    ! reverse angle for CICE
    angle = -angle

    do k = 1,nblocks
       write(ck,'(i1.1)')k
       dmax = 1.0e-30
       do j = 1,nj
          do i = 1,ni
             if (abs(angle(i,j,k)) .gt. 0.0) then
                diff = abs(angle(i,j,k) - cice6_angle(i,j,k))
                dmax = max(diff,dmax)
             end if
          end do
       end do
       call passfail(dmax, 1.0e-4, 'CICE6 angle, block '//trim(ck))
    end do

    if (debug) then
       print *,'angle'
       print *,'left'
       do j = 1,nj
          print '(5f15.6)',(angle(i,j,1), i = 1,ni)
       end do
       print *,'right'
       do j = 1,nj
          print '(5f15.6)',(angle(i,j,2), i = 1,ni)
       end do
       print *,'quad'
       do j = 1,nj
          print '(5f15.6)',(angle(i,j,3), i = 1,ni)
       end do
    end if

    !-----------------------------------------------------------------------
    ! find anglet calculated by CICE and test against cice6 and mom6 values
    !-----------------------------------------------------------------------
    call find_angchk((/2,4/),(/3,4/),angle(:,:,1),angchk(:,:,1))
    call find_angchk((/2,4/),(/3,4/),angle(:,:,2),angchk(:,:,2))
    call find_angchk((/2,4/),(/3,4/),angle(:,:,3),angchk(:,:,3))
    ! reverse angle for MOM6
    angchk = -angchk

    do k = 1,nblocks
       write(ck,'(i1.1)')k
       dmax = 1.0e-30
       do j = 1,nj
          do i = 1,ni
             if (abs(angchk(i,j,k)) .gt. 0.0) then
                diff = abs(-angchk(i,j,k) - cice6_anglet(i,j,k))
                dmax = max(diff,dmax)
             end if
          end do
       end do
       call passfail(dmax, 1.0e-4, 'angchk vs CICE6 anglet, block '//trim(ck))
    end do

    do k = 1,nblocks
       write(ck,'(i1.1)')k
       dmax = 1.0e-30
       do j = 1,nj
          do i = 1,ni
             if (abs(angchk(i,j,k)) .gt. 0.0) then
                diff = abs(angchk(i,j,k) - asin(mom6_sinrot(i,j,k)))
                dmax = max(diff,dmax)
             end if
          end do
       end do
       call passfail(dmax, 1.2e-2, 'angchk vs MOM6 anglet, block '//trim(ck))
    end do

    if (debug) then
       print *,'angchk'
       print *,'left'
       do j = 1,nj
          print '(5f15.6)',(angchk(i,j,1), i = 1,ni)
       end do
       print *,'right'
       do j = 1,nj
          print '(5f15.6)',(angchk(i,j,2), i = 1,ni)
       end do
       print *,'quad'
       do j = 1,nj
          print '(5f15.6)',(angchk(i,j,3), i = 1,ni)
       end do
    end if
  end program

  subroutine passfail(dmax,tolerance,msg)

    use gengrid_kinds, only : dbl_kind

    implicit none

    real(dbl_kind),    intent(in) :: dmax
    real(dbl_kind),    intent(in) :: tolerance
    character(len=*),  intent(in) :: msg

    if (dmax .le. tolerance) then
       print '(a,2f12.8)','SUCCESS! '//trim(msg)//'  ',dmax,tolerance
    else
       print '(a,2f12.8)','FAIL! '//trim(msg)//'  ',dmax,tolerance
       stop 1
    endif
  end subroutine passfail
