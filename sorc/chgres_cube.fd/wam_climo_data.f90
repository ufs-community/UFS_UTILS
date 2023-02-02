! This software incorporates the MSIS empirical atmospheric model software
! designed and provided by NRL. Use is governed by the Open Source Academic
! research License Agreement contained in the file msis2.1/nrlmsis2.1_license..txt
  subroutine gettemp(iday,nday,xlat,nlat,pr,np,pf,temp,n_o,n_o2,n_n2)
    use msis_init,      only: msisinit
    use msis_constants, only: rp

    implicit none

    integer(kind=rp), intent(in) :: iday(nday)  ! calender day
    integer, intent(in) :: nday                 ! number of days
    real(kind=rp),    intent(in) :: xlat(nlat)  ! latitude in degree
    integer, intent(in) :: nlat                 ! number of latitudes
    real(kind=rp),    intent(in) :: pr(np)      ! pressure in mb
    integer, intent(in) :: np                   ! number of pressure layer
    character(*), intent(in) :: pf              ! path to parmfile for msisinit
    real,   intent(out) :: temp(np,nlat,nday)   ! temperature
    real,   intent(out) :: n_o(np,nlat,nday)    ! number density of o
    real,   intent(out) :: n_o2(np,nlat,nday)   ! number density of o2
    real,   intent(out) :: n_n2(np,nlat,nday)   ! number density of n2
    real                :: alt(np,nlat,nday)    ! altitude in km
    real(kind=rp)       :: d(9),t,ap(7),ut,xlong,f107,f107a,zkm
    integer             :: k,il,ip
    real(4)             :: switch_legacy(1:25)

! set swich 7,8,10,14 to zero to avoid diurnal changes in output tempe
! #7 is for diurnal, #8 for semidiurnal, #10 is for all ut/longitudinal
! effect, #14 is for terdiurnal
    switch_legacy(:)  = 1.
    switch_legacy(7)  = 0.
    switch_legacy(8)  = 0.
    switch_legacy(10) = 0.
    switch_legacy(14) = 0.

    call msisinit(parmfile=pf, switch_legacy=switch_legacy)
! set 10.7cm flux, Ap to an average value
    f107=150.
    f107a=150.
    ap=9.
! set longitude, ut, it should not make difference to output
    ut=0.
    xlong=0.
! default z0
    alt = 100.
! calculate temperature, species, for each lat,pres level,day
    do k=1,nday
      do il=1,nlat
        do ip=1,np
          call ghp8(iday(k),ut,alt(ip,il,k),xlat(il),xlong,f107a,f107,  &
                    ap,pr(ip),zkm,d,t)
          temp(ip,il,k)=t
          n_n2(ip,il,k)=d(2)
          n_o2(ip,il,k)=d(3)
          n_o( ip,il,k)=d(4)
        enddo
      enddo
    enddo

  end subroutine gettemp

  subroutine ghp8(day,utsec,z0,glat,glon,f107a,f107,ap,pres,alt,dn,tn)

    use msis_constants, only : kB, NA, g0, rp
    use msis_calc, only      : msiscalc
    use msis_utils, only     : alt2gph

    implicit none

    real(kind=rp),intent(in)    :: day
    real(kind=rp),intent(in)    :: utsec
    real(kind=rp),intent(in)    :: z0     !! first guess
    real(kind=rp),intent(in)    :: glat,glon
    real(kind=rp),intent(in)    :: f107a,f107
    real(kind=rp),intent(in)    :: ap(7)
    real(kind=rp),intent(in)    :: pres   !!! pressure in hPa
    real(kind=rp),intent(out)   :: alt
    real(kind=rp),intent(out)   :: dn(10)
    real(kind=rp),intent(out)   :: tn
! Local variables
    real, parameter    :: tol = 0.000043
    integer, parameter :: maxit = 30
    real               :: plog,delta
    real(kind=rp)      :: tex,zkm,pzkm
    real               :: xn,gz,xmbar,scl
    integer            :: n
    real(8)            :: xlat,alt0,alt1

    plog = log10(pres*100.0_rp)
    zkm = z0
    delta = 1.0_rp

    n = 0

    do while ( abs(delta) .ge. tol .and. n .le. maxit )
      n = n + 1

      call msiscalc(day,utsec,zkm,glat,glon,f107a,f107,ap,tn,dn,tex)

      xn = sum(dn(2:8))
      pzkm = kB * xn * tn
      delta = plog - log10(pzkm)
      xmbar = dn(1) / xn / 1.66E-24_rp
      xlat = dble(glat)
      alt0 = dble(zkm)
      alt1 = alt0 + 1.0d0
      gz = real((alt2gph(xlat,alt1) - alt2gph(xlat,alt0)) * g0)
      scl = Na * kB * tn / (xmbar * gz)

      ! difference
      zkm = zkm - scl * delta / 1000.0_rp
    end do
    alt = zkm

  end subroutine ghp8
