! This software incorporates the MSIS empirical atmospheric model software
! designed and provided by NRL. Use is governed by the Open Source Academic
! research License Agreement contained in the file msis2.1/nrlmsis2.1_license..txt
  subroutine gettemp(iday,xlat,pr,np,pf,temp,n_o,n_o2,n_n2)
    use msis_init,      only: msisinit
    use msis_constants, only: rp
    use esmf,           only: esmf_kind_r8

    implicit none

    integer, intent(in) :: iday                   ! calender day
    real(kind=esmf_kind_r8), intent(in) :: xlat   ! latitude (degrees)
    real(kind=esmf_kind_r8), intent(in) :: pr(np) ! pressure (hPa)
    integer, intent(in) :: np                     ! number of pressure layers
    character(*), intent(in) :: pf                ! path to parmfile for msisinit
    real(kind=esmf_kind_r8), intent(out) :: temp(np) ! temperature (K)
    real(kind=esmf_kind_r8), intent(out) :: n_o(np)  ! number density of o
    real(kind=esmf_kind_r8), intent(out) :: n_o2(np) ! number density of o2
    real(kind=esmf_kind_r8), intent(out) :: n_n2(np) ! number density of n2
! Local variables
    real(kind=rp), parameter :: alt=100, ut=0, f107=150, f107a=150, ap(7)=9, xlong=0
    real(kind=rp) :: t, d(10), zkm
    integer       :: k,il,ip
    real(4)       :: switch_legacy(1:25)

! set swich 7,8,10,14 to zero to avoid diurnal changes in output tempe
! #7 is for diurnal, #8 for semidiurnal, #10 is for all ut/longitudinal
! effect, #14 is for terdiurnal
    switch_legacy(:)  = 1.
    switch_legacy(7)  = 0.
    switch_legacy(8)  = 0.
    switch_legacy(10) = 0.
    switch_legacy(14) = 0.

    call msisinit(parmfile=pf, switch_legacy=switch_legacy)
! calculate temperature, species, for each pres level
    do ip=1,np
      call ghp8(real(iday, rp), ut, alt, real(xlat, rp),   &
                xlong , f107a, f107, ap, real(pr(ip), rp), &
                zkm, d, t)
      temp(ip)=real(t,    esmf_kind_r8)
      n_n2(ip)=real(d(2), esmf_kind_r8)
      n_o2(ip)=real(d(3), esmf_kind_r8)
      n_o( ip)=real(d(4), esmf_kind_r8)
    enddo

  end subroutine gettemp

  subroutine ghp8(day,utsec,z0,glat,glon,f107a,f107,ap,pres,alt,dn,tn)

    use msis_constants, only: kB, NA, g0, rp
    use msis_calc,      only: msiscalc
    use msis_utils,     only: alt2gph
    use esmf,           only: esmf_kind_r8

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
