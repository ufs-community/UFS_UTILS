!> @file
!! @brief Snippets of noah model from sflx.F needed for land DA updates 
!! 
!! @author Clara Draper
      module sflx_snippet

      private

      public frh2o

      contains 

!> Calculate the liquid water (slc) for a given total moisture content 
!! and soil temperature. Used here to update slc when DA update to stc 
!! crosses freezing. Note this is an approximation, since in sflx.F (noah) 
!! the change in slc estimated by this routine is often clipped,
!! depending on the net energy input to the soil layer. Also, this 
!! routine is being called using stc, but in sflx.F it is called using 
!! the soil temp at the midpoint of the layer. However, testing shows 
!! the affects of these approximations is small (O(0.001 m3/m3)).
!! @param[in] tkelv Soil temperature in K
!! @param[in] smc Soil moisture 
!! @param[in] sh2o Input liquid soil moisture
!! @param[in] smcmax Max soil moisture
!! @param[in] bexp B exponent 
!! @param[in] psis Saturated matric potential
!! @param[out]  liqwat Output liquid soil moisture
      subroutine frh2o                                                  &
!  ---  inputs:
     &     ( tkelv, smc, sh2o, smcmax, bexp, psis,                      &
!  ---  outputs:
     &       liqwat                                                     &
     &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  subroutine frh2o calculates amount of supercooled liquid soil water  !
!  content if temperature is below 273.15k (t0).  requires newton-type  !
!  iteration to solve the nonlinear implicit equation given in eqn 17   !
!  of koren et al (1999, jgr, vol 104(d16), 19569-19585).               !
!                                                                       !
!  new version (june 2001): much faster and more accurate newton        !
!  iteration achieved by first taking log of eqn cited above -- less    !
!  than 4 (typically 1 or 2) iterations achieves convergence.  also,    !
!  explicit 1-step solution option for special case of parameter ck=0,  !
!  which reduces the original implicit equation to a simpler explicit   !
!  form, known as the "flerchinger eqn". improved handling of solution  !
!  in the limit of freezing point temperature t0.                       !
!                                                                       !
!  subprogram called:  none                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     tkelv    - real, temperature (k)                             1    !
!     smc      - real, total soil moisture content (volumetric)    1    !
!     sh2o     - real, liquid soil moisture content (volumetric)   1    !
!     smcmax   - real, saturation soil moisture content            1    !
!     bexp     - real, soil type "b" parameter                     1    !
!     psis     - real, saturated soil matric potential             1    !
!                                                                       !
!  outputs:                                                             !
!     liqwat   - real, supercooled liquid water content            1    !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  constant parameters:

! this block added from physconst.f for snippet

      implicit none

      real, parameter :: gs2     = 9.81        !< con_g in snowpack, frh2o
      real, parameter :: tfreez  = 2.7315e+2
      real, parameter :: lsubf   = 3.335e5     !< con_hfus=3.3358e+5
! end block added for snippet

      real, parameter :: ck    = 8.0
      real, parameter :: blim  = 5.5
      real, parameter :: error = 0.005

!  ---  inputs:
      real, intent(in) :: tkelv, smc, sh2o, smcmax, bexp, psis

!  ---  outputs:
      real, intent(out) :: liqwat

!  ---  locals:
      real :: bx, denom, df, dswl, fk, swl, swlk

      integer :: nlog, kcount
!
!===> ...  begin here
!
!  --- ...  limits on parameter b: b < 5.5  (use parameter blim)
!           simulations showed if b > 5.5 unfrozen water content is
!           non-realistically high at very low temperatures.

      bx = bexp
      if (bexp > blim)  bx = blim

!  --- ...  initializing iterations counter and iterative solution flag.

      nlog  = 0
      kcount= 0

!  --- ...  if temperature not significantly below freezing (t0), sh2o = smc

      if (tkelv > (tfreez-1.e-3)) then

        liqwat = smc

      else

        if (ck /= 0.0) then

!  --- ...  option 1: iterated solution for nonzero ck
!                     in koren et al, jgr, 1999, eqn 17

!  --- ...  initial guess for swl (frozen content)

          swl = smc - sh2o

!  --- ...  keep within bounds.

          swl = max( min( swl, smc-0.02 ), 0.0 )

!  --- ...  start of iterations

          do while ( (nlog < 10) .and. (kcount == 0) )
            nlog = nlog + 1

            df = alog( (psis*gs2/lsubf) * ( (1.0 + ck*swl)**2.0 )      &
              * (smcmax/(smc-swl))**bx ) - alog(-(tkelv-tfreez)/tkelv)

            denom = 2.0*ck/(1.0 + ck*swl) + bx/(smc - swl)
            swlk  = swl - df/denom

!  --- ...  bounds useful for mathematical solution.

            swlk = max( min( swlk, smc-0.02 ), 0.0 )

!  --- ...  mathematical solution bounds applied.

            dswl = abs(swlk - swl)
            swl = swlk

!  --- ...  if more than 10 iterations, use explicit method (ck=0 approx.)
!           when dswl less or eq. error, no more iterations required.

            if ( dswl <= error )  then
              kcount = kcount + 1
            endif
          enddo   !  end do_while_loop

!  --- ...  bounds applied within do-block are valid for physical solution.

          liqwat = smc - swl

        endif   ! end if_ck_block

!  --- ...  option 2: explicit solution for flerchinger eq. i.e. ck=0
!                     in koren et al., jgr, 1999, eqn 17
!           apply physical bounds to flerchinger solution

        if (kcount == 0) then
          fk = ( ( (lsubf/(gs2*(-psis)))                   & 
            * ((tkelv-tfreez)/tkelv) )**(-1/bx) ) * smcmax

          fk = max( fk, 0.02 )

          liqwat = min( fk, smc )
        endif

      endif   ! end if_tkelv_block
!
      return
!...................................
      end subroutine frh2o
!-----------------------------------
      end module sflx_snippet
