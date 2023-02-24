!> @file
!! @brief Print debugging information
!! @author Denise.Worthen@noaa.gov
!!
!> Print values for debugging
!! @author Denise.Worthen@noaa.gov

module debugprint

  use esmf,    only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
  use grdvars, only : ni,nj,ipole,angle,angleT
  use grdvars, only : htn,latCt,lonCt,latCv,lonCv,latCu,lonCu,latBu,lonBu
  use grdvars, only : xlatCt,xlonCt,xlatCu,xlonCu
  use grdvars, only : latBu_vert,lonBu_vert,latCv_vert,lonCv_vert
  use grdvars, only : latCt_vert,lonCt_vert,latCu_vert,lonCu_vert

  implicit none
  private

  public :: checkseam, checkxlatlon, checkpoint

contains
  !> Print values across the tripole seam
  !!
  !! @author Denise.Worthen@noaa.gov

  subroutine checkseam

    ! local variables
    integer :: j,i1,i2

    j = nj
    i1 = ipole(1); i2 = ipole(2)+1

    !htn must be the same along seam
    j = nj
    i1 = ipole(1); i2 = ipole(2)+1
    print *,'HTN across seam '
    print *,htn(i1-2,j),htn(i2+2,j)
    print *,htn(i1-1,j),htn(i2+1,j)
    print *,htn(i1,  j),htn(i2,  j)
    print *,htn(i1+1,j),htn(i2-1,j)
    print *,htn(i1+2,j),htn(i2-2,j)

    print *,'latCv across seam '
    print *,latCv(i1-2,j),latCv(i2+2,j)
    print *,latCv(i1-1,j),latCv(i2+1,j)
    print *,latCv(i1,  j),latCv(i2,  j)
    print *,latCv(i1+1,j),latCv(i2-1,j)
    print *,latCv(i1+2,j),latCv(i2-2,j)

    print *,'lonCv across seam '
    print *,lonCv(i1-2,j),lonCv(i2+2,j)
    print *,lonCv(i1-1,j),lonCv(i2+1,j)
    print *,lonCv(i1,  j),lonCv(i2,  j)
    print *,lonCv(i1+1,j),lonCv(i2-1,j)
    print *,lonCv(i1+2,j),lonCv(i2-2,j)

    print *,'angleT across seam '
    print *,angleT(i1-2,j),angleT(i2+2,j)
    print *,angleT(i1-1,j),angleT(i2+1,j)
    print *,angleT(i1,  j),angleT(i2,  j)
    print *,angleT(i1+1,j),angleT(i2-1,j)
    print *,angleT(i1+2,j),angleT(i2-2,j)

    print *,'latCu across seam '
    print *,latCu(i1-3,j),latCu(i2+2,j),latCu(i1-3,j)-latCu(i2+2,j)
    print *,latCu(i1-2,j),latCu(i2+1,j)
    print *,latCu(i1-1,j),latCu(i2+0,j)
    print *,latCu(i1,  j),latCu(i2-1,j)
    print *,latCu(i1+1,j),latCu(i2-2,j)
    print *,latCu(i1+2,j),latCu(i2-3,j)
    print *,latCu(i1+3,j),latCu(i2-4,j)

    print *,'lonCu across seam '
    print *,lonCu(i1-3,j),lonCu(i2+2,j),lonCu(i1-3,j)+lonCu(i2+2,j)
    print *,lonCu(i1-2,j),lonCu(i2+1,j)
    print *,lonCu(i1-1,j),lonCu(i2+0,j)
    print *,lonCu(i1,  j),lonCu(i2-1,j)
    print *,lonCu(i1+1,j),lonCu(i2-2,j)
    print *,lonCu(i1+2,j),lonCu(i2-3,j)
    print *,lonCu(i1+3,j),lonCu(i2-4,j)

    print *,'latCt across seam '
    print *,latCt(i1-3,j),latCt(i2+3,j),latCt(i1-3,j)-latCt(i2+3,j)
    print *,latCt(i1-2,j),latCt(i2+2,j)
    print *,latCt(i1-1,j),latCt(i2+1,j)
    print *,latCt(i1,  j),latCt(i2,  j)
    print *,latCt(i1+1,j),latCt(i2-1,j)
    print *,latCt(i1+2,j),latCt(i2-2,j)
    print *,latCt(i1+3,j),latCt(i2-3,j)

    print *,'lonCt across seam '
    print *,lonCt(i1-3,j),lonCt(i2+3,j),lonCt(i1-3,j)+lonCt(i2+3,j)
    print *,lonCt(i1-2,j),lonCt(i2+2,j)
    print *,lonCt(i1-1,j),lonCt(i2+1,j)
    print *,lonCt(i1,  j),lonCt(i2,  j)
    print *,lonCt(i1+1,j),lonCt(i2-1,j)
    print *,lonCt(i1+2,j),lonCt(i2-2,j)
    print *,lonCt(i1+3,j),lonCt(i2-3,j)
    print *
  end subroutine checkseam

  !> Print values near the poles and along the domain edges
  !!
  !! @author Denise.Worthen@noaa.gov

  subroutine checkxlatlon

    ! local variables
    integer :: i

    print *,'============== Ct grid ==============='
    print *,'============== Left pole ============'
    do i = ipole(1)-3,ipole(1)+3
       print '(i5,6f12.5)',i,lonCt(i,nj),xlonCt(i),lonCt(i,nj)+xlonCt(i),latCt(i,nj),xlatCt(i),latCt(i,nj)-xlatCt(i)
    enddo
    print *

    print *,'============ Right pole ============'
    do i = ipole(2)-3,ipole(2)+3
       print '(i5,6f12.5)',i,lonCt(i,nj),xlonCt(i),lonCt(i,nj)+xlonCt(i),latCt(i,nj),xlatCt(i),latCt(i,nj)-xlatCt(i)
    enddo
    print *

    print *,'============== Ct grid ==============='
    print *,'============== Left edge ============'
    do i = 1,5
       print '(i5,6f12.5)',i,lonCt(i,nj),xlonCt(i),lonCt(i,nj)+xlonCt(i),latCt(i,nj),xlatCt(i),latCt(i,nj)-xlatCt(i)
    enddo
    print *
    print *,'============== Right edge ==========='
    do i = ni-4,ni
       print '(i5,6f12.5)',i,lonCt(i,nj),xlonCt(i),lonCt(i,nj)+xlonCt(i),latCt(i,nj),xlatCt(i),latCt(i,nj)-xlatCt(i)
    enddo
    print *


    print *,'============== Cu grid ==============='
    print *,'============== Left pole ============='
    do i = ipole(1)-3,ipole(1)+3
       print '(i5,6f12.5)',i,lonCu(i,nj),xlonCu(i),lonCu(i,nj)+xlonCu(i),latCu(i,nj),xlatCu(i),latCu(i,nj)-xlatCu(i)
    enddo
    print *

    print *,'============ Right pole ============'
    do i = ipole(2)-3,ipole(2)+3
       print '(i5,6f12.5)',i,lonCu(i,nj),xlonCu(i),lonCu(i,nj)+xlonCu(i),latCu(i,nj),xlatCu(i),latCu(i,nj)-xlatCu(i)
    enddo
    print *

    print *,'============== Cu grid ==============='
    print *,'============== Left edge ============'
    do i = 1,5
       print '(i5,6f12.5)',i,lonCu(i,nj),xlonCu(i),lonCu(i,nj)+xlonCu(i),latCu(i,nj),xlatCu(i),latCu(i,nj)-xlatCu(i)
    enddo
    print *
    print *,'============== Right edge ==========='
    do i = ni-4,ni
       print '(i5,6f12.5)',i,lonCu(i,nj),xlonCu(i),lonCu(i,nj)+xlonCu(i),latCu(i,nj),xlatCu(i),latCu(i,nj)-xlatCu(i)
    enddo
    print *

  end subroutine checkxlatlon

  !> Print values at specified point
  !!
  !! @author Denise.Worthen@noaa.gov

  subroutine checkpoint

    ! local variables
    integer :: i,j

    ! check
    i = 1; j = nj
    print '(f12.5,a,f12.5)',latBu_vert(i,j,2),'        ',latBu_vert(i,j,1)
    print '(a12,f12.5)','          ',latBu(i,j)
    print '(f12.5,a,f12.5)',latBu_vert(i,j,3),'        ',latBu_vert(i,j,4)
    print *
    print '(f12.5,a,f12.5)',lonBu_vert(i,j,2),'        ',lonBu_vert(i,j,1)
    print '(a12,f12.5)','          ',lonBu(i,j)
    print '(f12.5,a,f12.5)',lonBu_vert(i,j,3),'        ',lonBu_vert(i,j,4)
    print *
    print *
    ! check
    print '(f12.5,a,f12.5)',latCv_vert(i,j,2),'        ',latCv_vert(i,j,1)
    print '(a12,f12.5)','          ',latCv(i,j)
    print '(f12.5,a,f12.5)',latCv_vert(i,j,3),'        ',latCv_vert(i,j,4)
    print *
    print '(f12.5,a,f12.5)',lonCv_vert(i,j,2),'        ',lonCv_vert(i,j,1)
    print '(a12,f12.5)','          ',lonCv(i,j)
    print '(f12.5,a,f12.5)',lonCv_vert(i,j,3),'        ',lonCv_vert(i,j,4)

    print *
    print *

    i = 1; j = 10
    print '(f12.5,a,f12.5)',latCt_vert(i,j,2),'        ',latCt_vert(i,j,1)
    print '(a12,f12.5)','          ',latCt(i,j)
    print '(f12.5,a,f12.5)',latCt_vert(i,j,3),'        ',latCt_vert(i,j,4)
    print *
    print '(f12.5,a,f12.5)',lonCt_vert(i,j,2),'        ',lonCt_vert(i,j,1)
    print '(a12,f12.5)','          ',lonCt(i,j)
    print '(f12.5,a,f12.5)',lonCt_vert(i,j,3),'        ',lonCt_vert(i,j,4)
    print *
    print *
    ! check
    print '(f12.5,a,f12.5)',latCu_vert(i,j,2),'        ',latCu_vert(i,j,1)
    print '(a12,f12.5)','          ',latCu(i,j)
    print '(f12.5,a,f12.5)',latCu_vert(i,j,3),'        ',latCu_vert(i,j,4)
    print *
    print '(f12.5,a,f12.5)',lonCu_vert(i,j,2),'        ',lonCu_vert(i,j,1)
    print '(a12,f12.5)','          ',lonCu(i,j)
    print '(f12.5,a,f12.5)',lonCu_vert(i,j,3),'        ',lonCu_vert(i,j,4)


    i = ni; j = 10
    print '(f12.5,a,f12.5)',latCt_vert(i,j,2),'        ',latCt_vert(i,j,1)
    print '(a12,f12.5)','          ',latCt(i,j)
    print '(f12.5,a,f12.5)',latCt_vert(i,j,3),'        ',latCt_vert(i,j,4)
    print *
    print '(f12.5,a,f12.5)',lonCt_vert(i,j,2),'        ',lonCt_vert(i,j,1)
    print '(a12,f12.5)','          ',lonCt(i,j)
    print '(f12.5,a,f12.5)',lonCt_vert(i,j,3),'        ',lonCt_vert(i,j,4)
    print *
    print *
    ! check
    print '(f12.5,a,f12.5)',latCu_vert(i,j,2),'        ',latCu_vert(i,j,1)
    print '(a12,f12.5)','          ',latCu(i,j)
    print '(f12.5,a,f12.5)',latCu_vert(i,j,3),'        ',latCu_vert(i,j,4)
    print *
    print '(f12.5,a,f12.5)',lonCu_vert(i,j,2),'        ',lonCu_vert(i,j,1)
    print '(a12,f12.5)','          ',lonCu(i,j)
    print '(f12.5,a,f12.5)',lonCu_vert(i,j,3),'        ',lonCu_vert(i,j,4)

    print *,"latCt minmax ",minval(latCt),maxval(latCt)
    print *,"latCu minmax ",minval(latCu),maxval(latCu)
    print *,"latCv minmax ",minval(latCv),maxval(latCv)
    print *,"latBu minmax ",minval(latBu),maxval(latBu)

    !   print *,minval(latCt_vert),maxval(latCt_vert)
    !   print *,minval(lonCt_vert),maxval(lonCt_vert)
    !   print *,minval(latBu_vert),maxval(latBu_vert)
    !   print *,minval(lonBu_vert),maxval(lonBu_vert)
  end subroutine checkpoint
end module debugprint
