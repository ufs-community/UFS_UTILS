!> @file
!! @brief Perform calculations needed when using an ocean file
!! @author Denise.Worthen@noaa.gov
!!
!> This module calculates quantities needed only in the ocean case
!! @author Denise.Worthen@noaa.gov
module ocncalc_mod

  use netcdf
  use init_mod ,  only : nlevs, nxr, nyr
  use init_mod,   only : debug, logunit

  use arrays_mod, only : hmin, maskspval, b3d, nbilin3d, rgb3d, eta
  use utils_mod,  only : getfield

  implicit none

  private

  public calc_eta
  public vfill

contains
  !> Calculate interface heights locally
  !!
  !! @param[in]  fname  the file name containing needed varables
  !! @param[in]  dims   the dimensions of the source domain
  !! @param[in]  bathy  the bathymetry on the source domain
  !!
  !! @author Denise.Worthen@noaa.gov
  subroutine calc_eta(fname,dims,bathy)

    character(len=*), intent(in)  :: fname
    integer,          intent(in)  :: dims(:)
    real(kind=8),     intent(in)  :: bathy(:)

    ! local variables
    integer      :: i,k
    real(kind=8) :: denom
    real(kind=8), allocatable, dimension(:)   :: ssh,dilate
    real(kind=8), allocatable, dimension(:,:) :: h
    real(kind=8), allocatable, dimension(:,:) :: etmp
    character(len=20) :: subname = 'calc_eta'
    !----------------------------------------------------------------------------

    if (debug)write(logunit,'(a)')'enter '//trim(subname)

    allocate(ssh(dims(1)*dims(2))); ssh = 0.0
    allocate(dilate(dims(1)*dims(2))); dilate = 0.0
    allocate(h(dims(3),dims(1)*dims(2))); h = 0.0
    ! note nlevs+1 for local array
    allocate(etmp(dims(3)+1,dims(1)*dims(2))); etmp = 0.0

    call getfield(trim(fname), 'sfc', (/dims(1),dims(2)/), ssh)
    call getfield(trim(fname),   'h', (/dims(1),dims(2),dims(3)/), h)

    ! ref: MOM_interface_heights.F90, SR find_eta_3d.F90, Boussinesq
    etmp(dims(3)+1,:) = -bathy(:)
    do k=dims(3),1,-1
       etmp(k,:) = etmp(k+1,:) + h(k,:)
    enddo

    do i = 1,dims(1)*dims(2)
       denom = etmp(1,i) + bathy(i)
       if (denom .ne. 0.0) then
          dilate(i) = (ssh(i) + bathy(i)) / (etmp(1,i) + bathy(i))
       end if
    end do

    eta = 0.0
    do k = 1,dims(3)
       eta(k,:) = dilate(:)*(etmp(k,:) + bathy(:)) - bathy(:)
    end do

    if (debug)write(logunit,'(a)')'exit '//trim(subname)
  end subroutine calc_eta

  !>  Fill water column vertically on the destination grid
  !!
  !! @author Denise.Worthen@noaa.gov
  subroutine vfill()

    ! local variables
    integer :: n,i,k
    integer :: idx1, klast
    character(len=20) :: subname = 'vfill'
    !----------------------------------------------------------------------------

    if (debug)write(logunit,'(a)')'enter '//trim(subname)

    idx1 = 0
    do n = 1,nbilin3d
       if (trim(b3d(n)%var_name) .eq. 'h') idx1 = n
    end do

    do i = 1,nxr*nyr
       klast = nlevs
       do k = 1,nlevs
          if (rgb3d(idx1,k,i) .lt. maskspval)klast = k
       end do
       do n = 1,nbilin3d
          do k = klast+1,nlevs
             if (trim(b3d(n)%var_name) .eq. 'h') then
                rgb3d(n,k,i) = hmin
             else
                rgb3d(n,k,i) = rgb3d(n,klast,i)
             end if
          end do
       end do
    end do

    if (debug)write(logunit,'(a)')'exit '//trim(subname)
  end subroutine vfill
end module ocncalc_mod
