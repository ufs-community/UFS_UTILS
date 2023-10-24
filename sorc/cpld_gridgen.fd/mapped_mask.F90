!> @file
!! @brief Create the mapped ocean mask files
!! @author Denise.Worthen@noaa.gov
!!
!> This writes the mapped ocean mask on the FV3 tiles
!! @author Denise.Worthen@noaa.gov

module mapped_mask

  use gengrid_kinds, only : dbl_kind,int_kind,CL,CM,CS
  use grdvars,       only : ni,nj,npx
  use charstrings,   only : dirout,res,atmres,logmsg
  use netcdf

  implicit none

contains

  !> Use ESMF weights to map the ocean land mask to the FV3 tiles and write the mapped mask to 6 tile files
  !!
  !! @param[in]  src a SCRIP file containing the land mask for the ocean domain
  !! @param[in]  wgt a file containing the ESMF weights to regrid from the ocean domain to the FV3 tile domain
  !!
  !! @author Denise.Worthen@noaa.gov

  subroutine make_frac_land(src, wgt)

    character(len=*), intent(in) :: src, wgt

    ! local variables
    integer, parameter :: ntile = 6
    integer(int_kind) :: n_a, n_b, n_s

    integer(int_kind), allocatable, dimension(:) :: col, row
    real(dbl_kind), allocatable, dimension(:) :: S
    real(dbl_kind), allocatable, dimension(:) :: lat1d, lon1d

    integer(int_kind), allocatable, dimension(:) :: src_field
    real(dbl_kind), allocatable, dimension(:) :: dst_field

    real(dbl_kind), allocatable, dimension(:,:)   :: dst2d
    real(dbl_kind), allocatable, dimension(:,:)   :: lat2d,lon2d

    character(len=CS) :: ctile
    character(len=CL) :: fdst
    integer :: i,ii,jj,id,rc,ncid, dim2(2)
    integer :: istr,iend
    integer :: idimid,jdimid

    character(len=CM) :: vname
    !---------------------------------------------------------------------
    ! retrieve the weights
    !---------------------------------------------------------------------

    rc = nf90_open(trim(wgt), nf90_nowrite, ncid)
    rc = nf90_inq_dimid(ncid, 'n_s', id)
    rc = nf90_inquire_dimension(ncid, id, len=n_s)
    rc = nf90_inq_dimid(ncid, 'n_a', id)
    rc = nf90_inquire_dimension(ncid, id, len=n_a)
    rc = nf90_inq_dimid(ncid, 'n_b', id)
    rc = nf90_inquire_dimension(ncid, id, len=n_b)

    allocate(col(1:n_s))
    allocate(row(1:n_s))
    allocate(  S(1:n_s))

    allocate(lat1d(1:n_b))
    allocate(lon1d(1:n_b))

    rc = nf90_inq_varid(ncid, 'col', id)
    rc = nf90_get_var(ncid,     id, col)
    rc = nf90_inq_varid(ncid, 'row', id)
    rc = nf90_get_var(ncid,     id, row)
    rc = nf90_inq_varid(ncid,   'S', id)
    rc = nf90_get_var(ncid,      id,  S)

    ! 1d-tiled lat,lon
    rc = nf90_inq_varid(ncid, 'yc_b',     id)
    rc = nf90_get_var(ncid,       id,  lat1d)
    rc = nf90_inq_varid(ncid, 'xc_b',     id)
    rc = nf90_get_var(ncid,       id,  lon1d)
    rc = nf90_close(ncid)

    !---------------------------------------------------------------------
    ! retrieve 1-d land mask from the SCRIP file and map it
    !---------------------------------------------------------------------

    allocate(src_field(1:n_a))
    allocate(dst_field(1:n_b))

    rc = nf90_open(trim(src), nf90_nowrite, ncid)

    !1-d ocean mask (integer)
    rc = nf90_inq_varid(ncid, 'grid_imask', id)
    rc = nf90_get_var(ncid,     id,  src_field)
    rc = nf90_close(ncid)

    dst_field = 0.0
    do i = 1,n_s
       ii = row(i); jj = col(i)
       dst_field(ii) = dst_field(ii) + S(i)*real(src_field(jj),dbl_kind)
    enddo

    !---------------------------------------------------------------------
    !
    !---------------------------------------------------------------------

    allocate(dst2d(npx,npx))
    allocate(lon2d(npx,npx)); allocate(lat2d(npx,npx))

    do i = 0,ntile-1
       istr = i*npx*npx+1
       iend = istr+npx*npx-1
       !print *,i,istr,iend

       write(ctile,'(a5,i1)')'.tile',i+1
       fdst = trim(dirout)//'/'//trim(atmres)//'.mx'//trim(res)//trim(ctile)//'.nc'
       logmsg = 'creating mapped ocean mask file '//trim(fdst)
       print '(a)',trim(logmsg)

       dst2d(:,:) = reshape(dst_field(istr:iend), (/npx,npx/))
       lat2d(:,:) = reshape(    lat1d(istr:iend), (/npx,npx/))
       lon2d(:,:) = reshape(    lon1d(istr:iend), (/npx,npx/))

       rc = nf90_create(trim(fdst), nf90_64bit_offset, ncid)
       rc = nf90_def_dim(ncid, 'grid_xt', npx, idimid)
       rc = nf90_def_dim(ncid, 'grid_yt', npx, jdimid)

       dim2(:) =  (/idimid, jdimid/)
       vname = 'grid_xt'
       rc = nf90_def_var(ncid, vname, nf90_double, dim2, id)
       vname = 'grid_yt'
       rc = nf90_def_var(ncid, vname, nf90_double, dim2, id)
       vname = 'land_frac'
       rc = nf90_def_var(ncid, vname, nf90_double, dim2, id)
       rc = nf90_enddef(ncid)

       rc = nf90_inq_varid(ncid,    'grid_xt',      id)
       rc = nf90_put_var(ncid,             id,   lon2d)
       rc = nf90_inq_varid(ncid,    'grid_yt',      id)
       rc = nf90_put_var(ncid,             id,   lat2d)
       rc = nf90_inq_varid(ncid,  'land_frac',      id)
       rc = nf90_put_var(ncid,             id,   dst2d)
       rc = nf90_close(ncid)
    end do

    !---------------------------------------------------------------------
    ! clean up
    !---------------------------------------------------------------------

    deallocate(col, row, S, lat1d, lon1d, src_field, dst_field)
    deallocate(dst2d,lon2d,lat2d)

  end subroutine make_frac_land
end module mapped_mask
