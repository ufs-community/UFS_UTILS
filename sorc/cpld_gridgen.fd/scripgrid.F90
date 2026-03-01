!> @file
!! @brief Write a SCRIP format file
!! @author Denise.Worthen@noaa.gov
!!
!> This module writes a SCRIP format file
!! @author Denise.Worthen@noaa.gov

module scripgrid

  use gengrid_kinds, only: dbl_kind,int_kind,CM
  use gengrid_utils, only: reshape_staggers
  use grdvars,       only: nv
  use charstrings,   only: logmsg
  use vartypedefs,   only: maxvars, scripvars, scripvars_typedefine
  use netcdf

  implicit none
  private

  public write_staggers

contains
  !> Reshape center and corner grid points for a given stagger location and write a SCRIP file
  !! @param[in]  fname             the file name to write
  !! @param[in]  iind,jjind        the grid domain bounds
  !! @param[out] lon,lat           2D center global lon,lat for a given stagger
  !! @param[out] lonvert, latvert  3D corner (vertices) global lon and lat for a given stagger
  !! @param[in]  imask (optional)  the land mask values
  subroutine write_staggers(fname,iind,jind,lon,lat,lonvert,latvert,imask)

    character(len=*), intent(in)           :: fname
    integer         , intent(in)           :: iind(:), jind(:)
    real(dbl_kind)  , intent(in)           :: lon(:,:), lat(:,:), lonvert(:,:,:), latvert(:,:,:)
    integer(int_kind),intent(in), optional :: imask(:,:)

    integer(int_kind), allocatable, dimension(:,:) :: lmask
    real(dbl_kind),    allocatable, dimension(:)   :: cnlons, cnlats
    real(dbl_kind),    allocatable, dimension(:,:) :: crlons, crlats
    integer(int_kind), allocatable, dimension(:)   :: cnmask

    integer :: ib,ie,jb,je
    integer :: idim, jdim

    ib = iind(1) ; ie = iind(2)
    jb = jind(1) ; je = jind(2)
    idim = (ie - ib) + 1
    jdim = (je - jb) + 1

    allocate(lmask(idim,jdim), source = 1_int_kind)
    if(present(imask))then
       lmask = imask
    end if

    allocate(cnlons(idim*jdim), source = 0.0_dbl_kind)
    allocate(cnlats(idim*jdim), source = 0.0_dbl_kind)
    allocate(cnmask(idim*jdim), source = 1_int_kind)
    allocate(crlons(nv,idim*jdim), source = 0.0_dbl_kind)
    allocate(crlats(nv,idim*jdim), source = 0.0_dbl_kind)

    call reshape_staggers((/ib,ie/),(/jb,je/),lon,lat,lmask,lonvert,latvert,cnlons,cnlats,cnmask,crlons,crlats)
    logmsg = 'creating SCRIP file '//trim(fname)
    print '(a)',trim(logmsg)
    call write_scripgrid(trim(fname),idim,jdim,cnlons,cnlats,crlons,crlats,cnmask)

    deallocate(lmask, cnlons, cnlats, crlons, crlats, cnmask)

  end subroutine write_staggers
  !> Write a SCRIP grid file
  !!
  !! @param[in]  fname             the file name to write
  !! @param[in]  idim, jdim        the 2D dimensions
  !! @param[in]  cnlons, cnlats    1D center lons,lats
  !! @param[in]  crlons, crlats    2D corner lons/lats
  !! @param[in]  cnmask (optional) the land mask values
  !!
  !! @author Denise.Worthen@noaa.gov
  subroutine write_scripgrid(fname, idim, jdim, cnlons, cnlats, crlons, crlats, cnmask)
    character(len=*),  intent(in) :: fname
    integer(int_kind), intent(in) :: idim,jdim
    real(dbl_kind),    intent(in) :: cnlons(:), cnlats(:)
    real(dbl_kind),    intent(in) :: crlons(:,:), crlats(:,:)
    integer(int_kind), intent(in) :: cnmask(:)

    integer, parameter :: grid_rank = 2

    integer :: ii, id, rc, ncid, dim2(2), dim1(1)
    integer :: idimid, jdimid, kdimid

    integer, dimension(grid_rank) :: gdims
    character(len=2)  :: vtype
    character(len=CM) :: vname
    character(len=CM) :: vunit

    gdims(:) = (/idim, jdim/)

    !---------------------------------------------------------------------
    ! create the netcdf file
    !---------------------------------------------------------------------

    ! define the output variables and file name
    call scripvars_typedefine
    ! create the file
    ! 64_bit offset reqd for 008 grid
    ! produces b4b results for smaller grids
    rc = nf90_create(trim(fname), nf90_64bit_offset, ncid)
    !logmsg = '==> writing SCRIP grid to '//trim(fname)
    !print '(a)',trim(logmsg)
    if(rc .ne. 0)print '(a)', 'nf90_create = '//trim(nf90_strerror(rc))

    rc = nf90_def_dim(ncid, 'grid_size', idim*jdim, idimid)
    rc = nf90_def_dim(ncid, 'grid_corners',     nv, jdimid)
    rc = nf90_def_dim(ncid, 'grid_rank', grid_rank, kdimid)

    !grid_dims
    dim1(:) = (/kdimid/)
    rc = nf90_def_var(ncid, 'grid_dims', nf90_int, dim1, id)
    ! mask
    dim1(:) = (/idimid/)
    rc = nf90_def_var(ncid, 'grid_imask', nf90_int, dim1, id)
    rc = nf90_put_att(ncid, id,     'units',      'unitless')

    ! centers
    do ii = 1,2
       vname = trim(scripvars(ii)%var_name)
       vunit = trim(scripvars(ii)%unit_name)
       vtype = trim(scripvars(ii)%var_type)
       dim1(:) =  (/idimid/)
       if(vtype .eq. 'r8')rc = nf90_def_var(ncid, vname, nf90_double, dim1, id)
       if(vtype .eq. 'r4')rc = nf90_def_var(ncid, vname, nf90_float,  dim1, id)
       if(vtype .eq. 'i4')rc = nf90_def_var(ncid, vname, nf90_int,    dim1, id)
       rc = nf90_put_att(ncid, id,     'units', vunit)
    enddo

    ! corners
    do ii = 3,4
       vname = trim(scripvars(ii)%var_name)
       vunit = trim(scripvars(ii)%unit_name)
       vtype = trim(scripvars(ii)%var_type)
       dim2(:) =  (/jdimid,idimid/)
       if(vtype .eq. 'r8')rc = nf90_def_var(ncid, vname, nf90_double, dim2, id)
       if(vtype .eq. 'r4')rc = nf90_def_var(ncid, vname, nf90_float,  dim2, id)
       if(vtype .eq. 'i4')rc = nf90_def_var(ncid, vname, nf90_int,    dim2, id)
       rc = nf90_put_att(ncid, id,     'units', vunit)
    enddo
    rc = nf90_enddef(ncid)

    rc = nf90_inq_varid(ncid, 'grid_dims', id)
    rc = nf90_put_var(ncid,            id, gdims)
    rc = nf90_inq_varid(ncid, 'grid_imask', id)
    rc = nf90_put_var(ncid,            id, cnmask)

    rc = nf90_inq_varid(ncid, 'grid_center_lon', id)
    rc = nf90_put_var(ncid,                  id, cnlons)
    rc = nf90_inq_varid(ncid, 'grid_center_lat', id)
    rc = nf90_put_var(ncid,                  id, cnlats)

    rc = nf90_inq_varid(ncid, 'grid_corner_lon', id)
    rc = nf90_put_var(ncid,                  id, crlons)
    rc = nf90_inq_varid(ncid, 'grid_corner_lat', id)
    rc = nf90_put_var(ncid,                  id, crlats)

    rc = nf90_close(ncid)

  end subroutine write_scripgrid
end module scripgrid
