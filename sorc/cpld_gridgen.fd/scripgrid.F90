!> @file
!! @brief Write a SCRIP format file
!! @author Denise.Worthen@noaa.gov
!!
!> This module writes a SCRIP format file
!! @author Denise.Worthen@noaa.gov

module scripgrid

  use gengrid_kinds, only: dbl_kind,int_kind,CM
  use grdvars,       only: ni,nj,nv
  use grdvars,       only: lonCt,latCt,lonCt_vert,latCt_vert
  use grdvars,       only: lonCu,latCu,lonCu_vert,latCu_vert
  use grdvars,       only: lonCv,latCv,lonCv_vert,latCv_vert
  use grdvars,       only: lonBu,latBu,lonBu_vert,latBu_vert
  use charstrings,   only: logmsg
  use vartypedefs,   only: maxvars, scripvars, scripvars_typedefine
  use netcdf

  implicit none
  private

  public write_scripgrid

contains
  !> Write a SCRIP grid file
  !!
  !! @param[in]  fname  the file name to write
  !! @param[in]  cstagger  the name of the stagger location
  !! @param[in]  imask (optional)  the land mask values
  !!
  !! @author Denise.Worthen@noaa.gov

  subroutine write_scripgrid(fname,cstagger, imask)

    character(len=*) , intent(in) :: fname
    character(len=*) , intent(in) :: cstagger
    integer(int_kind), optional, intent(in) :: imask(:,:)

    ! local variables
    integer, parameter :: grid_rank = 2

    integer :: ii,n,id,rc, ncid, dim2(2),dim1(1)
    integer :: idimid,jdimid,kdimid

    integer, dimension(grid_rank) :: gdims
    integer(int_kind), dimension(ni*nj)    :: cnmask          !1-d mask
    real(dbl_kind),    dimension(ni*nj)    :: cnlons, cnlats  !1-d center lats,lons
    real(dbl_kind),    dimension(nv,ni*nj) :: crlons, crlats  !2-d corner lats,lons

    real(dbl_kind), dimension(ni,nj) :: tmp

    character(len=2)  :: vtype
    character(len=CM) :: vname
    character(len=CM) :: vunit

    !---------------------------------------------------------------------
    !
    !---------------------------------------------------------------------

    gdims(:) = (/ni,nj/)
    if(trim(cstagger) .eq. 'Ct')then
       cnlons = reshape(lonCt, (/ni*nj/))
       cnlats = reshape(latCt, (/ni*nj/))
       do n = 1,nv
          tmp(:,:) = lonCt_vert(:,:,n)
          crlons(n,:) = reshape(tmp, (/ni*nj/))
          tmp(:,:) = latCt_vert(:,:,n)
          crlats(n,:) = reshape(tmp, (/ni*nj/))
       end do
    end if

    if(trim(cstagger) .eq. 'Cu')then
       cnlons = reshape(lonCu, (/ni*nj/))
       cnlats = reshape(latCu, (/ni*nj/))
       do n = 1,nv
          tmp(:,:) = lonCu_vert(:,:,n)
          crlons(n,:) = reshape(tmp, (/ni*nj/))
          tmp(:,:) = latCu_vert(:,:,n)
          crlats(n,:) = reshape(tmp, (/ni*nj/))
       end do
    end if

    if(trim(cstagger) .eq. 'Cv')then
       cnlons = reshape(lonCv, (/ni*nj/))
       cnlats = reshape(latCv, (/ni*nj/))
       do n = 1,nv
          tmp(:,:) = lonCv_vert(:,:,n)
          crlons(n,:) = reshape(tmp, (/ni*nj/))
          tmp(:,:) = latCv_vert(:,:,n)
          crlats(n,:) = reshape(tmp, (/ni*nj/))
       end do
    end if

    if(trim(cstagger) .eq. 'Bu')then
       cnlons = reshape(lonBu, (/ni*nj/))
       cnlats = reshape(latBu, (/ni*nj/))
       do n = 1,nv
          tmp(:,:) = lonBu_vert(:,:,n)
          crlons(n,:) = reshape(tmp, (/ni*nj/))
          tmp(:,:) = latBu_vert(:,:,n)
          crlats(n,:) = reshape(tmp, (/ni*nj/))
       end do
    end if

    if(present(imask))then
       cnmask = reshape(imask, (/ni*nj/))
    else
       cnmask = 1
    end if

    !---------------------------------------------------------------------
    ! create the netcdf file
    !---------------------------------------------------------------------

    ! define the output variables and file name
    call scripvars_typedefine
    ! create the file
    ! 64_bit offset reqd for 008 grid
    ! produces b4b results for smaller grids
    rc = nf90_create(trim(fname), nf90_64bit_offset, ncid)
    logmsg = '==> writing SCRIP grid to '//trim(fname)
    print '(a)',trim(logmsg)
    if(rc .ne. 0)print '(a)', 'nf90_create = '//trim(nf90_strerror(rc))

    rc = nf90_def_dim(ncid, 'grid_size',     ni*nj, idimid)
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

    rc = nf90_inq_varid(ncid,  'grid_dims',        id)
    rc = nf90_put_var(ncid,             id,     gdims)
    rc = nf90_inq_varid(ncid, 'grid_imask',        id)
    rc = nf90_put_var(ncid,             id,    cnmask)

    rc = nf90_inq_varid(ncid,  'grid_center_lon',        id)
    rc = nf90_put_var(ncid,                   id,    cnlons)
    rc = nf90_inq_varid(ncid,  'grid_center_lat',        id)
    rc = nf90_put_var(ncid,                   id,    cnlats)

    rc = nf90_inq_varid(ncid,  'grid_corner_lon',        id)
    rc = nf90_put_var(ncid,                   id,    crlons)
    rc = nf90_inq_varid(ncid,  'grid_corner_lat',        id)
    rc = nf90_put_var(ncid,                   id,    crlats)

    rc = nf90_close(ncid)

  end subroutine write_scripgrid
end module scripgrid
