!> @file
!! @brief Write the tripole grid file
!! @author Denise.Worthen@noaa.gov
!!
!> This module writes the main tripole grid file
!! @author Denise.Worthen@noaa.gov
module tripolegrid

  use gengrid_kinds, only: dbl_kind,int_kind,CM
  use charstrings,   only: logmsg,history
  use vartypedefs,   only: maxvars, fixvars, fixvars_typedefine
  use netcdf

  implicit none
  private

  public write_tripolegrid

contains
  !> Write the tripole grid file
  !!
  !! @param[in]  fname       the name of the tripole grid file to write
  !! @param[in]  iind,jind   the grid domain bounds
  !! @param[in]  G           the domain grid type
  !! @author Denise.Worthen@noaa.gov
  subroutine write_tripolegrid(fname, iind, jind, G)

    use grdvars, only : grid_type, nv, ncoord, nverts

    character(len=*), intent(in) :: fname
    integer         , intent(in) :: iind(:), jind(:)
    type(grid_type) , intent(in) :: G

    ! local variables
    integer :: ii,id,rc, ncid, dim2(2),dim3(3)
    integer :: idimid,jdimid,kdimid
    integer :: ib,ie,jb,je
    integer :: idim, jdim

    ib = iind(1) ; ie = iind(2)
    jb = jind(1) ; je = jind(2)
    idim = (ie - ib) + 1
    jdim = (je - jb) + 1

    !---------------------------------------------------------------------
    ! create the netcdf file
    !---------------------------------------------------------------------

    ! define the output variables and file name
    call fixvars_typedefine

    ! create the file
    ! 64_bit offset reqd for 008 grid
    ! produces b4b results for smaller grids
    rc = nf90_create(trim(fname), nf90_64bit_offset, ncid)
    logmsg = '==> writing tripole grid to '//trim(fname)
    print '(a)', trim(logmsg)
    if(rc .ne. 0)print '(a)', 'nf90_create = '//trim(nf90_strerror(rc))

    rc = nf90_def_dim(ncid, 'ni', idim, idimid)
    rc = nf90_def_dim(ncid, 'nj', jdim, jdimid)
    rc = nf90_def_dim(ncid, 'nv',   nv, kdimid)

    !mask
    dim2(:) = (/idimid, jdimid/)
    rc = nf90_def_var(ncid, 'wet',     nf90_int,   dim2, id)
    rc = nf90_put_att(ncid, id,     'units',           'nd')
    !area
    rc = nf90_def_var(ncid, 'area', nf90_double,   dim2, id)
    rc = nf90_put_att(ncid, id,     'units',           'm2')
    !angleT
    rc = nf90_def_var(ncid, 'anglet', nf90_double, dim2, id)
    rc = nf90_put_att(ncid, id,     'units',      'radians')
    !angle (angBu)
    rc = nf90_def_var(ncid,  'angle', nf90_double, dim2, id)
    rc = nf90_put_att(ncid, id,     'units',      'radians')
    !angchk
    rc = nf90_def_var(ncid, 'angchk', nf90_double, dim2, id)
    rc = nf90_put_att(ncid, id,     'units',      'radians')
    !bathymetry
    rc = nf90_def_var(ncid,  'depth', nf90_float,  dim2, id)
    rc = nf90_put_att(ncid, id,     'units',            'm')

    dim2(:) = (/idimid, jdimid/)
    do ii = 1,ncoord
       rc = nf90_def_var(ncid, trim(fixvars(ii)%var_name), nf90_double, dim2, id)
       rc = nf90_put_att(ncid, id,     'units', trim(fixvars(ii)%unit_name))
       rc = nf90_put_att(ncid, id, 'long_name', trim(fixvars(ii)%long_name))
       if(trim(fixvars(ii)%var_name(1:3)) .eq. "lon")then
          rc = nf90_put_att(ncid, id,  'lon_bounds', trim(fixvars(ii)%vertices))
       else
          rc = nf90_put_att(ncid, id,  'lat_bounds', trim(fixvars(ii)%vertices))
       endif
    enddo

    dim3(:) = (/idimid, jdimid, kdimid/)
    do ii = ncoord+1,ncoord+nverts
       rc = nf90_def_var(ncid, trim(fixvars(ii)%var_name), nf90_double, dim3, id)
       rc = nf90_put_att(ncid, id,     'units', trim(fixvars(ii)%unit_name))
       rc = nf90_put_att(ncid, id, 'long_name', trim(fixvars(ii)%long_name))
    enddo

    rc = nf90_put_att(ncid, nf90_global, 'history', trim(history))
    rc = nf90_enddef(ncid)

    rc = nf90_inq_varid(ncid,   'wet', id)
    rc = nf90_put_var(ncid,        id, int(G%wet4(ib:ie,jb:je)))

    rc = nf90_inq_varid(ncid,  'area', id)
    rc = nf90_put_var(ncid,        id, G%areaCt(ib:ie,jb:je))

    rc = nf90_inq_varid(ncid,'anglet', id)
    rc = nf90_put_var(ncid,        id, G%anglet(ib:ie,jb:je))

    rc = nf90_inq_varid(ncid, 'angle', id)
    rc = nf90_put_var(ncid,        id, G%angle(ib:ie,jb:je))

    rc = nf90_inq_varid(ncid,'angchk', id)
    rc = nf90_put_var(ncid,        id, G%angchk(ib:ie,jb:je))

    rc = nf90_inq_varid(ncid, 'depth', id)
    rc = nf90_put_var(ncid,        id, G%dp4(ib:ie,jb:je))

    rc = nf90_inq_varid(ncid, 'lonCt', id)
    rc = nf90_put_var(ncid,        id, G%Ct%lon(ib:ie,jb:je))

    rc = nf90_inq_varid(ncid, 'latCt', id)
    rc = nf90_put_var(ncid,        id, G%Ct%lat(ib:ie,jb:je))

    rc = nf90_inq_varid(ncid, 'lonCv', id)
    rc = nf90_put_var(ncid,        id, G%Cv%lon(ib:ie,jb:je))

    rc = nf90_inq_varid(ncid, 'latCv', id)
    rc = nf90_put_var(ncid,        id, G%Cv%lat(ib:ie,jb:je))

    rc = nf90_inq_varid(ncid, 'lonCu', id)
    rc = nf90_put_var(ncid,        id, G%Cu%lon(ib:ie,jb:je))

    rc = nf90_inq_varid(ncid, 'latCu', id)
    rc = nf90_put_var(ncid,        id, G%Cu%lat(ib:ie,jb:je))

    rc = nf90_inq_varid(ncid, 'lonBu', id)
    rc = nf90_put_var(ncid,        id, G%Bu%lon(ib:ie,jb:je))

    rc = nf90_inq_varid(ncid, 'latBu', id)
    rc = nf90_put_var(ncid,        id, G%Bu%lat(ib:ie,jb:je))

    ! vertices
    rc = nf90_inq_varid(ncid, 'lonCt_vert', id)
    rc = nf90_put_var(ncid,             id, G%Ct%lonvert(ib:ie,jb:je,:))

    rc = nf90_inq_varid(ncid, 'latCt_vert', id)
    rc = nf90_put_var(ncid,             id, G%Ct%latvert(ib:ie,jb:je,:))

    rc = nf90_inq_varid(ncid, 'lonCv_vert', id)
    rc = nf90_put_var(ncid,             id, G%Cv%lonvert(ib:ie,jb:je,:))

    rc = nf90_inq_varid(ncid, 'latCv_vert', id)
    rc = nf90_put_var(ncid,             id, G%Cv%latvert(ib:ie,jb:je,:))

    rc = nf90_inq_varid(ncid, 'lonCu_vert', id)
    rc = nf90_put_var(ncid,             id, G%Cu%lonvert(ib:ie,jb:je,:))

    rc = nf90_inq_varid(ncid, 'latCu_vert', id)
    rc = nf90_put_var(ncid,             id, G%Cu%latvert(ib:ie,jb:je,:))

    rc = nf90_inq_varid(ncid, 'lonBu_vert', id)
    rc = nf90_put_var(ncid,             id, G%Bu%lonvert(ib:ie,jb:je,:))

    rc = nf90_inq_varid(ncid, 'latBu_vert', id)
    rc = nf90_put_var(ncid,             id, G%Bu%latvert(ib:ie,jb:je,:))

    rc = nf90_close(ncid)

  end subroutine write_tripolegrid
end module tripolegrid
