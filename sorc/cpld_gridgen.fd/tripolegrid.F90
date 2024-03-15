!> @file
!! @brief Write the tripole grid file
!! @author Denise.Worthen@noaa.gov
!!
!> This module writes the main tripole grid file
!! @author Denise.Worthen@noaa.gov

module tripolegrid

  use gengrid_kinds, only: dbl_kind,int_kind,CM
  use grdvars,       only: ni,nj,nv,nverts,ncoord
  use grdvars,       only: lonCt,latCt,lonCt_vert,latCt_vert
  use grdvars,       only: lonCu,latCu,lonCu_vert,latCu_vert
  use grdvars,       only: lonCv,latCv,lonCv_vert,latCv_vert
  use grdvars,       only: lonBu,latBu,lonBu_vert,latBu_vert
  use grdvars,       only: wet4,areaCt,angleT,dp4,angle,angchk
  use charstrings,   only: logmsg,history
  use vartypedefs,   only: maxvars, fixvars, fixvars_typedefine
  use netcdf

  implicit none
  private

  public write_tripolegrid

contains
  !> Write the tripole grid file
  !!
  !! @param[in]  fname  the name of the tripole grid file to write
  !!
  !! @author Denise.Worthen@noaa.gov

  subroutine write_tripolegrid(fname)

    character(len=*), intent(in) :: fname

    ! local variables
    integer :: ii,id,rc, ncid, dim2(2),dim3(3)
    integer :: idimid,jdimid,kdimid

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

    rc = nf90_def_dim(ncid, 'ni', ni, idimid)
    rc = nf90_def_dim(ncid, 'nj', nj, jdimid)
    rc = nf90_def_dim(ncid, 'nv', nv, kdimid)

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

    rc = nf90_inq_varid(ncid,   'wet',        id)
    rc = nf90_put_var(ncid,        id, int(wet4))

    rc = nf90_inq_varid(ncid,  'area',      id)
    rc = nf90_put_var(ncid,        id,  areaCt)

    rc = nf90_inq_varid(ncid,'anglet',      id)
    rc = nf90_put_var(ncid,        id,  anglet)

    rc = nf90_inq_varid(ncid, 'angle',      id)
    rc = nf90_put_var(ncid,        id,   angle)

    rc = nf90_inq_varid(ncid,'angchk',      id)
    rc = nf90_put_var(ncid,        id,  angchk)

    rc = nf90_inq_varid(ncid, 'depth',      id)
    rc = nf90_put_var(ncid,        id,     dp4)

    rc = nf90_inq_varid(ncid,  'lonCt',     id)
    rc = nf90_put_var(ncid,        id,   lonCt)

    rc = nf90_inq_varid(ncid,  'latCt',     id)
    rc = nf90_put_var(ncid,        id,   latCt)

    rc = nf90_inq_varid(ncid, 'lonCv',      id)
    rc = nf90_put_var(ncid,        id,   lonCv)

    rc = nf90_inq_varid(ncid, 'latCv',      id)
    rc = nf90_put_var(ncid,        id,   latCv)

    rc = nf90_inq_varid(ncid, 'lonCu',      id)
    rc = nf90_put_var(ncid,        id,   lonCu)

    rc = nf90_inq_varid(ncid, 'latCu',      id)
    rc = nf90_put_var(ncid,        id,   latCu)

    rc = nf90_inq_varid(ncid, 'lonBu',      id)
    rc = nf90_put_var(ncid,        id,   lonBu)

    rc = nf90_inq_varid(ncid, 'latBu',      id)
    rc = nf90_put_var(ncid,        id,   latBu)

    ! vertices
    rc = nf90_inq_varid(ncid,  'lonCt_vert',     id)
    rc = nf90_put_var(ncid,         id,  lonCt_vert)

    rc = nf90_inq_varid(ncid,  'latCt_vert',     id)
    rc = nf90_put_var(ncid,         id,  latCt_vert)

    rc = nf90_inq_varid(ncid, 'lonCv_vert',      id)
    rc = nf90_put_var(ncid,        id,   lonCv_vert)

    rc = nf90_inq_varid(ncid, 'latCv_vert',      id)
    rc = nf90_put_var(ncid,        id,   latCv_vert)

    rc = nf90_inq_varid(ncid, 'lonCu_vert',      id)
    rc = nf90_put_var(ncid,        id,   lonCu_vert)

    rc = nf90_inq_varid(ncid, 'latCu_vert',      id)
    rc = nf90_put_var(ncid,        id,   latCu_vert)

    rc = nf90_inq_varid(ncid, 'lonBu_vert',      id)
    rc = nf90_put_var(ncid,        id,   lonBu_vert)

    rc = nf90_inq_varid(ncid, 'latBu_vert',      id)
    rc = nf90_put_var(ncid,        id,   latBu_vert)

    rc = nf90_close(ncid)

  end subroutine write_tripolegrid
end module tripolegrid
