!> @file
!! @brief Write the CICE6 grid file
!! @author Denise.Worthen@noaa.gov
!!
!> Write the CICE6 grid file
!! @author Denise.Worthen@noaa.gov

module cicegrid

  use grdvars,       only: ni,nj,ulat,ulon,htn,hte,angle,wet4
  use charstrings,   only: history, logmsg
  use vartypedefs,   only: maxvars, cicevars, cicevars_typedefine
  use gengrid_kinds, only: CM
  use netcdf

  implicit none
  private

  public write_cicegrid

contains
  !> Write the CICE6 grid file
  !!
  !! @param[in]  fname  the name of the CICE6 grid file to write
  !!
  !! @author Denise.Worthen@noaa.gov

  subroutine write_cicegrid(fname)

    character(len=*), intent(in) :: fname

    ! local variables
    integer :: ii,id,rc, ncid, dim2(2)
    integer :: idimid,jdimid

    character(len=2)  :: vtype
    character(len=CM) :: vname
    character(len=CM) :: vlong
    character(len=CM) :: vunit

    !---------------------------------------------------------------------
    ! create the netcdf file
    !---------------------------------------------------------------------

    ! define the output variables and file name
    call cicevars_typedefine

    rc = nf90_create(fname, nf90_write, ncid)
    logmsg = '==> writing CICE grid to '//trim(fname)
    print '(a)', trim(logmsg)
    if(rc .ne. 0)print '(a)', 'nf90_create = '//trim(nf90_strerror(rc))

    rc = nf90_def_dim(ncid, 'ni', ni, idimid)
    rc = nf90_def_dim(ncid, 'nj', nj, jdimid)

    do ii = 1,maxvars
       if(len_trim(cicevars(ii)%var_name) .gt. 0)then
          vname = trim(cicevars(ii)%var_name)
          vlong = trim(cicevars(ii)%long_name)
          vunit = trim(cicevars(ii)%unit_name)
          vtype = trim(cicevars(ii)%var_type)

          dim2(:) =  (/idimid, jdimid/)
          if(vtype .eq. 'r8')rc = nf90_def_var(ncid, vname, nf90_double, dim2, id)
          if(vtype .eq. 'r4')rc = nf90_def_var(ncid, vname, nf90_float,  dim2, id)
          if(vtype .eq. 'i4')rc = nf90_def_var(ncid, vname, nf90_int,    dim2, id)
          rc = nf90_put_att(ncid, id,     'units', vunit)
          rc = nf90_put_att(ncid, id, 'long_name', vlong)
       end if
    enddo
    rc = nf90_put_att(ncid, nf90_global, 'history', trim(history))
    rc = nf90_enddef(ncid)

    rc = nf90_inq_varid(ncid,  'ulon',      id)
    rc = nf90_put_var(ncid,        id,    ulon)

    rc = nf90_inq_varid(ncid,  'ulat',      id)
    rc = nf90_put_var(ncid,        id,    ulat)

    rc = nf90_inq_varid(ncid,   'htn',      id)
    rc = nf90_put_var(ncid,        id,     htn)

    rc = nf90_inq_varid(ncid,   'hte',      id)
    rc = nf90_put_var(ncid,        id,     hte)

    rc = nf90_inq_varid(ncid,  'angle',     id)
    rc = nf90_put_var(ncid,         id,  angle)

    rc = nf90_inq_varid(ncid,    'kmt',        id)
    rc = nf90_put_var(ncid,         id, int(wet4))

    rc = nf90_close(ncid)

  end subroutine write_cicegrid
end module cicegrid
