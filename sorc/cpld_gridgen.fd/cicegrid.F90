!> @file
!! @brief Write the CICE6 grid file
!! @author Denise.Worthen@noaa.gov
!!
!> Write the CICE6 grid file
!! @author Denise.Worthen@noaa.gov
module cicegrid

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
  !! @param[in]  iind,jind   the grid domain bounds
  !! @param[in]  G           the domain grid type
  !!
  !! @author Denise.Worthen@noaa.gov
  subroutine write_cicegrid(fname,iind,jind,G)

    use grdvars, only : grid_type

    character(len=*), intent(in) :: fname
    integer         , intent(in) :: iind(:), jind(:)
    type(grid_type) , intent(in) :: G

    ! local variables
    integer :: ii,id,rc, ncid, dim2(2)
    integer :: idimid,jdimid
    integer :: ib,ie,jb,je
    integer :: idim, jdim

    character(len=2)  :: vtype
    character(len=CM) :: vname
    character(len=CM) :: vlong
    character(len=CM) :: vunit

    ib = iind(1) ; ie = iind(2)
    jb = jind(1) ; je = jind(2)
    idim = (ie - ib) + 1
    jdim = (je - jb) + 1

    !---------------------------------------------------------------------
    ! create the netcdf file
    !---------------------------------------------------------------------

    ! define the output variables and file name
    call cicevars_typedefine

    rc = nf90_create(fname, nf90_write, ncid)
    logmsg = '==> writing CICE grid to '//trim(fname)
    print '(a)', trim(logmsg)
    if(rc .ne. 0)print '(a)', 'nf90_create = '//trim(nf90_strerror(rc))

    rc = nf90_def_dim(ncid, 'ni', idim, idimid)
    rc = nf90_def_dim(ncid, 'nj', jdim, jdimid)

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

    rc = nf90_inq_varid(ncid,  'ulon', id)
    rc = nf90_put_var(ncid,        id, G%ulon(ib:ie,jb:je))

    rc = nf90_inq_varid(ncid,  'ulat', id)
    rc = nf90_put_var(ncid,        id, G%ulat(ib:ie,jb:je))

    rc = nf90_inq_varid(ncid,   'htn', id)
    rc = nf90_put_var(ncid,        id, G%htn(ib:ie,jb:je))

    rc = nf90_inq_varid(ncid,   'hte', id)
    rc = nf90_put_var(ncid,        id, G%hte(ib:ie,jb:je))

    rc = nf90_inq_varid(ncid,  'angle', id)
    rc = nf90_put_var(ncid,         id, G%angle(ib:ie,jb:je))

    rc = nf90_inq_varid(ncid,    'kmt', id)
    rc = nf90_put_var(ncid,         id, int(G%wet4(ib:ie,jb:je)))

    rc = nf90_close(ncid)

  end subroutine write_cicegrid
end module cicegrid
