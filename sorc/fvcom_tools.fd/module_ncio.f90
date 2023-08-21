!> @file
!! @brief Functions to read and write netcdf files.
!! @author Ming Hu @date 2017-11-01

!> Functions to read and write netcdf files.
!! @author Ming Hu @date 2017-11-01
module module_ncio

  use netcdf
  implicit none

  public :: ncio
  ! set default to private
  private
  !
  type :: ncio
     character(len=256) :: filename !< Name of data file.
     integer :: ncid !< File ID.
     integer :: status !< Return code.
     integer :: debug_level !< Debug level.

     integer :: nDims !< Number of dims.
     integer :: ends(4) !< Counts of dims.
     integer :: xtype !< Type of data.
     character(len=40) :: dimname(4) !< Name of dims.
   contains
     procedure :: open => open_nc !< Open netCDF file. @return
     procedure :: close => close_nc !< Close netCDF file. @return
     procedure :: get_dim => get_dim_nc !< read in dimension from the nc file @return
     generic   :: get_att => get_att_nc_int,get_att_nc_real,get_att_nc_string !< Get attribute. @return
     procedure :: get_att_nc_int !< Get attribute. @return
     procedure :: get_att_nc_real !< Get attribute. @return
     procedure :: get_att_nc_string !< Get attribute. @return
     generic   :: get_var => get_var_nc_double_1d, get_var_nc_double_2d, & 
          get_var_nc_double_3d,                       &
          get_var_nc_real_1d,get_var_nc_real_2d,      &
          get_var_nc_real_3d,                         &
          get_var_nc_short_1d,get_var_nc_short_2d,    &
          get_var_nc_int_1d,get_var_nc_int_2d,        &
          get_var_nc_int_3d,                          &
          get_var_nc_char_1d,get_var_nc_char_2d,      &
          get_var_nc_char_3d      !< Read in a 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: get_var_nc_short !< Read in a 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: get_var_nc_short_1d !< Read in a 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: get_var_nc_short_2d !< Read in a 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: get_var_nc_int !< Read in a 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: get_var_nc_int_1d !< Read in a 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: get_var_nc_int_2d !< Read in a 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: get_var_nc_int_3d !< Read in a 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: get_var_nc_real !< Read in a 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: get_var_nc_real_1d !< Read in a 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: get_var_nc_real_2d !< Read in a 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: get_var_nc_real_3d !< Read in a 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: get_var_nc_double !< Read in a 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: get_var_nc_double_1d !< Read in a 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: get_var_nc_double_2d !< Read in a 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: get_var_nc_double_3d !< Read in a 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: get_var_nc_char !< Read in a 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: get_var_nc_char_1d !< Read in a 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: get_var_nc_char_2d !< Read in a 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: get_var_nc_char_3d !< Read in a 1d, 2d, 3d, or 4d field from the nc file. @return
     generic   :: replace_var => replace_var_nc_double_1d, replace_var_nc_double_2d, &
          replace_var_nc_double_3d,                           &
          replace_var_nc_real_1d,replace_var_nc_real_2d,      &
          replace_var_nc_real_3d,                             &
          replace_var_nc_int_1d,replace_var_nc_int_2d,        &
          replace_var_nc_int_3d,                              &
          replace_var_nc_char_1d,replace_var_nc_char_2d,      &
          replace_var_nc_char_3d !< Replace 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: replace_var_nc_int !< Replace 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: replace_var_nc_int_1d !< Replace 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: replace_var_nc_int_2d !< Replace 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: replace_var_nc_int_3d !< Replace 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: replace_var_nc_real !< Replace 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: replace_var_nc_real_1d !< Replace 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: replace_var_nc_real_2d !< Replace 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: replace_var_nc_real_3d !< Replace 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: replace_var_nc_double !< Replace 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: replace_var_nc_double_1d !< Replace 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: replace_var_nc_double_2d !< Replace 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: replace_var_nc_double_3d !< Replace 1d, 2d, 3d, or 4d field from the nc file. @return
     procedure :: replace_var_nc_char !< Replace character type variable. @return
     procedure :: replace_var_nc_char_1d !< Replace character type variable. @return
     procedure :: replace_var_nc_char_2d !< Replace character type variable. @return
     procedure :: replace_var_nc_char_3d !< Replace 3D character type variable. @return
     procedure :: handle_err !< Handle netCDF errors. @return
     procedure :: convert_theta2t_2dgrid !< Convert theta T (Kelvin) to T (deg C). @return
     generic   :: add_new_var => add_new_var_2d, &
          add_new_var_3d !< Add a new 2d or 3d variable to ouput file. @return
     procedure :: add_new_var_2d !< Add a new 2d variable to output file. @return
     procedure :: add_new_var_3d !< Add a new 3d variable to output file. @return
  end type ncio

contains

  !> Open a netcdf file, set initial debug level.
  !!
  !! @param this instance of an ncio class
  !! @param filename the file to open
  !! @param action "r" for read, "w" for write
  !! @param debug_level set to non-zero for some verbose output
  !! @author Ming Hu @date 2017-11-01
  subroutine open_nc(this,filename,action,debug_level)

    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: filename
    character(len=1),intent(in) :: action
    integer,intent(in),optional :: debug_level

    integer :: ncid, status

    this%debug_level=20
    if(present(debug_level)) this%debug_level=debug_level

    this%filename=trim(filename)
    ! open existing netCDF dataset
    if(action=="r" .or. action=="R") then
       status = nf90_open(path = trim(filename), mode = nf90_nowrite, ncid = ncid)
    elseif(action=="w" .or. action=="W") then
       status = nf90_open(path = trim(filename), mode = nf90_write, ncid = ncid)
    else
       write(6,*) 'unknow action :', action
       stop 123
    endif
    if (status /= nf90_noerr) call this%handle_err(status)
    this%ncid=ncid

    if(this%debug_level>0) then
       write(6,*) '>>> open file: ',trim(this%filename)
    endif

  end subroutine open_nc

  !> Close a netcdf file.
  !!
  !! @param this instance of an ncio class
  !! @author Ming Hu org: GSD/AMB @date 2017-04-10
  subroutine close_nc(this)

    implicit none
    !
    class(ncio) :: this

    integer :: ncid, status

    ncid=this%ncid
    !
    ! close netCDF dataset
    status = nf90_close(ncid)
    if (status /= nf90_noerr) call this%handle_err(status)

  end subroutine close_nc

  !> Get attribute in wrf netcdf file.
  !!
  !! @param this instance of an ncio class
  !! @param attname name of the attribute to get
  !! @param rval return value
  !! @author Ming Hu org: GSD/AMB @date 2017-10-04
  subroutine get_att_nc_real(this,attname,rval)
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: attname
    real, intent(out) :: rval

    integer :: ncid, status

    ! open existing netCDF dataset
    ncid=this%ncid

    !  get date from exisiting NC file
    status = nf90_get_att(ncid, NF90_GLOBAL, trim(attname), rval)
    if (status /= nf90_noerr) call this%handle_err(status)
    !
  end subroutine get_att_nc_real

  !> Get integer attribute in wrf netcdf file.
  !!
  !! @param this instance of an ncio class
  !! @param attname name of the attribute to get
  !! @param ival value of attribute.
  !! @author Ming Hu org: GSD/AMB @date 2017-10-04
  subroutine get_att_nc_int(this,attname,ival)
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: attname
    integer, intent(out) :: ival

    integer :: ncid, status

    ! open existing netCDF dataset
    ncid=this%ncid

    !  get date from exisiting NC file
    status = nf90_get_att(ncid, NF90_GLOBAL, trim(attname), ival)
    if (status /= nf90_noerr) call this%handle_err(status)
    !
  end subroutine get_att_nc_int

  !> Get string attribute in wrf netcdf file.
  !!
  !! @param this instance of an ncio class
  !! @param attname name of the attribute to get
  !! @param string value of attribute.
  !! @author Ming Hu org: GSD/AMB @date 2017-10-04
  subroutine get_att_nc_string(this,attname,string)
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: attname
    character(len=*), intent(out) :: string

    integer :: ncid, status

    ! open existing netCDF dataset
    ncid=this%ncid

    !  get date from exisiting NC file
    status = nf90_get_att(ncid, NF90_GLOBAL, trim(attname), string)
    if (status /= nf90_noerr) call this%handle_err(status)
    !
  end subroutine get_att_nc_string


  !> Get dimensions in netcdf file.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] dimname name of the dimension
  !! @param[out] dimvalue length of the dimension
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine get_dim_nc(this,dimname,dimvalue)
    implicit none
    !
    class(ncio) :: this
    character(len=*), intent(in) :: dimname
    integer,intent(out) :: dimvalue

    integer :: ncid, status
    integer :: DimId

    ! open existing netCDF dataset
    ncid=this%ncid

    !  get dimension from exisiting NC file
    status = nf90_inq_dimid(ncid,trim(dimname), DimId)
    if (status /= nf90_noerr) call this%handle_err(status)
    status = nf90_Inquire_Dimension(ncid, DimId, len = dimvalue)
    if (status /= nf90_noerr) call this%handle_err(status)
    !
  end subroutine get_dim_nc

  !> Replace 1D character type variable
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[in] field replacement field
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine replace_var_nc_char_1d(this,varname,nd1,field)

    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1              !  size of array dval
    character, intent(in) :: field(nd1)     !  values of the field read in
    integer :: ilength
    !
    character*40,parameter :: thissubname='replace_var_nc_char_1d'
    !
    integer :: i
    !
    !
    ilength=nd1
    !
    if(this%debug_level>100) then
       write(6,*) trim(thissubname),' show samples:'
       write(6,*) (field(i),i=1,min(nd1,10))
    endif

    call this%replace_var_nc_char(varname,ilength,field)
    !
  end subroutine replace_var_nc_char_1d

  !> Replace 2D character type variable
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[in] nd2 length of second dimension
  !! @param[in] field replacement field
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine replace_var_nc_char_2d(this,varname,nd1,nd2,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1,nd2          !  size of array dval
    character, intent(in) :: field(nd1,nd2) !  values of the field read in
    integer :: ilength
    !
    character,allocatable :: temp(:) 
    !
    character*40,parameter :: thissubname='replace_var_nc_char_2d'
    !
    integer :: i,j,k
    integer :: istart,iend
    !
    !
    ilength=nd1*nd2
    allocate(temp(ilength))

    do j=1,nd2
       istart=(j-1)*nd1+1
       iend=(j-1)*nd1+nd1
       temp(istart:iend)=field(:,j)
    enddo
    !
    if(this%debug_level>100) then
       write(6,*) trim(thissubname),' show samples:'
       write(6,*) field(1,1)
    endif
    !
    call this%replace_var_nc_char(varname,ilength,temp)

    deallocate(temp)
    !
  end subroutine replace_var_nc_char_2d

  !> Replace 3D character type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[in] nd2 length of second dimension
  !! @param[in] nd3 length of third dimension
  !! @param[in] field replacement field
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine replace_var_nc_char_3d(this,varname,nd1,nd2,nd3,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1,nd2,nd3      !  size of array dval
    character, intent(in) :: field(nd1,nd2,nd3) !  values of the field read in
    integer :: ilength
    !
    character,allocatable :: temp(:) 
    !
    character*40,parameter :: thissubname='replace_var_nc_char_3d'
    !
    integer :: i,j,k
    integer :: length2d
    integer :: istart,iend
    !
    !
    length2d=nd1*nd2
    ilength=length2d*nd3
    allocate(temp(ilength))

    do k=1,nd3
       do j=1,nd2
          istart=(k-1)*length2d+(j-1)*nd1+1
          iend  =(k-1)*length2d+(j-1)*nd1+nd1
          temp(istart:iend)=field(:,j,k)
       enddo
    enddo
    !
    if(this%debug_level>100) then
       write(6,*) trim(thissubname),' show samples:'
       write(6,*) field(1,1,1)
    endif

    call this%replace_var_nc_char(varname,ilength,temp)

    deallocate(temp)
    !
  end subroutine replace_var_nc_char_3d

  !> Replace character type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] ilength length of array
  !! @param[in] field replacement field
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine replace_var_nc_char(this,varname,ilength,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: ilength          !  size of array dval
    character, intent(in) :: field(ilength)   !  values of the field read in
    !
    integer :: ncid
    ! 
    integer :: status
    integer :: varid
    integer :: ends(4),start(4)

    integer :: length4d,length3d,length2d
    integer :: nDims,ndim
    integer :: dimids(4)
    integer :: xtype
    character*40 :: dimname

    character*40,parameter :: thissubname='replace_var_nc_char'
    !
    integer :: i,k
    !
    !
    ncid=this%ncid

    ! get variable IDs
    status = nf90_inq_varid(ncid, trim(varname), VarId)
    if(status /= nf90_NoErr) call this%handle_err(status)

    !  get dimensions
    ends=1
    start=1
    this%ends=1

    this%dimname="                           "
    ! get variable type
    status = nf90_inquire_variable(ncid, VarId, xtype = xtype)
    if(status /= nf90_NoErr) call this%handle_err(status)
    if(xtype==NF90_CHAR) then
       this%xtype=xtype
    else
       write(6,*) trim(thissubname),' ERROR: wrong data type, expect ',NF90_INT,' but read in ',xtype
       stop 123
    endif

    ! get dimension size
    status = nf90_inquire_variable(ncid, VarId, ndims = nDims)
    if(status /= nf90_NoErr) call this%handle_err(status)
    this%ndims=nDims
    !
    status = nf90_inquire_variable(ncid, VarId, dimids = dimids(1:nDims))
    if(status /= nf90_NoErr) call this%handle_err(status)
    do i=1,nDims
       dimname="       "
       status = nf90_inquire_dimension(ncid, dimids(i), dimname, len = ndim)
       if (status /= nf90_noerr) call this%handle_err(status)
       ends(i)=ndim
       this%ends(i)=ends(i)
       this%dimname(i)=trim(dimname)
       if(this%ends(i) < 1) then
          write(6,*) trim(thissubname),' Error, ends dimension should larger than 0 :', ends(i)
          stop 1234
       endif
    enddo
    length2d=ends(1)*ends(2)
    length3d=length2d*ends(3)
    length4d=length3d*ends(4)
    if(ilength .ne. length4d) then
       write(6,*) trim(thissubname),'ERROR: ',ilength,' should equal to ',length4d
       stop 123
    endif
    !
    if(nDims <=4 ) then
       status = nf90_put_var(ncid, VarId, field, &
            start = start(1:4) , &
            count = ends(1:4))
       if(status /= nf90_NoErr) call this%handle_err(status)
    else
       write(6,*) trim(thissubname),'Error: too many dimensions:',nDims
       stop 1234
    endif
    !
    if(this%debug_level>0) then
       write(6,'(a,a)') '>>>replace variable: ',trim(varname)
    endif
    if(this%debug_level>10) then
       write(6,'(8x,a,I10)') 'data type : ',this%xtype
       write(6,'(8x,a,I10)') 'dimension size: ',this%nDims
       do i=1,this%nDims
          write(6,'(8x,a,I5,I10,2x,a)') 'rank, ends, name=',i,this%ends(i),trim(this%dimname(i))
       enddo
    endif
    !
  end subroutine replace_var_nc_char
  !--- replace_var_nc_char

  !> Replace 1D real type variable
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[in] field replacement field
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine replace_var_nc_real_1d(this,varname,nd1,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1              !  size of array dval
    real(4), intent(in) :: field(nd1)     !  values of the field read in
    integer :: ilength
    !
    character*40,parameter :: thissubname='replace_var_nc_real_1d'
    !
    integer :: i
    !
    !
    ilength=nd1
    !
    if(this%debug_level>100) then
       write(6,*) trim(thissubname),' show samples:'
       write(6,*) (field(i),i=1,min(nd1,10))
    endif
    !
    call this%replace_var_nc_real(varname,ilength,field)
    !
  end subroutine replace_var_nc_real_1d

  !> Replace 2D real type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[in] nd2 length of second dimension
  !! @param[in] field replacement field
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine replace_var_nc_real_2d(this,varname,nd1,nd2,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1,nd2          !  size of array dval
    real(4), intent(in) :: field(nd1,nd2) !  values of the field read in
    integer :: ilength
    !
    real(4),allocatable :: temp(:) 
    !
    character*40,parameter :: thissubname='replace_var_nc_real_2d'
    !
    integer :: i,j,k
    integer :: istart,iend
    !
    !
    ilength=nd1*nd2
    allocate(temp(ilength))

    do j=1,nd2
       istart=(j-1)*nd1+1
       iend=(j-1)*nd1+nd1
       temp(istart:iend)=field(:,j)
    enddo
    !
    if(this%debug_level>100) then
       write(6,*) trim(thissubname),' show samples:'
       write(6,*) 'max,min:',maxval(field(:,:)),minval(field(:,:))
    endif

    call this%replace_var_nc_real(varname,ilength,temp)

    deallocate(temp)
    !
  end subroutine replace_var_nc_real_2d

  !> Replace 3D real type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[in] nd2 length of second dimension
  !! @param[in] nd3 length of third dimension
  !! @param[in] field replacement field
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine replace_var_nc_real_3d(this,varname,nd1,nd2,nd3,field)
    !
    ! read in one field 
    !
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1,nd2,nd3      !  size of array dval
    real(4), intent(in) :: field(nd1,nd2,nd3) !  values of the field read in
    integer :: ilength
    !
    real(4),allocatable :: temp(:) 
    !
    character*40,parameter :: thissubname='replace_var_nc_real_3d'
    !
    integer :: i,j,k
    integer :: length2d
    integer :: istart,iend
    !
    !
    length2d=nd1*nd2
    ilength=length2d*nd3
    allocate(temp(ilength))

    do k=1,nd3
       do j=1,nd2
          istart=(k-1)*length2d+(j-1)*nd1+1
          iend  =(k-1)*length2d+(j-1)*nd1+nd1
          temp(istart:iend)=field(:,j,k)
       enddo
    enddo
    !
    if(this%debug_level>100) then
       write(6,*) trim(thissubname),' show samples:'
       do k=1,nd3
          write(6,*) 'k,max,min:',k,maxval(field(:,:,k)),minval(field(:,:,k))
       enddo
    endif

    call this%replace_var_nc_real(varname,ilength,temp)

    deallocate(temp)
    !
  end subroutine replace_var_nc_real_3d

  !> Replace real type variable
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] ilength length of array
  !! @param[in] field replacement field
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine replace_var_nc_real(this,varname,ilength,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: ilength          !  size of array dval
    real(4), intent(in) :: field(ilength)   !  values of the field read in
    !
    integer :: ncid
    ! 
    integer :: status
    integer :: varid
    integer :: ends(4),start(4)

    integer :: length4d,length3d,length2d
    integer :: nDims,ndim
    integer :: dimids(4)
    integer :: xtype
    character*40 :: dimname

    character*40,parameter :: thissubname='replace_var_nc_real'
    !
    integer :: i,k
    !
    !
    ncid=this%ncid

    ! get variable IDs
    status = nf90_inq_varid(ncid, trim(varname), VarId)
    if(status /= nf90_NoErr) call this%handle_err(status)

    !  get dimensions
    ends=1
    start=1
    this%ends=1

    this%dimname="                           "
    ! get variable type
    status = nf90_inquire_variable(ncid, VarId, xtype = xtype)
    if(status /= nf90_NoErr) call this%handle_err(status)
    if(xtype==NF90_FLOAT) then
       this%xtype=xtype
    else
       write(6,*) trim(thissubname),' ERROR: wrong data type, expect ',NF90_INT,' but read in ',xtype
       stop 123
    endif

    ! get dimension size
    status = nf90_inquire_variable(ncid, VarId, ndims = nDims)
    if(status /= nf90_NoErr) call this%handle_err(status)
    this%ndims=nDims
    !
    status = nf90_inquire_variable(ncid, VarId, dimids = dimids(1:nDims))
    if(status /= nf90_NoErr) call this%handle_err(status)
    do i=1,nDims
       dimname="       "
       status = nf90_inquire_dimension(ncid, dimids(i), dimname, len = ndim)
       if (status /= nf90_noerr) call this%handle_err(status)
       ends(i)=ndim
       this%ends(i)=ends(i)
       this%dimname(i)=trim(dimname)
       if(this%ends(i) < 1) then
          write(6,*) trim(thissubname),' Error, ends dimension should larger than 0 :', ends(i)
          stop 1234
       endif
    enddo
    length2d=ends(1)*ends(2)
    length3d=length2d*ends(3)
    length4d=length3d*ends(4)
    if(ilength .ne. length4d) then
       write(6,*) trim(thissubname),'ERROR: ',ilength,' should equal to ',length4d
       stop 123
    endif
    !
    if(nDims <=4 ) then
       status = nf90_put_var(ncid, VarId, field, &
            start = start(1:4) , &
            count = ends(1:4))
       if(status /= nf90_NoErr) call this%handle_err(status)
    else
       write(6,*) trim(thissubname),'Error: too many dimensions:',nDims
       stop 1234
    endif
    !
    if(this%debug_level>0) then
       write(6,'(a,a)') '>>>replace variable: ',trim(varname)
    endif
    if(this%debug_level>10) then
       write(6,'(8x,a,I10)') 'data type : ',this%xtype
       write(6,'(8x,a,I10)') 'dimension size: ',this%nDims
       do i=1,this%nDims
          write(6,'(8x,a,I5,I10,2x,a)') 'rank, ends, name=',i,this%ends(i),trim(this%dimname(i))
       enddo
    endif
    !
  end subroutine replace_var_nc_real

  !> Replace 1D double type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[in] field replacement field
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine replace_var_nc_double_1d(this,varname,nd1,field)
    !
    ! read in one field 
    !
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1              !  size of array dval
    real(8), intent(in) :: field(nd1)     !  values of the field read in
    integer :: ilength
    !
    character*40,parameter :: thissubname='replace_var_nc_double_1d'
    !
    integer :: i
    !
    !
    ilength=nd1
    !
    if(this%debug_level>100) then
       write(6,*) trim(thissubname),' show samples:'
       write(6,*) (field(i),i=1,min(nd1,10))
    endif
    !
    call this%replace_var_nc_double(varname,ilength,field)
    !
  end subroutine replace_var_nc_double_1d

  !> Replace 2D double type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[in] nd2 length of second dimension
  !! @param[in] field replacement field
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine replace_var_nc_double_2d(this,varname,nd1,nd2,field)
    !
    ! read in one field 
    !
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1,nd2          !  size of array dval
    real(8), intent(in) :: field(nd1,nd2) !  values of the field read in
    integer :: ilength
    !
    real(8),allocatable :: temp(:) 
    !
    character*40,parameter :: thissubname='replace_var_nc_double_2d'
    !
    integer :: i,j,k
    integer :: istart,iend
    !
    !
    ilength=nd1*nd2
    allocate(temp(ilength))

    do j=1,nd2
       istart=(j-1)*nd1+1
       iend=(j-1)*nd1+nd1
       temp(istart:iend)=field(:,j)
    enddo
    !
    if(this%debug_level>100) then
       write(6,*) trim(thissubname),' show samples:'
       write(6,*) 'max,min:',maxval(field(:,:)),minval(field(:,:))
    endif

    call this%replace_var_nc_double(varname,ilength,temp)

    deallocate(temp)
    !
  end subroutine replace_var_nc_double_2d

  !> Replace 3D double type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[in] nd2 length of second dimension
  !! @param[in] nd3 length of third dimension
  !! @param[in] field replacement field
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine replace_var_nc_double_3d(this,varname,nd1,nd2,nd3,field)
    !
    ! read in one field 
    !
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1,nd2,nd3      !  size of array dval
    real(8), intent(in) :: field(nd1,nd2,nd3) !  values of the field read in
    integer :: ilength
    !
    real(8),allocatable :: temp(:) 
    !
    character*40,parameter :: thissubname='replace_var_nc_double_3d'
    !
    integer :: i,j,k
    integer :: length2d
    integer :: istart,iend
    !
    !
    length2d=nd1*nd2
    ilength=length2d*nd3
    allocate(temp(ilength))

    do k=1,nd3
       do j=1,nd2
          istart=(k-1)*length2d+(j-1)*nd1+1
          iend  =(k-1)*length2d+(j-1)*nd1+nd1
          temp(istart:iend)=field(:,j,k)
       enddo
    enddo
    !
    if(this%debug_level>100) then
       write(6,*) trim(thissubname),' show samples:'
       do k=1,nd3
          write(6,*) 'k,max,min:',k,maxval(field(:,:,k)),minval(field(:,:,k))
       enddo
    endif

    call this%replace_var_nc_double(varname,ilength,temp)

    deallocate(temp)
    !
  end subroutine replace_var_nc_double_3d
  !

  !> Replace double type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] ilength size of array
  !! @param[in] field replacement field
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine replace_var_nc_double(this,varname,ilength,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: ilength          !  size of array dval
    real(8), intent(in) :: field(ilength)   !  values of the field read in
    !
    integer :: ncid
    ! 
    integer :: status
    integer :: varid
    integer :: ends(4),start(4)

    integer :: length4d,length3d,length2d
    integer :: nDims,ndim
    integer :: dimids(4)
    integer :: xtype
    character*40 :: dimname

    character*40,parameter :: thissubname='replace_var_nc_double'
    !
    integer :: i,k
    !
    !
    ncid=this%ncid

    ! get variable IDs
    status = nf90_inq_varid(ncid, trim(varname), VarId)
    if(status /= nf90_NoErr) call this%handle_err(status)

    !  get dimensions
    ends=1
    start=1
    this%ends=1

    this%dimname="                           "
    ! get variable type
    status = nf90_inquire_variable(ncid, VarId, xtype = xtype)
    if(status /= nf90_NoErr) call this%handle_err(status)
    if(xtype==NF90_DOUBLE) then
       this%xtype=xtype
    else
       write(6,*) trim(thissubname),' ERROR: wrong data type, expect ',NF90_INT,' but read in ',xtype
       stop 123
    endif

    ! get dimension size
    status = nf90_inquire_variable(ncid, VarId, ndims = nDims)
    if(status /= nf90_NoErr) call this%handle_err(status)
    this%ndims=nDims
    !
    status = nf90_inquire_variable(ncid, VarId, dimids = dimids(1:nDims))
    if(status /= nf90_NoErr) call this%handle_err(status)
    do i=1,nDims
       dimname="       "
       status = nf90_inquire_dimension(ncid, dimids(i), dimname, len = ndim)
       if (status /= nf90_noerr) call this%handle_err(status)
       ends(i)=ndim
       this%ends(i)=ends(i)
       this%dimname(i)=trim(dimname)
       if(this%ends(i) < 1) then
          write(6,*) trim(thissubname),' Error, ends dimension should larger than 0 :', ends(i)
          stop 1234
       endif
    enddo
    length2d=ends(1)*ends(2)
    length3d=length2d*ends(3)
    length4d=length3d*ends(4)
    if(ilength .ne. length4d) then
       write(6,*) trim(thissubname),'ERROR: ',ilength,' should equal to ',length4d
       stop 123
    endif
    !
    if(nDims <=4 ) then
       status = nf90_put_var(ncid, VarId, field, &
            start = start(1:4) , &
            count = ends(1:4))
       if(status /= nf90_NoErr) call this%handle_err(status)
    else
       write(6,*) trim(thissubname),'Error: too many dimensions:',nDims
       stop 1234
    endif
    !
    if(this%debug_level>0) then
       write(6,'(a,a)') '>>>replace variable: ',trim(varname)
    endif
    if(this%debug_level>10) then
       write(6,'(8x,a,I10)') 'data type : ',this%xtype
       write(6,'(8x,a,I10)') 'dimension size: ',this%nDims
       do i=1,this%nDims
          write(6,'(8x,a,I5,I10,2x,a)') 'rank, ends, name=',i,this%ends(i),trim(this%dimname(i))
       enddo
    endif
    !
  end subroutine replace_var_nc_double

  !> Replace 1D integer type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 lenth of first dimension
  !! @param[in] field replacement field
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine replace_var_nc_int_1d(this,varname,nd1,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1              !  size of array dval
    integer, intent(in) :: field(nd1)     !  values of the field read in
    integer :: ilength
    !
    character*40,parameter :: thissubname='get_var_nc_int_1d'
    !
    integer :: i
    !
    !
    ilength=nd1
    !
    if(this%debug_level>100) then
       write(6,*) trim(thissubname),' show samples:'
       write(6,*) (field(i),i=1,min(nd1,10))
    endif

    call this%replace_var_nc_int(varname,ilength,field)
    !
  end subroutine replace_var_nc_int_1d

  !> Replace 2D integer type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[in] nd2 length of second dimension
  !! @param[in] field replacement field
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine replace_var_nc_int_2d(this,varname,nd1,nd2,field)
    !
    ! read in one field 
    !
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1,nd2          !  size of array dval
    integer, intent(in) :: field(nd1,nd2) !  values of the field read in
    integer :: ilength
    !
    integer,allocatable :: temp(:) 
    !
    character*40,parameter :: thissubname='replace_var_nc_int_2d'
    !
    integer :: i,j,k
    integer :: istart,iend
    !
    !
    ilength=nd1*nd2
    allocate(temp(ilength))

    do j=1,nd2
       istart=(j-1)*nd1+1
       iend=(j-1)*nd1+nd1
       temp(istart:iend)=field(:,j)
    enddo
    !
    if(this%debug_level>100) then
       write(6,*) trim(thissubname),' show samples:'
       write(6,*) 'max,min:',maxval(field(:,:)),minval(field(:,:))
    endif

    call this%replace_var_nc_int(varname,ilength,temp)

    deallocate(temp)
    !
  end subroutine replace_var_nc_int_2d

  !> Replace 3D integer type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[in] nd2 length of second dimension
  !! @param[in] nd3 length of third dimension
  !! @param[in] field replacement field
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine replace_var_nc_int_3d(this,varname,nd1,nd2,nd3,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1,nd2,nd3      !  size of array dval
    integer, intent(in) :: field(nd1,nd2,nd3) !  values of the field read in
    integer :: ilength
    !
    integer,allocatable :: temp(:) 
    !
    character*40,parameter :: thissubname='replace_var_nc_int_3d'
    !
    integer :: i,j,k
    integer :: length2d
    integer :: istart,iend
    !
    !
    length2d=nd1*nd2
    ilength=length2d*nd3
    allocate(temp(ilength))

    do k=1,nd3
       do j=1,nd2
          istart=(k-1)*length2d+(j-1)*nd1+1
          iend  =(k-1)*length2d+(j-1)*nd1+nd1
          temp(istart:iend)=field(:,j,k)
       enddo
    enddo
    !
    if(this%debug_level>100) then
       write(6,*) trim(thissubname),' show samples:'
       do k=1,nd3
          write(6,*) 'k,max,min:',k,maxval(field(:,:,k)),minval(field(:,:,k))
       enddo
    endif

    call this%replace_var_nc_int(varname,ilength,temp)

    deallocate(temp)
    !
  end subroutine replace_var_nc_int_3d

  !> Replace integer type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] ilength size of array
  !! @param[in] field replacement field
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine replace_var_nc_int(this,varname,ilength,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: ilength          !  size of array dval
    integer, intent(in) :: field(ilength)   !  values of the field read in
    !
    integer :: ncid
    ! 
    integer :: status
    integer :: varid
    integer :: ends(4),start(4)

    integer :: length4d,length3d,length2d
    integer :: nDims,ndim
    integer :: dimids(4)
    integer :: xtype
    character*40 :: dimname

    character*40,parameter :: thissubname='replace_var_nc_int'
    !
    integer :: i,k
    !
    !
    ncid=this%ncid

    ! get variable IDs
    status = nf90_inq_varid(ncid, trim(varname), VarId)
    if(status /= nf90_NoErr) call this%handle_err(status)

    !  get dimensions
    ends=1
    start=1
    this%ends=1

    this%dimname="                           "
    ! get variable type
    status = nf90_inquire_variable(ncid, VarId, xtype = xtype)
    if(status /= nf90_NoErr) call this%handle_err(status)
    if(xtype==NF90_INT) then
       this%xtype=xtype
    else
       write(6,*) trim(thissubname),' ERROR: wrong data type, expect ',NF90_INT,' but read in ',xtype
       stop 123
    endif

    ! get dimension size
    status = nf90_inquire_variable(ncid, VarId, ndims = nDims)
    if(status /= nf90_NoErr) call this%handle_err(status)
    this%ndims=nDims
    !
    status = nf90_inquire_variable(ncid, VarId, dimids = dimids(1:nDims))
    if(status /= nf90_NoErr) call this%handle_err(status)
    do i=1,nDims
       dimname="       "
       status = nf90_inquire_dimension(ncid, dimids(i), dimname, len = ndim)
       if (status /= nf90_noerr) call this%handle_err(status)
       ends(i)=ndim
       this%ends(i)=ends(i)
       this%dimname(i)=trim(dimname)
       if(this%ends(i) < 1) then
          write(6,*) trim(thissubname),' Error, ends dimension should larger than 0 :', ends(i)
          stop 1234
       endif
    enddo
    length2d=ends(1)*ends(2)
    length3d=length2d*ends(3)
    length4d=length3d*ends(4)
    if(ilength .ne. length4d) then
       write(6,*) trim(thissubname),'ERROR: ',ilength,' should equal to ',length4d
       stop 123
    endif
    !
    if(nDims <=4 ) then
       status = nf90_put_var(ncid, VarId, field, &
            start = start(1:4) , &
            count = ends(1:4))
       if(status /= nf90_NoErr) call this%handle_err(status)
    else
       write(6,*) trim(thissubname),'Error: too many dimensions:',nDims
       stop 1634
    endif
    !
    if(this%debug_level>0) then
       write(6,'(a,a)') '>>>replace variable: ',trim(varname)
    endif
    if(this%debug_level>10) then
       write(6,'(8x,a,I10)') 'data type : ',this%xtype
       write(6,'(8x,a,I10)') 'dimension size: ',this%nDims
       do i=1,this%nDims
          write(6,'(8x,a,I5,I10,2x,a)') 'rank, ends, name=',i,this%ends(i),trim(this%dimname(i))
       enddo
    endif
    !
  end subroutine replace_var_nc_int

  !> Read in 1D double type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 lenth of first dimension
  !! @param[out] field output variable
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine get_var_nc_double_1d(this,varname,nd1,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1              !  size of array dval
    real(8), intent(out) :: field(nd1)     !  values of the field read in
    integer :: ilength
    !
    character*40,parameter :: thissubname='get_var_nc_double_1d'
    !
    integer :: i
    !
    !
    ilength=nd1
    call this%get_var_nc_double(varname,ilength,field)
    !
    if(nd1==this%ends(1)) then
       if(this%debug_level>100) then
          write(6,*) trim(thissubname),' show samples:'
          write(6,*) (field(i),i=1,min(nd1,10))
       endif
    else
       write(6,*) trim(thissubname),' ERROR: dimension does not match.'
    endif
    !
  end subroutine get_var_nc_double_1d

  !> Read in 2D double type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[in] nd2 length of second dimension
  !! @param[out] field output variable
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine get_var_nc_double_2d(this,varname,nd1,nd2,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1,nd2          !  size of array dval
    real(8), intent(out) :: field(nd1,nd2) !  values of the field read in
    integer :: ilength
    !
    real(8),allocatable :: temp(:) 
    !
    character*40,parameter :: thissubname='get_var_nc_double_2d'
    !
    integer :: i,j,k
    integer :: istart,iend
    !
    !
    ilength=nd1*nd2
    allocate(temp(ilength))

    call this%get_var_nc_double(varname,ilength,temp)

    if(nd1==this%ends(1) .and. nd2==this%ends(2)) then
       do j=1,nd2
          istart=(j-1)*nd1+1
          iend=(j-1)*nd1+nd1
          field(:,j)=temp(istart:iend)
       enddo
       !
!       if(this%debug_level>100) then
!          write(*,*) trim(thissubname),' show samples:'
!          write(*,*) 'max,min:',maxval(field(:,:)),minval(field(:,:))
!       endif
    else
       write(6,*) trim(thissubname),' ERROR: dimension does not match.'
       write(6,*) nd1,this%ends(1),nd2,this%ends(2)
    endif
    deallocate(temp)
    !
  end subroutine get_var_nc_double_2d

  !> Read in 3D double type field.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[in] nd2 length of second dimension
  !! @param[in] nd3 length of third dimension
  !! @param[out] field output variable
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine get_var_nc_double_3d(this,varname,nd1,nd2,nd3,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1,nd2,nd3      !  size of array dval
    real(8), intent(out) :: field(nd1,nd2,nd3) !  values of the field read in
    integer :: ilength
    !
    real(8),allocatable :: temp(:) 
    !
    character*40,parameter :: thissubname='get_var_nc_double_3d'
    !
    integer :: i,j,k
    integer :: length2d
    integer :: istart,iend
    !
    !
    length2d=nd1*nd2
    ilength=length2d*nd3
    allocate(temp(ilength))

    call this%get_var_nc_double(varname,ilength,temp)

    if(nd1==this%ends(1) .and. nd2==this%ends(2) .and. nd3==this%ends(3)) then
       do k=1,nd3
          do j=1,nd2
             istart=(k-1)*length2d+(j-1)*nd1+1
             iend  =(k-1)*length2d+(j-1)*nd1+nd1
             field(:,j,k)=temp(istart:iend)
          enddo
       enddo
       !
!       if(this%debug_level>100) then
!          write(*,*) trim(thissubname),' show samples:'
!          do k=1,nd3
!             write(*,*) 'k,max,min:',k,maxval(field(:,:,k)),minval(field(:,:,k))
!          enddo
!       endif
    else
       write(6,*) trim(thissubname),' ERROR: dimension does not match.'
       write(6,*) nd1,this%ends(1),nd2,this%ends(2),nd3,this%ends(3)
    endif
    deallocate(temp)
    !
  end subroutine get_var_nc_double_3d

  !> Read in double type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] ilength size of array
  !! @param[out] field output variable
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine get_var_nc_double(this,varname,ilength,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: ilength          !  size of array dval
    real(8), intent(out) :: field(ilength)   !  values of the field read in
    !
    integer :: ncid
    ! 
    integer :: status
    integer :: varid
    integer :: ends(4),start(4)

    integer :: length4d,length3d,length2d
    integer :: nDims,ndim
    integer :: dimids(4)
    integer :: xtype
    character*40 :: dimname

    character*40,parameter :: thissubname='get_var_nc_double'
    !
    integer :: i,k
    !
    !
    ncid=this%ncid

    ! get variable IDs
    status = nf90_inq_varid(ncid, trim(varname), VarId)
    if(status /= nf90_NoErr) call this%handle_err(status)

    !  get dimensions
    ends=1
    start=1
    this%ends=1

    this%dimname="                           "
    ! get variable type
    status = nf90_inquire_variable(ncid, VarId, xtype = xtype)
    if(status /= nf90_NoErr) call this%handle_err(status)
    if(xtype==NF90_DOUBLE) then
       this%xtype=xtype
    else
       write(6,*) trim(thissubname),' ERROR: wrong data type, expect ',NF90_DOUBLE,' but read in ',xtype
       stop 123
    endif

    ! get dimension size
    status = nf90_inquire_variable(ncid, VarId, ndims = nDims)
    if(status /= nf90_NoErr) call this%handle_err(status)
    this%ndims=nDims
    !
    status = nf90_inquire_variable(ncid, VarId, dimids = dimids(1:nDims))
    if(status /= nf90_NoErr) call this%handle_err(status)
    do i=1,nDims
       dimname="       "
       write(6,*) 'dimids(i) = ', dimids(i)
       status = nf90_inquire_dimension(ncid, dimids(i), dimname, len = ndim)
       if (status /= nf90_noerr) call this%handle_err(status)
       ends(i)=ndim
       this%ends(i)=ends(i)
       this%dimname(i)=trim(dimname)
       if(this%ends(i) < 1) then
          write(6,*) trim(thissubname),' Error, ends dimension should larger than 0 :', ends(i)
          stop 1234
       endif
    enddo
    length2d=ends(1)*ends(2)
    length3d=length2d*ends(3)
    length4d=length3d*ends(4)
    if(ilength .ne. length4d) then
       write(6,*) trim(thissubname),'ERROR: ',ilength,' should equal to ',length4d
       stop 123
    endif
    !
    if(nDims <=4 ) then
       status = nf90_get_var(ncid, VarId, field, &
            start = start(1:4) , &
            count = ends(1:4))
       if(status /= nf90_NoErr) call this%handle_err(status)
    else
       write(6,*) trim(thissubname),'Error: too many dimensions:',nDims
       stop 1234
    endif
    !
    if(this%debug_level>0) then
       write(6,'(a,a)') '>>>read in variable: ',trim(varname)
    endif
    if(this%debug_level>10) then
       write(6,'(a,I10)') '        data type : ',this%xtype
       write(6,'(a,I10)')'        dimension size: ',this%nDims
       do i=1,this%nDims
          write(6,'(a,I5,I10,2x,a)') '        rank, ends, name=',i,this%ends(i),trim(this%dimname(i))
       enddo
    endif
    !
  end subroutine get_var_nc_double

  !> Read in 1D real type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[out] field output variable
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine get_var_nc_real_1d(this,varname,nd1,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1              !  size of array dval
    real(4), intent(out) :: field(nd1)     !  values of the field read in
    integer :: ilength
    !
    character*40,parameter :: thissubname='get_var_nc_real_1d'
    !
    integer :: i
    !
    !
    ilength=nd1
    call this%get_var_nc_real(varname,ilength,field)
    !
    if(nd1==this%ends(1)) then
       if(this%debug_level>100) then
          write(6,*) trim(thissubname),' show samples:'
          write(6,*) (field(i),i=1,min(nd1,10))
       endif
    else
       write(6,*) trim(thissubname),' ERROR: dimension does not match.'
    endif
    !
  end subroutine get_var_nc_real_1d

  !> Read in 2D real type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[in] nd2 length of second dimension
  !! @param[out] field output variable
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine get_var_nc_real_2d(this,varname,nd1,nd2,field)
    !
    ! read in one field 
    !
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1,nd2          !  size of array dval
    real(4), intent(out) :: field(nd1,nd2) !  values of the field read in
    integer :: ilength
    !
    real(4),allocatable :: temp(:) 
    !
    character*40,parameter :: thissubname='get_var_nc_real_2d'
    !
    integer :: i,j,k
    integer :: istart,iend
    !
    !
    ilength=nd1*nd2
    allocate(temp(ilength))

    call this%get_var_nc_real(varname,ilength,temp)

    if(nd1==this%ends(1) .and. nd2==this%ends(2)) then
       do j=1,nd2
          istart=(j-1)*nd1+1
          iend=(j-1)*nd1+nd1
          field(:,j)=temp(istart:iend)
       enddo
       !
       if(this%debug_level>100) then
          write(6,*) trim(thissubname),' show samples:'
          write(6,*) 'max,min:',maxval(field(:,:)),minval(field(:,:))
       endif
    else
       write(6,*) trim(thissubname),' ERROR: dimension does not match.'
       write(6,*) nd1,this%ends(1),nd2,this%ends(2)
    endif
    deallocate(temp)
    !
  end subroutine get_var_nc_real_2d

  !> Read in 3D real type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[in] nd2 length of second dimension
  !! @param[in] nd3 length of third dimension
  !! @param[out] field output variable
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine get_var_nc_real_3d(this,varname,nd1,nd2,nd3,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1,nd2,nd3      !  size of array dval
    real(4), intent(out) :: field(nd1,nd2,nd3) !  values of the field read in
    integer :: ilength
    !
    real(4),allocatable :: temp(:) 
    !
    character*40,parameter :: thissubname='get_var_nc_real_3d'
    !
    integer :: i,j,k
    integer :: length2d
    integer :: istart,iend
    !
    !
    length2d=nd1*nd2
    ilength=length2d*nd3
    allocate(temp(ilength))

    call this%get_var_nc_real(varname,ilength,temp)

    if(nd1==this%ends(1) .and. nd2==this%ends(2) .and. nd3==this%ends(3)) then
       do k=1,nd3
          do j=1,nd2
             istart=(k-1)*length2d+(j-1)*nd1+1
             iend  =(k-1)*length2d+(j-1)*nd1+nd1
             field(:,j,k)=temp(istart:iend)
          enddo
       enddo
       !
       if(this%debug_level>100) then
          write(6,*) trim(thissubname),' show samples:'
          do k=1,nd3
             write(6,*) 'k,max,min:',k,maxval(field(:,:,k)),minval(field(:,:,k))
          enddo
       endif
    else
       write(6,*) trim(thissubname),' ERROR: dimension does not match.'
       write(6,*) nd1,this%ends(1),nd2,this%ends(2),nd3,this%ends(3)
    endif
    deallocate(temp)
    !
  end subroutine get_var_nc_real_3d

  !> Read in real type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] ilength size of array
  !! @param[out] field output variable
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine get_var_nc_real(this,varname,ilength,field)
    !
    ! read in one field 
    !
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: ilength          !  size of array dval
    real(4), intent(out) :: field(ilength)   !  values of the field read in
    !
    integer :: ncid
    ! 
    integer :: status
    integer :: varid
    integer :: ends(4),start(4)

    integer :: length4d,length3d,length2d
    integer :: nDims,ndim
    integer :: dimids(4)
    integer :: xtype
    character*40 :: dimname

    character*40,parameter :: thissubname='get_var_nc_real'
    !
    integer :: i,k
    !
    !
    ncid=this%ncid

    ! get variable IDs
    status = nf90_inq_varid(ncid, trim(varname), VarId)
    if(status /= nf90_NoErr) call this%handle_err(status)

    !  get dimensions
    ends=1
    start=1
    this%ends=1

    this%dimname="                           "
    ! get variable type
    status = nf90_inquire_variable(ncid, VarId, xtype = xtype)
    if(status /= nf90_NoErr) call this%handle_err(status)
    if(xtype==NF90_FLOAT) then
       this%xtype=xtype
    else
       write(6,*) trim(thissubname),' ERROR: wrong data type, expect ',NF90_FLOAT,' but read in ',xtype
       stop 123
    endif

    ! get dimension size
    status = nf90_inquire_variable(ncid, VarId, ndims = nDims)
    if(status /= nf90_NoErr) call this%handle_err(status)
    this%ndims=nDims
    !
    status = nf90_inquire_variable(ncid, VarId, dimids = dimids(1:nDims))
    if(status /= nf90_NoErr) call this%handle_err(status)
    do i=1,nDims
       dimname="       "
       status = nf90_inquire_dimension(ncid, dimids(i), dimname, len = ndim)
       if (status /= nf90_noerr) call this%handle_err(status)
       ends(i)=ndim
       this%ends(i)=ends(i)
       this%dimname(i)=trim(dimname)
       if(this%ends(i) < 1) then
          write(6,*) trim(thissubname),' Error, ends dimension should larger than 0 :', ends(i)
          stop 1234
       endif
    enddo
    length2d=ends(1)*ends(2)
    length3d=length2d*ends(3)
    length4d=length3d*ends(4)
    if(ilength .ne. length4d) then
       write(6,*) trim(thissubname),'ERROR: ',ilength,' should equal to ',length4d
       stop 123
    endif
    !
    if(nDims <=4 ) then
       status = nf90_get_var(ncid, VarId, field, &
            start = start(1:4) , &
            count = ends(1:4))
       if(status /= nf90_NoErr) call this%handle_err(status)
    else
       write(6,*) trim(thissubname),'Error: too many dimensions:',nDims
       stop 1234
    endif
    !
    if(this%debug_level>0) then
       write(6,'(a,a)') '>>>read in variable: ',trim(varname)
    endif
    if(this%debug_level>10) then
       write(6,'(8x,a,I10)') 'data type : ',this%xtype
       write(6,'(8x,a,I10)') 'dimension size: ',this%nDims
       do i=1,this%nDims
          write(6,'(8x,a,I5,I10,2x,a)') 'rank, ends, name=',i,this%ends(i),trim(this%dimname(i))
       enddo
    endif
    !
  end subroutine get_var_nc_real

  !> Read in 1D integer variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[out] field output variable
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine get_var_nc_int_1d(this,varname,nd1,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1              !  size of array dval
    integer, intent(out) :: field(nd1)     !  values of the field read in
    integer :: ilength
    !
    character*40,parameter :: thissubname='get_var_nc_int_1d'
    !
    integer :: i
    !
    !
    ilength=nd1
    call this%get_var_nc_int(varname,ilength,field)
    !
    if(nd1==this%ends(1)) then
       if(this%debug_level>100) then
          write(6,*) trim(thissubname),' show samples:'
          write(6,*) (field(i),i=1,min(nd1,10))
       endif
    else
       write(6,*) trim(thissubname),' ERROR: dimension does not match.'
    endif
    !
  end subroutine get_var_nc_int_1d

  !> Read in 2D integer type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[in] nd2 length of second dimension
  !! @param[out] field output variable
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine get_var_nc_int_2d(this,varname,nd1,nd2,field)
    !
    ! read in one field 
    !
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1,nd2          !  size of array dval
    integer, intent(out) :: field(nd1,nd2) !  values of the field read in
    integer :: ilength
    !
    integer,allocatable :: temp(:) 
    !
    character*40,parameter :: thissubname='get_var_nc_int_2d'
    !
    integer :: i,j,k
    integer :: istart,iend
    !
    !
    ilength=nd1*nd2
    allocate(temp(ilength))

    call this%get_var_nc_int(varname,ilength,temp)

    if(nd1==this%ends(1) .and. nd2==this%ends(2)) then
       do j=1,nd2
          istart=(j-1)*nd1+1
          iend=(j-1)*nd1+nd1
          field(:,j)=temp(istart:iend)
       enddo
       !
       if(this%debug_level>100) then
          write(6,*) trim(thissubname),' show samples:'
          write(6,*) 'max,min:',maxval(field(:,:)),minval(field(:,:))
       endif
    else
       write(6,*) trim(thissubname),' ERROR: dimension does not match.'
       write(6,*) nd1,this%ends(1),nd2,this%ends(2)
    endif
    deallocate(temp)
    !
  end subroutine get_var_nc_int_2d

  !> Read in 3D integer type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[in] nd2 length of second dimension
  !! @param[in] nd3 length of third dimension
  !! @param[in] field output variable
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine get_var_nc_int_3d(this,varname,nd1,nd2,nd3,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1,nd2,nd3      !  size of array dval
    integer, intent(out) :: field(nd1,nd2,nd3) !  values of the field read in
    integer :: ilength
    !
    integer,allocatable :: temp(:) 
    !
    character*40,parameter :: thissubname='get_var_nc_int_3d'
    !
    integer :: i,j,k
    integer :: length2d
    integer :: istart,iend
    !
    !
    length2d=nd1*nd2
    ilength=length2d*nd3
    allocate(temp(ilength))

    call this%get_var_nc_int(varname,ilength,temp)

    if(nd1==this%ends(1) .and. nd2==this%ends(2) .and. nd3==this%ends(3)) then
       do k=1,nd3
          do j=1,nd2
             istart=(k-1)*length2d+(j-1)*nd1+1
             iend  =(k-1)*length2d+(j-1)*nd1+nd1
             field(:,j,k)=temp(istart:iend)
          enddo
       enddo
       !
       if(this%debug_level>100) then
          write(6,*) trim(thissubname),' show samples:'
          do k=1,nd3
             write(6,*) 'k,max,min:',k,maxval(field(:,:,k)),minval(field(:,:,k))
          enddo
       endif
    else
       write(6,*) trim(thissubname),' ERROR: dimension does not match.'
       write(6,*) nd1,this%ends(1),nd2,this%ends(2),nd3,this%ends(3)
    endif
    deallocate(temp)
    !
  end subroutine get_var_nc_int_3d

  !> Read in integer type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] ilength size of array
  !! @param[in] field output variable
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine get_var_nc_int(this,varname,ilength,field)
    !
    ! read in one field 
    !
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: ilength          !  size of array dval
    integer, intent(out) :: field(ilength)   !  values of the field read in
    !
    integer :: ncid
    ! 
    integer :: status
    integer :: varid
    integer :: ends(4),start(4)

    integer :: length4d,length3d,length2d
    integer :: nDims,ndim
    integer :: dimids(4)
    integer :: xtype
    character*40 :: dimname

    character*40,parameter :: thissubname='get_var_nc_int'
    !
    integer :: i,k
    !
    !
    ncid=this%ncid

    ! get variable IDs
    status = nf90_inq_varid(ncid, trim(varname), VarId)
    if(status /= nf90_NoErr) call this%handle_err(status)

    !  get dimensions
    ends=1
    start=1
    this%ends=1

    this%dimname="                           "
    ! get variable type
    status = nf90_inquire_variable(ncid, VarId, xtype = xtype)
    if(status /= nf90_NoErr) call this%handle_err(status)
    if(xtype==NF90_INT) then
       this%xtype=xtype
    else
       write(6,*) trim(thissubname),' ERROR: wrong data type, expect ',NF90_INT,' but read in ',xtype
       stop 123
    endif

    ! get dimension size
    status = nf90_inquire_variable(ncid, VarId, ndims = nDims)
    if(status /= nf90_NoErr) call this%handle_err(status)
    this%ndims=nDims
    !
    status = nf90_inquire_variable(ncid, VarId, dimids = dimids(1:nDims))
    if(status /= nf90_NoErr) call this%handle_err(status)
    do i=1,nDims
       dimname="       "
       status = nf90_inquire_dimension(ncid, dimids(i), dimname, len = ndim)
       if (status /= nf90_noerr) call this%handle_err(status)
       ends(i)=ndim
       this%ends(i)=ends(i)
       this%dimname(i)=trim(dimname)
       if(this%ends(i) < 1) then
          write(6,*) trim(thissubname),' Error, ends dimension should larger than 0 :', ends(i)
          stop 1234
       endif
    enddo
    length2d=ends(1)*ends(2)
    length3d=length2d*ends(3)
    length4d=length3d*ends(4)
    if(ilength .ne. length4d) then
       write(6,*) trim(thissubname),'ERROR: ',ilength,' should equal to ',length4d
       stop 123
    endif
    !
    if(nDims <=4 ) then
       status = nf90_get_var(ncid, VarId, field, &
            start = start(1:4) , &
            count = ends(1:4))
       if(status /= nf90_NoErr) call this%handle_err(status)
    else
       write(6,*) trim(thissubname),'Error: too many dimensions:',nDims
       stop 1234
    endif
    !
    if(this%debug_level>0) then
       write(6,'(a,a)') '>>>read in variable: ',trim(varname)
    endif
    if(this%debug_level>10) then
       write(6,'(8x,a,I10)') 'data type : ',this%xtype
       write(6,'(8x,a,I10)') 'dimension size: ',this%nDims
       do i=1,this%nDims
          write(6,'(8x,a,I5,I10,2x,a)') 'rank, ends, name=',i,this%ends(i),trim(this%dimname(i))
       enddo
    endif
    !
  end subroutine get_var_nc_int

  !> Read in 1D short type variable
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[out] field output variable
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine get_var_nc_short_1d(this,varname,nd1,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1              !  size of array dval
    integer(2), intent(out) :: field(nd1)     !  values of the field read in
    integer :: ilength
    !
    character*40,parameter :: thissubname='get_var_nc_short_1d'
    !
    integer :: i
    !
    !
    ilength=nd1
    call this%get_var_nc_short(varname,ilength,field)
    !
    if(nd1==this%ends(1)) then
       if(this%debug_level>100) then
          write(6,*) trim(thissubname),' show samples:'
          write(6,*) (field(i),i=1,min(nd1,10))
       endif
    else
       write(6,*) trim(thissubname),' ERROR: dimension does not match.'
    endif
    !
  end subroutine get_var_nc_short_1d

  !> Read in 2D short type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[in] nd2 length of second dimension
  !! @param[out] field output variable
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine get_var_nc_short_2d(this,varname,nd1,nd2,field)
    !
    ! read in one field 
    !
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1,nd2          !  size of array dval
    integer(2), intent(out) :: field(nd1,nd2) !  values of the field read in
    integer :: ilength
    !
    integer(2),allocatable :: temp(:)
    !
    character*40,parameter :: thissubname='get_var_nc_short_2d'
    !
    integer :: i,j,k
    integer :: istart,iend
    !
    !
    ilength=nd1*nd2
    allocate(temp(ilength))

    call this%get_var_nc_short(varname,ilength,temp)

    if(nd1==this%ends(1) .and. nd2==this%ends(2)) then
       do j=1,nd2
          istart=(j-1)*nd1+1
          iend=(j-1)*nd1+nd1
          field(:,j)=temp(istart:iend)
       enddo
       !
       if(this%debug_level>100) then
          write(6,*) trim(thissubname),' show samples:'
          write(6,*) 'max,min:',maxval(field(:,:)),minval(field(:,:))
       endif
    else
       write(6,*) trim(thissubname),' ERROR: dimension does not match.'
       write(6,*) nd1,this%ends(1),nd2,this%ends(2)
    endif
    deallocate(temp)
    !
  end subroutine get_var_nc_short_2d
  !
  !> Read in short type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] ilength size of array
  !! @param[out] field output variable
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine get_var_nc_short(this,varname,ilength,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: ilength          !  size of array dval
    integer(2), intent(out) :: field(ilength)   !  values of the field read in
    !
    integer :: ncid
    ! 
    integer :: status
    integer :: varid
    integer :: ends(4),start(4)

    integer :: length4d,length3d,length2d
    integer :: nDims,ndim
    integer :: dimids(4)
    integer :: xtype
    character*40 :: dimname

    character*40,parameter :: thissubname='get_var_nc_short'
    !
    integer :: i,k
    !
    !
    ncid=this%ncid

    ! get variable IDs
    status = nf90_inq_varid(ncid, trim(varname), VarId)
    if(status /= nf90_NoErr) call this%handle_err(status)

    !  get dimensions
    ends=1
    start=1
    this%ends=1

    this%dimname="                           "
    ! get variable type
    status = nf90_inquire_variable(ncid, VarId, xtype = xtype)
    if(status /= nf90_NoErr) call this%handle_err(status)
    if(xtype==NF90_SHORT) then
       this%xtype=xtype
    else
       write(6,*) trim(thissubname),' ERROR: wrong data type, expect ',NF90_SHORT,' but read in ',xtype
       stop 123
    endif

    ! get dimension size
    status = nf90_inquire_variable(ncid, VarId, ndims = nDims)
    if(status /= nf90_NoErr) call this%handle_err(status)
    this%ndims=nDims
    !
    status = nf90_inquire_variable(ncid, VarId, dimids = dimids(1:nDims))
    if(status /= nf90_NoErr) call this%handle_err(status)
    do i=1,nDims
       dimname="       "
       status = nf90_inquire_dimension(ncid, dimids(i), dimname, len = ndim)
       if (status /= nf90_noerr) call this%handle_err(status)
       ends(i)=ndim
       this%ends(i)=ends(i)
       this%dimname(i)=trim(dimname)
       if(this%ends(i) < 1) then
          write(6,*) trim(thissubname),' Error, ends dimension should larger than 0 :', ends(i)
          stop 1234
       endif
    enddo
    length2d=ends(1)*ends(2)
    length3d=length2d*ends(3)
    length4d=length3d*ends(4)
    if(ilength .ne. length4d) then
       write(6,*) trim(thissubname),'ERROR: ',ilength,' should equal to ',length4d
       stop 123
    endif
    !
    if(nDims <=4 ) then
       status = nf90_get_var(ncid, VarId, field, &
            start = start(1:4) , &
            count = ends(1:4))
       if(status /= nf90_NoErr) call this%handle_err(status)
    else
       write(6,*) trim(thissubname),'Error: too many dimensions:',nDims
       stop 1234
    endif
    !
    if(this%debug_level>0) then
       write(6,'(a,a)') '>>>read in variable: ',trim(varname)
    endif
    if(this%debug_level>10) then
       write(6,'(8x,a,I10)') 'data type : ',this%xtype
       write(6,'(8x,a,I10)') 'dimension size: ',this%nDims
       do i=1,this%nDims
          write(6,'(8x,a,I5,I10,2x,a)') 'rank, ends, name=',i,this%ends(i),trim(this%dimname(i))
       enddo
    endif
    !
  end subroutine get_var_nc_short

  !> Read in 1D character type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[out] field output variable
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine get_var_nc_char_1d(this,varname,nd1,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1              !  size of array dval
    character, intent(out) :: field(nd1)     !  values of the field read in
    integer :: ilength
    !
    character*40,parameter :: thissubname='get_var_nc_char_1d'
    !
    integer :: i
    !
    !
    ilength=nd1
    call this%get_var_nc_char(varname,ilength,field)
    !
    if(nd1==this%ends(1)) then
       if(this%debug_level>100) then
          write(6,*) trim(thissubname),' show samples:'
          write(6,*) (field(i),i=1,min(nd1,10))
       endif
    else
       write(6,*) trim(thissubname),' ERROR: dimension does not match.'
    endif
    !
  end subroutine get_var_nc_char_1d

  !> Read in 2D character type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[in] nd2 length of second dimension
  !! @param[out] field output variable
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine get_var_nc_char_2d(this,varname,nd1,nd2,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1,nd2          !  size of array dval
    character, intent(out) :: field(nd1,nd2) !  values of the field read in
    integer :: ilength
    !
    character,allocatable :: temp(:) 
    !
    character*40,parameter :: thissubname='get_var_nc_char_2d'
    !
    integer :: i,j,k
    integer :: istart,iend
    !
    !
    ilength=nd1*nd2
    allocate(temp(ilength))

    call this%get_var_nc_char(varname,ilength,temp)

    if(nd1==this%ends(1) .and. nd2==this%ends(2)) then
       do j=1,nd2
          istart=(j-1)*nd1+1
          iend=(j-1)*nd1+nd1
          field(:,j)=temp(istart:iend)
       enddo
       !
!       if(this%debug_level>100) then
!          write(*,*) trim(thissubname),' show samples:'
!          write(*,*) field(1,1)
!       endif
    else
       write(6,*) trim(thissubname),' ERROR: dimension does not match.'
       write(6,*) nd1,this%ends(1),nd2,this%ends(2)
    endif
    deallocate(temp)
    !
  end subroutine get_var_nc_char_2d

  !> Read in 3D character type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] nd1 length of first dimension
  !! @param[in] nd2 length of second dimension
  !! @param[in] nd3 length of third dimension
  !! @param[out] field output variable
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine get_var_nc_char_3d(this,varname,nd1,nd2,nd3,field)
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: nd1,nd2,nd3      !  size of array dval
    character, intent(out) :: field(nd1,nd2,nd3) !  values of the field read in
    integer :: ilength
    !
    character,allocatable :: temp(:) 
    !
    character*40,parameter :: thissubname='get_var_nc_char_3d'
    !
    integer :: i,j,k
    integer :: length2d
    integer :: istart,iend
    !
    !
    length2d=nd1*nd2
    ilength=length2d*nd3
    allocate(temp(ilength))

    call this%get_var_nc_char(varname,ilength,temp)

    if(nd1==this%ends(1) .and. nd2==this%ends(2) .and. nd3==this%ends(3)) then
       do k=1,nd3
          do j=1,nd2
             istart=(k-1)*length2d+(j-1)*nd1+1
             iend  =(k-1)*length2d+(j-1)*nd1+nd1
             field(:,j,k)=temp(istart:iend)
          enddo
       enddo
       !
!       if(this%debug_level>100) then
!          write(*,*) trim(thissubname),' show samples:'
!          write(*,*) field(1,1,1)
!       endif
    else
       write(6,*) trim(thissubname),' ERROR: dimension does not match.'
       write(6,*) nd1,this%ends(1),nd2,this%ends(2),nd3,this%ends(3)
    endif
    deallocate(temp)
    !
  end subroutine get_var_nc_char_3d
  !
  !> Read in character type variable.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] varname name of the variable
  !! @param[in] ilength size of array
  !! @param[out] field output variable
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine get_var_nc_char(this,varname,ilength,field)
    !
    ! read in one field 
    !
    use netcdf
    !
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname  ! name of the field to read
    integer, intent(in) :: ilength          !  size of array dval
    character, intent(out) :: field(ilength)   !  values of the field read in
    !
    integer :: ncid
    ! 
    integer :: status
    integer :: varid 
    integer :: ends(4),start(4)

    integer :: length4d,length3d,length2d
    integer :: nDims,ndim
    integer :: dimids(4)
    integer :: xtype
    character*40 :: dimname

    character*40,parameter :: thissubname='get_var_nc_char'
    !
    integer :: i,k
    !
    !
    ncid=this%ncid

    ! get variable IDs
    status = nf90_inq_varid(ncid, trim(varname), VarId)
    if(status /= nf90_NoErr) call this%handle_err(status)

    !  get dimensions
    ends=1
    start=1
    this%ends=1

    this%dimname="                           "
    ! get variable type
    status = nf90_inquire_variable(ncid, VarId, xtype = xtype)
    if(status /= nf90_NoErr) call this%handle_err(status)
    if(xtype==NF90_CHAR) then
       this%xtype=xtype
    else
       write(6,*) trim(thissubname),' ERROR: wrong data type, expect ',NF90_CHAR,' but read in ',xtype
       stop 123
    endif

    ! get dimension size
    status = nf90_inquire_variable(ncid, VarId, ndims = nDims)
    if(status /= nf90_NoErr) call this%handle_err(status)
    this%ndims=nDims
    !
    status = nf90_inquire_variable(ncid, VarId, dimids = dimids(1:nDims))
    if(status /= nf90_NoErr) call this%handle_err(status)
    do i=1,nDims
       dimname="       "
       status = nf90_inquire_dimension(ncid, dimids(i), dimname, len = ndim)
       if (status /= nf90_noerr) call this%handle_err(status)
       ends(i)=ndim
       this%ends(i)=ends(i)
       this%dimname(i)=trim(dimname)
       if(this%ends(i) < 1) then
          write(6,*) trim(thissubname),' Error, ends dimension should larger than 0 :', ends(i)
          stop 1234
       endif
    enddo
    length2d=ends(1)*ends(2)
    length3d=length2d*ends(3)
    length4d=length3d*ends(4)
    if(ilength .ne. length4d) then
       write(6,*) trim(thissubname),'ERROR: ',ilength,' should equal to ',length4d
       stop 123
    endif
    !
    if(nDims <=4 ) then
       status = nf90_get_var(ncid, VarId, field, &
            start = start(1:4) , &
            count = ends(1:4))
       if(status /= nf90_NoErr) call this%handle_err(status)
    else
       write(6,*) trim(thissubname),'Error: too many dimensions:',nDims
       stop 1234
    endif
    !
    if(this%debug_level>0) then
       write(6,'(a,a)') '>>>read in variable: ',trim(varname)
    endif
    if(this%debug_level>10) then
       write(6,'(8x,a,I10)') 'data type : ',this%xtype
       write(6,'(8x,a,I10)') 'dimension size: ',this%nDims
       do i=1,this%nDims
          write(6,'(8x,a,I5,I10,2x,a)') 'rank, ends, name=',i,this%ends(i),trim(this%dimname(i))
       enddo
    endif
    !
  end subroutine get_var_nc_char

  !> Handle netCDF errors.
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] status return code from neCDF
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine handle_err(this,status)
    use netcdf
    implicit none
    class(ncio) :: this
    !
    integer, intent ( in) :: status
    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop "Stopped"
    end if
  end subroutine handle_err

  !> Convert theta T (Kelvin) to T (deg C).
  !!
  !! @param[in] this instance of an ncio class
  !! @param[in] nx number of grid points in x-dir
  !! @param[in] ny number of grid points in y-dir
  !! @param[in] ps Pressure (Pa)
  !! @param[inout] t2 Pot. Temperature (Kelvin)
  !! @author Ming Hu org: GSD/AMB @date 2017-11-01
  subroutine convert_theta2t_2dgrid(this,nx,ny,ps,t2)
    implicit none
    class(ncio) :: this

    integer :: nx,ny
    real, intent(in   ) :: ps(nx,ny)
    real, intent(inout) :: t2(nx,ny)

    integer :: i,j
    real(8) :: rd,cp,rd_over_cp


    rd     = 2.8705e+2_8
    cp     = 1.0046e+3_8  !  specific heat of air @pressure (J/kg/K)
    rd_over_cp = rd/cp

    do j=1,ny
       do i=1,nx
          t2(i,j)=t2(i,j)*(ps(i,j)/1000.0)**rd_over_cp - 273.15
       enddo
    enddo

  end subroutine convert_theta2t_2dgrid

  !> Add a new variable to sfc_data.nc with dimensions (Time, yaxis_1,
  !! xaxis_1).
  !!
  !! @param this instance of an ncio class
  !! @param[in] varname Name of variable to be created in netcdf file
  !! @param[in] dname1 1st dimension name
  !! @param[in] dname2 2nd dimension name
  !! @param[in] dname3 3rd dimension name
  !! @param[in] lname long name output for netcdf variable
  !! @param[in] units units to use in netcdf variable
  !!
  !! @author David.M.Wright org: UM/GLERL @date 2020-09-01
  subroutine add_new_var_3d(this,varname,dname1,dname2,dname3,lname,units,dtype)
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname,dname1,dname2,dname3 &
         ,lname,units
    integer :: status, ncid, dim1id, dim2id, dim3id, varid
    character(len=*),intent(in) :: dtype

    status = nf90_redef(this%ncid) !Enter Define Mode
    if (status /= nf90_noerr) call this%handle_err(status)

    status = nf90_inq_dimid(this%ncid, dname1, dim1id)
    if (status /= nf90_noerr) call this%handle_err(status)
    status = nf90_inq_dimid(this%ncid, dname2, dim2id)
    if (status /= nf90_noerr) call this%handle_err(status)
    status = nf90_inq_dimid(this%ncid, dname3, dim3id)
    if (status /= nf90_noerr) call this%handle_err(status)

    if(trim(dtype)=="double") then
        status = nf90_def_var(this%ncid, varname, nf90_double, &
         (/ dim1id, dim2id, dim3id /), varid)
    elseif(trim(dtype)=="float") then
        status = nf90_def_var(this%ncid, varname, nf90_float, &
         (/ dim1id, dim2id, dim3id /), varid)
    elseif(trim(dtype)=="int") then
        status = nf90_def_var(this%ncid, varname, nf90_int, &
         (/ dim1id, dim2id, dim3id /), varid)
    else 
      write(*,*) ' undefined data type ', trim(dtype)
      call this%handle_err(status)
    endif
    if (status /= nf90_noerr) call this%handle_err(status)

    status = nf90_put_att(this%ncid, varid, 'long_name', lname)
    if (status /= nf90_noerr) call this%handle_err(status)
    status = nf90_put_att(this%ncid, varid, 'units', units)
    if (status /= nf90_noerr) call this%handle_err(status)

    status = nf90_enddef(this%ncid) !Exit Define Mode and
    ! return to Data Mode
    if (status /= nf90_noerr) call this%handle_err(status)

  end subroutine add_new_var_3d

  !> Add a new variable to sfc_data.nc with dimensions (yaxis_1,
  !! xaxis_1).
  !!
  !! @param this instance of an ncio class
  !! @param[in] varname Name of variable to be created in netcdf file
  !! @param[in] dname1 1st dimension name
  !! @param[in] dname2 2nd dimension name
  !! @param[in] lname long name output for netcdf variable
  !! @param[in] units units to use in netcdf variable
  !!
  !! @author David.M.Wright org: UM/GLERL @date 2021-10-07
  subroutine add_new_var_2d(this,varname,dname1,dname2,lname,units,dtype)
    implicit none
    !
    class(ncio) :: this
    character(len=*),intent(in) :: varname,dname1,dname2  &
         ,lname,units
    integer :: status, ncid, dim1id, dim2id, varid
    character(len=*),intent(in) :: dtype

    status = nf90_redef(this%ncid) !Enter Define Mode
    if (status /= nf90_noerr) call this%handle_err(status)

    status = nf90_inq_dimid(this%ncid, dname1, dim1id)
    if (status /= nf90_noerr) call this%handle_err(status)
    status = nf90_inq_dimid(this%ncid, dname2, dim2id)
    if (status /= nf90_noerr) call this%handle_err(status)

    if(trim(dtype)=="double") then
      status = nf90_def_var(this%ncid, varname, nf90_double, &
         (/ dim1id, dim2id /), varid)
    elseif(trim(dtype)=="float") then
      status = nf90_def_var(this%ncid, varname, nf90_float, &
         (/ dim1id, dim2id /), varid)
    elseif(trim(dtype)=="int") then
      status = nf90_def_var(this%ncid, varname, nf90_int, &
         (/ dim1id, dim2id /), varid)
    else
      write(*,*) ' undefined data type ', trim(dtype)
      call this%handle_err(status)
    endif
    if (status /= nf90_noerr) call this%handle_err(status)

    status = nf90_put_att(this%ncid, varid, 'long_name', lname)
    if (status /= nf90_noerr) call this%handle_err(status)
    status = nf90_put_att(this%ncid, varid, 'units', units)
    if (status /= nf90_noerr) call this%handle_err(status)

    status = nf90_enddef(this%ncid) !Exit Define Mode and
    ! return to Data Mode
    if (status /= nf90_noerr) call this%handle_err(status)

  end subroutine add_new_var_2d


end module module_ncio
