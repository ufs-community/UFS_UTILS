!> @file
!! @brief Generate fixed grid files required for coupled model
!!
!! @author Denise.Worthen@noaa.gov

!> Generate fixed grid files required for coupled model using the MOM6 super grid file and ocean mask file. It creates
!! a main grid file which is then used to create subsequent files which are required to create the fix and IC
!! files required for the S2S or S2SW application.
!!
!! This executable created with this source code runs within the shell scrip cpld_gridgen.sh in ../../ush, which
!! utilizes both NCO (netCDF Operators) and ESMF command line functions. The shell script creates a run-time grid.nml
!! from grid.nml.IN
!!
!! @author Denise.Worthen@noaa.gov
!! @return 0 for success, error code otherwise.
program gen_fixgrid

  use ESMF
  use mpi_f08

  use grdvars,           only: ni,nj,nv,deg2rad,sg_maxlat,iVertNE,jVertNE,maximum_lat
  use grdvars,           only: grid
  use gengrid_utils,     only: allocate_all
  use inputnml
  use gengrid_kinds,     only: CL, CS, dbl_kind, real_kind, int_kind
  use angles,            only: find_ang, find_angq, find_angchk
  use vertices,          only: fill_vertices
  use mapped_mask,       only: make_frac_land
  use postwgts,          only: make_postwgts
  use tripolegrid,       only: write_tripolegrid
  use cicegrid,          only: write_cicegrid
  use scripgrid,         only: write_staggers
  use topoedits,         only: add_topoedits, apply_topoedits
  use charstrings,       only: logmsg, res, atmres, dirsrc, dirout, fv3dir, editsfile
  use charstrings,       only: maskfile, maskname, topofile, toponame, editsfile, staggerlocs, cdate, history
  use netcdf

  implicit none

  ! local variables
  type(MPI_Comm) :: mpic  ! mpi_f08
  real(dbl_kind) :: dxT, dyT

  integer(int_kind)  :: ipole(2)

  real(real_kind),   allocatable, dimension(:,:) :: ww3dpth
  integer(int_kind), allocatable, dimension(:,:) :: ww3mask

  character(len=CL) :: fsrc, fdst, fwgt
  character(len= 2) :: cstagger

  ! Super-grid source grid variables
  real(dbl_kind), allocatable, dimension(:,:)   :: x      !Supergrid lon
  real(dbl_kind), allocatable, dimension(:,:)   :: y      !Supergrid lat
  real(dbl_kind), allocatable, dimension(:,:)   :: dx     !Supergrid cell width (m)
  real(dbl_kind), allocatable, dimension(:,:)   :: dy     !Supergrid cell width (m)

  integer :: int_mpic
  integer :: rc,ncid,id,xtype
  integer :: i,j,k,n,i2,j2,nvalid
  integer :: ii
  integer :: ierr
  integer :: localPet, nPet
  logical :: fexist = .false.
  logical :: maintask

  type(ESMF_RegridMethod_Flag) :: method
  type(ESMF_VM) :: vm

  !WW3 file format for mod_def generation
  character(len= 6) :: i4fmt = '(i4.4)'
  character(len=CS) :: form1
  character(len=CS) :: form2
  character(len= 6) :: cnx

  !-------------------------------------------------------------------------
  ! Initialize esmf environment. Everything except the generation of the
  ! ESMF weights is done on the root PE.
  !-------------------------------------------------------------------------

  call ESMF_Initialize()
  call ESMF_VMGetGlobal(vm)
  call ESMF_VMGet(vm, localPet=localPet, peCount=nPet, mpiCommunicator=int_mpic, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  mpic%mpi_val = int_mpic
  maintask = .false.
  if (localPet == 0) maintask=.true.
  if (maintask) then
     print '(a,i4,a)','Running on = ',npet,' tasks'
     !---------------------------------------------------------------------
     !
     !---------------------------------------------------------------------

     call read_inputnml('grid.nml')

     print '(a,2i6)',' output grid requested ',ni,nj
     print '(a,2i6)',' supergrid size used ', nx,ny
     print '(a)',' output grid tag '//trim(res)
     print '(a)',' supergrid source directory '//trim(dirsrc)
     print '(a)',' output grid directory '//trim(dirout)
     print '(a)',' atm mosaic directory '//trim(fv3dir)
     print '(a)',' MOM6 topography file '//trim(topofile)
     print '(a)',' MOM6 edits file '//trim(editsfile)
     print *,'editmask flag ',editmask
     print *,'debug flag ',debug
     print *,'do_postwgts flag ',do_postwgts
     print *

     call allocate_all(ni,nj,grid)

     call ESMF_LogWrite("Starting gen_fixgrid", ESMF_LOGMSG_INFO)
     !---------------------------------------------------------------------
     ! set up the arrays to retrieve the vertices
     !---------------------------------------------------------------------

     grid%Ct%iVert = iVertNE
     grid%Ct%jVert = jVertNE
     grid%Cu%iVert = grid%Ct%iVert + 1; grid%Cu%jVert = grid%Ct%jvert + 0
     grid%Cv%iVert = grid%Ct%iVert + 0; grid%Cv%jVert = grid%Ct%jVert + 1
     grid%Bu%iVert = grid%Ct%iVert + 1; grid%Bu%jVert = grid%Ct%jVert + 1

     print '(a8,4i6)','iVertCt ',(grid%Ct%iVert(i),i=1,4)
     print '(a8,4i6)','jVertCt ',(grid%Ct%jVert(i),i=1,4)
     print *
     print '(a8,4i6)','iVertCu ',(grid%Cu%iVert(i),i=1,4)
     print '(a8,4i6)','jVertCu ',(grid%Cu%jVert(i),i=1,4)
     print *
     print '(a8,4i6)','iVertCv ',(grid%Cv%iVert(i),i=1,4)
     print '(a8,4i6)','jVertCv ',(grid%Cv%jVert(i),i=1,4)
     print *
     print '(a8,4i6)','iVertBu ',(grid%Bu%iVert(i),i=1,4)
     print '(a8,4i6)','jVertBu ',(grid%Bu%jVert(i),i=1,4)
     print *

     grid%Ct%latvert = -9999.0; grid%Ct%lonvert = -9999.0
     grid%Cu%latvert = -9999.0; grid%Cu%lonvert = -9999.0
     grid%Cv%latvert = -9999.0; grid%Cv%lonvert = -9999.0
     grid%Bu%latvert = -9999.0; grid%Bu%lonvert = -9999.0

     !---------------------------------------------------------------------
     ! read the MOM6 land mask
     !---------------------------------------------------------------------

     fsrc = trim(dirsrc)//trim(maskfile)

     rc = nf90_open(fsrc, nf90_nowrite, ncid)
     print '(a)', 'reading ocean mask from '//trim(fsrc)
     if(rc .ne. 0)print '(a)', 'nf90_open = '//trim(nf90_strerror(rc))

     grid%wet4 = 0.0; grid%wet8 = 0.0
     rc = nf90_inq_varid(ncid,  trim(maskname), id)
     rc = nf90_inquire_variable(ncid, id, xtype=xtype)
     if(xtype .eq. 5)rc = nf90_get_var(ncid,      id,  grid%wet4)
     if(xtype .eq. 6)rc = nf90_get_var(ncid,      id,  grid%wet8)
     rc = nf90_close(ncid)

     if(xtype.eq. 6)grid%wet4 = real(grid%wet8,4)

     !---------------------------------------------------------------------
     ! read the MOM6 depth file
     !---------------------------------------------------------------------

     fsrc = trim(dirsrc)//trim(topofile)

     rc = nf90_open(fsrc, nf90_nowrite, ncid)
     print '(a)', 'reading ocean topography from '//trim(fsrc)
     if(rc .ne. 0)print '(a)', 'nf90_open = '//trim(nf90_strerror(rc))

     grid%dp4 = 0.0; grid%dp8 = 0.0
     rc = nf90_inq_varid(ncid,  trim(toponame), id)
     rc = nf90_inquire_variable(ncid, id, xtype=xtype)
     if(xtype .eq. 5)rc = nf90_get_var(ncid,      id,  grid%dp4)
     if(xtype .eq. 6)rc = nf90_get_var(ncid,      id,  grid%dp8)
     rc = nf90_close(ncid)

     if(xtype.eq. 6)grid%dp4 = real(grid%dp8,4)

     if(editmask)then
        !---------------------------------------------------------------------
        ! apply topoedits run time mask changes if required for this config
        ! this will create a modified topoedits file which accounts for any
        ! land mask changes created at run time by MOM6
        !---------------------------------------------------------------------

        if(trim(editsfile)  == 'none')then
           print '(a)', 'Need a valid editsfile to make mask edits '
           call abort()
        end if
        inquire(file=trim(dirsrc)//trim(editsfile),exist=fexist)
        if (.not. fexist) then
           print '(a)', 'Required topoedits file '//trim(editsfile) &
                //'for land mask changes is missing '
           call abort()
        end if

        fsrc = trim(dirsrc)//trim(editsfile)
        fdst = trim(dirout)//'ufs.'//trim(editsfile)
        call add_topoedits(fsrc,fdst,grid%wet4)
     endif

     !---------------------------------------------------------------------
     ! MOM6 reads the depth file, applies the topo edits and then adjusts
     ! depth using masking_depth and min/max depth. This call mimics
     ! MOM6 routines apply_topography_edits_from_file and limit_topography
     ! If the the topoedits file has been modified to account for MOM6 run
     ! time land mask changes (above), then the depth will be created using
     ! this modified topoedits file
     !---------------------------------------------------------------------

     fsrc = trim(dirsrc)//trim(editsfile)
     if(editmask)fsrc = trim(dirout)//'ufs.'//trim(editsfile)

     if (trim(editsfile) /= 'none') then
        inquire(file=trim(fsrc),exist=fexist)
        if (.not. fexist) then
           print '(a)', 'Required topoedits file '//trim(fsrc)//' is missing '
           call abort()
        end if
     end if
     call apply_topoedits(fsrc,grid%dp4)

     !---------------------------------------------------------------------
     ! read MOM6 supergrid file
     !---------------------------------------------------------------------

     allocate( x(0:nx,0:ny),  y(0:nx,0:ny) )
     allocate(  dx(nx,0:ny), dy(0:nx,ny) )

     fsrc = trim(dirsrc)//'ocean_hgrid.nc'

     rc = nf90_open(fsrc, nf90_nowrite, ncid)
     print '(a)', 'reading supergrid from '//trim(fsrc)
     if(rc .ne. 0)print '(a)', 'nf90_open = '//trim(nf90_strerror(rc))

     rc = nf90_inq_varid(ncid, 'x', id)  !lon
     rc = nf90_get_var(ncid,    id,  x)

     rc = nf90_inq_varid(ncid, 'y', id)  !lat
     rc = nf90_get_var(ncid,    id,  y)

     rc = nf90_inq_varid(ncid, 'dx', id)
     rc = nf90_get_var(ncid,     id, dx)

     rc = nf90_inq_varid(ncid, 'dy', id)
     rc = nf90_get_var(ncid,     id, dy)

     rc = nf90_close(ncid)
     sg_maxlat = maxval(y)
     write(logmsg,'(a,f12.2)')'max lat in super grid ',maxval(y)
     print '(a)',trim(logmsg)

     !---------------------------------------------------------------------
     ! fill grid variables
     !---------------------------------------------------------------------

     do j = 1,nj
        do i = 1,ni
           i2 = 2*i ; j2 = 2*j
           !deg->rad
           grid%ulon(i,j) =     x(i2,j2)*deg2rad
           grid%ulat(i,j) =     y(i2,j2)*deg2rad
           !m->cm
           grid%htn(i,j) = (dx(i2-1,j2) + dx(i2,j2))*100._dbl_kind
           grid%hte(i,j) = (dy(i2,j2-1) + dy(i2,j2))*100._dbl_kind
           !deg
           grid%Bu%lon(i,j) =     x(i2,j2)
           grid%Bu%lat(i,j) =     y(i2,j2)
           !deg
           grid%Ct%lon(i,j) =     x(i2-1,j2-1)
           grid%Cu%lon(i,j) =     x(i2,  j2-1)
           grid%Cv%lon(i,j) =     x(i2-1,j2  )
           !deg
           grid%Ct%lat(i,j) =     y(i2-1,j2-1)
           grid%Cu%lat(i,j) =     y(i2,  j2-1)
           grid%Cv%lat(i,j) =     y(i2-1,j2  )
           !m2
           dxT = dx(i2-1,j2-1) + dx(i2,j2-1)
           dyT = dy(i2-1,j2-1) + dy(i2-1,j2)
           grid%areaCt(i,j) = dxT*dyT
        enddo
     enddo

     !---------------------------------------------------------------------
     ! locate the ith index of the two poles on j=nj
     ! the corner points must lie on the pole
     !---------------------------------------------------------------------

     ipole = -1
     j = nj
     do i = 1,ni/2
        if(grid%Bu%lat(i,j) .eq. sg_maxlat)ipole(1) = i
     enddo
     do i = ni/2+1,ni
        if(grid%Bu%lat(i,j) .eq. sg_maxlat)ipole(2) = i
     enddo
     write(logmsg,'(a,2i6,2f12.2)')'poles found at i = ',ipole, &
          grid%Bu%lat(ipole(1),nj), grid%Bu%lat(ipole(2),nj)
     print '(a)',trim(logmsg)

     !---------------------------------------------------------------------
     ! find the angle on centers using the same procedure as MOM6
     !---------------------------------------------------------------------

     call find_ang((/1,ni/),(/1,nj/),grid%Bu%lon,grid%Bu%lat,grid%Ct%lon,grid%anglet)
     write(logmsg,'(a,2f12.2)')'ANGLET min,max: ',minval(grid%anglet),maxval(grid%anglet)
     print '(a)',trim(logmsg)
     write(logmsg,'(a,2f12.2)')'ANGLET edges i=1,i=ni: ',grid%anglet(1,nj),grid%anglet(ni,nj)
     print '(a)',trim(logmsg)

     grid%xangCt(:) = 0.0
     do i = 1,ni
        i2 = ipole(2)+(ipole(1)-i)+1
        grid%xangCt(i) = -grid%anglet(i2,nj)       ! angle changes sign across seam
     end do

     !---------------------------------------------------------------------
     ! find the angle on corners using the same procedure as CICE6
     !---------------------------------------------------------------------

     call find_angq((/1,ni/),(/1,nj/),grid%xangCt,grid%anglet,grid%angle)
     grid%angle(ni,:) = -grid%angle(1,:)
     ! reverse angle for CICE
     grid%angle = -grid%angle
     write(logmsg,'(a,2f12.2)')'ANGLE min,max: ',minval(grid%angle),maxval(grid%angle)
     print '(a)',trim(logmsg)
     write(logmsg,'(a,2f12.2)')'ANGLE edges i=1,i=ni: ',grid%angle(1,nj),grid%angle(ni,nj)
     print '(a)',trim(logmsg)

     !---------------------------------------------------------------------
     ! check the Bu angle
     !---------------------------------------------------------------------

     call find_angchk((/1,ni/),(/1,nj/),grid%angle,grid%angchk)
     grid%angchk(1,:) = -grid%angchk(ni,:)
     ! reverse angle for MOM6
     grid%angchk = -grid%angchk
     write(logmsg,'(a,2f12.2)')'ANGCHK min,max: ',minval(grid%angchk),maxval(grid%angchk)
     print '(a)',trim(logmsg)
     write(logmsg,'(a,2f12.2)')'ANGCHK edges i=1,i=ni: ',grid%angchk(1,nj),grid%angchk(ni,nj)
     print '(a)',trim(logmsg)

     !---------------------------------------------------------------------
     ! For the 1/4deg grid, hte at j=720 and j = 1440 is identically=0.0 for
     ! j > 840 (64.0N). These are land points, but since CICE uses hte to
     ! generate remaining variables, setting them to zero will cause problems
     ! For 1deg grid, hte at ni/2 and ni are very small O~10-12, so test for
     ! hte < 1.0
     !---------------------------------------------------------------------

     write(logmsg,'(a,2e12.5)')'min vals of hte at folds ', minval(grid%hte(ni/2,:)),minval(grid%hte(ni,:))
     print '(a)',trim(logmsg)
     do j = 1,nj
        ii = ni/2
        if(grid%hte(ii,j) .le. 1.0)grid%hte(ii,j) = 0.5*(grid%hte(ii-1,j) + grid%hte(ii+1,j))
        ii = ni
        if(grid%hte(ii,j) .le. 1.0)grid%hte(ii,j) = 0.5*(grid%hte(ii-1,j) + grid%hte(   1,j))
     enddo
     write(logmsg,'(a,2e12.5)')'min vals of grid%hte at folds ', minval(grid%hte(ni/2,:)),minval(grid%hte(ni,:))
     print '(a)',trim(logmsg)

     !---------------------------------------------------------------------
     ! find required extended values for setting all vertices
     !---------------------------------------------------------------------

     !if(debug)call checkseam

     do i = 1,ni
        i2 = ipole(2)+(ipole(1)-i)+1
        grid%Ct%xlon(i) = grid%Ct%lon(i2,nj)
        grid%Ct%xlat(i) = grid%Ct%lat(i2,nj)
     enddo

     do i = 1,ni
        i2 = ipole(2)+(ipole(1)-i)
        if(i2 .lt. 1)i2 = ni
        grid%Cu%xlon(i) = grid%Cu%lon(i2,nj)
        grid%Cu%xlat(i) = grid%Cu%lat(i2,nj)
     enddo

     !if(debug)call checkxlatlon

     ! values outside grid(j=0)
     do i = 1,ni
        grid%Bu%xlat(i) = grid%Bu%lat(i,1) + 2.0*(grid%Cu%lat(i,1) - grid%Bu%lat(i,1))
        grid%Cv%xlat(i) = grid%Ct%lat(i,1) + 2.0*(grid%Ct%lat(i,1) - grid%Cv%lat(i,1))
        grid%Bu%xlon(i) = grid%Bu%lon(i,1)
        grid%Cv%xlon(i) = grid%Cv%lon(i,1)
     enddo

     !---------------------------------------------------------------------
     ! fill grid vertices variables
     !---------------------------------------------------------------------

     call fill_vertices(grid%Ct%iVert, grid%Ct%jVert, grid%Bu%lat, grid%Bu%lon, grid%Bu%xlat, grid%Bu%xlon, grid%Ct%latvert, grid%Ct%lonvert, 0)
     call fill_vertices(grid%Cu%iVert, grid%Cu%jVert, grid%Cv%lat, grid%Cv%lon, grid%Cv%xlat, grid%Cv%xlon, grid%Cu%latvert, grid%Cu%lonvert, 0)
     call fill_vertices(grid%Cv%iVert, grid%Cv%jVert, grid%Cu%lat, grid%Cu%lon, grid%Cu%xlat, grid%Cu%xlon, grid%Cv%latvert, grid%Cv%lonvert)
     call fill_vertices(grid%Bu%iVert, grid%Bu%jVert, grid%Ct%lat, grid%Ct%lon, grid%Ct%xlat, grid%Ct%xlon, grid%Bu%latvert, grid%Bu%lonvert)

     !if(debug)call checkpoint

     if(minval(grid%Ct%latvert) .lt. -1.e3)stop
     if(minval(grid%Ct%lonvert) .lt. -1.e3)stop
     if(minval(grid%Cu%latvert) .lt. -1.e3)stop
     if(minval(grid%Cu%lonvert) .lt. -1.e3)stop
     if(minval(grid%Cv%latvert) .lt. -1.e3)stop
     if(minval(grid%Cv%lonvert) .lt. -1.e3)stop
     if(minval(grid%Bu%latvert) .lt. -1.e3)stop
     if(minval(grid%Bu%lonvert) .lt. -1.e3)stop
     deallocate(grid%Ct%xlon, grid%Ct%xlat, grid%Cu%xlon, grid%Cu%xlat, grid%Bu%xlat, grid%Bu%xlon, grid%Cv%xlat, grid%Cv%xlon)

     !---------------------------------------------------------------------
     ! write out grid file files
     !---------------------------------------------------------------------

     ! create a history attribute
     call date_and_time(date=cdate)
     history = 'created on '//trim(cdate)//' from '//trim(fsrc)

     ! write fix grid
     fdst = trim(dirout)//'tripole.mx'//trim(res)//'.nc'
     call write_tripolegrid(trim(fdst),(/1,ni/),(/1,nj/),grid)

     ! write cice grid
     fdst = trim(dirout)//'grid_cice_NEMS_mx'//trim(res)//'.nc'
     call write_cicegrid(trim(fdst),(/1,ni/),(/1,nj/),grid)
     deallocate(grid%ulon, grid%ulat, grid%htn, grid%hte)

     ! write SCRIP files for generation of positional weights
     cstagger = 'Ct'
     fdst = trim(dirout)//trim(cstagger)//'.mx'//trim(res)//'_SCRIP.nc'
     call write_staggers(trim(fdst),(/1,ni/),(/1,nj/),grid%Ct%lon,grid%Ct%lat,grid%Ct%lonvert,grid%Ct%latvert)
     fdst= trim(dirout)//trim(cstagger)//'.mx'//trim(res)//'_SCRIP_land.nc'
     call write_staggers(trim(fdst),(/1,ni/),(/1,nj/),grid%Ct%lon,grid%Ct%lat,grid%Ct%lonvert,grid%Ct%latvert,imask=int(grid%wet4))

     cstagger = 'Cu'
     fdst = trim(dirout)//trim(cstagger)//'.mx'//trim(res)//'_SCRIP.nc'
     call write_staggers(trim(fdst),(/1,ni/),(/1,nj/),grid%Cu%lon,grid%Cu%lat,grid%Cu%lonvert,grid%Cu%latvert)
     cstagger = 'Cv'
     fdst = trim(dirout)//trim(cstagger)//'.mx'//trim(res)//'_SCRIP.nc'
     call write_staggers(trim(fdst),(/1,ni/),(/1,nj/),grid%Cv%lon,grid%Cv%lat,grid%Cv%lonvert,grid%Cv%latvert)
     cstagger = 'Bu'
     fdst = trim(dirout)//trim(cstagger)//'.mx'//trim(res)//'_SCRIP.nc'
     call write_staggers(trim(fdst),(/1,ni/),(/1,nj/),grid%Bu%lon,grid%Bu%lat,grid%Bu%lonvert,grid%Bu%latvert)
     deallocate(grid%Ct%latvert, grid%Ct%lonvert)
     deallocate(grid%Cv%latvert, grid%Cv%lonvert)
     deallocate(grid%Cu%latvert, grid%Cu%lonvert)
     deallocate(grid%Bu%latvert, grid%Bu%lonvert)

     !---------------------------------------------------------------------
     ! write lat,lon,depth and mask arrays required by ww3 in creating
     ! mod_def file
     !---------------------------------------------------------------------

     write(cnx,i4fmt)nx
     write(form1,'(a)')'('//trim(cnx)//'f14.8)'
     write(form2,'(a)')'('//trim(cnx)//'i2)'

     allocate(ww3mask(1:ni,1:nj), source=0); ww3mask = int(grid%wet4)
     allocate(ww3dpth(1:ni,1:nj), source=0.0_real_kind); ww3dpth = grid%dp4

     where(grid%Ct%lat .ge. maximum_lat)ww3mask = 3
     !close last row
     ww3mask(:,nj) = 3

     open(unit=21,file=trim(dirout)//'ww3.mx'//trim(res)//'_x.inp',form='formatted')
     open(unit=22,file=trim(dirout)//'ww3.mx'//trim(res)//'_y.inp',form='formatted')
     open(unit=23,file=trim(dirout)//'ww3.mx'//trim(res)//'_bottom.inp',form='formatted')
     open(unit=24,file=trim(dirout)//'ww3.mx'//trim(res)//'_mapsta.inp',form='formatted')
     ! cice0 .ne. cicen requires obstruction map, should be initialized as zeros (w3grid,ln3032)
     open(unit=25,file=trim(dirout)//'ww3.mx'//trim(res)//'_obstr.inp',form='formatted')

     do j = 1,nj
        write( 21,trim(form1))grid%Ct%lon(:,j)
        write( 22,trim(form1))grid%Ct%lat(:,j)
     end do
     do j = 1,nj
        write( 23,trim(form1))ww3dpth(:,j)
        write( 24,trim(form2))ww3mask(:,j)
        !'obsx' and 'obsy' arrays ???
        write( 25,trim(form2))ww3mask(:,j)*0
        write( 25,trim(form2))ww3mask(:,j)*0
     end do

     close(21); close(22); close(23); close(24); close(25)
     deallocate(ww3mask); deallocate(ww3dpth)
     deallocate(grid%wet4, grid%wet8)

     nvalid = size(catm)
  end if ! if (maintask)
#ifdef test
  !---------------------------------------------------------------------
  ! set up for parallel work
  !---------------------------------------------------------------------
  call mpi_bcast(nvalid, 1, MPI_INTEGER, 0, mpic, ierr)
  if (ierr /= MPI_SUCCESS) then
     print *,' error in mpi broadcast for size(catm) '
     call mpi_abort(mpic, rc)
     call ESMF_Finalize(endflag=ESMF_END_ABORT)
  end if
  if (.not. maintask) then
     allocate(catm(nvalid))
  end if
  call mpi_bcast(catm, size(catm), MPI_INTEGER, 0, mpic, ierr)
  if (ierr /= MPI_SUCCESS) then
     print '(a)',' error in mpi broadcast for catm '
     call mpi_abort(mpic, rc)
     call ESMF_Finalize(endflag=ESMF_END_ABORT)
  end if
  call mpi_bcast(dirout, len(dirout), MPI_CHARACTER, 0, mpic, ierr)
  if (ierr /= MPI_SUCCESS) then
     print '(a)',' error in mpi broadcast for dirout '
     call mpi_abort(mpic, rc)
     call ESMF_Finalize(endflag=ESMF_END_ABORT)
  end if
  call mpi_bcast(res,    len(res),    MPI_CHARACTER, 0, mpic, ierr)
  if (ierr /= MPI_SUCCESS) then
     print '(a)',' error in mpi broadcast for res '
     call mpi_abort(mpic, rc)
     call ESMF_Finalize(endflag=ESMF_END_ABORT)
  end if
  call mpi_bcast(fv3dir, len(fv3dir), MPI_CHARACTER, 0, mpic, ierr)
  if (ierr /= MPI_SUCCESS) then
     print '(a)',' error in mpi broadcast for fv3dir '
     call mpi_abort(mpic, rc)
     call ESMF_Finalize(endflag=ESMF_END_ABORT)
  end if
  call mpi_bcast(do_postwgts, 1, MPI_LOGICAL, 0, mpic, ierr)
  if (ierr /= MPI_SUCCESS) then
     print '(a)',' error in mpi broadcast for do_postwgts '
     call mpi_abort(mpic, rc)
     call ESMF_Finalize(endflag=ESMF_END_ABORT)
  end if
  !---------------------------------------------------------------------
  ! use ESMF regridding to produce mapped ocean mask; first generate
  ! conservative regrid weights from ocean to tiles; then generate the
  ! tiled files containing the mapped ocean mask
  !---------------------------------------------------------------------

  do n = 1,size(catm)
     npx = catm(n)
     if (npx < 100) then
        write(atmres,'(a,i2)')'C',npx
     elseif (npx < 1000) then
        write(atmres,'(a,i3)')'C',npx
     else
        write(atmres,'(a,i4)')'C',npx
     end if

     method=ESMF_REGRIDMETHOD_CONSERVE
     fsrc = trim(dirout)//'Ct.mx'//trim(res)//'_SCRIP_land.nc'
     fdst = trim(fv3dir)//trim(atmres)//'/'//trim(atmres)//'_mosaic.nc'
     fwgt = trim(dirout)//'Ct.mx'//trim(res)//'.to.'//trim(atmres)//'.nc'
     logmsg = 'creating weight file '//trim(fwgt)
     if (maintask) print '(a)',trim(logmsg)

     call ESMF_RegridWeightGen(srcFile=trim(fsrc),dstFile=trim(fdst),         &
          weightFile=trim(fwgt), regridmethod=method,                         &
          unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, ignoreDegenerate=.true., &
          netcdf4fileFlag=.true., tileFilePath=trim(fv3dir)//trim(atmres)//'/', rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  end do

  !---------------------------------------------------------------------
  ! use ESMF to create positional weights for mapping a field from its
  ! native stagger location (Cu,Cv,Bu) onto the center (Ct) grid location
  ! these are used for both post and downscaling
  !---------------------------------------------------------------------

  method=ESMF_REGRIDMETHOD_BILINEAR
  fdst = trim(dirout)//'Ct.mx'//trim(res)//'_SCRIP.nc'
  do k = 2,nv
     cstagger = trim(staggerlocs(k))
     fsrc = trim(dirout)//trim(cstagger)//'.mx'//trim(res)//'_SCRIP.nc'
     fwgt = trim(dirout)//'tripole.mx'//trim(res)//'.'//trim(cstagger)//'.to.Ct.bilinear.nc'
     logmsg = 'creating weight file '//trim(fwgt)
     if (maintask) print '(a)',trim(logmsg)

     call ESMF_RegridWeightGen(srcFile=trim(fsrc),dstFile=trim(fdst), &
          weightFile=trim(fwgt), regridmethod=method,                 &
          ignoreDegenerate=.true., largeFileFlag=.true.,              &
          unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  end do

  !---------------------------------------------------------------------
  ! use ESMF to create positional weights for mapping a field from the
  ! center (Ct) grid location back to the native stagger location
  ! (Cu,Cv,Bu).
  !---------------------------------------------------------------------

  method=ESMF_REGRIDMETHOD_BILINEAR
  fsrc = trim(dirout)//'Ct.mx'//trim(res)//'_SCRIP.nc'
  do k = 2,nv
     cstagger = trim(staggerlocs(k))
     fdst = trim(dirout)//trim(cstagger)//'.mx'//trim(res)//'_SCRIP.nc'
     fwgt = trim(dirout)//'tripole.mx'//trim(res)//'.Ct.to.'//trim(cstagger)//'.bilinear.nc'
     logmsg = 'creating weight file '//trim(fwgt)
     if (maintask) print '(a)',trim(logmsg)

     call ESMF_RegridWeightGen(srcFile=trim(fsrc),dstFile=trim(fdst), &
          weightFile=trim(fwgt), regridmethod=method,                 &
          ignoreDegenerate=.true., largeFileFlag=.true.,              &
          unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  end do

  if(do_postwgts)call make_postwgts(maintask)
  if (maintask) then
     !---------------------------------------------------------------------
     ! make mapped ocean mask file and clean up
     !---------------------------------------------------------------------

     do n = 1,size(catm)
        npx = catm(n)
        if (npx < 100) then
           write(atmres,'(a,i2)')'C',npx
        elseif (npx < 1000) then
           write(atmres,'(a,i3)')'C',npx
        else
           write(atmres,'(a,i4)')'C',npx
        end if
        fsrc = trim(dirout)//'Ct.mx'//trim(res)//'_SCRIP_land.nc'
        fwgt = trim(dirout)//'Ct.mx'//trim(res)//'.to.'//trim(atmres)//'.nc'
        logmsg = 'creating mapped ocean mask for '//trim(atmres)
        print '(a)',trim(logmsg)
        call make_frac_land(trim(fsrc), trim(fwgt))
     end do

     !---------------------------------------------------------------------
     ! clean up
     !---------------------------------------------------------------------

     deallocate(x, y, dx, dy)
     deallocate(grid%areaCt, grid%anglet, grid%angle, grid%angchk)
     deallocate(grid%Ct%lat, grid%Ct%lon)
     deallocate(grid%Cv%lat, grid%Cv%lon)
     deallocate(grid%Cu%lat, grid%Cu%lon)
     deallocate(grid%Bu%lat, grid%Bu%lon)
  endif ! if (maintask)
#endif
end program gen_fixgrid
