!> Unit test for reshape_staggers routine
!!
!! This test checks the reshaping of staggered grid points from (i,j) or (i,j,nv)
!! to vector equivalents for both a global domain and an extracted subdomain
!!
!! @author Denise.Worthen@noaa.gov
program ftst_reshape_staggers

  use assertion_mod, only: assert_equal
  use gengrid_kinds, only: dbl_kind, int_kind, CL
  use gengrid_utils, only: reshape_staggers
  use grdvars      , only: nv

  implicit none

  integer, parameter :: nx = 5, ny = 4
  integer, parameter :: ngrids = 2
  integer, parameter :: maxtests = 50, nresults = ngrids*maxtests

  logical           :: ispassing(nresults)
  character(len=8)  :: gridname(ngrids) = (/'Global  ','Regional'/)
  character(len=CL) :: testmsg(nresults) = ' '

  integer           :: iind(2), jind(2)
  ! test data
  real(dbl_kind)    :: lon(nx,ny), lat(nx,ny)
  integer(int_kind) :: mask(nx,ny)
  real(dbl_kind)    :: lonvert(nx,ny,nv), latvert(nx,ny,nv)
  ! result data
  real(dbl_kind), allocatable    :: cnlons(:), cnlats(:)
  integer(int_kind), allocatable :: cnmask(:)
  real(dbl_kind), allocatable    :: crlons(:, :), crlats(:, :)

  integer :: ng, nt, ntests
  integer :: i, j, n
  integer :: ib, ie, jb, je, idim, jdim
  integer :: idx, jdx, idx1

  character(len=CL) :: msg, msg_out
  logical :: status

  ! Initialize global test data; coordinate encoded
  do j = 1, ny
     do i = 1, nx
        lon(i,j) =  10.0_dbl_kind * i + j
        lat(i,j) = -10.0_dbl_kind * i - j
        do n = 1, nv
           lonvert(i,j,n) = i*100.0 + j*10.0 + n
           latvert(i,j,n) = i*100.0 + j*10.0 + n
        end do
     end do
  end do
  do j = 1,ny
     mask(1:2,j) = 1_int_kind
     mask(:,3) = 0_int_kind
  end do

  nt = 0
  ispassing = .false.
  do ng = 1,ngrids
     if (ng .eq. 1) then
        iind = (/1, nx/)
        jind = (/1, ny/)
     else
        iind = (/3,5/)
        jind = (/1,2/)
     end if

     ib = iind(1) ; ie = iind(2)
     jb = jind(1) ; je = jind(2)
     idim = (ie - ib) + 1
     jdim = (je - jb) + 1

     ! point location in global domain also present in subdomain
     idx = 4; jdx = 2
     if (ng .eq. 1) then
        idx1 = (jdx - 1) * nx + idx
     else
        idx1 = 5
     end if

     allocate(cnlons(idim*jdim), source=0.0_dbl_kind)
     allocate(cnlats(idim*jdim), source=0.0_dbl_kind)
     allocate(cnmask(idim*jdim), source = 1_int_kind)
     allocate(crlons(nv,idim*jdim), source = 0.0_dbl_kind)
     allocate(crlats(nv,idim*jdim), source = 0.0_dbl_kind)

     call reshape_staggers((/ib,ie/),(/jb,je/), lon, lat, mask, lonvert, latvert, cnlons, cnlats, cnmask, crlons, crlats)

     ! compare start index
     nt = nt+1; msg = trim(gridname(ng))//' compare lon index (1,1) to index (1)'
     call assert_equal(cnlons(1),lon(ib,jb),0.0_dbl_kind,msg,status,msg_out)
     ispassing(nt) = status
     testmsg(nt) = trim(msg_out)

     nt = nt+1; msg = trim(gridname(ng))//' compare lat index (1,1) to index (1)'
     call assert_equal(cnlats(1),lat(ib,jb),0.0_dbl_kind,msg,status,msg_out)
     ispassing(nt) = status
     testmsg(nt) = trim(msg_out)

     nt = nt+1; msg = trim(gridname(ng))//' compare mask index (1,1) to index (1)'
     call assert_equal(cnmask(1),mask(ib,jb),msg,status,msg_out)
     ispassing(nt) = status
     testmsg(nt) = trim(msg_out)

     !  compare end index
     nt = nt+1; msg = trim(gridname(ng))//' compare lon index (ie,je) to index (idim*jdim)'
     call assert_equal(cnlons(idim*jdim),lon(ie,je),0.0_dbl_kind,msg,status,msg_out)
     ispassing(nt) = status
     testmsg(nt) = trim(msg_out)

     nt = nt+1; msg = trim(gridname(ng))//' compare lat index (ie,je) to index (idim*jdim)'
     call assert_equal(cnlats(idim*jdim),lat(ie,je),0.0_dbl_kind,msg,status,msg_out)
     ispassing(nt) = status
     testmsg(nt) = trim(msg_out)

     nt = nt+1; msg = trim(gridname(ng))//' compare mask index (ie,je) to index (idim*jdim)'
     call assert_equal(cnmask(idim*jdim),mask(ie,je),msg,status,msg_out)
     ispassing(nt) = status
     testmsg(nt) = trim(msg_out)

     !  compare reshaped corner arrays
     nt = nt+1; msg = trim(gridname(ng))//' compare crlons(:,1) to lonvert(iind(1),jind(1),:)'
     call assert_equal(crlons(:,1),lonvert(ib,jb,:),0.0_dbl_kind,msg,status,msg_out)
     ispassing(nt) = status
     testmsg(nt) = trim(msg_out)

     nt = nt+1; msg = trim(gridname(ng))//' compare crlats(:,1) to latvert(iind(1),jind(1),:)'
     call assert_equal(crlats(:,1),latvert(ib,jb,:),0.0_dbl_kind,msg,status,msg_out)
     ispassing(nt) = status
     testmsg(nt) = trim(msg_out)

     nt = nt+1; msg = trim(gridname(ng))//' compare crlons(:,idim*jdim) to lonvert(ie,je,:)'
     call assert_equal(crlons(:,idim*jdim),lonvert(ie,je,:),0.0_dbl_kind,msg,status,msg_out)
     ispassing(nt) = status
     testmsg(nt) = trim(msg_out)

     nt = nt+1; msg = trim(gridname(ng))//' compare crlats(:,idim*jdim) to latvert(ie,je,:)'
     call assert_equal(crlats(:,idim*jdim),latvert(ie,je,:),0.0_dbl_kind,msg,status,msg_out)
     ispassing(nt) = status
     testmsg(nt) = trim(msg_out)

     !  compare single point loc
     nt = nt+1; msg = trim(gridname(ng))//' compare lon index (idx,jdx) to index (idx1)'
     call assert_equal(cnlons(idx1),lon(idx,jdx),0.0_dbl_kind,msg,status,msg_out)
     ispassing(nt) = status
     testmsg(nt) = trim(msg_out)

     nt = nt+1; msg = trim(gridname(ng))//' compare lat index (idx,jdx) to index (idx1)'
     call assert_equal(cnlats(idx1),lat(idx,jdx),0.0_dbl_kind,msg,status,msg_out)
     ispassing(nt) = status
     testmsg(nt) = trim(msg_out)

     nt = nt+1; msg = trim(gridname(ng))//' compare mask index (idx,jdx) to index (idx1)'
     call assert_equal(cnmask(idx1),mask(idx,jdx),msg,status,msg_out)
     ispassing(nt) = status
     testmsg(nt) = trim(msg_out)

     !  compare single point loc, corners
     nt = nt+1; msg = trim(gridname(ng))//' compare crlons(:,idx1) to lonvert(idx,jdx,:)'
     call assert_equal(crlons(:,idx1),lonvert(idx,jdx,:),0.0_dbl_kind,msg,status,msg_out)
     ispassing(nt) = status
     testmsg(nt) = trim(msg_out)

     nt = nt+1; msg = trim(gridname(ng))//' compare crlats(:,idx1) to latvert(idx,jdx,:)'
     call assert_equal(crlats(:,idx1),latvert(idx,jdx,:),0.0_dbl_kind,msg,status,msg_out)
     ispassing(nt) = status
     testmsg(nt) = trim(msg_out)

     ! Mimic bad return values from reshape_arrays
     nt = nt+1; msg = trim(gridname(ng))//' compare lon index (1,1) to index (2)'
     cnlons(1) = cnlons(2)
     call assert_equal(cnlons(1),lon(ib,ib),0.0_dbl_kind,msg,status,msg_out)
     if (.not. status) ispassing(nt) = .true.    ! failure was caught
     testmsg(nt) = 'Expected '//trim(msg_out)

     nt = nt+1; msg = trim(gridname(ng))//' compare crlats(3,idim*jdim) to latvert(idim,jdim,2)'
     crlats(3,idim*jdim) = crlats(2,idim*jdim)
     call assert_equal(crlats(:,idim*jdim),latvert(idim,jdim,:),0.0_dbl_kind,msg,status,msg_out)
     if (.not. status) ispassing(nt) = .true.    ! failure was caught
     testmsg(nt) = 'Expected '//trim(msg_out)

     nt = nt+1; msg = trim(gridname(ng))//' compare lat index (idx,jdx) to index (idx1)+eps'
     cnlats(idx1) = cnlats(idx1)+1.0e-10
     call assert_equal(cnlats(idx1),lat(idx,jdx),0.0_dbl_kind,msg,status,msg_out)
     if (.not. status) ispassing(nt) = .true.    ! failure was caught
     testmsg(nt) = 'Expected '//trim(msg_out)

     deallocate(cnlons, cnlats, cnmask, crlats, crlons)
  end do

  ntests = nt
  if (all(ispassing(1:ntests))) then
     print '(a)', 'All unit tests passed '
  else
     print '(a)', 'FAIL: At least one unit test failed '
     stop 1
  end if

end program ftst_reshape_staggers
