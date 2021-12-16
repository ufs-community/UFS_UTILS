module topoedits

  use gengrid_kinds, only: real_kind,int_kind
  use grdvars,       only: ni,nj,nv,mastertask
  use grdvars,       only: wet4
  use charstrings,   only: logmsg,history
  use netcdf

  implicit none
  private

  public add_topoedits

  contains

  subroutine add_topoedits(fsrc,fdst)

   character(len=*), intent(in) :: fsrc, fdst

   integer :: rc,id,i,j,ii,jj,ncid,dimid,idimid,dim1(1)
   integer :: cnt1=0, cnt2=0, icnt
   integer(int_kind), allocatable, dimension(:) :: ieds1, jeds1, ieds2, jeds2
   real(real_kind), allocatable, dimension(:)   :: zeds1, zeds2

!---------------------------------------------------------------------
! read existing topo edits
!---------------------------------------------------------------------

  rc = nf90_open(trim(fsrc), nf90_nowrite, ncid)
  print '(a)','using topo edits file '//trim(fsrc)//' to edit land mask '

  rc = nf90_inq_dimid(ncid, 'nEdits', dimid)
  rc = nf90_inquire_dimension(ncid, dimid, len=cnt1)
  rc = nf90_close(ncid)

  ! return the existing values
  allocate(ieds1(cnt1)); ieds1 = 0
  allocate(jeds1(cnt1)); jeds2 = 0
  allocate(zeds1(cnt1)); zeds2 = 0.0

  rc = nf90_open(fsrc, nf90_nowrite, ncid)
  rc = nf90_inq_varid(ncid, 'iEdit', id)
  rc = nf90_get_var(ncid, id, ieds1)
  rc = nf90_inq_varid(ncid, 'jEdit', id)
  rc = nf90_get_var(ncid, id, jeds1)
  rc = nf90_inq_varid(ncid, 'zEdit', id)
  rc = nf90_get_var(ncid, id, zeds1)
  rc = nf90_close(ncid)

!---------------------------------------------------------------------
! determine the number of points to be added to topo-edits file
! check only j=1 now
!---------------------------------------------------------------------

  icnt = 0
     j = 1
  do i = 1,ni
   if(wet4(i,j) .eq. 1.0)icnt = icnt+1
  end do
  cnt2 = cnt2+icnt
  print '(a,i4,a,i4)', 'found ',icnt,' open water points at j=1 , cnt2 = ',cnt2

  cnt2 = cnt1 + cnt2
  ! allocate space for existing+new values and copy in original values
  allocate(ieds2(cnt2)); ieds2 = 0
  allocate(jeds2(cnt2)); jeds2 = 0
  allocate(zeds2(cnt2)); zeds2 = 0.0
 
  ieds2(1:cnt1) = ieds1(1:cnt1)
  jeds2(1:cnt1) = jeds1(1:cnt1)
  zeds2(1:cnt1) = zeds1(1:cnt1)

!---------------------------------------------------------------------
! fill in new values and write new topoedits file
!---------------------------------------------------------------------

  icnt = cnt1
   j = 1
  do i = 1,ni
   if(wet4(i,j) .eq. 1.0)then
     icnt = icnt+1
     ii = i-1; jj = j-1
     ieds2(icnt) = ii
     jeds2(icnt) = jj
     zeds2(icnt) = 0.0
   end if
  end do

  !do i = 1,cnt2
  ! print '(3i5,f12.4)',i,ieds2(i),jeds2(i),zeds2(i)
  !end do

  rc = nf90_create(fdst, nf90_write, ncid)
  print '(a)', 'writing new topo edits to '//trim(fdst)

  rc = nf90_def_dim(ncid, 'nEdits', cnt2, idimid)

  rc = nf90_def_var(ncid,    'ni', nf90_int,   id)
  rc = nf90_def_var(ncid,    'nj', nf90_int,   id)

  dim1(:) = (/idimid/)
  rc = nf90_def_var(ncid, 'iEdit', nf90_int,   dim1, id)
  rc = nf90_def_var(ncid, 'jEdit', nf90_int,   dim1, id)
  rc = nf90_def_var(ncid, 'zEdit', nf90_float, dim1, id)
  rc = nf90_put_att(ncid, nf90_global, 'history', trim(history))
  rc = nf90_enddef(ncid)

  rc = nf90_inq_varid(ncid,     'ni',     id)
  rc = nf90_put_var(ncid,         id,     ni)
  rc = nf90_inq_varid(ncid,     'nj',     id)
  rc = nf90_put_var(ncid,         id,     nj)

  rc = nf90_inq_varid(ncid,  'iEdit',     id)
  rc = nf90_put_var(ncid,         id,  ieds2)
  rc = nf90_inq_varid(ncid,  'jEdit',     id)
  rc = nf90_put_var(ncid,         id,  jeds2)
  rc = nf90_inq_varid(ncid,  'zEdit',     id)
  rc = nf90_put_var(ncid,         id,  zeds2)

  rc = nf90_close(ncid)

!---------------------------------------------------------------------
! adjust land mask by same edits used at run time
!---------------------------------------------------------------------

  do i = 1,cnt2
   ii = ieds2(i); jj = jeds2(i)
   if(wet4(ii+1,jj+1) .eq. 0.0 .and. zeds2(i) .gt. 0.0) then
     wet4(ii+1,jj+1) = 1.0
     print '(a,2i4,a)', 'switch point ',ii+1,jj+1,' from land->ocean at runtime'
   end if
   if(wet4(ii+1,jj+1) .eq. 1.0 .and. zeds2(i) .eq. 0.0) then
     wet4(ii+1,jj+1) = 0.0
     print '(a,2i4,a)', 'switch point ',ii+1,jj+1,' from ocean->land at runtime'
   end if
  end do

  end subroutine add_topoedits
end module topoedits
