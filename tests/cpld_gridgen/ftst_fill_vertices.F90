!> Unit test for fill_vertices routine
!!
!! This test checks the filling of vertex values for each stagger
!! location. This test relies on the known Arakawa C-grid (NE) staggering.
!! Vertices are ordered counter-clockwise from upper right
!!
!!          Bu----Cv(i,j)---Bu (i,j)
!!           |              |
!!          Cu    Ct(i,j)   Cu (i,j)
!!           |              |
!! (i-1,j-1) Bu----Cv-------Bu
!!
!! @author Denise.Worthen@noaa.gov
program ftst_fill_vertices

  use assertion_mod, only: assert_equal
  use gengrid_kinds, only: dbl_kind, int_kind, CL, CS
  use grdvars      , only: stagger_type, nv, iVertNE, jVertNE
  use gengrid_utils, only: allocate_staggers
  use vertices     , only: fill_vertices

  implicit none

  integer, parameter     :: nx = 5, ny = 4
  integer, parameter     :: maxtests = 50, nresults = maxtests

  logical                :: ispassing(nresults)
  character(len=CL)      :: testmsg(nresults) = ' '

  type(stagger_type) :: Ct, Cu, Cv, Bu

  real(dbl_kind) :: testverts(nv)

  integer :: nt, ntests
  integer :: i, j, n
  character(len=CS) :: stagger
  character(len=CL) :: msg, msg_out
  logical :: status

  ! Initialize global test data; coordinate encoded
  call allocate_staggers(nx,ny,Ct)
  call allocate_staggers(nx,ny,Cu)
  call allocate_staggers(nx,ny,Cv)
  call allocate_staggers(nx,ny,Bu)

  do j = 1,ny
     do i = 1,nx
        Ct%lon(i,j) = 10.0_dbl_kind * i + j
        Ct%lat(i,j) = 10.0_dbl_kind * i + j
     end do
  end do
  Bu%lon = Ct%lon + 0.5; Bu%lat = Ct%lat + 0.5
  Cu%lon = Ct%lon + 0.5; Cu%lat = Ct%lat
  Cv%lon = Ct%lon      ; Cv%lat = Ct%lat + 0.5

  do i = 1,nx
     Ct%xlon(i) = 100_dbl_kind * i
     Ct%xlat(i) = 100_dbl_kind * i
  end do
  Bu%xlon = Ct%xlon + 50.0; Bu%xlat = Ct%xlat + 50.0
  Cu%xlon = Ct%xlon + 50.0; Cu%xlat = Ct%xlat
  Cv%xlon = Ct%xlon       ; Cv%xlat = Ct%xlat + 50.0

  Ct%iVert = iVertNE
  Ct%jVert = jVertNE
  Cu%iVert = Ct%iVert + 1; Cu%jVert = Ct%jvert + 0
  Cv%iVert = Ct%iVert + 0; Cv%jVert = Ct%jVert + 1
  Bu%iVert = Ct%iVert + 1; Bu%jVert = Ct%jVert + 1

  nt = 0; testverts = 0.0_dbl_kind

  !------ Ct vertices from Bu grid ------!
  call fill_vertices(Ct%iVert, Ct%jVert, Bu%lat, Bu%lon, Bu%xlat, Bu%xlon, Ct%latvert, Ct%lonvert, 0)

  ! single loc
  i = 3; j = 2
  testverts(:) = (/Bu%lat(i,j),Bu%lat(i-1,j),Bu%lat(i-1,j-1),Bu%lat(i,j-1)/)
  nt = nt + 1; msg = 'Ct lat vertices from Bu grid at 3,2'
  call assert_equal(Ct%latvert(3,2,:),testverts,0.0_dbl_kind,msg,status,msg_out)
  ispassing(nt) = status
  testmsg(nt) = trim(msg_out)

  ! outer corners
  i = 1; j = 1
  testverts(:) = (/Bu%lat(i,j),Bu%lat(nx,j),Bu%xlat(nx),Bu%xlat(i)/)
  nt = nt + 1; msg = 'Ct lat vertices from Bu grid at 1,1'
  call assert_equal(Ct%latvert(1,1,:),testverts,0.0_dbl_kind,msg,status,msg_out)
  ispassing(nt) = status
  testmsg(nt) = trim(msg_out)

  i = nx; j = 1
  testverts(:) = (/Bu%lon(i,j),Bu%lon(i-1,j),Bu%xlon(i-1),Bu%xlon(i)/)
  nt = nt + 1; msg = 'Ct lon vertices from Bu grid at nx,1'
  call assert_equal(Ct%lonvert(nx,1,:),testverts,0.0_dbl_kind,msg,status,msg_out)
  ispassing(nt) = status
  testmsg(nt) = trim(msg_out)

  i = nx; j = ny
  testverts(:) = (/Bu%lat(i,j),Bu%lat(i-1,j),Bu%lat(i-1,j-1),Bu%lat(i,j-1)/)
  nt = nt + 1; msg = 'Ct lat vertices from Bu grid at nx,ny'
  call assert_equal(Ct%latvert(nx,ny,:),testverts,0.0_dbl_kind,msg,status,msg_out)
  ispassing(nt) = status
  testmsg(nt) = trim(msg_out)

  i = 1; j = ny
  testverts(:) = (/Bu%lon(i,j),Bu%lon(nx,j),Bu%lon(nx,j-1),Bu%lon(i,j-1)/)
  nt = nt + 1; msg = 'Ct lon vertices from Bu grid at 1,ny'
  call assert_equal(Ct%lonvert(1,ny,:),testverts,0.0_dbl_kind,msg,status,msg_out)
  ispassing(nt) = status
  testmsg(nt) = trim(msg_out)

  !------ Cu vertices from Cv grid ------!
  call fill_vertices(Cu%iVert, Cu%jVert, Cv%lat, Cv%lon, Cv%xlat, Cv%xlon, Cu%latvert, Cu%lonvert, 0)

  ! single loc
  i = 2; j = 3
  testverts(:) = (/Cv%lat(i+1,j),Cv%lat(i,j),Cv%lat(i,j-1),Cv%lat(i+1,j-1)/)
  nt = nt + 1; msg = 'Cu lat vertices from Cv grid at 2,3'
  call assert_equal(Cu%latvert(2,3,:),testverts,0.0_dbl_kind,msg,status,msg_out)
  ispassing(nt) = status
  testmsg(nt) = trim(msg_out)

  ! outer corners
  i = 1; j = 1
  testverts(:) = (/Cv%lat(i+1,j),Cv%lat(i,j),Cv%xlat(i),Cv%xlat(i+1)/)
  nt = nt + 1; msg = 'Cu lat vertices from Cv grid at 1,1'
  call assert_equal(Cu%latvert(1,1,:),testverts,0.0_dbl_kind,msg,status,msg_out)
  ispassing(nt) = status
  testmsg(nt) = trim(msg_out)

  i = nx; j = 1
  testverts(:) = (/Cv%lon(1,j),Cv%lon(nx,j),Cv%xlon(nx),Cv%xlon(1)/)
  nt = nt + 1; msg = 'Cu lon vertices from Cv grid at nx,1'
  call assert_equal(Cu%lonvert(nx,1,:),testverts,0.0_dbl_kind,msg,status,msg_out)
  ispassing(nt) = status
  testmsg(nt) = trim(msg_out)

  i = nx; j = ny
  testverts(:) = (/Cv%lat(1,j),Cv%lat(i,j),Cv%lat(i,j-1),Cv%lat(1,j-1)/)
  nt = nt + 1; msg = 'Cu lat vertices from Cv grid at nx,ny'
  call assert_equal(Cu%latvert(nx,ny,:),testverts,0.0_dbl_kind,msg,status,msg_out)
  ispassing(nt) = status
  testmsg(nt) = trim(msg_out)

  i = 1; j = ny
  testverts(:) = (/Cv%lon(i+1,j),Cv%lon(i,j),Cv%lon(i,j-1),Cv%lon(i+1,j-1)/)
  nt = nt + 1; msg = 'Cu lon vertices from Cv grid at 1,ny'
  call assert_equal(Cu%lonvert(1,ny,:),testverts,0.0_dbl_kind,msg,status,msg_out)
  ispassing(nt) = status
  testmsg(nt) = trim(msg_out)

  !------ Cv vertices from Cu grid ------!
  call fill_vertices(Cv%iVert, Cv%jVert, Cu%lat, Cu%lon, Cu%xlat, Cu%xlon, Cv%latvert, Cv%lonvert)

  ! single loc
  i = 2; j = 2
  testverts(:) = (/Cu%lat(i,j+1),Cu%lat(i-1,j+1),Cu%lat(i-1,j),Cu%lat(i,j)/)
  nt = nt + 1; msg = 'Cv lat vertices from Cu grid at 3,3'
  call assert_equal(Cv%latvert(2,2,:),testverts,0.0_dbl_kind,msg,status,msg_out)
  ispassing(nt) = status
  testmsg(nt) = trim(msg_out)

  ! outer corners
  i = 1; j = 1
  testverts(:) = (/Cu%lat(i,j+1),Cu%lat(nx,j+1),Cu%lat(nx,j),Cu%lat(i,j)/)
  nt = nt + 1; msg = 'Cv lat vertices from Cu grid at 1,1'
  call assert_equal(Cv%latvert(1,1,:),testverts,0.0_dbl_kind,msg,status,msg_out)
  ispassing(nt) = status
  testmsg(nt) = trim(msg_out)

  i = nx; j = 1
  testverts(:) = (/Cu%lon(i,j+1),Cu%lon(i-1,j+1),Cu%lon(i-1,j),Cu%lon(i,j)/)
  nt = nt + 1; msg = 'Cv lon vertices from Cu grid at nx,1'
  call assert_equal(Cv%lonvert(nx,1,:),testverts,0.0_dbl_kind,msg,status,msg_out)
  ispassing(nt) = status
  testmsg(nt) = trim(msg_out)

  i = nx; j = ny
  testverts(:) = (/Cu%xlat(nx),Cu%xlat(i-1),Cu%lat(i-1,j),Cu%lat(i,j)/)
  nt = nt + 1; msg = 'Cv lat vertices from Cu grid at nx,ny'
  call assert_equal(Cv%latvert(nx,ny,:),testverts,0.0_dbl_kind,msg,status,msg_out)
  ispassing(nt) = status
  testmsg(nt) = trim(msg_out)

  i = 1; j = ny
  testverts(:) = (/Cu%xlon(i),Cu%xlon(nx),Cu%lon(nx,j),Cu%lon(i,j)/)
  nt = nt + 1; msg = 'Cv lon vertices from Cu grid at 1,ny'
  call assert_equal(Cv%lonvert(1,ny,:),testverts,0.0_dbl_kind,msg,status,msg_out)
  ispassing(nt) = status
  testmsg(nt) = trim(msg_out)

  ! force fail
  i = 2; j = ny
  testverts(:) = (/Cu%xlat(nx),Cu%xlat(i-1)+1.0e-8,Cu%lat(i-1,j),Cu%lat(i,j)/)
  nt = nt + 1; msg = 'Cv lat vertices from Cu grid at 2,ny'
  call assert_equal(Cv%latvert(2,ny,:),testverts,0.0_dbl_kind,msg,status,msg_out)
  if (.not. status) ispassing(nt) = .true.    ! failure was caught
  testmsg(nt) = 'Expected '//trim(msg_out)

  !------ Bu vertices from Ct grid ------!
  ! Bu vertics from Ct grid
  call fill_vertices(Bu%iVert, Bu%jVert, Ct%lat, Ct%lon, Ct%xlat, Ct%xlon, Bu%latvert, Bu%lonvert)

  ! single loc
  i = 3; j = 3
  testverts(:) = (/Ct%lat(i+1,j+1),Ct%lat(i,j+1),Ct%lat(i,j),Ct%lat(i+1,j)/)
  nt = nt + 1; msg = 'Bu lat vertices from Ct grid at 3,3'
  call assert_equal(Bu%latvert(3,3,:),testverts,0.0_dbl_kind,msg,status,msg_out)
  ispassing(nt) = status
  testmsg(nt) = trim(msg_out)

  ! outer corners
  i = 1; j = 1
  testverts(:) = (/Ct%lat(i+1,j+1),Ct%lat(i,j+1),Ct%lat(i,j),Ct%lat(i+1,j)/)
  nt = nt + 1; msg = 'Bu lat vertices from Ct grid at 1,1'
  call assert_equal(Bu%latvert(1,1,:),testverts,0.0_dbl_kind,msg,status,msg_out)
  ispassing(nt) = status
  testmsg(nt) = trim(msg_out)

  i = nx; j = 1
  testverts(:) = (/Ct%lon(1,j+1),Ct%lon(i,j+1),Ct%lon(i,j),Ct%lon(1,j)/)
  nt = nt + 1; msg = 'Bu lon vertices from Ct grid at nx,1'
  call assert_equal(Bu%lonvert(nx,1,:),testverts,0.0_dbl_kind,msg,status,msg_out)
  ispassing(nt) = status
  testmsg(nt) = trim(msg_out)

  i = nx; j = ny
  testverts(:) = (/Ct%xlat(1),Ct%xlat(i),Ct%lat(i,j),Ct%lat(1,j)/)
  nt = nt + 1; msg = 'Bu lat vertices from Ct grid at nx,ny'
  call assert_equal(Bu%latvert(nx,ny,:),testverts,0.0_dbl_kind,msg,status,msg_out)
  ispassing(nt) = status
  testmsg(nt) = trim(msg_out)

  i = 1; j = ny
  testverts(:) = (/Ct%xlon(i+1),Ct%xlon(i),Ct%lon(i,j),Ct%lon(i+1,j)/)
  nt = nt + 1; msg = 'Bu lon vertices from Ct grid at 1,ny'
  call assert_equal(Bu%lonvert(1,ny,:),testverts,0.0_dbl_kind,msg,status,msg_out)
  ispassing(nt) = status
  testmsg(nt) = trim(msg_out)

  ! force fail
  i = 3; j = ny
  testverts(:) = (/Ct%xlat(1)+1.0e-8,Ct%xlat(i),Ct%lat(i,j),Ct%lat(1,j)/)
  nt = nt + 1; msg = 'Bu lat vertices from Ct grid at nx,ny'
  call assert_equal(Bu%latvert(nx,ny,:),testverts,0.0_dbl_kind,msg,status,msg_out)
  if (.not. status) ispassing(nt) = .true.    ! failure was caught
  testmsg(nt) = 'Expected '//trim(msg_out)

  ntests = nt
  if (all(ispassing(1:ntests))) then
     print '(a)', 'All unit tests passed '
  else
     print '(a)', 'FAIL: At least one unit test failed '
     stop 1
  end if

end program ftst_fill_vertices
