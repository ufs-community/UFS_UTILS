!> @file
!! @brief Program and associated routines for creating inland mask.
!! @author Ning Wang @date July 1, 2020

!> This program creates the inland mask and writes it to the orography data files.
!!
!! @author Ning Wang
!! @return 0 for success, error otherwise.
PROGRAM inland_mask
  USE cs_nb
  IMPLICIT NONE

  INTEGER :: tile, i, j
  TYPE(nb_gp_idx) :: nbs

  REAL, ALLOCATABLE :: inland(:,:,:)
  REAL, ALLOCATABLE :: land_frac(:,:,:)
  INTEGER :: i_ctr, j_ctr, tile_beg, tile_end
  INTEGER :: cs_res, x_res, y_res
  CHARACTER(len=32) :: arg
  INTEGER :: stat
  INTEGER :: max_rd 
  REAL :: cutoff
  CHARACTER(len=1) :: reg

  LOGICAL, ALLOCATABLE :: done(:,:,:)

  CALL getarg(0, arg) ! get the program name
  IF (iargc() /= 3 .AND. iargc() /= 4) THEN
    PRINT*, 'Usage: '
    PRINT*, trim(arg), ' [cres(48,96, ...)][nonland cutoff][max recursive depth][regional(r),global(g)]'  
    STOP
  ENDIF

  CALL getarg(1, arg)
  READ(arg,*,iostat=stat) cs_res
  CALL getarg(2, arg)
  READ(arg,*,iostat=stat) cutoff
  CALL getarg(3, arg)
  READ(arg,*,iostat=stat) max_rd
  CALL getarg(4, reg)

  IF (reg /= 'r') THEN
    tile_beg = 1; tile_end = 6
! read in orography data (land_frac)
    CALL read_orog(cs_res)

! init inter-panel neighbor index
    CALL idx_init(cs_res)

! create a inland mask 
    CALL mark_global_inland(cs_res)

! write back to the orography data files 
    CALL write_inland(cs_res)
  ELSE
! read in orography data (land_frac) for the SAR domain
    CALL read_orog_reg(cs_res)

! init inter-panel neighbor index
    CALL idx_init_reg(x_res, y_res)

! create a inland mask 
    CALL mark_inland_reg(cs_res)

! write back to the orography data files 
    CALL write_inland_reg(cs_res)
  ENDIF
  
  CALL free_mem()

CONTAINS

!> Create inland mask for global grid.
!!
!! @param[in] cs_res cubed sphere resolution (48, 96...)
!!
!! @author Ning Wang  
SUBROUTINE mark_global_inland(cs_res)
  INTEGER, INTENT(IN) :: cs_res
 
  ALLOCATE(done(cs_res,cs_res,6))
  ALLOCATE(inland(cs_res,cs_res,6))
  done = .false.
  inland = 1.0
  i_ctr = cs_res/2; j_ctr = cs_res/2

  CALL mark_global_inland_rec_d(i_ctr, j_ctr, 2, 0)
  
  DEALLOCATE(done)

END SUBROUTINE mark_global_inland

!> Create inland mask for regional grid.
!!
!! @param[in] cs_res cubed sphere resolution (48, 96...)
!!
!! @author Ning Wang
SUBROUTINE mark_inland_reg(cs_res)
  INTEGER, INTENT(IN) :: cs_res
  INTEGER :: i_seed, j_seed

  ALLOCATE(done(x_res,y_res,1))
  ALLOCATE(inland(x_res,y_res,1))

  done = .false.
  inland = 1.0

! set up three seeds, for the current domain
  i_seed = 1; j_seed = y_res/2 
  CALL mark_regional_inland_rec_d(i_seed, j_seed, 1, 0)

  i_seed = x_res; j_seed = y_res/2 
  CALL mark_regional_inland_rec_d(i_seed, j_seed, 1, 0)

  i_seed = x_res/3; j_seed = 1 
  CALL mark_regional_inland_rec_d(i_seed, j_seed, 1, 0)

  DEALLOCATE(done)

END SUBROUTINE mark_inland_reg

!> Recursively walk through neighbors marking inland points for global grid.
!!
!! @param[in] i cell index
!! @param[in] j cell index
!! @param[in] t tile number
!! @param[in] rd current recursive depth (starts at 0)
!!
!! @author Ning Wang
RECURSIVE SUBROUTINE mark_global_inland_rec_d(i, j, t, rd)
  INTEGER, INTENT(IN) :: i, j, t, rd

  TYPE(nb_gp_idx) :: nbs
  INTEGER :: k, nrd
  
  IF (t == 0) RETURN

  IF (land_frac(i,j,t) <= 0.15) THEN
    nrd = 1
  ELSE
    nrd = rd + 1
  ENDIF

  IF (nrd > max_rd) RETURN

  IF (done(i,j,t)) RETURN
  IF (land_frac(i,j,t) < cutoff) THEN
    done(i,j,t) = .true.
    inland(i,j,t) = 0.0
    CALL neighbors(t, i, j, nbs)
    ! recursively go through k neighbors
    DO k = 1, 4
      CALL mark_global_inland_rec_d(nbs%ijt(1,k),nbs%ijt(2,k),nbs%ijt(3,k),nrd)
    ENDDO
  ENDIF
  
END SUBROUTINE mark_global_inland_rec_d

!> Recursively walk through neighbors marking inland points for regional grid.
!!
!! @param[in] i cell index
!! @param[in] j cell index
!! @param[in] t tile face
!! @param[in] rd current recursive depth (starts at 0)
!!
!! @author Ning Wang
RECURSIVE SUBROUTINE mark_regional_inland_rec_d(i, j, t, rd)
  INTEGER, INTENT(IN) :: i, j, t, rd

  TYPE(nb_gp_idx) :: nbs
  INTEGER :: k, nrd
  
  IF (t == 0) RETURN

  IF (land_frac(i,j,t) <= 0.15) THEN
    nrd = 1
  ELSE
    nrd = rd + 1
  ENDIF

  IF (nrd > max_rd) RETURN

  IF (done(i,j,t)) RETURN
  IF (land_frac(i,j,t) < cutoff) THEN
    done(i,j,t) = .true.
    inland(i,j,t) = 0.0
    CALL neighbors_reg(i, j, nbs)
    ! recursively go through k neighbors
    DO k = 1, 4
      CALL mark_regional_inland_rec_d(nbs%ijt(1,k),nbs%ijt(2,k),nbs%ijt(3,k),nrd)
    ENDDO
  ENDIF
  
END SUBROUTINE mark_regional_inland_rec_d

!> Read in orography (land fraction) data.
!!
!! @param[in] cs_res cubed sphere resolution (48, 96...)
!!
!! @author Ning Wang
SUBROUTINE read_orog(cs_res)
  USE netcdf 
  INTEGER, INTENT(IN) :: cs_res

  INTEGER :: tile_sz, tile_num
  INTEGER :: stat, ncid, varid
  INTEGER :: land_frac_id, slmsk_id, geolon_id, geolat_id
  CHARACTER(len=256) :: filename,string
  CHARACTER(len=1) :: ich
  CHARACTER(len=4) res_ch
  CHARACTER(len=8) dimname

  INTEGER :: i, j
  REAL, ALLOCATABLE :: var_tmp(:,:)

  ALLOCATE(var_tmp(cs_res,cs_res))
  ALLOCATE(land_frac(cs_res,cs_res,6))

  WRITE(res_ch,'(I4)') cs_res
  DO tile_num = tile_beg, tile_end
    WRITE(ich, '(I1)') tile_num
    filename = "oro.C" // trim(adjustl(res_ch)) // ".tile" // ich // ".nc" 
    print *,'Read, update, and write ',trim(filename)
    stat = nf90_open(filename, NF90_NOWRITE, ncid)
    CALL nc_opchk(stat, "nf90_open oro_data.nc")

! original orodata netcdf file uses (y, x) order, so we made change to match it. 
    stat = nf90_inq_varid(ncid, "land_frac", land_frac_id)
    CALL nc_opchk(stat, "nf90_inq_varid: land_frac")
    stat = nf90_get_var(ncid, land_frac_id, var_tmp, &
           start = (/ 1, 1 /), count = (/ cs_res, cs_res /) )
    CALL nc_opchk(stat, "nf90_get_var: land_frac")
    land_frac(:,:,tile_num) = var_tmp(:,:)  
    stat = nf90_close(ncid)
    CALL nc_opchk(stat, "nf90_close oro_data.nc")
  ENDDO

  DEALLOCATE(var_tmp)

END SUBROUTINE read_orog

!> Read in orography (land fraction) data for regional grid.
!!
!! @param[in] cs_res cubed sphere resolution (48, 96...)
!!
!! @author Ning Wang
SUBROUTINE read_orog_reg(cs_res)
  USE netcdf 
  INTEGER, INTENT(IN) :: cs_res

  INTEGER :: tile_num
  INTEGER :: stat, ncid, x_dimid, y_dimid, varid
  INTEGER :: land_frac_id, slmsk_id, geolon_id, geolat_id
  CHARACTER(len=256) :: filename,string
  CHARACTER(len=1) :: ich
  CHARACTER(len=4) res_ch
  CHARACTER(len=8) dimname

  INTEGER :: i, j
  REAL, ALLOCATABLE :: var_tmp(:,:)

  WRITE(res_ch,'(I4)') cs_res
  tile_num = 7
  WRITE(ich, '(I1)') tile_num
  filename = "oro.C" // trim(adjustl(res_ch)) // ".tile" // ich // ".nc" 
  print *,'Read, update, and write ',trim(filename)
  stat = nf90_open(filename, NF90_NOWRITE, ncid)
  CALL nc_opchk(stat, "nf90_open oro_data.nc")
  stat = nf90_inq_dimid(ncid, "lon", x_dimid)
  CALL nc_opchk(stat, "nf90_inq_dim: x")
  stat = nf90_inq_dimid(ncid, "lat", y_dimid)
  CALL nc_opchk(stat, "nf90_inq_dim: y")
  stat = nf90_inquire_dimension(ncid,x_dimid,dimname,len=x_res)
  CALL nc_opchk(stat,'nf90_inquire_dimension: lon')
  stat = nf90_inquire_dimension(ncid,y_dimid,dimname,len=y_res)
  CALL nc_opchk(stat,'nf90_inquire_dimension: lon')

! original orodata netcdf file uses (y, x) order, so we made change to match it. 
  ALLOCATE(var_tmp(x_res,y_res))
  ALLOCATE(land_frac(x_res,y_res,1))
  stat = nf90_inq_varid(ncid, "land_frac", land_frac_id)
  CALL nc_opchk(stat, "nf90_inq_varid: land_frac")
  stat = nf90_get_var(ncid, land_frac_id, var_tmp, &
         start = (/ 1, 1 /), count = (/ x_res, y_res /) )
  CALL nc_opchk(stat, "nf90_get_var: land_frac")
  land_frac(:,:,1) = var_tmp(:,:)  
  stat = nf90_close(ncid)
  CALL nc_opchk(stat, "nf90_close oro_data.nc")

  DEALLOCATE(var_tmp)

END SUBROUTINE read_orog_reg

!> Write inland back to the orography data files for global grid
!!
!! @param[in] cs_res cubed sphere resolution (48, 96...)
!!
!! @author Ning Wang
SUBROUTINE write_inland(cs_res)
  USE netcdf 
  INTEGER, INTENT(IN) :: cs_res

  CHARACTER(len=256) :: filename
  CHARACTER(len=1) :: ich
  CHARACTER(len=4) res_ch

  INTEGER :: tile_num
  INTEGER :: stat, ncid, x_dimid, y_dimid, inland_id, dimids(2)
  REAL, ALLOCATABLE :: var_tmp(:,:)

  ALLOCATE(var_tmp(cs_res,cs_res))

  WRITE(res_ch,'(I4)') cs_res
  DO tile_num = tile_beg, tile_end
    WRITE(ich, '(I1)') tile_num
    filename = "oro.C" // trim(adjustl(res_ch)) // ".tile" // ich // ".nc" 
    print *,'write inland to ',trim(filename)
    stat = nf90_open(filename, NF90_WRITE, ncid)
    CALL nc_opchk(stat, "nf90_open oro_data.nc")
    stat = nf90_inq_dimid(ncid, "lon", x_dimid)
    CALL nc_opchk(stat, "nf90_inq_dim: x")
    stat = nf90_inq_dimid(ncid, "lat", y_dimid)
    CALL nc_opchk(stat, "nf90_inq_dim: y")

! original orodata netcdf file uses (y, x) order, so we made change to match it. 
    dimids = (/ x_dimid, y_dimid /)

! define a new variables 
    stat = nf90_inq_varid(ncid,"inland",inland_id) !safeguard if "inland" exists
    IF (stat /= nf90_noerr) THEN 
      stat = nf90_redef(ncid)
      CALL nc_opchk(stat, "nf90_redef")
      stat = nf90_def_var(ncid,"inland",NF90_FLOAT,dimids,inland_id)
      CALL nc_opchk(stat, "nf90_def_var: inland")
      stat = nf90_put_att(ncid, inland_id,'coordinates','geolon geolat')
      CALL nc_opchk(stat, "nf90_put_att: inland:coordinates") 
      stat = nf90_put_att(ncid, inland_id,'description', &
          'inland = 1 indicates grid cells away from coast')
      CALL nc_opchk(stat, "nf90_put_att: inland:description") 
      stat = nf90_enddef(ncid) 
      CALL nc_opchk(stat, "nf90_enddef")
    ENDIF

    var_tmp(:,:) = inland(:,:,tile_num)
    stat = nf90_put_var(ncid, inland_id, var_tmp, &
           start = (/ 1, 1 /), count = (/ cs_res, cs_res /) )
    CALL nc_opchk(stat, "nf90_put_var: inland") 

    stat = nf90_close(ncid)
    CALL nc_opchk(stat, "nf90_close oro_data.nc")
  ENDDO
  DEALLOCATE(var_tmp)

END SUBROUTINE write_inland

!> Write inland back to the orography data files for regional grid.
!!
!! @param[in] cs_res cubed sphere resolution (48, 96...)
!!
!! @author Ning Wang
SUBROUTINE write_inland_reg(cs_res)
  USE netcdf 
  INTEGER, INTENT(IN) :: cs_res

  CHARACTER(len=256) :: filename
  CHARACTER(len=1) :: ich
  CHARACTER(len=4) res_ch

  INTEGER :: tile_num
  INTEGER :: stat, ncid, x_dimid, y_dimid, inland_id, dimids(2)
  REAL, ALLOCATABLE :: var_tmp(:,:)
  CHARACTER(len=8) dimname

  ALLOCATE(var_tmp(x_res,y_res))

  WRITE(res_ch,'(I4)') cs_res
  tile_num = 7 
  WRITE(ich, '(I1)') tile_num
  filename = "oro.C" // trim(adjustl(res_ch)) // ".tile" // ich // ".nc" 
  PRINT*,'write inland to ',trim(filename)
  stat = nf90_open(filename, NF90_WRITE, ncid)
  CALL nc_opchk(stat, "nf90_open oro_data.nc")
  stat = nf90_inq_dimid(ncid, "lon", x_dimid)
  CALL nc_opchk(stat, "nf90_inq_dim: x")
  stat = nf90_inq_dimid(ncid, "lat", y_dimid)
  CALL nc_opchk(stat, "nf90_inq_dim: y")
  stat = nf90_inquire_dimension(ncid,x_dimid,dimname,len=x_res)
  CALL nc_opchk(stat,'nf90_inquire_dimension: lon')
  stat = nf90_inquire_dimension(ncid,y_dimid,dimname,len=y_res)
  CALL nc_opchk(stat,'nf90_inquire_dimension: lon')

! original orodata netcdf file uses (y, x) order, so we made change to match it. 
  dimids = (/ x_dimid, y_dimid /)

! define a new variables 
  stat = nf90_inq_varid(ncid,"inland",inland_id) !safeguard if "inland" exists
  IF (stat /= nf90_noerr) THEN 
    stat = nf90_redef(ncid)
    CALL nc_opchk(stat, "nf90_redef")
    stat = nf90_def_var(ncid,"inland",NF90_FLOAT,dimids,inland_id)
    CALL nc_opchk(stat, "nf90_def_var: inland")
    stat = nf90_put_att(ncid, inland_id,'coordinates','geolon geolat')
    CALL nc_opchk(stat, "nf90_put_att: inland:coordinates") 
    stat = nf90_put_att(ncid, inland_id,'description', &
        'inland = 1 indicates grid cells away from coast')
    CALL nc_opchk(stat, "nf90_put_att: inland:description") 
    stat = nf90_enddef(ncid) 
    CALL nc_opchk(stat, "nf90_enddef")
  ENDIF

  var_tmp(:,:) = inland(:,:,1)
  stat = nf90_put_var(ncid, inland_id, var_tmp, &
         start = (/ 1, 1 /), count = (/ x_res, y_res /) )
  CALL nc_opchk(stat, "nf90_put_var: inland") 

  stat = nf90_close(ncid)
  CALL nc_opchk(stat, "nf90_close oro_data.nc")

  DEALLOCATE(var_tmp)

END SUBROUTINE write_inland_reg

!> Deallocate module arrays
SUBROUTINE free_mem()

  DEALLOCATE(inland, land_frac)

END SUBROUTINE free_mem

!> Check NetCDF return code and print error message.
!!
!! @param[in] stat status code returned from NetCDF
!! @param[in] opname description of NetCDF operation called
!!
!! @author Ning Wang
SUBROUTINE nc_opchk(stat,opname)
   USE netcdf
   IMPLICIT NONE
   INTEGER stat
   CHARACTER(len=*) opname
   CHARACTER(64) msg

   IF (stat .NE.0)  THEN
     msg = trim(opname) // ' Error, status code and message:'
     PRINT*,trim(msg), stat, nf90_strerror(stat)
     STOP
   END IF

END SUBROUTINE nc_opchk

END PROGRAM inland_mask

