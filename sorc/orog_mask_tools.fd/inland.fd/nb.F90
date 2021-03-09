!> @file
!! This module contains the subroutines that find the adjacent neighbors of a
!! given cell, in a cubed sphere grid. Each neighbor is in the form of (i,j,tile).
!! 
!! @author Ning Wang @date July 1, 2020
!!
MODULE cs_nb
  IMPLICIT NONE

  !> Neighboring tile descriptor.
  TYPE nb_tile_idx
    INTEGER       :: nb_tile_num !< Neighbor tile number (1..6)
    CHARACTER (1) :: nb_tile_bndry !< Neighbor tile boundary (l)eft, (r)ight, (t)op, (b)ottom
  END TYPE nb_tile_idx

  !> Neighboring cell descriptor.
  TYPE nb_gp_idx
    INTEGER  :: gp_type !< Cell boundary type from bndry function
    INTEGER  :: ijt(3,8) !< Neighboring cell indices
  END TYPE nb_gp_idx

  TYPE(nb_tile_idx):: nb_tile(4,6) !< Descriptor for each edge of each tile on the cubed sphere.

  INTEGER :: cres !< Cubed sphere resolution
  integer :: xres !< x resolution of regional grid
  integer :: yres !< y resolution of regional grid
 
CONTAINS
  
  !> Initialize inter-panel neighbor index for global grid.
  !!
  !! @param[in] cres_in cubed sphere resolution (48, 96...)
  !!
  !! @verbatim  
  !!   _______1_______
  !!  |               |       1-upper, 2-bottom, 3-left, 4-right 
  !!  |               |
  !!  |               |
  !! 3|               |4
  !!  |               |
  !!  |               |
  !!  |_______________|
  !!          2
  !!      Figure 1. Boundary numbers
  !! @endverbatim
  !!
  !! @author Ning Wang
  SUBROUTINE idx_init(cres_in) 
    INTEGER :: cres_in

    nb_tile(1,1)%nb_tile_num = 3; nb_tile(1,1)%nb_tile_bndry ='l'
    nb_tile(2,1)%nb_tile_num = 6; nb_tile(2,1)%nb_tile_bndry ='t'
    nb_tile(3,1)%nb_tile_num = 5; nb_tile(3,1)%nb_tile_bndry ='t'
    nb_tile(4,1)%nb_tile_num = 2; nb_tile(4,1)%nb_tile_bndry ='l'

    nb_tile(1,2)%nb_tile_num = 3; nb_tile(1,2)%nb_tile_bndry ='b'
    nb_tile(2,2)%nb_tile_num = 6; nb_tile(2,2)%nb_tile_bndry ='r'
    nb_tile(3,2)%nb_tile_num = 1; nb_tile(3,2)%nb_tile_bndry ='r'
    nb_tile(4,2)%nb_tile_num = 4; nb_tile(4,2)%nb_tile_bndry ='b'

    nb_tile(1,3)%nb_tile_num = 5; nb_tile(1,3)%nb_tile_bndry ='l'
    nb_tile(2,3)%nb_tile_num = 2; nb_tile(2,3)%nb_tile_bndry ='t'
    nb_tile(3,3)%nb_tile_num = 1; nb_tile(3,3)%nb_tile_bndry ='t'
    nb_tile(4,3)%nb_tile_num = 4; nb_tile(4,3)%nb_tile_bndry ='l'

    nb_tile(1,4)%nb_tile_num = 5; nb_tile(1,4)%nb_tile_bndry ='b'
    nb_tile(2,4)%nb_tile_num = 2; nb_tile(2,4)%nb_tile_bndry ='r'
    nb_tile(3,4)%nb_tile_num = 3; nb_tile(3,4)%nb_tile_bndry ='r'
    nb_tile(4,4)%nb_tile_num = 6; nb_tile(4,4)%nb_tile_bndry ='b'

    nb_tile(1,5)%nb_tile_num = 1; nb_tile(1,5)%nb_tile_bndry ='l'
    nb_tile(2,5)%nb_tile_num = 4; nb_tile(2,5)%nb_tile_bndry ='t'
    nb_tile(3,5)%nb_tile_num = 3; nb_tile(3,5)%nb_tile_bndry ='t'
    nb_tile(4,5)%nb_tile_num = 6; nb_tile(4,5)%nb_tile_bndry ='l'

    nb_tile(1,6)%nb_tile_num = 1; nb_tile(1,6)%nb_tile_bndry ='b'
    nb_tile(2,6)%nb_tile_num = 4; nb_tile(2,6)%nb_tile_bndry ='r'
    nb_tile(3,6)%nb_tile_num = 5; nb_tile(3,6)%nb_tile_bndry ='r'
    nb_tile(4,6)%nb_tile_num = 2; nb_tile(4,6)%nb_tile_bndry ='b'
    
    cres = cres_in

  END SUBROUTINE idx_init

  !> Initialize resolution module variables for regional grid.
  !!
  !! @param[in] xres_in x resolution
  !! @param[in] yres_in y resolution
  !!
  !! @author Ning Wang
  SUBROUTINE idx_init_reg(xres_in, yres_in)
    INTEGER, INTENT(IN) :: xres_in, yres_in

    xres = xres_in
    yres = yres_in

  END SUBROUTINE idx_init_reg

  !> Get boundary type from indices for global grid.
  !!
  !! @param[in] i cell index
  !! @param[in] j cell index
  !! @return bndry cell boundary type
  !!
  !! @author Ning Wang
  INTEGER FUNCTION bndry(i, j)
    INTEGER :: i,j

    bndry = 0                ! no boundary

    IF (j == cres) THEN      ! upper boundary
      bndry = 1
      IF (i == 1) THEN
        bndry = 13
      ELSE IF (i == cres) THEN
        bndry = 14
      ENDIF  
    ELSE IF (j == 1) THEN    ! bottom boundary
      bndry = 2
      IF (i == 1) THEN
        bndry = 23
      ELSE IF (i == cres) THEN
        bndry = 24
      ENDIF  
    ELSE IF (i == 1) THEN    ! left boundary
      bndry = 3
    ELSE IF (i == cres) THEN ! right boundary
      bndry = 4
    ENDIF

  END FUNCTION bndry
  
  !> Get boundary type from indices for regional grid.
  !!
  !! @param[in] i cell index
  !! @param[in] j cell index
  !! @return cell boundary type
  !!
  !! @author Ning Wang
  INTEGER FUNCTION bndry_reg(i, j)
    INTEGER :: i,j

    bndry_reg = 0                ! no boundary

    IF (j == yres) THEN      ! upper boundary
      bndry_reg = 1
      IF (i == 1) THEN
        bndry_reg = 13
      ELSE IF (i == xres) THEN
        bndry_reg = 14
      ENDIF  
    ELSE IF (j == 1) THEN    ! bottom boundary
      bndry_reg = 2
      IF (i == 1) THEN
        bndry_reg = 23
      ELSE IF (i == xres) THEN
        bndry_reg = 24
      ENDIF  
    ELSE IF (i == 1) THEN    ! left boundary
      bndry_reg = 3
    ELSE IF (i == xres) THEN ! right boundary
      bndry_reg = 4
    ENDIF

  END FUNCTION bndry_reg

  !> Get neighbors of cell 'c' at (tile, i, j) for global grid.
  !!
  !!
  !! @verbatim
  !!     ______________
  !!    |    |    |    |              ________
  !!    | 5  | 1  | 6  |             /\ 1 \ 6 \
  !!    |____|____|____|            /  \___\___\
  !!    |    |    |    |           /\2 / c / 3 /
  !!    | 2  | c  | 3  |          /  \/___/___/
  !!    |____|____|____|          \7 / 4 / 8 /
  !!    |    |    |    |           \/___/___/
  !!    | 7  | 4  | 8  |       
  !!    |____|____|____|    
  !!  
  !! Figure 2.  Eight neighbors of cell 'c' and special cases at upper left
  !! cornner of the tile
  !! @endverbatim
  !!
  !! @param[in] tile tile face
  !! @param[in] i cell index
  !! @param[in] j cell index
  !! @param[out] nb neighbors
  !!
  !! @author Ning Wang
  SUBROUTINE neighbors(tile, i, j, nb)
    INTEGER :: tile, i, j
    TYPE(nb_gp_idx) :: nb  

    INTEGER :: bd, nb_t_num
    
    nb%gp_type = bndry(i,j)
    IF (nb%gp_type == 0) THEN  ! interior (non-boundary) cell
      ! top, bottom, left, and right
      nb%ijt(1,1) = i; nb%ijt(2,1) = j+1; nb%ijt(3,1) = tile 
      nb%ijt(1,2) = i-1; nb%ijt(2,2) = j; nb%ijt(3,2) = tile 
      nb%ijt(1,3) = i+1; nb%ijt(2,3) = j; nb%ijt(3,3) = tile 
      nb%ijt(1,4) = i; nb%ijt(2,4) = j-1; nb%ijt(3,4) = tile 
      ! top left, top right, bottom left, and bottom right
      nb%ijt(1,5) = i-1; nb%ijt(2,5) = j+1; nb%ijt(3,5) = tile 
      nb%ijt(1,6) = i+1; nb%ijt(2,6) = j+1; nb%ijt(3,6) = tile 
      nb%ijt(1,7) = i-1; nb%ijt(2,7) = j-1; nb%ijt(3,7) = tile 
      nb%ijt(1,8) = i+1; nb%ijt(2,8) = j-1; nb%ijt(3,8) = tile 
    ELSEIF (nb%gp_type == 1) THEN  ! top boundary cell
      bd = 1
      nb_t_num = nb_tile(nb%gp_type,tile)%nb_tile_num 
      nb%ijt(3,1)=nb_t_num; nb%ijt(3,5)=nb_t_num; nb%ijt(3,6)=nb_t_num
      IF (nb_tile(bd,tile)%nb_tile_bndry == 'l') THEN
        nb%ijt(1,1) = 1; nb%ijt(2,1) = cres+1-i; 
        nb%ijt(1,5) = 1; nb%ijt(2,5) = cres+1-(i-1); 
        nb%ijt(1,6) = 1; nb%ijt(2,6) = cres+1-(i+1); 
      ELSEIF (nb_tile(bd,tile)%nb_tile_bndry == 'b') THEN
        nb%ijt(1,1) = i; nb%ijt(2,1) = 1 
        nb%ijt(1,5) = i-1; nb%ijt(2,5) = 1 
        nb%ijt(1,6) = i+1; nb%ijt(2,6) = 1
      ENDIF
      nb%ijt(1,2) = i-1; nb%ijt(2,2) = j; nb%ijt(3,2) = tile 
      nb%ijt(1,3) = i+1; nb%ijt(2,3) = j; nb%ijt(3,3) = tile 
      nb%ijt(1,4) = i; nb%ijt(2,4) = j-1; nb%ijt(3,4) = tile 
      nb%ijt(1,7) = i-1; nb%ijt(2,7) = j-1; nb%ijt(3,7) = tile 
      nb%ijt(1,8) = i+1; nb%ijt(2,8) = j-1; nb%ijt(3,8) = tile 
    ELSEIF (nb%gp_type == 2) THEN ! bottom boundary cell
      bd = 2
      nb_t_num = nb_tile(nb%gp_type,tile)%nb_tile_num 
      nb%ijt(3,4)=nb_t_num; nb%ijt(3,7)=nb_t_num; nb%ijt(3,8)=nb_t_num
      IF (nb_tile(bd,tile)%nb_tile_bndry == 'r') THEN
        nb%ijt(1,4) = cres; nb%ijt(2,4) = cres+1-i; 
        nb%ijt(1,7) = cres; nb%ijt(2,7) = cres+1-(i-1); 
        nb%ijt(1,8) = cres; nb%ijt(2,8) = cres+1-(i+1); 
      ELSEIF (nb_tile(bd,tile)%nb_tile_bndry == 't') THEN
        nb%ijt(1,4) = i; nb%ijt(2,4) = cres 
        nb%ijt(1,7) = i-1; nb%ijt(2,7) = cres 
        nb%ijt(1,8) = i+1; nb%ijt(2,8) = cres
      ENDIF
      nb%ijt(1,1) = i; nb%ijt(2,1) = j+1; nb%ijt(3,1) = tile 
      nb%ijt(1,2) = i-1; nb%ijt(2,2) = j; nb%ijt(3,2) = tile 
      nb%ijt(1,3) = i+1; nb%ijt(2,3) = j; nb%ijt(3,3) = tile 
      nb%ijt(1,5) = i-1; nb%ijt(2,5) = j+1; nb%ijt(3,5) = tile 
      nb%ijt(1,6) = i+1; nb%ijt(2,6) = j+1; nb%ijt(3,6) = tile 
    ELSEIF (nb%gp_type == 3) THEN ! left boundary cell
      bd = 3
      nb_t_num = nb_tile(nb%gp_type,tile)%nb_tile_num 
      nb%ijt(3,2)=nb_t_num; nb%ijt(3,5)=nb_t_num; nb%ijt(3,7)=nb_t_num
      IF (nb_tile(bd,tile)%nb_tile_bndry == 'r') THEN
        nb%ijt(1,2) = cres; nb%ijt(2,2) = j; 
        nb%ijt(1,5) = cres; nb%ijt(2,5) = j+1; 
        nb%ijt(1,7) = cres; nb%ijt(2,7) = j-1; 
      ELSEIF (nb_tile(bd,tile)%nb_tile_bndry == 't') THEN
        nb%ijt(1,2) = cres+1-j; nb%ijt(2,2) = cres 
        nb%ijt(1,5) = cres+1-(j+1); nb%ijt(2,5) = cres 
        nb%ijt(1,7) = cres+1-(j-1); nb%ijt(2,7) = cres
      ENDIF
      nb%ijt(1,1) = i; nb%ijt(2,1) = j+1; nb%ijt(3,1) = tile 
      nb%ijt(1,3) = i+1; nb%ijt(2,3) = j; nb%ijt(3,3) = tile 
      nb%ijt(1,4) = i; nb%ijt(2,4) = j-1; nb%ijt(3,4) = tile 
      nb%ijt(1,6) = i+1; nb%ijt(2,6) = j+1; nb%ijt(3,6) = tile 
      nb%ijt(1,8) = i+1; nb%ijt(2,8) = j-1; nb%ijt(3,8) = tile 
    ELSEIF (nb%gp_type == 4) THEN ! right boundary cell
      bd = 4
      nb_t_num = nb_tile(nb%gp_type,tile)%nb_tile_num 
      nb%ijt(3,3)=nb_t_num; nb%ijt(3,6)=nb_t_num; nb%ijt(3,8)=nb_t_num
      IF (nb_tile(bd,tile)%nb_tile_bndry == 'l') THEN
        nb%ijt(1,3) = 1; nb%ijt(2,3) = j; 
        nb%ijt(1,6) = 1; nb%ijt(2,6) = j+1; 
        nb%ijt(1,8) = 1; nb%ijt(2,8) = j-1; 
      ELSEIF (nb_tile(bd,tile)%nb_tile_bndry == 'b') THEN
        nb%ijt(1,3) = cres+1-j; nb%ijt(2,3) = 1 
        nb%ijt(1,6) = cres+1-(j+1); nb%ijt(2,6) = 1 
        nb%ijt(1,8) = cres+1-(j-1); nb%ijt(2,8) = 1
      ENDIF
      nb%ijt(1,1) = i; nb%ijt(2,1) = j+1; nb%ijt(3,1) = tile 
      nb%ijt(1,2) = i-1; nb%ijt(2,2) = j; nb%ijt(3,2) = tile 
      nb%ijt(1,4) = i; nb%ijt(2,4) = j-1; nb%ijt(3,4) = tile 
      nb%ijt(1,5) = i-1; nb%ijt(2,5) = j+1; nb%ijt(3,5) = tile 
      nb%ijt(1,7) = i-1; nb%ijt(2,7) = j-1; nb%ijt(3,7) = tile 
    ELSEIF (nb%gp_type == 13) THEN ! upper left coner
      bd = 1
      nb_t_num = nb_tile(bd,tile)%nb_tile_num 
      nb%ijt(3,1)=nb_t_num; nb%ijt(3,6)=nb_t_num
      IF (nb_tile(bd,tile)%nb_tile_bndry == 'l') THEN
        nb%ijt(1,1) = 1; nb%ijt(2,1) = cres+1-i 
        nb%ijt(1,6) = 1; nb%ijt(2,6) = cres+1-(i+1) 
      ELSEIF (nb_tile(bd,tile)%nb_tile_bndry == 'b') THEN
        nb%ijt(1,1) = i; nb%ijt(2,1) = 1 
        nb%ijt(1,6) = i+1; nb%ijt(2,6) = 1
      ENDIF
      bd = 3
      nb_t_num = nb_tile(bd,tile)%nb_tile_num 
      nb%ijt(3,2)=nb_t_num; nb%ijt(3,7)=nb_t_num
      IF (nb_tile(bd,tile)%nb_tile_bndry == 'r') THEN
        nb%ijt(1,2) = cres; nb%ijt(2,2) = j
        nb%ijt(1,7) = cres; nb%ijt(2,7) = j-1
      ELSEIF (nb_tile(bd,tile)%nb_tile_bndry == 't') THEN
        nb%ijt(1,2) = cres+1-j; nb%ijt(2,2) = cres 
        nb%ijt(1,7) = cres+1-(j-1); nb%ijt(2,7) = cres
      ENDIF
      nb%ijt(3,5)=0 
      nb%ijt(1,3) = i+1; nb%ijt(2,3) = j; nb%ijt(3,3) = tile 
      nb%ijt(1,4) = i; nb%ijt(2,4) = j-1; nb%ijt(3,4) = tile 
      nb%ijt(1,8) = i+1; nb%ijt(2,8) = j-1; nb%ijt(3,8) = tile 
    ELSEIF (nb%gp_type == 14) THEN ! upper right coner
      bd = 1
      nb_t_num = nb_tile(bd,tile)%nb_tile_num 
      nb%ijt(3,1)=nb_t_num; nb%ijt(3,5)=nb_t_num
      IF (nb_tile(bd,tile)%nb_tile_bndry == 'l') THEN
        nb%ijt(1,1) = 1; nb%ijt(2,1) = cres+1-i 
        nb%ijt(1,5) = 1; nb%ijt(2,5) = cres+1-(i-1)
      ELSEIF (nb_tile(bd,tile)%nb_tile_bndry == 'b') THEN
        nb%ijt(1,1) = i; nb%ijt(2,1) = 1 
        nb%ijt(1,5) = i-1; nb%ijt(2,5) = 1 
      ENDIF
      bd = 4
      nb_t_num = nb_tile(bd,tile)%nb_tile_num 
      nb%ijt(3,3)=nb_t_num; nb%ijt(3,8)=nb_t_num
      IF (nb_tile(bd,tile)%nb_tile_bndry == 'l') THEN
        nb%ijt(1,3) = 1; nb%ijt(2,3) = j
        nb%ijt(1,8) = 1; nb%ijt(2,8) = j-1
      ELSEIF (nb_tile(bd,tile)%nb_tile_bndry == 'b') THEN
        nb%ijt(1,3) = cres+1-j; nb%ijt(2,3) = 1 
        nb%ijt(1,8) = cres+1-(j-1); nb%ijt(2,8) = 1
      ENDIF
      nb%ijt(3,6)=0
      nb%ijt(1,2) = i-1; nb%ijt(2,2) = j; nb%ijt(3,2) = tile 
      nb%ijt(1,4) = i; nb%ijt(2,4) = j-1; nb%ijt(3,4) = tile 
      nb%ijt(1,7) = i-1; nb%ijt(2,7) = j-1; nb%ijt(3,7) = tile 
    ELSEIF (nb%gp_type == 23) THEN ! lower left coner
      bd = 2
      nb_t_num = nb_tile(bd,tile)%nb_tile_num 
      nb%ijt(3,4)=nb_t_num; nb%ijt(3,8)=nb_t_num
      IF (nb_tile(bd,tile)%nb_tile_bndry == 'r') THEN
        nb%ijt(1,4) = cres; nb%ijt(2,4) = cres+1-i
        nb%ijt(1,8) = cres; nb%ijt(2,8) = cres+1-(i+1)
      ELSEIF (nb_tile(bd,tile)%nb_tile_bndry == 't') THEN
        nb%ijt(1,4) = i; nb%ijt(2,4) = cres 
        nb%ijt(1,8) = i+1; nb%ijt(2,8) = cres
      ENDIF
      bd = 3
      nb_t_num = nb_tile(bd,tile)%nb_tile_num 
      nb%ijt(3,2)=nb_t_num; nb%ijt(3,5)=nb_t_num
      IF (nb_tile(bd,tile)%nb_tile_bndry == 'r') THEN
        nb%ijt(1,2) = cres; nb%ijt(2,2) = j
        nb%ijt(1,5) = cres; nb%ijt(2,5) = j+1
      ELSEIF (nb_tile(bd,tile)%nb_tile_bndry == 't') THEN
        nb%ijt(1,2) = cres+1-j; nb%ijt(2,2) = cres 
        nb%ijt(1,5) = cres+1-(j+1); nb%ijt(2,5) = cres 
      ENDIF
      nb%ijt(3,7)=0
      nb%ijt(1,1) = i; nb%ijt(2,1) = j+1; nb%ijt(3,1) = tile 
      nb%ijt(1,3) = i+1; nb%ijt(2,3) = j; nb%ijt(3,3) = tile 
      nb%ijt(1,6) = i+1; nb%ijt(2,6) = j+1; nb%ijt(3,6) = tile 
    ELSEIF (nb%gp_type == 24) THEN ! lower right coner
      bd = 2
      nb_t_num = nb_tile(bd,tile)%nb_tile_num 
      nb%ijt(3,4)=nb_t_num; nb%ijt(3,7)=nb_t_num
      IF (nb_tile(bd,tile)%nb_tile_bndry == 'r') THEN
        nb%ijt(1,4) = cres; nb%ijt(2,4) = cres+1-i
        nb%ijt(1,7) = cres; nb%ijt(2,7) = cres+1-(i-1)
      ELSEIF (nb_tile(bd,tile)%nb_tile_bndry == 't') THEN
        nb%ijt(1,4) = i; nb%ijt(2,4) = cres 
        nb%ijt(1,7) = i-1; nb%ijt(2,7) = cres 
      ENDIF
      bd = 4
      nb_t_num = nb_tile(bd,tile)%nb_tile_num 
      nb%ijt(3,3)=nb_t_num; nb%ijt(3,6)=nb_t_num
      IF (nb_tile(bd,tile)%nb_tile_bndry == 'l') THEN
        nb%ijt(1,3) = 1; nb%ijt(2,3) = j
        nb%ijt(1,6) = 1; nb%ijt(2,6) = j+1
      ELSEIF (nb_tile(bd,tile)%nb_tile_bndry == 'b') THEN
        nb%ijt(1,3) = cres+1-j; nb%ijt(2,3) = 1 
        nb%ijt(1,6) = cres+1-(j+1); nb%ijt(2,6) = 1 
      ENDIF
      nb%ijt(3,8)=0
      nb%ijt(1,1) = i; nb%ijt(2,1) = j+1; nb%ijt(3,1) = tile 
      nb%ijt(1,2) = i-1; nb%ijt(2,2) = j; nb%ijt(3,2) = tile 
      nb%ijt(1,5) = i-1; nb%ijt(2,5) = j+1; nb%ijt(3,5) = tile 

    ENDIF
  END SUBROUTINE neighbors

  !> Get neighbors of cell 'c' at (tile, i, j) for regional grid.
  !!
  !! @param[in] i cell index
  !! @param[in] j cell index
  !! @param[out] nb neighbors
  !!
  !! @author Ning Wang
  SUBROUTINE neighbors_reg(i, j, nb)
    INTEGER :: i, j
    TYPE(nb_gp_idx) :: nb

! assign the standard interior cell neighbors as default values
    ! top, bottom, left, and right
    nb%ijt(1,1) = i; nb%ijt(2,1) = j+1; nb%ijt(3,1) = 1 
    nb%ijt(1,2) = i-1; nb%ijt(2,2) = j; nb%ijt(3,2) = 1 
    nb%ijt(1,3) = i+1; nb%ijt(2,3) = j; nb%ijt(3,3) = 1 
    nb%ijt(1,4) = i; nb%ijt(2,4) = j-1; nb%ijt(3,4) = 1 
    ! top left, top right, bottom left, and bottom right
    nb%ijt(1,5) = i-1; nb%ijt(2,5) = j+1; nb%ijt(3,5) = 1 
    nb%ijt(1,6) = i+1; nb%ijt(2,6) = j+1; nb%ijt(3,6) = 1 
    nb%ijt(1,7) = i-1; nb%ijt(2,7) = j-1; nb%ijt(3,7) = 1 
    nb%ijt(1,8) = i+1; nb%ijt(2,8) = j-1; nb%ijt(3,8) = 1 

    nb%gp_type = bndry_reg(i,j)
    IF (nb%gp_type == 1) THEN  !top boundary cell
      nb%ijt(3,1) = 0; nb%ijt(3,5) = 0; nb%ijt(3,6) = 0 
    ELSEIF (nb%gp_type == 2) THEN !bottom boundary cell
      nb%ijt(3,4) = 0; nb%ijt(3,7) = 0; nb%ijt(3,8) = 0 
    ELSEIF (nb%gp_type == 3) THEN !left boundary cell
      nb%ijt(3,2) = 0; nb%ijt(3,5) = 0; nb%ijt(3,7) = 0 
    ELSEIF (nb%gp_type == 4) THEN !right boundary cell
      nb%ijt(3,3) = 0; nb%ijt(3,6) = 0; nb%ijt(3,8) = 0 
    ELSEIF (nb%gp_type == 13) THEN ! upper left coner
      nb%ijt(3,1) = 0; nb%ijt(3,5) = 0; nb%ijt(3,6) = 0; 
      nb%ijt(3,2) = 0; nb%ijt(3,7) = 0 
    ELSEIF (nb%gp_type == 14) THEN ! upper right coner
      nb%ijt(3,1) = 0; nb%ijt(3,5) = 0; nb%ijt(3,6) = 0; 
      nb%ijt(3,3) = 0; nb%ijt(3,8) = 0 
    ELSEIF (nb%gp_type == 23) THEN ! lower left coner
      nb%ijt(3,4) = 0; nb%ijt(3,7) = 0; nb%ijt(3,8) = 0 
      nb%ijt(3,2) = 0; nb%ijt(3,5) = 0;  
    ELSEIF (nb%gp_type == 24) THEN ! lower right coner
      nb%ijt(3,4) = 0; nb%ijt(3,7) = 0; nb%ijt(3,8) = 0 
      nb%ijt(3,3) = 0; nb%ijt(3,6) = 0; 
    ENDIF

  END SUBROUTINE neighbors_reg
    
END MODULE cs_nb

#ifdef TEST_CS_NB
PROGRAM test_nb
  USE cs_nb

  INTEGER :: tile, i, j
  TYPE(nb_gp_idx) :: nbs

  INTEGER, PARAMETER :: res = 96

  CALL idx_init(res)
!  tile = 1; i = 1; j = 5
!  tile = 1; i = 5; j = 1
!  tile = 1; i = 96; j = 10
  tile = 2; i = 96; j = 96
  CALL neighbors(tile, i, j, nbs)
  print*, 'tile = ', tile, ' i = ',i, ' j = ', j
  print*, 'nbs%type: ', nbs%gp_type
  print*, 'nbs%ijt: '
  print*, nbs%ijt 
END PROGRAM test_nb
#endif
