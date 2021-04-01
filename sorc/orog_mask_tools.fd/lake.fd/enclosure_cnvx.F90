!> @file
!! @brief Test whether a point is inside a spherical convex polygon.
!! @author N Wang
#ifdef INCLUDE_TEST_DRIVER
    PROGRAM testenc
    IMPLICIT NONE
    REAL*8 :: v(2,4)
    REAL*8 :: p(2)

    REAL*8 :: d2r
    LOGICAL:: enclosure_cnvx, inside
    INTEGER :: co_gc
   
    d2r = acos(-1.0)/180.0D0

    v(1,1) = 10.0D0*d2r; v(2,1) = 20.0D0*d2r
    v(1,2) = 15.0D0*d2r; v(2,2) = 30.0D0*d2r
    v(1,3) = 17.7D0*d2r; v(2,3) = 25.0D0*d2r
    v(1,4) = 20.0D0*d2r; v(2,4) = 20.0D0*d2r

!    p(1) = 15.0D0*d2r; p(2) = 30.00000001D0*d2r
!    p(1) = 20.00000000D0*d2r; p(2) = 20.0D0*d2r
!    p(1) = 9.999999999D0*d2r; p(2) = 20.0D0*d2r
!    p(1) = 10.00000000*d2r; p(2) = 20.000000001D0*d2r
    p(1) = 17.7D0*d2r; p(2) = 25.000000001D0*d2r
 
    inside = enclosure_cnvx(v,4,p,co_gc)
    IF (inside) THEN
       PRINT*, 'inside ', co_gc
    ELSE 
       PRINT*, 'outside ', co_gc
    ENDIF
    
    END PROGRAM
#endif

!> Test whether a given point 'p' is inside a convex spherical polygon defined with a series of 'n' vertices.
!! Both the test point and the polygon are specified in latitude and longitude radians.
!! It is assumed that a side of a spherical polygon is a great-circle arc that connects 2 adjacent vertices of the polygon.
!!    
!! @param[in] v set of vertices
!! @param[in] n number of vertices in v
!! @param[in] p point to test
!! @param[out] co_gc 'i' if @p p is on or within epsilon to side 'i' of the given spherical polygon,
!! or 0 if point is not on any side of the given polygon.
!! @return enclosure_cnvx whether the point is inside the polygon
!!
!! @author N Wang   
LOGICAL FUNCTION enclosure_cnvx(v, n, p, co_gc)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: v(2,n), p(2)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(OUT) :: co_gc

    REAL*8 v_xy(2, n)
    REAL*8 cp_z(n), cos_d2c, eps
    
    INTEGER :: i, ip1

    
    eps = 0.000000000000001D0
    co_gc = 0
    DO i = 1, n
      cos_d2c = sin(p(1))*sin(v(1,i)) + cos(p(1))*cos(v(1,i))*cos(v(2,i)-p(2)) 
      v_xy(1,i) = (cos(v(1,i))*sin(v(2,i)-p(2)))/cos_d2c
      v_xy(2,i) = (cos(p(1))*sin(v(1,i))-sin(p(1))*cos(v(1,i))*cos(v(2,i)-p(2)))/cos_d2c

    ENDDO

    DO i = 1, n
      ip1 = mod(i,n)+1
      cp_z(i) = v_xy(1,i)*v_xy(2,ip1)-v_xy(2,i)*v_xy(1,ip1)
      IF (abs(cp_z(i)) < eps) co_gc = i
    ENDDO

    DO i = 1, n-1
      ip1 = mod(i,n)+1
      IF (cp_z(i)*cp_z(ip1) .LT. -eps) THEN
        enclosure_cnvx = .false.
        RETURN
      ENDIF
    ENDDO

    enclosure_cnvx = .true.
    RETURN

END FUNCTION enclosure_cnvx
