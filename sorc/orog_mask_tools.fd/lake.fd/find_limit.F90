!> @file
!! @brief Geo-reference utilities for a cubed-sphere grid.
!! @author Ning Wang

!> Given two points on a cubed-sphere grid, compute the
!! maximum and minimum latitudinal extent of the 
!! resulting great circle.
!!
!! @param[in] p1_in Latitude and longitude of point 1.
!! @param[in] p2_in Latitude and longitude of point 2.
!! @param[out] latmin Minimum latitudinal extent.
!! @param[out] latmax Maximum latitudinal extent.
!! @author Ning Wang
!#define DIAG
SUBROUTINE find_limit (p1_in, p2_in, latmin, latmax)
    REAL*8, INTENT(IN) :: p1_in(2), p2_in(2)
    REAL*8, INTENT(OUT) :: latmin, latmax

    REAL*8 :: p1(2),p2(2), pm(2)
    REAL*8 :: r2d = 180.0/acos(-1.0)

    p1 = p1_in/r2d; p2 = p2_in/r2d
    latmin = min(p1(1), p2(1))
    latmax = max(p1(1), p2(1))

    CALL middle(p1, p2, pm)
#ifdef DIAG
    PRINT*, 'before loop', p1(1)*r2d,p2(1)*r2d,pm(1)*r2d
#endif
    
    DO WHILE (abs(p1(1)-p2(1)) > 0.00001 .AND. &
              abs(p1(2)-p2(2)) > 0.00001 )
      IF (abs(p1(1)-pm(1)) < abs(p2(1)-pm(1))) THEN
        p2 = pm
      ELSE
        p1 = pm
      ENDIF
      CALL middle(p1, p2, pm)
#ifdef DIAG
      PRINT*, 'in loop', p1(1)*r2d,p2(1)*r2d, pm(1)*r2d
#endif
    ENDDO

    latmin = min(latmin, pm(1))
    latmax = max(latmax, pm(1))

    latmin = latmin *r2d
    latmax = latmax *r2d

END SUBROUTINE find_limit

!> Compute the latitude and longitude of the middle
!! point between two given points.
!!
!! There are two formulae available to compute it.
!!  
!! One derived from a more general m-sect formula:
!!  <pre>
!!  xyz = sin((1-f)*theta) / sin(theta) * xyz1 +
!!        sin(f*theta) /sin(theta) * xyz2 ;
!!  where theta is the angle of xyz1, and xyz2.
!!  </pre>
!!
!!  <pre>
!!  xyz = 0.5 / sqrt[(1+dot(xyz1,xyz2))/2] * (xyz1+xyz2)
!!  </pre>
!!
!! and the other one is the normalized middle point of
!! the two end points:
!!
!!  <pre>
!!  xyz = 0.5 * (xyz1+xyz2), xyz = xyz / sqrt(dot(xyz,xyz))
!!  </pre>
!!
!! @param[in] p1 Latitude/longitude of first end point.
!! @param[in] p2 Latitude/longitude of second end point
!! @param[out] p Latitude/longitude of the mid-point.
!! @author Ning Wang @date March, 2006
SUBROUTINE middle(p1,p2,p)
     IMPLICIT NONE

     ! Two given points in lat/lon:
     REAL*8, INTENT(IN) :: p1(2),p2(2)
     REAL*8, INTENT(OUT) :: p(2)
     REAL*8 :: pi

     REAL*8 :: xyz1(3),xyz2(3),xyz(3)

     pi = acos(-1.0) 

     ! Convert them into Cardesian coor:
     xyz1(1) = cos(p1(1)) * cos(p1(2))
     xyz1(2) = cos(p1(1)) * sin(p1(2))
     xyz1(3) = sin(p1(1))

     xyz2(1) = cos(p2(1)) * cos(p2(2))
     xyz2(2) = cos(p2(1)) * sin(p2(2))
     xyz2(3) = sin(p2(1))

     ! middle point:

!    coeff = 0.5 / sqrt((1.0 + dot_product(xyz1,xyz2)) / 2) 
!    xyz = coeff * (xyz1 + xyz2)

     xyz = 0.5 * (xyz1 + xyz2)

     xyz = xyz / sqrt(dot_product(xyz,xyz))

     ! Convert the middle point to lat/lon coor:
     p(1) = atan2(xyz(3), sqrt(xyz(1) * xyz(1) + xyz(2) * xyz(2))) 
     p(2) = atan2(xyz(2), xyz(1)) 

     IF (p(2) < -pi / 2.0) THEN 
       p(2) = p(2) + 2 * pi
     END IF

END SUBROUTINE middle

