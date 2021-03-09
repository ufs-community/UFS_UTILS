!> @file
!! @brief Return number of threads.
!! @author Mark Iredell NCEP/EMC

!> Return the number of omp threads.
!!
!! @return num_parthds The number of threads.
!! @author Mark Iredell NCEP/EMC
integer function num_parthds()
 use omp_lib
!$OMP PARALLEL
 num_parthds=omp_get_num_threads()
!$OMP END PARALLEL
 return
 end function num_parthds
