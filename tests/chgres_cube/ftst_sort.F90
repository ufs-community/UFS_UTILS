 program sort

 use input_data

 implicit none

 real*8 :: a(5)
 integer :: first, last

 data a /300.0, 600.0, 200.0, 50.0, 1000.0/

 print*,'hello world'

 first = 1
 last = 5
 call quicksort(a,first,last)

 print*,'after sort 1 ',a,first,last

 end program sort
