 program ftst_rh2spfh_gfs 

 use grib2_util 

 implicit none

 real, parameter    :: eps=1E-4

 real*4, dimension(3,2) :: rh_500hpa, rh_900hpa, rh_tmp                                                    
 real*4, dimension(3,2) :: spfh_gfs_900hpa_expected, spfh_gfs_500hpa_expected                                            
 real*8, dimension(3,2) :: temp_500hpa, temp_900hpa
 real*8 :: pres

data temp_900hpa /280.9, 282.6, 281.4, 289.9, 286.5, 288.7/
data rh_900hpa /17.6, 75.2, 99.9, 95.4, 86.6, 89.1/

data temp_500hpa /254.1, 259.7, 259.7, 267.2, 262.0, 269.5/
data rh_500hpa /10.7, 9.1, 0.5, 10.1, 4.2, 0.9/

data spfh_gfs_900hpa_expected /1.286259, 6.169457, 7.554718, 12.64650, 9.210890, 10.93505/                       
data spfh_gfs_500hpa_expected /0.1519615, 0.2257435, 1.2403490E-02, 0.4856212, 0.1288972, 5.2061662E-02/                               

i_input = 3
j_input = 2

! rh to spfh at 500 hPa with GFSv15/16 conversion
rh_tmp=rh_500hpa
pres=50000.0
call rh2spfh_gfs(rh_tmp,pres,temp_500hpa)
if (any((abs(rh_tmp*1000.0-spfh_gfs_500hpa_expected)) .gt. eps)) stop 500 

! rh to spfh at 900 hPa with GFSv15/16 conversion
rh_tmp=rh_900hpa
pres=90000.0
call rh2spfh_gfs(rh_tmp,pres,temp_900hpa)
if (any((abs(rh_tmp*1000.0-spfh_gfs_900hpa_expected)) .gt. eps)) stop 900

 print*, "OK"
 print*, "SUCCESS!"

 end program ftst_rh2spfh_gfs 

