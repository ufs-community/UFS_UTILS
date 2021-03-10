!> @file
!! @brief Driver for vcoord_gen.
!! @author Mark Iredell @date 2008-08-01

!> Driver for vcoord_gen.
!! @return 0 for success, error code otherwise.
!! @author Mark Iredell @date 2008-08-01
program driver
  implicit none
  integer levs,lupp,k
  real pbot,psig,ppre,pupp,ptop,dpbot,dpsig,dppre,dpupp,dptop,pmin
  real,allocatable:: ak(:),bk(:)
  write(0,*) 'Enter levs,lupp,pbot,psig,ppre,pupp,ptop,dpbot,dpsig,dppre,dpupp,dptop'
  read *,levs,lupp,pbot,psig,ppre,pupp,ptop,dpbot,dpsig,dppre,dpupp,dptop
  allocate(ak(0:levs),bk(0:levs))
  call vcoord_gen(levs,lupp,pbot,psig,ppre,pupp,ptop,dpbot,dpsig,dppre,dpupp,dptop,pmin,ak,bk)
  write(0,*) 'pmin=',pmin
  print '(2i6)',2,levs
  print '(f12.3,f12.8)',(ak(k),bk(k),k=0,levs)
end program
