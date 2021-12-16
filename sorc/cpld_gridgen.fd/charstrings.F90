module charstrings

  use gengrid_kinds, only : CL,CM,CS

  implicit none

  character(len=CL) :: dirsrc, dirout, fv3dir
  character(len=CS) :: res, atmres
  character(len=CL) :: logmsg

  character(len=CL) :: maskfile  = 'ocean_mask.nc'
  character(len=CS) :: maskname  = 'mask'
  character(len=CL) :: editsfile

  character(len=CL) :: topofile
  character(len=CS) :: toponame  = 'depth'

  character(len=CL) :: history
  character(len=CS) :: cdate

  character(len= 2), dimension(4) :: staggerlocs = (/'Ct','Cu','Cv','Bu'/)

end module charstrings
