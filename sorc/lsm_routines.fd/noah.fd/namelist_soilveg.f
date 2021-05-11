C> @file
C> @brief Module holding variables for Noah LSM
C> @author ? 

      module namelist_soilveg
      implicit none
      save

!> Routine copied from noah LSM, small edits for Doxygen compilation 
!! Clara Draper, May, 2021.

! Draper, hard-coded to git testing
      !INTEGER, PARAMETER ::  MAX_SLOPETYP = 30 
      !INTEGER, PARAMETER :: MAX_SOILTYP = 30 
      !INTEGER, PARAMETER ::  MAX_VEGTYP = 30 

      REAL SLOPE_DATA(30)
      REAL RSMTBL(30)
      REAL RGLTBL(30)
      REAL HSTBL(30)
      REAL SNUPX(30)
      REAL BB(30)
      REAL DRYSMC(30)
      REAL F11(30)
      REAL MAXSMC(30)
      REAL REFSMC(30)
      REAL SATPSI(30)
      REAL SATDK(30)
      REAL SATDW(30)
      REAL WLTSMC(30)
      REAL QTZ(30)
      LOGICAL LPARAM
      REAL ZBOT_DATA
      REAL SALP_DATA
      REAL CFACTR_DATA
      REAL CMCMAX_DATA
      REAL SBETA_DATA
      REAL RSMAX_DATA
      REAL TOPT_DATA
      REAL REFDK_DATA
      REAL FRZK_DATA
      INTEGER BARE
      INTEGER DEFINED_VEG
      INTEGER DEFINED_SOIL
      INTEGER DEFINED_SLOPE
      REAL FXEXP_DATA
      INTEGER NROOT_DATA(30)
      REAL REFKDT_DATA
      REAL Z0_DATA(30)
      REAL CZIL_DATA
      REAL LAI_DATA(30)
      REAL CSOIL_DATA

      end module namelist_soilveg
