C> @file
C> @brief Computes cosines of colatitude and Gaussian weights
C> for sets of latitudes.
C>
C> ### Program History Log
C> Date | Programmer | Comments
C> -----|------------|---------
C> 96-02-20 | Iredell | Initial.
C> 97-10-20 | Iredell | Adjust precision.
C> 98-06-11 | Iredell | Generalize precision using FORTRAN 90 intrinsic.
C> 1998-12-03 | Iredell | Generalize precision further.
C> 1998-12-03 | Iredell | Uses AIX ESSL BLAS calls.
C> 2009-12-27 | D. Stark | Updated to switch between ESSL calls on an AIX platform, and Numerical Recipies calls elsewise.
C> 2010-12-30 | Slovacek | Update alignment so preprocessor does not cause compilation failure.
C> 2012-09-01 | E. Mirvis & M.Iredell | Merging & debugging linux errors of _d and _8 using generic LU factorization.
C> 2012-11-05 | E. Mirvis | Generic FFTPACK and LU lapack were removed.
C>
C> @author Iredell @date 96-02-20

C> Computes cosines of colatitude and Gaussian weights
C> for one of the following specific global sets of latitudes.
C> - Gaussian latitudes (IDRT=4)
C> - Equally-spaced latitudes including poles (IDRT=0)
C> - Equally-spaced latitudes excluding poles (IDRT=256)
C>
C> The Gaussian latitudes are located at the zeroes of the
C> Legendre polynomial of the given order. These latitudes
C> are efficient for reversible transforms from spectral space.
C> (About twice as many equally-spaced latitudes are needed.)
C> The weights for the equally-spaced latitudes are based on
C> Ellsaesser (JAM,1966). (No weight is given the pole point.)
C> Note that when analyzing grid to spectral in latitude pairs,
C> if an equator point exists, its weight should be halved.
C> This version invokes the ibm essl matrix solver.
C>
C> @param[in] IDRT grid identifier
C>  - 4 for Gaussian grid
C>  - 0 for equally-spaced grid including poles
C>  - 256 for equally-spaced grid excluding poles
C> @param[in] JMAX number of latitudes
C> @param[out] SLAT sines of latitude
C> @param[out] WLAT Gaussian weights
C>
C> @author Iredell @date 96-02-20
      SUBROUTINE SPLAT(IDRT,JMAX,SLAT,WLAT)
      REAL SLAT(JMAX),WLAT(JMAX)
      INTEGER,PARAMETER:: KD=SELECTED_REAL_KIND(15,45)
      REAL(KIND=KD):: PK(JMAX/2),PKM1(JMAX/2),PKM2(JMAX/2)
      REAL(KIND=KD):: SLATD(JMAX/2),SP,SPMAX,EPS=10.*EPSILON(SP)
      PARAMETER(JZ=50)
      REAL BZ(JZ)
      DATA BZ        / 2.4048255577,  5.5200781103,
     $  8.6537279129, 11.7915344391, 14.9309177086, 18.0710639679,
     $ 21.2116366299, 24.3524715308, 27.4934791320, 30.6346064684,
     $ 33.7758202136, 36.9170983537, 40.0584257646, 43.1997917132,
     $ 46.3411883717, 49.4826098974, 52.6240518411, 55.7655107550,
     $ 58.9069839261, 62.0484691902, 65.1899648002, 68.3314693299,
     $ 71.4729816036, 74.6145006437, 77.7560256304, 80.8975558711,
     $ 84.0390907769, 87.1806298436, 90.3221726372, 93.4637187819,
     $ 96.6052679510, 99.7468198587, 102.888374254, 106.029930916,
     $ 109.171489649, 112.313050280, 115.454612653, 118.596176630,
     $ 121.737742088, 124.879308913, 128.020877005, 131.162446275,
     $ 134.304016638, 137.445588020, 140.587160352, 143.728733573,
     $ 146.870307625, 150.011882457, 153.153458019, 156.295034268 /
      REAL:: DLT,D1=1.
      REAL AWORK((JMAX+1)/2,((JMAX+1)/2)),BWORK(((JMAX+1)/2))
      INTEGER:: JHE,JHO,J0=0
      INTEGER IPVT((JMAX+1)/2)
      PARAMETER(PI=3.14159265358979,C=(1.-(2./PI)**2)*0.25)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  GAUSSIAN LATITUDES
      IF(IDRT.EQ.4) THEN
        JH=JMAX/2
        JHE=(JMAX+1)/2
        R=1./SQRT((JMAX+0.5)**2+C)
        DO J=1,MIN(JH,JZ)
          SLATD(J)=COS(BZ(J)*R)
        ENDDO
        DO J=JZ+1,JH
          SLATD(J)=COS((BZ(JZ)+(J-JZ)*PI)*R)
        ENDDO
        SPMAX=1.
        DO WHILE(SPMAX.GT.EPS)
          SPMAX=0.
          DO J=1,JH
            PKM1(J)=1.
            PK(J)=SLATD(J)
          ENDDO
          DO N=2,JMAX
            DO J=1,JH
              PKM2(J)=PKM1(J)
              PKM1(J)=PK(J)
              PK(J)=((2*N-1)*SLATD(J)*PKM1(J)-(N-1)*PKM2(J))/N
            ENDDO
          ENDDO
          DO J=1,JH
            SP=PK(J)*(1.-SLATD(J)**2)/(JMAX*(PKM1(J)-SLATD(J)*PK(J)))
            SLATD(J)=SLATD(J)-SP
            SPMAX=MAX(SPMAX,ABS(SP))
          ENDDO
        ENDDO
CDIR$ IVDEP
        DO J=1,JH
          SLAT(J)=SLATD(J)
          WLAT(J)=(2.*(1.-SLATD(J)**2))/(JMAX*PKM1(J))**2
          SLAT(JMAX+1-J)=-SLAT(J)
          WLAT(JMAX+1-J)=WLAT(J)
        ENDDO
        IF(JHE.GT.JH) THEN
          SLAT(JHE)=0.
          WLAT(JHE)=2./JMAX**2
          DO N=2,JMAX,2
            WLAT(JHE)=WLAT(JHE)*N**2/(N-1)**2
          ENDDO
        ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  EQUALLY-SPACED LATITUDES INCLUDING POLES
      ELSEIF(IDRT.EQ.0) THEN
        JH=JMAX/2
        JHE=(JMAX+1)/2
        JHO=JHE-1
        DLT=PI/(JMAX-1)
        SLAT(1)=1.
        DO J=2,JH
          SLAT(J)=COS((J-1)*DLT)
        ENDDO
        DO JS=1,JHO
          DO J=1,JHO
            AWORK(JS,J)=COS(2*(JS-1)*J*DLT)
          ENDDO
        ENDDO
        DO JS=1,JHO
          BWORK(JS)=-D1/(4*(JS-1)**2-1)
        ENDDO

        call ludcmp(awork,jho,jhe,ipvt)
        call lubksb(awork,jho,jhe,ipvt,bwork)

        WLAT(1)=0.
        DO J=1,JHO
          WLAT(J+1)=BWORK(J)
        ENDDO
CDIR$ IVDEP
        DO J=1,JH
           print *, j, jmax, JMAX+1-J
          SLAT(JMAX+1-J)=-SLAT(J)
          WLAT(JMAX+1-J)=WLAT(J)
        ENDDO
        IF(JHE.GT.JH) THEN
          SLAT(JHE)=0.
          WLAT(JHE)=2.*WLAT(JHE)
        ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  EQUALLY-SPACED LATITUDES EXCLUDING POLES
      ELSEIF(IDRT.EQ.256) THEN
        JH=JMAX/2
        JHE=(JMAX+1)/2
        JHO=JHE
        DLT=PI/JMAX
        SLAT(1)=1.
        DO J=1,JH
          SLAT(J)=COS((J-0.5)*DLT)
        ENDDO
        DO JS=1,JHO
          DO J=1,JHO
            AWORK(JS,J)=COS(2*(JS-1)*(J-0.5)*DLT)
          ENDDO
        ENDDO
        DO JS=1,JHO
          BWORK(JS)=-D1/(4*(JS-1)**2-1)
        ENDDO

        call ludcmp(awork,jho,jhe,ipvt,d)
        call lubksb(awork,jho,jhe,ipvt,bwork)

        WLAT(1)=0.
        DO J=1,JHO
          WLAT(J)=BWORK(J)
        ENDDO
CDIR$ IVDEP
        DO J=1,JH
          SLAT(JMAX+1-J)=-SLAT(J)
          WLAT(JMAX+1-J)=WLAT(J)
        ENDDO
        IF(JHE.GT.JH) THEN
          SLAT(JHE)=0.
          WLAT(JHE)=2.*WLAT(JHE)
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
