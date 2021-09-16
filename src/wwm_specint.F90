#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SEMI_IMPLICIT_INTEGRATION(IP,DT,WACOLD,WACNEW)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IP
      REAL(rkind), INTENT(IN) :: DT
      REAL(rkind), INTENT(INOUT) :: WACOLD(NUMSIG,NUMDIR)
      REAL(rkind), INTENT(OUT)   :: WACNEW(NUMSIG,NUMDIR)

      REAL(rkind)              :: PHI(NUMSIG,NUMDIR), DPHIDN(NUMSIG,NUMDIR)
      REAL(rkind)              :: SSLIM(NUMSIG,NUMDIR), SSBRL(NUMSIG,NUMDIR)

      REAL(rkind)              :: NEWDAC, MAXDAC(NUMSIG)

      INTEGER       :: IS, ID

      CALL COMPUTE_PHI_DPHI(IP,WACOLD,PHI,DPHIDN)

!      WRITE(*,'(A20,10F20.10)') 'PHI,DPHIDN', SUM(PHI), SUM(DPHIDN) 
!      WRITE(*,'(A20,10F20.10)') 'BEFORE INTEGRATION', SUM(WACOLD), SUM(WACNEW)

      DO IS = 1, NUMSIG
        DO ID = 1, NUMDIR
          NEWDAC = PHI(IS,ID) * DT / (ONE-DT*MIN(ZERO,DPHIDN(IS,ID)))
          WACNEW(IS,ID) = MAX( ZERO, WACOLD(IS,ID) + NEWDAC )
        END DO
      END DO

!      WRITE(*,'(A20,10F20.10)') 'AFTER INTEGRATION', SUM(WACOLD), SUM(WACNEW)

      IF (LMAXETOT) CALL BREAKING_LIMITER_LOCAL(IP,WACNEW,SSBRL)

      CALL POST_INTEGRATION(IP,WACOLD,WACNEW)

!      WRITE(*,'(A20,10F20.10)') 'AFTER POST', SUM(WACOLD), SUM(WACNEW)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_PHI_DPHI(IP,WALOC,PHI,DPHIDN)
      USE DATAPOOL 
      USE W3UOSTMD
      IMPLICIT NONE
      INTEGER, INTENT(IN)      :: IP
      INTEGER                  :: IS, ID
      REAL(rkind), INTENT(IN)  :: WALOC(NUMSIG,NUMDIR)
      REAL(rkind), INTENT(OUT) :: PHI(NUMSIG,NUMDIR), DPHIDN(NUMSIG,NUMDIR)
 
      REAL(rkind)              :: HS,TM01,TM02,TM10,KLM,WLM
      REAL(rkind)              :: MAXDAC(NUMSIG), NEWDAC, RATIO

      REAL(rkind) :: SSINL(NUMSIG,NUMDIR)
      REAL(rkind) :: SSINE(NUMSIG,NUMDIR),DSSINE(NUMSIG,NUMDIR)
      REAL(rkind) :: SSDS(NUMSIG,NUMDIR),DSSDS(NUMSIG,NUMDIR)
      REAL(rkind) :: SSNL3(NUMSIG,NUMDIR),DSSNL3(NUMSIG,NUMDIR)
      REAL(rkind) :: SSNL4(NUMSIG,NUMDIR),DSSNL4(NUMSIG,NUMDIR)
      REAL(rkind) :: SSBR(NUMSIG,NUMDIR),DSSBR(NUMSIG,NUMDIR)
      REAL(rkind) :: SSBRL(NUMSIG,NUMDIR)
      REAL(rkind) :: SSBF(NUMSIG,NUMDIR),DSSBF(NUMSIG,NUMDIR)
      REAL(rkind) :: SSUO(NUMSIG,NUMDIR),SDUO(NUMSIG,NUMDIR)
      REAL(rkind) :: SSIC(NUMSIG,NUMDIR),SDIC(NUMSIG,NUMDIR)

      LOGICAL     :: LMERGE = .FALSE. 
      LOGICAL     :: LRELAX = .TRUE. 

#ifdef DEBUG
      REAL(rkind)              :: SSINL_WW3(NUMSIG,NUMDIR), SSINE_WW3(NUMSIG,NUMDIR)
      REAL(rkind)              :: SSBRL_WW3(NUMSIG,NUMDIR), SSDS_WW3(NUMSIG,NUMDIR)
      REAL(rkind)              :: SSNL4_WW3(NUMSIG,NUMDIR), SSBF_WW3(NUMSIG,NUMDIR)
      REAL(rkind)              :: SSBR_WW3(NUMSIG,NUMDIR)
#endif

      DPHIDN = ZERO; PHI    = ZERO
      SSINE  = ZERO; DSSINE = ZERO
      SSDS   = ZERO; DSSDS  = ZERO
      SSNL3  = ZERO; DSSNL3 = ZERO
      SSNL4  = ZERO; DSSNL4 = ZERO
      SSBR   = ZERO; DSSBR  = ZERO
      SSBF   = ZERO; DSSBF  = ZERO
      SSUO   = ZERO; SDUO   = ZERO
      SSINL  = ZERO

#ifdef DEBUG
      IF (IP .eq. TESTNODE) THEN
         WRITE(740+myrank,*) 'Before integration'
         WRITE(740+myrank,*) 'Before integration', IP, IPLG(IP) 
         WRITE(740+myrank,*) 'sum(WALOC)=', sum(WALOC)
         CALL MEAN_PARAMETER(IP,WALOC,NUMSIG,HS,TM01,TM02,TM10,KLM,WLM)
         WRITE(740+myrank,*) 'HS=', HS, ' TM01=', TM01
      END IF
#endif
!
      CALL DEEP_WATER(IP, WALOC, PHI, DPHIDN, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
!
#ifdef DEBUG
      IF (IP .eq. TESTNODE) THEN
         WRITE(740+myrank,*) 'sum(SSINE)=', sum(SSINE), ' sum(DSSINE)=', sum(DSSINE)
      END IF
#endif
!
      IF (MESUO .GT. 0) THEN
         CALL UOST_SRCTRM(IP, WALOC, CG(:,IP), SSUO, SDUO)     
         PHI = PHI + SSUO
         DPHIDN = DPHIDN + SDUO
      ENDIF

      IF (MESIC .GT. 0) THEN
         CALL ICEUODIS_SRCTRM(IP, WALOC, SSIC, SDIC)
         PHI = PHI + SSIC
         DPHIDN = DPHIDN + SDIC
      ENDIF

!
!      WRITE(*,*) 'SUMS', SUM(PHI), SUM(DPHIDN)
!
      IF (ISHALLOW(IP) .EQ. 1) THEN
        CALL SHALLOW_WATER(IP, WALOC, LMERGE, PHI, DPHIDN, SSBR, DSSBR, SSBF, DSSBF, SSNL3, DSSNL3)
      ENDIF

      IF (LRELAX .and. ICOMP .GE. 2) THEN
        DO ID = 1, NUMDIR
          DO IS = 1, NUMSIG
            IF (SSNL3(IS,ID) .gt. ZERO) THEN 
              SSNL3(IS,ID)  =   2 * SSNL3(IS,ID)
              DSSNL3(IS,ID) = - DSSNL3(IS,ID) 
            ELSE
              SSNL3(IS,ID) =       SSNL3(IS,ID)
              DSSNL3(IS,ID) = 2 * DSSNL3(IS,ID)
            ENDIF
          END DO
        END DO
        DO ID = 1, NUMDIR
          DO IS = 1, NUMSIG
            SSBR(IS,ID)  =      SSBR(IS,ID)
            DSSBR(IS,ID) = 2 * DSSBR(IS,ID)
          END DO
        END DO
      ENDIF

      IF (MELIM .GT. 0) THEN
        CALL GET_MAXDAC(IP,MAXDAC)
        DO ID = 1, NUMDIR
          DO IS = 1, NUMSIG
            NEWDAC = PHI(IS,ID) * DT4A / (ONE-DT4A*MIN(ZERO,DPHIDN(IS,ID)))
            IF (MAXDAC(IS) .GT. ZERO) THEN
              RATIO = ONE/MAX(ONE,ABS(NEWDAC/MAXDAC(IS)))
            ELSE
              RATIO = ONE
            ENDIF
!            WRITE(*,*) IS, ID, RATIO, MAXDAC(IS), NEWDAC, PHI(IS,ID), DPHIDN(IS,ID), SUM(PHI), SUM(DPHIDN)
            PHI(IS,ID) = RATIO * PHI(IS,ID)
          END DO
        END DO
        DO ID = 1, NUMDIR
          DO IS = 1, NUMSIG
            NEWDAC = SSNL3(IS,ID) * DT4A / (ONE-DT4A*MIN(ZERO,DSSNL3(IS,ID)))
            IF (MAXDAC(IS) .GT. ZERO) THEN
              RATIO = ONE/MAX(ONE,ABS(NEWDAC/(100*MAXDAC(IS))))
            ELSE
              RATIO = ONE
            ENDIF
            SSNL3(IS,ID)  = RATIO * SSNL3(IS,ID)
            DSSNL3(IS,ID) = RATIO * DSSNL3(IS,ID)
            !IF (ABS(DSSNL3(IS,ID)) .gt. zero .or. ABS(SSNL3(IS,ID)) .gt. zero)  THEN
            !  IF (RATIO .LT. ONE) THEN
            !    WRITE(*,'(A20, 2I10, 10F20.10)') 'RATIO AND MAXLIM', IS, ID, RATIO, MAXDAC(IS)
            !    WRITE(*,'(A20, 2I10, 10F20.10)') 'OLD and NEW WAVEACTION', IS, ID, WALOC(IS,ID), NEWDAC
            !    WRITE(*,'(A20, 2I10, 10E20.10)') 'RIGHT and DIAGONAL', IS, ID, SSNL3(IS,ID), DSSNL3(IS,ID)
            !  ENDIF
            !ENDIF
          END DO
        END DO
      ENDIF

      PHI    = PHI     + SSBR  + SSBF  + SSNL3
      DPHIDN = DPHIDN + DSSBR + DSSBF + DSSNL3


#ifdef DEBUG
      IF (IP .eq. TESTNODE) THEN
         WRITE(740+myrank,*) 'After integration ISOURCE=', ISOURCE
         WRITE(740+myrank,*) 'sum(WALOC)=', sum(WALOC)
         CALL MEAN_PARAMETER(IP,WALOC,NUMSIG,HS,TM01,TM02,TM10,KLM,WLM)
         WRITE(740+myrank,*) 'HS=', HS, ' TM01=', TM01
         CALL CONVERT_VS_WWM_TO_WW3(IP, SSINL, SSINL_WW3)
         WRITE(740+myrank,*) 'WW3 : LINEAR INPUT =', SUM(SSINL_WW3), MINVAL(SSINL_WW3), MAXVAL(SSINL_WW3)
         CALL CONVERT_VS_WWM_TO_WW3(IP, SSINE, SSINE_WW3)
         WRITE(740+myrank,*) 'WW3 : EXP. INPUT   =', SUM(SSINE_WW3), MINVAL(SSINE_WW3), MAXVAL(SSINE_WW3)
         CALL CONVERT_VS_WWM_TO_WW3(IP, SSDS, SSDS_WW3)
         WRITE(740+myrank,*) 'WW3 : WHITECAP     =', SUM(SSDS_WW3), MINVAL(SSDS_WW3), MAXVAL(SSDS_WW3)
         CALL CONVERT_VS_WWM_TO_WW3(IP, SSNL4, SSNL4_WW3)
         WRITE(740+myrank,*) 'WW3 : SNL4         =', SUM(SSNL4_WW3), MINVAL(SSNL4_WW3), MAXVAL(SSNL4_WW3)
         CALL CONVERT_VS_WWM_TO_WW3(IP, SSBF, SSBF_WW3)
         WRITE(740+myrank,*) 'WW3 : BOT. FRIC.   =', SUM(SSBF_WW3), MINVAL(SSBF_WW3), MAXVAL(SSBF_WW3)
         CALL CONVERT_VS_WWM_TO_WW3(IP, SSBR, SSBR_WW3)
         WRITE(740+myrank,*) 'WW3 : BREAKING     =', SUM(SSBR_WW3), MINVAL(SSBR_WW3), MAXVAL(SSBR_WW3)
         
         WRITE(740+myrank,*) 'WAVE ACTION      =', SUM(WALOC), MINVAL(WALOC), MAXVAL(WALOC)
         WRITE(740+myrank,*) 'LINEAR INPUT     =', SUM(SSINL), MINVAL(SSINL), MAXVAL(SSINL)
         WRITE(740+myrank,*) 'BREAKING LIMITER =', SUM(SSBRL), MINVAL(SSBRL), MAXVAL(SSBRL)
         WRITE(740+myrank,*) 'EXP INPUT(SSINE) =', SUM(SSINE),  MINVAL(SSINE),  MAXVAL(SSINE)
         WRITE(740+myrank,*) 'EXP INPUT(DSSINE)=', SUM(DSSINE), MINVAL(DSSINE), MAXVAL(DSSINE)
         WRITE(740+myrank,*) 'WHITECAP(SSDS)   =', SUM(SSDS),  MINVAL(SSDS),  MAXVAL(SSDS)
         WRITE(740+myrank,*) 'WHITECAP(DSSDS)  =', SUM(DSSDS), MINVAL(DSSDS), MAXVAL(DSSDS)
         WRITE(740+myrank,*) 'SNL4(SSNL4)      =', SUM(SSNL4), MINVAL(SSNL4), MAXVAL(SSNL4)
         WRITE(740+myrank,*) 'SNL4(DSSNL4)     =', SUM(DSSNL4), MINVAL(DSSNL4), MAXVAL(DSSNL4)
         WRITE(740+myrank,*) 'BOT. FRIC(SSBF)  =', SUM(SSBF),  MINVAL(SSBF),  MAXVAL(SSBF)
         WRITE(740+myrank,*) 'BOT. FRIC(DSSBF) =', SUM(DSSBF), MINVAL(DSSBF), MAXVAL(DSSBF)
         WRITE(740+myrank,*) 'BREAKING(SSBR)   =', SUM(SSBR),  MINVAL(SSBR),  MAXVAL(SSBR)
         WRITE(740+myrank,*) 'BREAKING(DSSBR)  =', SUM(DSSBR), MINVAL(DSSBR), MAXVAL(DSSBR)
         WRITE(740+myrank,*) 'TOTAL SRC(PHI)   =', SUM(PHI), MINVAL(PHI), MAXVAL(PHI)
         WRITE(740+myrank,*) 'TOTAL SRC(DPHIDN)=', SUM(DPHIDN), MINVAL(DPHIDN), MAXVAL(DPHIDN)
      END IF
#endif

#ifdef DEBUG_SOURCE_TERM
      WRITE(740+myrank,'(A20,6E20.10)') 'WAVE ACTION', SUM(WALOC), MINVAL(WALOC), MAXVAL(WALOC)
      WRITE(740+myrank,'(A20,6E20.10)') 'LINEAR INPUT', SUM(SSINL), MINVAL(SSINL), MAXVAL(SSINL)
      WRITE(740+myrank,'(A20,6E20.10)') 'BREAKING LIMITER', SUM(SSBRL), MINVAL(SSBRL), MAXVAL(SSBRL)
      WRITE(740+myrank,'(A20,6E20.10)') 'EXP INPUT', SUM(SSINE), SUM(DSSINE), MINVAL(SSINE), MAXVAL(SSINE), MINVAL(DSSINE), MAXVAL(DSSINE)
      WRITE(740+myrank,'(A20,6E20.10)') 'WHITECAP', SUM(SSDS), SUM(DSSDS), MINVAL(SSDS), MAXVAL(SSDS), MINVAL(DSSDS), MAXVAL(DSSDS)
      WRITE(740+myrank,'(A20,6E20.10)') 'SNL4', SUM(SSNL4), SUM(DSSNL4), MINVAL(SSNL4), MAXVAL(SSNL4), MINVAL(DSSNL4), MAXVAL(DSSNL4)
      WRITE(740+myrank,'(A20,6E20.10)') 'BOTTOM FRICTION', SUM(SSBF), SUM(DSSBF), MINVAL(SSBF), MAXVAL(SSBF), MINVAL(DSSBF), MAXVAL(DSSBF)
      WRITE(740+myrank,'(A20,6E20.10)') 'BREAKING', SUM(SSBR), SUM(DSSBR), MINVAL(SSBR), MAXVAL(SSBR), MINVAL(DSSBR), MAXVAL(DSSBR)
      WRITE(740+myrank,'(A20,6E20.10)') 'TOTAL SOURCE TERMS', SUM(PHI), SUM(DPHIDN), MINVAL(PHI), MAXVAL(PHI), MINVAL(DPHIDN), MAXVAL(DPHIDN)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CONVERT_VS_WWM_TO_WW3(IP, VS_WWM, VS_WW3)
      USE DATAPOOL
      IMPLICIT NONE
      integer IP
      REAL(rkind), intent(in) :: VS_WWM(NUMSIG,NUMDIR)
      REAL(rkind), intent(out) :: VS_WW3(NUMSIG,NUMDIR)
      INTEGER ID,IS
      DO ID=1,NUMDIR
        DO IS=1,NUMSIG
          VS_WW3(IS,ID) = CG(IS,IP) * VS_WWM(IS,ID)
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE POST_INTEGRATION(IP,WACOLD,WACNEW)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)       :: IP
         REAL(rkind),INTENT(IN)    :: WACOLD(NUMSIG,NUMDIR)
         REAL(rkind),INTENT(OUT)   :: WACNEW(NUMSIG,NUMDIR)

         IF (ISOURCE == 1) THEN
           !CALL ST4_POST(IP,WACOLD, WACNEW)
         ELSE IF (ISOURCE == 2) THEN
           CALL ECMWF_POST(IP, WACNEW)
         ELSE IF (ISOURCE == 3) THEN
         ELSE IF (ISOURCE == 4) THEN
           CALL ST6_POST(IP,WACOLD, WACNEW)
         ENDIF
      END SUBROUTINE POST_INTEGRATION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SOURCES_EXPLICIT
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP
         REAL(rkind) :: WACOLD(NUMSIG,NUMDIR), WACNEW(NUMSIG,NUMDIR)
         REAL(rkind) :: SSLIM(NUMSIG,NUMDIR), SSBRL(NUMSIG,NUMDIR)
         REAL(rkind) :: MAXDAC(NUMSIG)

         DO IP = 1, MNP
           IF (IOBDP(IP) .GT. 0) THEN ! H .gt. DMIN
             WACOLD = AC2(:,:,IP)
             WACNEW = ZERO
             IF (LSOUBOUND .AND. IOBP(IP) .NE. 2) THEN ! CALL ALWAYS 
               CALL SEMI_IMPLICIT_INTEGRATION(IP,DT4S,WACOLD,WACNEW)
             ELSE IF (IOBP(IP) .EQ. 0 .AND. .NOT. LSOUBOUND) THEN ! CALL ONLY FOR NON BOUNDARY POINTS 
               CALL SEMI_IMPLICIT_INTEGRATION(IP,DT4S,WACOLD,WACNEW)
             ELSE
               WACNEW = WACOLD
             ENDIF 
             AC2(:,:,IP) = WACNEW
           ELSE
             AC2(:,:,IP) = ZERO
           ENDIF
         ENDDO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SOURCES_IMPLICIT
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER     :: IP, IS, ID
         REAL(rkind) :: WACOLD(NUMSIG,NUMDIR),PHI(NUMSIG,NUMDIR),DPHIDN(NUMSIG,NUMDIR)
         REAL(rkind) :: NEWDAC, MAXDAC(NUMSIG), RATIO

         REAL(rkind) :: SSINL(NUMSIG,NUMDIR)
         REAL(rkind) :: SSINE(NUMSIG,NUMDIR),DSSINE(NUMSIG,NUMDIR)
         REAL(rkind) :: SSDS(NUMSIG,NUMDIR),DSSDS(NUMSIG,NUMDIR)
         REAL(rkind) :: SSNL3(NUMSIG,NUMDIR),DSSNL3(NUMSIG,NUMDIR)
         REAL(rkind) :: SSNL4(NUMSIG,NUMDIR),DSSNL4(NUMSIG,NUMDIR)
         REAL(rkind) :: SSBR(NUMSIG,NUMDIR),DSSBR(NUMSIG,NUMDIR)
         REAL(rkind) :: SSBRL(NUMSIG,NUMDIR)
         REAL(rkind) :: SSBF(NUMSIG,NUMDIR),DSSBF(NUMSIG,NUMDIR)

         DO IP = 1, MNP
           IF (IOBDP(IP) .GT. 0) THEN ! H .gt. DMIN
             PHI = ZERO
             DPHIDN = ZERO
             WACOLD = AC2(:,:,IP)
             IF (LSOUBOUND  .AND. IOBP(IP) .NE. 2) THEN ! CALL ALWAYS
               CALL COMPUTE_PHI_DPHI(IP,WACOLD,PHI,DPHIDN)
             ELSE IF (IOBP(IP) .EQ. 0 .AND. .NOT. LSOUBOUND) THEN ! CALL ONLY FOR NON BOUNDARY POINTS
               CALL COMPUTE_PHI_DPHI(IP,WACOLD,PHI,DPHIDN)
             ENDIF
             PHIA(:,:,IP)    = PHI   
             DPHIDNA(:,:,IP) = DPHIDN
           ELSE
             PHIA(:,:,IP)    = ZERO
             DPHIDNA(:,:,IP) = ZERO
           ENDIF
         ENDDO

#ifdef DEBUG_SOURCE_TERM
         WRITE(DBG%FHNDL,*) 'SOURCES_IMPLICIT', SUM(PHIA), SUM(DPHIDNA)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
         SUBROUTINE GET_MAXDAC(IP,MAXDAC)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: IP
         INTEGER                    :: IS, ID
         REAL(rkind), INTENT(OUT)   :: MAXDAC(NUMSIG) 
         REAL(rkind)                :: DELFL(NUMSIG)
         REAL(rkind)                :: USFM, PHILMAXDAC

         DELFL  = COFRM4*DT4S

         DO IS = 1, NUMSIG
           PHILMAXDAC = 0.0081*0.1/(TWO*SPSIG(IS)*WK(IS,IP)*WK(IS,IP)*WK(IS,IP)*CG(IS,IP)) ! Phillips limiter following Komen et al. 
           IF (MELIM == 1 ) THEN  ! Phillips 
             MAXDAC(IS) = PHILMAXDAC
           ELSE IF (MELIM == 2) THEN ! Hersbach & Janssen 
             USFM       = UFRIC(IP)*MAX(FMEANWS(IP),FMEAN(IP))
             MAXDAC(IS) = USFM*COFRM4(IS)*DT4S/PI2/SPSIG(IS)
           ELSE IF (MELIM == 3) THEN ! Roland, 2018
             USFM       = UFRIC(IP)*MAX(FMEANWS(IP),FMEAN(IP))
             MAXDAC(IS) = MAX(PHILMAXDAC,USFM*DELFL(IS)/PI2/SPSIG(IS))
           END IF
         END DO

         END SUBROUTINE GET_MAXDAC 
!**********************************************************************
!*                                                                    *
!**********************************************************************
         SUBROUTINE ACTION_LIMITER_LOCAL(MAXDAC,WAOLD,WANEW,SSLIM)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER                    :: IS, ID
         REAL(rkind), INTENT(INOUT) :: WANEW(NUMSIG,NUMDIR), SSLIM(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(IN)    :: WAOLD(NUMSIG,NUMDIR), MAXDAC(NUMSIG)
         REAL(rkind)                :: NEWDAC, OLDAC, NEWAC, NEWDACL

!         WRITE(*,*) 'BEFORE', SUM(WAOLD), SUM(WANEW), SUM(WAOLD)-SUM(WANEW)

         DO ID = 1, NUMDIR
           DO IS = 1, NUMSIG
             NEWAC  = WANEW(IS,ID)
             OLDAC  = WAOLD(IS,ID)
             NEWDAC = NEWAC - OLDAC
             !IF (OLDAC .GT. 0. .or. NEWAC .gt. 0) WRITE(*,*) IS, ID, OLDAC, NEWAC, NEWDAC, MAXDAC(IS)
             NEWDAC = SIGN(MIN(MAXDAC(IS),ABS(NEWDAC)),NEWDAC)
             WANEW(IS,ID) = OLDAC + NEWDAC 
           END DO
         END DO

!         WRITE(*,*) 'AFTER', SUM(WAOLD), SUM(WANEW), SUM(WAOLD)-SUM(WANEW)

         END SUBROUTINE ACTION_LIMITER_LOCAL
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ACTION_LIMITER_GLOBAL(ACOLD,ACNEW)
        USE DATAPOOL
        IMPLICIT NONE
        REAL(rkind), INTENT(IN)  :: ACOLD(NUMSIG,NUMDIR,MNP)
        REAL(rkind), INTENT(OUT) :: ACNEW(NUMSIG,NUMDIR,MNP)
        REAL(rkind)              :: SSLIM(NUMSIG,NUMDIR)

        REAL(rkind)              :: MAXDAC(NUMSIG)
        REAL(rkind)              :: ACNEWLOC(NUMSIG,NUMDIR),ACOLDLOC(NUMSIG,NUMDIR)
        INTEGER                  :: IP

        DO IP = 1, NP_RES
          ACOLDLOC = ACOLD(:,:,IP)
          ACNEWLOC = ACNEW(:,:,IP)
          CALL GET_MAXDAC(IP,MAXDAC)
          CALL ACTION_LIMITER_LOCAL(MAXDAC,ACOLDLOC,ACNEWLOC,SSLIM)
          ACNEW(:,:,IP) = ACNEWLOC
        ENDDO

      END SUBROUTINE ACTION_LIMITER_GLOBAL
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE BREAKING_LIMITER_GLOBAL
        USE DATAPOOL
        IMPLICIT NONE
        REAL(rkind)              :: SSBRL(NUMSIG,NUMDIR)
        REAL(rkind)              :: MAXDAC(NUMSIG)
        REAL(rkind)              :: ACNEWLOC(NUMSIG,NUMDIR)
        INTEGER                  :: IP

        DO IP = 1, MNP
          ACNEWLOC = AC2(:,:,IP)
          CALL BREAKING_LIMITER_LOCAL(IP,ACNEWLOC,SSBRL)
          AC2(:,:,IP) = ACNEWLOC
        ENDDO

      END SUBROUTINE BREAKING_LIMITER_GLOBAL
!**********************************************************************
!*                                                                    *
!**********************************************************************
         SUBROUTINE BREAKING_LIMITER_LOCAL(IP,WALOC,SSBRL)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)  :: IP

         REAL(rkind), INTENT(INOUT)  :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT)    :: SSBRL(NUMSIG,NUMDIR)

         REAL(rkind)                 :: HS
         REAL(rkind)                 :: EMAX, RATIO, ETOT
         REAL(rkind)                 :: DINTSPEC

!AR: here we can save time to ommit this calculation of the mean value!!!
         ETOT = DINTSPEC(WALOC)
         HS   = 4.*SQRT(ETOT)
         !IF (LMONO_IN) HMAX(IP) = HMAX(IP) * SQRT(2.)
         EMAX = 1./16. * (HMAX(IP))**2

         IF (ETOT .GT. EMAX .AND. ETOT .GT. THR) THEN
           RATIO = EMAX/ETOT
           SSBRL = WALOC - RATIO * WALOC
           WALOC = RATIO * WALOC
         END IF

         END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
         SUBROUTINE RESCALE_SPECTRUM
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP, IS, ID
         REAL(rkind)    :: ETOTAL, EPOSIT
         REAL(rkind)    :: FACTOR
         LOGICAL :: ENEG

         DO IP = 1, MNP
            DO IS = 1, NUMSIG
               ETOTAL = ZERO 
               EPOSIT = ZERO 
               ENEG   = .FALSE.
               DO ID = 1, NUMDIR
                  ETOTAL = ETOTAL + AC2(IS,ID,IP)
                  IF (AC2(IS,ID,IP) > ZERO) THEN
                     EPOSIT = EPOSIT + AC2(IS,ID,IP)
                  ELSE
                     ENEG = .TRUE.
                  END IF
               END DO
               IF (ENEG) THEN
                  IF (EPOSIT .GT. VERYSMALL) THEN
                    FACTOR = ETOTAL/EPOSIT
                  ELSE 
                    FACTOR = 0.
                  END IF
                  DO ID = 1, NUMDIR
                     IF (AC2(IS,ID,IP) < ZERO) AC2(IS,ID,IP) = ZERO 
                     IF (FACTOR >= ZERO)  AC2(IS,ID,IP) = AC2(IS,ID,IP)*FACTOR
                     AC2(IS,ID,IP) = MAX(zero,AC2(IS,ID,IP))
                  END DO
               END IF
            END DO
         END DO
         END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
         SUBROUTINE SETSHALLOW
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP

         DO IP = 1, MNP
           IF (WK(1,IP)*DEP(IP) .LT. PI) THEN
             ISHALLOW(IP) = 1
           ELSE
             ISHALLOW(IP) = 0
           END IF
         END DO
         END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
