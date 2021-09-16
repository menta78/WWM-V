#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ST4_PRE (IP, WALOC, PHI, DPHIDN, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
      USE DATAPOOL
      USE W3SRC4MD
      IMPLICIT NONE

      INTEGER, INTENT(IN)        :: IP
      REAL(rkind), INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)

      REAL(rkind), INTENT(OUT)   :: PHI(NUMSIG,NUMDIR), DPHIDN(NUMSIG,NUMDIR)
      REAL(rkind), INTENT(OUT)   :: SSINE(NUMSIG,NUMDIR), DSSINE(NUMSIG,NUMDIR) 
      REAL(rkind), INTENT(OUT)   :: SSDS(NUMSIG,NUMDIR),DSSDS(NUMSIG,NUMDIR)
      REAL(rkind), INTENT(OUT)   :: SSNL4(NUMSIG,NUMDIR),DSSNL4(NUMSIG,NUMDIR)
      REAL(rkind), INTENT(OUT)   :: SSINL(NUMSIG,NUMDIR)

      INTEGER      :: IS, ID, ITH, IK, IS0

      REAL(rkind)  :: AWW3(NSPEC)
      REAL(rkind)  :: VDDS(NSPEC), VSDS(NSPEC), BRLAMBDA(NSPEC)
      REAL(rkind)  :: WN2(NUMSIG*NUMDIR), WHITECAP(ONE4)
      REAL(rkind)  :: VSIN(NSPEC), VDIN(NSPEC)

      REAL(rkind)  :: ETOT, FAVG, FMEANONE WNMEAN, AS, SUMWALOC, FAVGWS
      REAL(rkind)  :: TAUWAX, TAUWAY, AMAX, FPM, WINDONE, WINDTH
      REAL(rkind)  :: HS,SME0ONESMEONE,KME01,KMWAM,KMWAM2

      DO IS = ONE NUMSIG
        DO ID = ONE NUMDIR
          AWW3(ID + (IS-ONE * NUMDIR) = WALOC(IS,ID) * CG(IS,IP)
        END DO
      END DO

      DO IK=ONE NK
        WN2(ONE(IK-ONE*NTH) = WK(IK,IP)
      END DO

      DO IK=ONE NK
        IS0    = (IK-ONE*NTH
        DO ITH=2, NTH
          WN2(ITH+IS0) = WN2(ONEIS0)
        END DO
      END DO
!
! wind input
!
      TAUWX(IP)  = ZERO
      TAUWY(IP)  = ZERO               
      SSINL      = ZERO
      NUMSIG_HF(IP) = NUMSIG
      AS         = 0.
      BRLAMBDA   = ZERO

      IF (MESIN .GT. 0) THEN

        CALL SET_WIND( IP, WINDONE, WINDTH )
        CALL SET_FRICTION( IP, WALOC, WINDONE, WINDTH, FPM )
        LLWS=.TRUE.
#ifdef DEBUGSRC
        WRITE(740+myrank,*) 'ONE input value USTAR=', UFRIC(IP), ' USTDIR=', USTDIR(IP)
#endif
        CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), EMEAN(IP), FMEAN(IP), FMEANONE WNMEAN, AMAX, WINDONE, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS(IP))
#ifdef DEBUGSRC
        WRITE(740+myrank,*) 'ONE out value USTAR=', UFRIC(IP), ' USTDIR=', USTDIR(IP)
        WRITE(740+myrank,*) 'ONE out value EMEAN=', EMEAN(IP), ' FMEAN=', FMEAN(IP)
        WRITE(740+myrank,*) 'ONE out value FMEANONE', FMEAN1, ' WNMEAN=', WNMEAN
        WRITE(740+myrank,*) 'ONE out value CD=', CD(IP), ' Z0=', Z0(IP)
        WRITE(740+myrank,*) 'ONE out value ALPHA=', ALPHA_CH(IP), ' FMEANWS=', FMEANWS(IP)
#endif
        IF (EMEAN(IP) .LT. THR .AND. WINDONE .GT. THR) CALL SIN_LIN_CAV(IP,WINDTH,FPM,SSINL)
        CALL W3SIN4 ( IP, AWW3, CG(:,IP), WN2, WINDONE, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, VSIN, VDIN, LLWS, BRLAMBDA)
#ifdef DEBUGSRC
        WRITE(740+myrank,*) 'ONE WINDTH=', WINDTH, ' Z0=', Z0(IP), ' CD=', CD(IP)
        WRITE(740+myrank,*) 'ONE UFRIC=', UFRIC(IP), 'WINDONE=', WIND10, ' RHOAW=', RHOAW
        WRITE(740+myrank,*) 'ONE TAUWX=', TAUWX(IP), ' TAUWY=', TAUWY(IP)
        WRITE(740+myrank,*) 'ONE TAUWAX=', TAUWAX, ' TAUWAY=', TAUWAY
        WRITE(740+myrank,*) 'ONE W3SIN4min/max/sum(VSIN)=', minval(VSIN), maxval(VSIN), sum(VSIN)
        WRITE(740+myrank,*) 'ONE W3SIN4min/max/sum(VDIN)=', minval(VDIN), maxval(VDIN), sum(VDIN)
#endif
        CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), EMEAN(IP), FMEAN(IP), FMEANONE WNMEAN, AMAX, WINDONE, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS(IP))  
        CALL W3SIN4 ( IP, AWW3, CG(:,IP), WN2, WINDONE, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, VSIN, VDIN, LLWS, BRLAMBDA)
#ifdef DEBUGSRC
        WRITE(740+myrank,*) '2: W3SIN4min/max/sum(VSIN)=', minval(VSIN), maxval(VSIN), sum(VSIN)
        WRITE(740+myrank,*) '2: W3SIN4min/max/sum(VDIN)=', minval(VDIN), maxval(VDIN), sum(VDIN)
#endif
        CALL CONVERT_VS_VD_WWM(IP, VSIN, VDIN, SSINE, DSSINE)
      ENDIF

      IF (MESNL .GT. 0) THEN
         CALL MEAN_WAVE_PARAMETER(IP,WALOC,HS,ETOT,SME0ONESMEONE,KME01,KMWAM,KMWAM2)
         CALL DIASNL4WW3(IP, KMWAM, WALOC, SSNL4, DSSNL4)
      END IF

      IF (MESDS .GT. 0) THEN
        CALL W3SDS4(AWW3,WK(:,IP),CG(:,IP),UFRIC(IP),USTDIR(IP),DEP(IP),VSDS,VDDS,BRLAMBDA,WHITECAP)

#ifdef DEBUGSRC
        WRITE(740+myrank,*) '2: W3SDS4min/max/sum(VSDS)=', minval(VSDS), maxval(VSDS), sum(VSDS)
        WRITE(740+myrank,*) '2: W3SDS4min/max/sum(VDDS)=', minval(VDDS), maxval(VDDS), sum(VDDS)
#endif
        CALL CONVERT_VS_VD_WWM(IP, VSDS, VDDS, SSDS, DSSDS)
      ENDIF
!
      PHI    = SSINL + SSINE  + SSNL4  + SSDS
      DPHIDN =         DSSINE + DSSNL4 + DSSDS
!
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ST4_POST (IP, WALOC, SSINE, DSSINE, SSDS, DSSDS, SSINL)
        USE DATAPOOL
        USE W3SRC4MD
        IMPLICIT NONE
       
        INTEGER, INTENT(IN)        :: IP
        REAL(rkind), INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)
       
        REAL(rkind), INTENT(OUT)   :: SSINE(NUMSIG,NUMDIR),DSSINE(NUMSIG,NUMDIR)
        REAL(rkind), INTENT(OUT)   :: SSDS(NUMSIG,NUMDIR),DSSDS(NUMSIG,NUMDIR)
        REAL(rkind), INTENT(OUT)   :: SSINL(NUMSIG,NUMDIR)

        INTEGER                    :: IS, ID, IK, ITH, ITH2, IS0

        REAL(rkind)                :: PHI(NUMSIG,NUMDIR), DPHIDN(NUMSIG,NUMDIR)
        REAL(rkind)                :: AWW3(NSPEC), WN2(NUMSIG*NUMDIR), BRLAMBDA(NSPEC)
        REAL(rkind)                :: DPHIDNONE(NSPEC), PHIONE(NSPEC), TMP_DS(NUMSIG)

        REAL(rkind)                :: ETOT, FAVG, FMEANONE WNMEAN, AS, FAVGWS
        REAL(rkind)                :: TAUWAX, TAUWAY, AMAX, WINDONE, WINDTH
        REAL(rkind)                :: WHITECAP(ONE4), SUMWALOC, FPM

        DO IS = ONE NUMSIG
          DO ID = ONE NUMDIR
            AWW3(ID + (IS-ONE * NUMDIR) = WALOC(IS,ID) * CG(IS,IP)
          END DO
        END DO

        DO IK=ONE NK
          WN2(ONE(IK-ONE*NTH) = WK(IK,IP)
        END DO
        DO IK=ONE NK
          IS0    = (IK-ONE*NTH
          DO ITH=2, NTH
            WN2(ITH+IS0) = WN2(ONEIS0)
          END DO
        END DO
!
! wind input
!
        AS      = 0.
        NUMSIG_HF(IP) = NUMSIG
        CALL SET_WIND( IP, WINDONE, WINDTH )
        CALL SET_FRICTION( IP, WALOC, WINDONE, WINDTH, FPM )
        CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), EMEAN(IP), FMEAN(IP), FMEANONE WNMEAN, AMAX, WINDONE, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS(IP))
        IF (EMEAN(IP) .LT. THR .AND. WINDONE .GT. THR) CALL SIN_LIN_CAV(IP,WINDTH,FPM,SSINL)
        CALL W3SIN4 ( IP, AWW3, CG(:,IP), WN2,  WINDONE, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, PHIONE, DPHIDN1D, LLWS, BRLAMBDA)
        CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), EMEAN(IP), FMEAN(IP), FMEANONE WNMEAN, AMAX, WINDONE, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS(IP))
        CALL CONVERT_VS_VD_WWM(IP, PHIONE, DPHIDNONE, SSINE, DSSINE)
!
! dissipation 
!
        CALL W3SDS4(AWW3,WK(:,IP),CG(:,IP),UFRIC(IP),USTDIR(IP),DEP(IP),PHIONE,DPHIDNONE,BRLAMBDA,WHITECAP)
        CALL CONVERT_VS_VD_WWM(IP,PHIONE,DPHIDNONE,SSDS,DSSDS)
!
! missing high freq. tail contribution -> 2do
!

!
! adding the fluxes from waves to ocean ...
!
        WHITECAP(3)=0.
        HSTOT=0.
        DO IK=IKSONE NK
          FACTOR = DDEN(IK)/CGONEIK)                    !Jacobian to get energy in band
          FACTOR2= FACTOR*GRAV*WNONEIK)/SIG(IK)         ! coefficient to get momentum
          ! Wave direction is "direction to"
          ! therefore there is a PLUS sign for the stress
            IS   = (IK-ONE*NTH + ITH
            COSI(ONE=ECOS(IS)
            COSI(2)=ESIN(IS)
            PHIAW = PHIAW  + VSIN(IS) * DT * FACTOR / MAX ( ONE , (ONE-HDT*VDIN(IS))) ! semi-implict integration scheme
            PIBBL = PHIBBL - VSBT(IS) * DT * FACTOR / MAX ( ONE , (ONE-HDT*VDBT(IS))) ! semi-implict integration scheme
            PHINL = PHINL  + VSNL(IS) * DT * FACTOR / MAX ( ONE , (ONE-HDT*VDNL(IS))) ! semi-implict integration scheme
            IF (VSIN(IS).GT.0.) THEN
              WHITECAP(3) = WHITECAP(3) + SPEC(IS)  * FACTOR
            ELSE
            ! computes the upward energy flux (counted > 0 upward)
              CHARN = CHARN - (VSIN(IS))* DT * FACTOR / MAX ( ONE , (ONE-HDT*VDIN(IS))) ! semi-implict integration scheme
            END IF
          HSTOT = HSTOT + SPEC(IS) * FACTOR
          END DO
        END DO
        WHITECAP(3)=4.*SQRT(WHITECAP(3))
        HSTOT=4.*SQRT(HSTOT)
        TAUWIX= TAUWIX+ TAUWX * DRAT *DT
        TAUWIY= TAUWIY+ TAUWY * DRAT *DT
        TAUWNX= TAUWNX+ TAUWAX * DRAT *DT
        TAUWNY= TAUWNY+ TAUWAY * DRAT *DT

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
