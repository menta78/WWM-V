      MODULE WWM_PARTMD
!     This module was originally implemented in WW3 by 
!     B. Tracy and H. Tolman, with contributions by M. Szyszka and C. Bunney
!     
!     It contains an implementation of a "watershed by flooding" algorithm for
!     the partition of the directional spectrum.
!
!     Lorenzo Mentaschi adapted it to WWM
!
      USE DATAPOOL, ONLY: rkind
      PUBLIC
!
      INTEGER, PUBLIC, PARAMETER   :: NPOUTMAX = 3 ! maximum number of stored peaks
      INTEGER, PUBLIC, PARAMETER   :: DIMP = 6 ! number of parameters in partition (setting a high value)
      REAL(rkind), PUBLIC                :: SPART_PARAMS(DIMP,NPOUTMAX) ! array used to store the output partition parameters

      INTEGER, PRIVATE              :: MK = -1, MTH = -1
      INTEGER, ALLOCATABLE, PRIVATE :: NEIGH(:,:)

      INTEGER, PRIVATE              :: IHMAX = 50 ! number of levels used in the watershed algorithm
      INTEGER, PRIVATE              :: HSPMIN = 0.05 ! minimum Hs of a partition (for values below this threshold the peak parameters are put to UNDEF)
      
      CONTAINS



      SUBROUTINE WWM_PART ( IP, SPEC, NP )
!     Parameter list
!     ----------------------------------------------------------------
!       IP             I   IP of the current node
!       SPEC    R.A.   I   2-D spectrum E(f,theta).
!       NP      Int.   O   Number of partitions.
!                           -1 : Spectrum without minumum energy.
!                            0 : Spectrum with minumum energy.
!                                but no partitions.
!     ----------------------------------------------------------------
!
      USE DATAPOOL, ONLY: NUMSIG, NUMDIR, NSPEC, rkind
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)           :: IP
      INTEGER, INTENT(OUT)          :: NP
      REAL(rkind), INTENT(IN)              :: SPEC(NUMSIG,NUMDIR)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: ITH, IMI(NSPEC), IMD(NSPEC),         &
                                 IMO(NSPEC), IND(NSPEC), NP_MAX,      &
                                 IT(1),        &
                                 IPW, IPT, ISP, NK, NTH
!/S      INTEGER, SAVE           :: IENT = 0
      REAL(rkind)                :: ZP(NSPEC), ZMIN, ZMAX, Z(NSPEC),     &
                                 FACT, WSMAX, HSMAX
      REAL(rkind)                :: TP(DIMP,NPOUTMAX)
!/
!/ ------------------------------------------------------------------- /
! 0.  Initializations
!
!/S      CALL STRACE (IENT, 'W3PART')
!
      NP     = 0
      SPART_PARAMS     = 0.

      NK = NUMSIG
      NTH = NUMDIR
!
! -------------------------------------------------------------------- /
! 1.  Process input spectrum
! 1.a 2-D to 1-D spectrum
!
      DO ITH=1, NTH
        ZP(1+(ITH-1)*NK:ITH*NK) = SPEC(:,ITH)
        END DO

!
! 1.b Invert spectrum and 'digitize'
!
      ZMIN   = MINVAL ( ZP ) 
      ZMAX   = MAXVAL ( ZP ) 
      IF ( ZMAX-ZMIN .LT. 1.E-9 ) RETURN
!
      Z      = ZMAX - ZP
!
      FACT   = REAL(IHMAX-1) / ( ZMAX - ZMIN )
      IMI    = MAX ( 1 , MIN ( IHMAX , NINT ( 1. + Z*FACT ) ) )
!
! 1.c Sort digitized image
!
      CALL PTSORT ( IMI, IND, IHMAX )
!
! -------------------------------------------------------------------- /
! 2.  Perform partitioning
! 2.a Update nearest neighbor info as needed.
!
      CALL PTNGHB
!
! 2.b Incremental flooding
!
      CALL PT_FLD ( IMI, IND, IMO, ZP, NP_MAX )
!
! 2.c Compute parameters per partition
!     NP and NX initialized inside routine.
!
      CALL PTMEAN ( NP_MAX, IMO, IP, SPEC, NP )
!
      RETURN

      END SUBROUTINE WWM_PART



      SUBROUTINE PTSORT ( IMI, IND, IHMAX )
!  1. Purpose :
!
!     This subroutine sorts the image data in ascending order.
!     This sort original to F.T.Tracy (2006)
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMI     I.A.   I   Input discretized spectrum.
!       IND     I.A.   O   Sorted data.
!       IHMAX   Int.   I   Number of integer levels.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Sur.  W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
!/S      USE W3SERVMD, ONLY: STRACE
!
      USE DATAPOOL, ONLY: NSPEC
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)      :: IHMAX, IMI(NSPEC)
      INTEGER, INTENT(OUT)     :: IND(NSPEC)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: I, IN, IV
!/S      INTEGER, SAVE           :: IENT = 0
      INTEGER                 :: NUMV(IHMAX), IADDR(IHMAX),           &
                                 IORDER(NSPEC)
!/
!/S      CALL STRACE (IENT, 'PTSORT')
!
! -------------------------------------------------------------------- /
! 1.  Occurences per height
!
      NUMV   = 0
      DO I=1, NSPEC
        NUMV(IMI(I)) = NUMV(IMI(I)) + 1
        END DO
!
! -------------------------------------------------------------------- /
! 2.  Starting address per height
!
      IADDR(1) = 1
      DO I=1, IHMAX-1
        IADDR(I+1) = IADDR(I) + NUMV(I)
      END DO
!
! -------------------------------------------------------------------- /
! 3.  Order points
!
      DO I=1, NSPEC
        IV        = IMI(I)
        IN        = IADDR(IV)
        IORDER(I) = IN
        IADDR(IV) = IN + 1
        END DO
!
! -------------------------------------------------------------------- /
! 4.  Sort points
!
      DO I=1, NSPEC
        IND(IORDER(I)) = I
        END DO
!
      RETURN
!/
!/ End of PTSORT ----------------------------------------------------- /
!/
      END SUBROUTINE PTSORT


      SUBROUTINE PTNGHB 
!/
!  1. Purpose :
!
!     This subroutine computes the nearest neighbors for each grid
!     point. Wrapping of directional distribution (0 to 360)is taken
!     care of using the nearest neighbor system
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMI     I.A.   I   Input discretized spectrum.
!       IMD     I.A.   O   Sorted data.
!       IHMAX   Int.   I   Number of integer levels.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Sur.  W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
!/S      USE W3SERVMD, ONLY: STRACE
!
      USE DATAPOOL, ONLY: NUMSIG, NUMDIR, NSPEC
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!     INTEGER, INTENT(IN)      :: IHMAX, IMI(NSPEC)
!     INTEGER, INTENT(IN)      :: IMD(NSPEC)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: N, J, I, K, NK, NTH
!/S      INTEGER, SAVE           :: IENT = 0
!/
!/S      CALL STRACE (IENT, 'PTNGHB')
!
! -------------------------------------------------------------------- /
! 1.  Check on need of processing
!
      NK = NUMSIG
      NTH = NUMDIR
      IF ( MK.EQ.NK .AND. MTH.EQ.NTH ) RETURN
!
      IF ( MK.GT.0 ) DEALLOCATE ( NEIGH )
      ALLOCATE ( NEIGH(9,NSPEC) )
      MK     = NK
      MTH    = NTH
!
! -------------------------------------------------------------------- /
! 2.  Build map
!
      NEIGH  = 0
!
! ... Base loop
!
      DO N = 1, NSPEC
!
        J      = (N-1) / NK + 1
        I      = N - (J-1) * NK
        K      = 0
!
! ... Point at the left(1)
!
        IF ( I .NE. 1 ) THEN
            K           = K + 1
            NEIGH(K, N) = N - 1
          END IF
!
! ... Point at the right (2)
!
        IF ( I .NE. NK ) THEN 
            K           = K + 1
            NEIGH(K, N) = N + 1
          END IF
!
! ... Point at the bottom(3)
!
        IF ( J .NE. 1 ) THEN
            K           = K + 1
            NEIGH(K, N) = N - NK
          END IF
!
! ... ADD Point at bottom_wrap to top
!
        IF ( J .EQ. 1 ) THEN
            K          = K + 1
            NEIGH(K,N) = NSPEC - (NK-I)
          END IF
!
! ... Point at the top(4)
!
        IF ( J .NE. NTH ) THEN
            K           = K + 1
            NEIGH(K, N) = N + NK
          END IF
!
! ... ADD Point to top_wrap to bottom
!
         IF ( J .EQ. NTH ) THEN
             K          = K + 1
             NEIGH(K,N) = N - (NTH-1) * NK
            END IF
!
! ... Point at the bottom, left(5)
!
        IF ( (I.NE.1) .AND. (J.NE.1) ) THEN
            K           = K + 1
            NEIGH(K, N) = N - NK - 1
          END IF
!
! ... Point at the bottom, left with wrap.
!
         IF ( (I.NE.1) .AND. (J.EQ.1) ) THEN
             K          = K + 1
             NEIGH(K,N) = N - 1 + NK * (NTH-1)
           END IF
!
! ... Point at the bottom, right(6)
!
        IF ( (I.NE.NK) .AND. (J.NE.1) ) THEN
            K           = K + 1
            NEIGH(K, N) = N - NK + 1
          END IF
!
! ... Point at the bottom, right with wrap
!
        IF ( (I.NE.NK) .AND. (J.EQ.1) ) THEN
            K           = K + 1
            NEIGH(K,N) = N + 1 + NK * (NTH - 1)
          END  IF
!
! ... Point at the top, left(7)
!
        IF ( (I.NE.1) .AND. (J.NE.NTH) ) THEN
            K           = K + 1
            NEIGH(K, N) = N + NK - 1
          END IF
!
! ... Point at the top, left with wrap
!
         IF ( (I.NE.1) .AND. (J.EQ.NTH) ) THEN
             K           = K + 1
             NEIGH(K,N) = N - 1 - (NK) * (NTH-1)
           END IF
!
! ... Point at the top, right(8)
!
        IF ( (I.NE.NK) .AND. (J.NE.NTH) ) THEN
            K           = K + 1
            NEIGH(K, N) = N + NK + 1
          END IF
!
! ... Point at top, right with wrap
!
!
        IF ( (I.NE.NK) .AND. (J.EQ.NTH) ) THEN
            K           = K + 1
            NEIGH(K,N) = N + 1 - (NK) * (NTH-1)
          END IF
!
        NEIGH(9,N) = K
!
        END DO
!
      RETURN
!/
!/ End of PTNGHB ----------------------------------------------------- /
!/
      END SUBROUTINE PTNGHB



      SUBROUTINE PT_FLD ( IMI, IND, IMO, ZP, NPART )
!  1. Purpose :
!
!     This subroutine does incremental flooding of the image to
!     determine the watershed image.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMI     I.A.   I   Input discretized spectrum.
!       IND     I.A.   I   Sorted addresses.
!       IMO     I.A.   O   Output partitioned spectrum.
!       ZP      R.A.   I   Spectral array.
!       NPART   Int.   O   Number of partitions found.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Sur.  W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
!/S      USE W3SERVMD, ONLY: STRACE
!
      USE DATAPOOL, ONLY: NSPEC
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: IMI(NSPEC), IND(NSPEC)
      INTEGER, INTENT(OUT)    :: IMO(NSPEC), NPART
      REAL(rkind), INTENT(IN) :: ZP(NSPEC)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: MASK, INIT, IWSHED, IMD(NSPEC),      &
                                 IC_LABEL, IFICT_PIXEL, M, IH, MSAVE, &
                                 IP, I, IPP, IC_DIST, IEMPTY, IPPP,   &
                                 JL, JN, IPT, J
      INTEGER                 :: IQ(NSPEC), IQ_START, IQ_END
!/S      INTEGER, SAVE           :: IENT = 0
      REAL(rkind)             :: ZPMAX, EP1, DIFF
!/
!/S      CALL STRACE (IENT, 'PT_FLD')
!
! -------------------------------------------------------------------- /
! 0.  Initializations
!
      MASK        = -2
      INIT        = -1
      IWSHED      =  0
      IMO         = INIT
      IC_LABEL    =  0
      IMD         =  0
      IFICT_PIXEL = -100
!
      IQ_START    =  1
      IQ_END      =  1
!
      ZPMAX       = MAXVAL ( ZP )
!
! -------------------------------------------------------------------- /
! 1.  Loop over levels
!
      M      =  1
!
      DO IH=1, IHMAX
        MSAVE  = M
!
! 1.a Pixels at level IH
!
        DO
          IP     = IND(M)
          IF ( IMI(IP) .NE. IH ) EXIT
!
!     Flag the point, if it stays flagge, it is a separate minimum.
!
          IMO(IP) = MASK
!
!     Consider neighbors. If there is neighbor, set distance and add
!     to queue.
!
          DO I=1, NEIGH(9,IP)
            IPP    = NEIGH(I,IP)
            IF ( (IMO(IPP).GT.0) .OR. (IMO(IPP).EQ.IWSHED) ) THEN
                IMD(IP) = 1
                CALL FIFO_ADD (IP)
                EXIT
              END IF
            END DO
!
          IF ( M+1 .GT. NSPEC ) THEN
              EXIT
            ELSE
              M = M + 1
            END IF
!
          END DO
!
! 1.b Process the queue
!
        IC_DIST = 1
        CALL FIFO_ADD (IFICT_PIXEL)
!
        DO
          CALL FIFO_FIRST (IP)
!
!     Check for end of processing
!
          IF ( IP .EQ. IFICT_PIXEL ) THEN
              CALL FIFO_EMPTY (IEMPTY)
              IF ( IEMPTY .EQ. 1 ) THEN
                  EXIT
                ELSE
                  CALL FIFO_ADD (IFICT_PIXEL)
                  IC_DIST = IC_DIST + 1
                  CALL FIFO_FIRST (IP)
                END IF
            END IF
!
!     Process queue
!
          DO I=1, NEIGH(9,IP)
            IPP = NEIGH(I,IP)
!
!     Check for labeled watersheds or basins
!
            IF ( (IMD(IPP).LT.IC_DIST) .AND. ( (IMO(IPP).GT.0) .OR.  &
                 (IMO(IPP).EQ.IWSHED))) THEN
!
                IF ( IMO(IPP) .GT. 0 ) THEN
!
                    IF ((IMO(IP) .EQ. MASK) .OR. (IMO(IP) .EQ. &
                        IWSHED)) THEN
                        IMO(IP) = IMO(IPP)
                      ELSE IF (IMO(IP) .NE. IMO(IPP)) THEN
                        IMO(IP) = IWSHED
                      END IF
!
                  ELSE IF (IMO(IP) .EQ. MASK) THEN
!
                    IMO(IP) = IWSHED
!
                  END IF
!
              ELSE IF ( (IMO(IPP).EQ.MASK) .AND. (IMD(IPP).EQ.0) ) THEN
!
                 IMD(IPP) = IC_DIST + 1
                 CALL FIFO_ADD (IPP)
!
              END IF
!
            END DO
!
          END DO
!
! 1.c Check for mask values in IMO to identify new basins
!
        M = MSAVE
!
        DO
          IP     = IND(M)
          IF ( IMI(IP) .NE. IH ) EXIT
          IMD(IP) = 0
!
          IF (IMO(IP) .EQ. MASK) THEN
!
! ... New label for pixel
!
              IC_LABEL = IC_LABEL + 1
              CALL FIFO_ADD (IP)
              IMO(IP) = IC_LABEL
!
! ... and all connected to it ...
!
              DO
                CALL FIFO_EMPTY (IEMPTY)
                IF ( IEMPTY .EQ. 1 ) EXIT
                CALL FIFO_FIRST (IPP)
!
                DO I=1, NEIGH(9,IPP)
                  IPPP   = NEIGH(I,IPP)
                  IF ( IMO(IPPP) .EQ. MASK ) THEN
                      CALL FIFO_ADD (IPPP)
                      IMO(IPPP) = IC_LABEL
                    END IF
                  END DO
!
                END DO
!
            END IF
!
          IF ( M + 1 .GT. NSPEC ) THEN
              EXIT
            ELSE
              M = M + 1
            END IF
!
          END DO
!
        END DO
!
! -------------------------------------------------------------------- /
! 2.  Find nearest neighbor of 0 watershed points and replace
!     use original input to check which group to affiliate with 0
!     Soring changes first in IMD to assure symetry in adjustment.
!
      DO J=1, 5
        IMD    = IMO
        DO JL=1 , NSPEC
          IPT    = -1
          IF ( IMO(JL) .EQ. 0 ) THEN
              EP1    = ZPMAX
              DO JN=1, NEIGH (9,JL)
                DIFF   = ABS ( ZP(JL) - ZP(NEIGH(JN,JL)))
                IF ( (DIFF.LE.EP1) .AND. (IMO(NEIGH(JN,JL)).NE.0) ) THEN
                    EP1    = DIFF
                    IPT    = JN
                  END IF
                END DO
              IF ( IPT .GT. 0 ) IMD(JL) = IMO(NEIGH(IPT,JL))
            END IF
          END DO
        IMO    = IMD
        IF ( MINVAL(IMO) .GT. 0 ) EXIT
        END DO
!
      NPART = IC_LABEL
!
      RETURN
!
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE FIFO_ADD ( IV )
!
!     Add point to FIFO queue.
!
      INTEGER, INTENT(IN)      :: IV
!
      IQ(IQ_END) = IV
!
      IQ_END = IQ_END + 1
      IF ( IQ_END .GT. NSPEC ) IQ_END = 1
!
      RETURN
      END SUBROUTINE
!/ ------------------------------------------------------------------- /
      SUBROUTINE FIFO_EMPTY ( IEMPTY )
!
!     Check if queue is empty.
!
      INTEGER, INTENT(OUT)     :: IEMPTY
!
      IF ( IQ_START .NE. IQ_END ) THEN
        IEMPTY = 0
      ELSE
        IEMPTY = 1
      END IF
!
      RETURN
      END SUBROUTINE
!/ ------------------------------------------------------------------- /
      SUBROUTINE FIFO_FIRST ( IV )
!
!     Get point out of queue.
!
      INTEGER, INTENT(OUT)     :: IV
!
      IV = IQ(IQ_START)
!
      IQ_START = IQ_START + 1
      IF ( IQ_START .GT. NSPEC ) IQ_START = 1
!
      RETURN
      END SUBROUTINE
!/
!/ End of PT_FLD ----------------------------------------------------- /
!/
      END SUBROUTINE PT_FLD



      SUBROUTINE PTMEAN ( NPI, IMO, IP, SPEC, NPO )
!  1. Purpose :
!
!     Compute mean parameters per partition.
!     The output is stored in the array SPART_PARAMS
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       NPI     Int.   I   Number of partitions found.
!       IMO     I.A.   I   Partition map.
!       IP      Int.   I   Id of the current node
!       SPEC      R.A.   I   Input spectrum.
!       NPO     Int.   O   Number of partitions with mean parameters.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Sur.  W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
!
      USE DATAPOOL, ONLY: NUMSIG, NUMDIR, NSPEC, rkind
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: NPI, IMO(NSPEC), IP
      INTEGER, INTENT(OUT)    :: NPO
      REAL(rkind), INTENT(IN)        :: SPEC(NUMSIG,NUMDIR)

      INTEGER  :: IC_LABEL, MSK(NUMSIG,NUMDIR)
      REAL(rkind)     :: IMO2D(NUMSIG,NUMDIR), PSPEC(NUMSIG,NUMDIR)
      INTEGER                 :: IK, ITH, ISP, IPART, IFPMAX(0:NPI), INDX(1)
      REAL(rkind)                    :: HS,TM01,TM02,KLM,WLM,TM10
      REAL(rkind)                    :: ETOTS,ETOTC,DM,DSPR
      REAL(rkind)                    :: XPALL(DIMP,NPI)


      NPO = MIN(NPI, NPOUTMAX)
      SPART_PARAMS = 0

      IF ( NPO .EQ. 0 ) RETURN
!
      DO ITH=1, NUMDIR
        IMO2D(:,ITH) = IMO(1+(ITH-1)*NUMSIG:ITH*NUMSIG)
      END DO

! computing the parameters for all the peaks
      XPALL = 0
      DO IC_LABEL=1, NPI
        MSK = MERGE(1,0,IMO2D.EQ.IC_LABEL)
        PSPEC = SPEC*MSK
        CALL MEAN_PARAMETER(IP,PSPEC,NUMSIG,HS,TM01,TM02,TM10,KLM,WLM)
        CALL MEAN_DIRECTION_AND_SPREAD(IP,PSPEC,NUMSIG,ETOTS,ETOTC,DM,DSPR)
        XPALL(1,IC_LABEL) = HS
        XPALL(2,IC_LABEL) = TM01
        XPALL(3,IC_LABEL) = TM02
        XPALL(4,IC_LABEL) = WLM
        XPALL(5,IC_LABEL) = DM
        XPALL(6,IC_LABEL) = DSPR
      ENDDO
     ! substituting nan with 0
      WHERE(XPALL /= XPALL)
        XPALL = 0
      END WHERE

      SPART_PARAMS = 0
! sorting by HS and getting the first NPO peaks
      DO IPART=1, NPO
        INDX          = MAXLOC(XPALL(1,1:NPI))
        SPART_PARAMS(:,IPART)    = XPALL(:,INDX(1))
        XPALL(:,INDX) = -1
      END DO

      END SUBROUTINE PTMEAN

      END MODULE WWM_PARTMD
