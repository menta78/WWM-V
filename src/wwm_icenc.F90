#include "wwm_functions.h"
! author: Lorenzo Mentaschi
! for now, the ice concentration is loaded from a netcdf on a rectangular grid (regular or rotated or curvilinear)
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_ICE_INPUT
      USE NETCDF
      USE DATAPOOL

      IMPLICIT NONE
      INTEGER           :: fid, varid, ndims, dimids(nf90_max_var_dims)
      REAL(rkind), ALLOCATABLE :: CF_LON(:,:), CF_LAT(:,:), CF_LON_1D(:), CF_LAT_1D(:)
      character (len = *), parameter :: CallFct="INIT_NETCDF_CF"
      character (len=100) :: Xname, Yname
      type(FD_FORCING_GRID) TheInfo
      integer IX, IY
      REAL(rkind) :: cf_w1, cf_w2

      IF (PrintLOG) THEN
        WRITE(STAT%FHNDL,*) 'INIT_ICE_INPUT start'
        FLUSH(STAT%FHNDL)
      ENDIF

      MNP_ICE=np_total
      IF (myrank .eq. 0) THEN
        allocate(XP_ICE(MNP_ICE), YP_ICE(MNP_ICE), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_icenc, allocate error 1')
        XP_ICE=XPtotal
        YP_ICE=YPtotal

      ENDIF

      allocate(ICECONC(MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_icenc, allocate error 1')
      ICECONC=0
      CALL INIT_DIRECT_NETCDF_CF_ICE(eVAR_ICE, .FALSE., ICEIN%FNAME, ICENCVAR)

      IF (myrank .eq. 0) THEN
        ISTAT = nf90_open(ICEIN%FNAME, nf90_nowrite, fid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, ISTAT)

        ISTAT = nf90_inq_varid(fid, ICENCVAR, varid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, ISTAT)

        ISTAT = nf90_inq_varid(fid, ICEXNCVAR, varid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, ISTAT)

        ISTAT = nf90_inq_varid(fid, ICEYNCVAR, varid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, ISTAT)

        Xname=ICEXNCVAR
        Yname=ICEYNCVAR
        IF (PrintLOG) THEN
          WRITE(STAT%FHNDL,*) 'Xname=', TRIM(Xname)
          WRITE(STAT%FHNDL,*) 'Yname=', TRIM(Yname)
          FLUSH(STAT%FHNDL)
        END IF

        ! Reading lontitude/latitude array

        ISTAT = nf90_inq_varid(fid, Xname, varid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, ISTAT)

        ISTAT = nf90_inquire_variable(fid, varid, ndims=ndims, dimids=dimids)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, ISTAT)

        IF (ndims.EQ.2) THEN
          ISTAT = nf90_inquire_dimension(fid, dimids(1), len=NDX_ICE_FD)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, ISTAT)

          ISTAT = nf90_inquire_dimension(fid, dimids(2), len=NDY_ICE_FD)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, ISTAT)

          IF (PrintLOG) THEN
            WRITE(STAT%FHNDL,*) 'NDX_ICE_FD=', NDX_ICE_FD
            WRITE(STAT%FHNDL,*) 'NYX_ICE_FD=', NDY_ICE_FD
            FLUSH(STAT%FHNDL)
          END IF

          allocate(CF_LON(NDX_ICE_FD, NDY_ICE_FD), CF_LAT(NDX_ICE_FD, NDY_ICE_FD), ICECONC_FD(NDX_ICE_FD, NDY_ICE_FD), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_icenc, allocate error 47')

          ISTAT = nf90_inq_varid(fid, Xname, varid)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, ISTAT)

          ISTAT = nf90_get_var(fid, varid, CF_LON)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, ISTAT)

          ISTAT = nf90_inq_varid(fid, Yname, varid)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, ISTAT)

          ISTAT = nf90_get_var(fid, varid, CF_LAT)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, ISTAT)
        ELSE IF (ndims.EQ.1) THEN
          ISTAT = nf90_inquire_dimension(fid, dimids(1), len=NDX_ICE_FD)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, ISTAT)

          ALLOCATE(CF_LON_1D(NDX_ICE_FD))
          IF (istat/=0) CALL WWM_ABORT('wwm_icenc, allocate error 47')
          ISTAT = nf90_get_var(fid, varid, CF_LON_1D)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, ISTAT)

          ISTAT = nf90_inq_varid(fid, Yname, varid)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 15, ISTAT)

          ISTAT = nf90_inquire_variable(fid, varid, ndims=ndims, dimids=dimids)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 16, ISTAT)
          IF (ndims.NE.1) CALL WWM_ABORT('wwm_icenc, x and y must be both 2d or both 1d')

          ISTAT = nf90_inquire_dimension(fid, dimids(1), len=NDY_ICE_FD)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 17, ISTAT)

          ALLOCATE(CF_LAT_1D(NDY_ICE_FD))
          IF (istat/=0) CALL WWM_ABORT('wwm_icenc, allocate error 47')
          ISTAT = nf90_get_var(fid, varid, CF_LAT_1D)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 18, ISTAT)

          allocate(CF_LON(NDX_ICE_FD, NDY_ICE_FD), CF_LAT(NDX_ICE_FD, NDY_ICE_FD), ICECONC_FD(NDX_ICE_FD, NDY_ICE_FD), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_icenc, allocate error 47')

          DO IX=1,NDX_ICE_FD
            DO IY=1,NDY_ICE_FD
              CF_LON(IX,IY) = CF_LON_1D(IX)
              CF_LAT(IX,IY) = CF_LAT_1D(IY)
            END DO
          END DO
        ELSE
          ! x and y must be 1-d or 2-d
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 555, ISTAT)
        END IF
        TheInfo % nx_dim = NDX_ICE_FD
        TheInfo % ny_dim = NDY_ICE_FD
        allocate(TheInfo%LON(NDX_ICE_FD, NDY_ICE_FD), TheInfo%LAT(NDX_ICE_FD, NDY_ICE_FD), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_icenc, allocate error 47')
        DO IX=1,NDX_ICE_FD
          DO IY=1,NDY_ICE_FD
            TheInfo % LON(IX,IY) = CF_LON(IX,IY)
            TheInfo % LAT(IX,IY) = CF_LAT(IX,IY)
          END DO
        END DO
        DEALLOCATE(CF_LON, CF_LAT)
        CALL COMPUTE_CF_COEFFICIENTS_ICE(TheInfo)
        Deallocate(TheInfo % LON, TheInfo % LAT)
      END IF

      ALLOCATE(tmp_ice1(MNP),tmp_ice2(MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_icenc, allocate error 1')
      CALL GET_CF_TIME_INDEX_ICE(eVAR_ICE, REC1_ice_new,REC2_ice_new,cf_w1,cf_w2)
      CALL READ_INTERP_NETCDF_CF_WWM_ICE(REC1_ice_new,tmp_ice1)
      IF (cf_w1.NE.1) THEN
        CALL READ_INTERP_NETCDF_CF_WWM_ICE(REC2_ice_new,tmp_ice2)
        ICECONC(:) = cf_w1*tmp_ice1(:)+cf_w2*tmp_ice2(:)
      ELSE
        ICECONC(:) = cf_w1*tmp_ice1(:)
      END IF

      END SUBROUTINE INIT_ICE_INPUT
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE UPDATE_ICE(K)
      USE NETCDF
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind)             :: TMP(MNP,2)
      REAL(rkind)             :: cf_w1, cf_w2
      INTEGER                 :: IT, IFILE
      INTEGER, intent(in)     :: K
!AR: All crap ... defining K without using means that nobody has ever checked the results or anything else, so why coding at all?
!AR: Mathieu can you please fix this !!!
!LM: keeping this comment by AR, taken from UPDATE_WIND (hopefully the problem is now solved?)
      IF (PrintLOG) THEN
        WRITE(STAT%FHNDL,*) 'MAIN%TMJD=', MAIN%TMJD
        WRITE(STAT%FHNDL,*) 'SEWI(TMJD,EMJD)=', SEWI%TMJD, SEWI%EMJD
      END IF
      IF ( (MAIN%TMJD .ge. SEWI%TMJD-1.E-8) .AND. (MAIN%TMJD .le. SEWI%EMJD+1.e-8) ) THEN
        IF (K.EQ.1) THEN
          REC1_ice_old = 0
          REC2_ice_old = 0
        END IF
        CALL GET_CF_TIME_INDEX_ICE(eVAR_ICE, REC1_ice_new,REC2_ice_new,cf_w1,cf_w2)
        IF (REC1_ice_new.NE.REC1_ice_old) THEN
          CALL READ_INTERP_NETCDF_CF_WWM_ICE(REC1_ice_new,tmp_ice1)
        END IF
        IF (REC2_ice_new.NE.REC2_ice_old) THEN
          CALL READ_INTERP_NETCDF_CF_WWM_ICE(REC2_ice_new,tmp_ice2)
        END IF
        IF (cf_w1.NE.1) THEN
          ICECONC(:) = cf_w1*tmp_ice1(:)+cf_w2*tmp_ice2(:)
        ELSE
          ICECONC(:) = cf_w1*tmp_ice1(:)
        END IF
        REC1_ice_old = REC1_ice_new
        REC2_ice_old = REC2_ice_new
      END IF
      IF (PrintLOG) THEN
        WRITE(STAT%FHNDL,*) 'ICECONC, min/max=', minval(ICECONC), maxval(ICECONC)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_CF_COEFFICIENTS_ICE(TheInfo)
      USE DATAPOOL
      IMPLICIT NONE
      type(FD_FORCING_GRID), intent(in) :: TheInfo
      integer I
      integer eCF_IX, eCF_IY
      real(rkind) eCF_COEFF(4)
      integer :: nbExtrapolation = 0
      character(len=256) :: FileSave = "wwm_filesave_interp_array_ice.nc"
      logical success
      logical EXTRAPO_OUT
      real(rkind) eX, eY
      IF (PrintLOG) THEN
        WRITE(STAT%FHNDL,*) 'Starting node loop for calcs of coefs (ice)'
      END IF
      allocate(CF_IX_ICE(MNP_ICE), CF_IY_ICE(MNP_ICE), CF_COEFF_ICE(4,MNP_ICE), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_icenc, allocate error 52')
      IF (PrintLOG) THEN
        WRITE(STAT%FHNDL,*) 'LSAVE_INTERP_ARRAY_ICE=', LSAVE_INTERP_ARRAY_ICE
      END IF

      IF (LSAVE_INTERP_ARRAY_ICE) THEN
        CALL LOAD_INTERP_ARRAY_ICE(FileSave, success)
        IF (PrintLOG) THEN
          WRITE(STAT%FHNDL,*) 'success=', success
        END IF
        IF (success .eqv. .TRUE.) RETURN
      END IF

      CF_IX_ICE=0
      CF_IY_ICE=0
      CF_COEFF_ICE=0
      IF (PrintLOG) THEN
        WRITE(STAT%FHNDL,*) 'min(lon)=', minval(TheInfo % LON)
        WRITE(STAT%FHNDL,*) 'max(lon)=', maxval(TheInfo % LON)
        WRITE(STAT%FHNDL,*) 'min(lat)=', minval(TheInfo % LAT)
        WRITE(STAT%FHNDL,*) 'max(lat)=', maxval(TheInfo % LAT)
      END IF
      DO I = 1, MNP_ICE
        eX=XP_ICE(I)
        eY=YP_ICE(I)
        CALL COMPUTE_SINGLE_INTERPOLATION_INFO(TheInfo, EXTRAPOLATION_ALLOWED_ICE, eX, eY, eCF_IX, eCF_IY, eCF_COEFF, EXTRAPO_OUT)
        IF (PrintLOG) THEN
           WRITE(STAT%FHNDL,'(4I10,10F20.10)') I, MNP_ICE, eCF_IX, eCF_IY, eCF_COEFF
        END IF
        CF_IX_ICE(I) = eCF_IX
        CF_IY_ICE(I) = eCF_IY
        CF_COEFF_ICE(:,I) = eCF_COEFF
        IF (EXTRAPO_OUT) nbExtrapolation = nbExtrapolation + 1
      END DO

      IF (LSAVE_INTERP_ARRAY_ICE) THEN
        CALL SAVE_INTERP_ARRAY_ICE(FileSave)
      END IF

      IF (PrintLOG) THEN
        IF (EXTRAPOLATION_ALLOWED_ICE) THEN
          WRITE(STAT%FHNDL,*) ' ICE interp. nbExtrapolation=', nbExtrapolation
        END IF
        WRITE(STAT%FHNDL,*) ' ICE: done interp calcs'
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE LOAD_INTERP_ARRAY_ICE(FileSave, success)
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
      logical, intent(out) :: success
      character(len=256), intent(in) :: FileSave
      character (len = *), parameter :: CallFct = "LOAD_INTERP_ARRAY_ICE"
      integer, allocatable :: CF_IX_GLOBAL(:), CF_IY_GLOBAL(:)
      real(rkind), allocatable :: CF_COEFF_GLOBAL(:,:)
      integer iret, ncid, varid
      integer IP, IPglob, iPROC, NPloc, IPloc
      INQUIRE(FILE=TRIM(FileSave), EXIST=LPRECOMP_EXIST)
      IF (LPRECOMP_EXIST .eqv. .FALSE.) THEN
        success=.FALSE.
        RETURN
      END IF
      success=.TRUE.
      IF (myrank .eq. 0) THEN
        allocate(CF_IX_GLOBAL(np_total), CF_IY_GLOBAL(np_total), CF_COEFF_GLOBAL(4,np_total), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_icenc, allocate error 52')
        !
        iret=nf90_open(TRIM(FileSave), NF90_NOWRITE, ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)
        !
        iret=nf90_inq_varid(ncid, "CF_IX", varid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, ISTAT)
        iret=NF90_GET_VAR(ncid, varid, CF_IX_GLOBAL)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, ISTAT)
        !
        iret=nf90_inq_varid(ncid, "CF_IY", varid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, ISTAT)
        iret=NF90_GET_VAR(ncid, varid, CF_IY_GLOBAL)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, ISTAT)
        !
        iret=nf90_inq_varid(ncid, "CF_COEFF", varid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, ISTAT)
        iret=NF90_GET_VAR(ncid, varid, CF_COEFF_GLOBAL)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, ISTAT)
        !
        iret=nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 27, iret)
        !
        CF_IX_ICE=CF_IX_GLOBAL
        CF_IY_ICE=CF_IY_GLOBAL
        CF_COEFF_ICE=CF_COEFF_GLOBAL
        deallocate(CF_IX_GLOBAL, CF_IY_GLOBAL, CF_COEFF_GLOBAL)
      ENDIF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SAVE_INTERP_ARRAY_ICE(FileSave)
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
      character(len=256), intent(in) :: FileSave
      character (len = *), parameter :: CallFct = "SAVE_INTERP_ARRAY"
      integer, allocatable :: CF_IX_GLOBAL(:), CF_IY_GLOBAL(:)
      real(rkind), allocatable :: CF_COEFF_GLOBAL(:,:)
      integer, allocatable :: CF_IX_loc(:), CF_IY_loc(:)
      real(rkind), allocatable :: CF_COEFF_loc(:,:)
      integer, allocatable :: ListFirstMNP(:)
      integer ncid, iret, var_id
      integer mnp_dims, four_dims
      integer IP, IPglob, iPROC, NP_RESloc, IPloc
      IF (PrintLOG) THEN
        WRITE(STAT%FHNDL,*) 'minval(CF_IX_ICE)=', minval(CF_IX_ICE)
        WRITE(STAT%FHNDL,*) 'minval(CF_IY_ICE)=', minval(CF_IY_ICE)
        WRITE(STAT%FHNDL,*) 'minval(CF_COEFF_ICE)=', minval(CF_COEFF_ICE)
      END IF

      allocate(CF_IX_GLOBAL(np_total), CF_IY_GLOBAL(np_total), CF_COEFF_GLOBAL(4,np_total), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_icenc, allocate error 52')
      CF_IX_GLOBAL=CF_IX_ICE
      CF_IY_GLOBAL=CF_IY_ICE
      CF_COEFF_GLOBAL=CF_COEFF_ICE
      !
      ! Now writing up
      !
      IF (myrank .eq. 0) THEN
        iret=nf90_create(TRIM(FileSave), nf90_CLOBBER, ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)
        !
        iret = nf90_def_dim(ncid, 'mnp', np_total, mnp_dims)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, iret)
        !
        iret = nf90_def_dim(ncid, 'four', 4, four_dims)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, iret)
        !
        iret=nf90_def_var(ncid,'CF_IX',NF90_INT,(/ mnp_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, iret)
        !
        iret=nf90_def_var(ncid,'CF_IY',NF90_INT,(/ mnp_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)
        !
        iret=nf90_def_var(ncid,'CF_COEFF',NF90_RUNTYPE,(/ four_dims, mnp_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
        !
        iret=nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, iret)
        !
        ! Now writing the data
        !
        IF (PrintLOG) THEN
          WRITE(STAT%FHNDL,*) 'ICE: minval(CF_IX_GLOBAL)=', minval(CF_IX_GLOBAL)
          WRITE(STAT%FHNDL,*) 'ICE: minval(CF_IY_GLOBAL)=', minval(CF_IY_GLOBAL)
          WRITE(STAT%FHNDL,*) 'ICE: minval(CF_COEFF_GLOBAL)=', minval(CF_COEFF_GLOBAL)
        END IF
        !
        iret=nf90_open(TRIM(FileSave), NF90_WRITE, ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, iret)
        !
        iret=nf90_inq_varid(ncid,'CF_IX',var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, iret)
        !
        iret=nf90_put_var(ncid,var_id,CF_IX_GLOBAL,start = (/ 1 /), count=(/np_total/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, iret)
        !
        iret=nf90_inq_varid(ncid,'CF_IY',var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, iret)
        !
        iret=nf90_put_var(ncid,var_id,CF_IY_GLOBAL,start = (/ 1 /), count=(/np_total/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, iret)
        !
        iret=nf90_inq_varid(ncid,'CF_COEFF',var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, iret)
        !
        iret=nf90_put_var(ncid,var_id,CF_COEFF_GLOBAL,start = (/ 1, 1 /), count=(/4, np_total/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, iret)
        !
        iret=nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 15, iret)
        !
        deallocate(CF_IX_GLOBAL, CF_IY_GLOBAL, CF_COEFF_GLOBAL)
      END IF
      END SUBROUTINE SAVE_INTERP_ARRAY_ICE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_INTERP_NETCDF_CF_WWM_ICE(RECORD_IN, outice)
      USE NETCDF
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(in)                :: RECORD_IN
      REAL(rkind), INTENT(out)           :: outice(MNP)
      REAL(rkind) :: varTotal(MNP_ICE), Vlocal(MNP)
      character (len = *), parameter :: CallFct="READ_INTERP_NETCDF_CF_WWM_ICE"
      INTEGER                            :: FID, ID
      IF (myrank .eq. 0) THEN
        CALL TEST_FILE_EXIST_DIE("Missing ice file : ", TRIM(ICEIN%FNAME))
        ISTAT = NF90_OPEN(ICEIN%FNAME, NF90_NOWRITE, FID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, ISTAT)

        ISTAT = NF90_inq_varid(FID, ICENCVAR, ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, ISTAT)

        ISTAT = NF90_GET_VAR(FID, ID, ICECONC_FD, start = (/ 1, 1, RECORD_IN /), count = (/ NDX_ICE_FD, NDY_ICE_FD, 1 /))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, ISTAT)

        ISTAT = NF90_CLOSE(FID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, ISTAT)
        CALL KERNEL_INTERP_ICEFD(varTotal)
      END IF
      CALL SCATTER_ONED_ARRAY(varTotal, Vlocal)
      outice(:)=Vlocal
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE KERNEL_INTERP_ICEFD(outice)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER I, J
      REAL(rkind), INTENT(out)           :: outice(MNP_ICE)
      REAL(rkind) :: ic
      INTEGER IX, IY
      REAL(rkind) :: cf_scale_factor, cf_add_offset
      INTEGER SHIFTXY(4,2)
      SHIFTXY(1,1)=0
      SHIFTXY(1,2)=0
      SHIFTXY(2,1)=1
      SHIFTXY(2,2)=0
      SHIFTXY(3,1)=0
      SHIFTXY(3,2)=1
      SHIFTXY(4,1)=1
      SHIFTXY(4,2)=1
      cf_scale_factor = eVAR_ICE % cf_scale_factor
      cf_add_offset = eVAR_ICE % cf_add_offset
      DO I = 1, MNP_ICE
        ic=ZERO
        IX=CF_IX_ICE(I)
        IY=CF_IY_ICE(I)
        DO J=1,4
          ic=ic + CF_COEFF_ICE(J,I)*ICECONC_FD(IX+SHIFTXY(J,1),IY+SHIFTXY(J,2))
        END DO
        outice(I)=ic*cf_scale_factor + cf_add_offset
        IF (outice(I).GT.1) outice(I) = 1
        IF (outice(I).LT.0) outice(I) = 0
      END DO
      IF (PrintLOG) THEN
        WRITE(STAT%FHNDL,*) 'KERNEL_INTERP_ICEFD'
        WRITE(STAT%FHNDL,*) 'ICECONC_FD, min/max=', minval(ICECONC_FD), maxval(ICECONC_FD)
        WRITE(STAT%FHNDL,*) 'ice interp, min/max=', minval(outice), maxval(outice)
        FLUSH(STAT%FHNDL)
      END IF
      END SUBROUTINE
!****************************************************************************
!*                                                                          *
!****************************************************************************
      SUBROUTINE INIT_DIRECT_NETCDF_CF_ICE(eVAR, MULTIPLE_IN, eFileName, eString)
      USE NETCDF
      USE DATAPOOL
      IMPLICIT NONE
      TYPE(VAR_NETCDF_CF), intent(inout) :: eVAR
      logical, intent(in) :: MULTIPLE_IN
      character(len=100), intent(in) :: eFileName
      character(len=*), intent(in) :: eString
!
      INTEGER           :: fid, varid
      INTEGER           :: dimidsB(nf90_max_var_dims), dimids(nf90_max_var_dims)
      integer nbChar
      character (len=20) :: eTimeStr
      character(len=100) :: CHRERR
      character (len = *), parameter :: CallFct="INIT_DIRECT_NETCDF_CF"
      character (len=100) :: eStrUnitTime
      real(rkind) :: ConvertToDay
      real(rkind) :: eTimeStart
      real(rkind) :: cf_scale_factor, cf_add_offset
      integer nbtime_mjd
      real(rkind), allocatable :: ice_time_mjd(:)
      integer eInt(1)
      real(rkind) :: eReal(2)
      integer IPROC
      integer np_file, NumberDim
      eVAR % MULTIPLE_IN = .FALSE.
      IF (myrank .eq. 0) THEN
        ISTAT = nf90_open(TRIM(eFileName), nf90_nowrite, fid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, ISTAT)

        ! Reading ice attributes

        WRITE(STAT%FHNDL,*), 'eString=', TRIM(eString)
        WRITE(STAT%FHNDL,*), 'FNAME=', TRIM(eFileName)
        ISTAT = nf90_inq_varid(fid, TRIM(eString), varid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, ISTAT)

        ISTAT = nf90_inquire_variable(fid, varid, dimids=dimidsB)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, ISTAT)

        ISTAT = nf90_inquire_variable(fid, varid, ndims=NumberDim)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, ISTAT)



        ISTAT = nf90_inquire_dimension(fid, dimidsB(1), len=np_file)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, ISTAT)
        IF (NumberDim .eq. 2) THEN
          IF (np_file .ne. np_total) THEN
            Print *, 'We have np_file=', np_file
            Print *, '   But np_total=', np_total
            CALL WWM_ABORT("Direct forcing file seems to have wrong dimension")
          END IF
        END IF

        ISTAT = nf90_inquire_dimension(fid, dimidsB(NumberDim), name=eTimeStr)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, ISTAT)
        IF (PrintLOG) THEN
          WRITE(STAT%FHNDL,*) 'ice: variable used for time=', TRIM(eTimeStr)
          FLUSH(STAT%FHNDL)
        END IF

        ISTAT = nf90_get_att(fid, varid, "scale_factor", cf_scale_factor)
        IF (ISTAT /= 0) THEN
          CHRERR = nf90_strerror(ISTAT)
          IF (PrintLOG) THEN
            WRITE(STAT%FHNDL,*) 'CHRERR=', TRIM(CHRERR)
          END IF
          cf_scale_factor=ONE
        ENDIF
        IF (PrintLOG) THEN
          WRITE(STAT%FHNDL,*) 'cf_scale_factor=', cf_scale_factor
          FLUSH(STAT%FHNDL)
        END IF

        ISTAT = nf90_get_att(fid, varid, "add_offset", cf_add_offset)
        IF (ISTAT /= 0) THEN
          CHRERR = nf90_strerror(ISTAT)
          IF (PrintLOG) THEN
            WRITE(STAT%FHNDL,*) 'CHRERR=', TRIM(CHRERR)
          END IF
          cf_add_offset=ZERO
        ENDIF
        IF (PrintLOG) THEN
          WRITE(STAT%FHNDL,*) 'cf_add_offset=', cf_add_offset
          FLUSH(STAT%FHNDL)
        END IF

        ! Reading time

        ISTAT = nf90_inq_varid(fid, eTimeStr, varid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, ISTAT)

        ISTAT = nf90_inquire_attribute(fid, varid, "units", len=nbChar)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, ISTAT)

        ISTAT = nf90_get_att(fid, varid, "units", eStrUnitTime)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, ISTAT)
        CALL CF_EXTRACT_TIME(eStrUnitTime, ConvertToDay, eTimeStart)
        IF (PrintLOG) THEN
          WRITE(STAT%FHNDL,*) 'eTimeStart=', eTimeStart
          FLUSH(STAT%FHNDL)
        END IF

        ISTAT = nf90_inquire_variable(fid, varid, dimids=dimids)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, ISTAT)

        ISTAT = nf90_inquire_dimension(fid, dimids(1), len=nbtime_mjd)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, ISTAT)

        allocate(ice_time_mjd(nbtime_mjd), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_icenc, allocate error 48')

        ISTAT = nf90_get_var(fid, varid, ice_time_mjd)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, ISTAT)

        ISTAT = nf90_close(fid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, ISTAT)

        ice_time_mjd(:) = ice_time_mjd(:)*ConvertToDay + eTimeStart
        CALL CHECK_ICE_TIME(nbtime_mjd, ICE_TIME_MJD)
      END IF

      IF (myrank .eq. 0) THEN
        eInt(1)=nbtime_mjd
        DO IPROC=2,nproc
          CALL MPI_SEND(eInt,1,itype, iProc-1, 811, comm, ierr)
        END DO
        eReal(1)=cf_scale_factor
        eReal(2)=cf_add_offset
        DO IPROC=2,nproc
          CALL MPI_SEND(eReal,2,rtype, iProc-1, 813, comm, ierr)
        END DO
        DO IPROC=2,nproc
          CALL MPI_SEND(ICE_TIME_MJD,nbtime_mjd,rtype, iProc-1, 812, comm, ierr)
        END DO
      ELSE
        CALL MPI_RECV(eInt,1,itype, 0, 811, comm, istatus, ierr)
        nbtime_mjd=eInt(1)
        CALL MPI_RECV(eReal,2,rtype, 0, 813, comm, istatus, ierr)
        cf_scale_factor=eReal(1)
        cf_add_offset=eReal(2)
        ALLOCATE(ICE_TIME_MJD(nbtime_mjd), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_icenc, allocate error 3')
        CALL MPI_RECV(ICE_TIME_MJD,nbtime_mjd,rtype, 0, 812, comm, istatus, ierr)
      END IF

      eVAR % cf_scale_factor = cf_scale_factor
      eVAR % cf_add_offset = cf_add_offset
      eVAR % nbTime = nbtime_mjd
      allocate(eVAR % ListTime(nbtime_mjd))
      eVAR % ListTime = ICE_TIME_MJD
      eVAR % eFileName = eFileName
      eVAR % eString = eString
      eVAR % idVar = 1
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CHECK_ICE_TIME(nbtime_mjd, ICE_TIME_MJD)
      USE DATAPOOL, only : SEWI, STAT, rkind, THR, wwmerr
      IMPLICIT NONE
      integer, intent(in) :: nbtime_mjd
      real(rkind), intent(in) :: ICE_TIME_MJD(nbtime_mjd)
      CHARACTER(LEN=15) :: eTimeStr
      real(rkind) DiffTime
      real(rkind) :: tolDay = 0.00001
      IF (SEWI%BMJD .LT. minval(ICE_TIME_MJD) - tolDay) THEN
        WRITE(STAT%FHNDL,*) 'END OF RUN'
        WRITE(STAT%FHNDL,*) 'ICE START TIME is outside CF ice_time range!'
        CALL MJD2CT(SEWI%BMJD,eTimeStr)
        WRITE(STAT%FHNDL,*) 'SEWI%BMJD=', SEWI%BMJD, ' date=', eTimeStr
        CALL MJD2CT(SEWI%EMJD,eTimeStr)
        WRITE(STAT%FHNDL,*) 'SEWI%EMJD=', SEWI%EMJD, ' date=', eTimeStr
        CALL MJD2CT(minval(ICE_TIME_MJD),eTimeStr)
        WRITE(STAT%FHNDL,*) 'min(ICE_TIME_MJD)=', minval(ICE_TIME_MJD), ' date=', eTimeStr
        CALL MJD2CT(maxval(ICE_TIME_MJD),eTimeStr)
        WRITE(STAT%FHNDL,*) 'max(ICE_TIME_MJD)=', maxval(ICE_TIME_MJD), ' date=', eTimeStr
        FLUSH(STAT%FHNDL)
        WRITE(wwmerr, *) 'Error in ICE_TIME_MJD 1, read ', TRIM(STAT%FNAME)
        CALL WWM_ABORT(wwmerr)
      END IF
      IF (SEWI%EMJD .GT. maxval(ICE_TIME_MJD) + tolDay) THEN
        WRITE(STAT%FHNDL,*) 'END OF RUN'
        WRITE(STAT%FHNDL,*) 'ICE END TIME is outside CF ice_time range!'
        DiffTime=maxval(ICE_TIME_MJD) + THR - SEWI%EMJD
        WRITE(STAT%FHNDL,*) 'DiffTime=', DiffTime
        CALL MJD2CT(SEWI%BMJD,eTimeStr)
        WRITE(STAT%FHNDL,*) 'SEWI%BMJD=', SEWI%BMJD, ' date=', eTimeStr
        CALL MJD2CT(SEWI%EMJD,eTimeStr)
        WRITE(STAT%FHNDL,*) 'SEWI%EMJD=', SEWI%EMJD, ' date=', eTimeStr
        CALL MJD2CT(minval(ICE_TIME_MJD),eTimeStr)
        WRITE(STAT%FHNDL,*) 'min(ICE_TIME_MJD)=', minval(ICE_TIME_MJD), ' date=', eTimeStr
        CALL MJD2CT(maxval(ICE_TIME_MJD),eTimeStr)
        WRITE(STAT%FHNDL,*) 'max(ICE_TIME_MJD)=', maxval(ICE_TIME_MJD), ' date=', eTimeStr
        FLUSH(STAT%FHNDL)
        WRITE(wwmerr, *) 'Error in ICE_TIME_MJD 2, read ', TRIM(STAT%FNAME)
        CALL WWM_ABORT(wwmerr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GET_CF_TIME_INDEX_ICE(eVAR, REC1, REC2, w1, w2)
      ! For given wwm_time and ice_time return records to get and weights for time
      ! interpolation F(wwm_time)=F(rec1)*w1 + F(rec2)*w2
      !
      USE DATAPOOL
      IMPLICIT NONE
      TYPE(VAR_NETCDF_CF), intent(in)           :: eVAR
      REAL(rkind), INTENT(OUT)            :: w1, w2
      INTEGER, INTENT(OUT)                :: REC1, REC2
      REAL(rkind) :: eTime1, eTime2
      INTEGER  :: iTime

      DO iTime=2,eVAR % nbTime
        eTime1=eVAR % ListTime(iTime-1)
        eTime2=eVAR % ListTime(iTime)
        IF ((eTime1 .le. MAIN%TMJD).and.(MAIN%TMJD .le. eTime2)) THEN
          REC2=iTime
          REC1=iTime-1
          w2=(MAIN % TMJD - eTime1)/(eTime2-eTime1)
          w1=(eTime2 - MAIN % TMJD)/(eTime2-eTime1)
          RETURN
        END IF
      END DO
      IF (PrintLOG) THEN
        WRITE(STAT%FHNDL,*) 'Time error in wind for CF'
        WRITE(STAT%FHNDL,*) 'MAIN % TMJD=', MAIN%TMJD
        WRITE(STAT%FHNDL,*) 'min(ice_time_mjd)=', minval(eVAR % ListTime)
        WRITE(STAT%FHNDL,*) 'max(ice_time_mjd)=', maxval(eVAR % ListTime)
        FLUSH(STAT%FHNDL)
      END IF
      CALL WWM_ABORT('Error in CF ice forcing time setup')
      END SUBROUTINE
