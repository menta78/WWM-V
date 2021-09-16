!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    filename = "m0_sboltz_7r3_1_1e"      in ~/msm/TSA/subrs7r3_99_m/ !
!!    from       "sbtsa-1-norm-Dec15-08.f-7r3"  in  same  dir.         !
!!    for        "w3snlxmd.ftn-NLX-7r3-TSA2-99-1-1e"                   !
!!    used for both FBI and TSA                                        !
!!                                              Bash Toulany           !
!!    -----------------------------------------------------------------#
!!
!!op2
!!
!!              Define  nsep  as the min # of bins that separates npk & npk2
!!              and     nbins as the min # of bins > npk or npk2
!!                    to guarantee 1 bin in equi. range
!!              then    nbins1 as the actual # of bins > npk or npk2
!!
!!              nsep   = 1
!!              nbins  = 1
!!              nbins1 = nfs - npk  (or = nrng - npk2)
!!              if ( nbins1.gt.14 ) nbins1 = 14  !* to limit equi. range to 1.98*fp
!!
!!op2
!!
!!    ------------------------------------------------------------------
!!
!!-
!a    program sboltz
!!
      MODULE W3SNLXMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    20-Dec-2004 : Origination.                        ( version 3.06 )
!/    23-Jun-2006 : Formatted for submitting code for   ( version 3.09 )
!/                  inclusion in WAVEWATCH III.
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!!-=
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    Generic shallow-water Boltzmann integral for                     !
!a    "read-in" log-spaced E(f,theta) plus TSA estimate                !
!b    passed in from ww3 log-spaced action spectrum A(theta,k)         !
!!                                                                     !
!!    This is a new program for wave-wave interactions (TSA or FBI)    !
!!    This is the time-stepping version.                               !
!!                                                                     !
!!    -----------------------------------------------------------------!
!!                                                                     !
!!    Nov. 2010; It was blended inside WAVEWATCH III ver. 3.14         !
!!                                                   by Bash Toulany   !
!!    -----------------------------------------------------------------!
!!                                                                     !
!!    March/April 2005 added two-scale analysis                        !
!!                                                                     !
!!    February/March 2004 - fiddled with by Redwing                    !
!!                                                                     !
!!    March 2003 version                                               !
!!                                                                     !
!!    This is Resio version of the Boltzmann code as of March 2003.    !
!!    It should contain all changes/mods from Tracy's version.         !
!!                                                                     !
!!    Computes over the full-circle, and includes fluxes.              !
!!                                                                     !
!!    Runs best on Thursdays;                                          !
!!                 usually doesn't work too well on Fridays            !
!!                                                                     !
!!                                                                     !
!!    S h a l l o w  -  W a t e r    B o l t z m a n n    C o d e      !
!!                                                                     !
!!                                                                     !
!!    -----------------------------------------------------------------!
!!                                                                     !
!!    The included file 'parfile_2.f' (below) contains just            !
!!    parameters that sets the dimensions nrng, nang, npts and NZZ     !
!!    for some arrays (below)                                          !
!!                                                                     !
!!    parameters in 'parfile_2.f' are:                                 !
!!    nrng  =   67 = maximum number of rings                           !
!!    nang  =   36 = maximum number of angles                          !
!!    npts  =   30 = maximum number of points around locus             !
!!    NZZ   = 2278 =  NZZ = (nrng*(nrng+1))/2                          !
!!                                                                     !
!a    include 'parfile_2.f'
!!                                                                     !
!!    -----------------------------------------------------------------#
!!
!!------------------------------------------------------------------------------
!!--------------------------------------------------------------------72------80
!!==============================================================================
!!
!!+
!!
!! 1. Purpose :
!!
!!    Dummy slot for nonlinear interaction source term.
!!     This module can be used to include user-defined source term(s),
!!    or to submit a source term to be included into the WAVEWATCH III
!!    distribution. See section 5 for requirements for submissions.
!!       Codes may be included in WAVEWATCH III as a standard option
!!    with it's own dedicated compile switch, or as a version of this
!!    module that can be plugged in by a user.
!!+
!!    Here is being used as
!!    Interface module for TSA type nonlinear interactions.
!!    Based on Resio and Perrie (2008) and Perrie and Resio (2009)
!!
!! 2. Variables and types :
!!
!!     Name      Type  Scope    Description
!!    ------------------------------------------------------------------
!!    ------------------------------------------------------------------
!!
!! 3. Subroutines and functions :
!!
!!     Name      Type  Scope    Description
!!    ------------------------------------------------------------------
!!     W3SNLX    Subr. Public   Main interface for TSA subroutines.
!!                              Replaces main program "sboltz" in
!!                              "sbtsa-1-norm-Dec15-08.f" with
!!                              initialization done in subr. INSNLX
!!     INSNLX    Subr. Public   Corresponding initialization routine.
!!
!!
!!                              TSA subroutines
!!                              ---------------
!!     gridsetr  Subr. Public   Setup geometric integration grid
!!     shlocr    Subr. Public   General locus solution
!!     shloxr    Subr. Public   Locus solving routine - must converges
!!     cplshr    Subr. Public   Computes Boltzmann coupling coeff.
!!     ------
!!     optsa     Subr. Public   Converts Cart. Energy density (f,theta)
!!                              to Polar Action density (k,theta)
!!                              then split it into large and small scale
!!     snlr      Subr. Public   Computes dN(k,theta)/dt due to
!!     ------                   wave-wave inter.
!!
!!     cgf       fnc.  Public   Calculate group velocity "cgf" (m/s)
!!                              from frequency "f" (Hz), phase speed
!!                              "c" (m/s) and water depth "d" (m)
!!     wkfnc     fnc.  Public   Compute wave number "k" for given
!!     -----                    freq "f" (Hz) and water depth "d" (m)
!!    ------------------------------------------------------------------
!!
!! 4. Subroutines and functions used :
!!
!!     Name      Type  Module   Description
!!    ------------------------------------------------------------------
!!     STRACE    Subr. W3SERVMD Subroutine tracing.
!!    ------------------------------------------------------------------
!!
!! 5. Remarks :
!
!     WAVEWATCH III is designed as a highly plug-compatible code.
!     Source term modules can be included as self-contained modules,
!     with limited changes needed to the interface of routine calls
!     in W3SRCE, and in the point postprocessing programs only.
!     Codes submitted for inclusion in WAVEWATCH III should be
!     self-contained in the way described below, and might be
!     provided with distributions fully integrated in the data
!     structure, or as an optional version of this module to be
!     included by the user.
!
!     Rules for preparing a module to be included in or distributed
!     with WAVEWATCH III :
!
!      - Fully document the code following the outline given in this
!        file, and according to all other WAVEWATCH III routines.
!      - Provide a file with necessary modifications to W3SRCE and
!        all other routines that require modification.
!      - Provide a test case with expected results.
!      - It is strongly recommended that the programming style used
!        in WAVEWATCH III is followed, in particular
!          a) for readability, write as if in fixed FORTRAN format
!             regarding column use, even though all files are F90
!             free format.
!          b) I prefer upper case programming for permanent code,
!             as I use lower case in debugging and temporary code.
!
!     This module needs to be self-contained in the following way.
!
!      a) All saved variables connected with this source term need
!         to be declared in the module header. Upon acceptance as
!         permanent code, they will be converted to the WAVEWATCH III
!         dynamic data structure.
!      b) Provide a separate computation and initialization routine.
!         In the submission, the initialization should be called
!         from the computation routine upon the first call to the
!         routine. Upon acceptance as permanent code, the
!         initialization routine will be moved to a more appropriate
!         location in the code (i.e., being absorbed in ww3_grid or
!         being moved to W3IOGR).
!
!     See notes in the file below where to add these elements.
!
!!
!! 6. Switches :
!!
!!      !/S      Enable subroutine tracing.
!!
!! 7. Source code :
!/
!/ ------------------------------------------------------------------- /
!/
!     *****************************************
!     ***    Declare saved variables here   ***
!     ***  public or private as appropriate ***
!     *****************************************
!
!!
!!    ------------------------------------------------------------------
!!    ------------------------------------------------------------------
!!
!!
      PUBLIC
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!pub
!!    Declare  all gridsetr 11 returned arrays as PUBLIC arrays
      integer, allocatable, dimension(:,:,:) :: kref2, kref4     !* /fr/
      integer, allocatable, dimension(:,:,:) :: jref2, jref4     !* /fr/
      real,    allocatable, dimension(:,:,:) :: wtk2,  wtk4      !* /fr/
      real,    allocatable, dimension(:,:,:) :: wta2,  wta4      !* /fr/
      real,    allocatable, dimension(:,:,:) :: tfac2, tfac4     !* /fr/
      real,    allocatable, dimension(:,:,:) :: grad             !* /dn/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!pub
!!    Declare  all shloxr/shlocr 5 returned arrays as PUBLIC arrays
      real,    allocatable, dimension(:)     :: wk2x, wk2y       !* /a/
      real,    allocatable, dimension(:)     :: wk4x, wk4y, ds   !* /a/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!pub
!!    Declare  the only cplshr 1 returned variable as PUBLIC
!xx   real,    allocatable                   :: csq
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!pub
!!    Declare  the optsa 2 returned arrays as PUBLIC arrays
      real,    allocatable, dimension(:,:)   :: dens, dens2   !* /c/ /z/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!pub
!!    Declare  all snlr  4 returned arrays as PUBLIC arrays
!!    tsa, diag  used for -tsa
!!    fbi, diag2 used for -fbi
      real,    allocatable, dimension(:,:)   :: tsa, diag
      real,    allocatable, dimension(:,:)   :: fbi, diag2
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      CONTAINS
!!
!!------------------------------------------------------------------------------
!!--------------------------------------------------------------------72------80
!!==============================================================================
!!
!!
!!
!!
!!    ------------------------------------------------------------------
      SUBROUTINE INSNLX
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         23-Jun-2006 |
!/                  +-----------------------------------+
!/
!/    20-Dec-2004 : Origination.                        ( version 3.06 )
!/    23-Jun-2006 : Origination.                        ( version 3.09 )
!/
!  1. Purpose :
!
!     Initialization for source term routine.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SNLX    Subr. W3SNLXMD Corresponding source term.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3ODATMD, ONLY: NDSE
      USE W3SERVMD, ONLY: EXTCDE
!/S      USE W3SERVMD, ONLY: STRACE
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
!/S      INTEGER, SAVE           :: IENT = 0
!/
!/ ------------------------------------------------------------------- /
!/
!/S      CALL STRACE (IENT, 'INSNLX')
!
! 1.  .... ----------------------------------------------------------- *
!
!/
!/ End of INSNLX ----------------------------------------------------- /
!/
      END SUBROUTINE INSNLX
!!
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!
!!
!!    --------------------------------------------------------------- &
      SUBROUTINE W3SNLX ( A, CG, WN, DEPTH, NZZ,  S, D )
!!    --------------------------------------------------------------- &
!!
!!
!! 1. Purpose :
!!
!!    Interface module for TSA type nonlinear interactions.
!!    Based on Resio and Perrie (2008) and Perrie and Resio (2009)
!!
!! 2. Method :
!!
!! 3. Parameters :
!!
!!    Parameter list
!!    ------------------------------------------------------------------
!!      A       R.A.  I   Action spectrum A(NTH,NK) as a function of
!!                        direction (rad)  and wavenumber.
!!      CG      R.A.  I   Group velocities (dimension NK=nrng).
!!      WN      R.A.  I   Wavenumbers      (dimension NK=nrng).
!!      DEPTH   Real  I   Water depth in meters.
!!      NZZ     Int.  I   NZZ = (NK*(NK+1))/2
!!      S       R.A.  O   Source term.
!!      D       R.A.  O   Diagonal term of derivative.
!!    ------------------------------------------------------------------
!!
!! 4. Subroutines used :
!!
!!     Name      Type  Module   Description
!!    ------------------------------------------------------------------
!!     STRACE    Subr. W3SERVMD Subroutine tracing.
!!    ------------------------------------------------------------------
!!
!! 5. Called by :
!!
!!     Name      Type  Module   Description
!!    ------------------------------------------------------------------
!!     W3SRCE    Subr. W3SRCEMD Source term integration.
!!     W3EXPO    Subr.   N/A    Point output post-processor.
!!     GXEXPO    Subr.   N/A    GrADS point output post-processor.
!!    ------------------------------------------------------------------
!!
!! 6. Error messages :
!!
!!      None.
!!
!! 7. Remarks :
!!
!! 8. Structure :
!!
!!    See source code.
!!
!! 9. Switches :
!!
!!    !/S   Enable subroutine tracing.
!!
!!10. Source code :
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ----------------------------------------------------------------72
!!    ==================================================================
!!
!!
!b    USE CONSTANTS, ONLY: TPI, GRAV
      USE CONSTANTS, ONLY: TPI
      USE W3GDATMD,  ONLY: NK,  NTH,  XFR, DTH,  SIG, TH, ECOS, ESIN
!!    dimension: SIG(0:NK+1),TH(NTH), ECOS(NSPEC+NTH), ESIN(NSPEC+NTH)
!!
      USE W3SERVMD, ONLY: EXTCDE
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!+=
!!
!!
      IMPLICIT NONE
!!
!!+
!!    Parameter list
!!    --------------
      INTEGER, INTENT(IN)  :: NZZ
      REAL,    INTENT(IN)  :: A(NTH,NK), CG(NK), WN(NK), DEPTH
      REAL,    INTENT(OUT) :: S(NTH,NK), D(NTH,NK)
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      LOGICAL, SAVE           :: FIRST_TSA = .TRUE.
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!+=
!! 
!!    Local Parameters
!!    ----------------
!b    integer           :: nrng, nang, npts, NZZ !* where in parfile_2.f
!b    parameter          ( nrng=67, nang=36, npts=30, NZZ=2278 )
!!
      integer           :: nrng, nang, npts
!!
!!-
!b    Bash; for sbtsa run; NK, NTH & NZZ must be set here (below)
!b          nrng, nang will be set later and 
!b          npts will be hardwired to 30 as it is not likely to change
!a    integer           :: NK,   NTH
!a    parameter          ( NK = 67, NTH =36, NZZ=2278 )
!!
!b    Bash; for ww3 run; NK, NTH are coming in via USE W3GDATMD and 
!b          NZZ  is passed via the call (see above)
!b          npts will be hardwired to 30 as it is not likely to change
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!-=
!!
!!-
!!    Input variables (from Tape9 = densin)
!!    where  densin = test.out output from logspc-2d-plot-new.f
!!
!a    integer           :: npk0             !* read in from tape9
      integer           :: npk              !* re-calculated npk
!a    real              :: fpk0             !* read in from tape9
      real              :: fpk              !* re-calculated fpk
      real              :: dep              !* get it from WW3 DEPTH
      real              :: dfrq             !* get it from WW3 XFR
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      real              :: ef2(NK,NTH)      !* Cartesian Energy(f,theta)
                                            !* get it from WW3 A(NTH,NK)
!!                                          !* formely /enrgy/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!
!!    Local variables
!!
!!    to avoid cross between ierr and IERR used in WW3 change ierr to ierr_gr
      integer           :: irng,iang, ierr_gr
!a    integer           :: ifr
!wrt
!wrt  integer           :: ipt, iz
!!
      real              :: twopi            !* get it from WW3 TPI
!a    real              :: deg              !* 360./twopi
!x    real              :: pi, pi2
!x    real              :: g, gsq           !* get it from WW3 GRAV
!!
      real              :: f0               !* get it from WW3 oma(1)
      real              :: delfx            !* multiplier for bandwidth
!!                                          !* calc.  from WW3 XFR
      real              :: ainc             !* get it from WW3 DTH
!b    real              :: dwka             !* store as an array
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!a    real              :: cgf                          !* real function
!a    real              :: wkfnc                        !* real function
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      real              :: fac
      real              :: dwka(NK)                     !* k*dk*dtheta
      real              :: ef1(NK)                      !* 1D Energy
      real              :: e1max, e1sum, h1sig
!a    real              :: e2max
!!-=
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!op2
      integer           :: npk2             !* bin# of second peak freq.
      integer           :: nsep             !* # of bins that separates npk & npk2
      integer           :: nbins            !* min # of bins > npk or npk2
!!                                          !* to guarantee 1 bin in equi. range
      integer           :: nbins1           !* actual # of bins > npk or npk2
      integer           :: npeaks           !* # of peaks (=0, 1, or 2)
      integer           :: nfs, nfs2        !* bin# of freq. separation
      integer           :: nfrq             !* # of bin per freq. regime
      real              :: fpk2             !* second peak freq.
      real              :: e1max2           !* 1D energy at fpk2
      real              :: sumd1            !* sum dens+dens2 at nfs
      real              :: sumd2            !* sum dens+dens2 at nfs+1
      real              :: densat1          !* averaged dens  at nfs
      real              :: densat2          !* averaged dens  at nfs+1
      real              :: dens2sum         !* dbl-sum  dens2
      real              :: dens1sum         !* dbl-sum  dens
      real              :: dens2ov1         !* ratio of dens2sum/dens1sum
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!+
!!    Bash; Add to print
!prt  real              :: fbia1(NK,NTH),  tsaa1(NK,NTH),  afac1
!prt  real              :: fbia2(NK,NTH),  tsaa2(NK,NTH),  afac2
!prt  real              :: fbia3(NK,NTH),  tsaa3(NK,NTH),  afac3
!prt  real              :: fbia4(NK,NTH),  tsaa4(NK,NTH),  afac4
!prt  real              :: fbie1(NK,NTH),  tsae1(NK,NTH),  efac1
!prt  real              :: fbie2(NK,NTH),  tsae2(NK,NTH),  efac2
!prt  real              :: fbie3(NK,NTH),  tsae3(NK,NTH),  efac3
!prt  real              :: fbie4(NK,NTH),  tsae4(NK,NTH),  efac4
!!
!prt  real              :: fa1max, fa2max, fa3max, fa4max
!prt  real              :: fe1max, fe2max, fe3max, fe4max
!prt  real              :: ta1max, ta2max, ta3max, ta4max
!prt  real              :: te1max, te2max, te3max, te4max
!!
!prt  real              :: fa1sum, fa2sum, fa3sum, fa4sum
!prt  real              :: fe1sum, fe2sum, fe3sum, fe4sum
!prt  real              :: ta1sum, ta2sum, ta3sum, ta4sum
!prt  real              :: te1sum, te2sum, te3sum, te4sum
!!
!prt  real              :: deltaf                       !* df
!!+=
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------::--------------------------------------------72
!!    --------------------------------------------------------------- &
!!    ==================================================================
!!
      real              :: wka(NK),    oma(NK),     frqa(NK)      !* /b/
      real              :: angl(NTH),  sinan(NTH),  cosan(NTH)    !* /b/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      real              :: pha(NK)                                !* /c/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      real              :: sumint(NK,NTH)                        !* /dn/
      real              :: cga(NK),           cgnrng             !* /dn/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!pub
!!    Declare all gridsetr 11 returned arrays as PUBLIC (global) arrays
!pub  integer           :: kref2(30,NTH,NZZ), kref4(30,NTH,NZZ)  !* /fr/
!pub  integer           :: jref2(30,NTH,NZZ), jref4(30,NTH,NZZ)  !* /fr/
!pub  real              ::  wtk2(30,NTH,NZZ),  wtk4(30,NTH,NZZ)  !* /fr/
!pub  real              ::  wta2(30,NTH,NZZ),  wta4(30,NTH,NZZ)  !* /fr/
!pub  real              :: tfac2(30,NTH,NZZ), tfac4(30,NTH,NZZ)  !* /fr/
!pub  real              :: grad(30,NTH,NZZ)                      !* /dn/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!pub
!pub  real              :: dens(NK,NTH)                           !* /c/
!pub  real              :: dens2(NK,NTH)                          !* /z/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      real              :: sumintp(NK,NTH)                        !* /z/
      real              :: sumintx(NK,NTH),    sumintsa(NK,NTH)   !* /z/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!    for -tsa
!pub
!pub  real              :: tsa(NK,NTH),    diag(NK,NTH)
!!
!!    for -fbi
!pub
!pub  real              :: fbi(NK,NTH),    diag2(NK,NTH)
!!    ------------------::--------------------------------------------72
!!    ##################################################################
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!
!!-
!!    Initial constants
!!    -----------------
!!
      nrng  = NK
      nang  = NTH
      npts  = 30
!prt
!x    pi    = 3.141592654
!x    pi2   = 1.570796327
!a    twopi = 6.283185308                        !* 2pi = 6.283,185,308 rad
!a    deg   = 360./twopi
!x    g     = 9.8
!x    gsq   = 96.04
!!
      twopi = TPI                                !* twopi = 8.*atan(1.)
!x    pi    = twopi/2.
!x    pi2   = pi/2.
!x    g     = GRAV
!x    gsq   = GRAV**2
!!-=
!!    ------------------------------------------------------------------
!!
!!+
!!    Initialization of the output arrays
!!    before calling TSA subroutines.
!!    -----------------------------------
      do 12 iang=1,nang
        do 11 irng=1,nrng
          S(iang,irng) = 0.0
          D(iang,irng) = 0.0
  11    continue
  12  continue
!!+=
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!-
!!
!!-----------------------------------------------------------------------------#
!!                                                                             !
!!    dep  = [m] depth of water                                                !
!!    dfrq = frequency multiplier for logarithmic frequency spacing = 1.05     !
!!    npk  = the magic number 14 ???  (as read from tape9)                     !
!!    fpk  = [Hz] peak frequency of initial frequency spectrum and             !
!!              a parameter of the solution grid                               !
!!                                                                             !
!!-----------------------------------------------------------------------------#
!!
!!    open & read spectral ef2(f,theta) Energy fit parameters (fit large-scale)
!!    & separate ef2(f,theta) into large-scale (dens1) and small-scale (dens2)
!!
!!    open & read in tape9 = e(f,theta)
!a    open(9,file='densin',status='old')
!!
!a    read (9,*) fpk0,dep,dfrq,npk0
!a    do 110 irng=1,nrng
!a       read(9,'(36e11.3)') (ef2(irng,iang), iang=1,nang)
!110  continue
!!
!a    close(9)
!!-=
!!
      dep  = DEPTH
      dfrq = XFR
!!
!!    still have to find npk and fpk of the 1D Energy spectrum (see below)
!!    --------------------------------------------------------------------------
!!
!!-
!!    check input value of dfrq
      if (dfrq .le. 1.) then
         print *, ' dfrq </= 1;  no longer supported '
         print *, ' use geometric spacing - dfrq > 1 '
!a       stop
         CALL EXTCDE ( 110 )
      end if
!!-=
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!-
!!    open all output files
!a    open(80, file='2d-energy.dat', status='unknown')
!a    open(81, file='2d-Snl-e.dat',  status='unknown')
!a    open(84, file='2d-TSA-e.dat',  status='unknown')
!a    open(85, file='2d-energy-norm.dat',status='unknown')
!!
!a    open(11, file='freq.dat', status='unknown')
!a    open(12, file='direc.dat', status='unknown')
!a    open(40, file='2D_Energy.dat', status='unknown')
!a    open(41, file='2D_FBI_Snl.dat',  status='unknown')
!a    open(44, file='2D_TSA_Snl.dat',  status='unknown')
!!
!!    Bash; Add to print
!b    open(16, file='dbl_sum.dat',   status='unknown')      !* use Tape6
!!-=
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!-
!!    Initialize frequency arrays and related parameters
!!    --------------------------------------------------
!a    f0    = fpk0 * (dfrq**(-(npk0-1)))          !* npk=14 is magic - no more
      delfx = 0.5 * (dfrq-1./dfrq)                !* multiplier for bandwidth
!!                                    !* needed for  deltaf = delfx*frqa(irng)
!a    do 115 irng=1,nrng                          !* on each computational ring:
!a       frqa(irng) = f0 * (dfrq**(irng-1))       !* frequency f [Hz}
!a       wka(irng)  = wkfnc ( frqa(irng), dep )   !* wavenumber modulus k[rad/m]
!a       oma(irng)  = twopi * frqa(irng)          !* frequency omega [rad/s]
!a       cga(irng)  = cgf ( frqa(irng), dep, oma(irng)/wka(irng) )
!!                                                !* oma/wka = phase speed
!115  continue
!!    --------------------------------------------------------------------------
!!
      oma    = SIG(1:NK)
      frqa   = oma / twopi
      f0     = frqa(1)
!!
      wka    = WN(1:NK)
      cga    = CG(1:NK)
      cgnrng = cga(nrng)                          !* group velocity, outer ring
!!-=
!!------------------------------------------------------------------------------
!!
!!
!!-
!!    Initialize direction arrays and related parameters
!!    --------------------------------------------------
!a    ainc  = twopi/nang                          !* angle increment (radians)
!a    do 2 iang=1,nang
!a       angl(iang)  = (iang-1)*ainc              !* angle array (radians)
!a       cosan(iang) = cos(angl(iang))            !* array of cosine(angle)
!a       sinan(iang) = sin(angl(iang))            !* array of sine(angle)
!a 2  continue
!!    --------------------------------------------------------------------------
!!
      ainc   = DTH
      angl   = TH(1:NTH)
      cosan  = ECOS(1:NTH)
      sinan  = ESIN(1:NTH)
!!-=
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    initialize array pha = k*dk*dtheta, the base area at a grid      !
!!    intersection for use in integration of 2-d density functions;    !
!!    dk at a given ring here is the spacing between adjacent cell     !
!!    centers, with edge points consistent with the geometric          !
!!    spacing of frequency rings                                       !
!!                                                                     !
!!    -----------------------------------------------------------------#
!!
!!-
!!    In the line below "wkfnc(frqa(1)/dfrq,dep)" is like wka(0)
      dwka(1) = (wka(2) - wkfnc(frqa(1)/dfrq,dep)) / 2.
!!                                              !* dk at ring 1
      pha(1)  = wka(1)*dwka(1)*ainc             !* k*dk*dtheta at ring 1
!!
      do 21 irng=2,nrng-1
        dwka(irng) = (wka(irng+1)-wka(irng-1)) / 2. !* dk at irng
        pha(irng)  = wka(irng)*dwka(irng)*ainc    !* k*dk*dtheta at irng
  21  continue
!!
!!    In the line below "wkfnc(dfrq*frqa(nrng),dep)" is like wka(nrng+1)
      dwka(nrng) = (wkfnc(dfrq*frqa(nrng),dep) - wka(nrng-1)) / 2.
!!                                                !* dk at nrng
      pha(nrng)  = wka(nrng)*dwka(nrng)*ainc      !* k*dk*dtheta at nrng
!!-=
!!------------------------------------------------------------------------------
!!
!!-
!!    -----------------------------------------------------------------#
!!    Since this run is for a rectangular ocean with uniform depth     |
!!    call grid geometry subr "gridsetr" only once (from here) for     |
!!    all spectra at all grid points and for all the time steps.       |
!!                                                                     |
!!-   We can also do cga and wka (both fnc of frqa and dep) once (here)|
!!-   and since they are used in gridsetr we have to do them first.    |
!!-   However, for now ignore this option.                             |
!!    -----------------------------------------------------------------#
!!
!!
      IF ( FIRST_TSA ) THEN
!!
!pub
!!     allocate all gridsetr 11 returned arrays declared above as PUBLIC
       allocate(kref2(npts, nang, NZZ))
       allocate(kref4(npts, nang, NZZ))

       allocate(jref2(npts, nang, NZZ))
       allocate(jref4(npts, nang, NZZ))

       allocate(wtk2(npts, nang, NZZ))
       allocate(wtk4(npts, nang, NZZ))

       allocate(wta2(npts, nang, NZZ))
       allocate(wta4(npts, nang, NZZ))

       allocate(tfac2(npts, nang, NZZ))
       allocate(tfac4(npts, nang, NZZ))

       allocate(grad(npts, nang, NZZ))
!!
!pub
!!     allocate all shloxr/shlocr 5 returned arrays as PUBLIC arrays
       allocate(wk2x(npts))
       allocate(wk2y(npts))
       allocate(wk4x(npts))
       allocate(wk4y(npts))
       allocate(ds(npts))
!!
!pub
!!     allocate the only cplshr 1 returned variable as PUBLIC
!xx    can't declare as allocatable sibgle variable ?!
!xx    allocate(csq)
!!
!pub
!!     allocate the optsa 2 returned arrays as PUBLIC arrays
       allocate(dens(nrng, nang))
       allocate(dens2(nrng, nang))
!!
!pub
!!     allocate all snlr  4 returned arrays as PUBLIC arrays
!!     tsa, diag  used for -tsa
!!     fbi, diag2 used for -fbi
       allocate(tsa(nrng, nang))
       allocate(diag(nrng, nang))
       allocate(fbi(nrng, nang))
       allocate(diag2(nrng, nang))
!!
!pub
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!    set up integration grid
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!
!!    --------------------------------------------------------------- &
!pub  call gridsetr ( nrng,nang,npts,NZZ, dep,dfrq,ierr_gr,           &
!pub                  frqa,wka, cgnrng, sinan,cosan,                  &
!pub                  kref2,kref4,jref2,jref4,wtk2,wtk4,              &
!pub                  wta2,wta4,tfac2,tfac4,grad )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call gridsetr ( nrng,nang,npts,NZZ, dep,dfrq,ierr_gr,           &
                      frqa,wka, cgnrng, sinan,cosan )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!
!!                    wka,    frqa,      sinan,cosan,    !* from /b/
!!                    grad,           cgnrng,            !* from /dn/
!!                    kref2,kref4,jref2,jref4,wtk2,wtk4, !* =    /fr/
!!                    wta2,wta4,tfac2,tfac4 )            !* =    /fr/
!!                    from  /b/, /dn/, and /fr/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
      if (ierr_gr .ne. 0) then
         write(*,'('' error '',i10,'' from gridset; bail'')') ierr_gr
!a       stop
         CALL EXTCDE ( 120 )
      end if
!!
      FIRST_TSA  = .FALSE.
      print *, ' Done calling gridsetr OK '
!!-=
!!    ------------------------------------------------------------------
!!    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!
!!
!!
!!wrt   Bash; Add test write output kref2()
!wrt    open(103, file='kref2_at00.dat', status='unknown')    !* Tape103
!wrt    write(103,903) 'kref2(npts,NTH,NZZ)'
!wrt    do 614 ipt =1,npts
!wrt    do 613 iang=1,nang
!wrt       write(103,904) ipt, iang
!wrt       write(103,940) (kref2(ipt,iang,iz), iz=1,NZZ)
!613    continue
!614    continue
!wrt    close(103)
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!wrt   Bash; Add test write output kref4()
!wrt    open(104, file='kref4_at00.dat', status='unknown')    !* Tape104
!wrt    write(104,903) 'kref4(npts,NTH,NZZ)'
!wrt    do 616 ipt =1,npts
!wrt    do 615 iang=1,nang
!wrt       write(104,904) ipt, iang
!wrt       write(104,940) (kref4(ipt,iang,iz), iz=1,NZZ)
!615    continue
!616    continue
!wrt    close(104)
!!      == = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!!      == = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!!
!!
!!wrt   Bash; Add test write output jref2()
!wrt    open(105, file='jref2_at00.dat', status='unknown')    !* Tape105
!wrt    write(105,903) 'jref2(npts,NTH,NZZ)'
!wrt    do 618 ipt =1,npts
!wrt    do 617 iang=1,nang
!wrt       write(105,904) ipt, iang
!wrt       write(105,940) (jref2(ipt,iang,iz), iz=1,NZZ)
!617    continue
!618    continue
!wrt    close(105)
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!wrt   Bash; Add test write output jref4()
!wrt    open(106, file='jref4_at00.dat', status='unknown')    !* Tape106
!wrt    write(106,903) 'jref4(npts,NTH,NZZ)'
!wrt    do 620 ipt =1,npts
!wrt    do 619 iang=1,nang
!wrt       write(106,904) ipt, iang
!wrt       write(106,940) (jref4(ipt,iang,iz), iz=1,NZZ)
!619    continue
!620    continue
!wrt    close(106)
!!      == = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!!      == = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!!
!!
!!wrt   Bash; Add test write output wtk2()
!wrt    open(107, file='wtk2_at00.dat', status='unknown')     !* Tape107
!wrt    write(107,903) 'wtk2(npts,NTH,NZZ)'
!wrt    do 622 ipt =1,npts
!wrt    do 621 iang=1,nang
!wrt       write(107,904) ipt, iang
!wrt       write(107,950) (wtk2(ipt,iang,iz), iz=1,NZZ)
!621    continue
!622    continue
!wrt    close(107)
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!wrt   Bash; Add test write output wtk4()
!wrt    open(108, file='wtk4_at00.dat', status='unknown')     !* Tape108
!wrt    write(108,903) 'wtk4(npts,NTH,NZZ)'
!wrt    do 624 ipt =1,npts
!wrt    do 623 iang=1,nang
!wrt       write(108,904) ipt, iang
!wrt       write(108,950) (wtk4(ipt,iang,iz), iz=1,NZZ)
!623    continue
!624    continue
!wrt    close(108)
!!      == = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!!      == = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!!
!!
!!wrt   Bash; Add test write output wta2()
!wrt    open(109, file='wta2_at00.dat', status='unknown')     !* Tape109
!wrt    write(109,903) 'wta2(npts,NTH,NZZ)'
!wrt    do 626 ipt =1,npts
!wrt    do 625 iang=1,nang
!wrt       write(109,904) ipt, iang
!wrt       write(109,950) (wta2(ipt,iang,iz), iz=1,NZZ)
!625    continue
!626    continue
!wrt    close(109)
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!wrt   Bash; Add test write output wta4()
!wrt    open(110, file='wta4_at00.dat', status='unknown')     !* Tape110
!wrt    write(110,903) 'wta4(npts,NTH,NZZ)'
!wrt    do 628 ipt =1,npts
!wrt    do 627 iang=1,nang
!wrt       write(110,904) ipt, iang
!wrt       write(110,950) (wta4(ipt,iang,iz), iz=1,NZZ)
!627    continue
!628    continue
!wrt    close(110)
!!      == = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!!      == = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!!
!!
!!wrt   Bash; Add test write output tfac2()
!wrt    open(111, file='tfac2_at00.dat', status='unknown')    !* Tape111
!wrt    write(111,903) 'tfac2(npts,NTH,NZZ)'
!wrt    do 630 ipt =1,npts
!wrt    do 629 iang=1,nang
!wrt       write(111,904) ipt, iang
!wrt       write(111,950) (tfac2(ipt,iang,iz), iz=1,NZZ)
!629    continue
!630    continue
!wrt    close(111)
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!wrt   Bash; Add test write output tfac4()
!wrt    open(112, file='tfac4_at00.dat', status='unknown')    !* Tape112
!wrt    write(112,903) 'tfac4(npts,NTH,NZZ)'
!wrt    do 632 ipt =1,npts
!wrt    do 631 iang=1,nang
!wrt       write(112,904) ipt, iang
!wrt       write(112,950) (tfac4(ipt,iang,iz), iz=1,NZZ)
!631    continue
!632    continue
!wrt    close(112)
!!      == = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!!      == = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!!
!!
!!wrt   Bash; Add test write output grad()
!wrt    open(113, file='grad_at00.dat', status='unknown')     !* Tape113
!wrt    write(113,903) 'grad(npts,NTH,NZZ)'
!wrt    do 634 ipt =1,npts
!wrt    do 633 iang=1,nang
!wrt       write(113,904) ipt, iang
!wrt       write(113,950) (grad(ipt,iang,iz), iz=1,NZZ)
!633    continue
!634    continue
!wrt    close(113)
!!      == = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!!      == = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!!
!!
!903    format(2x,A)
!904    format(2x,'ipt=',i4,2x,'iang=',i4)
!940    format(2x,10I5)
!950    format(2x,10E12.5)
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!      == = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!!      == = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!!
!!
!!dens  Bash; Add test write output ef2()
!b      open(123, file='ef2_2pks.dat', status='unknown')      !* Tape123
!b      write(123,903) 'ef2(NK,NTH) + 2 peaks'
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!dens  Bash; Add test write output dens()
!b      open(124, file='dens_2pks.dat', status='unknown')     !* Tape124
!b      write(124,903) 'dens(NK,NTH) + 2 peaks'
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!      == = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!!      == = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!!
!!
      END IF
!!    END IF ( FIRST_TSA )
!!-=
!!    ------------------------------------------------------------------
!!    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!
!!+
!!    Convert input WW3 Cartesian 2D Action density spectrum   A(theta,k)
!!                   to Cartesian 2D Energy density spectrum ef2(theta,f)
!!                                   and  reverse indices to ef2(f,theta)
!!    ==>  ef2(f,theta) = A(theta,k) * 2*pi*oma(k)/cga(k)
!!         It's this ef2(f,theta) that I send to optsa
!!         But first I used it to calc. 1D Energy ef1(f)
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      do 32 irng=1,nrng
        fac = twopi*oma(irng)/cga(irng)
        do 31 iang=1,nang
          ef2(irng,iang) = A(iang,irng) * fac
  31    continue
  32  continue
!!
!!+
!!    Convert input WW3 Cartesian 2D Action density spectrum   A(theta,k)
!!            2)     to Polar     2D Action density spectrum af2(theta,k)
!!                                   and  reverse indices to af2(k,theta)
!!    ==>  af2(k,theta) = A(theta,k) * 1.0*/wka(k)
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!b    do 29 irng=1,nrng
!b      fac = 1.0 / wka(irng)
!b      do 28 iang=1,nang
!b        af2(irng,iang) = A(iang,irng) * fac
!b28    continue
!b29  continue
!!
!!
!!    First e2max = max ef2(f,theta)
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!a    e2max = 0.0
!a    do 34 irng=1,nrng
!b    do 34 irng=1,28       !* over first 28 freq. only
!a      do 33 iang=1,nang
!a        if ( ef2(irng,iang).gt.e2max ) e2max = ef2(irng,iang)
!a33    continue
!a34  continue
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!
!!
!!    First calculte the 1D Energy density "ef1(f)"
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      do 42 irng=1,nrng
        ef1(irng) = 0.0
        do 41 iang=1,nang
          ef1(irng) = ef1(irng) + ef2(irng,iang)
  41    continue
        ef1(irng) = ef1(irng) * ainc           !* sum(ef2(f,theta))*ainc
  42  continue
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!!
!!                                                   op2-99c
!!
!!op2
!!
!!*   Bash;
!!*   Find 1 or 2 peaks that satisfy TSA min condition (below) ------- *
!!*   before calling TSA subrs. otherwise bailout (return) ----------- *
!!*   otherwise return with init. values of S & D = 0.0 -------------- *
!!
!!*   nsep   = # of bins that separates between npk & npk2
!!*   nbins  = min # of bins (incl nrng) > npk
!!*    or      min # of bins (incl  nfs) > npk2 to guarantee 1 bin in equi. range
!!*   nbins1 = actual # of bins > npk or nrng; nbins1 = nfs - npk (or nrng - npk2)
!!
!!*   ===>  TSA min condition relative to nrng is satisfied when ----- * 
!!*   ===>  npk.le.nrng-nbins guarantee min nbins (incl nrng) > fpk -- *  <<<<<
!!*   ===>  of which we only use 1 in optsa2 equi. range defined ----- *  <<<<<
!!*   ===>  to be between [dfrq**(nbins) -/+ 0.025] * fp  ------------ *  <<<<<
!!*   ===>  skip if condition is not met ie if  npk.gt.nrng-nbins ---- *  <<<<<
!!
!!*   ===>  In case of 2 peaks we have 2nd condition to satifies ----- * 
!!*   ===>  the TSA min condition for the low freq peak -------------- * 
!!
!!*   ===>  TSA min condition relative to nfs  is satisfied when ----- * 
!!*   ===>  npk.le.nfs-nbins guarantee min nbins (incl nfs) > fpk2 --- *  <<<<<
!!*   ===>  of which we only use 1 in optsa2 equi. range defined ----- *  <<<<<
!!*   ===>  to be between dfrq**(nbins-1)*fp & dfrq**(nbins+1)*fp ---- *  <<<<<
!!*   ===>  skip if condition is not met ie if  npk2.gt.nfs-nbins ---- *  <<<<<
!!*
!!*   ===>  with:
!!*         iabs(npk-npk2) > nsep (min nsep bins separating the 2 peaks)  <<<<<
!!
!!    With dfrq = 1.05, for a combination of (nsep, nbins);
!!                      the equi. range used in optsa2 is set to
!!*   ===>  nsep = 1, nbins = 1  ==>  1.02*fp  &  1.07*fp
!!*   ===>  nsep = 3, nbins = 2  ==>  1.07*fp  &  1.13*fp
!!*   ===>  nsep = 5, nbins = 3  ==>  1.13*fp  &  1.18*fp
!!*         o o o
!!*   ===>  nsep = 3, nbins =14  ==>  1.95*fp  &  2.00*fp
!!*   ===>  nsep = 5, nbins =14  ==>  1.95*fp  &  2.00*fp
!!
!!    The orig equi. range was  1.55*fp  &  2.45*fp;
!!    The new  equi. range is   Narrow and close to fp ???
!!    ------------------------------------------------------------------
!!    ==================================================================
!!    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!    ==================================================================
!!
!!
!!op2
!!    Define nsep & nbins
      nsep  = 1
      nbins = 1
!!
!!
!!op2
!!    First find the overall peak with e1max must be > 0.000001
!!    Starting from low freq. find the Energy max "e1max" and
!!    corresp. peak freq. "fpk" and its freq. number "npk".
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    sum energy over the 1st bin + the last "nbins-1"  excluded from search
!!    The 1st bin is before 1st  possible peak loc. at 2
!!    last nbins are after  last possible peak loc. at nrng-nbins         <<<<<
      e1sum = ef1(1)*delfx*frqa(1)
      do 40 irng=nrng-nbins+1,nrng                   !* nrng-nbins+1,nrng <<<<<
        e1sum = e1sum + ef1(irng) * delfx*frqa(irng) !* sum E(f)*df
  40  continue
!!
      e1max  = 0.0
      fpk    = 0.0
      npk    = 0
      npeaks = 0
!!
!!    look in the freq range that works for TSA call (see condition below)
      do 43 irng=2,nrng-nbins          !* last peak loc. is at nrng-nbins <<<<<
!!
        e1sum = e1sum + ef1(irng) * delfx*frqa(irng) !* sum E(f)*df
!!
!b      if ( ef1(irng).ge.ef1(irng-1) .and. ef1(irng).ge.ef1(irng+1)  &
!b                                    .and. ef1(irng).ge.e1max ) then
!!      pick the abs. local peak, if 2 equal peaks are found take the lower freq.
!!      take the lower freq. and update npk, fpk & e1max
        if ( ef1(irng).gt.ef1(irng-1) .and. ef1(irng).gt.ef1(irng+1)  &
                                      .and. ef1(irng).gt.e1max ) then
          npk    = irng                      !* update npk
          fpk    = frqa(npk)                 !* update fpk
          e1max  = ef1(npk)                  !* update e1max
          npeaks = 1
        endif
!!
  43  continue
!!B   if a peak is not found then of course  e1max = 0.0 < eps
!!    ------------------------------------------------------------------
!!
!!
!!B   if a peak is not found (npeaks=0 & e1max=0.0<eps)
!!B   if a peak is found with a tiny peak   (e1max<eps) or
!!B   if a peak is found with a tiny energy (e1sum<eps) or
!!B   if TSA min condition is not met rel. to nrng (npk.gt.nrng-nbins)    <<<<<
!!B   drop everything & return with init. values of S() and D() arrays=0.0
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      if ( e1max.lt.0.000001 .or. e1sum.lt.0.000001 .or.              &
                                  npk.gt.nrng-nbins )   then          !*  <<<<<
!!      --------------------------------------------------------------!*  <<<<<
!!
!!*     Not suitable spectrum for tsa; don't call tsa routines ------- *
!!*     print info line, skip and return (don't stop) ---------------- *
!!      --------------------------------------------------------------!*  <<<<<
        write(6,206) nbins, npk,fpk,e1max,e1sum                       !*  <<<<<
 206    format(' *** W3SNLX: e1max,e1sum< or npk>nrng-',I1,           &
                                                       ' skip tsa;'   &
               ' npk, fpk, e1max,e1sum =', I4,f9.6,E14.6,f10.4)
!!
        return     !* return from here
!!                 !* with init. values of S() and D() arrays = 0.0
!!
      else
!!
        h1sig  = 4.0 * sqrt(e1sum)
!!
      endif
!!    ------------------------------------------------------------------
!!    ==================================================================
!!    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!    ==================================================================
!!
!!
!!op2
!!    Bash; if we are here (i.e. we did not return) then we must
!!    have found the 1st good peak (= overall peak) with e1max>eps
!!    then look for new 2nd peak at least 'nsep' bins away from the 1st   <<<<<
!!                                          iabs(irng-npk).gt.nsep        <<<<<
!!    before calling new "optsa2" (thw 2nd peak will have e1max2 < e1max)
!!    Again look in the freq range that is in line with TSA min condition
!!    and find the 2nd highest peak with  eps < e1max2 < e1max
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!b    if ( npeaks.eq.1 .and. e1max.gt.0.000001 ) then
        e1max2 = 0.0
        fpk2   = 0.0
        npk2   = 0
        do 45 irng=2,nrng-nbins        !* last peak loc. is at nrng-nbins <<<<<
!!
!b        if ( ef1(irng).ge.ef1(irng-1).and.ef1(irng).ge.ef1(irng+1)  &
!b        .and.ef1(irng).ge.e1max2.and.iabs(irng-npk).gt.nsep ) then  !*  <<<<<
!!
          if ( ef1(irng).gt.ef1(irng-1).and.ef1(irng).gt.ef1(irng+1)  &
          .and.ef1(irng).ge.e1max2.and.iabs(irng-npk).gt.nsep ) then  !*  <<<<<
!!        ------------------------------------------------------------!*  <<<<<
!!
!!          pick new abs. local peak, if 2 equal peaks are found 
!!          take the higher freq. and update npk2, fpk2 & e1max2
!!          then continue the loop to find the biggest e1max2  peak
            npk2   = irng                      !* update npk2
            fpk2   = frqa(npk2)                !* update fpk2
            e1max2 = ef1(npk2)                 !* update e1max2
            npeaks = 2
          endif
  45    continue
!b    endif
!!B   if 2nd peak is not found then we have single peak spectrum (npeaks=1)
!!    ------------------------------------------------------------------
!!
!!op2
!!B   if 2nd peak is found then make sure:
!!B   (1) its e1max2>eps (not tiny) and
!!B   (2) it satisfies TSA min condition relative to nrng
!!B   if not drop this 2nd peak & stick to the overall one peak spectrum
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (  npeaks.eq.2 .and.                                         &
           (e1max2.lt.0.000001 .or. npk2.gt.nrng-nbins) )  npeaks = 1 !*  <<<<<
!!    ----------------------------------------------------------------!*  <<<<<
!!
!!op2 op2 op2
!!B   if we still have 2 peaks (npeaks=2)
!!B   then find the bin in the middle to divide the freq. regime into two
      if ( npeaks.eq.2 ) then
        nfs = int ( (npk+npk2+1) / 2.0 )  !* give the higher bin # to nfs
      endif
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!op2 op2 op2
!!B   if we still have 2 peaks (npeaks=2)
!!B   then check and make sure that:
!!B   the lower freq peak satisfies TSA min condition relative to nfs     <<<<<
!!B                                                   npk2.le.nfs-nbins   <<<<<
!!B   if not drop the idea of 2 freq regimes & go back to 1 overall peak
!!B   which we already know it satisfies TSA min condition rel. to nrng
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      if ( npeaks.eq.2 ) then
!!
!!      we have 2 peaks with npk the lager peak & npk2 the smaller peak
!!
!!      --------------------------------------------------------------!*  <<<<<
        if ( npk2.lt.npk  .and. npk2.gt.nfs-nbins ) then              !*  <<<<<
!!B       low freq (small) peak npk2 failed TSA min condition relative to nfs
          write(6,216) npk2,fpk2,e1max2,nfs
 216      format(' *** W3SNLX: low  freq sml peak; failed nfs test ;' &
                 ' npk2,fpk2,e1max2,nfs  =', I4,f9.6,E14.6,6x,I4)
          write(6,217) npk,fpk,e1max,h1sig
 217      format(' ***         high freq lrg peak   is at ;         ' &
                 ' npk, fpk, e1max,h1sig =', I4,f9.6,E14.6,f10.4)
          npeaks = 1
        endif
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!      --------------------------------------------------------------!*  <<<<<
        if ( npk .lt.npk2 .and. npk .gt.nfs-nbins ) then              !*  <<<<<
!!B       low freq (big) peak npk failed TSA min condition relative to nfs
          write(6,218) npk,fpk,e1max,nfs
 218      format(' *** W3SNLX: low  freq lrg peak; failed nfs test ;' &
                 ' npk, fpk, e1max,nfs   =', I4,f9.6,E14.6,6x,I4)
          write(6,219) npk2,fpk2,e1max2,h1sig
 219      format(' ***         high freq sml peak   is at ;         ' &
                 ' npk2,fpk2,e1max2,h1sig=', I4,f9.6,E14.6,f10.4)
          npeaks = 1
        endif
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      endif
!!    ------------------------------------------------------------------
!!
!!op2
!!B   Here if we still have 2 peaks the 1st is the dominant peak,
!!B   re-order the 2 peaks forcing npk to be always < npk2
!!B   and  re-alligne corresp. fpk,e1max & fpk2,e1max2
!!B   This says nothing about which peak is the dominant peak after the shift
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if ( npeaks.eq.2 .and. npk2.lt.npk ) then
        nfs2   = npk2
        npk2   = npk
        npk    = nfs2
        fpk    = frqa(npk)
        fpk2   = frqa(npk2)
        e1max  = ef1(npk)
        e1max2 = ef1(npk2)
      endif                         !*  this way  npk < npk2  always
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!    ==================================================================
!!
!!
!!op2
!!    Bash; Using the new "optsa2" you are allowed
!!          a single call or possibly 2 calls (if 2 peaks)
!!          to account for spectra with double peaks.
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      if ( npeaks.eq.1 ) then
!!
!!      Do this this  call for the only one peak  (if 1 peak)
!!
!!-0    print info line about the 1st peak and continue -------------- *
!!-0    this condition guarantees a min. of 4 bins in the equi. range
        write(6,207) npk,fpk,e1max,h1sig
 207    format(' *** W3SNLX: continue - one peak; call tsa subrs ;'   &
               ' npk, fpk, e1max,h1sig =', I4,f9.6,E14.6,f10.4)
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!-0    one call to optsa2 for the whole freq. regime ( 1 --> nrng )
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        nfrq   = nrng                    !* = nrng - 1 +1
        nbins1 = nrng - npk
        if ( nbins1.gt.14 ) nbins1 = 14  !* to limit equi. range to 1.98*fp
        call optsa2 ( 1,nrng,     nrng,nang,   npk, fpk,   frqa,      &
                      oma,wka,cga, angl,cosan, ef1,ef2, dfrq,nbins1 )
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      endif
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!!op2
      if ( npeaks.eq.2 ) then
!!
!!      Now make two calls to new "optsa2" one for each freq regime.
!!
!!-1    1st call to optsa2 for the lower freq. regime ( 1 --> nfs )
!!-1    print info line about the 1st peak and continue -------------- *
!!-1    this condition guarantees a min. of 4 bins in the equi. range
        write(6,208) npk,fpk,e1max,h1sig
 208    format(' *** W3SNLX: low  freq  1st peak; call tsa subrs ;'   &
               ' npk, fpk, e1max,h1sig =', I4,f9.6,E14.6,f10.4)
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        nfrq   = nfs                     !* = nfs - 1 +1
        nbins1 = nfs - npk
        if ( nbins1.gt.14 ) nbins1 = 14  !* to limit equi. range to 1.98*fp
        call optsa2 ( 1,nfs,      nrng,nang,   npk, fpk,   frqa,      &
                      oma,wka,cga, angl,cosan, ef1,ef2, dfrq,nbins1 )
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!
!!
!!-2    2nd call to optsa2 for the higher freq. regime ( nfs+1 --> nrng )
!!-2    print info line about the 2nd peak and continue -------------- *
!!-2    this condition guarantees a min. of 4 bins in the equi. range
        write(6,209) npk2,fpk2,e1max2,nfs
 209    format(' ***         high freq  2nd peak; call tsa subrs ;'   &
               ' npk2,fpk2,e1max2, nfs =', I4,f9.6,E14.6,6x,I4)
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        nfrq   = nrng - nfs              !* = nrng - (nfs+1) + 1
        nbins1 = nrng - npk2
        if ( nbins1.gt.14 ) nbins1 = 14  !* to limit equi. range to 1.98*fp
        call optsa2 ( nfs+1,nrng, nrng,nang,   npk2,fpk2,  frqa,      &
                      oma,wka,cga, angl,cosan, ef1,ef2, dfrq,nbins1 )
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!
!!
!!-3    Remove the step like jump (if exist) in dens() between nfs & nfs+1
        do 440 iang=1,nang
          sumd1 = dens(nfs,iang)   + dens2(nfs,iang)    !* sum at nfs
          sumd2 = dens(nfs+1,iang) + dens2(nfs+1,iang)  !* sum at nfs+1
!!
!!        do 3 bin average dens() at nfs   and store in densat1
          densat1 = ( dens(nfs-1,iang) + dens(nfs,iang) +             &
                                         dens(nfs+1,iang) ) / 3.
!!        do 3 bin average dens() at nfs+1 and store in densat2
          densat2 = ( dens(nfs,iang)   + dens(nfs+1,iang) +           &
                                         dens(nfs+2,iang) ) / 3.
!!
!!        subtitute back into dens(nfs,iang) & dens(nfs+1,iang)
          dens(nfs,iang)    = densat1             ! dens at nfs
          dens(nfs+1,iang)  = densat2             ! dens at nfs+1
!!
!!        recalculate dens2(nfs,iang) & dens2(nfs+1,iang)
          dens2(nfs,iang)   = sumd1 - densat1     ! dens2 at nfs
          dens2(nfs+1,iang) = sumd2 - densat2     ! dens2 at nfs+1
 440    continue
!!
      endif
!!    ------------------------------------------------------------------
!!    ==================================================================
!!    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!    ==================================================================
!!
!!
!!    Calculate dbl-sum of small-scale Action density dens2(k,theta) 
!!    and       dbl-sum of large-scale Action density dens(k,theta) 
!!    then      calculate their ratio  dens2ov1 and write it out
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      dens2sum = 0.0
      dens1sum = 0.0
      dens2ov1 = 0.0
      do 442 irng=1,nrng
        do 441 iang=1,nang
          dens2sum = dens2sum + abs(dens2(irng,iang))
          dens1sum = dens1sum +      dens(irng,iang)
 441    continue
        dens2sum = dens2sum * ainc * delfx*frqa(irng)
        dens1sum = dens1sum * ainc * delfx*frqa(irng)
 442  continue
      dens2ov1 = dens2sum / dens1sum
!!
      write(6,418) dens2ov1
 418  format(' *** W3SNLX: dbl-sum(abs(dens2)) / dbl-sum(dens) ;'     &
      ' ------------ dens2ov1 =', f11.4,4x,'----------------------')
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!
!!
!!op2
!!
!b      if ( npeaks.eq.2 ) then
!!
!!dens    Now after making 2 calls to "optsa2" (if 2 peaks)
!!dens    write the 2D energy ef2() which have 2 peaks to Tape123
!b        write(123,210) npk,fpk,e1max, npk2,fpk2,e1max2, h1sig
!b        do 46 iang=1,nang
!b          write(123,211) (ef2(irng,iang), irng=1,nrng)
!b46      continue
!!        -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!dens    write the 2D broad scale 'dens() which have 2 peaks to Tape124
!b        write(124,210) npk,fpk,e1max, npk2,fpk2,e1max2, h1sig
!b        do 47 iang=1,nang
!b          write(124,211) (dens(irng,iang), irng=1,nrng)
!b47      continue
!!        -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!210      format(1x,/,' npk, fpk, e1max       =',I4,f9.6,E14.6,/,     &
!b                    ' npk2,fpk2,e1max2,h1sig=',I4,f9.6,E14.6,f10.4)
!211      format(1x,10E11.4)
!!        -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!b      endif
!!    ------------------------------------------------------------------
!!    ==================================================================
!!    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!    ==================================================================
!!
!!
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    Get Snl source term and its diagonal term  from "snlr"           !
!!                                                                     !
!!    -----------------------------------------------------------------#
!!
!!
!!
!!    ==================================================================
!!    ------------------------------------------------------------------
!!                                                                     !
!!    -----------------------------------------------------------------#
!!
!!
!!    --------------------------------------------------------------- &
!pub  call snlr  ( nrng,nang,npts,NZZ,                                &
!pub               frqa,pha, dens,dens2,                              &
!pub               kref2,kref4,jref2,jref4,wtk2,wtk4,                 &
!pub               wta2,wta4,tfac2,tfac4,grad,                        &
!pub               sumint,sumintsa,sumintp,sumintx,                   &
!pub               tsa,diag,  fbi,diag2 )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call snlr  ( nrng,nang,npts,NZZ, frqa,pha,                      &
                   sumint,sumintsa,sumintp,sumintx )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!
!!                 frqa,                                !* from /b/
!!                 pha, dens,                           !* from /c/
!!                 grad,sumint,                         !* from /dn/
!!                 kref2,kref4,jref2,jref4,wtk2,wtk4,   !* from /fr/
!!                 wta2,wta4,tfac2,tfac4,               !* from /fr/
!!                 dens2,sumintp,sumintx,sumintsa,      !* from /z/
!!                 tsa,diag,  fbi,diag2 )               !* Added by Bash
!!                 from  /b/, /c/, /dn/, /fr/, and /z/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!!    Pack results in proper format ---------------------------------- *
!!    S() & D() arrays are to be returned to WW3 in Cart. (k,theta) space
!!
      do 52 irng=1,nrng
        do 51 iang=1,nang
!!
!!        Bash; use tsa() & diag()
!!        Convert Polar tsa(k,theta) to Cartesian S(theta,k)
          S(iang,irng) = tsa(irng,iang) * wka(irng)   !* <=============
          D(iang,irng) = diag(irng,iang)
!!        --------------------------------
!!
!!        Bash; use fbi() & diag2()
!!        Convert Polar fbi(k,theta) to Cartesian S(theta,k)
!b        S(iang,irng) = fbi(irng,iang) * wka(irng)   !* <=============
!b        D(iang,irng) = diag2(irng,iang)
!!        --------------------------------
!!
  51    continue
  52  continue
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!
!!
!!    o o o o o o o o o
!!
!!    Bash; Add to print
!!    ------------------------------------------------------------------
!!
!!    After subroutine "snlr" have returned:
!!    Polar fbi(k,theta) and Polar tsa((k,theta) Snl in Action Density Units
!!    where fbi() & tsa() were calc. in "snlr" as follow:
!!    fbi(k,theta) = sumint(k,tetha)+sumintp(k,tetha)+sumintx(k,tetha)
!!    tsa(k,theta) = sumint(k,tetha)+sumintsa(k,tetha)
!!
!!    We calculate 8 diff flavors of fbi(k,theta) & fbi(f,theta)
!!    by multiplying it by 8 diff factors "fac??"
!!-1     1          * Polar     fbia1(k,theta) in Action Density Units
!!-2     k          * Cartesian fbia2(k,theta) in Action Density Units
!!-3     2pi/Cg     * Polar     fbia3(f,theta) in Action Density Units
!!-4     k*2pi/Cg   * Cartesian fbia4(f,theta) in Action Density Units
!!-5     w          * Polar     fbie1(k,theta) in Energy Density Units
!!-6     w*k        * Cartesian fbie2(k,theta) in Energy Density Units
!!-7     w*2pi/Cg   * Polar     fbie3(f,theta) in Energy Density Units
!!-8     w*k*2pi/Cg * Cartesian fbie4(f,theta) in Energy Density Units
!!    ------------------------------------------------------------------
!!
!!    and
!!    We calculate 8 diff flavors of tsa(k,theta) & tsa(f,theta)
!!    by multiplying it by 8 diff factors "fac??"
!!-1     1          * Polar     tsaa1(k,theta) in Action Density Units
!!-2     k          * Cartesian tsaa2(k,theta) in Action Density Units
!!-3     2pi/Cg     * Polar     tsaa3(f,theta) in Action Density Units
!!-4     k*2pi/Cg   * Cartesian tsaa4(f,theta) in Action Density Units
!!-5     w          * Polar     tsae1(k,theta) in Energy Density Units
!!-6     w*k        * Cartesian tsae2(k,theta) in Energy Density Units
!!-7     w*2pi/Cg   * Polar     tsae3(f,theta) in Energy Density Units
!!-8     w*k*2pi/Cg * Cartesian tsae4(f,theta) in Energy Density Units
!!    ------------------------------------------------------------------
!!
!!    for fbi dbl sums
!prt  fa1sum =  0.0         !* for dbl sum fbia1(k,theta)
!prt  fa2sum =  0.0         !* for dbl sum fbia2(k,theta)
!prt  fa3sum =  0.0         !* for dbl sum fbia3(f,theta)
!prt  fa4sum =  0.0         !* for dbl sum fbia4(f,theta)
!prt  fe1sum =  0.0         !* for dbl sum fbie1(k,theta)
!prt  fe2sum =  0.0         !* for dbl sum fbie2(k,theta)
!prt  fe3sum =  0.0         !* for dbl sum fbie3(f,theta)
!prt  fe4sum =  0.0         !* for dbl sum fbie4(f,theta)
!!    ------------------------------------------------------------------
!!
!!    for tsa dbl sums
!prt  ta1sum =  0.0         !* for dbl sum tsaa1(k,theta)
!prt  ta2sum =  0.0         !* for dbl sum tsaa2(k,theta)
!prt  ta3sum =  0.0         !* for dbl sum tsaa3(f,theta)
!prt  ta4sum =  0.0         !* for dbl sum tsaa4(f,theta)
!prt  te1sum =  0.0         !* for dbl sum tsae1(k,theta)
!prt  te2sum =  0.0         !* for dbl sum tsae2(k,theta)
!prt  te3sum =  0.0         !* for dbl sum tsae3(f,theta)
!prt  te4sum =  0.0         !* for dbl sum tsae4(f,theta)
!!    ------------------------------------------------------------------
!!
!!
!prt  do 62 irng=1,nrng
!!
!prt    afac1 = 1.0                                 !* afac1 = 1.0
!prt    afac2 = wka(irng)                           !* afac2 = k
!prt    afac3 = twopi/cga(irng)                     !* afac3 = 2pi/Cg
!prt    afac4 = wka(irng)*twopi/cga(irng)           !* afac4 = k*2pi/Cg
!prt    efac1 = oma(irng)*1.0                       !* efac1 = w*1.0
!prt    efac2 = oma(irng)*wka(irng)                 !* efac2 = w*k
!prt    efac3 = oma(irng)*twopi/cga(irng)           !* efac3 = w*2pi/Cg
!prt    efac4 = oma(irng)*wka(irng)*twopi/cga(irng) !* efac4 = w*k*2pi/Cg
!!
!prt    deltaf = delfx*frqa(irng)
!!      ----------------------------------------------------------------
!!
!!
!prt    do 61 iang=1,nang
!!
!!        calculate the 8 diff flavors of fbi(k,theta) & fbi(f,theta)
!prt      fbia1(irng,iang) = fbi(irng,iang) * afac1  !* fbia1(k,theta)
!prt      fbia2(irng,iang) = fbi(irng,iang) * afac2  !* fbia2(k,theta)
!prt      fbia3(irng,iang) = fbi(irng,iang) * afac3  !* fbia3(f,theta)
!prt      fbia4(irng,iang) = fbi(irng,iang) * afac4  !* fbia4(f,theta)
!prt      fbie1(irng,iang) = fbi(irng,iang) * efac1  !* fbie1(k,theta)
!prt      fbie2(irng,iang) = fbi(irng,iang) * efac2  !* fbie2(k,theta)
!prt      fbie3(irng,iang) = fbi(irng,iang) * efac3  !* fbie3(f,theta)
!prt      fbie4(irng,iang) = fbi(irng,iang) * efac4  !* fbie4(f,theta)
!!        --------------------------------------------------------------
!!
!!        calculate the 8 diff flavors of tsa(k,theta) & tsa(f,theta)
!prt      tsaa1(irng,iang) = tsa(irng,iang) * afac1  !* tsaa1(k,theta)
!prt      tsaa2(irng,iang) = tsa(irng,iang) * afac2  !* tsaa2(k,theta)
!prt      tsaa3(irng,iang) = tsa(irng,iang) * afac3  !* tsaa3(f,theta)
!prt      tsaa4(irng,iang) = tsa(irng,iang) * afac4  !* tsaa4(f,theta)
!prt      tsae1(irng,iang) = tsa(irng,iang) * efac1  !* tsae1(k,theta)
!prt      tsae2(irng,iang) = tsa(irng,iang) * efac2  !* tsae2(k,theta)
!prt      tsae3(irng,iang) = tsa(irng,iang) * efac3  !* tsae3(f,theta)
!prt      tsae4(irng,iang) = tsa(irng,iang) * efac4  !* tsae4(f,theta)
!!        --------------------------------------------------------------
!!
!!        calc. sums: dbl sum fbi??(k,theta) or dbl sum fbi??(f,theta)
!prt      fa1sum = fa1sum + fbia1(irng,iang)
!prt      fa2sum = fa2sum + fbia2(irng,iang)
!prt      fa3sum = fa3sum + fbia3(irng,iang)
!prt      fa4sum = fa4sum + fbia4(irng,iang)
!prt      fe1sum = fe1sum + fbie1(irng,iang)
!prt      fe2sum = fe2sum + fbie2(irng,iang)
!prt      fe3sum = fe3sum + fbie3(irng,iang)
!prt      fe4sum = fe4sum + fbie4(irng,iang)
!!        --------------------------------------------------------------
!!
!!        calc. sums: dbl sum tsa??(k,theta) or dbl sum tsa??(f,theta)
!prt      ta1sum = ta1sum + tsaa1(irng,iang)
!prt      ta2sum = ta2sum + tsaa2(irng,iang)
!prt      ta3sum = ta3sum + tsaa3(irng,iang)
!prt      ta4sum = ta4sum + tsaa4(irng,iang)
!prt      te1sum = te1sum + tsae1(irng,iang)
!prt      te2sum = te2sum + tsae2(irng,iang)
!prt      te3sum = te3sum + tsae3(irng,iang)
!prt      te4sum = te4sum + tsae4(irng,iang)
!!        --------------------------------------------------------------
!!
!p61    continue
!!
!!
!!      calc. dbl sums:  dbl sum fbi??(k,theta) * dk * dtetha
!!            or         dbl sum fbi??(f,theta) * df * dtetha
!prt    fa1sum = fa1sum * ainc*dwka(irng)  !* for dbl sum fbia1(k,theta)
!prt    fa2sum = fa2sum * ainc*dwka(irng)  !* for dbl sum fbia2(k,theta)
!prt    fa3sum = fa3sum * ainc*deltaf      !* for dbl sum fbia3(f,theta)
!prt    fa4sum = fa4sum * ainc*deltaf      !* for dbl sum fbia4(f,theta)
!prt    fe1sum = fe1sum * ainc*dwka(irng)  !* for dbl sum fbie1(k,theta)
!prt    fe2sum = fe2sum * ainc*dwka(irng)  !* for dbl sum fbie2(k,theta)
!prt    fe3sum = fe3sum * ainc*deltaf      !* for dbl sum fbie3(f,theta)
!prt    fe4sum = fe4sum * ainc*deltaf      !* for dbl sum fbie4(f,theta)
!!      ----------------------------------------------------------------
!!
!!      calc. dbl sums:  dbl sum tsa??(k,theta) * dk * dtetha
!!            or         dbl sum tsa??(f,theta) * df * dtetha
!prt    ta1sum = ta1sum * ainc*dwka(irng)  !* for dbl sum tsaa1(k,theta)
!prt    ta2sum = ta2sum * ainc*dwka(irng)  !* for dbl sum tsaa2(k,theta)
!prt    ta3sum = ta3sum * ainc*deltaf      !* for dbl sum tsaa3(f,theta)
!prt    ta4sum = ta4sum * ainc*deltaf      !* for dbl sum tsaa4(f,theta)
!prt    te1sum = te1sum * ainc*dwka(irng)  !* for dbl sum tsae1(k,theta)
!prt    te2sum = te2sum * ainc*dwka(irng)  !* for dbl sum tsae2(k,theta)
!prt    te3sum = te3sum * ainc*deltaf      !* for dbl sum tsae3(f,theta)
!prt    te4sum = te4sum * ainc*deltaf      !* for dbl sum tsae4(f,theta)
!!      ----------------------------------------------------------------
!!
!p62  continue
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!!    for fbi max
!prt  fa1max =  0.0         !* for     max fbia1(k,theta)
!prt  fa2max =  0.0         !* for     max fbia2(k,theta)
!prt  fa3max =  0.0         !* for     max fbia3(f,theta)
!prt  fa4max =  0.0         !* for     max fbia4(f,theta)
!prt  fe1max =  0.0         !* for     max fbie1(k,theta)
!prt  fe2max =  0.0         !* for     max fbie2(k,theta)
!prt  fe3max =  0.0         !* for     max fbie3(f,theta)
!prt  fe4max =  0.0         !* for     max fbie4(f,theta)
!!
!!    for tsa max
!prt  ta1max =  0.0         !* for     max tsaa1(k,theta)
!prt  ta2max =  0.0         !* for     max tsaa2(k,theta)
!prt  ta3max =  0.0         !* for     max tsaa3(f,theta)
!prt  ta4max =  0.0         !* for     max tsaa4(f,theta)
!prt  te1max =  0.0         !* for     max tsae1(k,theta)
!prt  te2max =  0.0         !* for     max tsae2(k,theta)
!prt  te3max =  0.0         !* for     max tsae3(f,theta)
!prt  te4max =  0.0         !* for     max tsae4(f,theta)
!!
!prt  do 75 irng=1,nrng
!b    do 75 irng=1,28       !* over first 28 freq. only
!!
!prt    do 73 iang=1,nang
!!
!!        calc. fbi??(?,theta) max:
!prt      if ( fbia1(irng,iang).gt.fa1max ) fa1max = fbia1(irng,iang)
!prt      if ( fbia2(irng,iang).gt.fa2max ) fa2max = fbia2(irng,iang)
!prt      if ( fbia3(irng,iang).gt.fa3max ) fa3max = fbia3(irng,iang)
!prt      if ( fbia4(irng,iang).gt.fa4max ) fa4max = fbia4(irng,iang)
!prt      if ( fbie1(irng,iang).gt.fe1max ) fe1max = fbie1(irng,iang)
!prt      if ( fbie2(irng,iang).gt.fe2max ) fe2max = fbie2(irng,iang)
!prt      if ( fbie3(irng,iang).gt.fe3max ) fe3max = fbie3(irng,iang)
!prt      if ( fbie4(irng,iang).gt.fe4max ) fe4max = fbie4(irng,iang)
!!
!!        calc. tsa??(?,theta) max:
!prt      if ( tsaa1(irng,iang).gt.ta1max ) ta1max = tsaa1(irng,iang)
!prt      if ( tsaa2(irng,iang).gt.ta2max ) ta2max = tsaa2(irng,iang)
!prt      if ( tsaa3(irng,iang).gt.ta3max ) ta3max = tsaa3(irng,iang)
!prt      if ( tsaa4(irng,iang).gt.ta4max ) ta4max = tsaa4(irng,iang)
!prt      if ( tsae1(irng,iang).gt.te1max ) te1max = tsae1(irng,iang)
!prt      if ( tsae2(irng,iang).gt.te2max ) te2max = tsae2(irng,iang)
!prt      if ( tsae3(irng,iang).gt.te3max ) te3max = tsae3(irng,iang)
!prt      if ( tsae4(irng,iang).gt.te4max ) te4max = tsae4(irng,iang)
!!
!p73    continue
!!
!p75  continue
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!prt  write(6,500)
!prt  write(6,501)
!prt  write(6,701) fa1max, fa1sum
!prt  write(6,702) fa2max, fa2sum
!prt  write(6,703) fa3max, fa3sum
!prt  write(6,704) fa4max, fa4sum
!prt  write(6,500)
!prt  write(6,705) fe1max, fe1sum
!prt  write(6,706) fe2max, fe2sum
!prt  write(6,707) fe3max, fe3sum
!prt  write(6,708) fe4max, fe4sum
!prt  write(6,511)
!!
!prt  write(6,500)
!prt  write(6,502)
!prt  write(6,701) ta1max, ta1sum
!prt  write(6,702) ta2max, ta2sum
!prt  write(6,703) ta3max, ta3sum
!prt  write(6,704) ta4max, ta4sum
!prt  write(6,500)
!prt  write(6,705) te1max, te1sum
!prt  write(6,706) te2max, te2sum
!prt  write(6,707) te3max, te3sum
!prt  write(6,708) te4max, te4sum
!prt  write(6,511)
!!
!500  format ('  ')
!501  format ('  8-dbl sums for FBI 2D Snl ')
!502  format ('  8-dbl sums for TSA 2D Snl ')
!511  format (' ======================================='              &
!             '========================================')
!!             1234567890123456789012345678901234567890
!!
!701  format (' max & dbl sum Polar     (k,theta) in Action =',2E14.6)
!702  format (' max & dbl sum Cartesian (k,theta) in Action =',2E14.6)
!703  format (' max & dbl sum Polar     (f,theta) in Action =',2E14.6)
!704  format (' max & dbl sum Cartesian (f,theta) in Action =',2E14.6)
!705  format (' max & dbl sum Polar     (k,theta) in Energy =',2E14.6)
!706  format (' max & dbl sum Cartesian (k,theta) in Energy =',2E14.6)
!707  format (' max & dbl sum Polar     (f,theta) in Energy =',2E14.6)
!708  format (' max & dbl sum Cartesian (f,theta) in Energy =',2E14.6)
!!    ------------------------------------------------------------------
!!
!!    Bash; End / Add to print
!!    ------------------------------------------------------------------
!!
!!    o o o o o o o o o
!!
!!    skip writing to Tape80, Tape81, Tape84 & Tape85  and
!!    skip writing to Tape11, Tape12, Tape40, Tape41 & Tape44
!!
!!
!!    Tape80 = '2d-energy.dat' is the 2D Energy density vs. (f,theta)
!!    ------------------------------------------------------------------
!a    write(80,'(''ZONE,T="Energy density",I=35,J=36,F=point'')')
!a    do 380 iang=1,nang
!a    do 380 ifr=1,nrng 
!a       write(80,179)          frqa(ifr),angl(iang)*deg,ef2(ifr,iang)
!179     format(3e11.3)
!380  continue
!!    ==================================================================
!!
!!
!!    write it for Matlab
!!    Tape11='freq.dat'
!a    write(11,111)
!111  format('% freq array (1-35) in Hz',/,'%  i    f(i)')
!111  format('% freq array (1-28) in Hz',/,'%  i    f(i)')
!a    do 311 ifr=1,nrng
!b    do 311 ifr=1,28
!a      write(11,112) ifr, frqa(ifr)
!112    format(i4,f10.5)
!311  continue
!!    ==================================================================
!!
!!
!!    write it for Matlab
!!    Tape12='direc.dat'
!a    write(12,121)
!121  format('% direc array (1-37) in Deg',/,'%  j    th(j)')
!a    do 312 iang=1,nang
!a      write(12,122) iang, angl(iang)*deg
!122    format(i4,f10.5)
!312  continue
!a    write(12,122) 37, 360.0
!!    ==================================================================
!!
!!
!!    write it for Matlab
!!    Tape40='2D_Energy.dat'
!!    ef2(f,theta)= Cartesian 2D Energy Density in (f,theta)
!a    write(40,*)'% Cartesian 2D Energy Density in (f,theta), I=35,J=37'
!b    write(40,*)'% Cartesian 2D Energy Density in (f,theta), I=28,J=37'
!a    write(40,*)'% 2D Enenrgy max = ', e2max
!a    do 340 iang=1,nang
!a      write(40,182) ( ef2(ifr,iang), ifr=1,nrng )
!b      write(40,182) ( ef2(ifr,iang), ifr=1,28 )
!182    format(35e14.6)
!340  continue
!a    write(40,182) ( ef2(ifr,1), ifr=1,nrng )
!b    write(40,182) ( ef2(ifr,1), ifr=1,28 )
!!    ==================================================================
!!
!!
!!    Tape81 = '2d-Snl-e.dat' is the 2D Snl (FBI) for Energy density
!!    ------------------------------------------------------------------
!a    write(81,'(''ZONE,T="Snl for energy",I=35,J=36,F=POINT'')')
!a    do 381 iang=1,nang
!a    do 381 ifr=1,nrng 
!a       write(81,179)                                             !* &
!a   &     frqa(ifr), angl(iang)*deg, twopi*frqa(ifr)*             !* &
!a   &     (sumint(ifr,iang) + sumintp(ifr,iang) + sumintx(ifr,iang))
!381  continue
!!    ==================================================================
!!
!!
!!    write it for Matlab
!!    Tape41='FBI_2D_Snl.dat'
!!    fbie1(k,theta)= Polar FBI 2D Snl in (k,theta) Energy Density Units
!a    write(41,*)'% Polar FBI 2D Snl in (k,theta) Energy Unit,I=35,J=37'
!b    write(41,*)'% Polar FBI 2D Snl in (k,theta) Energy Unit,I=28,J=37'
!a    write(41,*)'% FBI 2D Snl max = ', fe1max
!a    do 341 iang=1,nang
!a      write(41,182) ( fbie1(ifr,iang), ifr=1,nrng )
!b      write(41,182) ( fbie1(ifr,iang), ifr=1,28 )
!341  continue
!a    write(41,182) ( fbie1(ifr,1), ifr=1,nrng )
!b    write(41,182) ( fbie1(ifr,1), ifr=1,28 )
!!    ==================================================================
!!
!!
!!    Tape84 = '2d-TSA-e.dat' is the 2D Snl (TSA) for Energy density
!!    ------------------------------------------------------------------
!a    write(84,'(''ZONE,T="TSA for energy",I=35,J=36,F=POINT'')')
!a    do 384 iang=1,nang
!a    do 384 ifr=1,35  
!a       write(84,179) frqa(ifr),angl(iang)*deg, twopi*frqa(ifr)*  !* &
!a   &                (sumint(ifr,iang) + sumintsa(ifr,iang))
!384  continue
!!    ==================================================================
!!
!!
!!    write it for Matlab
!!    Tape44='TSA_2D_Snl.dat'
!!    tsae1(k,theta)= Polar TSA 2D Snl in (k,theta) Energy Density Units
!a    write(44,*)'% Polar TSA 2D Snl in (k,theta) Energy Unit,I=35,J=37'
!b    write(44,*)'% Polar TSA 2D Snl in (k,theta) Energy Unit,I=28,J=37'
!a    write(44,*)'% TSA 2D Snl max = ', te1max
!a    do 344 iang=1,nang
!a      write(44,182) ( tsae1(ifr,iang), ifr=1,nrng )
!b      write(44,182) ( tsae1(ifr,iang), ifr=1,28 )
!344  continue
!a    write(44,182) ( tsae1(ifr,1), ifr=1,nrng )
!b    write(44,182) ( tsae1(ifr,1), ifr=1,28 )
!!    ==================================================================
!!
!!
!!    Tape85='2d-energy-norm.dat' the 2D Energy density vs (f/f14,theta)
!!    (like tape80 ='2d-energy.dat' the 2D Energy density vs (f,theta))
!!    ------------------------------------------------------------------
!a    write(85,'(''ZONE,T="Energy density-norm",I=35,J=36,F=point'')')
!a    do 385 iang=1,nang
!a    do 385 ifr=1,nrng 
!a       write(85,179) frqa(ifr)/frqa(14),angl(iang)*deg,ef2(ifr,iang)
!385  continue
!!    ==================================================================
!!
!!    o o o o o o o o o
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!-----------------------------------------------------------------------------#
!!
!!
      RETURN
!!
!!
!!b   Bash; changed it to include internal function wkfnc
!!b         moved END SUBR. to the end and added CONTAINS + real fuction wkfnc
!!
!a    end                                         !* End of Main program sboltz
!!b   END SUBROUTINE W3SNLX
!!
      CONTAINS
!!
!!------------------------------------------------------------------------------
!!==============================================================================
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!
!!
!!b   Bash; Add real function wkfnc as an internal function
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    filename = "m0_wkfnc_7r3_1_1e"       in ~/msm/TSA/subrs7r3_99_m/ !
!!    from       "sbtsa-1-norm-Dec15-08.f-7r3"  in  same  dir.         !
!!    for        "w3snlxmd.ftn-NLX-7r3-TSA2-99-1-1e"                   !
!!    used for both FBI and TSA                                        !
!!                                              Bash Toulany           !
!!    -----------------------------------------------------------------#
!!
!!    --------------------------------------------------------------- &
      real function wkfnc ( f, dep )
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    Looks like a "Pade approximation" of an inversion of the linear  !
!!    wave dispersion relation. f = frequency (Hz), dep = depth (m),   !
!!                                                                     !
!!         sigma^2 = gk*tanh(kd), sigma = 2*pi*f,                      !
!!                                                                     !
!!         k = wavenumber (rad/m) = wkfnc.                             !
!!                                                                     !
!!    -----------------------------------------------------------------#
!!
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!!
      IMPLICIT NONE
!!
!!
      real, intent(in)  :: f, dep
!!
!!
!!    Local variables
!!
      real              :: twopi, g, y, x
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------::--------------------------------------------72
!!    ==================================================================
!!
!!
!!
!b    twopi = 6.283185308                       !* 2pi = 6.283,185,308 rad
!b    g     = 9.800
!!
      twopi = 6.283185400                       !* set = TPI  as in CONSTANTS
      g     = 9.806                             !* set = GRAV as in CONSTANTS
!!
      y     = ( (twopi*f)**2 ) * dep / g        !* sigma^2 d/g
!!
!!    --------------------------------------------------------------- &
      x     = y * ( y +                                               &
            1./(1.00000+y*(0.66667+y*(0.35550+y*(0.16084+y*(0.06320   &
            +y*(0.02174+y*(0.00654+y*(0.00171+y*(0.00039+y*0.00011)   &
            )))))))))
!!    --------------------------------------------------------------- &
!!
      x     = sqrt(x)                             !* kd
!!
      wkfnc = x / dep                             !* k
!!
      return
!!
!b    end                                         !* end function wkfnc
!!
      end function wkfnc
!!
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!b    end                                         !* End of Main program sboltz
      END SUBROUTINE W3SNLX
!!
!!------------------------------------------------------------------------------
!!==============================================================================
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    filename = "m1_gridsetr_7r3_1_1e"    in ~/msm/TSA/subrs7r3_99_m/ !
!!    from       "sbtsa-1-norm-Dec15-08.f-7r3"  in  same  dir.         !
!!    for        "w3snlxmd.ftn-NLX-7r3-TSA2-99-1-1e"                   !
!!    used for both FBI and TSA                                        !
!!                                              Bash Toulany           !
!!    -----------------------------------------------------------------#
!!
!!    --------------------------------------------------------------- &
!pub  SUBROUTINE gridsetr ( nrng,nang,npts,NZZ, dep,dfrq,ierr_gr,     &
!pub                  frqa,wka, cgnrng, sinan,cosan )        !* ,     &
!pub                  kref2,kref4,jref2,jref4,wtk2,wtk4,              &
!pub                  wta2,wta4,tfac2,tfac4,grad )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUBROUTINE gridsetr ( nrng,nang,npts,NZZ, dep,dfrq,ierr_gr,     &
                      frqa,wka, cgnrng, sinan,cosan )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!
!!                    wka,    frqa,      sinan,cosan,      !* from /b/
!!                    grad,           cgnrng,              !* from /dn/
!!                    kref2,kref4,jref2,jref4,wtk2,wtk4,   !* =    /fr/
!!                    wta2,wta4,tfac2,tfac4 )              !* =    /fr/
!!                    from  /b/, /dn/, and /fr/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!!-----------------------------------------------------------------------------#
!!                                                                             !
!!    This routine sets up the geometric part of the Boltzmann integral        !
!!    based on a grid of wave frequencies and directions, with wave-           !
!!    numbers related to frequency and depth by linear dispersion.  It         !
!!    is adapted from Don's original code with changes to modify the           !
!!    indexing so there are fewer unused elements, and a number of algo-       !
!!    rithmic changes that are mathematically equivalent to Don's but          !
!!    take advantage of intrinsic functions to form smooth results with        !
!!    less reliance on if statements.                                          !
!!    It calls locus-solving routines shloxr and shlocr and coupling           !
!!    coefficient routine cplshr.  If shlocr does not converge, ierr_gr        !
!!    will be something other than 0 and the routine will terminate,           !
!!    returning ierr_gr to the calling program.                                !
!!                                                                             !
!!    It returns array grad(,,), which is an estimate of the product           !
!!    C(k1,k2,k3,k4)*H(k1,k3,k4)*ds/|dW/dn| (where n and the k's are all       !
!!    vectors) as given, for example, by Eq.(7) of 'Nonlinear energy           !
!!    fluxes and the finite depth equilibrium range in wave spectra,'          !
!!    by Resio, Pihl, Tracy and Vincent (2001, JGR, 106(C4), p. 6985),         !
!!    as well as arrays for indexing, interpolating and weighting locus-       !
!!    based wavenumber vectors within the discrete solution grid.              !
!!-----------------------------------------------------------------------------!
!!                                                                             !
!!    The included file 'parfile_2.f' (below) contains just                    !
!!    parameters that sets the dimensions nrng, nang, npts and NZZ             !
!!    for some arrays (below)                                                  !
!!                                                                             !
!!    parameters in 'parfile_2.f' are:                                         !
!!    nrng  =   67 = maximum number of rings                                   !
!!    nang  =   36 = maximum number of angles                                  !
!!    npts  =   30 = maximum number of points around locus                     !
!!    NZZ   = 2278 =  NZZ = (nrng*(nrng+1))/2                                  !
!!                                                                             !
!b    include 'parfile_2.f'
!!                                                                             !
!!-----------------------------------------------------------------------------#
!!
!!------------------------------------------------------------------------------
!!--------------------------------------------------------------------72------80
!!==============================================================================
!!
!!
!!
      USE W3SERVMD, ONLY: EXTCDE
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
      IMPLICIT NONE
!!
!!
      integer, intent(in)  :: nrng, nang, npts, NZZ
      integer, intent(out) :: ierr_gr
!!
      real,    intent(in)  :: dep, dfrq
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      real,    intent(in)  :: wka(nrng),   frqa(nrng),                &
                              sinan(nang), cosan(nang)            !* /b/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      real,    intent(in)  :: cgnrng                             !* /dn/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!pub  integer, intent(out) :: kref2(npts,nang,NZZ),                   &
!pub                          kref4(npts,nang,NZZ),                   &
!pub                          jref2(npts,nang,NZZ),                   &
!pub                          jref4(npts,nang,NZZ)
!pub  real,    intent(out) ::  wtk2(npts,nang,NZZ),                   &
!pub                           wtk4(npts,nang,NZZ),                   &
!pub                           wta2(npts,nang,NZZ),                   &
!pub                           wta4(npts,nang,NZZ),                   &
!pub                          tfac2(npts,nang,NZZ),                   &
!pub                          tfac4(npts,nang,NZZ)               !* /fr/
!pub  real,    intent(out) ::  grad(npts,nang,NZZ)               !* /dn/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ---------------------::-----------------------------------------72
!!    ==================================================================
!!
!!
!!
!!    Local variables
!!
      integer           :: irng, krng, iang, kang, ipt, npmid  !*, na2p1
      integer           :: iizz, izz, ir, i
!!
      real              :: twopi, g, gsq, ainc
      real              :: f0,alf0,aldfrq, wk1x,wk1y, wk3x,wk3y
!!
      real              :: wn2,th2, wn2d,tnh2, om2,f2,cg2, er,tt2,w2
      real              :: wn4,th4, wn4d,tnh4, om4,f4,cg4,    tt4,w4
      real              :: dWdnsq, dWdn, dif13, dif14, Heaviside
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!pub
!!    Declare all shloxr/shlocr 5 returned arrays as PUBLIC arrays
!pub  real              :: wk2x(npts), wk2y(npts),                    &
!pub                       wk4x(npts), wk4y(npts), ds(npts)       !* /a/
!!
!pub
!!    Declare the only cplshr 1 returned variable as PUBLIC
!xx   can't declare as allocatable sibgle variable ?!
      real              :: csq
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------::--------------------------------------------72
!!    ##################################################################
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!
!!    initial constants
!!
!b    twopi  = 6.283185308                        !* 2pi = 6.283,185,308 rad
!b    g      = 9.800
!b    gsq    = 96.04
!!
      twopi  = 6.283185400                        !* set = TPI  as in CONSTANTS
      g      = 9.806                              !* set = GRAV as in CONSTANTS
      gsq    = 96.157636                          !* set = GRAV**2
!!
!!
      npmid  = npts/2 + 1                         !* mid-index of locus array
!b    na2p1  = nang/2 + 1                         !* the angle opposite to 1
      ainc   = twopi/nang                         !* angle increment (radians)
      f0     = frqa(1)                            !* lowest freq in grid
      alf0   = alog(frqa(1))                      !* ln(f0) for ir calc. below
      aldfrq = alog(dfrq)                         !* ln(dfrq)       "
!!
!!    initialize array grad
!!    Bash; add initialization of other returned arrays.
!!          Take it back, keep it as in original code
!!
      do 10 irng=1,nrng
         iizz = (nrng-1)*(irng-1)-((irng-2)*(irng-1))/2
         do 11 krng=irng,nrng
            izz = krng+iizz
            do 12 kang=1,nang
               do 13 ipt=1,npts
!!
                 grad(ipt,kang,izz)  = 0.0
!!
                 kref2(ipt,kang,izz) = 0
                 kref4(ipt,kang,izz) = 0
                 jref2(ipt,kang,izz) = 0
                 jref4(ipt,kang,izz) = 0
!!
                 wtk2(ipt,kang,izz)  = 0.0
                 wtk4(ipt,kang,izz)  = 0.0
                 wta2(ipt,kang,izz)  = 0.0
                 wta4(ipt,kang,izz)  = 0.0
!!
                 tfac2(ipt,kang,izz) = 0.0
                 tfac4(ipt,kang,izz) = 0.0
!!
  13           continue
  12        continue
  11     continue
  10  continue
!!------------------------------------------------------------------------------
!!
!!
!!wrt   Bash; Add test write output som1, som2, som3   in cplshr
!!      or test write output domsq23= denominator of   t1    in cplshr
!!      or test write output sumom  = denominator of csqhatd in cplshr
!!
!wrt    open(125, file='som123_at00.dat', status='unknown')   !* Tape125
!wrt    write(125,902) 'som1,som2,som3 from cplshr at irng,krng,izz,kang,ipt'
!!
!wrt    open(126, file='domsq23_at00.dat', status='unknown')  !* Tape125
!wrt    write(126,902) 'domsq23 from cplshr at irng,krng,izz,kang,ipt'
!!
!wrt    open(127, file='sumom_at00.dat', status='unknown')    !* Tape125
!wrt    write(127,902) 'sumom from cplshr at irng,krng,izz,kang,ipt'
!!
!902    format(1x,A)
!!------------------------------------------------------------------------------
!!
!!
!!    irng and iang are k1 parameters; krng and kang are k3 parameters
!!
      iang = 1                               !* set = 1 and will remain = 1
!!
!!20
      do 20 irng=1,nrng
        wk1x   = wka(irng)
        wk1y   = 0.0                         !* set = 0.0 and will remain = 0.0
        iizz = (nrng-1)*(irng-1)-((irng-2)*(irng-1))/2
!!
!!30
        do 30 krng=irng,nrng
!!
!!        Bash; check1
!!        Bash; change this ratio from > 4 to > 3  to make it
!!        consistent with similar test done in subr. snlr
          if ( frqa(krng)/frqa(irng) .gt. 4. ) go to 30   !* original
!b        if ( frqa(krng)/frqa(irng) .gt. 3. ) go to 30   !* Bash; use .gt. 3
!!
          izz = krng+iizz
!!
!!40
          do 40 kang=1,nang
            wk3x = wka(krng)*cosan(kang)
            wk3y = wka(krng)*sinan(kang)
!!
            if ( krng .eq. irng ) then            !* wn3 = wn1
!!---
!!            Bash; keep the opposite angles to k1 - original setting
              if ( kang .eq. iang ) go to 40      !* inf. if th3 = th1
!!                                                !* Note:   iang = 1
!!---
!!            Bash; eliminate the opposite angles to k1 - did not help
!b            if ( kang.eq.iang .or. kang.eq.na2p1 ) go to 40
!!                                                !* inf. if th3=+/-th1
!!                                                !* Note: iang=1,na2p1
!!---
!!            ------------------------------------------------------- &
!pub          call shloxr ( dep, wk1x,wk1y, wk3x,wk3y, npts,          &
!pub                             wk2x,wk2y, wk4x,wk4y, ds )
!!            -- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              call shloxr ( dep, wk1x,wk1y, wk3x,wk3y, npts )
!!            -- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!            ----------------------------------------------------------
!!
!!
            else                                  !* wn3 > wn1
!!
!!            ------------------------------------------------------- &
!pub          call shlocr ( dep, wk1x,wk1y, wk3x,wk3y, npts,          &
!pub                             wk2x,wk2y, wk4x,wk4y, ds, ierr_gr )
!!            -- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              call shlocr ( dep, wk1x,wk1y, wk3x,wk3y, npts,ierr_gr )
!!            -- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!            ----------------------------------------------------------
!!
              if ( ierr_gr .ne. 0 ) go to 99
!!
            end if
!!
!!          set the Heaviside coefficient
            dif13 = (wk1x-wk3x)**2+(wk1y-wk3y)**2       !* Note: wk1y = 0.0
!!
!!50
            do 50 ipt=1,npts
!!---
!!            Bash; keep the opposite angles to k1 - original setting
              if ( kang.eq.1 ) then                      !* if th3 = th1
                if (ipt.eq.1 .or. ipt.eq.npmid) go to 50 !* skip x-axis loci
              end if
!!---
!!            Bash; eliminate the opposite angles to k1 - did not help
!b            if ( kang.eq.1 .or. kang.eq.na2p1 ) then   !* if th3=+/- th1
!b              if (ipt.eq.1 .or. ipt.eq.npmid) go to 50 !* skip x-axis loci
!b            end if
!!---
!!
!!            set the coupling coefficient for ipt'th locus position
!!
!!            ----------------------------------------------------------
!!
!!wrt         Bash; Add test write output som1, som2, som3   in cplshr
!!            or test write output domsq23= denominator of   t1    in cplshr
!!            or test write output sumom  = denominator of csqhatd in cplshr
!!            -- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!            ------------------------------------------------------- &
!wrt
!wrt          call cplshr ( wk4x(ipt),wk4y(ipt), wk3x,wk3y,           &
!wrt                        wk2x(ipt),wk2y(ipt), dep,csq,             &
!wrt                        irng, krng, izz, kang, ipt )  !* for wrt
!!            -- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!---
              call cplshr ( wk4x(ipt),wk4y(ipt), wk3x,wk3y,           &
                            wk2x(ipt),wk2y(ipt), dep,csq)
!!            -- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!            ----------------------------------------------------------
!!
!!wn2         set parameters related to ipt'th locus wavenumber vector k2
!!            ----------------------------------------------------------
              wn2  = sqrt(wk2x(ipt)**2+wk2y(ipt)**2)  !* k2
              th2  = atan2(wk2y(ipt),wk2x(ipt))       !* k2 direction
              if ( th2 .lt. 0. ) th2 = th2 + twopi
              wn2d = wn2*dep                      !* k2*depth
              tnh2 = tanh(wn2d)                   !* tanh(k2*depth)
              om2  = sqrt(g*wn2*tnh2)             !* omega2 (rad)
              cg2  = 0.5*(om2/wn2)*(1.+wn2d*(1.-tnh2**2)/tnh2) !* group velocity
              f2   = om2/twopi                    !* f2 (Hz)
!!            ----------------------------------------------------------
!!
!!wn4         set parameters related to ipt'th locus wavenumber vector k4
!!            ----------------------------------------------------------
              wn4  = sqrt(wk4x(ipt)**2+wk4y(ipt)**2)
              th4  = atan2(wk4y(ipt),wk4x(ipt))
              if ( th4 .lt. 0. ) th4 = th4 + twopi
              wn4d = wn4*dep
              tnh4 = tanh(wn4d)
              om4  = sqrt(g*wn4*tnh4)
              cg4  = 0.5*(om4/wn4)*(1.+wn4d*(1.-tnh4**2)/tnh4)
              f4   = om4/twopi
!!            ----------------------------------------------------------
!!
!!
!!            ----------------------------------------------------------
!!
!!            set the Heaviside coefficient
              dif14 = (wk1x-wk4x(ipt))**2+(wk1y-wk4y(ipt))**2
!!                                                !* Note: wk1y = 0.0
              if ( dif13 .gt. dif14 ) then
                 Heaviside = 0.                   !* Eq(12) of RPTV
              else
                 Heaviside = 1.                   !* Eq(11) of RPTV
              end if
!!            ----------------------------------------------------------
!!
!!
!!            dWdn is the same as sqrt(zzsum) in Don's code, here reduced to a
!!            simpler but mathematically equivalent form that should vary
!!            smoothly between deep and intermediate water owing to identities
!!            using the computer's tanh() function
!!            ----------------------------------------------------------
!!
!!            set grad(,,);
!!            looks like the g^2 goes with csq (Webb'1978, eq. A2)
!!            ----------------------------------------------------------
!!
!sq
!b            dWdn = sqrt(cg2**2-2.*cg2*cg4*cos(th2-th4)+cg4**2)
!b            grad(ipt,kang,izz) = Heaviside*ds(ipt)*csq*gsq/dWdn
!sq
!!            Bash; check argument of dWdn sqrt
              dWdnsq = cg2**2 - 2.*cg2*cg4*cos(th2-th4) + cg4**2
!!
              if ( dWdnsq .le. 0.0 ) then
                print *, ' *** gridsetr; dWdn sqrt argument goes < 0.0 '
!b              print *, ' *** gridsetr; Bash set it to zero, didn't exit'
                dWdn               = 0.0
                grad(ipt,kang,izz) = 0.0
              else
                dWdn               = sqrt(dWdnsq)
                grad(ipt,kang,izz) = Heaviside*ds(ipt)*csq*gsq/dWdn
              endif
!!            --------------------------------------------------------72
!sq
!!
!!            ----------------------------------------------------------
!!
!!            set interpolation, indexing and weight parameters for
!!            computations along wavenumber radials
!!            ----------------------------------------------------------
!!
!!f2          --------------------
              if ( f2 .lt. f0 ) then
                 wtk2(ipt,kang,izz)  = 1.
                 tfac2(ipt,kang,izz) = 0.
                 kref2(ipt,kang,izz) = 1
              else
                 ir = 1+int((alog(f2)-alf0)/aldfrq)
!!
                 if ( ir+1 .gt. nrng ) then
                   wtk2(ipt,kang,izz)  = 0.
                   er = (wka(nrng)/wn2)**(2.5)
                   tt2= er*(cg2/cgnrng)*(frqa(nrng)/f2)*(wka(nrng)/wn2)
                   tfac2(ipt,kang,izz) = tt2
                   kref2(ipt,kang,izz) = nrng - 1
!!-
!!               Bash; add this check to make sure 'ir' doesn't go < 1
                 elseif ( ir .lt. 1 ) then
                   print *,' gridsetr: stop220-failed ir; ir < 1 ', ir
!b                 stop 220
                   CALL EXTCDE ( 220 )
!!-
                 else
                   w2 = (f2-frqa(ir))/(frqa(ir+1)-frqa(ir))
                   wtk2(ipt,kang,izz)  = 1. - w2
                   tfac2(ipt,kang,izz) = 1.
                   kref2(ipt,kang,izz) = ir
                 end if
              end if
!!            ----------------------------------------------------------
!!
!!f4          --------------------
              if ( f4 .lt. f0 ) then
                 wtk4(ipt,kang,izz)  = 1.
                 tfac4(ipt,kang,izz) = 0.
                 kref4(ipt,kang,izz) = 1
              else
                 ir = 1+int((alog(f4)-alf0)/aldfrq)
                 if ( ir+1 .gt. nrng ) then
                   wtk4(ipt,kang,izz)  = 0.
                   er = (wka(nrng)/wn4)**2.5
                   tt4= er*(cg4/cgnrng)*(frqa(nrng)/f4)*(wka(nrng)/wn4)
                   tfac4(ipt,kang,izz) = tt4
                   kref4(ipt,kang,izz) = nrng - 1
!!-
!!               Bash; add this check to make sure 'ir' doesn't go < 1
                 elseif ( ir .lt. 1 ) then
                   print *,' gridsetr: stop240-failed ir; ir < 1 ', ir
!b                 stop 240
                   CALL EXTCDE ( 240 )
!!-
                 else
                   w2 = (f4-frqa(ir))/(frqa(ir+1)-frqa(ir))
                   wtk4(ipt,kang,izz)  = 1. - w2
                   tfac4(ipt,kang,izz) = 1.
                   kref4(ipt,kang,izz) = ir
                 end if
              end if
!!            ----------------------------------------------------------
!!
!!
!!            set indexing and weight parameters for computations around
!!            azimuths; it appears that jref2 and jref4 should be bounded
!!            between 0 and nang-1 so that when iang (=1,nang) is added in
!!            the integration section, the proper bin index will arise;
!!            the weights wta2 and wta4 seem to be the fractional bin
!!            widths between th2 or th4 and the next increasing
!!            directional bin boundary
!!            ----------------------------------------------------------
!!
              i = int(th2/ainc)
              wta2(ipt,kang,izz)  = 1. - abs(th2-i*ainc)/ainc
              if ( i .ge. nang )  i = i - nang
              jref2(ipt,kang,izz) = i
!!
              i = int(th4/ainc)
              wta4(ipt,kang,izz)  = 1. - abs(th4-i*ainc)/ainc
              if ( i .ge. nang )  i = i - nang
              jref4(ipt,kang,izz) = i
!!
  50        continue                              !* end of ipt loop
!!
  40      continue                                !* end of kang loop
!!
  30    continue                                  !* end of krng loop
!!
  20  continue                                    !* end of irng loop
!!
  99  continue                                    !* ierr_gr bail
!!
!!wrt Bash; Add test write output som1, som2, som3   in cplshr
!!    or test write output domsq23= denominator of   t1    in cplshr
!!    or test write output sumom  = denominator of csqhatd in cplshr
!wrt  close(125)
!wrt  close(126)
!wrt  close(127)
!!
      return
!!
!b    end                                         !* end SUBROUTINE gridsetr
!!
      END SUBROUTINE gridsetr
!!
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    filename = "m2_shloxr_7r3_1_1e"      in ~/msm/TSA/subrs7r3_99_m/ !
!!    from       "sbtsa-1-norm-Dec15-08.f-7r3"  in  same  dir.         !
!!    for        "w3snlxmd.ftn-NLX-7r3-TSA2-99-1-1e"                   !
!!    used for both FBI and TSA                                        !
!!                                              Bash Toulany           !
!!    -----------------------------------------------------------------#
!!
!!    --------------------------------------------------------------- &
!pub  SUBROUTINE shloxr ( dep, wk1x,wk1y, wk3x,wk3y, npts,            &
!pub                           wk2x,wk2y, wk4x,wk4y, ds )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUBROUTINE shloxr ( dep, wk1x,wk1y, wk3x,wk3y, npts )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!!-----------------------------------------------------------------------------#
!!                                                                             !
!!    General locus solution for input vectors (wk1x,wk1y) and                 !
!!    (wk3x,wk3y) of the same magnitude but NOT in the same direction          !
!!    (or singularness will occur), output vectors (wk2x,wk2y) and             !
!!    (wk4x,wk4y), and element length ds along locus curve:                    !
!!                                                                             !
!!    With wavenumber vector n identified by (wknx,wkny), its magnitude        !
!!    given by wkn = sqrt(wknx**2+wkny**2) and its associated radian           !
!!    frequency given by sign = sqrt[g*wkn*tanh(wkn*dep)], where g is          !
!!    gravitational acceleration and dep is water depth, the four-wave         !
!!    resonance condition is satisfied along a locus of pts defined by         !
!!                                                                             !
!!    [1]  (wk1x,wk1y) + (wk2x,wk2y) - (wk3x,wk3y) - (wk4x,wk4y) = 0           !
!!                                                                             !
!!    [2]  sig1 + sig2 - sig3 - sig4 = 0                                       !
!!                                                                             !
!!    In the case where k1 [= sqrt(wk1x**2+wk1y**2)] is equal to k3            !
!!    [= sqrt(wk3x**2+wk3y**2)], we have by dispersion,                        !
!!                                                                             !
!!    [3]  sig1 = sqrt[g*k1*tanh(k1*h)] = sqrt[g*k3*tanh(k3*h)] = sig3         !
!!                                                                             !
!!    so sig1 - sig3 = 0 and [2] becomes sig2 = sig4, where, again by          !
!!    dispersion,                                                              !
!!                                                                             !
!!    [4]  sig2 = sqrt(g*k2*tanh(k2*h)] = sig4 = sqrt(g*k4*tanh(k4*h)]         !
!!                                                                             !
!!    and consequently k2 = k4.  This simplifies the locus solution            !
!!    considerably, and it can be shown that the (wk2x,wk2y) locus is          !
!!    along the perpendicular bisector of the (px,py) vector given by          !
!!                                                                             !
!!    [5]  (px,py) = (wk3x-wk1x,wk3y-wk1y)                                     !
!!                                                                             !
!!    and thereby from [1]                                                     !
!!    [6]  (wk4x,wk4y) = (wk2x,wk2y) - (px,py)                                 !
!!                                                                             !
!!    We note that these loci are independent of depth, although depth         !
!!    is used to set the length of the locus line by requiring that its        !
!!    range on either side of the p vector correspond to a wave with a         !
!!    freq four times that of the k1 vector (the locus line is made            !
!!    up of npts segments of length ds; the outer edges of the terminal        !
!!    segments satisfy the length constraint; vectors k2 and k4 extend         !
!!    to segment centers and will sufficiently approximate the length          !
!!    constraint).  As compared to srshlocr.f, we can do all                   !
!!    calculations here in dimensional space.                                  !
!!                                                                             !
!!-----------------------------------------------------------------------------|
!!                                                                             !
!!    The included file 'parfile_2.f' (below) contains just                    !
!!    parameters that sets the dimensions nrng, nang, npts and NZZ             !
!!    for some arrays (below)                                                  !
!!                                                                             !
!!    parameters in 'parfile_2.f' are:                                         !
!!    nrng  =   67 = maximum number of rings                                   !
!!    nang  =   36 = maximum number of angles                                  !
!!    npts  =   30 = maximum number of points around locus                     !
!!    NZZ   = 2278 =  NZZ = (nrng*(nrng+1))/2                                  !
!!                                                                             !
!b    include 'parfile_2.f'
!!                                                                             !
!!-----------------------------------------------------------------------------#
!!
!!------------------------------------------------------------------------------
!!--------------------------------------------------------------------72------80
!!==============================================================================
!!
!!
!/
!sq
!b    USE W3SERVMD, ONLY: EXTCDE
!/
      IMPLICIT NONE
!!
!!
      integer, intent(in)   :: npts
!!
      real,    intent(in)   :: dep
      real,    intent(in)   :: wk1x,wk1y, wk3x,wk3y
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!
!pub  real,    intent(out)  :: wk2x(npts), wk2y(npts),                &
!pub                           wk4x(npts), wk4y(npts), ds(npts)   !* /a/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ----------------------::----------------------------------------72
!!    ==================================================================
!!
!!
!!
!!    Local variables
!!
      integer           :: m, n, np2p1
!!
      real              :: twopi, g
      real              :: wk1, f1, wkx, fx, db, px, py, p,           &
                           thp, halfp, a, b, dth
!a    real              :: wkfnc                        !* real function
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------::--------------------------------------------72
!!    ##################################################################
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!
!!    initial constants
!b    twopi = 6.283185308                       !* 2pi = 6.283,185,308 rad
!b    g     = 9.800
!!
      twopi = 6.283185400                       !* set = TPI  as in CONSTANTS
      g     = 9.806                             !* set = GRAV as in CONSTANTS
!!
!!
!!    initial all returned arrays before they are computed
      do 5 n=1,npts
         wk2x(n) = 0.0
         wk2y(n) = 0.0
         wk4x(n) = 0.0
         wk4y(n) = 0.0
         ds(n)   = 0.0 
   5  continue
!!
      np2p1 = npts/2 + 1                       !* npts=30 ==> np2p1=16
!!
      wk1   = sqrt(wk1x**2+wk1y**2)            !* k1=wk1x since wk1y=0.0
      f1    = sqrt(g*wk1*tanh(wk1*dep))/twopi  !* f1=f(k1,dep)
!!
      fx    = 4. * f1                          !* fx = 4*f1
      wkx   = wkfnc(fx,dep)                    !* kx = k(4*f1,dep)
      db    = wkx/float(np2p1-1)               !* locus length increment
!!
!!
      px    = wk3x - wk1x
      py    = wk3y - wk1y                      !* Note: wk1y=0.0
!sq
!!    Bash; check argument of p sqrt
      p     = sqrt(px**2+py**2)
!sq
!b    if ( px**2+py**2 .le. 0.0 ) then
!b      print *, ' *** shloxr; sqrt argument goes < 0.0 "
!b      print *, ' *** shloxr; px**2+py**2 < 0.0; px, py =', px,py
!b      print *, ' *** shloxr; Bash called EXTCDE (111) to exit here'
!b      CALL EXTCDE ( 111 )
!b    else
!b      p   = sqrt(px**2+py**2)
!b    endif
!!    ----------------------------------------------------------------72
!sq
      thp   = atan2(py,px)
      halfp = 0.5*p
!!
!!
      do 10 n=np2p1,npts                       !* n = 16 --> 30
!!
         b   = 0.5 * db + float(n-np2p1) * db
         a   = sqrt(1.+(2.*b/p)**2)
         dth = acos(1./a)
!!
         wk2x(n) = a * halfp * cos(thp+dth)
         wk2y(n) = a * halfp * sin(thp+dth)
         wk4x(n) = wk2x(n) - px
         wk4y(n) = wk2y(n) - py
         ds(n)   = db
!!
         m       = npts - n + 1                !* m = 15 --> 1
         wk2x(m) = a * halfp * cos(thp-dth)
         wk2y(m) = a * halfp * sin(thp-dth)
         wk4x(m) = wk2x(m) - px
         wk4y(m) = wk2y(m) - py
         ds(m)   = db
!!
  10  continue
!!
      return
!!
!!
!!b   Bash; changed it to include internal function wlfnc
!!b         moved END SUBR. to the end and added CONTAINS + real fuction wkfnc
!!
!b    end                                         !* end SUBROUTINE shloxr
!b    END SUBROUTINE shloxr
!!
      CONTAINS
!!
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    filename = "m2_wkfnc_7r3_1_1e"       in ~/msm/TSA/subrs7r3_99_m/ !
!!    from       "sbtsa-1-norm-Dec15-08.f-7r3"  in  same  dir.         !
!!    for        "w3snlxmd.ftn-NLX-7r3-TSA2-99-1-1e"                   !
!!    used for both FBI and TSA                                        !
!!                                              Bash Toulany           !
!!    -----------------------------------------------------------------#
!!
!!    --------------------------------------------------------------- &
      real function wkfnc ( f, dep )
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    Looks like a "Pade approximation" of an inversion of the linear  !
!!    wave dispersion relation. f = frequency (Hz), dep = depth (m),   !
!!                                                                     !
!!         sigma^2 = gk*tanh(kd), sigma = 2*pi*f,                      !
!!                                                                     !
!!         k = wavenumber (rad/m) = wkfnc.                             !
!!                                                                     !
!!    -----------------------------------------------------------------#
!!
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!!
      IMPLICIT NONE
!!
!!
      real, intent(in)  :: f, dep
!!
!!
!!    Local variables
!!
      real              :: twopi, g, y, x
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------::--------------------------------------------72
!!    ==================================================================
!!
!!
!!
!b    twopi = 6.283185308                       !* 2pi = 6.283,185,308 rad
!b    g     = 9.800
!!
      twopi = 6.283185400                       !* set = TPI  as in CONSTANTS
      g     = 9.806                             !* set = GRAV as in CONSTANTS
!!
      y     = ( (twopi*f)**2 ) * dep / g        !* sigma^2 d/g
!!
!!    --------------------------------------------------------------- &
      x     = y * ( y +                                               &
            1./(1.00000+y*(0.66667+y*(0.35550+y*(0.16084+y*(0.06320   &
            +y*(0.02174+y*(0.00654+y*(0.00171+y*(0.00039+y*0.00011)   &
            )))))))))
!!    --------------------------------------------------------------- &
!!
      x     = sqrt(x)                             !* kd
!!
      wkfnc = x / dep                             !* k
!!
      return
!!
!b    end                                         !* end function wkfnc
!!
      end function wkfnc
!!
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!a    end                                         !* end SUBROUTINE shloxr
      END SUBROUTINE shloxr
!!
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    filename = "m3_shlocr_7r3_1_1e"      in ~/msm/TSA/subrs7r3_99_m/ !
!!    from       "sbtsa-1-norm-Dec15-08.f-7r3"  in  same  dir.         !
!!    for        "w3snlxmd.ftn-NLX-7r3-TSA2-99-1-1e"                   !
!!    used for both FBI and TSA                                        !
!!                                              Bash Toulany           !
!!    -----------------------------------------------------------------#
!!
!!    --------------------------------------------------------------- &
!pub  SUBROUTINE shlocr ( dep, wk1x,wk1y, wk3x,wk3y, npts,            &
!pub                           wk2x,wk2y, wk4x,wk4y, ds, ierr_gr )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUBROUTINE shlocr ( dep, wk1x,wk1y, wk3x,wk3y, npts,ierr_gr )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    With wavenumber vector n identified by (wknx,wkny), its magnitude!
!!    given by wkn = sqrt(wknx**2+wkny**2) and its associated radian   !
!!    frequency given by sign = sqrt[g*wkn*tanh(wkn*dep)], where g is  !
!!    gravitational acceleration and dep is water depth, the four-wave !
!!    resonance condition is satisfied along a locus of pts defined by !
!!                                                                     !
!!    [1]  (wk1x,wk1y) + (wk2x,wk2y) - (wk3x,wk3y) - (wk4x,wk4y) = 0   !
!!                                                                     !
!!    [2]  sig1 + sig2 - sig3 - sig4 = 0                               !
!!                                                                     !
!!    Because of the influence of depth, it is convenient to define new!
!!    vectors (wnx,wny) = (wknx*dep,wkny*dep) with magnitudes wn =     !
!!    sqrt(wnx**2+wny**2) = wkn*dep such that a dimensionless frequency!
!!    is sign*sqrt(dep/g) = sqrt[wkn*dep*tanh(wkn*dep)]                !
!!                        = sqrt[wn*tanh(wn)].                         !
!!    With these definitions and vectors (wk1x,wk1y) and (wk3x,wk3y)   !
!!    given as input, we can write (with some rearrangement            !
!!    of [1] and [2])                                                  !
!!                                                                     !
!!    [3]  w3x - w1x = px = w2x - w4x                                  !
!!    [4]  w3y - w1y = py = w2y - w4y                                  !
!!    [5]  sqrt[w3*tanh(w3)] - sqrt[w1*tanh(w1)] = q                   !
!!         = sqrt[w2*tanh(w2)] - sqrt[w4*tanh(w4)]                     !
!!                                                                     !
!!    With dimensionless vector (px,py) = (w3x-w1x,w3y-w1y) [magnitude !
!!    p = sqrt(px**2+py**2), direction atan2(py,px)] and dimensionless !
!!    frequency difference q = sqrt(w3*tanh(w3)] - sqrt(w1*tanh(w1)]   !
!!    defined by input parameters, we see from [3] and [4] that        !
!!    (w4x,w4y) = (w2x-px,w2y-py) [magnitude w4 = sqrt((w2x-px)**2 +   !
!!    (w2y-py)**2)] and thus from [5] we must basically find elements  !
!!    w2x and w2y that satisfy                                         !
!!                                                                     !
!!    [6] sqrt[sqrt(w2x**2+w2y**2)*tanh(sqrt(w2x**2+w2y**2))] -        !
!!        sqrt[sqrt((w2x-px)**2+(w2y-py)**2) *                         !
!!             tanh(sqrt((w2x-px)**2+(w2y-py)**2] = q                  !
!!                                                                     !
!!    The locus curve defined by the set of pts (w2x,w2y) crosses the  !
!!    p-vector axis at two points; one with magnitude w2=rmin*p with   !
!!    0.5 < rmin < 1.0 and one with magnitude w2=rmax*p with rmax > 1. !
!!    We first isolate rmin, rmax using various iterative algorithms,  !
!!    and then find locus pts that are not on the p-vector axis with   !
!!    another iterative scheme.  At the end, we un-normalize the w2    !
!!    and w4 vectors to find the wk2 and wk4 vectors.                  !
!!                                                                     !
!!    -----------------------------------------------------------------!
!!                                                                     !
!!    The included file 'parfile_2.f' (below) contains just            !
!!    parameters that sets the dimensions nrng, nang, npts and NZZ     !
!!    for some arrays (below)                                          !
!!                                                                     !
!!    parameters in 'parfile_2.f' are:                                 !
!!    nrng  =   67 = maximum number of rings                           !
!!    nang  =   36 = maximum number of angles                          !
!!    npts  =   30 = maximum number of points around locus             !
!!    NZZ   = 2278 =  NZZ = (nrng*(nrng+1))/2                          !
!!                                                                     !
!b    include 'parfile_2.f'
!!                                                                     !
!!    -----------------------------------------------------------------#
!!
!!------------------------------------------------------------------------------
!!--------------------------------------------------------------------72------80
!!==============================================================================
!!
!!
!/
!sq
!b    USE W3SERVMD, ONLY: EXTCDE
!/
      IMPLICIT NONE
!!
!!
      integer, intent(in)   :: npts
!!
      real,    intent(in)   :: dep
      real,    intent(in)   :: wk1x,wk1y, wk3x,wk3y
!!
      integer, intent(inout) :: ierr_gr
!!
!pub  real,    intent(out)  :: wk2x(npts), wk2y(npts),                &
!pub                           wk4x(npts), wk4y(npts), ds(npts)
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ----------------------::----------------------------------------72
!!    ==================================================================
!!
!!
!!
!!    Local variables
!!
      integer           :: n, np, nnp, np2p1, nplace
!!
      real              :: p,   px,   py,  q,    qrtp,  qsqp,         &
                           dth, thp,  dr,  dphi, cphi,                &
                           w1,  w1x, w1y,  wk1,  w3,  w3x, w3y,  wk3, &
                           rold,  rold1, rold2,  rnew,  rnew1, rnew2, &
                           pxod, pyod,  zpod,                         &
                           t, t1, t2, t3, tm, tp, ds1, ds2,           &
                           rmin, rmax, rcenter, rradius
!!
      double precision  :: dbt3, dbt4, dbt5, dbt6, dbz, dbp, dbqrtp,  &
                           cdthold, cdthnew, wate1, wate2
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------::--------------------------------------------72
!!    ##################################################################
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!
!!    initial all returned arrays before they are computed
      do 5 n=1,npts
        wk2x(n) = 0.0
        wk2y(n) = 0.0
        wk4x(n) = 0.0
        wk4y(n) = 0.0
        ds(n)   = 0.0 
   5  continue
!!
!!
      wk1   = sqrt(wk1x**2+wk1y**2)
      wk3   = sqrt(wk3x**2+wk3y**2)
!!
      w1    = wk1  * dep
      w1x   = wk1x * dep
      w1y   = wk1y * dep
      w3    = wk3  * dep
      w3x   = wk3x * dep
      w3y   = wk3y * dep
!!
      px    = w3x - w1x
      py    = w3y - w1y
!sq
!!    Bash; check argument of p sqrt
      p     = sqrt(px**2+py**2)
!sq
!b    if ( px**2+py**2 .le. 0.0 ) then
!b      print *, ' *** shlocr; sqrt argument goes < 0.0 "
!b      print *, ' *** shlocr; px**2+py**2 < 0.0; px, py =', px,py
!b      print *, ' *** shlocr; Bash called EXTCDE (112) to exit here'
!b      CALL EXTCDE ( 112 )
!b    else
!b      p   = sqrt(px**2+py**2)
!b    endif
!!    ----------------------------------------------------------------72
!sq
      thp  = atan2(py,px)
      q    = sqrt(w3*tanh(w3)) - sqrt(w1*tanh(w1))
      qrtp = q / sqrt(p)
      qsqp = qrtp**2
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    for (w2x,w2y) = rmin*(px,py) (locus crossing the p-vector axis   !
!!    nearest the origin), we have (w4x,w4y) = (w2x,w2y) - (px,py)     !
!!    = rmin*(px,py) - (px,py) = (rmin - 1)*(px,py); note that because !
!!    rmin < 1, the length of (w4x,w4y) is w4 = (1 - rmin)*p; then [6] !
!!    takes the simpler form                                           !
!!                                                                     !
!!    [7] sqrt(rmin*p*tanh(rmin*p)) - sqrt[(1-rmin)*p*tanh((1-rmin)*p)]!
!!        = q                                                          !
!!                                                                     !
!!    assuming the tanh() functions are slowly varying and can be      !
!!    treated as separate entities, [7] can be written as a quadratic  !
!!    in sqrt(rmin), i.e.,                                             !
!!                                                                     !
!!               2*qrtp*sqrt(rmin*p)                                   !
!!    [8] rmin - ------------------- sqrt(rmin) +                      !
!!                        T                                            !
!!                                                                     !
!!                                   qsqp-p*tanh((1-rmin)*p)           !
!!                                   -----------------------  =  0,    !
!!                                              T                      !
!!                                                                     !
!!    where   T    = tanh(rmin*p) + tanh((1-rmin)*p),                  !
!!            qrtp = q/sqrt(p)  and qsqp = qrtp**2                     !
!!                                                                     !
!!    the square of the most positive root of [8] can (with some       !
!!    algebra) be written                                              !
!!                                                                     !
!!    [9] rmin =                                                       !
!!   (1/T**2)*{qsqp*[tanh(rmin*p)-tanh((1-rmin)*p)]+T*tanh((1-rmin)*p) !
!!           +2*qrtp*sqrt[tanh(rmin*p)*tanh((1-rmin)*p)]*sqrt(T-qsqp)} !
!!                                                                     !
!!    setting rnew=rmin on the LHS and rold=rmin on the RHS in all     !
!!    instances allows the creation of an iterative algorithm for rmin;!
!!    convergence can be slow in general, so a coarse search for the   !
!!    crossing of [9] with the rnew=rold line is conducted first, then !
!!    a weighted iterative replacement loop is executed until the      !
!!    desired accuracy is achieved;                                    !
!!    Note that if p is sufficiently large, all tanh() -> 1            !
!!    and [9] becomes the analytic expression                          !
!!            rmin = 0.5*[1 + qrtp*sqrt(2-qsqp)]                       !
!!                                                                     !
!!    following is the coarse search using rold1 = 0.5,0.1,1.0,        !
!!                                         rold2 = rold1 + 0.1:        !
!!                                                                     !
!!    -----------------------------------------------------------------#
!!
!!
      ierr_gr = 0
!!
      rold1 = 0.5
      tp    = tanh(rold1 * p)
      tm    = tanh((1.-rold1) * p)
      t     = tp + tm
      t1    = qsqp * (tp-tm)
      t2    = t  * tm
      t3    = 2. * qrtp * sqrt(tp*tm) * sqrt(t-qsqp)
      rnew1 = (t1 + t2 + t3) / (t**2)
!!
!!
      do 10 n=1,4
        rold2 = rold1 + 0.1
        tp    = tanh(rold2 * p)
        tm    = tanh((1.-rold2) * p)
        t     = tp + tm
        t1    = qsqp * (tp-tm)
        t2    = t  * tm
        t3    = 2. * qrtp * sqrt(tp*tm) * sqrt(t-qsqp)
        rnew2 = (t1 + t2 + t3) / (t**2)
        if ( rnew2 .lt. rold2 ) then
          rold = (rold2*rnew1-rold1*rnew2)/(rold2-rold1-rnew2+rnew1)
          go to 11
        end if
        rold1 = rold2
        rnew1 = rnew2
  10  continue
      rold = 0.9                       !* default if not otherwise found
  11  continue
!!    ------------------------------------------------------------------
!!
!!
!!    iterative replacement search for rmin
      do 20 n=1,50
        tp   = tanh(rold * p)
        tm   = tanh((1.-rold) * p)
        t    = tp + tm
        t1   = qsqp * (tp-tm)
        t2   = t  * tm
        t3   = 2. * qrtp*sqrt(tp*tm)*sqrt(t-qsqp)
        rnew = (t1 + t2 + t3) / (t**2)
        if ( abs(rnew-rold) .lt. 0.00001 ) then
          rmin = rnew
          go to 21
        end if
        rold = 0.5 * (rold + rnew)
  20  continue
      ierr_gr = ierr_gr + 1  !* set 1's flag in ierr_gr if no convergence
      rmin = rnew
  21  continue
!!    ------------------------------------------------------------------
!!
!!
!!    set (dimensional) wavenumber components for this point on locus
!!
      wk2x(1) =  rmin * px / dep
      wk2y(1) =  rmin * py / dep
      wk4x(1) = (rmin-1.) * px / dep
      wk4y(1) = (rmin-1.) * py / dep
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    for (w2x,w2y) = rmax*(px,py) (locus crossing the p-vector axis   !
!!    farthest from the origin), we have (w4x,w4y)=(w2x,w2y) - (px,py) !
!!    = rmax*(px,py) - (px,py) = (rmax - 1)*(px,py);                   !
!!    here, because rmax > 1, the length of (w4x,w4y) is               !
!!    w4 = (rmax - 1)*p; then [6] takes the form                       !
!!                                                                     !
!!    [10] sqrt(rmax*p*tanh(rmax*p))-sqrt[(rmax-1)*p*tanh((rmax-1)*p)] !
!!         = q                                                         !
!!                                                                     !
!!    rearranging terms, squaring both sides and again rearranging     !
!!    terms yields                                                     !
!!                                                                     !
!!    [11] rmax*p*[tanh(rmax*p) - tanh((rmax-1)*p)]                    !
!!         = 2*q*sqrt(tanh(rmax*p))*sqrt(rmax*p) +                     !
!!           q**2 + p*tanh((rmax-1)*p)                                 !
!!                                                                     !
!!    because the difference of the two tanh()'s on the LHS tend to    !
!!    make the whole term small, we solve for rmax from the rapidly    !
!!    varying part of the first term on the RHS, i.e.,                 !
!!                                                                     !
!!                         [tanh((rmax-1)*p)+qsqp + rmax*T]**2         !
!!    [12]          rmax = ----------------------------------- ,       !
!!                                 4*qsqp*tanh(rmax*p)                 !
!!                                                                     !
!!    where, in this algorithm, T = tanh(rmax*p) - tanh((rmax-1)*p);   !
!!    as for rmin in [9], setting rnew=rmax on the LHS and rold=rmax   !
!!    in all instances on the RHS allows the formation of an iterative !
!!    algorithm; initially, we only know rmax > 1 so we do a coarse    !
!!    search in the 10's place out to some reasonably big number to    !
!!    try to find the place where [12] crosses the rnew = rold line    !
!!    (if this fails, we set an error flag); in refining the estimate, !
!!    it appears that [12] can get a little squirrely, so we do a      ! 
!!    brute force successive decimation search to nplace decimal places!
!!    to home in on the answer; note that if p is big enough for the   !
!!    tanh()'s to reach unity, [12] becomes exact and                  !
!!    rmax = [(1 + qsqp)**2]/(4*qsqp)                                  !
!!                                                                     !
!!    following is the coarse search with rold = 1,10,2001:            !
!!                                                                     !
!!    -----------------------------------------------------------------#
!!
!!
!!
      rold = 1.0
!!
      do 30 n=1,200
        rold = rold + 10.
        tp   = tanh(rold * p)
        tm   = tanh((rold-1.) * p)
        t    = tp - tm
        t1   = tm + qsqp
        t2   = 4. * tp * qsqp
        rnew = ((t1+rold*t)**2) / t2
        if ( rnew .lt. rold ) then
          rold = rold - 10.
          go to 31
        end if
  30  continue
      ierr_gr = ierr_gr + 10    !* set 10's place in ierr_gr if no sol'n
  31  continue
!!    ------------------------------------------------------------------
!!
!!
!!
!!    successive decimation search to refine rmax
      dr = 10.
      do 40 nplace=1,6
        dr = dr/10.
        do 50 n=1,10
          rold = rold + dr
          tp   = tanh(rold * p)
          tm   = tanh((rold-1.) * p)
          t    = tp - tm
          t1   = tm + qsqp
          t2   = 4. * tp * qsqp
          rnew = ((t1+rold*t)**2) / t2
          if ( rnew .lt. rold ) then
            rold = rold - dr
            go to 51
          end if
  50    continue
  51    continue
  40  continue
!!
      rmax = rold
!!
!!    set (dimensional) wavenumber components for this locus point
!!
      np2p1 = npts/2 + 1                       !* npts=30 ==> np2p1=16
!!
      wk2x(np2p1) =  rmax * px / dep
      wk2y(np2p1) =  rmax * py / dep
      wk4x(np2p1) = (rmax-1.) * px / dep
      wk4y(np2p1) = (rmax-1.) * py / dep
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    search for cos(dth) for off-p-vector solutions; use a circle     !
!!    centered on the p-vector axis at a distance                      ! 
!!    rcenter = 0.5*(rmax+rmin)  from the  origin with a               !
!!    radius  = 0.5*(rmax-rmin); radii from the center of the circle   !
!!    at successive angle increments np*dphi intersect the circle at   !
!!    distances r*p from the origin of the p vector such that          !
!!                                                                     !
!!    [13]  r**2 = rradius**2 + rcenter**2 -                           !
!!                 2*rcenter*rradius*cos(np*dphi)                      !
!!                                                                     !
!!    and makes an angle dth with the p vector satisfying              !
!!                                                                     !
!!    [14]   cdth = cos(dth) = (rcenter/r) - (rradius/r)*cos(np*dphi)  !
!!                                                                     !
!!    we then rotate this vector, holding its length=r*p constant and  !
!!    successively estimating cdth (using the above equation as an     !
!!    initial guess) until it intersects the locus curve; some algebra !
!!    yields the estimation equation as                                !
!!                                                                     !
!!                  r**2 + 1      [sqrt(r*tanh(rp)) - q/sqrt(p)]**4    !
!!    [15] cdthnew= ------- - ---------------------------------------- !
!!                    2*r     2*r*[tanh(p*sqrt(r**2-2*r*cdthold+1))]**2!
!!                                                                     !
!!    we use a weighted new estimate of cdthold with the weights based !
!!    on the argument of the tanh() function in the denominator        !
!!    (if the argument is big, tanh() -> 1 and cdthnew is found in one !
!!    pass; for small arguments, convergence is faster with equal      !
!!    weighting of old and new estimates; all this is empirical to try !
!!    to increase speed); double precision is used to gain enough      !
!!    accuracy when the arccos is taken; note that if p is big enough  !
!!    for all tanh()'s -> 1, [15] is exact and                         !
!!    cdthnew = cdth = [r**2 + 1 - (sqrt(r) - qrtp)**4]/(2*r)          !
!!                                                                     !
!!    -----------------------------------------------------------------#
!!
!!
      rcenter = 0.5 * (rmax + rmin)
      rradius = 0.5 * (rmax - rmin)
      t1      = rradius**2 + rcenter**2
      t2      = 2. * rradius * rcenter
      dphi    = 6.283185308 / float(npts)
      pxod    = px / dep
      pyod    = py / dep
!!
      dbp     = dble(p)
      dbqrtp  = dble(qrtp)
!!
!!
      do 60 np=2,npts/2                        !* np = 2 --> 15
!!
        cphi    = cos(float(np-1)*dphi)
        dbz     = dsqrt(dble(t1-t2*cphi))
        cdthold = dble(rcenter-rradius*cphi) / dbz
        dbt3    = (dbz**2) + 1.d0
        dbt4    = dbt3 / (2.d0*dbz)
        dbt5    = ((dsqrt(dbz*dtanh(dbz*dbp))-dbqrtp)**4)/(2.d0*dbz)
        dbt6    = dbp * dsqrt(dbt3-2.d0*dbz*cdthold)
!!
        if ( dbt6 .gt. 0.55d0 ) then
          wate1 = dtanh(dbt6)
          wate2 = 1.d0 - wate1
        else
          wate1 = 0.5d0
          wate2 = 0.5d0
        end if
!!
        do 70 n=1,25
          cdthnew = dbt4 - dbt5 / ((dtanh(dbt6))**2)
          if ( dabs(cdthnew-cdthold) .lt. 0.0000001d0 ) go to 71
          cdthold = wate1 * cdthnew + wate2 * cdthold
          dbt6    = dbp * dsqrt(dbt3-2.d0*dbz*cdthold)
  70    continue
        ierr_gr = ierr_gr + 100   !* add to 100's place for every failure
  71    continue
!!
        dth  = sngl(dacos(cdthnew))
        zpod = sngl(dbz) * p / dep
!!
        wk2x(np) = zpod * cos(thp+dth)
        wk2y(np) = zpod * sin(thp+dth)
        wk4x(np) = wk2x(np) - pxod
        wk4y(np) = wk2y(np) - pyod
!!
        nnp = npts-np+2                        !* nnp = 30 --> 17
!!
        wk2x(nnp) = zpod * cos(thp-dth)
        wk2y(nnp) = zpod * sin(thp-dth)
        wk4x(nnp) = wk2x(nnp) - pxod
        wk4y(nnp) = wk2y(nnp) - pyod
!!
  60  continue
!!
!!
!!    set arc length ds as the sum of half the segment lengths on either
!!    side of a given point
!!
      ds1   = sqrt((wk2x(2)-wk2x(1))**2+(wk2y(2)-wk2y(1))**2)
      ds(1) = ds1
      do 80 np=3,npts/2+1
        ds2   = sqrt((wk2x(np)-wk2x(np-1))**2+(wk2y(np)-wk2y(np-1))**2)
        ds(np-1)      = 0.5*(ds1+ds2)
        ds(npts-np+3) = ds(np-1)
        ds1           = ds2
  80  continue
      ds(npts/2+1)    = ds2
!!    ------------------------------------------------------------------
!!
      return
!!
!b    end                                      !* end SUBROUTINE shlocr
!!
      END SUBROUTINE shlocr
!!
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    filename = "m4_cplshr_7r3_1_1e"      in ~/msm/TSA/subrs7r3_99_m/ !
!!    from       "sbtsa-1-norm-Dec15-08.f-7r3"  in  same  dir.         !
!!    for        "w3snlxmd.ftn-NLX-7r3-TSA2-99-1-1e"                   !
!!    used for both FBI and TSA                                        !
!!                                              Bash Toulany           !
!!    -----------------------------------------------------------------#
!!
!!    ------------------------------------------------------------------
!!
!!wrt Bash; Add test write output som1, som2, som3   in cplshr
!!    or test write output domsq23= denominator of   t1    in cplshr
!!    or test write output sumom  = denominator of csqhatd in cplshr
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    --------------------------------------------------------------- &
!wrt
!wrt  SUBROUTINE cplshr ( w1x0,w1y0, w2x0,w2y0, w3x0,w3y0, h,csq,     &
!wrt                        irng, krng, izz, kang, ipt )
!---
      SUBROUTINE cplshr ( w1x0,w1y0, w2x0,w2y0, w3x0,w3y0, h,csq)
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    Calculates four-wave Boltzmann coupling coefficient in shallow   !
!!    water given k1,k2,k3 and following at least Hasselmann (1962)    !
!!    and probably Herterich and Hasselmann (1982).  Dimensional       !
!!    wavenumbers are (wnx0,wny0), n = 1,3, h = depth, csq = coupling  !
!!    coefficient.  This is the same as Don's cplesh, except within    !
!!    the algorithm, wavenumbers are made dimensionless with h and     !
!!    frequencies with sqrt(h/g), g = gravitational acceleration (the  !
!!    idea is to simplify and speed up the calculations while keeping  !
!!    a reasonable machine resolution of the result).  At the end,     !
!!    dimensionless csqhat is redimensioned as csq = csqhat/(h**6)     !
!!    so it is returned as a dimensional entity.                       !
!!                                                                     !
!!    This calculation can be a touchy bird, so we use double precision!
!!    for internal calculations, using single precision for input and  !
!!    output.                                                          !
!!                                                                     !
!!    -----------------------------------------------------------------#
!!
!!------------------------------------------------------------------------------
!!--------------------------------------------------------------------72------80
!!==============================================================================
!!
!!
!!
      IMPLICIT NONE
!!
!!
!!wrt Bash; Add test write output som1, som2, som3   in cplshr
!!    or test write output domsq23= denominator of   t1    in cplshr
!!    or test write output sumom  = denominator of csqhatd in cplshr
!wrt  integer, intent(in)   :: irng, krng, izz, kang, ipt
!!
      real,    intent(in)   :: w1x0,w1y0, w2x0,w2y0, w3x0,w3y0
      real,    intent(in)   :: h
      real,    intent(out)  :: csq
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ----------------------::----------------------------------------72
!!    ==================================================================
!!
!!
!!
!!    Local variables
!!
      integer           :: ipass
!!
      double precision  ::  pi4,   hh,  scple,                        &
                             t1,   t2,   t3,    t4,    t5,            &
                             s1,   s2,   s3,    k1,    k2,    k3,     &
                            k1x,  k2x,  k3x,   k1y,   k2y,   k3y,     &
                            om1,  om2,  om3, om1sq, om2sq, om3sq,     &
                           som1, som2, som3,   k23,  k23x,  k23y,     &
                           dot23,    dot123,  omsq23
!!
      double precision  :: k1x0, k2x0, k3x0,  k1zx,                   &
                           k1y0, k2y0, k3y0,  k1zy,                   &
                             p1,   p2,   p3,    p4,                   &
                             di,    e, csqd, csqhatd
!!
!!wrt Bash; test write output domsq23= denominator of   t1    in cplshr
      double precision  :: domsq23
!!
!!wrt Bash; test write output sumom  = denominator of csqhatd in cplshr
      double precision  :: sumom
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------::--------------------------------------------72
!!    ##################################################################
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!
!b    pi4   = 0.785398164d0               !* pi/4
      pi4   = 0.785398175d0               !* Set = PI/4 as in CONSTANTS
      hh    = dble(h)                     !* single to dbl precision
      scple = 0.d0                        !* initialize accumulator
      csq   = 0.d0                        !* initialize returned var.
!!
      do 10 ipass=1,3
!p1
        if (ipass .eq. 1) then            !* initial pass (+1,+1,-1)
          s1   =  1.d0
          s2   =  1.d0
          s3   = -1.d0
          k1x0 = dble(w1x0) * hh          !* norm. k elements with h
          k1y0 = dble(w1y0) * hh
          k2x0 = dble(w2x0) * hh
          k2y0 = dble(w2y0) * hh
          k3x0 = dble(w3x0) * hh
          k3y0 = dble(w3y0) * hh
!p1
!p2
        else if (ipass .eq. 2) then       !* 1st permutation (+1,-1,+1)
          s1   =  1.d0
          s2   = -1.d0
          s3   =  1.d0
          k1zx = k1x0
          k1zy = k1y0
          k1x0 = k2x0
          k1y0 = k2y0
          k2x0 = k3x0
          k2y0 = k3y0
          k3x0 = k1zx
          k3y0 = k1zy
!p2
!p3
        else                              !* 2nd permutation (-1,+1,+1)
          s1   = -1.d0
          s2   =  1.d0
          s3   =  1.d0
          k1zx = k1x0
          k1zy = k1y0
          k1x0 = k2x0
          k1y0 = k2y0
          k2x0 = k3x0
          k2y0 = k3y0
          k3x0 = k1zx
          k3y0 = k1zy
!p3
        end if
!!
!!
        k1x = s1 * k1x0                   !* sign the norm'ed k parts
        k1y = s1 * k1y0
        k2x = s2 * k2x0
        k2y = s2 * k2y0
        k3x = s3 * k3x0
        k3y = s3 * k3y0
!!
        k1  = dsqrt(k1x**2 + k1y**2)      !* normalized |k|
        k2  = dsqrt(k2x**2 + k2y**2)
        k3  = dsqrt(k3x**2 + k3y**2)
!!
        om1 = dsqrt(k1*dtanh(k1))         !* norm. omega (by sqrt(h/g))
        om2 = dsqrt(k2*dtanh(k2))
        om3 = dsqrt(k3*dtanh(k3))
!!
        om1sq = om1**2
        om2sq = om2**2
        om3sq = om3**2
!!
        som1 = s1 * om1                   !* sign the norm'ed omega's
        som2 = s2 * om2
        som3 = s3 * om3
!!
!!wrt   Bash; Add test write output som1, som2, som3   in cplshr
!wrt    write(125,909) irng,krng,izz,kang,ipt,                        &
!wrt                                 sngl(som1),sngl(som2),sngl(som3)
!909    format(2i3,i4,2i3,2x,3E14.6)
!!      == = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!!      ----------------------------------------------------------------
!!
        dot23  = k2x*k3x + k2y*k3y        !*  vector k2 dot vector k3
        k23x   = k2x + k3x                !* (vector k2  +  vector k3)_x
        k23y   = k2y + k3y                !* (vector k2  +  vector k3)_y
        k23    = dsqrt(k23x**2+k23y**2)   !* |vector k2  +  vector k3|
        omsq23 = k23 * dtanh(k23)         !* norm sq frq of v.k2+v.k3
        dot123 = k1x*k23x + k1y*k23y      !* v.k1 dot (v.k2 + v.k3)
!!      ----------------------------------------------------------------
!!
!!      note: the "i**2" factor from some reference is included in this term
!!
        di = -(som2+som3)*(om2sq*om3sq-dot23)+0.5d0 *                 &
              (som2*(k3**2-om3sq**2)+som3*(k2**2-om2sq**2))
!!
        e  = 0.5d0*(dot23-som2*som3*(om2sq+om3sq+som2*som3))
!!
        p1 = 2.d0 * (som1 + som2 + som3) * (om1sq*omsq23 - dot123)
        p2 = -som1 * (k23**2 - omsq23**2)
        p3 = -(som2 + som3) * (k1**2 - om1sq**2)
        p4 = k1**2 - om1sq**2
!!
!!
!!wrt   Bash; test write output domsq23= denominator of   t1    in cplshr
        domsq23 = omsq23 - ((som2+som3)**2)     !* Bash; needed for test
!!
!wrt    write(126,919) irng,krng,izz,kang,ipt, sngl(domsq23)
!919    format(2i3,i4,2i3,2x,E14.6)
!!
!!      Bash; discoverd that the denominator of t1 taking 0.0 values
!!            making t1 Not a Number "NaN" and the run will terminate.
!!      ----------------------------------------------------------------
!!---
!b2     t1 = di * (p1+p2+p3) / (omsq23 - ((som2+som3)**2))
!!---
!b      Bash; defined domsq23 = the denominator of t1 above and
!b            added a test for it to avoid dividing by 0.0
!!
        if ( domsq23 .eq. 0.d0 ) then             !* Bash; needed  test
!!        only domsq23 = 0.0 is bad               !* to avoid div. by 0.
!wrt      print *, ' domsq23=0. at irng,krng,izz, kang,ipt, omsq23,'  &
!wrt               ' som2,som3,-(som2+som3)**2=', irng,krng,izz,      &
!wrt                 kang,ipt, omsq23,som2,som3, -(som2+som3)**2
          t1 = 0.d0
        else
!!        domsq23 not 0.0 is OK
          t1 = di * (p1+p2+p3) / domsq23
        endif
!!---
!!      == = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!!      ----------------------------------------------------------------
!!
!!
        t2 = -di * som1 * (om1sq+omsq23)
        t3 = e * ((som1**3) * (som2+som3) - dot123 - p4)
!!
        t4 = 0.5d0 * som1 * dot23 *                                   &
          ((som1+som2+som3) * (om2sq+om3sq) + som2*som3*(som2+som3))
!!
        t5 = -0.5d0 * som1 *                                          &
            (om2sq * (k3**2) * (som1 + som2 + 2.d0 * som3)  +         &
             om3sq * (k2**2) * (som1 + som3 + 2.d0 * som2))
!!
        scple = scple + t1 + t2 + t3 + t4 + t5
!!
  10  continue
!!
!!
!!wrt Bash; test write output sumom  = denominator of csqhatd in cplshr
      sumom = om1*om2*om3*(om2+om3-om1)
!wrt  write(127,929) irng,krng,izz,kang,ipt, sngl(sumom)
!929  format(2i3,i4,2i3,2x,E14.6)
!!    ------------------------------------------------------------------
!!
!!
!!    Bash; also check the denominator of csqhatd from taking 0.0 values
!!          making csqhatd Not a Number "NaN" and the run will terminate.
!!      ----------------------------------------------------------------
!!---
!b2   csqhatd = scple*scple*pi4/(om1*om2*om3*(om2+om3-om1))  !* Bash; ok
!!---
!!wrt Bash; test write output sumom  = denominator of csqhatd in cplshr
      if ( sumom .eq. 0.d0 ) then     !* Bash; this test was Not needed
!!      only sumom = 0.0 is bad       !* div. by 0. causes NaN
!wrt    print *,' sumom=om1*om2*om3*(om2+om3-om1)=0.; om1,2,3,...= ', &
!wrt                            om1,om2,om3, (om2+om3-om1), sumom
        csqhatd = 0.d0                !* sumom was never 0.0
      else
        csqhatd = scple*scple*pi4 / sumom
      endif
!!---
!!    == = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!!    ------------------------------------------------------------------
!!
      csqd    = csqhatd / (hh**6)
!!
      csq     = sngl(csqd)                 !* from dbl to single precision
!b    csq     =      csqd
!!
      return
!!
!b    end                                         !* end SUBROUTINE cplshr
!!
      END SUBROUTINE cplshr
!!
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    filename = "m5_optsa2_7r3_1_1e"      in ~/msm/TSA/subrs7r3_99_m/ !
!!    from       "sbtsa-1-norm-Dec15-08.f-7r3"  in  same  dir.         !
!!    for        "w3snlxmd.ftn-NLX-7r3-TSA2-99-1-1e"                   !
!!    used for both FBI and TSA                                        !
!!                                              Bash Toulany           !
!!    -----------------------------------------------------------------#
!!
!!    ------------------------------------------------------------------
!!
!!    --------------------------------------------------------------- &
!!
!b    SUBROUTINE optsa2 ( nrng,nang,npk, f0,fpk,dfrq,dep, frqa,oma,   &
!b                        wka,cga, angl,cosan, ef1,ef2,  dens,dens2 )
!!
!!    --------------------------------------------------------------- &
!!
!pub  SUBROUTINE optsa2 ( nrng,nang,npk, fpk,dep, frqa,oma,           &
!pub                      wka,cga, angl,cosan, ef1,ef2,  dens,dens2 )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!
!!op2 Bash; new for optsa2, add # of 1st & last bins & total # of bins.
!!          follow !op2 for all th simple changes.
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!op2  SUBROUTINE optsa2 ( nrng,nang,npk, fpk,dep, frqa,               &
!op2                     oma,wka,cga, angl,cosan, ef1,ef2 )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      SUBROUTINE optsa2 (nrmn,nrmx,   nrng,nang, npk,fpk, frqa,      &
                         oma,wka,cga, angl,cosan, ef1,ef2,dfrq,nbins1)
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!
!!                 wka,oma,frqa, angl,   cosan,        !* from /b/
!!                 dens,                               !* from /c/
!!                 cga,                                !* from /dn/
!!                 ef2,                                !* from /enrgy/
!!                 dens2 )                             !* from /z/
!!                 from  /b/, /c/, /dn/, /enrgy/, and /z/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!!    ------------------------------------------------------------------
!!
!!    It returns variables dens(nrng,nang) and dens2(nrng,nang)
!!    ------------------------------------------------------------------
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    Using the input energy density ef2(nrng,nang)                    !
!!    calc. the polar (k,theta) action density                         !
!!    act2d(nrng,nang) = ef2*fac  where  fac = Cg/(2*pi*w*k)           !
!!                                                                     !
!!    Calculates and returns action density separated into two parts:  !
!!    (1) large-scale part into dens(nrng,nang)  and                   !
!!    (2) small-scale part into dens2(nrng,nang)                       !
!!    Both dens & dens2 are in Polar Action density (k,theta) space    !
!!                                                                     !
!!    NOTE: current version expects maxang to be located at angle=pi   !
!!                                                                     !
!!    -----------------------------------------------------------------#
!!
!!    3. Parameters :
!!
!!    Parameter list
!!    ------------------------------------------------------------------
!!op2   nrmn     Int.  I   number of first freq. bin  1 =< nrmn  < nrng
!!op2   nrmx     Int.  I   number of last  freq. bin  1  < nrmx =< nrng
!!op2   nbins1   Int.  I   actual # of bins > npk or npk2 (a min of ? bins is
!!                         guaranteed in equi. range, see nbins in subr. W3SNLX)
!!
!!      nrng     Int.  I   number of frequencies (sigma)
!!      nang     Int.  I   number of directions
!!      npk      Int.  I   number of peak frequency ( npk = 14 )
!!
!!      ef2      R.A.  I   Cartesian 2D Energy density ef2(f,theta)
!!                                             ------ dim=(nrng,nang)
!!
!!      ef1      R.A.  I   1D Energy density ef1(f) from ef2(f,theta)
!!                                             ------ dim=(nrng)
!!
!!      dfrq     Real  I   freq mult. for log freq spacing
!!      frqa     R.A.  I   radian frequencies (Hz) ------- dim=(nrng)
!!      fpk      Real  I   peak freq. [Hz] of initial freq spectrum
!!b     dep      Real  I   water depth [m]
!!
!!      angl     R.A.  I   dir. array (rad) (full circle); dim=(nrng)
!!      cosan    R.A.  I   cosine angles array ----------- dim=(nang)
!!
!!      oma      R.A.  I   rel. freq. array (rad*Hz) ----- dim=(nrng)
!!                         = twopi*frqa
!!      wka      R.A.  I   wavenumbers array [1/m] ------- dim=(nrng)
!!      cga      R.A.  I   group velocities array [m/s] -- dim=(nrng)
!!
!!      dens     R.A.  O   large-scale Action density (k,theta)
!!      dens2    R.A.  O   Small-scale Action density (k,theta)
!!                                               both dim=(nrng,nang)
!!    ------------------------------------------------------------------
!!    ----------------------------------------------------------------72
!!    ==================================================================
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    The included file 'parfile_2.f' (below) contains just            !
!!    parameters that sets the dimensions nrng, nang, npts and NZZ     !
!!    for some arrays (below)                                          !
!!                                                                     !
!!    parameters in 'parfile_2.f' are:                                 !
!!    nrng  =   67 = maximum number of rings                           !
!!    nang  =   36 = maximum number of angles                          !
!!    npts  =   30 = maximum number of points around locus             !
!!    NZZ   = 2278 =  NZZ = (nrng*(nrng+1))/2                          !
!!                                                                     !
!b    include 'parfile_2.f'
!!                                                                     !
!!    -----------------------------------------------------------------#
!!
!!------------------------------------------------------------------------------
!!--------------------------------------------------------------------72------80
!!==============================================================================
!!
!!
!!
      IMPLICIT NONE
!!
!!
!!op2 Bash; new for optsa2
      integer, intent(in)  :: nrmn, nrmx, nbins1
!!
      integer, intent(in)  :: nrng, nang, npk
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!b    real,    intent(in)  :: f0, fpk, dfrq, dep
      real,    intent(in)  ::     fpk, dfrq
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      real,    intent(in)  :: ef2(nrng,nang)                  !* /enrgy/
      real,    intent(in)  :: ef1(nrng)
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      real,    intent(in)  :: frqa(nrng), oma(nrng), wka(nrng),       &
                              angl(nang), cosan(nang)            !* /b/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      real,    intent(in)  :: cga(nrng)                          !* /dn/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!
!pub  real,    intent(out) :: dens(nrng,nang)                    !* /c/
!pub  real,    intent(out) :: dens2(nrng,nang)                   !* /z/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ---------------------::-----------------------------------------72
!!    ==================================================================
!!
!!
!!
!!    Local variables
!!
      integer              :: irng, iang
      integer              :: n1, n2, nn1, nn2, m, mm, ii
      integer              :: neq,  maxang,   idif, igam
!!
      real                 :: twopi, ainc, sum, fac
      real                 :: fovfp, beta, gam, emax, y, qmin
      real                 :: adif,  fdenp, fr, ratio, z, ddd
!b    real                 :: sigz
!!
      real                 :: fk(nrng),    fknrm(nrng)
      real                 :: bscl1(nrng), fkscl1(nrng)
!b    real                 :: act1d(nrng)
      real                 :: psi2(nang)
!!
      real                 :: q(16)
      real                 :: act2d(nrng,nang)
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ---------------------::-----------------------------------------72
!!    ##################################################################
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!
!!    NOTE:  current version expects maxang to be located at angle=PI
!!
!b    twopi= 8. * atan(1.0)                       !* 2pi rad
      twopi = 6.283185400                         !* set = TPI  as in CONSTANTS
!!
!!    Bash; use the passed in 'npk' parameter instead of 14
!b    f0   = fpk * (dfrq**(-13))                  !* 13=14-1 is magic
!x    f0   = fpk * (dfrq**(-(npk-1)))             !* npk=14  is magic
!!                                                !* f0    now  input
!!
!b    ainc is input, but keep this one for print accuracy
      ainc = twopi/nang                           !* angle increment (radians)
!!
!!
!!    solve for normalizing coefficient for INT (cos**m) **(-1)
!b    n1 = -nang/4 + 1
!b    n2 =  nang/4 + 1
!!
!!    Bash redefined n1 & n2 
!!    the prev. n1, n2 spans half circle (from 270. to 90. going through 0.)
!!    with n1 always falls on 270. & n2 always falls on 90. (after correction)
!!    and with cosine 270 & 90 are 0.0, there is no contribution to "sum"
      n1 = -nang/4 + 2         !* Bash redefined n1,  was  n1 = -nang/4 + 1
      n2 =  nang/4             !* Bash redefined n2,  was  n2 =  nang/4 + 1
!!
!!    and here the values of q(16) almost the same for any nang
!!    m, q(m) =  1    0.500317
!!    m, q(m) =  2    0.636620
!!    m, q(m) =  3    0.749999
!!    m, q(m) =  4    0.848826
!!    m, q(m) =  5    0.937500
!!    m, q(m) =  6    1.018591
!!    m, q(m) =  7    1.093750
!!    m, q(m) =  8    1.164104
!!    m, q(m) =  9    1.230468
!!    m, q(m) = 10    1.293449
!!    m, q(m) = 11    1.353515
!!    m, q(m) = 12    1.411035
!!    m, q(m) = 13    1.466308
!!    m, q(m) = 14    1.519576
!!    m, q(m) = 15    1.571044
!!    m, q(m) = 16    1.620882
!!
!prt  print *, '  '
!prt  print *, '  *** optsa2 ***  '
      do 15 m=1,16
        sum = 0.
        do 16 iang=n1,n2
          ii = iang
          if ( iang .lt. 1 ) ii = iang + nang
          sum = sum + cosan(ii)**m
  16    continue
        q(m) = 1./(sum*ainc)
!b      write(6,105) m, q(m)
!105    format(2x,'  m, q(m) =',I3,F12.6)
  15  continue
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!
!prt  print *, '  '
!op2  do 24 irng=1,nrng
      do 24 irng=nrmn,nrmx
        fac   = cga(irng)/(twopi*oma(irng)*wka(irng))
        do 25 iang=1,nang
!!        Convert 2d Energy Density ef2(f,theta)   Cartesian
!!             to 2d Action Density act2d(k,theta) Polar.
          act2d(irng,iang) = ef2(irng,iang) * fac
  25    continue
!!
!!      Convert ef1(f) to fk(k); both are 1d Energy Density Cartesian
        fk(irng)    = cga(irng)*ef1(irng)/twopi  !* fk(k)    energy
!!
!!      Normalize the 1d cartesian Energy Density fk(k) to give fknrm(k)
        fknrm(irng) = fk(irng)*wka(irng)**2.5    !* fknrm(k)=norm. fk(k)
!!
!!      Convert 1d Energy Density fk(k) to 1d Action Density act1d(k);
!!      both are Cartesian
!x      act1d(irng) = fk(irng)/oma(irng)         !* act1d(k) action
!!
!!      don't print ef2(f,19)*ainc/ef1(f) at all f,
!!                            as one ef1(f) may be 0.0
!b      print 601, ' f,e(f), fknrm(f), ef2(f,19)*ainc/ef1(f) = ',     &
!b      frqa(irng),ef1(irng),fknrm(irng),ef2(irng,19)*ainc/ef1(irng)
!601    format(4E14.6)
  24  continue
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!
!!    fit parameters to spectrum
!!    --------------------------
      sum = 0.
      neq = 0
!op2  do 26 irng=1,nrng
      do 26 irng=nrmn,nrmx
        fovfp = frqa(irng)/fpk
!!      Bash; check2 test equilibrium range
!b      if ( fovfp.ge.1.55.and.fovfp.le.2.45 ) then !* original   equi range
!b      if ( fovfp.ge.1.20.and.fovfp.le.2.20 ) then !* Bash; wide equi range
!!      Bash; select the narrow equi. range in line  with TSA min condition
!b      if ( fovfp.ge.1.90.and.fovfp.le.2.20 ) then !* narrow equi range
!!op2
!!      Bash; select new narrow equi. range in line  with TSA min condition
!!            that takes into account the low freq peak for double peaks
!b      if ( fovfp.ge.1.20.and.fovfp.le.1.40 ) then !* narrow equi range
!b      if ( fovfp.ge.1.13.and.fovfp.le.1.18 ) then !* narrow equi range  <<<<<
!!      --------------------------------------------------------------!*  <<<<<
!!
        if ( fovfp.ge.(dfrq**(nbins1))-0.025 .and.                    &
             fovfp.le.(dfrq**(nbins1))+0.025 ) then  !* narrow equi range <<<<<
          sum = sum + fknrm(irng)
          neq = neq + 1
        endif
!!
  26  continue
      beta = sum / neq
!!    scale beta down to adjust the equi range to a propre range assuming f**-4
!b    beta = 0.22 * beta
!!
!!    Bash; use the passed in 'npk' parameter instead of 14
!b    gam  = fknrm(14) / beta              !*  14      is magic
      gam  = fknrm(npk) / beta             !*  npk=14  is magic
!b+   print *, ' estimated beta, gamma     = ', beta, gam
!!
!op2  do 226 irng=1,nrng
      do 226 irng=nrmn,nrmx
         fknrm(irng) = fknrm(irng) / beta
 226  continue
!!
!!
!!   compare value of peak with q-array for closest fit to cos**m at peak
!!   find peak direction at spectral peak
!!
      emax   = 0.
      maxang = 0
      do 27 iang=1,nang
!!
!!       Bash; use the passed in 'npk' parameter instead of 14
!b       if ( ef2(14,iang).gt.emax ) then  !*  14      is magic
         if ( ef2(npk,iang).gt.emax ) then !*  npk=14  is magic
!b           emax   = ef2(14,iang)         !*  14      is magic
             emax   = ef2(npk,iang)        !*  npk=14  is magic
             maxang = iang                 !* 1 <= maxang <= nang
         endif
  27  continue
!!
!x    print *, ' maxang = ', maxang
!!
!!
!!    Bash; use the passed in 'npk' parameter instead of 14
!b    y  = ef2(14,maxang) / ef1(14)        !*  14      is magic
      y  = ef2(npk,maxang) / ef1(npk)      !*  npk=14  is magic
!x    print *, ' y = ', y
!!
      mm   = 1
      qmin = abs(q(1)-y)
      do 28 m=2,16
         adif = abs(q(m)-y)
         if ( adif.lt.qmin ) then
            qmin = adif
            mm   = m
         endif
  28  continue
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!
!b+   print *, ' estimated spreading power = ', mm
!prt  print *, ' estimated beta, gamma, sprd pwr = ', beta,gam,mm

!!    construct directional spreading function
!!    ----------------------------------------
!!
!!    Bash; initialize psi2() array here
      do 128 iang=1,nang
         psi2(iang) = 0.0
 128  continue
!!
!!    Bash; check3 - to make sure nn1, nn2 don't go out of range
      nn1 = maxang - nang/4
      nn2 = maxang + nang/4
!!
      do 29 iang=nn1,nn2
!!
!!       Bash; add this check on iang so it won't get out of range.
!!       reset iang into ii so it doesn't take value outside its range
!!       warning: iang can take values outside the 1 to nang range
!!       iang can take -ve values 1-1=0 to 1-6=-5 or nang+1 to nang+6
!!       use ii rotated iang that reamins in range
!!       ---------------------------------------------------------------
         ii = iang
         if ( ii .lt. 1 ) then
            ii = ii + nang
         elseif ( ii .gt. nang ) then
            ii = ii - nang
         endif
!!       ---------------------------------------------------------------
!!
!!       Bash; check4 - check idif
         idif = iabs(maxang-iang) + 1   !* Bash; use  iang (original)
!b       idif = iabs(maxang-ii)   + 1   !* Bash; don't use ii  xx
!!       ---------------------------------------------------------------
!!
         if ( idif.gt.nang/4 ) then
!b          psi2(iang) = 0.0            !* Bash; iang sometimes goes < 0
            psi2(ii)   = 0.0            !* Bash; here use ii
         else                           !* since 1 < ii < 36
!!                                                 =    =
!b          psi2(iang) = q(mm) * cos(angl(idif))**mm    !* ditto (above)
            psi2(ii)   = q(mm) * cos(angl(idif))**mm
!!                                      !* Bash; here use ii
         endif                          !* since 1 < ii < 36
!!                                                 =    =
!!
  29  continue
!!    ------------------------------------------------------------------
!!
!!
!!    estimate parametric spectrum and deviation from parametric spectrum
!!    ------------------------------------------------------------------
!!
      igam = (gam-0.4)*10 + 0.5
!b-   print *, ' gamma, igam     = ', gam, igam    !* gam is redundent
!!
!b    sigz  = 0.109
      gam   = igam/10. + 0.4
!!
!!    Bash; use the passed in 'npk' parameter instead of 14
!b    fdenp = gam * beta / wka(14)**2.5       !*  14      is magic
      fdenp = gam * beta / wka(npk)**2.5      !*  npk=14  is magic
!b-   print *, ' estimated gam from igam   = ', gam
!prt  print *, ' igam, estimated gam from igam   = ', igam,gam
!!
!op2  do 40 irng=1,nrng
      do 40 irng=nrmn,nrmx
!x       print *, 'inside 40 loop, f =  ', frqa(irng)
         fr = frqa(irng) / fpk
!x       print *, 'f,fpk,fr = ', frqa(irng),fpk,fr
!!
         if ( fr.le.1.0001 ) then
            if ( fr.ge.0.85 ) then
               ratio = 1.-(1.-fr)*0.7/0.15
            else
              ratio = 0.3*exp(-17.3*(0.85-fr))
            endif
!prt        print *, ' fr, ratio   = ', fr, ratio
            fkscl1(irng) = fdenp*ratio
            bscl1(irng)  = fkscl1(irng)/oma(irng)
         else
!b          z = 0.5*((fr-1.)/sigz)**1.2
            z = 0.5*((fr-1.)/0.109)**1.2          !* sigz = 0.109
            if ( z.gt.6. ) z = 6.
            ratio = 1.+exp(-z)*(gam-1.)
            fkscl1(irng) = beta*ratio/wka(irng)**2.5
            bscl1(irng)  = fkscl1(irng)/oma(irng)
         endif
!!
!x       print *, 'fpk,fr,fknrm = ',                                  &
!x                 fpk,fr,fkscl1(irng)*wka(irng)**2.5/beta
!!
!prt     print *, ' fk, fkscl   = ', fk(irng),fkscl1(irng)
!!
         do 41 iang=1,nang
            ddd = bscl1(irng) * psi2(iang) / wka(irng)    !* large-scale
            dens(irng,iang)  = ddd                        !* large-scale
            dens2(irng,iang) = act2d(irng,iang) - ddd     !* small-scale
  41     continue
!!
!x       print *, irng
!x       print 4949,  (dens(irng,iang), iang=1,nang)
!x       print 4949, (dens2(irng,iang), iang=1,nang)
!x       print 4949, (act2d(irng,iang), iang=1,nang)
!4949     format(36e10.3)

  40  continue
!!    ==================================================================
!!
!!
      return
!!
!!
!b    end                                         !* end SUBROUTINE optsa2
!!
      END SUBROUTINE optsa2
!!
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    filename = "m6_snlr_7r3_1_1e"        in ~/msm/TSA/subrs7r3_99_m/ !
!!    from       "sbtsa-1-norm-Dec15-08.f-7r3"  in  same  dir.         !
!!    for        "w3snlxmd.ftn-NLX-7r3-TSA2-99-1-1e"                   !
!!    used for both FBI and TSA                                        !
!!                                              Bash Toulany           !
!!    -----------------------------------------------------------------#
!!
!!    --------------------------------------------------------------- &
!pub  SUBROUTINE snlr ( nrng,nang,npts,NZZ,                           &
!pub               frqa,pha, dens,dens2,                              &
!pub               kref2,kref4,jref2,jref4,wtk2,wtk4,                 &
!pub               wta2,wta4,tfac2,tfac4,grad,                        &
!pub               sumint,sumintsa,sumintp,sumintx,                   &
!pub               tsa,diag,  fbi,diag2 )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUBROUTINE snlr ( nrng,nang,npts,NZZ, frqa,pha,                 &
                   sumint,sumintsa,sumintp,sumintx )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!
!!                 frqa,                                !* from /b/
!!                 pha, dens,                           !* from /c/
!!                 grad,sumint,                         !* from /dn/
!!                 kref2,kref4,jref2,jref4,wtk2,wtk4,   !* from /fr/
!!                 wta2,wta4,tfac2,tfac4,               !* from /fr/
!!                 dens2,sumintp,sumintx,sumintsa,      !* from /z/
!!                 tsa,diag,  fbi,diag2 )               !* Added by Bash
!!                 from  /b/, /c/, /dn/, /fr/, and /z/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!    it returns: sumint, sumintsa, sumintp, sumintx and
!!                tsa, diag, fbi, diag2  all with dim=(nrng,nang)
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    For a given action density array dens(k,theta), computes the     !
!!    rate-of-change array sumint(k,theta) = dN(k,theta)/dt owing to   !
!!    wave-wave interaction, as well as some ancillary arrays          !
!!    relating to positive and negative fluxes and their integrals.    !
!!                                                                     !
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    Compute:                                                         !
!!    --------                                                         !
!!    1 sumint    contains scale 1             " -tsa " and " -fbi "   !
!!    2 sumintsa  contains tsa approximation     " -tsa "              !
!!    3 sumintp   contains scale 2               " -fbi "              !
!!    4 sumintx   contains cross interactions    " -fbi "              !
!!                     between scales 1 and 2                          !
!!                                                                     !
!!    comparison to full integral can be made between                  !
!!    tsa=sumint+sumintsa  versus  fbi=sumint+sumintp+sumintx          !
!!    -----------------------------------------------------------------#
!!
!!
!!    3. Parameters :
!!
!!    Parameter list  ( all I/O arrays are now Public)
!!    ------------------------------------------------------------------
!!      nrng     Int.  I   number of wavenumbers
!!      nang     Int.  I   number of directions
!!      NZZ      Int.  I   NZZ = nrng * (nrng+1) / 2
!!      npts     Int.  I   number of points around locus
!!
!!      frqa     R.A.  I   radian frequencies (Hz);   dim=(nrng)
!!      pha      R.A.  I   = k*dk*dtheta          ;   dim=(nrng)
!!      dens     R.A.  I   lrg-scl Action density (k,theta);
!!      dens2    R.A.  I   Sml-scl Action density (k,theta);
!!                         dens and dens2 dimension = (nrng,nang)
!!
!!      kref2    I.A.  I   integration grid geometry variable;
!!      kref4    I.A.  I   ditto
!!      jref2    I.A.  I   ditto
!!      jref4    I.A.  I   ditto
!!      wtk2     R.A.  I   ditto
!!      wtk4     R.A.  I   ditto
!!      wta2     R.A.  I   ditto
!!      wta4     R.A.  I   ditto
!!      tfac2    R.A.  I   ditto
!!      tfac4    R.A.  I   ditto
!!                         dimension = (npts,nang,NZZ)
!!
!!      grad     R.A.  I   = C * H * g**2 * ds / |dW/dn|;
!!                         dimension = (npts,nang,NZZ)
!!
!!      sumint   R.A.  O   contains scale 1 contribution to Snl -tsa and -fbi
!!      sumintsa R.A.  O   contains tsa approximation "  "   "  -tsa
!!      sumintp  R.A.  O   contains scale 2 contribution to Snl -fbi
!!      sumintx  R.A.  O   contains cross interactions " "   "  -fbi
!!                         all  dimension = (nrng,nang)
!!
!!      for -tsa
!!      tsa      R.A.  O   Snl-tsa = sumint + sumintsa
!!      diag     R.A.  O   Snl-tsa diagonal term = [dN/dn1]
!!
!!      for -fbi
!!      fbi      R.A.  O   Snl-fbi = sumint + sumintp  + sumintx
!!      diag2    R.A.  O   Snl-fbi diagonal term = [dN/dn1]
!!    ------------------------------------------------------------------
!!    ----------------------------------------------------------------72
!!    ==================================================================
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    The included file 'parfile_2.f' (below) contains just            !
!!    parameters that sets the dimensions nrng, nang, npts and NZZ     !
!!    for some arrays (below)                                          !
!!                                                                     !
!!    parameters in 'parfile_2.f' are:                                 !
!!    nrng  =   67 = maximum number of rings                           !
!!    nang  =   36 = maximum number of angles                          !
!!    npts  =   30 = maximum number of points around locus             !
!!    NZZ   = 2278 =  NZZ = (nrng*(nrng+1))/2                          !
!!                                                                     !
!b    include 'parfile_2.f'
!!                                                                     !
!!    -----------------------------------------------------------------#
!!
!!------------------------------------------------------------------------------
!!--------------------------------------------------------------------72------80
!!==============================================================================
!!
!!
!!
      IMPLICIT NONE
!!
!!
      integer, intent(in)  :: nrng, nang, npts, NZZ
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      real,    intent(in)  :: frqa(nrng)                         !* /b/
      real,    intent(in)  :: pha(nrng)                          !* /c/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!pub  integer, intent(in)  :: kref2(npts,nang,NZZ),kref4(npts,nang,NZZ)
!pub  integer, intent(in)  :: jref2(npts,nang,NZZ),jref4(npts,nang,NZZ)
!pub  real,    intent(in)  ::  wtk2(npts,nang,NZZ), wtk4(npts,nang,NZZ)
!pub  real,    intent(in)  ::  wta2(npts,nang,NZZ), wta4(npts,nang,NZZ)
!pub  real,    intent(in)  :: tfac2(npts,nang,NZZ),tfac4(npts,nang,NZZ)
!!                                                               !* /fr/
!pub  real,    intent(in)  :: grad(npts,nang,NZZ)                !* /dn/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!pub  real,    intent(in)  :: dens(nrng,nang)                    !* /c/
!pub  real,    intent(in)  :: dens2(nrng,nang)                   !* /z/
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!
!!    for -tsa
      real,    intent(out) :: sumint(nrng,nang)                  !* /dn/
      real,    intent(out) :: sumintsa(nrng,nang)                !* /z/
!pub  real,    intent(out) :: tsa(nrng,nang), diag(nrng,nang)
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!    for -fbi
      real,    intent(out) :: sumintp(nrng,nang)                 !* /z/
      real,    intent(out) :: sumintx(nrng,nang)                 !* /z/
!pub  real,    intent(out) :: fbi(nrng,nang), diag2(nrng,nang)
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ---------------------::-----------------------------------------72
!!    ==================================================================
!!
!!
!!
!!    Local variables
!!
      integer              :: irng,krng, iang,kang, ipt, iizz,izz
      integer              :: iamin                      ! index of min angle
      integer              :: ia2, ia2p, k2, k2p
      integer              :: ia4, ia4p, k4, k4p
      integer              :: nref, nklimit, nalimit
!b    integer              :: mat1(nrng,nang), mat2(nrng,nang)
!!
!!    for both -tsa and -fbi
      real                 :: d1, d3,  d2, d4,  dp1, dp3, dz4, dz5
      real                 :: dx13, ds13,  dxp13, dsp13
      real                 :: dgm,  t31,   tr31
      real                 :: w2, w2p, wa2, wa2p, d2a, d2b, tt2
      real                 :: w4, w4p, wa4, wa4p, d4a, d4b, tt4
!!
!!    for -tsa
      real                 :: dz2a, dz3a,      ttsa, trtsa
      real                 :: diagk1, diagk3,  ddn1, ddn3
!!
!!    for -fbi
      real                 :: dp2, dp4, dz1, dz2, dz3, dz6, dz7, dz8
      real                 :: d2pa, d4pa, d2pb, d4pb
      real                 :: dgmp, tp31,trp31, dzsum, txp31,trx31
      real                 :: diag2k1, diag2k3, dd2n1, dd2n3
!!
!!    for -fbi
!!    Bash added 4 new terms for a full expression of diagonal terms:
!!    ddpi = di + dpi for i = 1,4
      real                 :: ddp1, ddp2,  ddp3, ddp4
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ---------------------::-----------------------------------------72
!!    ##################################################################
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!    Bash; hardwire these two parameters
      nklimit = 6
      nalimit = 6
!!
!!
!!    Bash; initialize returned arrays here instead of below
!!
      do 22 irng=1,nrng
        do 21 iang=1,nang
!!
!!        for both -tsa and -fbi
!!        sumint is now initialized here instead of below!
          sumint(irng,iang)   = 0.0
!!
!!        for -tsa
!!        sumintsa are now initialized here instead of below!
          sumintsa(irng,iang) = 0.0
          tsa(irng,iang)      = 0.0
          diag(irng,iang)     = 0.0
!!
!!        for -fbi
!!        sumintp and sumintx are now initialized here instead of below!
          sumintp(irng,iang)  = 0.0
          sumintx(irng,iang)  = 0.0
          fbi(irng,iang)      = 0.0
          diag2(irng,iang)    = 0.0
!!
  21    continue
  22  continue
!!    ------------------------------------------------------------------
!!    ------------------------------------------------------------------
!!
!!
!!    Bash; Test matrices mat1 and mat2 are to monitor where
!!    the 2 pairs of indices (irng,iang) and (krng,kang) have been
!!          = 0 means not touched (remains as initialized)
!!          = n means # of times being touched (touched or filled)
!b    do 24 irng=1,nrng
!b      do 23 iang=1,nang
!b        mat1(irng,iang) = 0
!b        mat2(irng,iang) = 0
!b23    continue
!b24  continue
!!
!!    Bash; wrote these matrices to 2 tapes
!b    open(1, file='mat1.plt', status='unknown')
!b    open(2, file='mat2.plt', status='unknown')
!!    ------------------------------------------------------------------
!!    ##################################################################
!!
!!
!!
      ddn1    = 0.0                           !* for -tsa diag [dN/dn1]
      ddn3    = 0.0                           !* for -tsa diag [dN/dn3]
!!
      dd2n1   = 0.0                           !* for -fbi diag [dN/dn1]
      dd2n3   = 0.0                           !* for -fbi diag [dN/dn3]
!!
!!
!!    50
      do 50 irng=1,nrng
!!
        iizz = (nrng-1)*(irng-1) - ((irng-2)*(irng-1))/2
!!      ----------------------------------------------------------------
!!
!!      60
        do 60 iang=1,nang
!!
          d1   = dens(irng,iang)
          dp1  = dens2(irng,iang)
          ddp1 = d1 + dp1           !! for full expression of diag. term
!!
!!        Bash; move the initialization up (see above)
!!        for -tsa
!x        sumint(irng,iang)   = 0.0
!x        sumintsa(irng,iang) = 0.0
!x        tsa(irng,iang)      = 0.0
!x        diag(irng,iang)     = 0.0
!!
!!        for -fbi
!x        sumintp(irng,iang)  = 0.0
!x        sumintx(irng,iang)  = 0.0
!x        fbi(irng,iang)      = 0.0
!x        diag2(irng,iang)    = 0.0
!!        --------------------------------------------------------------
!!
!!
!!        70
          do 70 krng=irng,nrng
!!
!!          Bash; check5  be consistent with gridsetr
!!          moved here from below (was after do 80 kang=1,nang)
!!          and changed go to 80 into go to 70 (i.e. go to next krng)
            if ( frqa(krng)/frqa(irng) .gt. 4. ) go to 70  !* gridsetr
!b          if ( frqa(krng)/frqa(irng) .gt. 3. ) go to 70  !* original
!!
            izz = krng + iizz
!!          ------------------------------------------------------------
!!
!!          Bash; Grab index of the min angle to remove self interaction
            iamin = 1
            if ( krng.eq.irng ) iamin = iang + 1
!!          ------------------------------------------------------------
!!
!!
!!          80
!b          do 80 kang=1,nang
            do 80 kang=iamin,nang
!!
              d3   = dens(krng,kang)
              dp3  = dens2(krng,kang)
              ddp3 = d3 + dp3       !! for full expression of diag. term
!!
              nref = kang - iang + 1
              if ( nref .lt. 1 ) nref = nref + nang
!!
!!            Bash; check5  be consistent with gridsetr
!!                  and move this test above right after do 70 krng=irng,nrng
!x            if ( frqa(krng)/frqa(irng) .gt. 4. ) go to 80  !* gridsetr
!b            if ( frqa(krng)/frqa(irng) .gt. 3. ) go to 80  !* original
!!
!!            for both -tsa and -fbi
              dgm     = 0.0             !* not necessary
              t31     = 0.0             !* must be reset to 0.0
              tr31    = 0.0             !* not necessary
!!
!!            for -tsa 
              ttsa    = 0.0             !* must be reset to 0.0
              diagk1  = 0.0             !* must be reset to 0.0
              diagk3  = 0.0             !* must be reset to 0.0
!!
!!            for -fbi 
              dgmp    = 0.0             !* not necessary
              tp31    = 0.0             !* must be reset to 0.0
              trp31   = 0.0             !* not necessary
              dzsum   = 0.0             !* not necessary
              txp31   = 0.0             !* must be reset to 0.0
              trx31   = 0.0             !* not necessary
              diag2k1 = 0.0             !* must be reset to 0.0
              diag2k3 = 0.0             !* must be reset to 0.0
!!
!!            for both -tsa and -fbi
              dx13  = d1  * d3
              ds13  = d3  - d1
              dxp13 = dp1 * dp3
              dsp13 = dp3 - dp1
!!            ----------------------------------------------------------
!!
!!
!!            90
              do 90 ipt=1,npts                       !* begin locus loop
!!
!!              save time by skipping insignificant contributions
                if ( grad(ipt,nref,izz) .lt. 1.e-30 ) go to 90
!!
!!2             estimation of density for wave #2
                k2  = kref2(ipt,nref,izz)
                k2p = k2 + 1
                w2  = wtk2(ipt,nref,izz)
                w2p = 1. - w2
!!
                ia2 = iang + jref2(ipt,nref,izz)
                if ( ia2 .gt. nang )  ia2  = ia2  - nang
!!
                ia2p = ia2 + 1
                if ( ia2p .gt. nang ) ia2p = ia2p - nang
!!
                wa2  = wta2(ipt,nref,izz)
                wa2p = 1. - wa2
                d2a  = w2 * dens(k2,ia2)   + w2p * dens(k2p,ia2)
                d2pa = w2 * dens2(k2,ia2)  + w2p * dens2(k2p,ia2)
                d2b  = w2 * dens(k2,ia2p)  + w2p * dens(k2p,ia2p)
                d2pb = w2 * dens2(k2,ia2p) + w2p * dens2(k2p,ia2p)
                tt2  = tfac2(ipt,nref,izz)
                d2   = (wa2*d2a  + wa2p*d2b)  * tt2
                dp2  = (wa2*d2pa + wa2p*d2pb) * tt2   !* for -fbi
                ddp2 = d2 + dp2     !! for full expression of diag. term
!!              --------------------------------------------------------
!!
!!4             estimation of density for wave #4
                k4  = kref4(ipt,nref,izz)
                k4p = k4 + 1
                w4  = wtk4(ipt,nref,izz)
                w4p = 1. - w4
!!
                ia4 = iang + jref4(ipt,nref,izz)
                if ( ia4 .gt. nang )  ia4  = ia4  - nang
!!
                ia4p= ia4 + 1
                if ( ia4p .gt. nang ) ia4p = ia4p - nang
!!
                wa4  = wta4(ipt,nref,izz)
                wa4p = 1. - wa4
                d4a  = w4*dens(k4,ia4)   + w4p*dens(k4p,ia4)
                d4pa = w4*dens2(k4,ia4)  + w4p*dens2(k4p,ia4)
                d4b  = w4*dens(k4,ia4p)  + w4p*dens(k4p,ia4p)
                d4pb = w4*dens2(k4,ia4p) + w4p*dens2(k4p,ia4p)
                tt4  = tfac4(ipt,nref,izz)
                d4   = (wa4*d4a  + wa4p*d4b)  * tt4
                dp4  = (wa4*d4pa + wa4p*d4pb) * tt4   !* for -fbi
                ddp4 = d4 + dp4     !! for full expression of diag. term
!!              --------------------------------------------------------
!!
!!
!!              for both -tsa and -fbi
                dgm  = dx13*(d4-d2) + ds13*d4*d2
                t31  = t31  + dgm  * grad(ipt,nref,izz)
!!
!!              for -fbi
                dgmp = dxp13*(dp4-dp2) + dsp13*dp4*dp2
                tp31 = tp31 + dgmp * grad(ipt,nref,izz)
!!
!!
!!              for -tsa : -diag
!!              use this expression for the diagonal term
!!              whose derivation neglect "dp2" & "dp4"
                ddn1 = (d3+dp3)*(d4-d2) - d4*d2              !* dN/dn1
                ddn3 = (d1+dp1)*(d4-d2) + d4*d2              !* dN/dn3
!!
!!              or use this expresiion for the diagonal term
!!              whose derivation neglects all the "dp#'s" and 
!!              keeps only the large scale terms the "d#'s"
!b              ddn1 = d3*(d4-d2) - d4*d2                    !* dN/dn1
!b              ddn3 = d1*(d4-d2) + d4*d2                    !* dN/dn3
!!              --------------------------------------------------------
!!
                diagk1  = diagk1 + ddn1 * grad(ipt,nref,izz)
                diagk3  = diagk3 + ddn3 * grad(ipt,nref,izz)
!!              ========================================================
!!
!!
!!              diagonal term for -fbi : -diag2
!!
!!              use the full expression for the diagonal terms
!!              whose derivation keeps all large + small scale
                dd2n1 = ddp3*(ddp4-ddp2) - ddp4*ddp2          !* dN/dn1
                dd2n3 = ddp1*(ddp4-ddp2) + ddp4*ddp2          !* dN/dn3
!!              --------------------------------------------------------
!!
                diag2k1 = diag2k1 + dd2n1 * grad(ipt,nref,izz)
                diag2k3 = diag2k3 + dd2n3 * grad(ipt,nref,izz)
!!              ========================================================
!!
!!
!!              for -fbi
                dz1  = dx13    * (dp4-dp2)
                dz2  = d1*dp3  * ((d4-d2)+(dp4-dp2))
                dz3  = d3*dp1  * ((d4-d2)+(dp4-dp2))
!!
!!              for both -tsa and -fbi
                dz4  = dxp13   * (d4-d2)
                dz5  = d2*d4   *  dsp13
!!
!!              for -fbi
                dz6   = d2*dp4  * (ds13+dsp13)
                dz7   = d4*dp2  * (ds13+dsp13)
                dz8   = dp2*dp4 *  ds13
                dzsum = dz1 + dz2 + dz3 + dz4 + dz5 + dz6 + dz7 + dz8
                txp31 = txp31 + dzsum * grad(ipt,nref,izz)
!!              --------------------------------------------------------
!!
!!
!!              Cross-interactions between parametric and perturbation
!!              Note: These occur only if k3 is close enough to k1
!!              for -tsa
!!
!!              Bash; check6
!!              Bash; add extra check on (nang-nalimit)
!!                    and change .lt. to .le.           -  NO
!!              --------------------------------------------------------
!!
!!              for -tsa
!b              if ( iabs(irng-krng).lt.nklimit .and.                 &
!b                   iabs(iang-kang).lt.nalimit )    then   !* original
!!
                if (     (krng-irng).lt.nklimit .and.                 &
                   ( iabs(iang-kang).lt.nalimit .or.                  &
                     iabs(iang-kang).gt.(nang-nalimit) ) )  then
!!
                  dz2a = d1*dp3 * (d4-d2)
                  dz3a = d3*dp1 * (d4-d2)
!!
                  ttsa = ttsa + (dz4+dz5+dz2a+dz3a)*grad(ipt,nref,izz)
!!
                endif
!!                                                        !* Bash
!!              --------------------------------------------------------
!!
!!
  90          continue                        !* end of ipt (locus) loop
!!
!!            ----------------------------------------------------------
!!
!!
!!            multiply the following components by factor 2. in here
              tr31  = 2. * t31                !* for both -tsa and -fbi
              trtsa = 2. * ttsa               !* for -tsa
              trp31 = 2. * tp31               !* for -fbi
              trx31 = 2. * txp31              !* for -fbi
!!
              diagk1  = 2. * diagk1           !* for Snl-TSA diag
              diagk3  = 2. * diagk3           !* for Snl-TSA diag
              diag2k1 = 2. * diag2k1          !* for Snl-fbi diag2
              diag2k3 = 2. * diag2k3          !* for Snl-fbi diag2
!!            ----------------------------------------------------------
!!
!!
!!            for both -tsa and -fbi
              sumint(irng,iang)  = sumint(irng,iang)  + tr31*pha(krng)
              sumint(krng,kang)  = sumint(krng,kang)  - tr31*pha(irng)
!!
!!            for -tsa
              sumintsa(irng,iang)= sumintsa(irng,iang)+ trtsa*pha(krng)
              sumintsa(krng,kang)= sumintsa(krng,kang)- trtsa*pha(irng)
!!
!!            diagonal term for -tsa : -diag
              diag(irng,iang) = diag(irng,iang)  + diagk1*pha(krng)
              diag(krng,kang) = diag(krng,kang)  - diagk3*pha(irng)
!!
!!            for -fbi
              sumintp(irng,iang) = sumintp(irng,iang) + trp31*pha(krng)
              sumintp(krng,kang) = sumintp(krng,kang) - trp31*pha(irng)
!!
              sumintx(irng,iang) = sumintx(irng,iang) + trx31*pha(krng)
              sumintx(krng,kang) = sumintx(krng,kang) - trx31*pha(irng)
!!
!!            diagonal term for -fbi : -diag2
              diag2(irng,iang) = diag2(irng,iang) + diag2k1*pha(krng)
              diag2(krng,kang) = diag2(krng,kang) - diag2k3*pha(irng)
!!            ----------------------------------------------------------
!!
!b            mat1(irng,iang) = mat1(irng,iang) + 1 !* touched index loc.
!b            mat2(krng,kang) = mat2(krng,kang) + 1 !* touched index loc.
!!            ----------------------------------------------------------
!!
  80        continue                              !* end of kang loop
!!
  70      continue                                !* end of krng loop
!!
  60    continue                                  !* end of iang loop
!!
  50  continue                                    !* end of irng loop
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!    Calc. tsa() and fbi() to be returned. diag() and diag2() should be ok
      do 52 irng=1,nrng
      do 51 iang=1,nang
         fbi(irng,iang) = sumint(irng,iang) + sumintp(irng,iang) +    &
                                              sumintx(irng,iang)
         tsa(irng,iang) = sumint(irng,iang) + sumintsa(irng,iang)
  51  continue
  52  continue
!!
!!
!!    Bash; wrote these matrices to 2 tapes
!!    The results of this test,  all the indices of mat1 and mat2 
!!    have been touched => all the elements of  mat1 and mat2 are 1's
!!
!b    do 55 iang=1,nang
!b      write(1,901)  (mat1(irng,iang), irng=1,nrng)
!b      write(2,901)  (mat2(irng,iang), irng=1,nrng)
!901    format(2x,51I1)
!b55  continue
!!
!b    close(1)
!b    close(2)
!!    ------------------------------------------------------------------
!!
      return
!!
!b    end                                         !* end SUBROUTINE snlr
!!
      END SUBROUTINE snlr
!!
!!------------------------------------------------------------------------------
!!==============================================================================
!!
      END MODULE W3SNLXMD
!!
!!------------------------------------------------------------------------------
!!==============================================================================
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
