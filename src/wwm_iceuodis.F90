      ! author: Lorenzo Mentaschi
      !  Zero-order approximation of the wave attenuation due to ice concentration.
      !  ICEUODIS stands for Ice Unresolved Obstacles Dissipation
      !  Based on UOST (Mentaschi et al. 2015, 2018, 2020)
      !  
      !  From the ice concentration an isotropic transparency coefficient is estimated.
      !  To simplify the problem, for now it is assumed that beta==alpha
      !  which means that all the energy is dissipated in the current cell.
      !  This gets rid of the shadow problems.
      !  This is inaccurate for the current cell, but should dissipate all the energy needed.
      !  Furthermore, local wave growth (which reduces the effect of unresolved obstacles) 
      !  is neglected
      !
      !  In the future, if this approximation will not be enough, 
      !  some improvement could be introduced:
      !  - a beta different from alpha could be estimated by assuming 
      !  a uniform distribution of the ice in the cell, or by loading
      !  the distribution of the ice,
      !  and the shadow could be estimated for the neighboring cells.
      !  - local wave growth could be taken into account (see the psi function in UOST)
      !  - if there is information on the size of the ice flows, the transparency coeff. 
      !  could be made frequency-dependent.
      SUBROUTINE ICEUODIS_SRCTRM(IP, SPEC, S, D)
         USE DATAPOOL

         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP ! this is the node id in the local partition
         REAL(rkind), INTENT(IN) :: SPEC(NUMSIG, NUMDIR)
         REAL(rkind), INTENT(OUT) :: S(NUMSIG, NUMDIR), D(NUMSIG, NUMDIR)
         REAL(rkind) :: ICEC, OBSTSECTION, BETA, CELLAREA, CELLSIZE,& 
             DEG2M, CGI, GAM
         REAL(rkind) :: GAMMAUP = 200
         INTEGER  :: IK
   
         S = 0
         D = 0
         ICEC = MAX(MIN(ICECONC(IP), 1.), 0.)

         IF (ICEC .LE. THR) RETURN

         ! computing the transparency coefficient
         ! the total obstruction coefficient is given by sqrt(concentration)
         ! the total transparency alpha is given by 1-obstruction
         ! here we assume that beta==alpha
         OBSTSECTION = SQRT(ICEC)
         BETA = 1-OBSTSECTION

         CELLAREA = SI(IP)
         IF (LSPHE) THEN
           DEG2M = REARTH*PI/180._rkind
           CELLAREA = CELLAREA*DEG2M*DEG2M ! converting cell area to meters
         END IF
         ! cellsize computed as the ray of the equivalent circle
         CELLSIZE = SQRT(CELLAREA/PI) 

         GAM = (1 - BETA)/BETA
         GAM = MIN(GAM, GAMMAUP)
   
         DO IK = 1,NUMSIG
           CGI = CG(IK,IP)
           D(IK, :) = - CGI/CELLSIZE * GAM
         END DO
         S = D*SPEC
         
      END SUBROUTINE ICEUODIS_SRCTRM

     
