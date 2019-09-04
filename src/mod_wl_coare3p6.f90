! AeroBulk / 2016 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
! LOLO: if SST measured at a relatively shallow depth, say 1 m, then it should already include the warm-layer effect???
!      => correction makes only sense if SST is taken relatively deep and that heat flux is so big and wind so weak that the warm layer is becomes thinner ??? 
!
!
MODULE mod_wl_coare3p6
   !!====================================================================================
   !!       Warm-Layer correction of SST (if needed)
   !!
   !!       Routine "wl_coare3p6" maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, 2019
   !!
   !!====================================================================================
   USE mod_const   !: physical and othe constants
   USE mod_phymbl  !: thermodynamics

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: WL_COARE3P6

   REAL(wp), PARAMETER :: rich   = 0.65_wp   !: critical Richardson number
   REAL(wp), PARAMETER :: z_sst  = 1._wp     !: depth at which bulk SST is taken...
   REAL(wp), PARAMETER :: dz_max = 19._wp    !: maximum depth of warm layer (adjustable)
   REAL(wp), PARAMETER :: Qabs_thr = 50._wp  !: threshold for heat flux absorbed in WL
   REAL(wp), PARAMETER :: zfs0   = 0.5_wp    !: initial value of solar flux absorption

   !LOGICAL, PUBLIC, SAVE :: l_wl_c36_never_called

CONTAINS

   SUBROUTINE WL_COARE3P6( pQsw, pQnsol, pTau, pSST, pdT, pTau_ac, pQ_ac, plon, isd, rdt, &
      &                    Hwl )
      !!---------------------------------------------------------------------
      !!
      !!  Cool-Skin Warm-Layer scheme according to COARE 3.6 (Fairall et al, 2019)
      !!     ------------------------------------------------------------------
      !!
      !!  **   INPUT:
      !!
      !!     *pQsw*       surface net solar radiation into the ocean     [W/m^2] => >= 0 !
      !!     *pQnsol*     surface net non-solar heat flux into the ocean [W/m^2] => normally < 0 !
      !!     *pTau*       surface wind stress                            [N/m^2]
      !!     *pSST*       bulk SST at depth z_sst                        [K]
      !!     *plon*       longitude                                      [deg.E]
      !!     *isd*        current UTC time, counted in second since 00h of the current day
      !!     *rdt*        physical time step between two successive calls to this routine [s]
      !!
      !!   **  INPUT/OUTPUT:
      !!     *pdT*        dT due to warming at depth of pSST such that SST_actual = pSST + pdT
      !!
      !!   ** OPTIONAL OUTPUT:
      !!     *Hwl*        depth of warm layer [m]
      !!
      !!------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pQsw     ! surface net solar radiation into the ocean [W/m^2]     => >= 0 !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pQnsol   ! surface net non-solar heat flux into the ocean [W/m^2] => normally < 0 !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pTau     ! wind stress [N/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pSST     ! bulk SST at depth z_sst [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pdT      ! dT due to warming at depth of pSST such that pSST_actual = pSST + pdT
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pTau_ac  ! time integral / accumulated momentum Tauxdt => [N.s/m^2] (reset to zero every midnight)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pQ_ac    ! time integral / accumulated heat stored by the warm layer Qxdt => [J/m^2] (reset to zero every midnight)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: plon     ! longitude ! lolo
      INTEGER ,                     INTENT(in)    :: isd      ! current UTC time, counted in second since 00h of the current day
      REAL(wp),                     INTENT(in)    :: rdt      ! physical time step between two successive call to this routine [s]
      !!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out), OPTIONAL :: Hwl    ! depth of warm layer [m]
      !
      !
      INTEGER :: ji,jj,iflg
      !
      REAL(wp) :: dT_wl, dz_wl, zQabs, zfs
      REAL(wp) :: zqac, ztac
      REAL(wp) :: zAl, zcd1, zcd2

      REAL(wp) :: rlag_gw_h  ! local solar time lag in hours   / Greenwich meridian (lon==0) => ex: ~ -10.47 hours for Hawai

      INTEGER  :: ilag_gw_s, &  ! local solar time LAG in seconds / Greenwich meridian (lon==0) => ex: ~ INT( -10.47*3600. ) seconds for Hawai
         &        isd_sol,   &  ! local solar time in number of seconds since local solar midnight
         &        jl

      !! INITIALIZATION:

      pdT = 0._wp         ! dT initially set to 0._wp
      dT_wl = 0._wp       ! total warming (amplitude) in warm layer
      dz_wl = dz_max      ! initial depth set to max value
      zQabs  = 0._wp       ! total heat absorped in warm layer
      zfs   = zfs0        ! initial value of solar flux absorption

      DO jj = 1, jpj
         DO ji = 1, jpi

            !! Need to know the local solar time from longitude and isd
            !! *********************************************************
            rlag_gw_h = -1._wp * MODULO( ( 360._wp - MODULO(plon(ji,jj),360._wp) ) / 15._wp , 24._wp )
            rlag_gw_h = -1._wp * SIGN( MIN(ABS(rlag_gw_h) , ABS(MODULO(rlag_gw_h,24._wp))), rlag_gw_h + 12._wp )
            ilag_gw_s = INT( rlag_gw_h*3600._wp )
            isd_sol = MODULO( isd + ilag_gw_s , 24*3600 )
            !PRINT *, ' Lag in hours / Greenwich for local solar time =', rlag_gw_h
            !PRINT *, '     UTC     time in seconds:', isd
            !PRINT *, ' Local solar time in seconds:', isd_sol
            !PRINT *, '     UTC     time in hours:',   REAL(isd    ,wp)/3600._wp
            PRINT *, '  [WL_COARE3P6] Local solar time in hours:',   REAL(isd_sol,wp)/3600._wp
            !**********************************************************

            !*****  variables for warm layer  ***
            zAl   = 2.1e-5*(pSST(ji,jj)-rt0 + 3.2_wp)**0.79
            zcd1 = SQRT(2._wp*rich*rCp0_w/(zAl*grav*rho0_w))        !mess-o-constants 1
            zcd2 = SQRT(2._wp*zAl*grav/(rich*rho0_w))/(rCp0_w**1.5) !mess-o-constants 2

            !********************************************************
            !****  Compute apply warm layer  correction *************
            !********************************************************

            IF (isd_sol < rdt ) THEN    !re-zero at midnight ! LOLO improve: risky if real midnight (00:00:00) is not a time in vtime...
               PRINT *, '  [WL_COARE3P6] MIDNIGHT RESET !!!!, isd_sol =>', isd_sol
               zfs            = zfs0
               dz_wl          = dz_max
               pTau_ac(ji,jj) = 0._wp
               ztac           = 0._wp
               pQ_ac(ji,jj)   = 0._wp
               zqac           = 0._wp
               dT_wl          = 0._wp
            END IF

            IF ( isd_sol >= 21600 ) THEN  ! (21600 == 6am)

               PRINT *, '  [WL_COARE3P6] WE DO WL !'
               PRINT *, '  [WL_COARE3P6] isd_sol, pTau, pSST, pdT =', isd_sol, REAL(pTau(ji,jj),4), REAL(pSST(ji,jj),4), REAL(pdT(ji,jj),4)

               !************************************
               !****   set warm layer constants  ***
               !************************************

               zQabs = zfs*pQsw(ji,jj) + pQnsol(ji,jj)       ! tot heat absorbed in warm layer

               PRINT *, '  [WL_COARE3P6] rdt,  pQsw, pQnsol, zQabs =', rdt, REAL(pQsw(ji,jj),4), REAL(pQnsol(ji,jj),4), REAL(zQabs,4)

               IF ( zQabs >= Qabs_thr ) THEN         ! Check for threshold

                  PRINT *, '  [WL_COARE3P6] pTau_ac, pQ_ac =', REAL(pTau_ac(ji,jj),4), REAL(pQ_ac(ji,jj),4)

                  !pTau_ac(ji,jj) = pTau_ac(ji,jj) + MAX(.002_wp , pTau(ji,jj))*rdt      ! momentum integral
                  ztac = pTau_ac(ji,jj) + MAX(.002_wp , pTau(ji,jj))*rdt      ! updated momentum integral

                  IF ( pQ_ac(ji,jj) + zQabs*rdt > 0._wp ) THEN         !check threshold for warm layer existence
                     !******************************************
                     ! Compute the absorption profile
                     !******************************************
                     DO jl = 1, 5                           !loop 5 times for zfs
                        zfs = 1. - ( 0.28*0.014*(1. - EXP(-dz_wl/0.014)) + 0.27*0.357*(1. - EXP(-dz_wl/0.357)) &
                           &        + 0.45*12.82*(1-EXP(-dz_wl/12.82)) ) / dz_wl
                        zqac = pQ_ac(ji,jj) + (zfs*pQsw(ji,jj) + pQnsol(ji,jj))*rdt ! updated heat absorbed
                        IF (zqac <= 0._wp) STOP'ERROR: zqac <= 0 !!! #1'
                        dz_wl = MIN( dz_max , zcd1*ztac/SQRT(zqac)) ! Warm-layer depth (normally: zqac > 0 !)
                     END DO

                  ELSE
                     !***********************
                     ! Warm layer wiped out
                     !***********************
                     zfs    = 0.75
                     dz_wl  = dz_max
                     zqac   = pQ_ac(ji,jj) + (zfs*pQsw(ji,jj) + pQnsol(ji,jj))*rdt ! updated heat absorbed
                  END IF

                  ! normally: zqac > 0 !

                  !*******  compute dt_wl  ******
                  !IF (zqac > 0._wp) dT_wl = zcd2*zqac**1.5/ztac  ! dT_wl remains = 0 otherwize...
                  dT_wl = zcd2*zqac**1.5/ztac * MAX(zqac/ABS(zqac),0._wp)  !! => IF(zqac>0._wp): dT_wl=zcd2*zqac**1.5/ztac ; ELSE: dT_wl=0.



               END IF ! IF ( zQabs >= Qabs_thr )

               ! Warm layer correction
               iflg = INT( 0.5 + SIGN( 0.5 , z_sst-dz_wl ) ) ! => 1 when z_sst>dz_wl (pdT(ji,jj) = dT_wl) | 0 when z_sst<dz_wl (pdT(ji,jj) = dT_wl*z_sst/dz_wl)
               pdT(ji,jj) = dT_wl * ( iflg + (1-iflg)*z_sst/dz_wl )

            END IF ! IF ( isd_sol >= 21600 ) THEN  ! (21600 == 6am)

            IF ( (zQabs >= Qabs_thr).AND.(isd_sol >= 21600) ) THEN
               pQ_ac(ji,jj)   = zqac ! Updating pQ_ac, heat integral
               pTau_ac(ji,jj) = ztac !
            END IF

            IF ( PRESENT(Hwl) ) Hwl(ji,jj) = dz_wl

         END DO
      END DO

   END SUBROUTINE WL_COARE3P6

END MODULE mod_wl_coare3p6
