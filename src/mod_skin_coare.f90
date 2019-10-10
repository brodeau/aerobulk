! AeroBulk / 2016 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_skin_coare
   !!====================================================================================
   !!       Cool-Skin and warm-layer correction of SST
   !!    Cool-skin and warm-layer parametrization (Fairall et al. 1996)
   !!
   !!       Routine "cs_coare" and "wl_coare" maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, 2019
   !!
   !!====================================================================================
   USE mod_const   !: physical and othe constants
   USE mod_phymbl  !: thermodynamics

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: CS_COARE, WL_COARE

   !! Cool-skin related parameters:
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: delta_vl  !: thickness of the surface viscous layer (in the water) right below the air-sea interface [m]

   !! Warm-layer related parameters:
   REAL(wp), PARAMETER, PUBLIC :: H_wl_max = 20._wp    !: maximum depth of warm layer (adjustable)
   !
   REAL(wp), PARAMETER :: rich   = 0.65_wp   !: critical Richardson number
   !
   !REAL(wp), PARAMETER :: Qabs_thr = 50._wp  !: threshold for heat flux absorbed in WL !DO NOT GET IT !!!!

   REAL(wp), PARAMETER :: zfr0   = 0.5_wp    !: initial value of solar flux absorption
   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: Tau_ac  ! time integral / accumulated momentum Tauxdt => [N.s/m^2] (reset to zero every midnight)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: Qnt_ac    ! time integral / accumulated heat stored by the warm layer Qxdt => [J/m^2] (reset to zero every midnight)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: H_wl     ! depth of warm-layer [m]
   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE CS_COARE( pQsw, pQnsol, pustar, pSST, pQlat,  pdT )
      !!---------------------------------------------------------------------
      !!
      !!  Cool-Skin scheme according to Fairall et al. 1996, revisited for COARE 3.6 (Fairall et al., 2019)
      !!     ------------------------------------------------------------------
      !!
      !!  **   INPUT:
      !!     *pQsw*       surface net solar radiation into the ocean     [W/m^2] => >= 0 !
      !!     *pQnsol*     surface net non-solar heat flux into the ocean [W/m^2] => normally < 0 !
      !!     *pustar*     friction velocity u*                           [m/s]
      !!     *pSST*       bulk SST (taken at depth gdept_1d(1))          [K]
      !!     *pQlat*      surface latent heat flux                       [K]
      !!
      !!  **  INPUT/OUTPUT:
      !!     *pdT*  : as input =>  previous estimate of dT cool-skin
      !!              as output =>  new estimate of dT cool-skin
      !!
      !!------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pQsw   ! net solar a.k.a shortwave radiation into the ocean (after albedo) [W/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pQnsol ! non-solar heat flux to the ocean [W/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pustar  ! friction velocity, temperature and humidity (u*,t*,q*)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pSST ! bulk SST [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pQlat  ! latent heat flux [W/m^2]
      !!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pdT    ! dT due to cool-skin effect
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) :: zQnet, zQnsol, zlamb, zdelta, zalpha_w, zfr, &
         &        zz1, zz2, zus, &
         &        ztf
      !!---------------------------------------------------------------------

      DO jj = 1, jpj
         DO ji = 1, jpi

            zalpha_w = alpha_sw( pSST(ji,jj) ) ! thermal expansion coefficient of sea-water (SST accurate enough!)

            zQnsol = MAX( 1._wp , - pQnsol(ji,jj) ) ! Non-solar heat loss to the atmosphere

            zdelta = delta_vl(ji,jj)   ! using last value of delta

            !! Fraction of the shortwave flux absorbed by the cool-skin sublayer:
            zfr   = MAX( 0.137_wp + 11._wp*zdelta - 6.6E-5_wp/zdelta*(1._wp - EXP(-zdelta/8.E-4_wp)) , 0.01_wp ) ! Eq.16 (Fairall al. 1996b) /  !LB: why 0.065 and not 0.137 like in the paper??? Beljaars & Zeng use 0.065, not 0.137 !

            zQnet = MAX( 1._wp , zQnsol - zfr*pQsw(ji,jj) )  ! Total cooling at the interface

            ztf = 0.5 + SIGN(0.5_wp, zQnet) ! Qt > 0 => cooling of the layer => ztf = 1
            !                                 Qt < 0 => warming of the layer => ztf = 0

            !! Term alpha*Qb (Qb is the virtual surface cooling inc. buoyancy effect of salinity due to evap):
            zz1 = zalpha_w*zQnet - 0.026*MIN(pQlat(ji,jj),0._wp)*rCp0_w/rLevap  ! alpha*(Eq.8) == alpha*Qb "-" because Qlat < 0
            !! LB: this terms only makes sense if > 0 i.e. in the cooling case
            !! so similar to what's done in ECMWF:
            zz1 = MAX(0._wp , zz1)    ! 1. instead of 0.1 though ZQ = MAX(1.0,-pQlw(ji,jj) - pQsen(ji,jj) - pQlat(ji,jj))

            zus = MAX(pustar(ji,jj), 1.E-4_wp) ! Laurent: too low wind (u*) might cause problem in stable cases:
            zz2 = zus*zus * roadrw
            zz2 = zz2*zz2
            zlamb =  6._wp*( 1._wp + (rcst_cs*zz1/zz2)**0.75 )**(-1./3.) ! Lambda (Eq.14) (Saunders)

            ! Updating molecular sublayer thickness (delta):
            zz2    = rnu0_w/(SQRT(roadrw)*zus)
            zdelta =      ztf    *          zlamb*zz2   &  ! Eq.12 (when alpha*Qb>0 / cooling of layer)
               &    + (1._wp - ztf) * MIN(0.007_wp , 6._wp*zz2 )    ! Eq.12 (when alpha*Qb<0 / warming of layer)
            !LB: changed 0.01 to 0.007

            !! Once again with the new zdelta:
            zfr   = MAX( 0.137_wp + 11._wp*zdelta - 6.6E-5_wp/zdelta*(1._wp - EXP(-zdelta/8.E-4_wp)) , 0.01_wp ) ! Solar absorption / Eq.16 (Fairall al. 1996b)
            zQnet = MAX( 1._wp , zQnsol - zfr*pQsw(ji,jj) ) ! Total cooling at the interface

            !! Update!
            pdT(ji,jj) =  MIN( - zQnet*zdelta/rk0_w , 0._wp )   ! temperature increment
            delta_vl(ji,jj) = zdelta

         END DO
      END DO

   END SUBROUTINE CS_COARE




   SUBROUTINE WL_COARE( kt,  pQsw, pQnsol, pTau, pSST, plon, isd, iwait,  pdT, &
      &                    Hwl, mask_wl )
      !!---------------------------------------------------------------------
      !!
      !!  Warm-Layer scheme according to COARE 3.6 (Fairall et al, 2019)
      !!     ------------------------------------------------------------------
      !!
      !!  **   INPUT:
      !!     *pQsw*       surface net solar radiation into the ocean     [W/m^2] => >= 0 !
      !!     *pQnsol*     surface net non-solar heat flux into the ocean [W/m^2] => normally < 0 !
      !!     *pTau*       surface wind stress                            [N/m^2]
      !!     *pSST*       bulk SST  (taken at depth gdept_1d(1))         [K]
      !!     *plon*       longitude                                      [deg.E]
      !!     *isd*        current UTC time, counted in second since 00h of the current day
      !!     *iwait*      if /= 0 then wait before updating accumulated fluxes, we are within a converging itteration loop...
      !!
      !!  **   OUTPUT:
      !!     *pdT*   dT due to warm-layer effect => difference between "almost surface (right below viscous layer) and depth of bulk SST
      !!---------------------------------------------------------------------
      !!
      !!   ** OPTIONAL OUTPUT:
      !!     *Hwl*        depth of warm layer [m]
      !!     *mask_wl*    mask for possible existence of a warm-layer (1) or not (0)
      !!
      !!------------------------------------------------------------------
      INTEGER ,                     INTENT(in)  :: kt
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pQsw     ! surface net solar radiation into the ocean [W/m^2]     => >= 0 !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pQnsol   ! surface net non-solar heat flux into the ocean [W/m^2] => normally < 0 !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pTau     ! wind stress [N/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pSST     ! bulk SST at depth gdept_1d(1) [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: plon     ! longitude [deg.E]
      INTEGER ,                     INTENT(in)  :: isd      ! current UTC time, counted in second since 00h of the current day
      INTEGER ,                     INTENT(in)  :: iwait    ! if /= 0 then wait before updating accumulated fluxes
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pdT      ! dT due to warm-layer effect => difference between "almost surface (right below viscous layer) and depth of bulk SST
      !!
      REAL(wp),   DIMENSION(jpi,jpj), INTENT(out), OPTIONAL :: Hwl     ! depth of warm layer [m]
      INTEGER(1), DIMENSION(jpi,jpj), INTENT(out), OPTIONAL :: mask_wl ! mask for possible existence of a warm-layer (1) or not (0)
      !
      !
      INTEGER :: ji,jj
      !
      REAL(wp) :: zdT_wl, zQabs, zfr, zdz
      REAL(wp) :: zqac, ztac
      REAL(wp) :: zalpha_w, zcd1, zcd2, flg
      !!---------------------------------------------------------------------

      REAL(wp) :: rlag_gw_h, &  ! local solar time lag in hours   / Greenwich meridian (lon==0) => ex: ~ -10.47 hours for Hawai
         &        rhr_sol       ! local solar time in hours since midnight

      INTEGER  :: ilag_gw_s, &  ! local solar time LAG in seconds / Greenwich meridian (lon==0) => ex: ~ INT( -10.47*3600. ) seconds for Hawai
         &        isd_sol,   &  ! local solar time in number of seconds since local solar midnight
         &        jl

      !! INITIALIZATION:
      pdT(:,:) = 0._wp    ! dT initially set to 0._wp
      zdT_wl = 0._wp       ! total warming (amplitude) in warm layer
      zQabs  = 0._wp       ! total heat flux absorped in warm layer
      zfr    = zfr0        ! initial value of solar flux absorption
      ztac   = 0._wp
      zqac   = 0._wp
      IF ( PRESENT(mask_wl) ) mask_wl(:,:) = 0

      DO jj = 1, jpj
         DO ji = 1, jpi

            zdz   = MAX( MIN(H_wl(ji,jj),H_wl_max) , 0.1_wp) ! depth of warm layer

            !! Need to know the local solar time from longitude and isd:
            rlag_gw_h = -1._wp * MODULO( ( 360._wp - MODULO(plon(ji,jj),360._wp) ) / 15._wp , 24._wp )
            rlag_gw_h = -1._wp * SIGN( MIN(ABS(rlag_gw_h) , ABS(MODULO(rlag_gw_h,24._wp))), rlag_gw_h + 12._wp )
            ilag_gw_s = INT( rlag_gw_h*3600._wp )
            isd_sol = MODULO( isd + ilag_gw_s , 24*3600 )
            rhr_sol = REAL( isd_sol , wp) / 3600._wp
            !PRINT *, ' Lag in hours / Greenwich for local solar time =', rlag_gw_h
            !PRINT *, '     UTC     time in seconds:', isd
            !PRINT *, ' Local solar time in seconds:', isd_sol
            !PRINT *, '     UTC     time in hours:',   REAL(isd    ,wp)/3600._wp
            !PRINT *, '  [WL_COARE] Local solar time in hours:',   REAL(isd_sol,wp)/3600._wp

            !*****  variables for warm layer  ***
            zalpha_w = alpha_sw( pSST(ji,jj) ) ! thermal expansion coefficient of sea-water (SST accurate enough!)

            zcd1 = SQRT(2._wp*rich*rCp0_w/(zalpha_w*grav*rho0_w))        !mess-o-constants 1
            zcd2 = SQRT(2._wp*zalpha_w*grav/(rich*rho0_w))/(rCp0_w**1.5) !mess-o-constants 2

            !IF (isd_sol <= rdt ) THEN    !re-zero at midnight ! LOLO improve: risky if real midnight (00:00:00) is not a time in vtime...
            IF ( (rhr_sol > 23.5_wp).OR.(rhr_sol < 4._wp) ) THEN
               !PRINT *, '  [WL_COARE] MIDNIGHT RESET !!!!, isd_sol =>', isd_sol
               zdz           = H_wl_max
               Tau_ac(ji,jj) = 0._wp
               Qnt_ac(ji,jj) = 0._wp
            END IF

            IF ( rhr_sol > 5._wp ) THEN  ! ( 5am)

               !PRINT *, '  [WL_COARE] WE DO WL !'
               !PRINT *, '  [WL_COARE] isd_sol, pTau, pSST, pdT =', isd_sol, REAL(pTau(ji,jj),4), REAL(pSST(ji,jj),4), REAL(pdT(ji,jj),4)


               zQabs = zfr*pQsw(ji,jj) + pQnsol(ji,jj)       ! tot. heat flux absorbed in warm layer

               !PRINT *, '  [WL_COARE] rdt,  pQsw, pQnsol, zQabs =', rdt, REAL(pQsw(ji,jj),4), REAL(pQnsol(ji,jj),4), REAL(zQabs,4)

               !IF ( zQabs >= Qabs_thr ) THEN         ! Check for threshold

               ztac = Tau_ac(ji,jj) + MAX(.002_wp , pTau(ji,jj))*rdt      ! updated momentum integral

               IF ( Qnt_ac(ji,jj) + zQabs*rdt > 0._wp ) THEN         !check threshold for warm layer existence
                  !******************************************
                  ! Compute the absorption profile
                  !******************************************
                  DO jl = 1, 5                           !loop 5 times for zfr
                     zfr = 1. - ( 0.28*0.014*(1. - EXP(-zdz/0.014)) + 0.27*0.357*(1. - EXP(-zdz/0.357)) &
                        &        + 0.45*12.82*(1-EXP(-zdz/12.82)) ) / zdz
                     zqac = Qnt_ac(ji,jj) + (zfr*pQsw(ji,jj) + pQnsol(ji,jj))*rdt ! updated heat absorbed
                     IF ( zqac > 1._wp ) zdz = MAX( MIN( H_wl_max , zcd1*ztac/SQRT(zqac)) , 0.1_wp ) ! Warm-layer depth
                  END DO

               ELSE
                  !***********************
                  ! Warm layer wiped out
                  !***********************
                  zfr  = 0.75
                  zdz  = H_wl_max
                  zqac = Qnt_ac(ji,jj) + (zfr*pQsw(ji,jj) + pQnsol(ji,jj))*rdt ! updated heat absorbed

               END IF !IF ( Qnt_ac(ji,jj) + zQabs*rdt > 0._wp )

               IF ( zqac > 1._wp ) zdT_wl = zcd2*zqac**1.5/ztac * MAX(zqac/ABS(zqac),0._wp)  !! => IF(zqac>0._wp): zdT_wl=zcd2*zqac**1.5/ztac ; ELSE: zdT_wl=0. / ! normally: zqac > 0 !

               ! Warm layer correction
               flg = 0.5_wp + SIGN( 0.5_wp , gdept_1d(1)-zdz )               ! => 1 when gdept_1d(1)>zdz (pdT(ji,jj) = zdT_wl) | 0 when gdept_1d(1)<zdz (pdT(ji,jj) = zdT_wl*gdept_1d(1)/zdz)
               pdT(ji,jj) = zdT_wl * ( flg + (1._wp-flg)*gdept_1d(1)/zdz )

               !END IF ! IF ( zQabs >= Qabs_thr )

            END IF ! IF ( isd_sol >= 21600 ) THEN  ! (21600 == 6am)

            IF ( iwait == 0 ) THEN
               !IF ( (zQabs >= Qabs_thr).AND.(rhr_sol >= 5._wp) ) THEN
               IF ( rhr_sol >= 5._wp ) THEN
                  !PRINT *, '  [WL_COARE] WE UPDATE ACCUMULATED FLUXES !!!'
                  Qnt_ac(ji,jj) = zqac ! Updating Qnt_ac, heat integral
                  Tau_ac(ji,jj) = ztac !
                  IF ( PRESENT(mask_wl) ) mask_wl(ji,jj) = 1
               END IF
            END IF

            H_wl(ji,jj) = zdz

            IF ( PRESENT(Hwl) ) Hwl(ji,jj) = H_wl(ji,jj)

         END DO
      END DO

   END SUBROUTINE WL_COARE


   !!======================================================================
END MODULE mod_skin_coare
