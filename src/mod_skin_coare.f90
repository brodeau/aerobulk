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
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC :: delta_vl  !: thickness of the surface viscous layer (in the water) right below the air-sea interface [m]

   !! Warm-layer related parameters:
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC :: &
      &                        dT_wl,     &  ! dT due to warm-layer effect => difference between "almost surface (right below viscous layer) and depth of bulk SST
      &                        Hz_wl,      &  ! depth of warm-layer [m]
      &                        Qnt_ac,    &  ! time integral / accumulated heat stored by the warm layer Qxdt => [J/m^2] (reset to zero every midnight)
      &                        Tau_ac        ! time integral / accumulated momentum Tauxdt => [N.s/m^2] (reset to zero every midnight)

   REAL(wp), PARAMETER, PUBLIC :: Hwl_max = 20._wp    !: maximum depth of warm layer (adjustable)
   !
   REAL(wp), PARAMETER :: rich   = 0.65_wp   !: critical Richardson number
   !
   !REAL(wp), PARAMETER :: Qabs_thr = 50._wp  !: threshold for heat flux absorbed in WL !DO NOT GET IT !!!!

   REAL(wp), PARAMETER :: zfr0   = 0.5_wp    !: initial value of solar flux absorption
   !
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
      REAL(wp) :: zQabs, zQnsol, zlamb, zdelta, zalpha_w, zfr, &
         &        zz1, zz2, zusw, &
         &        ztf
      !!---------------------------------------------------------------------

      DO jj = 1, jpj
         DO ji = 1, jpi

            zalpha_w = alpha_sw( pSST(ji,jj) ) ! thermal expansion coefficient of sea-water (SST accurate enough!)

            zQnsol = MAX( 1._wp , - pQnsol(ji,jj) ) ! Non-solar heat loss to the atmosphere

            zdelta = delta_vl(ji,jj)   ! using last value of delta

            !! Fraction of the shortwave flux absorbed by the cool-skin sublayer:
            zfr   = MAX( 0.137_wp + 11._wp*zdelta - 6.6E-5_wp/zdelta*(1._wp - EXP(-zdelta/8.E-4_wp)) , 0.01_wp ) ! Eq.16 (Fairall al. 1996b) /  !LB: why 0.065 and not 0.137 like in the paper??? Beljaars & Zeng use 0.065, not 0.137 !

            zQabs = MAX( 1._wp , zQnsol - zfr*pQsw(ji,jj) )  ! Total cooling at the interface

            ztf = 0.5 + SIGN(0.5_wp, zQabs) ! Qt > 0 => cooling of the layer => ztf = 1
            !                                 Qt < 0 => warming of the layer => ztf = 0

            !! Term alpha*Qb (Qb is the virtual surface cooling inc. buoyancy effect of salinity due to evap):
            zz1 = zalpha_w*zQabs - 0.026*MIN(pQlat(ji,jj),0._wp)*rCp0_w/rLevap  ! alpha*(Eq.8) == alpha*Qb "-" because Qlat < 0
            !! LB: this terms only makes sense if > 0 i.e. in the cooling case
            !! so similar to what's done in ECMWF:
            zz1 = MAX(0._wp , zz1)    ! 1. instead of 0.1 though ZQ = MAX(1.0,-pQlw(ji,jj) - pQsen(ji,jj) - pQlat(ji,jj))

            zusw = MAX(pustar(ji,jj), 1.E-4_wp) * sq_radrw   ! u* in water / ! Laurent: too low wind (u*) might cause problem in stable cases:
            zz2 = zusw*zusw
            zz2 = zz2*zz2
            
            zlamb =  6._wp*( 1._wp + (rcst_cs*zz1/zz2)**0.75 )**(-1./3.) ! Lambda (Eq.14) (Saunders)

            ! Updating molecular sublayer thickness (delta):
            zz2    = rnu0_w/zusw
            zdelta =      ztf    *          zlamb*zz2   &  ! Eq.12 (when alpha*Qb>0 / cooling of layer)
               &    + (1._wp - ztf) * MIN(0.007_wp , 6._wp*zz2 )    ! Eq.12 (when alpha*Qb<0 / warming of layer)
            !LB: changed 0.01 to 0.007

            !! Once again with the new zdelta:
            zfr   = MAX( 0.137_wp + 11._wp*zdelta - 6.6E-5_wp/zdelta*(1._wp - EXP(-zdelta/8.E-4_wp)) , 0.01_wp ) ! Solar absorption / Eq.16 (Fairall al. 1996b)
            zQabs = MAX( 1._wp , zQnsol - zfr*pQsw(ji,jj) ) ! Total cooling at the interface

            !! Update!
            pdT(ji,jj) =  MIN( - zQabs*zdelta/rk0_w , 0._wp )   ! temperature increment
            delta_vl(ji,jj) = zdelta

         END DO
      END DO

   END SUBROUTINE CS_COARE


   SUBROUTINE WL_COARE( pQsw, pQnsol, pTau, pSST, plon, isd, iwait )
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
      !!---------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pQsw     ! surface net solar radiation into the ocean [W/m^2]     => >= 0 !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pQnsol   ! surface net non-solar heat flux into the ocean [W/m^2] => normally < 0 !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pTau     ! wind stress [N/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pSST     ! bulk SST at depth gdept_1d(1) [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: plon     ! longitude [deg.E]
      INTEGER ,                     INTENT(in)    :: isd      ! current UTC time, counted in second since 00h of the current day
      INTEGER ,                     INTENT(in  )  :: iwait    ! if /= 0 then wait before updating accumulated fluxes
      !!
      INTEGER :: ji,jj
      !
      REAL(wp) :: zdTwl, zHwl, zQabs, zfr
      REAL(wp) :: zqac, ztac
      REAL(wp) :: zalpha_w, zcd1, zcd2, flg
      !!---------------------------------------------------------------------

      REAL(wp) :: rlag_gw_h, &  ! local solar time lag in hours   / Greenwich meridian (lon==0) => ex: ~ -10.47 hours for Hawai
         &        rhr_sol       ! local solar time in hours since midnight

      INTEGER  :: ilag_gw_s, &  ! local solar time LAG in seconds / Greenwich meridian (lon==0) => ex: ~ INT( -10.47*3600. ) seconds for Hawai
         &        isd_sol,   &  ! local solar time in number of seconds since local solar midnight
         &        jl

      LOGICAL :: l_exit, l_destroy_wl

      !! INITIALIZATION:
      zQabs  = 0._wp       ! total heat flux absorped in warm layer
      zfr    = zfr0        ! initial value of solar flux absorption !LOLO: save it and use previous value !!!

      DO jj = 1, jpj
         DO ji = 1, jpi

            l_exit       = .FALSE.
            l_destroy_wl = .FALSE.
            
            zdTwl =  dT_wl(ji,jj)                          ! value of previous time step as first guess
            zHwl  = MAX( MIN(Hz_wl(ji,jj),Hwl_max),0.1_wp) !   "                  "           "
            
            zqac = Qnt_ac(ji,jj) ! previous time step Qnt_ac
            ztac = Tau_ac(ji,jj)

            !! Need to know the local solar time from longitude and isd:
            rlag_gw_h = -1._wp * MODULO( ( 360._wp - MODULO(plon(ji,jj),360._wp) ) / 15._wp , 24._wp )
            rlag_gw_h = -1._wp * SIGN( MIN(ABS(rlag_gw_h) , ABS(MODULO(rlag_gw_h,24._wp))), rlag_gw_h + 12._wp )
            ilag_gw_s = INT( rlag_gw_h*3600._wp )
            isd_sol = MODULO( isd + ilag_gw_s , 24*3600 )
            rhr_sol = REAL( isd_sol , wp) / 3600._wp

            PRINT *, '#LBD:'
            PRINT *, '#LBD: dT_wl at former ts:', zdTwl
            PRINT *, '#LBD:  Hz_wl at former ts:', zHwl
            
            !*****  variables for warm layer  ***
            zalpha_w = alpha_sw( pSST(ji,jj) ) ! thermal expansion coefficient of sea-water (SST accurate enough!)

            zcd1 = SQRT(2._wp*rich*rCp0_w/(zalpha_w*grav*rho0_w))        !mess-o-constants 1
            zcd2 = SQRT(2._wp*zalpha_w*grav/(rich*rho0_w))/(rCp0_w**1.5) !mess-o-constants 2


            !IF ( rhr_sol > 5._wp ) THEN  ! ( 5am)

            IF ( (rhr_sol > 4._wp).AND.(rhr_sol <= 6.5_wp) ) THEN
               ! Dawn reset to 0!
               l_exit       = .TRUE.
               l_destroy_wl = .TRUE.
            END IF
            
            IF ( .NOT. l_exit ) THEN
               !! Initial test on initial guess of absorbed heat flux in warm-layer:
               zfr = 1._wp - ( 0.28*0.014*(1. - EXP(-zHwl/0.014)) + 0.27*0.357*(1. - EXP(-zHwl/0.357)) &
                  &        + 0.45*12.82*(1-EXP(-zHwl/12.82)) ) / zHwl
               zQabs = zfr*pQsw(ji,jj) + pQnsol(ji,jj) ! first guess of tot. heat flux absorbed in warm layer !LOLO: depends of zfr, which is wild guess... Wrong!!!
               PRINT *, '#LBD:  Initial Qsw & Qnsol:', NINT(pQsw(ji,jj)), NINT(pQnsol(ji,jj))
               PRINT *, '#LBD:       =>Qabs:', zQabs,' zfr=', zfr
               
               IF ( (ABS(zdTwl) < 1.E-6_wp) .AND. (zQabs <= 0._wp) ) THEN
                  ! We have not started to build a WL yet (dT==0) and there's no way it can occur now
                  ! since zQabs <= 0._wp
                  ! => no need to go further
                  PRINT *, '#LBD: we have not started to to build a WL yet (dT==0)'
                  PRINT *, '#LBD: and theres no way it can occur now since zQabs=', zQabs
                  PRINT *, '#LBD: => leaving without changing anything...'
                  l_exit = .TRUE.
               END IF

            END IF

            ! Okay test on updated absorbed flux:
            !LOLO: remove??? has a strong influence !!!
            IF ( (.NOT.(l_exit)) .AND. (Qnt_ac(ji,jj) + zQabs*rdt <= 0._wp) ) THEN
               PRINT *, '#LBD: Oh boy! Next Qnt_ac looking weak! =>', Qnt_ac(ji,jj) + zQabs*rdt
               PRINT *, '#LBD:  => time to destroy the warm-layer!'
               l_exit       = .TRUE.
               l_destroy_wl = .TRUE.
            END IF
            

            IF ( .NOT. l_exit) THEN
               
               ! Two possibilities at this point:
               ! 1/ A warm layer already exists (dT>0) but it is cooling down because Qabs<0
               ! 2/ Regardless of WL formed (dT==0 or dT>0), we are in the process to initiate one or warm further it !

               PRINT *, '#LBD:======================================================'
               PRINT *, '#LBD: WL action makes sense now! => zQabs,dT_wl=', REAL(zQabs,4), REAL(zdTwl,4)
               PRINT *, '#LBD:======================================================'
               PRINT *, '#LBD: current values for Qac and Tac=', REAL(Qnt_ac(ji,jj),4), REAL(Tau_ac(ji,jj),4)

               ztac = Tau_ac(ji,jj) + MAX(.002_wp , pTau(ji,jj))*rdt      ! updated momentum integral
               PRINT *, '#LBD: updated value for Tac=',  REAL(ztac,4)

               !! We update the value of absorbtion and zQabs:
               !! some part is useless if Qsw=0 !!!
               DO jl = 1, 5
                  zfr = 1. - ( 0.28*0.014*(1. - EXP(-zHwl/0.014)) + 0.27*0.357*(1. - EXP(-zHwl/0.357)) &
                     &        + 0.45*12.82*(1-EXP(-zHwl/12.82)) ) / zHwl
                  zQabs = zfr*pQsw(ji,jj) + pQnsol(ji,jj)
                  zqac  = Qnt_ac(ji,jj) + zQabs*rdt ! updated heat absorbed
                  IF ( zqac <= 0._wp ) EXIT
                  zHwl = MAX( MIN( Hwl_max , zcd1*ztac/SQRT(zqac)) , 0.1_wp ) ! Warm-layer depth
               END DO
               PRINT *, '#LBD: updated absorption and WL depth=',  REAL(zfr,4), REAL(zHwl,4)
               PRINT *, '#LBD: updated value for Qabs=',  REAL(zQabs,4), 'W/m2'
               PRINT *, '#LBD: updated value for Qac =',  REAL(zqac,4), 'J'

               IF ( zqac <= 0._wp ) THEN
                  l_destroy_wl = .TRUE.
                  l_exit       = .TRUE.
               ELSE
                  zdTwl = zcd2*zqac**1.5/ztac * MAX(zqac/ABS(zqac),0._wp)  !! => IF(zqac>0._wp): zdTwl=zcd2*zqac**1.5/ztac ; ELSE: zdTwl=0. / ! normally: zqac > 0 !
                  PRINT *, '#LBD: updated preliminary value for dT_wl=',  REAL(zdTwl,4)
                  ! Warm layer correction
                  flg = 0.5_wp + SIGN( 0.5_wp , gdept_1d(1)-zHwl )               ! => 1 when gdept_1d(1)>zHwl (zdTwl = zdTwl) | 0 when gdept_1d(1)<zHwl (zdTwl = zdTwl*gdept_1d(1)/zHwl)
                  zdTwl = zdTwl * ( flg + (1._wp-flg)*gdept_1d(1)/zHwl )
               END IF
               
            END IF !IF ( .NOT. l_exit)
            
            IF ( l_destroy_wl ) THEN
               zdTwl = 0._wp
               zfr   = 0.75_wp
               zHwl  = Hwl_max
               zqac  = 0._wp
               ztac  = 0._wp
            END IF
            
            PRINT *, '#LBD: exit values for Qac & Tac:', REAL(zqac,4), REAL(ztac,4)

            IF ( iwait == 0 ) THEN
               !! Iteration loop within bulk algo is over, time to update what needs to be updated:
               dT_wl(ji,jj)  = zdTwl
               Hz_wl(ji,jj)  = zHwl
               PRINT *, '#LBD: FINAL EXIT values for dT_wl & Hz_wl:', REAL(dT_wl(ji,jj),4), REAL(Hz_wl(ji,jj),4)
               Qnt_ac(ji,jj) = zqac ! Updating Qnt_ac, heat integral
               Tau_ac(ji,jj) = ztac
               PRINT *, '#LBD: FINAL EXIT values for Qac & Tac:', REAL(Qnt_ac(ji,jj),4), REAL(Tau_ac(ji,jj),4)
               PRINT *, '#LBD'
            END IF

         END DO
      END DO


      !STOP 'lolo1'
   END SUBROUTINE WL_COARE


   !!======================================================================
END MODULE mod_skin_coare
