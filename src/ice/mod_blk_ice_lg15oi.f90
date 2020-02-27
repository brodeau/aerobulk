! AeroBulk / 2020 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
! NEMO actually only needs the flux over sea-ice: F_ice
!  * F_ice is then used by SI3 which returns the basal ice-ocean flux F_b.
!  * As its SBC "Fs", OPA then uses Fs = (1.-ifr)*F_oce + ifr*F_b
!
! So we do not need the present routine to return a flux over open water
! At least when ifr is rather small...
!
!
!
MODULE mod_blk_ice_lg15oi
   !!====================================================================================
   !!       Computes turbulent components of surface fluxes over sea-ice
   !!       Following Lupkes & Gryanik, 2015
   !!
   !!              => case when 100 % sea-ice
   !!
   !!       Routine turb_ice_lg15oi maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, January 2020
   !!
   !!====================================================================================
   USE mod_const       !: physical and other constants
   USE mod_phymbl      !: misc. physical functions

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: turb_ice_lg15oi

   REAL(wp), PARAMETER ::   rce10_i_0 = 3.46e-3_wp ! (Eq.48) MIZ
   REAL(wp), PARAMETER ::   rbeta_0   = 1.4_wp     ! (Eq.47) MIZ
   REAL(wp), PARAMETER ::   ralpha_0  = 0.2_wp     ! (Eq.12) (ECHAM6 value)

   !! To be namelist parameters in NEMO:
   REAL(wp), PARAMETER ::   rz0_s_0  = 0.69e-3_wp  ! Eq. 43 [m]
   REAL(wp), PARAMETER ::   rz0_f_0  = 4.54e-4_wp  ! bottom p.562 MIZ [m]
   LOGICAL,  PARAMETER :: l_use_form_drag = .TRUE.
   LOGICAL,  PARAMETER :: l_use_pond_info = .FALSE.
   LOGICAL,  PARAMETER :: l_dbg_print     = .FALSE.


   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE turb_ice_lg15oi( kt, zt, zu, Ts_i, Ts_w, t_zt, qs_i, qs_w, q_zt, U_zu, frice, &
      &                        Cd_i, Ch_i, Ce_i, t_zu_i, q_zu_i,                            &
      &                        Cd_w, Ch_w, Ce_w, t_zu_w, q_zu_w,                            &
      &                        xz0, xu_star, xL, xUN10 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_ice_lg15oi  ***
      !!
      !! ** Purpose :   Computestransfert coefficients of turbulent surface
      !!                fluxes according
      !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!                Returns the effective bulk wind speed at zu to be used in the bulk formulas
      !!
      !! INPUT :
      !! -------
      !!    *  kt   : current time step (starts at 1)
      !!    *  zt   : height for temperature and spec. hum. of air            [m]
      !!    *  zu   : height for wind speed (usually 10m)                     [m]
      !!    *  Ts_i  : surface temperature of sea-ice                         [K]
      !!    *  Ts_w  : surface temperature of water (sea)                     [K]
      !!    *  t_zt : potential air temperature at zt                         [K]
      !!    *  qs_i  : saturation specific humidity at temp. Ts_i over ice    [kg/kg]
      !!    *  qs_w  : saturation specific humidity at temp. Ts_w over water  [kg/kg]
      !!    *  q_zt : specific humidity of air at zt                          [kg/kg]
      !!    *  U_zu : scalar wind speed at zu                                 [m/s]
      !!    * frice : sea-ice concentration        (fraction)
      !!
      !! OUTPUT :
      !! --------
      !!    *  Cd_i   : drag coefficient over sea-ice
      !!    *  Ch_i   : sensible heat coefficient over sea-ice
      !!    *  Ce_i   : sublimation coefficient over sea-ice
      !!    *  t_zu_i : pot. air temp. adjusted at zu over sea-ice             [K]
      !!    *  q_zu_i : spec. hum. of air adjusted at zu over sea-ice          [kg/kg]
      !!
      !! OPTIONAL OUTPUT:
      !! ----------------
      !!    *  Cd_w   : drag coefficient over water
      !!    *  Ch_w   : sensible heat coefficient over water
      !!    *  Ce_w   : sublimation coefficient over water
      !!    *  t_zu_w : pot. air temp. adjusted at zu over water             [K]
      !!    *  q_zu_w : spec. hum. of air adjusted at zu over water          [kg/kg]
      !!    * xz0     : return the aerodynamic roughness length (integration constant for wind stress) [m]
      !!    * xu_star : return u* the friction velocity                    [m/s]
      !!    * xL      : return the Obukhov length                          [m]
      !!    * xUN10   : neutral wind speed at 10m                          [m/s]
      !!
      !! ** Author: L. Brodeau, January 2020 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      INTEGER,  INTENT(in )                     :: kt    ! current time step
      REAL(wp), INTENT(in )                     :: zt    ! height for t_zt and q_zt                    [m]
      REAL(wp), INTENT(in )                     :: zu    ! height for U_zu                             [m]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: Ts_i  ! ice surface temperature                [Kelvin]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: Ts_w  ! water surface temperature              [Kelvin]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: t_zt  ! potential air temperature              [Kelvin]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: qs_i  ! sat. spec. hum. at ice/air interface    [kg/kg]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: qs_w  ! sat. spec. hum. at water/air interface  [kg/kg]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: q_zt  ! spec. air humidity at zt               [kg/kg]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: U_zu  ! relative wind module at zu                [m/s]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: frice ! sea-ice concentration        (fraction)
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Cd_i  ! drag coefficient over sea-ice
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Ch_i  ! transfert coefficient for heat over ice
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Ce_i  ! transfert coefficient for sublimation over ice
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: t_zu_i  ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: q_zu_i  ! spec. humidity adjusted at zu           [kg/kg]
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(jpi,jpj) :: Cd_w    ! drag coefficient over sea-ice
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(jpi,jpj) :: Ch_w    ! transfert coefficient for heat over ice
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(jpi,jpj) :: Ce_w    ! transfert coefficient for sublimation over ice
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(jpi,jpj) :: t_zu_w  ! pot. air temp. adjusted at zu over water    [K]
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(jpi,jpj) :: q_zu_w  ! spec. humidity adjusted at zu over water [kg/kg]
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(jpi,jpj) :: xz0     ! Aerodynamic roughness length                [m]
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(jpi,jpj) :: xu_star ! u*, friction velocity
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(jpi,jpj) :: xL      ! Obukhov length                              [m]
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(jpi,jpj) :: xUN10   ! Neutral wind at zu
      !!
      INTEGER :: j_itt
      LOGICAL :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !!
      REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: Ub, xtmp1, xtmp2        !! bulk wind speed at zu + temporary arrays
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: dt_zu, dq_zu, zt_zu, zq_zu  ! third dimension
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zz0_s, zz0_f, RiB ! third dimensions (size=2):
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zCd, zCh, zCdN_s, zChN_s, zCdN_f, zChN_f

      LOGICAL :: lreturn_z0=.FALSE., lreturn_ustar=.FALSE., lreturn_L=.FALSE., lreturn_UN10=.FALSE.
      LOGICAL :: lreturn_o_water=.FALSE.
      CHARACTER(len=40), PARAMETER :: crtnm = 'turb_ice_lg15oi@mod_blk_ice_lg15oi.f90'
      !!----------------------------------------------------------------------------------
      ALLOCATE ( Ub(jpi,jpj), xtmp1(jpi,jpj), xtmp2(jpi,jpj) )
      ALLOCATE ( dt_zu(jpi,jpj,2), dq_zu(jpi,jpj,2), zt_zu(jpi,jpj,2), zq_zu(jpi,jpj,2) )
      ALLOCATE ( zz0_s(jpi,jpj,2),  zz0_f(jpi,jpj,2),    RiB(jpi,jpj,2), &
         &      zCdN_s(jpi,jpj,2), zChN_s(jpi,jpj,2), zCdN_f(jpi,jpj,2), zChN_f(jpi,jpj,2) )
      ALLOCATE ( zCd(jpi,jpj,2), zCh(jpi,jpj,2) )

      lreturn_o_water =  PRESENT(Cd_w) .AND. PRESENT(Ch_w) .AND. PRESENT(Ce_w) .AND. PRESENT(t_zu_w) .AND. PRESENT(q_zu_w) 

      IF( PRESENT(Cd_w) ) lreturn_o_water = .TRUE.
      IF( PRESENT(xz0) )     lreturn_z0    = .TRUE.
      IF( PRESENT(xu_star) ) lreturn_ustar = .TRUE.
      IF( PRESENT(xL) )      lreturn_L     = .TRUE.
      IF( PRESENT(xUN10) )   lreturn_UN10  = .TRUE.

      l_zt_equal_zu = .FALSE.
      IF( ABS(zu - zt) < 0.01_wp )   l_zt_equal_zu = .TRUE. ! testing "zu == zt" is risky with double precision

      !! Scalar wind speed cannot be below 0.2 m/s
      Ub = MAX( U_zu, wspd_thrshld_ice )

      !! First guess of temperature and humidity at height zu:
      zt_zu(:,:,1) = MAX( t_zt(:,:) ,   100._wp )   ! who knows what's given on masked-continental regions...
      zq_zu(:,:,1) = MAX( q_zt(:,:) , 0.1e-6_wp )   !               "
      IF( lreturn_o_water ) THEN
         zt_zu(:,:,2) = MAX( t_zt(:,:) ,   100._wp )   ! who knows what's given on masked-continental regions...
         zq_zu(:,:,2) = MAX( q_zt(:,:) , 0.1e-6_wp )   !               "
      END IF

      !! Air-Ice & Air-Sea differences (and we don't want them to be 0!)
      dt_zu(:,:,1) = zt_zu(:,:,1) - Ts_i
      dq_zu(:,:,1) = zq_zu(:,:,1) - qs_i
      IF( lreturn_o_water ) THEN
         dt_zu(:,:,2) = zt_zu(:,:,2) - Ts_w
         dq_zu(:,:,2) = zq_zu(:,:,2) - qs_w
      END IF
      dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6_wp), dt_zu )
      dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9_wp), dq_zu )

      !! Very crude first guess:
      Cd_i(:,:) = rCd_ice
      Ch_i(:,:) = rCd_ice
      Ce_i(:,:) = rCd_ice
      IF( lreturn_o_water ) THEN
         Cd_w(:,:) = 0.001_wp
         Ch_w(:,:) = 0.001_wp
         Ce_w(:,:) = 0.001_wp
      END IF

      !! For skin drag :
      zz0_s(:,:,1) = rz0_s_0        !LOLO/RFI! ! Room for improvement. We use the same z0_skin everywhere (= rz0_s_0)...
      xtmp1(:,:) = LOG( zu / zz0_s(:,:,1) )
      zCdN_s(:,:,1) = vkarmn2 / ( xtmp1(:,:) * xtmp1(:,:) )                          ! (Eq.7)   [ index 1 is for ice, 2 for water ]
      zChN_s(:,:,1) = vkarmn2 / ( xtmp1(:,:) * LOG( zu / (ralpha_0*zz0_s(:,:,1)) ) )     ! (Eq.11,12)  [ "" ]

      !! For form drag in MIZ:
      zz0_f(:,:,:)  = 0._wp
      zCdN_f(:,:,:) = 0._wp
      zChN_f(:,:,:) = 0._wp
      IF ( l_use_form_drag ) THEN
         zz0_f(:,:,1) = rz0_f_0        !LOLO/RFI! ! Room for improvement. We use the same z0_form everywhere !!!
         xtmp1(:,:) = 1._wp / zz0_f(:,:,1)
         xtmp2(:,:) = rce10_i_0 * ( LOG( 10._wp * xtmp1(:,:) ) / LOG( zu * xtmp1(:,:) ) )**2      ! part of (Eq.46)
         zCdN_f(:,:,1) = xtmp2(:,:) * frice(:,:) * (1._wp - frice(:,:))**rbeta_0                  ! (Eq.46)  [ index 1 is for ice, 2 for water ]
         zChN_f(:,:,1) = zCdN_f(:,:,1)/( 1._wp + LOG(1._wp/ralpha_0)/vkarmn*SQRT(zCdN_f(:,:,1)) ) ! (Eq.60,61)   [ "" ]
      END IF

      !! Some other first guess values, needed to compute wind at zt:
      zCd(:,:,1) = zCdN_s(:,:,1) + zCdN_f(:,:,1)
      zCh(:,:,1) = zChN_s(:,:,1) + zChN_f(:,:,1)
      RiB(:,:,1) = Ri_bulk( zt, Ts_i(:,:), t_zt(:,:), qs_i(:,:), q_zt(:,:), Ub(:,:) )  ! over ice (index=1)


      !! ITERATION BLOCK
      DO j_itt = 1, nb_itt

         IF(l_dbg_print) PRINT *, 'LOLO: LOOP #', INT(j_itt,1)
         IF(l_dbg_print) PRINT *, 'LOLO: theta_zu, Ts_i, Ub =', REAL(zt_zu(:,:,1),4), REAL(Ts_i(:,:),4), REAL(Ub(:,:),4)
         IF(l_dbg_print) PRINT *, 'LOLO:     q_zu =', REAL(zq_zu(:,:,1),4)
         IF(l_dbg_print) PRINT *, 'LOLO:  CdN_s, zCdN_f   =', REAL(zCdN_s(:,:,1),4), REAL(zCdN_f(:,:,1),4)


         !! Bulk Richardson Number
         !! ======================
         !! PROBLEM: when computed at z=zu, with adjusted theta and q, it is numerically unstable in some rare events (unstable)
         !!          => fix: compute RiB at zt, with ajusted wind at zt... => seems to be more stable
         IF( .NOT. l_zt_equal_zu ) THEN
            ! U_zt = U_zu + u_star/vkarmn*(LOG(zt/zu) + psi_m_coare(zu/L) - psi_m_coare(zt/L))
            xtmp1(:,:) = zCdN_s(:,:,1) + zCdN_f(:,:,1)    ! total neutral drag coeff!
            xtmp2(:,:) = zz0_s(:,:,1) + zz0_f(:,:,1)      ! total roughness length z0
            xtmp1 = LOG(zt/zu) + f_h_louis( zu, RiB(:,:,1), xtmp1(:,:), xtmp2(:,:) ) &
               &               - f_h_louis( zt, RiB(:,:,1), xtmp1(:,:), xtmp2(:,:) )
            xtmp2(:,:) = MAX( Ub(:,:) + (SQRT(zCd(:,:,1))*Ub)*xtmp1 , wspd_thrshld_ice ) ! wind at zt ( SQRT(zCd(:,:,1))*Ub == u* !)
            xtmp2(:,:) = MIN( xtmp2(:,:) , Ub(:,:) )
            IF(l_dbg_print) PRINT *, 'LOLO: ADJUSTED WIND AT ZT =', xtmp2
         ELSE
            xtmp2(:,:) = Ub(:,:)
         END IF
         RiB(:,:,1) = Ri_bulk( zt, Ts_i(:,:), t_zt(:,:), qs_i(:,:), q_zt(:,:), xtmp2(:,:) )  ! over ice (index=1)
         IF(l_dbg_print) PRINT *, 'LOLO: RiB_zt =', RiB(:,:,1)

         !RiB(:,:,1) = Ri_bulk( zu, Ts_i(:,:), zt_zu(:,:,1), qs_i(:,:), zq_zu(:,:,1), Ub(:,:) )
         !IF(l_dbg_print) PRINT *, 'LOLO: RiB_zu =', RiB(:,:,1)

         IF( lreturn_o_water ) THEN
            RiB(:,:,2) = Ri_bulk( zu, Ts_w(:,:), zt_zu(:,:,2), qs_w(:,:), zq_zu(:,:,2), Ub(:,:) )  ! over water (index=2)
            IF(l_dbg_print) PRINT *, 'LOLO: over water RiB_zt =', RiB(:,:,2)
         END IF

         ! Momentum and Heat transfer coefficients WITHOUT FORM DRAG / (Eq.6) and (Eq.10):
         zCd(:,:,1) = zCdN_s(:,:,1) * f_m_louis( zu, RiB(:,:,1), zCdN_s(:,:,1), zz0_s(:,:,1) ) ! (Eq.6)
         zCh(:,:,1) = zChN_s(:,:,1) * f_h_louis( zu, RiB(:,:,1), zCdN_s(:,:,1), zz0_s(:,:,1) ) ! (Eq.10) / LOLO: why "zCdN_s" (xtmp1) and not "zChn" ???
         IF(l_dbg_print) PRINT *, 'LOLO: f_m_louis_s =', f_m_louis( zu, RiB(:,:,1), zCdN_s(:,:,1), zz0_s(:,:,1) )
         IF(l_dbg_print) PRINT *, 'LOLO: f_h_louis_s =', f_h_louis( zu, RiB(:,:,1), zCdN_s(:,:,1), zz0_s(:,:,1) )
         IF(l_dbg_print) PRINT *, 'LOLO: Cd / skin only / ice   =', REAL(zCd(:,:,1),4)

         IF( lreturn_o_water ) THEN
            zCd(:,:,2) = zCdN_s(:,:,2) * f_m_louis( zu, RiB(:,:,2), zCdN_s(:,:,2), zz0_s(:,:,2) ) ! (Eq.6)
            zCh(:,:,2) = zChN_s(:,:,2) * f_h_louis( zu, RiB(:,:,2), zCdN_s(:,:,2), zz0_s(:,:,2) ) ! (Eq.10) / LOLO: why "zCdN_s" (xtmp1) and not "zChn" ???
            IF(l_dbg_print) PRINT *, 'LOLO: Cd / skin only / water =', REAL(zCd(:,:,2),4)
         END IF


         IF ( l_use_form_drag ) THEN
            !! Form-drag-related NEUTRAL momentum and Heat transfer coefficients:
            !!   MIZ:
            zCd(:,:,1) = zCd(:,:,1) + zCdN_f(:,:,1) * f_m_louis( zu, RiB(:,:,1), zCdN_f(:,:,1), zz0_f(:,:,1) ) ! (Eq.6)
            zCh(:,:,1) = zCh(:,:,1) + zChN_f(:,:,1) * f_h_louis( zu, RiB(:,:,1), zCdN_f(:,:,1), zz0_f(:,:,1) ) ! (Eq.10) / LOLO: why "zCdN_f" and not "zChn" ???
            IF(l_dbg_print) PRINT *, 'LOLO: f_m_louis_f =', f_m_louis( zu, RiB(:,:,1), zCdN_f(:,:,1), zz0_f(:,:,1) )
            IF(l_dbg_print) PRINT *, 'LOLO: f_h_louis_f =', f_h_louis( zu, RiB(:,:,1), zCdN_f(:,:,1), zz0_f(:,:,1) )

            IF(l_dbg_print) PRINT *, 'LOLO: Cd / form only / ice   =', REAL(zCdN_f(:,:,1) * f_m_louis( zu, RiB(:,:,1), zCdN_f(:,:,1), zz0_f(:,:,1) ),4)
            !zCd(:,:,2) = ???
            !zCh(:,:,2) = ???

         END IF

         IF(l_dbg_print) PRINT *, 'LOLO: Cd, Ch / TOTAL / ice   =', REAL(zCd(:,:,1),4), REAL(zCh(:,:,1),4)


         !! Adjusting temperature and humidity from zt to zu:
         IF( .NOT. l_zt_equal_zu ) THEN

            !! Over ice:
            xtmp1(:,:) = zCdN_s(:,:,1) + zCdN_f(:,:,1)    ! total neutral drag coeff!
            xtmp2(:,:) = zz0_s(:,:,1) + zz0_f(:,:,1)      ! total roughness length z0
            xtmp1 = LOG(zt/zu) + f_h_louis( zu, RiB(:,:,1), xtmp1(:,:), xtmp2(:,:) ) &
               &               - f_h_louis( zt, RiB(:,:,1), xtmp1(:,:), xtmp2(:,:) )
            xtmp2 = 1._wp/SQRT(zCd(:,:,1))

            zt_zu(:,:,1) = t_zt - (zCh(:,:,1) * dt_zu(:,:,1) * xtmp2) / vkarmn * xtmp1   ! t_star = Ch * dt_zu / SQRT(Cd)
            zq_zu(:,:,1) = q_zt - (zCh(:,:,1) * dq_zu(:,:,1) * xtmp2) / vkarmn * xtmp1   ! q_star = Ce * dq_zu / SQRT(Cd)
            zq_zu(:,:,1) = MAX(0._wp, zq_zu(:,:,1))

            dt_zu(:,:,1) = zt_zu(:,:,1) - Ts_i
            dq_zu(:,:,1) = zq_zu(:,:,1) - qs_i

            !! Over water:
            IF( lreturn_o_water ) THEN
               xtmp1(:,:) = zCdN_s(:,:,2) + zCdN_f(:,:,2)    ! total neutral drag coeff!
               xtmp2(:,:) = zz0_s(:,:,2) + zz0_f(:,:,2)      ! total roughness length z0
               xtmp1 = LOG(zt/zu) + f_h_louis( zu, RiB(:,:,2), xtmp1(:,:), xtmp2(:,:) ) &
                  &               - f_h_louis( zt, RiB(:,:,2), xtmp1(:,:), xtmp2(:,:) )
               xtmp2 = 1._wp/SQRT(zCd(:,:,2))
               zt_zu(:,:,2) = t_zt - (zCh(:,:,2) * dt_zu(:,:,2) * xtmp2) / vkarmn * xtmp1   ! t_star = Ch * dt_zu / SQRT(Cd)
               zq_zu(:,:,2) = q_zt - (zCh(:,:,2) * dq_zu(:,:,2) * xtmp2) / vkarmn * xtmp1   ! q_star = Ce * dq_zu / SQRT(Cd)
               zq_zu(:,:,2) = MAX(0._wp, zq_zu(:,:,2))
               dt_zu(:,:,2) = zt_zu(:,:,2) - Ts_w
               dq_zu(:,:,2) = zq_zu(:,:,2) - qs_w
            END IF
            
            dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6_wp), dt_zu )
            dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9_wp), dq_zu )
         END IF

         IF(l_dbg_print) PRINT *, ''!LOLO

      END DO !DO j_itt = 1, nb_itt
      IF(l_dbg_print) PRINT *, ''!LOLO

      !! Result is ice + ocean:
      !t_zu(:,:) = mix_val_msh(zt_zu, frice)
      !q_zu(:,:) = mix_val_msh(zq_zu, frice)
      !Cd(:,:) = mix_val_msh(zCd, frice)
      !Ch(:,:) = mix_val_msh(zCh, frice)

      !! Result is over ice only:
      t_zu_i(:,:) = zt_zu(:,:,1)
      q_zu_i(:,:) = zq_zu(:,:,1)
      Cd_i(:,:)   =   zCd(:,:,1)
      Ch_i(:,:)   =   zCh(:,:,1)
      Ce_i(:,:)   =  Ch_i(:,:)

      IF( lreturn_o_water ) THEN
         t_zu_w(:,:) = zt_zu(:,:,2)
         q_zu_w(:,:) = zq_zu(:,:,2)
         Cd_w(:,:)   =   zCd(:,:,2)
         Ch_w(:,:)   =   zCh(:,:,2)
         Ce_w(:,:)   =  Ch_w(:,:)
      END IF
               
      IF( lreturn_z0 ) xz0   = z0_from_Cd( zu, zCdN_s(:,:,1)+zCdN_f(:,:,1) )

      IF( lreturn_ustar ) xu_star = SQRT(Cd_i) * Ub
      IF( lreturn_L ) THEN
         xtmp1 = SQRT(Cd_i)
         !xL    = 1./One_on_L(t_zu_i, q_zu, xtmp1*Ub, Ch_i*mix_val_msh(dt_zu,frice)/xtmp1, Ce_i*mix_val_msh(dq_zu,frice)/xtmp1)
         xL    = 1./One_on_L( t_zu_i, q_zu_i, xtmp1*Ub, Ch_i*dt_zu(:,:,1)/xtmp1, Ce_i*dq_zu(:,:,1)/xtmp1 )
      END IF

      IF( lreturn_UN10 ) THEN
         xtmp1 = zCdN_s(:,:,1) + zCdN_f(:,:,1)  ! => CdN
         xUN10 = SQRT(Cd_i) * Ub/vkarmn * LOG( 10._wp / z0_from_Cd(zu, xtmp1) )
         !xtmp2 = f_m_louis( zu, RiB(:,:,1), xtmp1, z0_from_Cd(zu, xtmp1) ) ! => f_m
         !xUN10 = UN10_from_CD( zu, Ub, Cd_i, ppsi=xtmp2 )
      END IF



      DEALLOCATE ( Ub, xtmp1, xtmp2 )
      DEALLOCATE ( dt_zu, dq_zu, zt_zu, zq_zu )
      DEALLOCATE ( zz0_s, zz0_f, RiB, zCdN_s, zChN_s, zCdN_f, zChN_f )
      DEALLOCATE ( zCd, zCh )

   END SUBROUTINE turb_ice_lg15oi

   !!======================================================================



   FUNCTION mix_val_msh( pfld, pfri )
      REAL(wp), DIMENSION(jpi,jpj) :: mix_val_msh
      REAL(wp), DIMENSION(jpi,jpj,2), INTENT(in) :: pfld  ! field array to "water/ice average" on the mesh [ over ice=>(:,:,1), over water=>(:,:,2) ]
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in) :: pfri  ! sea-ice concentration (fraction)
      mix_val_msh(:,:) = pfri(:,:) * pfld(:,:,1) + (1._wp - pfri(:,:)) * pfld(:,:,2)
   END FUNCTION mix_val_msh

END MODULE mod_blk_ice_lg15oi
