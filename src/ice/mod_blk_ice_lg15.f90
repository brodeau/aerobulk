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
! In "ice.F90" of SI3 (NEMO), the following pond information is available (also the same per ice cat.):
!  at_ip      !: total melt pond concentration
!  hm_ip      !: mean melt pond depth                     [m]
!  vt_ip      !: total melt pond volume per gridcell area [m]
!
!
!
MODULE mod_blk_ice_lg15
   !!====================================================================================
   !!       Computes turbulent components of surface fluxes over sea-ice
   !!       Following Lupkes & Gryanik, 2015
   !!
   !!
   !!  Lüpkes, C., and Gryanik, V. M. ( 2015), A stability‐dependent parametrization
   !!  of transfer coefficients for momentum and heat over polar sea ice to be used in climate models,
   !!  J. Geophys. Res. Atmos., 120, 552– 581, doi:10.1002/2014JD022418.
   !!
   !!       => Sespite the fact that the sea-ice concentration (frice) must be provided,
   !!          only transfer coefficients, and air temp. + hum. height adjustement
   !!          over ice are returned/performed.
   !!        ==> 'frice' is only here to estimate the form drag caused by sea-ice...
   !!
   !!       Routine turb_ice_lg15 maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, January 2020
   !!
   !!====================================================================================
   USE mod_const       !: physical and other constants
   USE mod_phymbl      !: misc. physical functions
   USE mod_cdn_form_ice

   IMPLICIT NONE
   PRIVATE

   REAL(wp), PARAMETER ::   ralpha_0  = 0.2_wp     ! (Eq.12) (ECHAM6 value)

   !! To be namelist parameters in NEMO:
   REAL(wp), PARAMETER :: rz0_i_s_0  = 0.69e-3_wp  !           Eq. 43 [m]
   REAL(wp), PARAMETER :: rz0_i_f_0  = 4.54e-4_wp  ! bottom p.562 MIZ [m]

   LOGICAL,  PARAMETER :: l_add_form_drag = .TRUE.
   LOGICAL,  PARAMETER :: l_use_pond_info = .FALSE.
   LOGICAL,  PARAMETER :: l_dbg_print     = .FALSE.

   PUBLIC :: turb_ice_lg15

   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE turb_ice_lg15( kt, zt, zu, Ts_i, t_zt, qs_i, q_zt, U_zu, frice, &
      &                      Cd_i, Ch_i, Ce_i, t_zu_i, q_zu_i, Ubzu,            &
      &                      CdN_frm, &
      &                      CdN, ChN, CeN, xz0, xu_star, xL, xUN10 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_ice_lg15  ***
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
      !!    *  t_zt : potential air temperature at zt                         [K]
      !!    *  qs_i  : saturation specific humidity at temp. Ts_i over ice    [kg/kg]
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
      !!    *  Ubzu  : bulk wind speed at zu that was used                    [m/s]
      !!
      !! OPTIONAL OUTPUT:
      !! ----------------
      !!    * CdN_frm : sea-ice-related neutral FORM drag coefficient caused by floe egges etc...
      !!    * CdN     : neutral-stability drag coefficient
      !!    * ChN     : neutral-stability sensible heat coefficient
      !!    * CeN     : neutral-stability evaporation coefficient
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
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: t_zt  ! potential air temperature              [Kelvin]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: qs_i  ! sat. spec. hum. at ice/air interface    [kg/kg]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: q_zt  ! spec. air humidity at zt               [kg/kg]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: U_zu  ! relative wind module at zu                [m/s]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: frice ! sea-ice concentration        (fraction)
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Cd_i  ! drag coefficient over sea-ice
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Ch_i  ! transfert coefficient for heat over ice
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Ce_i  ! transfert coefficient for sublimation over ice
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: t_zu_i ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: q_zu_i ! spec. humidity adjusted at zu           [kg/kg]
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Ubzu     ! bulk wind speed at zu                     [m/s]
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj), OPTIONAL :: CdN_frm  ! form drag
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj), OPTIONAL :: CdN
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj), OPTIONAL :: ChN
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj), OPTIONAL :: CeN
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj), OPTIONAL :: xz0  ! Aerodynamic roughness length   [m]
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj), OPTIONAL :: xu_star  ! u*, friction velocity
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj), OPTIONAL :: xL  ! zeta (zu/L)
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj), OPTIONAL :: xUN10  ! Neutral wind at zu
      !!
      INTEGER :: j_itt
      LOGICAL :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !!
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: xtmp1, xtmp2      ! temporary stuff
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: dt_zu, dq_zu, zt_zu, zq_zu  ! third dimension
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zz0_s, zz0_f, RiB ! third dimensions (size=2):
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zCd, zCh, zCdN_s, zChN_s, zCdN_f, zChN_f
      
      LOGICAL ::  lreturn_cdfrm=.FALSE., lreturn_cdn=.FALSE., lreturn_chn=.FALSE., lreturn_cen=.FALSE., &
         &        lreturn_z0=.FALSE., lreturn_ustar=.FALSE., lreturn_L=.FALSE., lreturn_UN10=.FALSE.
      !!
      CHARACTER(len=40), PARAMETER :: crtnm = 'turb_ice_lg15@mod_blk_ice_lg15.f90'
      !!----------------------------------------------------------------------------------
      ALLOCATE ( xtmp1(jpi,jpj), xtmp2(jpi,jpj) )
      ALLOCATE ( dt_zu(jpi,jpj), dq_zu(jpi,jpj), zt_zu(jpi,jpj), zq_zu(jpi,jpj) )
      ALLOCATE ( zz0_s(jpi,jpj),  zz0_f(jpi,jpj),    RiB(jpi,jpj), &
         &      zCdN_s(jpi,jpj), zChN_s(jpi,jpj), zCdN_f(jpi,jpj), zChN_f(jpi,jpj) )
      ALLOCATE ( zCd(jpi,jpj), zCh(jpi,jpj) )
      
      lreturn_cdfrm = PRESENT(CdN_frm)
      lreturn_cdn   = PRESENT(CdN)
      lreturn_chn   = PRESENT(ChN)
      lreturn_cen   = PRESENT(CeN)
      lreturn_z0    = PRESENT(xz0)
      lreturn_ustar = PRESENT(xu_star)
      lreturn_L     = PRESENT(xL)
      lreturn_UN10  = PRESENT(xUN10)

      l_zt_equal_zu = ( ABS(zu - zt) < 0.01_wp )

      !! Scalar wind speed cannot be below 0.2 m/s
      Ubzu = MAX( U_zu, wspd_thrshld_ice )

      !! First guess of temperature and humidity at height zu:
      zt_zu(:,:) = MAX( t_zt(:,:) ,   100._wp )   ! who knows what's given on masked-continental regions...
      zq_zu(:,:) = MAX( q_zt(:,:) , 0.1e-6_wp )   !               "

      !! Air-Ice & Air-Sea differences (and we don't want them to be 0!)
      dt_zu(:,:) = zt_zu(:,:) - Ts_i
      dq_zu(:,:) = zq_zu(:,:) - qs_i
      dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6_wp), dt_zu )
      dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9_wp), dq_zu )

      !! Very crude first guess:
      Cd_i(:,:) = rCd_ice
      Ch_i(:,:) = rCd_ice
      Ce_i(:,:) = rCd_ice

      !! For skin drag :
      zz0_s(:,:)  = rz0_i_s_0        !LOLO/RFI! ! Room for improvement. We use the same z0_skin everywhere (= rz0_i_s_0)...
      zCdN_s(:,:) = Cd_from_z0( zu, zz0_s(:,:) )
      zChN_s(:,:) = vkarmn2 / ( LOG( zu / zz0_s(:,:) ) * LOG( zu / (ralpha_0*zz0_s(:,:)) ) )     ! (Eq.11,12)  [ "" ]

      !! For form drag in MIZ:
      zz0_f(:,:)  = 0._wp
      zCdN_f(:,:) = 0._wp
      zChN_f(:,:) = 0._wp
      IF ( l_add_form_drag ) THEN
         zz0_f(:,:)  = rz0_i_f_0        !LOLO/RFI! ! Room for improvement. We use the same z0_form everywhere !!!
         zCdN_f(:,:) = CdN_f_LG15_light( zu, frice(:,:), zz0_f(:,:) )
         zChN_f(:,:) = zCdN_f(:,:)/( 1._wp + LOG(1._wp/ralpha_0)/vkarmn*SQRT(zCdN_f(:,:)) ) ! (Eq.60,61)   [ "" ]
      END IF

      !! Some other first guess values, needed to compute wind at zt:
      zCd(:,:) = zCdN_s(:,:) + zCdN_f(:,:)
      zCh(:,:) = zChN_s(:,:) + zChN_f(:,:)
      RiB(:,:) = Ri_bulk( zt, Ts_i(:,:), t_zt(:,:), qs_i(:,:), q_zt(:,:), Ubzu(:,:) )  ! over ice (index=1)


      !! ITERATION BLOCK
      DO j_itt = 1, nb_itt

         IF(l_dbg_print) PRINT *, 'LOLO: LOOP #', INT(j_itt,1)
         IF(l_dbg_print) PRINT *, 'LOLO: theta_zu, Ts_i, Ubzu =', REAL(zt_zu(:,:),4), REAL(Ts_i(:,:),4), REAL(Ubzu(:,:),4)
         IF(l_dbg_print) PRINT *, 'LOLO:     q_zu =', REAL(zq_zu(:,:),4)
         IF(l_dbg_print) PRINT *, 'LOLO:  CdN_s, zCdN_f   =', REAL(zCdN_s(:,:),4), REAL(zCdN_f(:,:),4)


         !! Bulk Richardson Number
         !! ======================
         !! PROBLEM: when computed at z=zu, with adjusted theta and q, it is numerically unstable in some rare events (unstable)
         !!          => fix: compute RiB at zt, with ajusted wind at zt... => seems to be more stable
         IF( .NOT. l_zt_equal_zu ) THEN
            ! U_zt = U_zu + u_star/vkarmn*(LOG(zt/zu) + psi_m_coare(zu/L) - psi_m_coare(zt/L))
            xtmp1(:,:) = zCdN_s(:,:) + zCdN_f(:,:)    ! total neutral drag coeff!
            xtmp2(:,:) = zz0_s(:,:) + zz0_f(:,:)      ! total roughness length z0
            xtmp1 = LOG(zt/zu) + f_h_louis( zu, RiB(:,:), xtmp1(:,:), xtmp2(:,:) ) &
               &               - f_h_louis( zt, RiB(:,:), xtmp1(:,:), xtmp2(:,:) )
            xtmp2(:,:) = MAX( Ubzu(:,:) + (SQRT(zCd(:,:))*Ubzu)*xtmp1 , wspd_thrshld_ice ) ! wind at zt ( SQRT(zCd(:,:))*Ubzu == u* !)
            xtmp2(:,:) = MIN( xtmp2(:,:) , Ubzu(:,:) )
            IF(l_dbg_print) PRINT *, 'LOLO: ADJUSTED WIND AT ZT =', xtmp2
         ELSE
            xtmp2(:,:) = Ubzu(:,:)
         END IF
         RiB(:,:) = Ri_bulk( zt, Ts_i(:,:), t_zt(:,:), qs_i(:,:), q_zt(:,:), xtmp2(:,:) )  ! over ice (index=1)
         IF(l_dbg_print) PRINT *, 'LOLO: RiB_zt =', RiB(:,:)


         ! Momentum and Heat transfer coefficients WITHOUT FORM DRAG / (Eq.6) and (Eq.10):
         zCd(:,:) = zCdN_s(:,:) * f_m_louis( zu, RiB(:,:), zCdN_s(:,:), zz0_s(:,:) ) ! (Eq.6)
         zCh(:,:) = zChN_s(:,:) * f_h_louis( zu, RiB(:,:), zCdN_s(:,:), zz0_s(:,:) ) ! (Eq.10) / LOLO: why "zCdN_s" (xtmp1) and not "zChn" ???
         IF(l_dbg_print) PRINT *, 'LOLO: f_m_louis_s =', f_m_louis( zu, RiB(:,:), zCdN_s(:,:), zz0_s(:,:) )
         IF(l_dbg_print) PRINT *, 'LOLO: f_h_louis_s =', f_h_louis( zu, RiB(:,:), zCdN_s(:,:), zz0_s(:,:) )
         IF(l_dbg_print) PRINT *, 'LOLO: Cd / skin only / ice   =', REAL(zCd(:,:),4)



         IF ( l_add_form_drag ) THEN
            !! Form-drag-related NEUTRAL momentum and Heat transfer coefficients:
            !!   MIZ:
            zCd(:,:) = zCd(:,:) + zCdN_f(:,:) * f_m_louis( zu, RiB(:,:), zCdN_f(:,:), zz0_f(:,:) ) ! (Eq.6)
            zCh(:,:) = zCh(:,:) + zChN_f(:,:) * f_h_louis( zu, RiB(:,:), zCdN_f(:,:), zz0_f(:,:) ) ! (Eq.10) / LOLO: why "zCdN_f" and not "zChn" ???
            IF(l_dbg_print) PRINT *, 'LOLO: f_m_louis_f =', f_m_louis( zu, RiB(:,:), zCdN_f(:,:), zz0_f(:,:) )
            IF(l_dbg_print) PRINT *, 'LOLO: f_h_louis_f =', f_h_louis( zu, RiB(:,:), zCdN_f(:,:), zz0_f(:,:) )

            IF(l_dbg_print) PRINT *, 'LOLO: Cd / form only / ice   =', REAL(zCdN_f(:,:) * f_m_louis( zu, RiB(:,:), zCdN_f(:,:), zz0_f(:,:) ),4)

         END IF

         IF(l_dbg_print) PRINT *, 'LOLO: Cd, Ch / TOTAL / ice   =', REAL(zCd(:,:),4), REAL(zCh(:,:),4)


         !! Adjusting temperature and humidity from zt to zu:
         IF( .NOT. l_zt_equal_zu ) THEN

            !! Over ice:
            xtmp1(:,:) = zCdN_s(:,:) + zCdN_f(:,:)    ! total neutral drag coeff!
            xtmp2(:,:) = zz0_s(:,:) + zz0_f(:,:)      ! total roughness length z0
            xtmp1 = LOG(zt/zu) + f_h_louis( zu, RiB(:,:), xtmp1(:,:), xtmp2(:,:) ) &
               &               - f_h_louis( zt, RiB(:,:), xtmp1(:,:), xtmp2(:,:) )
            xtmp2 = 1._wp/SQRT(zCd(:,:))

            zt_zu(:,:) = t_zt - (zCh(:,:) * dt_zu(:,:) * xtmp2) / vkarmn * xtmp1   ! t_star = Ch * dt_zu / SQRT(Cd)
            zq_zu(:,:) = q_zt - (zCh(:,:) * dq_zu(:,:) * xtmp2) / vkarmn * xtmp1   ! q_star = Ce * dq_zu / SQRT(Cd)
            zq_zu(:,:) = MAX(0._wp, zq_zu(:,:))

            dt_zu(:,:) = zt_zu(:,:) - Ts_i
            dq_zu(:,:) = zq_zu(:,:) - qs_i


            dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6_wp), dt_zu )
            dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9_wp), dq_zu )
         END IF

         IF(l_dbg_print) PRINT *, ''!LOLO

      END DO !DO j_itt = 1, nb_itt
      IF(l_dbg_print) PRINT *, ''!LOLO


      IF( lreturn_cdfrm ) CdN_frm = zCdN_f(:,:)      
      IF( lreturn_cdn )   CdN = zCdN_s(:,:)+zCdN_f(:,:)
      IF( lreturn_chn )   ChN = zChN_s(:,:)+zChN_f(:,:)
      IF( lreturn_cen )   CeN = zChN_s(:,:)+zChN_f(:,:)

      !! Result is ice + ocean:
      !t_zu(:,:) = mix_val_msh(zt_zu, frice)
      !q_zu(:,:) = mix_val_msh(zq_zu, frice)
      !Cd(:,:) = mix_val_msh(zCd, frice)
      !Ch(:,:) = mix_val_msh(zCh, frice)

      !! Result is over ice only:
      t_zu_i(:,:) = zt_zu(:,:)
      q_zu_i(:,:) = zq_zu(:,:)
      Cd_i(:,:)   =   zCd(:,:)
      Ch_i(:,:)   =   zCh(:,:)
      Ce_i(:,:)   =  Ch_i(:,:)


      IF( lreturn_z0 ) xz0   = z0_from_Cd( zu, zCdN_s(:,:)+zCdN_f(:,:) )

      IF( lreturn_ustar ) xu_star = SQRT(Cd_i) * Ubzu
      IF( lreturn_L ) THEN
         xtmp1 = SQRT(Cd_i)
         xL    = 1./One_on_L( t_zu_i, q_zu_i, xtmp1*Ubzu, Ch_i*dt_zu(:,:)/xtmp1, Ce_i*dq_zu(:,:)/xtmp1 )
      END IF

      IF( lreturn_UN10 ) THEN
         xtmp1 = zCdN_s(:,:) + zCdN_f(:,:)  ! => CdN
         xUN10 = SQRT(Cd_i) * Ubzu/vkarmn * LOG( 10._wp / z0_from_Cd(zu, xtmp1) )
         !xtmp2 = f_m_louis( zu, RiB(:,:), xtmp1, z0_from_Cd(zu, xtmp1) ) ! => f_m
         !xUN10 = UN10_from_CD( zu, Ubzu, Cd_i, ppsi=xtmp2 )
      END IF



      DEALLOCATE ( xtmp1, xtmp2 )
      DEALLOCATE ( dt_zu, dq_zu, zt_zu, zq_zu )
      DEALLOCATE ( zz0_s, zz0_f, RiB, zCdN_s, zChN_s, zCdN_f, zChN_f )
      DEALLOCATE ( zCd, zCh )

   END SUBROUTINE turb_ice_lg15

   !!======================================================================

END MODULE mod_blk_ice_lg15
