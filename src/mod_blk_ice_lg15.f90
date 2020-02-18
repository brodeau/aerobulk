! AeroBulk / 2020 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_blk_ice_lg15
   !!====================================================================================
   !!       Computes turbulent components of surface fluxes over sea-ice
   !!       Following Lupkes & Gryanik, 2015
   !!
   !!              => case when 100 % sea-ice
   !!
   !!       Routine turb_ice_lg15 maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, January 2020
   !!
   !!====================================================================================
   USE mod_const       !: physical and other constants
   USE mod_phymbl      !: misc. physical functions 

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: turb_ice_lg15 !, Cx_LG15, Cx_skin_LG15

   ! ECHAM6 constants
   REAL(wp), PARAMETER ::   rz0_s_0  = 0.69e-3_wp  ! Eq. 43 [m]
   REAL(wp), PARAMETER ::   rz0_f_0  = 4.54e-4_wp  ! bottom p.562 MIZ [m]
   REAL(wp), PARAMETER ::   rz0_0       = 1.00e-3_wp  ! Eq. 15 [m]
   REAL(wp), PARAMETER ::   rce10_i_0   = 3.46e-3_wp  ! Eq. 48 MIZ
   REAL(wp), PARAMETER ::   rbeta_0     = 1.4_wp      ! Eq. 47 MIZ
   REAL(wp), PARAMETER ::   ralpha_0     = 0.2_wp      ! Eq. 12
   REAL(wp), PARAMETER ::   rbetah_0     = 1.e-3_wp    ! Eq. 26
   REAL(wp), PARAMETER ::   rgamma_0     = 1.25_wp     ! Eq. 26
   REAL(wp), PARAMETER ::   r1_gamma_0   = 1._wp / rgamma_0
   
   LOGICAL, PARAMETER :: l_use_form_drag = .TRUE.
   LOGICAL, PARAMETER :: l_use_pond_info = .FALSE.

   
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE turb_ice_lg15( kt, zt, zu, Ti_s, t_zt, qi_s, q_zt, U_zu, &
      &                     Cd, Ch, Ce, t_zu, q_zu, U_blk,             &
      &                     xz0, xu_star, xL, xUN10 )
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
      !!    *  Ti_s  : surface temperature of sea-ice                         [K]
      !!    *  t_zt : potential air temperature at zt                         [K]
      !!    *  qi_s  : SSQ aka saturation specific humidity at temp. Ti_s     [kg/kg]
      !!    *  q_zt : specific humidity of air at zt                          [kg/kg]
      !!    *  U_zu : scalar wind speed at zu                                 [m/s]
      !!
      !! OUTPUT :
      !! --------
      !!    *  Cd     : drag coefficient
      !!    *  Ch     : sensible heat coefficient
      !!    *  Ce     : evaporation coefficient
      !!    *  t_zu   : pot. air temperature adjusted at wind height zu       [K]
      !!    *  q_zu   : specific humidity of air        //                    [kg/kg]
      !!    *  U_blk  : bulk wind speed at zu                                 [m/s]
      !!
      !! OPTIONAL OUTPUT:
      !! ----------------
      !!    * xz0         : return the aerodynamic roughness length (integration constant for wind stress) [m]
      !!    * xu_star     : return u* the friction velocity                    [m/s]
      !!    * xL          : return the Obukhov length                          [m]
      !!    * xUN10       : neutral wind speed at 10m                          [m/s]
      !!
      !! ** Author: L. Brodeau, January 2020 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      INTEGER,  INTENT(in )                     ::   kt       ! current time step
      REAL(wp), INTENT(in )                     ::   zt       ! height for t_zt and q_zt                    [m]
      REAL(wp), INTENT(in )                     ::   zu       ! height for U_zu                             [m]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) ::   Ti_s     ! ice surface temperature                [Kelvin]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) ::   t_zt     ! potential air temperature              [Kelvin]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) ::   qi_s     ! specific humidity at ice/air interface  [kg/kg]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) ::   q_zt     ! specific air humidity at zt             [kg/kg]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) ::   U_zu     ! relative wind module at zu                [m/s]
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) ::   Cd       ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) ::   Ch       ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) ::   Ce       ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) ::   t_zu     ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) ::   q_zu     ! spec. humidity adjusted at zu           [kg/kg]
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) ::   U_blk    ! bulk wind speed at zu                     [m/s]
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(jpi,jpj) ::   xz0  ! Aerodynamic roughness length   [m]
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(jpi,jpj) ::   xu_star  ! u*, friction velocity
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(jpi,jpj) ::   xL  ! zeta (zu/L)
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(jpi,jpj) ::   xUN10  ! Neutral wind at zu
      !
      INTEGER :: j_itt
      LOGICAL :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !
      REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: xtmp1, xtmp2      ! temporary stuff
      REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: dt_zu, dq_zu
      REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: zfrac_ice      
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zz0_s, zz0_f, RiB, zCdN_s, zChN_s, zCdN_f, zChN_f ! third dimensions (size=2): 1 => ice, 2 => water

      LOGICAL :: lreturn_z0=.FALSE., lreturn_ustar=.FALSE., lreturn_L=.FALSE., lreturn_UN10=.FALSE.
      CHARACTER(len=40), PARAMETER :: crtnm = 'turb_ice_lg15@mod_blk_ice_lg15.f90'
      !!----------------------------------------------------------------------------------
      ALLOCATE ( xtmp1(jpi,jpj), xtmp2(jpi,jpj) )
      ALLOCATE ( dt_zu(jpi,jpj), dq_zu(jpi,jpj), zfrac_ice(jpi,jpj) )
      ALLOCATE ( zz0_s(jpi,jpj,2), zz0_f(jpi,jpj,2), RiB(jpi,jpj,2), zCdN_s(jpi,jpj,2), zChN_s(jpi,jpj,2), zCdN_f(jpi,jpj,2), zChN_f(jpi,jpj,2) )

      
      IF( PRESENT(xz0) )     lreturn_z0    = .TRUE.
      IF( PRESENT(xu_star) ) lreturn_ustar = .TRUE.
      IF( PRESENT(xL) )      lreturn_L     = .TRUE.
      IF( PRESENT(xUN10) )   lreturn_UN10  = .TRUE.

      l_zt_equal_zu = .FALSE.
      IF( ABS(zu - zt) < 0.01_wp )   l_zt_equal_zu = .TRUE.    ! testing "zu == zt" is risky with double precision

      !! Scalar wind speed cannot be below 0.2 m/s
      U_blk = MAX( U_zu, 0.2_wp )
           
      !! First guess of temperature and humidity at height zu:
      t_zu = MAX( t_zt ,   100._wp )   ! who knows what's given on masked-continental regions...
      q_zu = MAX( q_zt , 0.1e-6_wp )   !               "
      
      !! Air-Ice differences (and we don't want it to be 0!)
      dt_zu = t_zu - Ti_s ;   dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6_wp), dt_zu )
      dq_zu = q_zu - qi_s ;   dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9_wp), dq_zu )

      !! Very crude first guess:
      Cd(:,:) = rCd_ice
      Ch(:,:) = rCd_ice
      Ce(:,:) = rCd_ice

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
         zfrac_ice(:,:) = 0.9 !!! LOLO: fix, use model value !
         xtmp1(:,:) = 1._wp / zz0_f(:,:,1)
         xtmp2(:,:) = rce10_i_0 * ( LOG( 10._wp * xtmp1(:,:) ) / LOG( zu * xtmp1(:,:) ) )**2       ! part of (Eq.46)
         zCdN_f(:,:,1) = xtmp2(:,:) * zfrac_ice(:,:) * (1._wp - zfrac_ice(:,:))**rbeta_0           ! (Eq.46)  [ index 1 is for ice, 2 for water ]
         zChN_f(:,:,1) = zCdN_f(:,:,1) / ( 1._wp + LOG(1._wp/ralpha_0)/vkarmn * SQRT(zCdN_f(:,:,1)) ) ! (Eq.60,61)   [ "" ]
      END IF
      
      !! ITERATION BLOCK
      DO j_itt = 1, nb_itt
                  
         ! Bulk Richardson Number:
         RiB(:,:,1) = Ri_bulk( zu, Ti_s(:,:), t_zu(:,:), qi_s(:,:), q_zu(:,:), U_blk(:,:) )
         
         ! Momentum and Heat transfer coefficients WITHOUT FORM DRAG / (Eq.6) and (Eq.10):
         Cd(:,:) = zCdN_s(:,:,1) * f_m_louis( zu, RiB(:,:,1), zCdN_s(:,:,1), zz0_s(:,:,1) ) ! (Eq.6)
         Ch(:,:) = zChN_s(:,:,1) * f_h_louis( zu, RiB(:,:,1), zCdN_s(:,:,1), zz0_s(:,:,1) ) ! (Eq.10) / LOLO: why "zCdN_s" (xtmp1) and not "zChn" ???
         PRINT *, 'LOLO: Cd / skin only / ice =', REAL(Cd,4)


         IF ( l_use_form_drag ) THEN
            !! Form-drag-related NEUTRAL momentum and Heat transfer coefficients:
            !!   MIZ:
            Cd(:,:) = Cd(:,:) + zCdN_f(:,:,1) * f_m_louis( zu, RiB(:,:,1), zCdN_f(:,:,1), zz0_f(:,:,1) ) ! (Eq.6)
            Ch(:,:) = Ch(:,:) + zChN_f(:,:,1) * f_h_louis( zu, RiB(:,:,1), zCdN_f(:,:,1), zz0_f(:,:,1) ) ! (Eq.10) / LOLO: why "zCdN_f" (xtmp1) and not "zChn" ???
         END IF
         
         Ce(:,:) = Ch(:,:)
         PRINT *, 'LOLO: Cd / total / ice =', REAL(Cd,4)
         
         !! Adjusting temperature and humidity from zt to zu:
         IF( .NOT. l_zt_equal_zu ) THEN
            xtmp1(:,:) = zCdN_s(:,:,1) + zCdN_f(:,:,1)    ! total neutral drag coeff!
            xtmp2(:,:) = zz0_s(:,:,1) + zz0_f(:,:,1)      ! total roughness length z0
            xtmp1 = LOG(zt/zu) + f_h_louis( zu, RiB(:,:,1), xtmp1(:,:), xtmp2(:,:) ) &
               &               - f_h_louis( zt, RiB(:,:,1), xtmp1(:,:), xtmp2(:,:) )
            xtmp2 = 1._wp/SQRT(Cd)
            t_zu = t_zt - (Ch * dt_zu * xtmp2) / vkarmn * xtmp1   ! t_star = Ch * dt_zu / SQRT(Cd)
            q_zu = q_zt - (Ce * dq_zu * xtmp2) / vkarmn * xtmp1   ! q_star = Ce * dq_zu / SQRT(Cd)
            q_zu = MAX(0._wp, q_zu)
            PRINT *, 'LOLO: height adjustment (mod_blk_ice_lg15)! Cd=', REAL(Cd,4)
            dt_zu = t_zu - Ti_s ;   dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6_wp), dt_zu )
            dq_zu = q_zu - qi_s ;   dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9_wp), dq_zu )
         END IF

         PRINT *, ''!LOLO         
      END DO !DO j_itt = 1, nb_itt


      PRINT *, 'LOLO: CDN10_skin_ice =', REAL( zCdN_s(:,:,1), 4)
      
      IF( lreturn_z0 ) xz0   = z0_from_Cd( zu, zCdN_s(:,:,1)+zCdN_f(:,:,1) )
      
      IF( lreturn_ustar ) xu_star = SQRT(Cd) * U_blk
      IF( lreturn_L ) THEN
         xtmp1 = SQRT(Cd)
         xL    = 1./One_on_L(t_zu, q_zu, xtmp1*U_blk, Ch*dt_zu/xtmp1, Ce*dq_zu/xtmp1)
      END IF
      
      IF( lreturn_UN10 ) xUN10   = (SQRT(Cd) * U_blk) / vkarmn*LOG(10./(z0_from_Cd( zu, zCdN_s(:,:,1)+zCdN_f(:,:,1) )) )
      
      DEALLOCATE ( xtmp1, xtmp2 )
      DEALLOCATE ( dt_zu, dq_zu, zfrac_ice )      
      DEALLOCATE ( zz0_s, zz0_f, RiB, zCdN_s, zChN_s, zCdN_f, zChN_f )

   END SUBROUTINE turb_ice_lg15





   SUBROUTINE Cx_LG15( zu, t_zu, q_zu, Ui_zu, Ts_i, qs_i, pcd, pch )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  Cx_LG15  ***
      !!
      !!                         CASE 100 % sea-ice covered !!!
      !!
      !! ** Purpose :    Alternative turbulent transfer coefficients formulation
      !!                 between sea-ice and atmosphere with distinct momentum
      !!                 and heat coefficients depending on sea-ice concentration
      !!                 and atmospheric stability (no meltponds effect for now).
      !!
      !! ** Method :     The parameterization is adapted from Lupkes & Gryanik (2015)
      !!                 and ECHAM6 atmospheric model. Compared to Lupkes2012 scheme,
      !!                 it considers specific skin and form drags (Andreas et al. 2010)
      !!                 to compute neutral transfer coefficients for both heat and
      !!                 momemtum fluxes. Atmospheric stability effect on transfer
      !!                 coefficient is also taken into account following Louis (1979).
      !!
      !! ** References : Lupkes & Gryanik JGR 2015 (theory)
      !!                 Lupkes & Gryanik ECHAM6 documentation 2015 (implementation)
      !!
      !! ** History :
      !!              - G. Samson (2018,2019) original code
      !!              - L. Brodeau (2020) AeroBulk
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   )                   :: zu     ! reference height (height for Uo_zu)   [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: t_zu   ! potential air temperature              [Kelvin]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: q_zu   ! specific air humidity at zt             [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: Ui_zu  ! relative wind module at zu over ice    [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: Ts_i   ! sea-ice surface temperature               [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: qs_i   ! humidity at saturation over ice at T=Ts_i [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pcd    ! momentum transfer coefficient
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pch    ! heat transfer coefficient
      !
      REAL(wp) ::   zfi, zfo
      REAL(wp) ::   zrib_i, ztmp
      REAL(wp) ::   zCdn_skin, zCdn_form, zCdn
      REAL(wp) ::   zChn_skin, zChn_form
      REAL(wp) ::   z0i, zfmi, zfhi
      REAL(wp) ::   zCdn_form_tmp, zwndspd_i
      INTEGER  ::   ji, jj         ! dummy loop indices
      !!----------------------------------------------------------------------

      ! Momentum Neutral Transfer Coefficients (should be a constant)
      zCdn_form_tmp = rce10_i_0 * ( LOG( 10._wp / rz0_f_0 + 1._wp ) / LOG( zu / rz0_f_0 + 1._wp ) )**2   ! Eq. 46
      zCdn_skin = ( vkarmn                                      / LOG( zu / rz0_s_0 + 1._wp ) )**2   ! Eq. 7
      zCdn      = zCdn_skin   ! Eq. 7
      !zCdn     = 1.89e-3         ! old ECHAM5 value (cf Eq. 32)

      ! Heat Neutral Transfer Coefficients
      zChn_skin = vkarmn**2 / ( LOG( zu / rz0_0 + 1._wp ) * LOG( zu * ralpha_0 / rz0_s_0 + 1._wp ) )   ! Eq. 50 + Eq. 52

      DO jj = 1, jpj
         DO ji = 1, jpi

            zfi  = 1._wp  ! fraction of sea-ice
            zwndspd_i = MAX( 0.5, Ui_zu(ji,jj) )

            zfo  = 0._wp  ! fraction of open ocean

            ! Bulk Richardson Number:
            zrib_i = Ri_bulk( zu, Ts_i(ji,jj), t_zu(ji,jj), qs_i(ji,jj), q_zu(ji,jj), zwndspd_i )

            ! Momentum and Heat Neutral Transfer Coefficients
            zCdn_form = zCdn_form_tmp * zfi * zfo**rbeta_0                          ! Eq. 40 !LOLO: WHAT????? zfi * zfo is always 0 !!!
            zChn_form = zCdn_form / ( 1._wp + ( LOG( 1._wp/ralpha_0 ) / vkarmn ) * SQRT( zCdn_form ) )   ! Eq. 53

            ! Momentum and Heat Stability functions (possibility to use psi_m_ecmwf instead ?)
            z0i = rz0_s_0                                        ! over ice

            zfmi = f_m_louis( zu, zrib_i, zCdn, z0i )
            zfhi = f_h_louis( zu, zrib_i, zCdn, z0i )  !LOLO: why "zCdn" and not "zChn" ???

            ! Momentum and Heat transfer coefficients (Eq. 38) and (Eq. 49):
            ztmp       = 1._wp / MAX( 1.e-06, zfi )
            pcd(ji,jj) = zCdn_skin * zfmi  +  zCdn_form * ( zfmi*zfi ) * ztmp
            pch(ji,jj) = zChn_skin * zfhi  +  zChn_form * ( zfhi*zfi ) * ztmp

         END DO
      END DO
      !
   END SUBROUTINE Cx_LG15


   
   SUBROUTINE Cx_skin_LG15( zu, t_zu, q_zu, Ui_zu, Ts_i, qs_i, pcd, pch )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  Cx_skin_LG15  ***
      !!
      !!                         CASE 100 % sea-ice covered !!!
      !!
      !!                        Cd and Ch, WITHOUT FORM DRAG
      !!
      !! ** Purpose :    Alternative turbulent transfer coefficients formulation
      !!                 between sea-ice and atmosphere with distinct momentum
      !!                 and heat coefficients depending on sea-ice concentration
      !!                 and atmospheric stability (no meltponds effect for now).
      !!
      !! ** Method :     The parameterization is adapted from Lupkes & Gryanik (2015)
      !!                 and ECHAM6 atmospheric model. Compared to Lupkes2012 scheme,
      !!                 it considers specific skin and form drags (Andreas et al. 2010)
      !!                 to compute neutral transfer coefficients for both heat and
      !!                 momemtum fluxes. Atmospheric stability effect on transfer
      !!                 coefficient is also taken into account following Louis (1979).
      !!
      !! ** References : Lupkes & Gryanik JGR 2015 (theory)
      !!                 Lupkes & Gryanik ECHAM6 documentation 2015 (implementation)
      !!
      !! ** History :
      !!              - G. Samson (2018,2019) original code
      !!              - L. Brodeau (2020) AeroBulk
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   )                   :: zu     ! reference height (height for Uo_zu)   [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: t_zu   ! potential air temperature              [Kelvin]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: q_zu   ! specific air humidity at zt             [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: Ui_zu  ! relative wind module at zu over ice    [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: Ts_i   ! sea-ice surface temperature               [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: qs_i   ! humidity at saturation over ice at T=Ts_i [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pcd    ! momentum transfer coefficient
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pch    ! heat transfer coefficient
      !
      REAL(wp) ::   zrib_i
      REAL(wp) ::   zCdn_skin, zChn_skin
      REAL(wp) ::   ztmp
      REAL(wp) ::   zwndspd_i
      INTEGER  ::   ji, jj         ! dummy loop indices
      !!----------------------------------------------------------------------
      
      ! Momentum Neutral Transfer Coefficients
      !zCdn_skin = vkarmn2 / ( LOG( zu / rz0_s_0 + 1._wp ) )**2   ! Eq. 7
      !zChn_skin = vkarmn2 / ( LOG( zu / z0      + 1._wp ) * LOG( zu * ralpha_0 / rz0_s_0 + 1._wp ) )   ! Eq. 50 + Eq. 52
      ztmp = LOG( zu / rz0_s_0 )
      zCdn_skin = vkarmn2 / ( ztmp * ztmp )   ! Eq. 7
      zChn_skin = vkarmn2 / ( ztmp * LOG( zu / (ralpha_0*rz0_s_0) ) )   ! Eq. 11, 12
      
      DO jj = 1, jpj
         DO ji = 1, jpi
            
            zwndspd_i = MAX( 0.5, Ui_zu(ji,jj) )
            
            ! Bulk Richardson Number:
            zrib_i = Ri_bulk( zu, Ts_i(ji,jj), t_zu(ji,jj), qs_i(ji,jj), q_zu(ji,jj), zwndspd_i )
            
            ! Momentum and Heat transfer coefficients (Eq. 38) and (Eq. 49):
            pcd(ji,jj) = zCdn_skin * f_m_louis( zu, zrib_i, zCdn_skin, rz0_s_0 )
            pch(ji,jj) = zChn_skin * f_h_louis( zu, zrib_i, zCdn_skin, rz0_s_0 ) !LOLO: why "zCdn_skin" and not "zChn" ???

         END DO
      END DO
      !
   END SUBROUTINE Cx_skin_LG15







   
   !!======================================================================

END MODULE mod_blk_ice_lg15
