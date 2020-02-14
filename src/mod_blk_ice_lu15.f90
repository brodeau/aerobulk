! AeroBulk / 2020 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_blk_ice_lu15
   !!====================================================================================
   !!       Computes turbulent components of surface fluxes over sea-ice




   !!       Routine turb_ice_lu15 maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, January 2020
   !!
   !!====================================================================================
   USE mod_const       !: physical and other constants
   USE mod_phymbl      !: thermodynamics

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: turb_ice_lu15, Cx_Lupkes2015

   ! ECHAM6 constants
   REAL(wp), PARAMETER ::   z0_skin_ice  = 0.69e-3_wp  ! Eq. 43 [m]
   REAL(wp), PARAMETER ::   z0_form_ice  = 0.57e-3_wp  ! Eq. 42 [m]
   REAL(wp), PARAMETER ::   z0_ice       = 1.00e-3_wp  ! Eq. 15 [m]
   REAL(wp), PARAMETER ::   zce10        = 2.80e-3_wp  ! Eq. 41
   REAL(wp), PARAMETER ::   zbeta        = 1.1_wp      ! Eq. 41
   REAL(wp), PARAMETER ::   z1_alpha     = 1._wp / 0.2_wp  ! Eq. 51
   REAL(wp), PARAMETER ::   z1_alphaf    = z1_alpha    ! Eq. 56
   REAL(wp), PARAMETER ::   zbetah       = 1.e-3_wp    ! Eq. 26
   REAL(wp), PARAMETER ::   zgamma       = 1.25_wp     ! Eq. 26
   REAL(wp), PARAMETER ::   z1_gamma     = 1._wp / zgamma
   REAL(wp), PARAMETER ::   r1_3         = 1._wp / 3._wp

   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE turb_ice_lu15( kt, zt, zu, Ti_s, t_zt, qi_s, q_zt, U_zu, &
      &                     Cd, Ch, Ce, t_zu, q_zu, U_blk,             &
      &                     xz0, xu_star, xL, xUN10 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_ice_lu15  ***
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
      !
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(jpi,jpj) ::   xz0  ! Aerodynamic roughness length   [m]
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(jpi,jpj) ::   xu_star  ! u*, friction velocity
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(jpi,jpj) ::   xL  ! zeta (zu/L)
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(jpi,jpj) ::   xUN10  ! Neutral wind at zu
      !
      INTEGER :: j_itt
      LOGICAL :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !
      REAL(wp), DIMENSION(:,:), ALLOCATABLE  ::  &
         &  u_star, t_star, q_star, &
         &  dt_zu, dq_zu

      LOGICAL :: lreturn_z0=.FALSE., lreturn_ustar=.FALSE., lreturn_L=.FALSE., lreturn_UN10=.FALSE.
      CHARACTER(len=40), PARAMETER :: crtnm = 'turb_ice_lu15@mod_blk_ice_lu15.f90'
      !!----------------------------------------------------------------------------------
      ALLOCATE ( u_star(jpi,jpj), t_star(jpi,jpj), q_star(jpi,jpj),  &
         &       dt_zu(jpi,jpj),  dq_zu(jpi,jpj) )

      IF( PRESENT(xz0) )     lreturn_z0    = .TRUE.
      IF( PRESENT(xu_star) ) lreturn_ustar = .TRUE.
      IF( PRESENT(xL) )      lreturn_L     = .TRUE.
      IF( PRESENT(xUN10) )   lreturn_UN10  = .TRUE.

      l_zt_equal_zu = .FALSE.
      IF( ABS(zu - zt) < 0.01_wp )   l_zt_equal_zu = .TRUE.    ! testing "zu == zt" is risky with double precision

      !! Scalar wind speed cannot be below 0.2 m/s
      U_blk = MAX( U_zu, 0.2_wp )

      !! Initializing values at z_u with z_t values:
      t_zu = t_zt
      q_zu = q_zt
      
      !! Pot. temp. difference (and we don't want it to be 0!)
      dt_zu = t_zu - Ti_s ;   dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6_wp), dt_zu )
      dq_zu = q_zu - qi_s ;   dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9_wp), dq_zu )

      CALL Cx_Lupkes2015( zu, t_zu, q_zu, U_blk, Ti_s, qi_s, Cd, Ch )
      Ce = Ch

      u_star = SQRT(Cd) * U_blk
      t_star = Cd * U_blk * dt_zu / u_star
      q_star = Cd * U_blk * dq_zu / u_star
      
      !! ITERATION BLOCK
      DO j_itt = 1, nb_itt
         !
         dt_zu = t_zu - Ti_s ;   dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6_wp), dt_zu )
         dq_zu = q_zu - qi_s ;   dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9_wp), dq_zu )
         
         ! Updating turbulent scales :   (L&Y 2004 Eq. (7))
         u_star = SQRT(Cd) * U_blk
         t_star = Cd * U_blk * dt_zu / u_star
         q_star = Cd * U_blk * dq_zu / u_star

         !! Shifting temperature and humidity at zu (L&Y 2004 Eq. (9b-9c))
         IF( .NOT. l_zt_equal_zu ) THEN
            !!t_zu = t_zt - t_star/vkarmn * LOG(zt/zu) + f_h_louis( zu, Rib, Chn, z0 ) - f_h_louis( zt, Rib, Chn, z0 )
            !!q_zu = q_zt - q_star/vkarmn * LOG(zt/zu) + f_h_louis( zu, Rib, Chn, z0 ) - f_h_louis( zt, Rib, Chn, z0 )
            !!q_zu = MAX(0._wp, q_zu)
            PRINT *, 'LOLO: fix me height adjustment into mod_blk_ice_lu15.f90 !!!'
            !!
            !! PROBLEM HERE IS THAT WE DO NOT USE STABILITY FUNCTIONS PSI !!!
            !! => find out the way to adjust at zu based on Louis functions !!!
            !ztmp0 = zt*ztmp0 ! zeta_t !
            !ztmp0 = SIGN( MIN(ABS(ztmp0),10._wp), ztmp0 )  ! Temporaty array ztmp0 == zeta_t !!!
            !ztmp0 = LOG(zt/zu) + psi_h_ice(zeta_u) - psi_h_ice(ztmp0)                   ! ztmp0 just used as temp array again!
            !t_zu = t_zt - ztmp1/vkarmn*ztmp0    ! ztmp1 is still theta*  L&Y 2004 Eq. (9b)
            !q_zu = q_zt - ztmp2/vkarmn*ztmp0    ! ztmp2 is still q*      L&Y 2004 Eq. (9c)
            !q_zu = MAX(0._wp, q_zu)
         END IF

      END DO

      IF( lreturn_z0 )    xz0     = z0_from_Cd( zu, Cd )
      IF( lreturn_ustar ) xu_star = u_star
      IF( lreturn_L )     xL      = 1./One_on_L(t_zu, q_zu, u_star, t_star, q_star)
      !IF( lreturn_UN10 )  xUN10   = u_star/vkarmn*LOG(10./z0) !LOLO: fix me needs z0  !!!
      
      DEALLOCATE ( u_star, t_star, q_star, dt_zu, dq_zu )

   END SUBROUTINE turb_ice_lu15





   SUBROUTINE Cx_Lupkes2015( zu, t_zu, q_zu, Ui_zu, Ts_i, qs_i, pcd, pch )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  Cx_Lupkes2015  ***
      !!
      !!                         CASE 100 % sea-ice covered !!!
      !!
      !! ** Purpose :    Alternative turbulent transfer coefficients formulation
      !!                 between sea-ice and atmosphere with distinct momentum
      !!                 and heat coefficients depending on sea-ice concentration
      !!                 and atmospheric stability (no meltponds effect for now).
      !!
      !! ** Method :     The parameterization is adapted from Lupkes et al. (2015)
      !!                 and ECHAM6 atmospheric model. Compared to Lupkes2012 scheme,
      !!                 it considers specific skin and form drags (Andreas et al. 2010)
      !!                 to compute neutral transfer coefficients for both heat and
      !!                 momemtum fluxes. Atmospheric stability effect on transfer
      !!                 coefficient is also taken into account following Louis (1979).
      !!
      !! ** References : Lupkes et al. JGR 2015 (theory)
      !!                 Lupkes et al. ECHAM6 documentation 2015 (implementation)
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
      REAL(wp) ::   zCdn_skin_ice, zCdn_form_ice, zCdn_ice
      REAL(wp) ::   zChn_skin_ice, zChn_form_ice
      REAL(wp) ::   z0i, zfmi, zfhi
      REAL(wp) ::   zCdn_form_tmp, zwndspd_i
      INTEGER  ::   ji, jj         ! dummy loop indices
      !!----------------------------------------------------------------------

      ! Momentum Neutral Transfer Coefficients (should be a constant)
      zCdn_form_tmp = zce10 * ( LOG( 10._wp / z0_form_ice + 1._wp ) / LOG( zu / z0_form_ice + 1._wp ) )**2   ! Eq. 46
      zCdn_skin_ice = ( vkarmn                                      / LOG( zu / z0_skin_ice + 1._wp ) )**2   ! Eq. 7
      zCdn_ice      = zCdn_skin_ice   ! Eq. 7
      !zCdn_ice     = 1.89e-3         ! old ECHAM5 value (cf Eq. 32)

      ! Heat Neutral Transfer Coefficients
      zChn_skin_ice = vkarmn**2 / ( LOG( zu / z0_ice + 1._wp ) * LOG( zu * z1_alpha / z0_skin_ice + 1._wp ) )   ! Eq. 50 + Eq. 52

      DO jj = 1, jpj
         DO ji = 1, jpi

            zfi  = 1._wp  ! fraction of sea-ice
            zwndspd_i = MAX( 0.5, Ui_zu(ji,jj) )

            zfo  = 0._wp  ! fraction of open ocean

            ! Bulk Richardson Number:
            zrib_i = Ri_bulk( zu, Ts_i(ji,jj), t_zu(ji,jj), qs_i(ji,jj), q_zu(ji,jj), zwndspd_i )

            ! Momentum and Heat Neutral Transfer Coefficients
            zCdn_form_ice = zCdn_form_tmp * zfi * zfo**zbeta                          ! Eq. 40 !LOLO: WHAT????? zfi * zfo is always 0 !!!
            zChn_form_ice = zCdn_form_ice / ( 1._wp + ( LOG( z1_alphaf ) / vkarmn ) * SQRT( zCdn_form_ice ) )   ! Eq. 53

            ! Momentum and Heat Stability functions (possibility to use psi_m_ecmwf instead ?)
            z0i = z0_skin_ice                                        ! over ice

            zfmi = f_m_louis( zu, zrib_i, zCdn_ice, z0i )
            zfhi = f_h_louis( zu, zrib_i, zCdn_ice, z0i )  !LOLO: why "zCdn_ice" and not "zChn_ice" ???

            ! Momentum and Heat transfer coefficients (Eq. 38) and (Eq. 49):
            ztmp       = 1._wp / MAX( 1.e-06, zfi )
            pcd(ji,jj) = zCdn_skin_ice * zfmi + zCdn_form_ice * ( zfmi*zfi ) * ztmp
            pch(ji,jj) = zChn_skin_ice * zfhi + zChn_form_ice * ( zfhi*zfi ) * ztmp

         END DO
      END DO
      !
   END SUBROUTINE Cx_Lupkes2015

   !!======================================================================

END MODULE mod_blk_ice_lu15
