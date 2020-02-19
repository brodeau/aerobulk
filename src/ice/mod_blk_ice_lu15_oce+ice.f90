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

   PUBLIC :: turb_ice_lu15, Cdn10_Lupkes2015

   ! ECHAM6 constants
   REAL(wp), PARAMETER ::   z0_skin_ice  = 0.69e-3_wp  ! Eq. 43 [m]
   REAL(wp), PARAMETER ::   z0_form_ice  = 0.57e-3_wp  ! Eq. 42 [m]
   REAL(wp), PARAMETER ::   z0_ice       = 1.00e-3_wp  ! Eq. 15 [m]
   REAL(wp), PARAMETER ::   zce10        = 2.80e-3_wp  ! Eq. 41
   REAL(wp), PARAMETER ::   zbeta        = 1.1_wp      ! Eq. 41
   REAL(wp), PARAMETER ::   zc           = 5._wp       ! Eq. 13
   REAL(wp), PARAMETER ::   zc2          = zc * zc
   REAL(wp), PARAMETER ::   zam          = 2. * zc     ! Eq. 14
   REAL(wp), PARAMETER ::   zah          = 3. * zc     ! Eq. 30
   REAL(wp), PARAMETER ::   z1_alpha     = 1._wp / 0.2_wp  ! Eq. 51
   REAL(wp), PARAMETER ::   z1_alphaf    = z1_alpha    ! Eq. 56
   REAL(wp), PARAMETER ::   zbetah       = 1.e-3_wp    ! Eq. 26
   REAL(wp), PARAMETER ::   zgamma       = 1.25_wp     ! Eq. 26
   REAL(wp), PARAMETER ::   z1_gamma     = 1._wp / zgamma
   REAL(wp), PARAMETER ::   r1_3         = 1._wp / 3._wp

   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE turb_ice_lu15( kt, zt, zu, Ti_s, t_zt, qi_s, q_zt, U_zu,   &
      &                     Cd, Ch, Ce, t_zu, q_zu, U_blk,        &
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
      !!
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

      t_zu = MAX( t_zt ,   100._wp )   ! who knows what's given on masked-continental regions...
      q_zu = MAX( q_zt , 0.1e-6_wp )   !               "

      !! Scalar wind speed cannot be below 0.2 m/s
      U_blk = MAX( U_zu, 0.2_wp )

      !! Pot. temp. difference (and we don't want it to be 0!)
      dt_zu = t_zu - Ti_s ;   dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6_wp), dt_zu )
      dq_zu = q_zu - qi_s ;   dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9_wp), dq_zu )

      Cd = rCd_ice
      Ch = rCd_ice
      Ce = rCd_ice

      u_star = SQRT(rCd_ice) * U_blk
      t_star = rCd_ice * U_blk * dt_zu / u_star
      q_star = rCd_ice * U_blk * dq_zu / u_star


      IF( lreturn_z0 )    xz0     = EXP( LOG(10._wp) - vkarmn/SQRT(rCd_ice) )
      IF( lreturn_ustar ) xu_star = u_star
      IF( lreturn_L )     xL      = 1./One_on_L(t_zu, q_zu, u_star, t_star, q_star)
      !IF( lreturn_UN10 )  xUN10   = u_star/vkarmn*LOG(10./z0)
      IF( lreturn_UN10 )  xUN10   = U_blk

      DEALLOCATE ( u_star, t_star, q_star, dt_zu, dq_zu )

   END SUBROUTINE turb_ice_lu15





   SUBROUTINE Cdn10_Lupkes2015( zu, Ts_o, t_zu, q_zu, Uo_zu, Ui_zu, Cdn_o, Chn_o, &
      &                         pic, Ts_i, pslp, pcd, pch )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  Cdn10_Lupkes2015  ***
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
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: Ts_o   ! sea surface temperature  [Kelvin]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: t_zu   ! potential air temperature              [Kelvin]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: q_zu   ! specific air humidity at zt             [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: Uo_zu  ! relative wind module at zu over ocean  [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: Ui_zu  ! relative wind module at zu over ice    [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: Cdn_o  ! neutral drag coefficient over ocean
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: Chn_o  ! neutral heat coefficient over ocean
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pic    ! ice concentration [fraction]  => at_i_b
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: Ts_i ! sea-ice surface temperature [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pslp   ! sea-level pressure [Pa]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pcd    ! momentum transfer coefficient
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pch    ! heat transfer coefficient
      !
      REAL(wp) ::   zfi, zfo, zqsat_o, zqsat_i
      REAL(wp) ::   zrib_o, zrib_i, ztmp
      REAL(wp) ::   zCdn_skin_ice, zCdn_form_ice, zCdn_ice
      REAL(wp) ::   zChn_skin_ice, zChn_form_ice
      REAL(wp) ::   z0w, z0i, zfmi, zfmw, zfhi, zfhw
      REAL(wp) ::   zCdn_form_tmp, zwndspd_o, zwndspd_i
      INTEGER  ::   ji, jj         ! dummy loop indices
      !!----------------------------------------------------------------------

      ! Momentum Neutral Transfer Coefficients (should be a constant)
      zCdn_form_tmp = zce10 * ( LOG( 10._wp / z0_form_ice + 1._wp ) / LOG( zu / z0_form_ice + 1._wp ) )**2   ! Eq. 40
      zCdn_skin_ice = ( vkarmn                                      / LOG( zu / z0_skin_ice + 1._wp ) )**2   ! Eq. 7
      zCdn_ice      = zCdn_skin_ice   ! Eq. 7
      !zCdn_ice     = 1.89e-3         ! old ECHAM5 value (cf Eq. 32)

      ! Heat Neutral Transfer Coefficients
      zChn_skin_ice = vkarmn**2 / ( LOG( zu / z0_ice + 1._wp ) * LOG( zu * z1_alpha / z0_skin_ice + 1._wp ) )   ! Eq. 50 + Eq. 52

      DO jj = 1, jpj
         DO ji = 1, jpi

            zfi = pic(ji,jj)    ! fraction of sea-ice
            zfo = 1._wp - zfi    ! fraction of open ocean

            zwndspd_o = MAX( 0.5, Uo_zu(ji,jj) )
            zwndspd_i = MAX( 0.5, Ui_zu(ji,jj) )

            ! Specific humidities at saturation at air-sea and air-ice interface [kg/kg]:
            zqsat_o = rdct_qsat_salt*q_sat( Ts_o(ji,jj), pslp(ji,jj) )
            zqsat_i =                q_sat( Ts_i(ji,jj), pslp(ji,jj) ) !LOLO!!! No!!! Must use a special 1 over ice!

            ! Bulk Richardson Number:
            zrib_o = Ri_bulk( zu, Ts_o(ji,jj), t_zu(ji,jj), zqsat_o, q_zu(ji,jj), zwndspd_o ) !LOLO: is t_zu potential T? It should!!!!
            zrib_i = Ri_bulk( zu, Ts_i(ji,jj), t_zu(ji,jj), zqsat_i, q_zu(ji,jj), zwndspd_i )

            ! Momentum and Heat Neutral Transfer Coefficients
            zCdn_form_ice = zCdn_form_tmp * zfi * zfo**zbeta                          ! Eq. 40
            zChn_form_ice = zCdn_form_ice / ( 1._wp + ( LOG( z1_alphaf ) / vkarmn ) * SQRT( zCdn_form_ice ) )   ! Eq. 53

            ! Momentum and Heat Stability functions (possibility to use psi_m_ecmwf instead ?)
            z0w = zu * EXP( -1._wp * vkarmn / SQRT( Cdn_o(ji,jj) ) ) ! over water
            z0i = z0_skin_ice                                        ! over ice

            IF( zrib_o <= 0._wp ) THEN
               ztmp = virt_temp( Ts_o(ji,jj), zqsat_o ) - virt_temp( t_zu(ji,jj), q_zu(ji,jj) ) ! difference of potential virtual temperature
               zfmw = 1._wp - zam*zrib_o / ( 1._wp + 3._wp*zc2*Cdn_o(ji,jj)*SQRT( -zrib_o*( zu / z0w + 1._wp ) ) )  ! Eq. 10
               zfhw = ( 1._wp + ( zbetah*ztmp**r1_3 / ( Chn_o(ji,jj) * zwndspd_o ) )**zgamma )**z1_gamma      ! Eq. 26
            ELSE
               ztmp = zrib_o / SQRT( 1._wp + zrib_o )
               zfmw = 1._wp / ( 1._wp + zam * ztmp )   ! Eq. 12
               zfhw = 1._wp / ( 1._wp + zah * ztmp )   ! Eq. 28
            ENDIF

            IF( zrib_i <= 0._wp ) THEN
               ztmp = zrib_i / ( 1._wp + 3._wp * zc2 * zCdn_ice * SQRT( -zrib_i * ( zu / z0i + 1._wp) ) )
               zfmi = 1._wp - zam * ztmp   ! Eq.  9
               zfhi = 1._wp - zah * ztmp   ! Eq. 25
            ELSE
               ztmp = zrib_i / SQRT( 1._wp + zrib_i )
               zfmi = 1._wp / ( 1._wp + zam * ztmp )   ! Eq. 11
               zfhi = 1._wp / ( 1._wp + zah * ztmp )   ! Eq. 27
            ENDIF



            ! Momentum and Heat transfer coefficients (Eq. 38) and (Eq. 49):
            ztmp       = 1._wp / MAX( 1.e-06, zfi )
            pcd(ji,jj) = zCdn_skin_ice * zfmi + zCdn_form_ice * ( zfmi*zfi + zfmw*zfo ) * ztmp
            pch(ji,jj) = zChn_skin_ice * zfhi + zChn_form_ice * ( zfhi*zfi + zfhw*zfo ) * ztmp

         END DO
      END DO
      !
      !CALL lbc_lnk_multi( 'sbcblk', pcd, 'T',  1., pch, 'T', 1. )
      !
   END SUBROUTINE Cdn10_Lupkes2015

   !!======================================================================

END MODULE mod_blk_ice_lu15
