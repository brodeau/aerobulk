! AeroBulk / 2020 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_blk_ice_best
   !!====================================================================================
   !!       Computes turbulent components of surface fluxes over sea-ice
   !!
   !!
   !!       Routine turb_ice_best maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, February 2020
   !!
   !!====================================================================================
   USE mod_const       !: physical and other constants
   USE mod_phymbl      !: thermodynamics

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: TURB_ICE_BEST

   ! ECHAM6 constants
   REAL(wp), PARAMETER ::   z0_skin_ice  = 0.69e-3_wp  ! Eq. 43 [m]
   REAL(wp), PARAMETER ::   z0_form_ice  = 0.57e-3_wp  ! Eq. 42 [m]
   REAL(wp), PARAMETER ::   z0_ice       = 1.00e-3_wp  ! Eq. 15 [m]
   REAL(wp), PARAMETER ::   zce10        = 2.80e-3_wp  ! Eq. 41
   REAL(wp), PARAMETER ::   zbeta        = 1.1_wp      ! Eq. 41
   REAL(wp), PARAMETER ::   zc           = 5._wp       ! Eq. 13
   REAL(wp), PARAMETER ::   z1_alpha     = 1._wp / 0.2_wp  ! Eq. 51
   REAL(wp), PARAMETER ::   z1_alphaf    = z1_alpha    ! Eq. 56
   REAL(wp), PARAMETER ::   zbetah       = 1.e-3_wp    ! Eq. 26
   REAL(wp), PARAMETER ::   zgamma       = 1.25_wp     ! Eq. 26
   REAL(wp), PARAMETER ::   z1_gamma     = 1._wp / zgamma
   REAL(wp), PARAMETER ::   r1_3         = 1._wp / 3._wp

   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE turb_ice_best( kt, zt, zu, Ts_i, t_zt, qi_s, q_zt, U_zu, &
      &                     Cd, Ch, Ce, t_zu, q_zu, U_blk,             &
      &                     xz0, xu_star, xL, xUN10 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_ice_best  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to
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
      !!    *  qi_s  : SSQ aka saturation specific humidity at temp. Ts_i     [kg/kg]
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
      !! ** Author: L. Brodeau, February 2020 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      INTEGER,  INTENT(in )                     ::   kt       ! current time step
      REAL(wp), INTENT(in )                     ::   zt       ! height for t_zt and q_zt                    [m]
      REAL(wp), INTENT(in )                     ::   zu       ! height for U_zu                             [m]
      REAL(wp), INTENT(in ), DIMENSION(:,:) ::   Ts_i     ! ice surface temperature                [Kelvin]
      REAL(wp), INTENT(in ), DIMENSION(:,:) ::   t_zt     ! potential air temperature              [Kelvin]
      REAL(wp), INTENT(in ), DIMENSION(:,:) ::   qi_s     ! specific humidity at ice/air interface  [kg/kg]
      REAL(wp), INTENT(in ), DIMENSION(:,:) ::   q_zt     ! specific air humidity at zt             [kg/kg]
      REAL(wp), INTENT(in ), DIMENSION(:,:) ::   U_zu     ! relative wind module at zu                [m/s]
      REAL(wp), INTENT(out), DIMENSION(:,:) ::   Cd       ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(out), DIMENSION(:,:) ::   Ch       ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(out), DIMENSION(:,:) ::   Ce       ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), INTENT(out), DIMENSION(:,:) ::   t_zu     ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(out), DIMENSION(:,:) ::   q_zu     ! spec. humidity adjusted at zu           [kg/kg]
      REAL(wp), INTENT(out), DIMENSION(:,:) ::   U_blk    ! bulk wind speed at zu                     [m/s]
      !
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(:,:) ::   xz0  ! Aerodynamic roughness length   [m]
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(:,:) ::   xu_star  ! u*, friction velocity
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(:,:) ::   xL  ! zeta (zu/L)
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(:,:) ::   xUN10  ! Neutral wind at zu
      !
      INTEGER :: Ni, Nj, jit
      LOGICAL :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   Cx_n10        ! 10m neutral latent/sensible coefficient
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   sqrtCdn10   ! square root of Cd_n10
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   zeta_u        ! stability parameter at height zu
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   ztmp0, ztmp1, ztmp2
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   sqrtCd       ! square root of Cd
      !
      LOGICAL :: lreturn_z0=.FALSE., lreturn_ustar=.FALSE., lreturn_L=.FALSE., lreturn_UN10=.FALSE.
      CHARACTER(len=40), PARAMETER :: crtnm = 'turb_ice_best@mod_blk_ice_best.f90'
      !!----------------------------------------------------------------------------------
      Ni = SIZE(Ts_i,1)
      Nj = SIZE(Ts_i,2)
      ALLOCATE( Cx_n10(Ni,Nj), sqrtCdn10(Ni,Nj), &
         &    zeta_u(Ni,Nj), sqrtCd(Ni,Nj),      &
         &    ztmp0(Ni,Nj),  ztmp1(Ni,Nj), ztmp2(Ni,Nj) )

      IF( PRESENT(xz0) )     lreturn_z0    = .TRUE.
      IF( PRESENT(xu_star) ) lreturn_ustar = .TRUE.
      IF( PRESENT(xL) )      lreturn_L     = .TRUE.
      IF( PRESENT(xUN10) )   lreturn_UN10  = .TRUE.

      l_zt_equal_zu = .FALSE.
      IF( ABS(zu - zt) < 0.01_wp )   l_zt_equal_zu = .TRUE.    ! testing "zu == zt" is risky with double precision

      !! Scalar wind speed cannot be below 0.2 m/s
      U_blk = MAX( U_zu, wspd_thrshld_ice )

      !! Initializing values at z_u with z_t values:
      t_zu = t_zt
      q_zu = q_zt

      CALL Cx_Lupkes2015( zu, t_zu, q_zu, U_blk, Ts_i, qi_s, Cd, Ch )
      Ce = Ch
      sqrtCd = SQRT( Cd )
      !LOLO:STOP


      !! ITERATION BLOCK
      DO jit = 1, nb_iter
         !
         ztmp1 = t_zu - Ts_i   ! Updating air/sea differences
         ztmp2 = q_zu - qi_s

         ! Updating turbulent scales :   (L&Y 2004 Eq. (7))
         ztmp0 = sqrtCd*U_blk       ! u*
         ztmp1 = Ch/sqrtCd*ztmp1    ! theta*
         ztmp2 = Ce/sqrtCd*ztmp2    ! q*

         ! Estimate the inverse of Obukov length (1/L) at height zu:
         ztmp0 = One_on_L( t_zu, q_zu, ztmp0, ztmp1, ztmp2 )
         !ztmp0 = One_on_L( 0.5*(t_zu + Ts_i), 0.5*(q_zu + qi_s), ztmp0, ztmp1, ztmp2 ) ! using an approximation of mean
         !                                                                ! theta & q in surface layer rather than values at zu...

         !! Stability parameters :
         zeta_u   = zu*ztmp0
         zeta_u   = sign( min(abs(zeta_u),10._wp), zeta_u )

         !! Shifting temperature and humidity at zu (L&Y 2004 Eq. (9b-9c))
         IF( .NOT. l_zt_equal_zu ) THEN
            ztmp0 = zt*ztmp0 ! zeta_t !
            ztmp0 = SIGN( MIN(ABS(ztmp0),10._wp), ztmp0 )  ! Temporaty array ztmp0 == zeta_t !!!
            ztmp0 = LOG(zt/zu) + psi_h_ice(zeta_u) - psi_h_ice(ztmp0)                   ! ztmp0 just used as temp array again!
            t_zu = t_zt - ztmp1/vkarmn*ztmp0    ! ztmp1 is still theta*  L&Y 2004 Eq. (9b)
            q_zu = q_zt - ztmp2/vkarmn*ztmp0    ! ztmp2 is still q*      L&Y 2004 Eq. (9c)
            q_zu = MAX(0._wp, q_zu)
         END IF

         ! Update neutral wind speed at 10m and neutral Cd at 10m (L&Y 2004 Eq. 9a)...
         !   In very rare low-wind conditions, the old way of estimating the
         !   neutral wind speed at 10m leads to a negative value that causes the code
         !   to crash. To prevent this a threshold of 0.25m/s is imposed.
         ztmp2 = psi_m_ice(zeta_u)
         ztmp0 = MAX( wspd_thrshld_ice , U_blk/(1._wp + sqrtCdn10/vkarmn*(LOG(zu/10._wp) - ztmp2)) ) ! U_n10 (ztmp2 == psi_m_ice(zeta_u))

         CALL Cx_Lupkes2015( zu, t_zu, q_zu, ztmp0, Ts_i, qi_s, Cd, Cx_n10 )
         sqrtCdn10 = sqrt(Cd)

         !! Update of transfer coefficients:
         ztmp1  = 1._wp + sqrtCdn10/vkarmn*(LOG(zu/10._wp) - ztmp2)   ! L&Y 2004 Eq. (10a) (ztmp2 == psi_m_ice(zeta_u))
         Cd     = Cd / ( ztmp1*ztmp1 )
         sqrtCd = SQRT( Cd )

         ztmp0  = ( LOG(zu/10._wp) - psi_h_ice(zeta_u) ) / vkarmn / sqrtCdn10
         ztmp2  = sqrtCd / sqrtCdn10
         ztmp1  = 1._wp + Cx_n10*ztmp0
         Ch     = Cx_n10 * ztmp2 / ztmp1
         Ce     = Ch

      END DO !DO jit = 1, nb_iter

      IF( lreturn_z0 )    xz0     = z0_from_Cd( zu, Cd,  ppsi=psi_m_ice(zeta_u) )
      IF( lreturn_ustar ) xu_star = SQRT( Cd )*U_blk
      IF( lreturn_L )     xL      = zu/zeta_u
      IF( lreturn_UN10 )  xUN10   = U_blk/(1. + sqrtCdn10/vkarmn*(LOG(zu/10.) - psi_m_ice(zeta_u)))

      DEALLOCATE( Cx_n10, sqrtCdn10, zeta_u, sqrtCd, ztmp0, ztmp1, ztmp2 ) !

   END SUBROUTINE turb_ice_best



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
      REAL(wp),                 INTENT(in)  :: zu     ! reference height (height for Uo_zu)   [m]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: t_zu   ! potential air temperature              [Kelvin]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: q_zu   ! specific air humidity at zt             [kg/kg]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: Ui_zu  ! relative wind module at zu over ice    [m/s]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: Ts_i   ! sea-ice surface temperature               [K]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: qs_i   ! humidity at saturation over ice at T=Ts_i [kg/kg]
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pcd    ! momentum transfer coefficient
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pch    ! heat transfer coefficient
      !
      REAL(wp) ::   zfi, zfo
      REAL(wp) ::   zrib_i, ztmp
      REAL(wp) ::   zCdn_skin_ice, zCdn_form_ice, zCdn_ice
      REAL(wp) ::   zChn_skin_ice, zChn_form_ice
      REAL(wp) ::   z0i, zfmi, zfhi
      REAL(wp) ::   zCdn_form_tmp, zwndspd_i
      INTEGER  ::   Ni, Nj, ji, jj         ! dummy loop indices
      !!----------------------------------------------------------------------
      Ni = SIZE(Ts_i,1)
      Nj = SIZE(Ts_i,2)
      
      ! Momentum Neutral Transfer Coefficients (should be a constant)
      zCdn_form_tmp = zce10 * ( LOG( 10._wp / z0_form_ice + 1._wp ) / LOG( zu / z0_form_ice + 1._wp ) )**2   ! Eq. 46
      zCdn_skin_ice = ( vkarmn                                      / LOG( zu / z0_skin_ice + 1._wp ) )**2   ! Eq. 7
      zCdn_ice      = zCdn_skin_ice   ! Eq. 7
      !zCdn_ice     = 1.89e-3         ! old ECHAM5 value (cf Eq. 32)

      ! Heat Neutral Transfer Coefficients
      zChn_skin_ice = vkarmn**2 / ( LOG( zu / z0_ice + 1._wp ) * LOG( zu * z1_alpha / z0_skin_ice + 1._wp ) )   ! Eq. 50 + Eq. 52

      DO jj = 1, Nj
         DO ji = 1, Ni

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


   FUNCTION psi_m_ice( pzeta )
      !!----------------------------------------------------------------------------------
      !! ** Purpose: compute the universal profile stability function for momentum
      !!
      !!
      !!     Andreas et al 2005 == Jordan et al. 1999
      !!
      !!     Psi:
      !!     Unstable => Paulson 1970
      !!     Stable   => Holtslag & De Bruin 1988
      !!
      !!             pzeta : stability paramenter, z/L where z is altitude
      !!                     measurement and L is M-O length
      !!
      !! ** Author: L. Brodeau, 2020 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in) :: pzeta
      REAL(wp), DIMENSION(SIZE(pzeta,1),SIZE(pzeta,2)) :: psi_m_ice
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) :: zta, zx, zpsi_u, zpsi_s, zstab
      !!----------------------------------------------------------------------------------
      DO jj = 1, SIZE(pzeta,2)
         DO ji = 1, SIZE(pzeta,1)
            !
            zta = pzeta(ji,jj)
            !
            ! Unstable stratification:
            zx = ABS(1._wp - 16._wp*zta)**.25              !  (16 here, not 15!)

            zpsi_u =      LOG( (1._wp + zx*zx)/2. ) &  ! Eq.(30) Jordan et al. 1999
               &     + 2.*LOG( (1._wp + zx       )/2. ) &
               &    - 2.*ATAN( zx ) + 0.5*rpi

            ! Stable stratification:
            zpsi_s = - ( 0.7_wp*zta + 0.75_wp*(zta - 14.3_wp)*EXP( -0.35*zta) + 10.7_wp )  ! Eq.(33) Jordan et al. 1999

            !! Combine:
            zstab = 0.5 + SIGN(0.5_wp, zta)
            psi_m_ice(ji,jj) = (1._wp - zstab) * zpsi_u   & ! Unstable (zta < 0)
               &                   + zstab     * zpsi_s     ! Stable (zta > 0)
            !
         END DO
      END DO
   END FUNCTION psi_m_ice


   FUNCTION psi_h_ice( pzeta )
      !!----------------------------------------------------------------------------------
      !! ** Purpose: compute the universal profile stability function for
      !!             temperature and humidity
      !!
      !!
      !!     Andreas et al 2005 == Jordan et al. 1999
      !!
      !!     Psi:
      !!     Unstable => Paulson 1970
      !!     Stable   => Holtslag & De Bruin 1988
      !!
      !!             pzeta : stability paramenter, z/L where z is altitude
      !!                     measurement and L is M-O length
      !!
      !! ** Author: L. Brodeau, 2020 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in) :: pzeta
      REAL(wp), DIMENSION(SIZE(pzeta,1),SIZE(pzeta,2)) :: psi_h_ice
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) :: zta, zx, zpsi_u, zpsi_s, zstab
      !!----------------------------------------------------------------------------------
      DO jj = 1, SIZE(pzeta,2)
         DO ji = 1, SIZE(pzeta,1)
            !
            zta = pzeta(ji,jj)
            !
            ! Unstable stratification:
            zx = ABS(1._wp - 16._wp*zta)**.25              !  (16 here, not 15!)

            zpsi_u =   2.*LOG( (1._wp + zx*zx)/2. )  ! Eq.(31) Jordan et al. 1999

            ! Stable stratification (identical to psi_m_ice!):
            zpsi_s = - ( 0.7_wp*zta + 0.75_wp*(zta - 14.3_wp)*EXP( -0.35*zta) + 10.7_wp )  ! Eq.(33) Jordan et al. 1999

            !! Combine:
            zstab = 0.5 + SIGN(0.5_wp, zta)
            psi_h_ice(ji,jj) = (1._wp - zstab) * zpsi_u   & ! Unstable (zta < 0)
               &                   + zstab     * zpsi_s     ! Stable (zta > 0)
            !
         END DO
      END DO
   END FUNCTION psi_h_ice

   !!======================================================================
END MODULE mod_blk_ice_best
