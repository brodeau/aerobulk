! AeroBulk / 2020 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_blk_ice_an05
   !!====================================================================================
   !!       Computes turbulent components of surface fluxes over sea-ice
   !!
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (usually 2m) to zu (usually 10m) if needed
   !!   * the "effective" bulk wind speed at zu: Ub (including gustiness contribution in unstable conditions)
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!       Routine turb_ice_an05 maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, January 2020
   !!
   !!====================================================================================
   USE mod_const       !: physical and other constants
   USE mod_phymbl      !: misc. physical functions

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: TURB_ICE_AN05
   PUBLIC :: rough_leng_m, rough_leng_tq

   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE turb_ice_an05( kt, zt, zu, Ts_i, t_zt, qs_i, q_zt, U_zu,         &
      &                      Cd, Ch, Ce, t_zu, q_zu, Ub,                       &
      &                      CdN, ChN, CeN, xz0, xu_star, xL, xUN10 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_ice_an05  ***
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
      !!
      !! OUTPUT :
      !! --------
      !!    *  Cd     : drag coefficient
      !!    *  Ch     : sensible heat coefficient
      !!    *  Ce     : evaporation coefficient
      !!    *  t_zu   : pot. air temperature adjusted at wind height zu       [K]
      !!    *  q_zu   : specific humidity of air        //                    [kg/kg]
      !!    *  Ub     : bulk wind speed at zu that we used                    [m/s]
      !!
      !! OPTIONAL OUTPUT:
      !! ----------------
      !!    * CdN      : neutral-stability drag coefficient
      !!    * ChN      : neutral-stability sensible heat coefficient
      !!    * CeN      : neutral-stability evaporation coefficient
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
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Cd    ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Ch    ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Ce    ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: t_zu  ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: q_zu  ! spec. humidity adjusted at zu           [kg/kg]
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Ub ! bulk wind speed at zu                     [m/s]
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   CdN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   ChN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   CeN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xz0  ! Aerodynamic roughness length   [m]
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xu_star  ! u*, friction velocity
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xL  ! zeta (zu/L)
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xUN10  ! Neutral wind at zu
      !!
      INTEGER :: j_itt
      LOGICAL :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !
      REAL(wp), DIMENSION(:,:), ALLOCATABLE  ::  &
         &  u_star, t_star, q_star, &
         &  dt_zu, dq_zu,    &
         &  znu_a,           & !: Nu_air = kinematic viscosity of air
         &  z0
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE  :: z0tq

      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zeta_u        ! stability parameter at height zu
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zeta_t        ! stability parameter at height zt
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: ztmp0, ztmp1, ztmp2
      !
      LOGICAL ::  lreturn_cdn=.FALSE., lreturn_chn=.FALSE., lreturn_cen=.FALSE., &
         &        lreturn_z0=.FALSE., lreturn_ustar=.FALSE., lreturn_L=.FALSE., lreturn_UN10=.FALSE.
      !!
      CHARACTER(len=40), PARAMETER :: crtnm = 'turb_ice_an05@mod_blk_ice_an05.f90'
      !!----------------------------------------------------------------------------------
      ALLOCATE ( u_star(jpi,jpj), t_star(jpi,jpj),  q_star(jpi,jpj),  &
         &       zeta_u(jpi,jpj),  dt_zu(jpi,jpj),   dq_zu(jpi,jpj),  &
         &        znu_a(jpi,jpj),  ztmp1(jpi,jpj),   ztmp2(jpi,jpj),  &
         &           z0(jpi,jpj),   z0tq(jpi,jpj,2), ztmp0(jpi,jpj)   )

      IF( PRESENT(CdN) )     lreturn_cdn   = .TRUE.
      IF( PRESENT(ChN) )     lreturn_chn   = .TRUE.
      IF( PRESENT(CeN) )     lreturn_cen   = .TRUE.
      IF( PRESENT(xz0) )     lreturn_z0    = .TRUE.
      IF( PRESENT(xu_star) ) lreturn_ustar = .TRUE.
      IF( PRESENT(xL) )      lreturn_L     = .TRUE.
      IF( PRESENT(xUN10) )   lreturn_UN10  = .TRUE.

      l_zt_equal_zu = .FALSE.
      IF( ABS(zu - zt) < 0.01_wp )   l_zt_equal_zu = .TRUE.    ! testing "zu == zt" is risky with double precision
      IF( .NOT. l_zt_equal_zu )  ALLOCATE( zeta_t(jpi,jpj) )


      !! Scalar wind speed cannot be below 0.2 m/s
      Ub = MAX( U_zu, wspd_thrshld_ice )

      !! First guess of temperature and humidity at height zu:
      t_zu = MAX( t_zt ,   100._wp )   ! who knows what's given on masked-continental regions...
      q_zu = MAX( q_zt , 0.1e-6_wp )   !               "

      !! Air-Ice differences (and we don't want it to be 0!)
      dt_zu = t_zu - Ts_i ;   dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6_wp), dt_zu )
      dq_zu = q_zu - qs_i ;   dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9_wp), dq_zu )

      znu_a = visc_air(t_zu) ! Air viscosity (m^2/s) at zt given from temperature in (K)

      !! Very crude first guesses of z0:
      z0 = 8.0E-4_wp

      !! Crude first guess of turbulent scales
      u_star = 0.035_wp*Ub*LOG( 10._wp/z0 )/LOG( zu/z0 )
      z0 = rough_leng_m( u_star , znu_a )

      DO j_itt = 1, 2
         u_star = MAX ( Ub*vkarmn/(LOG(zu) - LOG(z0)) , 1.E-9 )
         z0 = rough_leng_m( u_star , znu_a )
      END DO

      z0tq = rough_leng_tq( z0, u_star , znu_a )
      t_star = dt_zu*vkarmn/(LOG(zu/z0tq(:,:,1)))
      q_star = dq_zu*vkarmn/(LOG(zu/z0tq(:,:,2)))



      !! ITERATION BLOCK
      DO j_itt = 1, nb_itt

         !!Inverse of Obukov length (1/L) :
         ztmp0 = One_on_L(t_zu, q_zu, u_star, t_star, q_star)  ! 1/L == 1/[Obukhov length]
         ztmp0 = SIGN( MIN(ABS(ztmp0),200._wp), ztmp0 ) ! (prevents FPE from stupid values from masked region later on...)

         !! Stability parameters "zeta" :
         zeta_u = zu*ztmp0
         zeta_u = SIGN( MIN(ABS(zeta_u),50.0_wp), zeta_u )
         IF( .NOT. l_zt_equal_zu ) THEN
            zeta_t = zt*ztmp0
            zeta_t = SIGN( MIN(ABS(zeta_t),50.0_wp), zeta_t )
         END IF

         !! Roughness lengthes z0, z0t, & z0q :
         z0   = rough_leng_m (     u_star , znu_a )
         z0tq = rough_leng_tq( z0, u_star , znu_a )

         !! Turbulent scales at zu :
         ztmp0   = psi_h_ice(zeta_u)
         t_star =      dt_zu*vkarmn/(LOG(zu) - LOG(z0tq(:,:,1)) - ztmp0)
         q_star =      dq_zu*vkarmn/(LOG(zu) - LOG(z0tq(:,:,2)) - ztmp0)
         u_star = MAX( Ub*vkarmn/(LOG(zu) - LOG(z0(:,:)) - psi_m_ice(zeta_u)) , 1.E-9 )

         IF( .NOT. l_zt_equal_zu ) THEN
            !! Re-updating temperature and humidity at zu if zt /= zu :
            ztmp1 = LOG(zt/zu) + ztmp0 - psi_h_ice(zeta_t)
            t_zu = t_zt - t_star/vkarmn*ztmp1
            q_zu = q_zt - q_star/vkarmn*ztmp1
            dt_zu = t_zu - Ts_i ;  dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6_wp), dt_zu )
            dq_zu = q_zu - qs_i ;  dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9_wp), dq_zu )
         END IF

      END DO !DO j_itt = 1, nb_itt

      ! compute transfer coefficients at zu :
      ztmp0 = u_star/Ub
      Cd   = ztmp0*ztmp0
      Ch   = ztmp0*t_star/dt_zu
      Ce   = ztmp0*q_star/dq_zu

      IF( lreturn_cdn .OR. lreturn_chn .OR. lreturn_cen ) ztmp0 = 1._wp/LOG( zu/z0(:,:) )
      IF( lreturn_cdn )   CdN = vkarmn2*ztmp0*ztmp0
      IF( lreturn_chn )   ChN = vkarmn2*ztmp0/LOG(zu/z0tq(:,:,1))
      IF( lreturn_cen )   CeN = vkarmn2*ztmp0/LOG(zu/z0tq(:,:,2))

      IF( lreturn_z0 )    xz0     = z0
      IF( lreturn_ustar ) xu_star = u_star
      IF( lreturn_L )     xL      = 1./One_on_L(t_zu, q_zu, u_star, t_star, q_star)
      IF( lreturn_UN10 )  xUN10   = u_star/vkarmn*LOG(10./z0)

      DEALLOCATE ( u_star, t_star, q_star, zeta_u, dt_zu, dq_zu, z0, z0tq, znu_a, ztmp0, ztmp1, ztmp2 )
      IF( .NOT. l_zt_equal_zu ) DEALLOCATE( zeta_t )

   END SUBROUTINE turb_ice_an05



   FUNCTION rough_leng_m( pus , pnua )
      !!----------------------------------------------------------------------------------
      !! Computes the roughness length of sea-ice according to Andreas et al. 2005, (eq. 19)
      !!
      !! Author: L. Brodeau, January 2020 / AeroBulk  (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: rough_leng_m      ! roughness length over sea-ice [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pus   ! u* = friction velocity    [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pnua  ! kinematic viscosity of air [m^2/s]
      !!
      INTEGER  :: ji, jj    ! dummy loop indices
      REAL(wp) :: zus, zz
      !!----------------------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi

            zus = MAX( pus(ji,jj) , 1.E-9_wp )

            zz = (zus - 0.18_wp) / 0.1_wp

            rough_leng_m(ji,jj) = 0.135*pnua(ji,jj)/zus + 0.035*zus*zus/grav*( 5.*EXP(-zz*zz) + 1._wp ) ! Eq.(19) Andreas et al., 2005

         END DO
      END DO
      !!
   END FUNCTION rough_leng_m

   FUNCTION rough_leng_tq( pz0, pus , pnua )
      !!----------------------------------------------------------------------------------
      !! Computes the roughness length of sea-ice according to Andreas et al. 2005, (eq. 22)
      !!    => which still relies on Andreas 1987 !
      !!
      !! Author: L. Brodeau, January 2020 / AeroBulk  (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,2)           :: rough_leng_tq     ! temp.,hum. roughness lengthes over sea-ice [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pz0   ! roughness length            [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pus   ! u* = friction velocity    [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pnua  ! kinematic viscosity of air [m^2/s]
      !!
      INTEGER  :: ji, jj    ! dummy loop indices
      REAL(wp) :: zz0, zus, zre, zsmoot, ztrans, zrough
      REAL(wp) :: zb0, zb1, zb2, zlog, zlog2, zlog_z0s_on_z0
      !!----------------------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi

            zz0 = pz0(ji,jj)
            zus = MAX( pus(ji,jj) , 1.E-9_wp )
            zre = MAX( zus*zz0/pnua(ji,jj) , 0._wp ) ! Roughness Reynolds number

            !! *** TABLE 1 of Andreas et al. 2005 ***
            !! Smooth flow condition (R* <= 0.135):
            zsmoot = 0.5_wp + SIGN( 0.5_wp, (0.135_wp   - zre) ) ! zre <= 0.135: zsmoot==1 ; otherwize: zsmoot==0
            !! Transition (0.135 < R* < 2.5):
            ztrans = 0.5_wp + SIGN( 0.5_wp, (2.49999_wp - zre) ) - zsmoot
            !! Rough ( R* > 2.5):
            zrough = 0.5_wp + SIGN( 0.5_wp, (zre - 2.5_wp) )

            IF( (zsmoot+ztrans+zrough > 1.001_wp).OR.(zsmoot+ztrans+zrough < 0.999_wp) ) &
               CALL ctl_stop( ' rough_leng_tq@mod_blk_ice_an05.f90 => something wrong with zsmoot, ztrans, zrough!' )

            zlog  = LOG(zre)
            zlog2 = zlog*zlog

            !! z0t:
            zb0 = zsmoot*1.25_wp + ztrans*0.149_wp + zrough*0.317_wp
            zb1 =                - ztrans*0.550_wp - zrough*0.565_wp
            zb2 =                                  - zrough*0.183_wp
            zlog_z0s_on_z0 = zb0 + zb1*zlog + zb2*zlog2
            rough_leng_tq(ji,jj,1) = zz0 * EXP( zlog_z0s_on_z0 )

            !! z0q:
            zb0 = zsmoot*1.61_wp + ztrans*0.351_wp + zrough*0.396_wp
            zb1 =                - ztrans*0.628_wp - zrough*0.512_wp
            zb2 =                                  - zrough*0.180_wp
            zlog = LOG(zre)
            zlog_z0s_on_z0 = zb0 + zb1*zlog + zb2*zlog2
            rough_leng_tq(ji,jj,2) = zz0 * EXP( zlog_z0s_on_z0 )

         END DO
      END DO
      !!
   END FUNCTION rough_leng_tq









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
      REAL(wp), DIMENSION(jpi,jpj) :: psi_m_ice
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) :: zta, zx, zpsi_u, zpsi_s, zstab
      !!----------------------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
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
      REAL(wp), DIMENSION(jpi,jpj) :: psi_h_ice
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) :: zta, zx, zpsi_u, zpsi_s, zstab
      !!----------------------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zta = pzeta(ji,jj)
            !
            ! Unstable stratification:
            zx = ABS(1._wp - 16._wp*zta)**.25              !  (16 here, not 15!)

            zpsi_u =   2.*LOG( (1._wp + zx*zx)/2. )  ! Eq.(31) Jordan et al. 1999

            ! Stable stratification (identical to Psi_m!):
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
END MODULE mod_blk_ice_an05
