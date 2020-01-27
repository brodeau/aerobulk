! AeroBulk / 2020 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_blk_ice
   !!====================================================================================
   !!       Computes turbulent components of surface fluxes over sea-ice
   !!
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (usually 2m) to zu (usually 10m) if needed
   !!   * the "effective" bulk wind speed at zu: U_blk (including gustiness contribution in unstable conditions)
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!       Routine turb_ice maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, January 2020
   !!
   !!====================================================================================
   USE mod_const       !: physical and othe constants
   USE mod_phymbl      !: thermodynamics

   IMPLICIT NONE
   PRIVATE

   !PUBLIC :: ICE_INIT, TURB_ICE
   PUBLIC :: TURB_ICE

   !!
   !REAL(wp), PARAMETER :: zi0   = 600._wp     ! scale height of the atmospheric boundary layer...
   !REAL(wp), PARAMETER :: Beta0 =  1.2_wp     ! gustiness parameter
   !!----------------------------------------------------------------------
CONTAINS

   !SUBROUTINE ice_init(l_use_cs, l_use_wl)
   !   !!---------------------------------------------------------------------
   !   !!                  ***  FUNCTION ice_init  ***
   !   !!
   !   !! INPUT :
   !   !! -------
   !   !!    * l_use_cs : use the cool-skin parameterization
   !   !!    * l_use_wl : use the warm-layer parameterization
   !   !!---------------------------------------------------------------------
   !   LOGICAL , INTENT(in) ::   l_use_cs ! use the cool-skin parameterization
   !   LOGICAL , INTENT(in) ::   l_use_wl ! use the warm-layer parameterization
   !   INTEGER :: ierr
   !   !!---------------------------------------------------------------------
   !   IF ( l_use_wl ) THEN
   !      ierr = 0
   !      ALLOCATE ( Tau_ac(jpi,jpj) , Qnt_ac(jpi,jpj), dT_wl(jpi,jpj), Hz_wl(jpi,jpj), STAT=ierr )
   !      IF( ierr > 0 ) CALL ctl_stop( ' ICE_INIT => allocation of Tau_ac, Qnt_ac, dT_wl & Hz_wl failed!' )
   !      Tau_ac(:,:) = 0._wp
   !      Qnt_ac(:,:) = 0._wp
   !      dT_wl(:,:)  = 0._wp
   !      Hz_wl(:,:)  = Hwl_max
   !   END IF
   !   IF ( l_use_cs ) THEN
   !      ierr = 0
   !      ALLOCATE ( dT_cs(jpi,jpj), STAT=ierr )
   !      IF( ierr > 0 ) CALL ctl_stop( ' ICE_INIT => allocation of dT_cs failed!' )
   !      dT_cs(:,:) = -0.25_wp  ! First guess of skin correction
   !   END IF
   !END SUBROUTINE ice_init



   SUBROUTINE turb_ice( kt, zt, zu, Ti_s, t_zt, qi_s, q_zt, U_zu,   &
      &                     Cd, Ch, Ce, t_zu, q_zu, U_blk,        &
      &                     xz0, xu_star, xL, xUN10 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_ice  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according
      !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!                Returns the effective bulk wind speed at zu to be used in the bulk formulas
      !!
      !! INPUT :
      !! -------
      !!    *  kt   : current time step (starts at 1)
      !!    *  zt   : height for temperature and spec. hum. of air            [m]
      !!    *  zu   : height for wind speed (usually 10m)                     [m]
      !!    *  t_zt : potential air temperature at zt                         [K]
      !!    *  q_zt : specific humidity of air at zt                          [kg/kg]
      !!    *  U_zu : scalar wind speed at zu                                 [m/s]
      !!
      !! INPUT/OUTPUT:
      !! -------------
      !!    *  Ti_s  : surface temperature of sea-ice                         [K]
      !!    *  qi_s  : SSQ aka saturation specific humidity at temp. Ti_s     [kg/kg]
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
      !!    * xL          : return the Monin-Obukhov length                    [m]
      !!    * xUN10       : neutral wind speed at 10m                          [m/s]
      !!
      !! ** Author: L. Brodeau, January 2020 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      INTEGER,  INTENT(in   )                     ::   kt       ! current time step
      REAL(wp), INTENT(in   )                     ::   zt       ! height for t_zt and q_zt                    [m]
      REAL(wp), INTENT(in   )                     ::   zu       ! height for U_zu                             [m]
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj) ::   Ti_s      ! sea surface temperature                [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   t_zt     ! potential air temperature              [Kelvin]
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj) ::   qi_s      ! sea surface specific humidity           [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   q_zt     ! specific air humidity at zt             [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   U_zu     ! relative wind module at zu                [m/s]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Cd       ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ch       ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ce       ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   t_zu     ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   q_zu     ! spec. humidity adjusted at zu           [kg/kg]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   U_blk    ! bulk wind speed at zu                     [m/s]
      !
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xz0  ! Aerodynamic roughness length   [m]
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xu_star  ! u*, friction velocity
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xL  ! zeta (zu/L)
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xUN10  ! Neutral wind at zu
      !
      INTEGER :: j_itt
      LOGICAL :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !
      REAL(wp), DIMENSION(:,:), ALLOCATABLE  ::  &
         &  u_star, t_star, q_star, &
         &  dt_zu, dq_zu,    &
         &  znu_a,           & !: Nu_air = kinematic viscosity of air
         &  z0, z0t
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zeta_u        ! stability parameter at height zu
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zeta_t        ! stability parameter at height zt
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: ztmp0, ztmp1, ztmp2
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zsst     ! to back up the initial bulk SST
      !
      LOGICAL :: lreturn_z0=.FALSE., lreturn_ustar=.FALSE., lreturn_L=.FALSE., lreturn_UN10=.FALSE.
      CHARACTER(len=40), PARAMETER :: crtnm = 'turb_ice@mod_blk_ice.f90'
      !!----------------------------------------------------------------------------------
      ALLOCATE ( u_star(jpi,jpj), t_star(jpi,jpj), q_star(jpi,jpj),  &
         &       zeta_u(jpi,jpj),  dt_zu(jpi,jpj),  dq_zu(jpi,jpj),  &
         &        znu_a(jpi,jpj),     z0(jpi,jpj),    z0t(jpi,jpj),  &
         &        ztmp0(jpi,jpj),  ztmp1(jpi,jpj),  ztmp2(jpi,jpj) )

      !IF ( kt == nit000 ) CALL ICE_INIT(l_use_cs, l_use_wl)

      IF( PRESENT(xz0) )     lreturn_z0    = .TRUE.
      IF( PRESENT(xu_star) ) lreturn_ustar = .TRUE.
      IF( PRESENT(xL) )      lreturn_L     = .TRUE.
      IF( PRESENT(xUN10) )   lreturn_UN10  = .TRUE.

      l_zt_equal_zu = .FALSE.
      IF( ABS(zu - zt) < 0.01_wp )   l_zt_equal_zu = .TRUE.    ! testing "zu == zt" is risky with double precision
      IF( .NOT. l_zt_equal_zu )  ALLOCATE( zeta_t(jpi,jpj) )

      !! First guess of temperature and humidity at height zu:
      t_zu = MAX( t_zt ,  180._wp )   ! who knows what's given on masked-continental regions...
      q_zu = MAX( q_zt , 1.e-6_wp )   !               "

      !! Pot. temp. difference (and we don't want it to be 0!)
      dt_zu = t_zu - Ti_s ;   dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6_wp), dt_zu )
      dq_zu = q_zu - qi_s ;   dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9_wp), dq_zu )

      znu_a = visc_air(t_zu) ! Air viscosity (m^2/s) at zt given from temperature in (K)

      U_blk = SQRT(U_zu*U_zu + 0.5_wp*0.5_wp) ! initial guess for wind gustiness contribution

      ztmp0   = LOG(    zu*10000._wp) ! optimization: 10000. == 1/z0 (with z0 first guess == 0.0001)
      ztmp1   = LOG(10._wp*10000._wp) !       "                    "               "
      u_star = 0.035_wp*U_blk*ztmp1/ztmp0       ! (u* = 0.035*Un10)

      !z0     = rough_leng(U_zu)*u_star*u_star/grav + 0.11_wp*znu_a/u_star
      z0     = MIN( MAX(ABS(z0), 1.E-9) , 1._wp )                      ! (prevents FPE from stupid values from masked region later on)

      z0t    = 1._wp / ( 0.1_wp*EXP(vkarmn/(0.00115/(vkarmn/ztmp1))) )
      z0t    = MIN( MAX(ABS(z0t), 1.E-9) , 1._wp )                      ! (prevents FPE from stupid values from masked region later on)

      Cd     = (vkarmn/ztmp0)**2    ! first guess of Cd

      ztmp0 = vkarmn*vkarmn/LOG(zt/z0t)/Cd

      ztmp2 = Ri_bulk( zu, Ti_s, t_zu, qi_s, q_zu, U_blk ) ! Bulk Richardson Number (BRN)

      !! First estimate of zeta_u, depending on the stability, ie sign of BRN (ztmp2):
      ztmp1 = 0.5 + SIGN( 0.5_wp , ztmp2 )
      ztmp0 = ztmp0*ztmp2
      !zeta_u = (1._wp-ztmp1) * (ztmp0/(1._wp+ztmp2/(-zu/(zi0*0.004_wp*Beta0**3)))) & !  BRN < 0
      !   &  +     ztmp1   * (ztmp0*(1._wp + 27._wp/9._wp*ztmp2/ztmp0))               !  BRN > 0
      !#LB: should make sure that the "ztmp0" of "27./9.*ztmp2/ztmp0" is "ztmp0*ztmp2" and not "ztmp0==vkarmn*vkarmn/LOG(zt/z0t)/Cd" !

      !! First guess M-O stability dependent scaling params.(u*,t*,q*) to estimate z0 and z/L
      ztmp0  = vkarmn/(LOG(zu/z0t) - psi_h_ice(zeta_u))

      u_star = MAX ( U_blk*vkarmn/(LOG(zu) - LOG(z0)  - psi_m_ice(zeta_u)) , 1.E-9 )  !  (MAX => prevents FPE from stupid values from masked region later on)
      t_star = dt_zu*ztmp0
      q_star = dq_zu*ztmp0

      ! What needs to be done if zt /= zu:
      IF( .NOT. l_zt_equal_zu ) THEN
         !! First update of values at zu (or zt for wind)
         zeta_t = zt*zeta_u/zu
         ztmp0 = psi_h_ice(zeta_u) - psi_h_ice(zeta_t)
         ztmp1 = LOG(zt/zu) + ztmp0
         t_zu = t_zt - t_star/vkarmn*ztmp1
         q_zu = q_zt - q_star/vkarmn*ztmp1
         q_zu = (0.5_wp + SIGN(0.5_wp,q_zu))*q_zu !Makes it impossible to have negative humidity :
         !
         dt_zu = t_zu - Ti_s  ; dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6_wp), dt_zu )
         dq_zu = q_zu - qi_s  ; dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9_wp), dq_zu )
      END IF

      !! ITERATION BLOCK
      DO j_itt = 1, nb_itt

         !!Inverse of Monin-Obukov length (1/L) :
         ztmp0 = One_on_L(t_zu, q_zu, u_star, t_star, q_star)  ! 1/L == 1/[Monin-Obukhov length]
         ztmp0 = SIGN( MIN(ABS(ztmp0),200._wp), ztmp0 ) ! (prevents FPE from stupid values from masked region later on...)

         ztmp1 = u_star*u_star   ! u*^2

         !! Update wind at zu with convection-related wind gustiness in unstable conditions (Fairall et al. 2003, Eq.8):
         !ztmp2 = Beta0*Beta0*ztmp1*(MAX(-zi0*ztmp0/vkarmn,0._wp))**(2._wp/3._wp) ! square of wind gustiness contribution, ztmp2 == Ug^2
         !!   ! Only true when unstable (L<0) => when ztmp0 < 0 => explains "-" before zi0
         U_blk = MAX(SQRT(U_zu*U_zu + ztmp2), 0.2_wp)        ! include gustiness in bulk wind speed
         ! => 0.2 prevents U_blk to be 0 in stable case when U_zu=0.

         !! Stability parameters:
         zeta_u = zu*ztmp0
         zeta_u = SIGN( MIN(ABS(zeta_u),50.0_wp), zeta_u )
         IF( .NOT. l_zt_equal_zu ) THEN
            zeta_t = zt*ztmp0
            zeta_t = SIGN( MIN(ABS(zeta_t),50.0_wp), zeta_t )
         END IF

         !! Adjustment the wind at 10m (not needed in the current algo form):
         !IF ( zu \= 10._wp ) U10 = U_zu + u_star/vkarmn*(LOG(10._wp/zu) - psi_m_ice(10._wp*ztmp0) + psi_m_ice(zeta_u))

         !! Roughness lengthes z0, z0t (z0q = z0t) :
         ztmp2 = u_star/vkarmn*LOG(10./z0)                                 ! Neutral wind speed at 10m
         !z0    = rough_leng(ztmp2)*ztmp1/grav + 0.11_wp*znu_a/u_star   ! Roughness length (eq.6) [ ztmp1==u*^2 ]
         z0     = MIN( MAX(ABS(z0), 1.E-9) , 1._wp )                      ! (prevents FPE from stupid values from masked region later on)

         ztmp1 = ( znu_a / (z0*u_star) )**0.72_wp
         z0t   = MIN( 1.6E-4_wp , 5.8E-5_wp*ztmp1 )
         z0t   = MIN( MAX(ABS(z0t), 1.E-9) , 1._wp )                      ! (prevents FPE from stupid values from masked region later on)

         !! Turbulent scales at zu :
         ztmp0   = psi_h_ice(zeta_u)
         ztmp1   = vkarmn/(LOG(zu) - LOG(z0t) - ztmp0) ! #LB: in ztmp0, some use psi_h_ice(zeta_t) rather than psi_h_ice(zeta_t) ???

         t_star = dt_zu*ztmp1
         q_star = dq_zu*ztmp1
         u_star = MAX( U_blk*vkarmn/(LOG(zu) - LOG(z0) - psi_m_ice(zeta_u)) , 1.E-9 )  !  (MAX => prevents FPE from stupid values from masked region later on)

         IF( .NOT. l_zt_equal_zu ) THEN
            !! Re-updating temperature and humidity at zu if zt /= zu :
            ztmp1 = LOG(zt/zu) + ztmp0 - psi_h_ice(zeta_t)
            t_zu = t_zt - t_star/vkarmn*ztmp1
            q_zu = q_zt - q_star/vkarmn*ztmp1
         END IF


         IF( .NOT. l_zt_equal_zu ) THEN
            dt_zu = t_zu - Ti_s ;  dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6_wp), dt_zu )
            dq_zu = q_zu - qi_s ;  dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9_wp), dq_zu )
         END IF

      END DO !DO j_itt = 1, nb_itt

      ! compute transfer coefficients at zu :
      ztmp0 = u_star/U_blk
      Cd   = ztmp0*ztmp0
      Ch   = ztmp0*t_star/dt_zu
      Ce   = ztmp0*q_star/dq_zu

      IF( lreturn_z0 )    xz0     = z0
      IF( lreturn_ustar ) xu_star = u_star
      IF( lreturn_L )     xL      = 1./One_on_L(t_zu, q_zu, u_star, t_star, q_star)
      IF( lreturn_UN10 )  xUN10   = u_star/vkarmn*LOG(10./z0)

      DEALLOCATE ( u_star, t_star, q_star, zeta_u, dt_zu, dq_zu, z0, z0t, znu_a, ztmp0, ztmp1, ztmp2 )
      IF( .NOT. l_zt_equal_zu ) DEALLOCATE( zeta_t )

   END SUBROUTINE turb_ice


   FUNCTION rough_leng_m( pus , pta )
      !!----------------------------------------------------------------------------------
      !! Computes the roughness length of sea-ice according to Andreas et al. 2005, (eq. 19)
      !!
      !! Author: L. Brodeau, January 2020 / AeroBulk  (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: rough_leng_m      ! roughness length of sea-ice [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pus ! u* = friction velocity    [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pta ! air temperature (only used to estimate kinematic viscosity of air) [K]
      !!
      INTEGER  :: ji, jj    ! dummy loop indices
      REAL(wp) :: znu, zus, zz
      !!----------------------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
            
            znu = visc_air( t_zu(ji,jj) ) ! Air viscosity (m^2/s) from temperature in (K)
            zus = pus(ji,jj)

            zz = (zus - 0.18_wp) / 0.1_wp
            
            rough_leng_m(ji,jj) = 0.135*znu/zus + 0.035*us*us/grav*( 5.*EXP(-zz*zz) + 1._wp )
            
         END DO
      END DO
      !!
   END FUNCTION rough_leng_m


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
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
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
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
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
END MODULE mod_blk_ice
