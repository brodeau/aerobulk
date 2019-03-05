! AeroBulk / 2016 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_blk_coare
   !!====================================================================================
   !!       Computes turbulent components of surface fluxes
   !!         according to Fairall et al. 2003 (COARE v3)
   !!
   !!       With Cool-Skin and Warm-Layer correction of SST (if needed)
   !!
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (2m) to zu (10m) if needed
   !!   * the effective bulk wind speed at 10m U_blk
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!    Using the bulk formulation/param. of COARE v3, Fairall et al. 2003
   !!      + consideration of cool-skin warm layer parametrization (Fairall et al. 1996)
   !!
   !!       Routine turb_coare maintained and developed in AeroBulk
   !!                     (http://aerobulk.sourceforge.net/)
   !!
   !!            Author: Laurent Brodeau, 2016
   !!
   !!====================================================================================
   USE mod_const   !: physical and othe constants
   USE mod_thermo  !: thermodynamics

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: TURB_COARE

   !! COARE own values for given constants:
   REAL(wp), PARAMETER ::  &
      &   zi0     = 600.,  & !: scale height of the atmospheric boundary layer...1
      &  Beta0    = 1.25,  & !: gustiness parameter
      &  charn0_max = 0.028  !: for COARE 3.5:
   !                         !:  -> VALUE above which the Charnock paramter levels off for winds > 18

CONTAINS

   SUBROUTINE turb_coare( cver, zt, zu, T_s, t_zt, q_s, q_zt, U_zu, &
      &                   Cd, Ch, Ce, t_zu, q_zu, U_blk,            &
      &                   rad_sw, rad_lw, slp,                      &
      &                   xz0, xu_star, xL, xUN10 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_coare  ***
      !!
      !!            2015: L. Brodeau
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to Fairall et al. (2003)
      !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!
      !!                Applies the cool-skin warm-layer correction of the SST to T_s
      !!                if the downwelling radiative fluxes at the surface (rad_sw & rad_lw)
      !!                and the SLP are provided as arguments!
      !!
      !! ** Method : Monin Obukhov Similarity Theory
      !!======================================================================================
      !!
      !! INPUT :
      !! -------
      !!    *  cver : version of COARE to use => '3.0' or '3.5'
      !!    *  zt   : height for temperature and spec. hum. of air            [m]
      !!    *  zu   : height for wind speed (generally 10m)                   [m]
      !!    *  U_zu : scalar wind speed at 10m                                [m/s]
      !!    *  t_zt : potential air temperature at zt                         [K]
      !!    *  q_zt : specific humidity of air at zt                          [kg/kg]
      !!
      !! INPUT/OUTPUT:
      !! -------------
      !!    *  T_s  : SST or skin temperature                                 [K]
      !!    *  q_s  : SSQ aka saturation specific humidity (at temp. T_s)     [kg/kg]
      !!              -> doesn't need to be given a value if skin temp computed (in case l_use_skin=True)
      !!              -> MUST be given the correct value if not computing skint temp. (in case l_use_skin=False)
      !!
      !! OPTIONAL INPUT (will trigger l_use_skin=TRUE if present!):
      !! ---------------
      !!    *  rad_sw : downwelling shortwave radiation at the surface (>0)   [W/m^2]
      !!    *  rad_lw : downwelling longwave radiation at the surface  (>0)   [W/m^2]
      !!    *  slp    : sea-level pressure                                    [Pa]
      !!
      !! OUTPUT :
      !! --------
      !!    *  Cd     : drag coefficient
      !!    *  Ch     : sensible heat coefficient
      !!    *  Ce     : evaporation coefficient
      !!    *  t_zu   : pot. air temperature adjusted at wind height zu       [K]
      !!    *  q_zu   : specific humidity of air        //                    [kg/kg]
      !!    *  U_blk  : bulk wind at 10m                                      [m/s]
      !!
      !! OPTIONAL OUTPUT:
      !! ----------------
      !!    * xz0         : return the aerodynamic roughness length (integration constant for wind stress) [m]
      !!    * xu_star     : return u* the friction velocity                    [m/s]
      !!    * xL          : return the Monin-Obukhov length                    [m]
      !!    * xUN10       : return the Monin-Obukhov length                    [m/s]
      !!
      !!============================================================================
      CHARACTER(len=3), INTENT(in   )             ::   cver     ! version of COARE to use (use '3.0' if u don't know what to chose!'
      REAL(wp), INTENT(in   )                     ::   zt       ! height for t_zt and q_zt                   [m]
      REAL(wp), INTENT(in   )                     ::   zu       ! height for U_zu                              [m]
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj) ::   T_s      ! sea surface temperature              [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   t_zt     ! potential air temperature            [Kelvin]
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj) ::   q_s      ! saturation sea surface spec. hum.     [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   q_zt     ! specific air humidity                 [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   U_zu     ! relative wind module at zu            [m/s]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Cd       ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ch       ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ce       ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   t_zu     ! pot. air temp. adjusted at zu             [K]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   q_zu     ! spec. humidity adjusted at zu             [kg/kg]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   U_blk    ! bulk wind at 10m                          [m/s]
      !
      REAL(wp), INTENT(in   ), OPTIONAL, DIMENSION(jpi,jpj) ::   rad_sw   !             [W/m^2]
      REAL(wp), INTENT(in   ), OPTIONAL, DIMENSION(jpi,jpj) ::   rad_lw   !             [W/m^2]
      REAL(wp), INTENT(in   ), OPTIONAL, DIMENSION(jpi,jpj) ::   slp      !             [Pa]
      !
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xz0  ! Aerodynamic roughness length   [m]
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xu_star  ! u*, friction velocity
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xL  ! zeta (zu/L)
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xUN10  ! Neutral wind at zu
      !
      INTEGER :: j_itt
      LOGICAL ::   l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !
      REAL(wp), DIMENSION(:,:), ALLOCATABLE  ::  &
         &  u_star, t_star, q_star, &
         &  dt_zu, dq_zu,    &
         &  zalpha,          & !: Charnock parameter
         &  znu_a,           & !: Nu_air, Viscosity of air
         &  z0, z0t
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   zeta_u        ! stability parameter at height zu
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   zeta_t        ! stability parameter at height zt
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   ztmp0, ztmp1, ztmp2
      !
      ! Cool skin:
      LOGICAL :: l_use_skin = .FALSE.
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: &
         &                zsst,   &  ! to back up the initial bulk SST
         &                zrhoa,  &  ! densitty of air
         &                zQsw,   &  ! net solar flux to the ocean (after albedo)
         &                zdelta     ! thickness of the viscous (skin) layer
      !
      LOGICAL :: lreturn_z0=.FALSE., lreturn_ustar=.FALSE., lreturn_L=.FALSE., lreturn_UN10=.FALSE.
      !!============================================================================

      ALLOCATE ( u_star(jpi,jpj), t_star(jpi,jpj), q_star(jpi,jpj), &
         &     zeta_u(jpi,jpj), zalpha(jpi,jpj),  &
         &     dt_zu(jpi,jpj), dq_zu(jpi,jpj),    &
         &     znu_a(jpi,jpj),   &
         &     z0(jpi,jpj), z0t(jpi,jpj),         &
         &     ztmp0(jpi,jpj), ztmp1(jpi,jpj), ztmp2(jpi,jpj) )

      ! Cool skin ?
      IF( PRESENT(rad_sw) .AND. PRESENT(rad_lw) .AND. PRESENT(slp) ) THEN
         l_use_skin = .TRUE.
         ALLOCATE ( zsst(jpi,jpj) , zrhoa(jpi,jpj), zQsw(jpi,jpj), zdelta(jpi,jpj) )
      END IF

      IF( PRESENT(xz0) )     lreturn_z0    = .TRUE.
      IF( PRESENT(xu_star) ) lreturn_ustar = .TRUE.
      IF( PRESENT(xL) )      lreturn_L     = .TRUE.
      IF( PRESENT(xUN10) )   lreturn_UN10  = .TRUE.


      l_zt_equal_zu = .FALSE.
      IF( ABS(zu - zt) < 0.01 ) l_zt_equal_zu = .TRUE.    ! testing "zu == zt" is risky with double precision

      IF( .NOT. l_zt_equal_zu )   ALLOCATE ( zeta_t(jpi,jpj) )

      !! Initialization for cool skin:
      IF( l_use_skin ) THEN
         zsst   = T_s    ! save the bulk SST
         zQsw   = (1. - oce_alb0)*rad_sw   ! Solar flux available for the ocean:
         zrhoa  = MAX(rho_air(t_zt, q_zt, slp), 1._wp) ! No updat needed! Fine enough!! For some reason seems to be negative sometimes
         T_s    = T_s - 0.25                      ! First guess of correction
         q_s    = 0.98*q_sat(MAX(T_s, 200._wp), slp) ! First guess of q_s
         zdelta = 0.001                    ! First guess of zdelta
      END IF

      !! First guess of temperature and humidity at height zu:
      t_zu = MAX(t_zt , 0.0_wp)    ! who knows what's given on masked-continental regions...
      q_zu = MAX(q_zt , 1.e-6_wp)  !               "

      !! Pot. temp. difference (and we don't want it to be 0!)
      dt_zu = t_zu - T_s ;  dt_zu = SIGN( MAX(ABS(dt_zu),1e-6_wp), dt_zu )
      dq_zu = q_zu - q_s ;  dq_zu = SIGN( MAX(ABS(dq_zu),1e-9_wp), dq_zu )

      znu_a = visc_air(t_zt) ! Air viscosity (m^2/s) at zt given from temperature in (K)

      ztmp2 = 0.5*0.5  ! initial guess for wind gustiness contribution
      U_blk = SQRT(U_zu*U_zu + ztmp2)

      ztmp2   = 10000.     ! optimization: ztmp2 == 1/z0 (with z0 first guess == 0.0001)
      ztmp0   = LOG(zu*ztmp2)
      ztmp1   = LOG(10.*ztmp2)
      u_star = 0.035*U_blk*ztmp1/ztmp0       ! (u* = 0.035*Un10)

      ! Charnock Parameter
      SELECT CASE (cver)
      CASE('3.0')
         zalpha = alfa_charn_3p0(U_zu)
      CASE('3.5')
         zalpha = MIN( 0.0017_wp*U_zu - 0.005_wp , charn0_max) !: alpha Charnock parameter (Eq. 13 Edson al. 2013)
         zalpha = MAX( zalpha , 0._wp )
      CASE DEFAULT
         PRINT *, 'Unknown version for COARE algorithm: ',cver ; PRINT *, ''
         STOP
      END SELECT

      z0     = zalpha*u_star*u_star/grav + 0.11*znu_a/u_star
      z0t    = 1. / ( 0.1*EXP(vkarmn/(0.00115/(vkarmn/ztmp1))) )

      ztmp2  = vkarmn/ztmp0
      Cd     = ztmp2*ztmp2    ! first guess of Cd

      ztmp0 = vkarmn*vkarmn/LOG(zt/z0t)/Cd

      !Ribcu = -zu/(zi0*0.004*Beta0**3) !! Saturation Rib, zi0 = tropicalbound. layer depth
      ztmp2  = grav*zu*(dt_zu + rctv0*t_zu*dq_zu)/(t_zu*U_blk*U_blk)  !! Ribu Bulk Richardson number
      ztmp1 = 0.5 + SIGN(0.5_wp , ztmp2)
      ztmp0 = ztmp0*ztmp2
      !!             Ribu < 0                                 Ribu > 0   Beta = 1.25
      zeta_u = (1.-ztmp1) * (ztmp0/(1.+ztmp2/(-zu/(zi0*0.004*Beta0**3)))) &
         &  +     ztmp1   * (ztmp0*(1. + 27./9.*ztmp2/ztmp0))

      !! First guess M-O stability dependent scaling params.(u*,t*,q*) to estimate z0 and z/L
      ztmp0  = vkarmn/(LOG(zu/z0t) - psi_h_coare(zeta_u))

      u_star = U_blk*vkarmn/(LOG(zu) - LOG(z0)  - psi_m_coare(zeta_u))
      t_star = dt_zu*ztmp0
      q_star = dq_zu*ztmp0

      ! What's need to be done if zt /= zu:
      IF( .NOT. l_zt_equal_zu ) THEN

         zeta_t = zt*zeta_u/zu

         !! First update of values at zu (or zt for wind)
         ztmp0 = psi_h_coare(zeta_u) - psi_h_coare(zeta_t)
         ztmp1 = LOG(zt/zu) + ztmp0
         t_zu = t_zt - t_star/vkarmn*ztmp1
         q_zu = q_zt - q_star/vkarmn*ztmp1
         q_zu = (0.5 + SIGN(0.5_wp,q_zu))*q_zu !Makes it impossible to have negative humidity :

         dt_zu = t_zu - T_s  ; dt_zu = SIGN( MAX(ABS(dt_zu),1e-6_wp), dt_zu )
         dq_zu = q_zu - q_s  ; dq_zu = SIGN( MAX(ABS(dq_zu),1e-9_wp), dq_zu )

      END IF

      !! ITERATION BLOCK
      DO j_itt = 1, nb_itt

         !!Inverse of Monin-Obukov length (1/L) :
         ztmp0 = One_on_L(t_zu, q_zu, u_star, t_star, q_star)  ! 1/L == 1/[Monin-Obukhov length]

         ztmp1 = u_star*u_star   ! u*^2

         !! Update wind at 10m taking into acount convection-related wind gustiness:
         ! Ug = Beta*w*  (Beta = 1.25, Fairall et al. 2003, Eq.8):
         ztmp2 = Beta0*Beta0*ztmp1*(MAX(-zi0*ztmp0/vkarmn,0._wp))**(2./3.)   ! => ztmp2 == Ug^2
         !!   ! Only true when unstable (L<0) => when ztmp0 < 0 => explains "-" before 600.
         U_blk = MAX(sqrt(U_zu*U_zu + ztmp2), 0.2_wp)        ! include gustiness in bulk wind speed
         ! => 0.2 prevents U_blk to be 0 in stable case when U_zu=0.

         IF( cver == '3.5' ) THEN
            !! Need to update Charnock parameter from neutral wind speed!
            ztmp2 = u_star/vkarmn*LOG(10./z0)   ! UN10 Neutral wind at 10m!
            zalpha = MIN( 0.0017_wp*ztmp2 - 0.005_wp , charn0_max)  ! alpha Charnock parameter (Eq. 13 Edson al. 2013)
            zalpha = MAX( zalpha , 0._wp )
         END IF

         !! Roughness lengthes z0, z0t (z0q = z0t) :
         z0    = zalpha*ztmp1/grav + 0.11*znu_a/u_star  ! Roughness length (eq.6)
         ztmp1 = z0*u_star/znu_a                        ! Re_r: roughness Reynolds number

         SELECT CASE (cver)
         CASE('3.0')
            z0t   = MIN( 1.1E-4_wp , 5.5E-5_wp*ztmp1**(-0.6_wp) ) ! Scalar roughness for both theta and q (eq.28)
         CASE('3.5')
            ! Chris Fairall and Jim Edsson, private communication, March 2016 / COARE 3.5 :
            !  -> these thermal roughness lengths give CE and CH that closely approximate COARE3.0
            z0t   = MIN( 1.6e-4_wp , 5.8E-5_wp*ztmp1**(-0.72_wp))
            !                                            !
         END SELECT

         !! Stability parameters:
         zeta_u = zu*ztmp0 ; zeta_u = sign( min(abs(zeta_u),50.0_wp), zeta_u )
         IF( .NOT. l_zt_equal_zu ) THEN
            zeta_t = zt*ztmp0 ;  zeta_t = sign( min(abs(zeta_t),50.0_wp), zeta_t )
         END IF

         !! Turbulent scales at zu=10m :
         ztmp0   = psi_h_coare(zeta_u)
         ztmp1   = vkarmn/(LOG(zu) - LOG(z0t) - ztmp0)

         t_star = dt_zu*ztmp1
         q_star = dq_zu*ztmp1
         u_star = U_blk*vkarmn/(LOG(zu) - LOG(z0) - psi_m_coare(zeta_u))

         IF( .NOT. l_zt_equal_zu ) THEN
            ! What's need to be done if zt /= zu
            !! Re-updating temperature and humidity at zu :
            ztmp2 = ztmp0 - psi_h_coare(zeta_t)
            ztmp1 = log(zt/zu) + ztmp2
            t_zu = t_zt - t_star/vkarmn*ztmp1
            q_zu = q_zt - q_star/vkarmn*ztmp1
         END IF

         !! SKIN related part
         IF( l_use_skin ) THEN
            CALL CSWL_COARE( t_zu, q_zu, zsst, slp, U_blk, u_star, t_star, q_star, &
               &             zrhoa, rad_lw, zQsw, zdelta, T_s, q_s )
         END IF

         dt_zu = t_zu - T_s ;  dt_zu = SIGN( MAX(ABS(dt_zu),1e-6_wp), dt_zu )
         dq_zu = q_zu - q_s ;  dq_zu = SIGN( MAX(ABS(dq_zu),1e-9_wp), dq_zu )

      END DO
      !
      ! compute transfer coefficients at zu :
      ztmp0 = u_star/U_blk
      Cd   = ztmp0*ztmp0
      Ch   = ztmp0*t_star/dt_zu
      Ce   = ztmp0*q_star/dq_zu

      IF( lreturn_z0 )    xz0     = z0
      IF( lreturn_ustar ) xu_star = u_star
      IF( lreturn_L )     xL      = 1./One_on_L(t_zu, q_zu, u_star, t_star, q_star)
      IF( lreturn_UN10 )  xUN10   = u_star/vkarmn*LOG(10./z0)

      DEALLOCATE ( u_star, t_star, q_star, zeta_u, zalpha, dt_zu, dq_zu, z0, z0t, znu_a, ztmp0, ztmp1, ztmp2 )
      IF( .NOT. l_zt_equal_zu )   DEALLOCATE ( zeta_t )

      IF( l_use_skin ) THEN
         DEALLOCATE ( zsst, zrhoa, zQsw, zdelta )
      END IF

   END SUBROUTINE turb_coare


   FUNCTION alfa_charn_3p0(pwnd)
      !!-------------------------------------------------------------------
      !! Compute the Charnock parameter as a function of the wind speed
      !!
      !! (Fairall et al., 2003 p.577-578)
      !!
      !! Wind below 10 m/s :  alfa = 0.011
      !! Wind between 10 and 18 m/s : linear increase from 0.011 to 0.018
      !! Wind greater than 18 m/s :  alfa = 0.018
      !!
      !! Author: L. Brodeau, june 2016 / AeroBulk  (https://sourceforge.net/p/aerobulk)
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: alfa_charn_3p0
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pwnd   ! wind speed
      !
      INTEGER  ::   ji, jj         ! dummy loop indices
      REAL(wp) :: zw, zgt10, zgt18
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zw = pwnd(ji,jj)   ! wind speed
            !
            ! Charnock's constant, increases with the wind :
            zgt10 = 0.5 + SIGN(0.5_wp,(zw - 10)) ! If zw<10. --> 0, else --> 1
            zgt18 = 0.5 + SIGN(0.5_wp,(zw - 18.)) ! If zw<18. --> 0, else --> 1
            !
            alfa_charn_3p0(ji,jj) =  (1. - zgt10)*0.011    &    ! wind is lower than 10 m/s
               &     + zgt10*((1. - zgt18)*(0.011 + (0.018 - 0.011) &
               &      *(zw - 10.)/(18. - 10.)) + zgt18*( 0.018 ) )    ! Hare et al. (1999)
            !
         END DO
      END DO
      !
   END FUNCTION alfa_charn_3p0


   FUNCTION One_on_L( ptha, pqa, pus, pts, pqs )
      !!------------------------------------------------------------------------
      !!
      !! Evaluates the 1./(Monin Obukhov length) from air temperature and
      !!  specific humidity, and frictional scales u*, t* and q*
      !!
      !! Author: L. Brodeau, june 2016 / AeroBulk
      !!         (https://sourceforge.net/p/aerobulk)
      !!------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: One_on_L         !: 1./(Monin Obukhov length) [m^-1]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: ptha,  &  !: average potetntial air temperature [K]
         &                                        pqa,   &  !: average specific humidity of air   [kg/kg]
         &                                      pus, pts, pqs   !: frictional velocity, temperature and humidity
      !
      INTEGER  ::   ji, jj         ! dummy loop indices
      REAL(wp) ::     zqa          ! local scalar
      !!-------------------------------------------------------------------
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zqa = (1. + rctv0*pqa(ji,jj))
            !
            One_on_L(ji,jj) =  grav*vkarmn*(pts(ji,jj)*zqa + rctv0*ptha(ji,jj)*pqs(ji,jj)) &
               &                      / ( pus(ji,jj)*pus(ji,jj) * ptha(ji,jj)*zqa )
            !
         END DO
      END DO
      !
   END FUNCTION One_on_L


   FUNCTION psi_m_coare( pzeta )
      !!----------------------------------------------------------------------------------
      !! ** Purpose: compute the universal profile stability function for momentum
      !!             COARE 3.0, Fairall et al. 2003
      !!             pzeta : stability paramenter, z/L where z is altitude
      !!                     measurement and L is M-O length
      !!       Stability function for wind speed and scalars matching Kansas and free
      !!       convection forms with weighting f convective form, follows Fairall et
      !!       al (1996) with profile constants from Grachev et al (2000) BLM stable
      !!       form from Beljaars and Holtslag (1991)
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: psi_m_coare
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) :: zta, zphi_m, zphi_c, zpsi_k, zpsi_c, zf, zc, zstab
      !!----------------------------------------------------------------------------------
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zta = pzeta(ji,jj)
            !
            zphi_m = ABS(1. - 15.*zta)**.25    !!Kansas unstable
            !
            zpsi_k = 2.*LOG((1. + zphi_m)/2.) + LOG((1. + zphi_m*zphi_m)/2.)   &
               & - 2.*ATAN(zphi_m) + 0.5*rpi
            !
            zphi_c = ABS(1. - 10.15*zta)**.3333                   !!Convective
            !
            zpsi_c = 1.5*LOG((1. + zphi_c + zphi_c*zphi_c)/3.) &
               &     - 1.7320508*ATAN((1. + 2.*zphi_c)/1.7320508) + 1.813799447
            !
            zf = zta*zta
            zf = zf/(1. + zf)
            zc = MIN(50._wp, 0.35_wp*zta)
            zstab = 0.5 + SIGN(0.5_wp, zta)
            !
            psi_m_coare(ji,jj) = (1. - zstab) * ( (1. - zf)*zpsi_k + zf*zpsi_c ) & ! (zta < 0)
               &                -   zstab     * ( 1. + 1.*zta     &                ! (zta > 0)
               &                         + 0.6667*(zta - 14.28)/EXP(zc) + 8.525 )   !     "
            !
         END DO
      END DO
      !
   END FUNCTION psi_m_coare


   FUNCTION psi_h_coare( pzeta )
      !!---------------------------------------------------------------------
      !! Universal profile stability function for temperature and humidity
      !! COARE 3.0, Fairall et al. 2003
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! Stability function for wind speed and scalars matching Kansas and free
      !! convection forms with weighting f convective form, follows Fairall et
      !! al (1996) with profile constants from Grachev et al (2000) BLM stable
      !! form from Beljaars and Holtslag (1991)
      !!
      !! Author: L. Brodeau, june 2016 / AeroBulk
      !!         (https://sourceforge.net/p/aerobulk)
      !!----------------------------------------------------------------
      !!
      REAL(wp), DIMENSION(jpi,jpj) :: psi_h_coare
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) :: zta, zphi_h, zphi_c, zpsi_k, zpsi_c, zf, zc, zstab
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zta = pzeta(ji,jj)
            !
            zphi_h = (ABS(1. - 15.*zta))**.5  !! Kansas unstable   (zphi_h = zphi_m**2 when unstable, zphi_m when stable)
            !
            zpsi_k = 2.*LOG((1. + zphi_h)/2.)
            !
            zphi_c = (ABS(1. - 34.15*zta))**.3333   !! Convective
            !
            zpsi_c = 1.5*LOG((1. + zphi_c + zphi_c*zphi_c)/3.) &
               &    -1.7320508*ATAN((1. + 2.*zphi_c)/1.7320508) + 1.813799447
            !
            zf = zta*zta
            zf = zf/(1. + zf)
            zc = MIN(50._wp,0.35_wp*zta)
            zstab = 0.5 + SIGN(0.5_wp, zta)
            !
            psi_h_coare(ji,jj) = (1. - zstab) * ( (1. - zf)*zpsi_k + zf*zpsi_c ) &
               &                -   zstab     * ( (ABS(1. + 2.*zta/3.))**1.5     &
               &                           + .6667*(zta - 14.28)/EXP(zc) + 8.525 )
            !
         END DO
      END DO
      !
   END FUNCTION psi_h_coare


   SUBROUTINE CSWL_COARE( pTzu, pqzu, psst, pslp, pU10, pus, pts, pqs, &
      &                   prhoa, pRlw, pQsw, pdelta, pT_s, pq_s )
      !!---------------------------------------------------------------------
      !!
      !!  Cool-Skin Warm-Layer scheme according to Fairall et al. 1996
      !!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pTzu, pqzu, psst, pslp, &
         &                                           pU10, pus, pts, pqs
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: prhoa, pRlw, pQsw
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pdelta, pT_s, pq_s
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) :: zz0, zz1, zz2, zCe, zCh, zus, zQsen, zQlat, zQlw, zfr, &
         &        zdt, zdq, ztf, &
         &        zdelta, zlamb, zalpha, zQt
      !!---------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi

            zdelta = pdelta(ji,jj)

            zdt = pTzu(ji,jj) - pT_s(ji,jj)  ; zdt = SIGN( MAX(ABS(zdt),1e-6_wp), zdt )
            zdq = pqzu(ji,jj) - pq_s(ji,jj)  ; zdq = SIGN( MAX(ABS(zdq),1e-9_wp), zdq )

            !! compute transfer coefficients at zu :
            zz0 = pus(ji,jj)/pU10(ji,jj)    !LB! not needed inside the loop !
            zCh = zz0*pts(ji,jj)/zdt
            zCe = zz0*pqs(ji,jj)/zdq

            ! Turbulent heat fluxes:
            zz1 = prhoa(ji,jj)*pU10(ji,jj)
            zQlat = MIN( L0vap*zCe*zz1*(pqzu(ji,jj) - pq_s(ji,jj)) , 0._wp )
            zQsen =      Cp0_a*zCh*zz1*(pTzu(ji,jj) - pT_s(ji,jj))

            ! Net longwave flux:
            zz1  = pT_s(ji,jj)*pT_s(ji,jj)
            zQlw = 0.97*(pRlw(ji,jj) - sigma0*zz1*zz1)

            !! Fraction of the shortwave flux absorbed by the cool-skin sublayer:
            !zQsw_f = 0.065 + 11.*zdelta - 6.6e-5/zdelta*(1. - EXP(-zdelta/8.e-4)) ! Eq.16 (Fairall al. 1996b)
            zfr = 0.137 + 11.*zdelta - 6.6e-5/zdelta*(1. - EXP(-zdelta/8.e-4)) ! Eq.16 (Fairall al. 1996b)
            !LB: why 0.065 and not 0.137 like in the paper??? Beljaars & Zeng use 0.065
            !LB: maybe comes from Wick et al 2005 ...

            zQt = -(zQlw + zQsen + zQlat + zfr*pQsw(ji,jj))  ! Total cooling at the interface

            ztf    = 0.5 + SIGN(0.5_wp, zQt) ! Qt > 0 => cooling of the layer => ztf = 1
            !                               Qt < 0 => warming of the layer => ztf = 0

            zalpha = 2.1e-5*MAX(pT_s(ji,jj)-rt0 + 3.2_wp, 0._wp)**0.79  ! alpha = thermal expansion of water (~2.5E-4) LB: remove from loop, sst accurate enough!

            !! Term alpha*Qb (Qb is the virtual surface cooling inc. buoyancy effect of salinity due to evap):
            zz1 = zalpha*zQt - 0.026*zQlat*Cp0_w/L0vap  ! alpha*(Eq.8) == alpha*Qb "-" because Qlat < 0
            !! LB: this terms only makes sense if > 0 i.e. in the cooling case
            !! so similar to what's donce in ECMWF:
            zz1 = MAX(0._wp , zz1)    ! 1. instead of 0.1 though ZQ = MAX(1.0,-pQlw(ji,jj) - pQsen(ji,jj) - pQlat(ji,jj))

            !! Laurent: too low wind (u*) might cause problem in stable cases:
            zus = MAX(pus(ji,jj), 1.E-4_wp)

            ! Lambda (=> zz0, empirical coeff.) (Eq.14):
            zz0 = 16. * zz1 * grav * rho0_w * Cp0_w * nu0_w*nu0_w*nu0_w  ! (numerateur) zz1 == alpha*Q
            zz2 = zus*zus * prhoa(ji,jj) / rho0_w * k0_w
            zz2 =  zz2*zz2                                             ! denominateur
            !LB:  zz0 has the sign of zz1 and therefore of Qb !
            zlamb =  6.*( 1. + (zz0/zz2)**(3./4.) )**(-1./3.) !  Eq.14   (Saunders)

            ! Updating molecular sublayer thickness (delta):
            zz2    = nu0_w/(SQRT(prhoa(ji,jj)/rho0_w)*zus)
            zdelta =      ztf    *          zlamb*zz2   &  ! Eq.12 (when alpha*Qb>0 / cooling of layer)
               &    + (1. - ztf) * MIN(0.007_wp , 6._wp*zz2 )    ! Eq.12 (when alpha*Qb<0 / warming of layer)
            !LB: changed 0.01 to 0.007
            pdelta(ji,jj) = zdelta

            ! Updating pT_s and q_s ...
            zz2 =  - zQt*zdelta/k0_w   ! temperature increment !  Eq.13 Cool skin
            !
            pT_s(ji,jj) = psst(ji,jj) + zz2
            !
         END DO
      END DO

      pq_s = 0.98*q_sat(MAX(pT_s, 200._wp), pslp)   !skin !LB: just to avoid problem on masked regions

   END SUBROUTINE CSWL_COARE


END MODULE mod_blk_coare
