! AeroBulk / 2016 / L. Brodeau (brodeau@gmail.com)
! https://sourceforge.net/p/aerobulk
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_blk_ecmwf
   !!====================================================================================
   !!       Computes turbulent components of surface fluxes
   !!         according to IFS of the ECMWF
   !!
   !!       With Cool-Skin and Warm-Layer correction of SST (if needed)
   !!
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (2m) to zu (10m) if needed
   !!   * the effective bulk wind speed at 10m U_blk
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!    Using the bulk formulation/param. of ECMWF
   !!      + consideration of cool-skin warm layer parametrization (CS: Fairall et al. 1996; WL: Zeng & Beljaars, 2005 )
   !!
   !!       Routine turb_ecmwf maintained and developed in AeroBulk
   !!                     (http://aerobulk.sourceforge.net/)
   !!
   !!            Author: Laurent Brodeau, 2016, brodeau@gmail.com
   !!
   !!====================================================================================
   USE mod_const   !: physical and othe constants
   USE mod_thermo  !: thermodynamics

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: TURB_ECMWF

   !! ECMWF own values for given constants, taken form IFS documentation...
   REAL(wp),    PARAMETER :: &
      & charn0 = 0.018, &      !:  Charnock constant (pretty high value here!!!
                                !!                       !:   =>  Usually 0.011 for moderate winds)
                                !!  Note that in the ECMWF system, when coupled to tje wave model, the
                                !!  the Charnock parameter is provided by the wave model!
      &   zi0     = 1000.,  &  !: scale height of the atmospheric boundary layer...1
      &  Beta0    = 1. ,    &  !: gustiness parameter ( = 1.25 in COAREv3)
      &   alpha_M = 0.11,   &  !: For roughness length (smooth surface term)
      &   alpha_H = 0.40,   &  !: (Chapter 3, p.34, IFS doc Cy31r1)
      &   alpha_Q = 0.62

   !! Cool-Skin / Warm-Layer related parameters:
   REAL(wp),    PARAMETER :: &
      &  rdt0    = 3600.*1.5, & !: time step
      &  rd0     = 3. ,      &  !: Depth scale [m], "d" in Eq.11 (Zeng & Beljaars 2005)
      &  rNu0    = 0.5          !: Nu (exponent of temperature profile) Eq.11
   !                            !: (Zeng & Beljaars 2005) !: set to 0.5 instead of
   !                            !: 0.3 to respect a warming of +3 K in calm
   !                            !: condition for the insolation peak of +1000W/m^2
   INTEGER,    PARAMETER :: &
      &  nb_itt_wl = 10         !: number of sub-itterations for solving the differential equation in warm-layer part
   !                            !:  => use "nb_itt_wl = 1" for No itteration! => way cheaper !!!
   !                            !:    => assumes balance between the last 2 terms of Eq.11 (lhs of eq.11 = 0)
   !                            !:    => in that case no need for sub-itterations !
   !                            !:    => ACCEPTABLE IN MOST CONDITIONS ! (UNLESS: sunny + very calm/low-wind conditions)
   !                            !:  => Otherwize use "nb_itt_wl = 10"

CONTAINS

   SUBROUTINE TURB_ECMWF( zt, zu, T_s, t_zt, q_s, q_zt, U_zu,       &
      &                   Cd, Ch, Ce, t_zu, q_zu, U_blk,            &
      &                   rad_sw, rad_lw, slp,                      &
      &                   xz0, xu_star, xL, xUN10 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_ecmwf  ***
      !!
      !!            2015: L. Brodeau (brodeau@gmail.com)
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to IFS doc. (cycle 40)
      !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!
      !!                Applies the cool-skin warm-layer correction of the SST to T_s
      !!                if the downwelling radiative fluxes at the surface (rad_sw & rad_lw)
      !!                and the SLP are provided as arguments!
      !!
      !! ** Method : Monin Obukhov Similarity Theory
      !!
      !! INPUT :
      !! -------
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
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
      !!----------------------------------------------------------------------------------
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
      LOGICAL :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !
      REAL(wp), DIMENSION(:,:), ALLOCATABLE  ::  &
         &  u_star, t_star, q_star, &
         &  dt_zu, dq_zu,    &
         &  znu_a,           & !: Nu_air, Viscosity of air
         &  Linv,            & !: 1/L (inverse of Monin Obukhov length...
         &  z0, z0t, z0q
      !
      ! Cool skin:
      LOGICAL :: l_use_skin = .FALSE.
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: &
         &                zsst,   &  ! to back up the initial bulk SST
         &                zrhoa,  &  ! densitty of air
         &                zQsw       ! thickness of the viscous (skin) layer
      !
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   func_m, func_h
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   ztmp0, ztmp1, ztmp2
      !
      LOGICAL :: lreturn_z0=.FALSE., lreturn_ustar=.FALSE., lreturn_L=.FALSE., lreturn_UN10=.FALSE.
      !!----------------------------------------------------------------------------------

      ALLOCATE ( u_star(jpi,jpj), t_star(jpi,jpj), q_star(jpi,jpj), &
         &     func_m(jpi,jpj), func_h(jpi,jpj),  &
         &     dt_zu(jpi,jpj), dq_zu(jpi,jpj),    &
         &     znu_a(jpi,jpj), Linv(jpi,jpj),   &
         &     z0(jpi,jpj), z0t(jpi,jpj), z0q(jpi,jpj), &
         &     ztmp0(jpi,jpj), ztmp1(jpi,jpj), ztmp2(jpi,jpj) )

      ! Cool skin ?
      IF( PRESENT(rad_sw) .AND. PRESENT(rad_lw) .AND. PRESENT(slp) ) THEN
         l_use_skin = .TRUE.
         ALLOCATE ( zsst(jpi,jpj) , zrhoa(jpi,jpj), zQsw(jpi,jpj) )
      END IF

      IF( PRESENT(xz0) )     lreturn_z0    = .TRUE.
      IF( PRESENT(xu_star) ) lreturn_ustar = .TRUE.
      IF( PRESENT(xL) )      lreturn_L     = .TRUE.
      IF( PRESENT(xUN10) )   lreturn_UN10  = .TRUE.
      !
      ! Identical first gess as in COARE, with IFS parameter values though
      !
      l_zt_equal_zu = .FALSE.
      IF( ABS(zu - zt) < 0.01 )   l_zt_equal_zu = .TRUE.    ! testing "zu == zt" is risky with double precision

      !! Initialization for cool skin:
      IF( l_use_skin ) THEN
         zsst   = T_s    ! save the bulk SST
         zQsw   = (1. - oce_alb0)*rad_sw   ! Solar flux available for the ocean:
         zrhoa  = MAX(rho_air(t_zt, q_zt, slp), 1.) ! No updat needed! Fine enough!! For some reason seems to be negative sometimes
         T_s    = T_s - 0.25                      ! First guess of correction
         q_s    = 0.98*q_sat(MAX(T_s, 200.), slp) ! First guess of q_s
      END IF

      !! First guess of temperature and humidity at height zu:
      t_zu = MAX( t_zt , 0.0  )   ! who knows what's given on masked-continental regions...
      q_zu = MAX( q_zt , 1.e-6)   !               "

      !! Pot. temp. difference (and we don't want it to be 0!)
      dt_zu = t_zu - T_s ;   dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6), dt_zu )
      dq_zu = q_zu - q_s ;   dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9), dq_zu )

      znu_a = visc_air(t_zt) ! Air viscosity (m^2/s) at zt given from temperature in (K)

      ztmp2 = 0.5*0.5  ! initial guess for wind gustiness contribution
      U_blk = SQRT(U_zu*U_zu + ztmp2)

      ztmp2   = 10000.     ! optimization: ztmp2 == 1/z0 (with z0 first guess == 0.0001)
      ztmp0   = LOG(zu*ztmp2)
      ztmp1   = LOG(10.*ztmp2)
      u_star = 0.035*U_blk*ztmp1/ztmp0       ! (u* = 0.035*Un10)

      z0     = charn0*u_star*u_star/grav + 0.11*znu_a/u_star
      z0t    = 1. / ( 0.1*EXP(vkarmn/(0.00115/(vkarmn/ztmp1))) )

      ztmp2  = vkarmn/ztmp0
      Cd     = ztmp2*ztmp2    ! first guess of Cd

      ztmp0 = vkarmn*vkarmn/LOG(zt/z0t)/Cd

      ztmp2 = Ri_bulk( zu, t_zu, dt_zu, q_zu, dq_zu, U_blk )   ! Ribu = Bulk Richardson number

      !! First estimate of zeta_u, depending on the stability, ie sign of Ribu (ztmp2):
      ztmp1 = 0.5 + SIGN( 0.5 , ztmp2 )
      ztmp0 = ztmp0*ztmp2
      !!             Ribu < 0                                 Ribu > 0   Beta = 1.25
      func_h = (1.-ztmp1) * (ztmp0/(1.+ztmp2/(-zu/(zi0*0.004*Beta0**3)))) &  ! temporary array !!! func_h == zeta_u
         &  +     ztmp1   * (ztmp0*(1. + 27./9.*ztmp2/ztmp0))

      !! First guess M-O stability dependent scaling params.(u*,t*,q*) to estimate z0 and z/L
      ztmp0  = vkarmn/(LOG(zu/z0t) - psi_h_ecmwf(func_h))

      u_star = U_blk*vkarmn/(LOG(zu) - LOG(z0)  - psi_m_ecmwf(func_h))
      t_star = dt_zu*ztmp0
      q_star = dq_zu*ztmp0

      ! What's need to be done if zt /= zu:
      IF( .NOT. l_zt_equal_zu ) THEN
         !
         !! First update of values at zu (or zt for wind)
         ztmp0 = psi_h_ecmwf(func_h) - psi_h_ecmwf(zt*func_h/zu)    ! zt*func_h/zu == zeta_t
         ztmp1 = LOG(zt/zu) + ztmp0
         t_zu = t_zt - t_star/vkarmn*ztmp1
         q_zu = q_zt - q_star/vkarmn*ztmp1
         q_zu = (0.5 + SIGN(0.5,q_zu))*q_zu !Makes it impossible to have negative humidity :

         dt_zu = t_zu - T_s  ; dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6), dt_zu )
         dq_zu = q_zu - q_s  ; dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9), dq_zu )
         !
      ENDIF

      !! First guess of inverse of Monin-Obukov length (1/L) :
      ztmp0 = (1. + rctv0*q_zu)  ! the factor to apply to temp. to get virt. temp...
      Linv  =  grav*vkarmn*(t_star*ztmp0 + rctv0*t_zu*q_star) / ( u_star*u_star * t_zu*ztmp0 )

      !! Functions such as  u* = U_blk*vkarmn/func_m
      ztmp0 = zu*Linv
      func_m = LOG(zu) - LOG(z0)  - psi_m_ecmwf(ztmp0) + psi_m_ecmwf( z0*Linv)
      func_h = LOG(zu) - LOG(z0t) - psi_h_ecmwf(ztmp0) + psi_h_ecmwf(z0t*Linv)

      !! ITERATION BLOCK
      DO j_itt = 1, nb_itt

         !! Bulk Richardson Number at z=zu (Eq. 3.25)
         ztmp0 = Ri_bulk(zu, t_zu, dt_zu, q_zu, dq_zu, U_blk)

         !! New estimate of the inverse of the Monin-Obukhon length (Linv == zeta/zu) :
         Linv = ztmp0*func_m*func_m/func_h / zu ! From Eq. 3.23, Chap.3.2.3, IFS doc - Cy40r1
         !! Note: it is slightly different that the L we would get with the usual
         !! expression, as in coare algorithm or in 'mod_thermo.f90' (One_on_L_MO())

         !! Update func_m with new Linv:
         func_m = LOG(zu) -LOG(z0) - psi_m_ecmwf(zu*Linv) + psi_m_ecmwf(z0*Linv)

         !! Need to update roughness lengthes:
         u_star = U_blk*vkarmn/func_m
         ztmp1  = u_star*u_star
         ztmp2  = znu_a/u_star
         z0     = alpha_M*ztmp2 + charn0*ztmp1/grav
         z0t    = alpha_H*ztmp2                              ! eq.3.26, Chap.3, p.34, IFS doc - Cy31r1
         z0q    = alpha_Q*ztmp2

         !! Update wind at 10m taking into acount convection-related wind gustiness:
         !! => Chap. 3.2, IFS doc - Cy40r1, Eq.3.17 and Eq.3.18 + Eq.3.8
         ztmp1 = ztmp1 * MAX( -zi0*Linv/vkarmn ,0. )**(2./3.) ! => w*^2
         !! => equivalent using Beta=1 (gustiness parameter, 1.25 for COARE, also zi0=600 in COARE..)
         U_blk = MAX(SQRT(U_zu*U_zu + ztmp1), 0.2)    ! => 0.2 prevents U_blk to be 0 in stable case when U_zu=0.

         !! Need to update "theta" and "q" at zu in case they are given at different heights
         !! as well the air-sea differences:
         IF( .NOT. l_zt_equal_zu ) THEN

            !! Arrays func_m and func_h are free for a while so using them as temporary arrays...
            func_h = psi_h_ecmwf(zu*Linv) ! temporary array !!!
            func_m = psi_h_ecmwf(zt*Linv) ! temporary array !!!

            ztmp2  = psi_h_ecmwf(z0t*Linv)
            ztmp0  = func_h - ztmp2
            ztmp1  = vkarmn/(LOG(zu) - LOG(z0t) - ztmp0)
            t_star = dt_zu*ztmp1
            ztmp2  = ztmp0 - func_m + ztmp2
            ztmp1  = LOG(zt/zu) + ztmp2
            t_zu   = t_zt - t_star/vkarmn*ztmp1

            ztmp2  = psi_h_ecmwf(z0q*Linv)
            ztmp0  = func_h - ztmp2
            ztmp1  = vkarmn/(LOG(zu) - LOG(z0q) - ztmp0)
            q_star = dq_zu*ztmp1
            ztmp2  = ztmp0 - func_m + ztmp2
            ztmp1  = LOG(zt/zu) + ztmp2
            q_zu   = q_zt - q_star/vkarmn*ztmp1

         END IF

         !! Updating because of updated z0 and z0t and new Linv...
         ztmp0  = zu*Linv
         func_m = LOG(zu) - LOG(z0 ) - psi_m_ecmwf(ztmp0) + psi_m_ecmwf( z0*Linv)
         func_h = LOG(zu) - LOG(z0t) - psi_h_ecmwf(ztmp0) + psi_h_ecmwf(z0t*Linv)

         !! SKIN related part
         !! -----------------
         IF( l_use_skin ) THEN
            !! compute transfer coefficients at zu : lolo: verifier...
            Cd = vkarmn*vkarmn/(func_m*func_m)
            Ch = vkarmn*vkarmn/(func_m*func_h)
            ztmp1 = LOG(zu) - LOG(z0q) - psi_h_ecmwf(ztmp0) + psi_h_ecmwf(z0q*Linv)   ! func_q
            Ce = vkarmn*vkarmn/(func_m*ztmp1)

            ! Non-Solar heat flux to the ocean:
            ztmp1 = U_blk*zrhoa     ! rho*U10
            ztmp2 = T_s*T_s
            ztmp1 = ztmp1 * ( Ce*L0vap*(q_zu - q_s) + Ch*Cp_dry*(t_zu - T_s) ) & ! Total turb. heat flux
               &     + 0.97*(rad_lw - sigma0*ztmp2*ztmp2)                        ! Net longwave flux

            !! Updating the values of the skin temperature T_s and q_s :
            CALL CSWL_ECMWF( zQsw, ztmp1, u_star, zsst, T_s )

            q_s = 0.98*q_sat(MAX(T_s, 200.), slp)  ! 200 -> just to avoid numerics problem on masked regions if silly values are given

         END IF

         dt_zu = t_zu - T_s ;  dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6), dt_zu )
         dq_zu = q_zu - q_s ;  dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9), dq_zu )

      END DO

      Cd = vkarmn*vkarmn/(func_m*func_m)
      Ch = vkarmn*vkarmn/(func_m*func_h)
      ztmp1 = LOG(zu) - LOG(z0q) - psi_h_ecmwf(zu*Linv) + psi_h_ecmwf(z0q*Linv)   ! func_q
      Ce = vkarmn*vkarmn/(func_m*ztmp1)

      IF( lreturn_z0 )    xz0     = z0
      IF( lreturn_ustar ) xu_star = u_star
      IF( lreturn_L )     xL      = 1./Linv
      IF( lreturn_UN10 )  xUN10   = u_star/vkarmn*LOG(10./z0)

      DEALLOCATE ( u_star, t_star, q_star, func_m, func_h, &
         &       dt_zu, dq_zu, z0, z0t, z0q, znu_a, Linv, ztmp0, ztmp1, ztmp2 )

      IF( l_use_skin ) THEN
         DEALLOCATE ( zsst, zrhoa, zQsw ) ! Cool skin
      END IF

   END SUBROUTINE turb_ecmwf


   FUNCTION psi_m_ecmwf( pzeta )
      !!----------------------------------------------------------------------------------
      !! Universal profile stability function for momentum
      !!     ECMWF / as in IFS cy31r1 documentation, available online
      !!     at ecmwf.int
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: psi_m_ecmwf
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) :: zzeta, zx, ztmp, psi_unst, psi_stab, stab
      !!----------------------------------------------------------------------------------
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zzeta = MIN( pzeta(ji,jj) , 5. ) !! Very stable conditions (L positif and big!):
            !
            ! Unstable (Paulson 1970):
            !   eq.3.20, Chap.3, p.33, IFS doc - Cy31r1
            zx = SQRT(ABS(1. - 16.*zzeta))
            ztmp = 1. + SQRT(zx)
            ztmp = ztmp*ztmp
            psi_unst = LOG( 0.125*ztmp*(1. + zx) )   &
               &       -2.*ATAN( SQRT(zx) ) + 0.5*rpi
            !
            ! Unstable:
            ! eq.3.22, Chap.3, p.33, IFS doc - Cy31r1
            psi_stab = -2./3.*(zzeta - 5./0.35)*EXP(-0.35*zzeta) &
               &       - zzeta - 2./3.*5./0.35
            !
            ! Combining:
            stab = 0.5 + SIGN(0.5, zzeta) ! zzeta > 0 => stab = 1
            !
            psi_m_ecmwf(ji,jj) = (1. - stab) * psi_unst & ! (zzeta < 0) Unstable
               &                +      stab  * psi_stab   ! (zzeta > 0) Stable
            !
         END DO
      END DO
      !
   END FUNCTION psi_m_ecmwf


   FUNCTION psi_h_ecmwf( pzeta )
      !!----------------------------------------------------------------------------------
      !! Universal profile stability function for temperature and humidity
      !!     ECMWF / as in IFS cy31r1 documentation, available online
      !!     at ecmwf.int
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: psi_h_ecmwf
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) ::  zzeta, zx, psi_unst, psi_stab, stab
      !!----------------------------------------------------------------------------------
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zzeta = MIN(pzeta(ji,jj) , 5.)   ! Very stable conditions (L positif and big!):
            !
            zx  = ABS(1. - 16.*zzeta)**.25        ! this is actually (1/phi_m)**2  !!!
            !                                     ! eq.3.19, Chap.3, p.33, IFS doc - Cy31r1
            ! Unstable (Paulson 1970) :
            psi_unst = 2.*LOG(0.5*(1. + zx*zx))   ! eq.3.20, Chap.3, p.33, IFS doc - Cy31r1
            !
            ! Stable:
            psi_stab = -2./3.*(zzeta - 5./0.35)*EXP(-0.35*zzeta) & ! eq.3.22, Chap.3, p.33, IFS doc - Cy31r1
               &       - ABS(1. + 2./3.*zzeta)**1.5 - 2./3.*5./0.35 + 1.
            ! LB: added ABS() to avoid NaN values when unstable, which contaminates the unstable solution...
            !
            stab = 0.5 + SIGN(0.5, zzeta) ! zzeta > 0 => stab = 1
            !
            !
            psi_h_ecmwf(ji,jj) = (1. - stab) * psi_unst &   ! (zzeta < 0) Unstable
               &                +    stab    * psi_stab     ! (zzeta > 0) Stable
            !
         END DO
      END DO
      !
   END FUNCTION psi_h_ecmwf


   FUNCTION Ri_bulk( pz, ptz, pdt, pqz, pdq, pub )
      !!----------------------------------------------------------------------------------
      !! Bulk Richardson number (Eq. 3.25 IFS doc)
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: Ri_bulk
      REAL(wp),                     INTENT(in) :: pz       !: height above the sea [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: ptz, &   !: air temperature at pz m [K]
         &                                        pdt, &   !: ptz - sst               [K]
         &                                        pqz, &   !: air temperature at pz m [kg/kg]
         &                                        pdq, &   !: pqz - ssq               [kg/kg]
         &                                        pub      !: bulk wind speed         [m/s]
      !!-------------------------------------------------------------------------------
      !
      Ri_bulk =   grav*pz/(pub*pub)   &
         &      * ( pdt/(ptz - 0.5*(pdt + grav*pz/(Cp_dry + Cp_vap*pqz))) &
         &          + rctv0*pdq )
      !
   END FUNCTION Ri_bulk


   SUBROUTINE CSWL_ECMWF( pQsw, pQnsol, pustar, pSST, pTs )
      !!---------------------------------------------------------------------
      !!
      !!  Cool-Skin Warm-Layer scheme according to Zeng & Beljaars, 2005 (GRL)
      !!  " A prognostic scheme of sea surface skin temperature for modeling and data assimilation "
      !!
      !!    As included in IFS Cy40   /  E.C.M.W.F.
      !!     ------------------------------------------------------------------
      !!
      !!  **   INPUT:
      !!
      !!     *pQsw*       net solar radiative flux to the ocean
      !!     *pQnsol*     net non-solar heat flux to the ocean
      !!     *pustar*     friction velocity u*
      !!     *pSST*       SST
      !!
      !!   **  INPUT/OUTPUT:
      !!     *pTs*  : as input  =>  previous estimate of skin temperature
      !!             as output =>  new estimate of skin temperature
      !!
      !!------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(IN)    :: pQsw
      REAL(wp), DIMENSION(jpi,jpj), INTENT(IN)    :: pQnsol
      REAL(wp), DIMENSION(jpi,jpj), INTENT(IN)    :: pustar
      REAL(wp), DIMENSION(jpi,jpj), INTENT(IN)    :: pSST
      !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(INOUT) :: pTs
      !
      INTEGER :: ji,jj
      !
      REAL(wp) :: rmult, &
         & zRhoCp_w, &
         & ZCON2,ZCON3,ZCON4,ZCON5, zQnsol ,zQnet, zlamb, zdelta,&
         & ZSRD,ZDSST,ZZ,ZEPDU2,&
         & ZFI,zdL,zdL2, ztmp, &
         & ZPHI,ZROADRW, &
         & zus_a, &
         & zdt, zfs, zsgn
      !
      REAL(wp), DIMENSION(jpi,jpj) :: &
         &  zalpha_w, &       !: thermal expansion coefficient of seawater
         & zus_w, zus_w2, &   !: u* and u*^2 in water
         &       zdT_c,   &   !: cool skin temperature increment !lolo rm!
         &       zdT_w        !: warm skin temperature increment
      !
      INTEGER :: nbi, jwl
      !!------------------------------------------------------------------
      !
      !     1. Initialize constants for ocean warm layer and cool skin
      !
      !     1.1 General
      !
      ZEPDU2  = 0.01   !    security constant for velocity**2   (m2/s2)
      ZROADRW = rho0_a/rho0_w          ! Density ratio                      (-)
      zRhoCp_w = rho0_w*Cp0_w
      !
      !     1.2C Warm layer parametrization constants
      !
      !    ZFI = Fraction of solar radiation absorbed in warm layer (-)
      ZFI = 1. -0.28*EXP(-71.5*rd0) -0.27*EXP(-2.8*rd0) - 0.45*EXP(-0.07*rd0)  !: Eq. 8.135
      !
      ZCON3 = rd0*vkarmn*grav/(ZROADRW)**1.5
      ZCON4 = (rNu0 + 1.0)*vkarmn/rd0
      ZCON5 = (rNu0 + 1.0)/(rNu0*rd0)
      !
      !     1.3 Cool skin parametrization constants
      ZCON2 = 16.*grav*zRhoCp_w*nu0_w**3/(k0_w**2)
      !
      ! Friction velocities
      ! "MAX( pustar(:,:), 1.E-4)" is u* in the air !
      zus_w(:,:)  = MAX( pustar(:,:), 1.E-4)*SQRT(ZROADRW)       ! u* in the water
      zus_w2(:,:) = zus_w(:,:)*zus_w(:,:)
      !
      ! Ocean buoyancy
      zalpha_w(:,:) = MAX( 1.E-5 , 1.E-5*(pTs(:,:) - rt0) ) ! thermal expansion coefficient of water
      !
      zdT_c = 0.0
      zdT_w = 0.0
      !
      !  3. Cool skin (Fairall et al. 1996)
      !------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            ! Non-solar heat loss to the atmosphere:
            zQnsol = MAX( 1. , - pQnsol(ji,jj) )

            zlamb = 6.*(1. + (zQnsol*zalpha_w(ji,jj)*ZCON2/(zus_w2(ji,jj)*zus_w2(ji,jj)))**0.75)**(-1./3.)

            zdelta = zlamb*nu0_w/zus_w(ji,jj)

            !    Solar absorption
            zfs   = 0.065 + 11.*zdelta - (6.6E-5/zdelta)*(1. - EXP(-zdelta/8.E-4))  ! Eq. 8.131 / IFS cy40r1, doc, Part IV,
            zfs   = MAX(zfs , 0.01)
            zQnet = MAX( 1. , -zfs*pQsw(ji,jj) + zQnsol )
            zdT_c(ji,jj) = -zdelta*zQnet/k0_w

         END DO
      END DO

      ! 2.2 Warm layer; formulation C (Xubin Zeng)
      !!--------------------------------------------
      nbi = MAX( 1 , nb_itt_wl)

      IF( nbi > 1 ) THEN
         !! Itterating for warm-layer solution
         zdt   = rdt0/REAL(nbi)
         rmult = 1.
      ELSE
         !! No itteration! Way cheaper !!!
         !! Assuming balance between the last 2 terms of Eq. 11 (lhs of 11 = 0):
         !!   => in that case no need for sub-itterations !
         !! Acceptable in most conditions (UNLESS: sunny + no-wind conditions)
         !!
         zdt   = 1.
         rmult = 0
      END IF

      DO jwl = 1, nbi  ! itteration to solve implicitely equation for warm layer

         DO jj = 1, jpj
            DO ji = 1, jpi

               ZDSST = pTs(ji,jj) - pSST(ji,jj) - zdT_c(ji,jj)

               !! Buoyancy flux and stability parameter (zdl = -z/L) in water
               !
               !! Qt/(rho_w*Cpw):
               ZSRD = ( pQsw(ji,jj)*ZFI + pQnsol(ji,jj) )/zRhoCp_w
               !
               zsgn = 0.5 + SIGN(0.5, ZSRD)  ! ZSRD > 0. => 1.  / ZSRD < 0. => 0.
               ztmp = MAX(ZDSST,0.)
               zdl = (zsgn + 1.)*( zus_w2(ji,jj) * SQRT(ztmp/(5.*rd0*grav*zalpha_w(ji,jj)/rNu0)) ) & ! (ZDSST > 0.0 .AND. ZSRD < 0.0)
                  &  +   zsgn   * ZSRD                                                  !   otherwize
               !
               zus_a = MAX( pustar(ji,jj), 1.E-4 )
               zdL = ZCON3*zalpha_w(ji,jj)*zdL/(zus_a*zus_a*zus_a)

               !! Stability function Phi_t(-z/L) (zdL is -z/L) :
               zsgn = 0.5 + SIGN(0.5, zdL)  ! zdl > 0. => 1.  / zdl < 0. => 0.
               zdL2 = zdL*zdL
               ZPHI =     zsgn     * (1. + (5.*zdL + 4.*zdL2)/(1. + 3.*zdL + 0.25*zdL2) ) &  ! (zdL > 0) Takaya et al.
                  &  + (1. + zsgn) * ( 1./SQRT(1. - 16.*(-ABS(zdL))) )        ! (zdl < 0) Eq. 8.136
               !
               !! FOR zdL > 0.0, old relations:
               !         ZPHI = 1.+5.*zdL                                ! Eq. 8.136 (Large et al. 1994)
               !         ZPHI = 1.+5.0*(zdL+zdL**2)/(1.0+3.0*zdL+zdL**2) ! SHEBA, Grachev et al. 2007

               !! Solving 11 by itteration with time step of zdt...
               ZZ = rmult*1. + ZCON4*zdt*zus_w(ji,jj)/ZPHI
               ZZ = SIGN( MAX(ABS(ZZ) , 1.E-4), ZZ )
               zdT_w(ji,jj) = MAX( 0. , (rmult*ZDSST + ZCON5*ZSRD*zdt)/ZZ )

            END DO
         END DO

         ! 3. Apply warm layer and cool skin effects
         !------------------------------------------
         pTs(:,:) = pSST(:,:) + zdT_w(:,:) + zdT_c(:,:)

      END DO  !: sub-itteration

   END SUBROUTINE CSWL_ECMWF


END MODULE mod_blk_ecmwf
