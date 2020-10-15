! AeroBulk / 2019 / L. Brodeau
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
   !!   * air temp. and spec. hum. adjusted from zt (usually 2m) to zu (usually 10m) if needed
   !!   * the "effective" bulk wind speed at zu: Ub
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!    Using the bulk formulation/param. of ECMWF
   !!      + consideration of cool-skin warm layer parametrization (CS: Fairall et al. 1996; WL: Zeng & Beljaars, 2005 )
   !!
   !!       Routine turb_ecmwf maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, July 2019
   !!
   !!====================================================================================
   USE mod_const       !: physical and othe constants
   USE mod_phymbl      !: thermodynamics
   USE mod_skin_ecmwf  !: cool-skin & warm-layer parameterizations of
   !                   !: Zeng and Beljaars, 1995 WITH update from Takaya et al. 2010...

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: ECMWF_INIT, TURB_ECMWF, psi_m_ecmwf, psi_h_ecmwf

   !! ECMWF own values for given constants, taken form IFS documentation...
   REAL(wp), PARAMETER, PUBLIC :: charn0_ecmwf = 0.018_wp    ! Charnock constant (pretty high value here !!!
   !                                          !    =>  Usually 0.011 for moderate winds)
   REAL(wp), PARAMETER ::   zi0     = 1000.   ! scale height of the atmospheric boundary layer...1
   REAL(wp), PARAMETER ::   Beta0    = 1.     ! gustiness parameter ( = 1.25 in COAREv3)
   REAL(wp), PARAMETER ::   alpha_M = 0.11    ! For roughness length (smooth surface term)
   REAL(wp), PARAMETER ::   alpha_H = 0.40    ! (Chapter 3, p.34, IFS doc Cy31r1)
   REAL(wp), PARAMETER ::   alpha_Q = 0.62    !


   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE ecmwf_init(l_use_cs, l_use_wl)
      !!---------------------------------------------------------------------
      !!                  ***  FUNCTION ecmwf_init  ***
      !!
      !! INPUT :
      !! -------
      !!    * l_use_cs : use the cool-skin parameterization
      !!    * l_use_wl : use the warm-layer parameterization
      !!---------------------------------------------------------------------
      LOGICAL , INTENT(in) ::   l_use_cs ! use the cool-skin parameterization
      LOGICAL , INTENT(in) ::   l_use_wl ! use the warm-layer parameterization
      INTEGER :: ierr
      !!---------------------------------------------------------------------
      IF( l_use_wl ) THEN
         ierr = 0
         ALLOCATE ( dT_wl(jpi,jpj), Hz_wl(jpi,jpj), STAT=ierr )
         IF( ierr > 0 ) CALL ctl_stop( ' ECMWF_INIT => allocation of dT_wl & Hz_wl failed!' )
         dT_wl(:,:)  = 0._wp
         Hz_wl(:,:)  = rd0 ! (rd0, constant, = 3m is default for Zeng & Beljaars)
      ENDIF
      IF( l_use_cs ) THEN
         ierr = 0
         ALLOCATE ( dT_cs(jpi,jpj), STAT=ierr )
         IF( ierr > 0 ) CALL ctl_stop( ' ECMWF_INIT => allocation of dT_cs failed!' )
         dT_cs(:,:) = -0.25_wp  ! First guess of skin correction
      ENDIF
   END SUBROUTINE ecmwf_init



   SUBROUTINE turb_ecmwf( kt, zt, zu, T_s, t_zt, q_s, q_zt, U_zu, l_use_cs, l_use_wl,    &
      &                      Cd, Ch, Ce, t_zu, q_zu, Ub,                                 &
      &                      Qsw, rad_lw, slp, pdT_cs,                                   & ! optionals for cool-skin (and warm-layer)
      &                      pdT_wl, pHz_wl,                                             & ! optionals for warm-layer only
      &                      CdN, ChN, CeN, xz0, xu_star, xL, xUN10 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_ecmwf  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to IFS doc. (cycle 45r1)
      !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!                Returns the effective bulk wind speed at zu to be used in the bulk formulas
      !!
      !!                Applies the cool-skin warm-layer correction of the SST to T_s
      !!                if the net shortwave flux at the surface (Qsw), the downwelling longwave
      !!                radiative fluxes at the surface (rad_lw), and the sea-leve pressure (slp)
      !!                are provided as (optional) arguments!
      !!
      !! INPUT :
      !! -------
      !!    *  kt   : current time step (starts at 1)
      !!    *  zt   : height for temperature and spec. hum. of air            [m]
      !!    *  zu   : height for wind speed (usually 10m)                     [m]
      !!    *  t_zt : potential air temperature at zt                         [K]
      !!    *  q_zt : specific humidity of air at zt                          [kg/kg]
      !!    *  U_zu : scalar wind speed at zu                                 [m/s]
      !!    * l_use_cs : use the cool-skin parameterization
      !!    * l_use_wl : use the warm-layer parameterization
      !!
      !! INPUT/OUTPUT:
      !! -------------
      !!    *  T_s  : always "bulk SST" as input                              [K]
      !!              -> unchanged "bulk SST" as output if CSWL not used      [K]
      !!              -> skin temperature as output if CSWL used              [K]
      !!
      !!    *  q_s  : SSQ aka saturation specific humidity at temp. T_s       [kg/kg]
      !!              -> doesn't need to be given a value if skin temp computed (in case l_use_cs=True or l_use_wl=True)
      !!              -> MUST be given the correct value if not computing skint temp. (in case l_use_cs=False or l_use_wl=False)
      !!
      !! OPTIONAL INPUT:
      !! ---------------
      !!    *  Qsw    : net solar flux (after albedo) at the surface (>0)     [W/m^2]
      !!    *  rad_lw : downwelling longwave radiation at the surface  (>0)   [W/m^2]
      !!    *  slp    : sea-level pressure                                    [Pa]
      !!
      !! OPTIONAL OUTPUT:
      !! ----------------
      !!    * pdT_cs  : SST increment "dT" for cool-skin correction           [K]
      !!    * pdT_wl  : SST increment "dT" for warm-layer correction          [K]
      !!    * pHz_wl  : thickness of warm-layer                               [m]
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
      !!    * xz0      : return the aerodynamic roughness length (integration constant for wind stress) [m]
      !!    * xu_star  : return u* the friction velocity                    [m/s]
      !!    * xL       : return the Obukhov length                          [m]
      !!    * xUN10    : neutral wind speed at 10m                          [m/s]
      !!
      !! ** Author: L. Brodeau, June 2019 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      INTEGER,  INTENT(in   )                     ::   kt       ! current time step
      REAL(wp), INTENT(in   )                     ::   zt       ! height for t_zt and q_zt                    [m]
      REAL(wp), INTENT(in   )                     ::   zu       ! height for U_zu                             [m]
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj) ::   T_s      ! sea surface temperature                [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   t_zt     ! potential air temperature              [Kelvin]
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj) ::   q_s      ! sea surface specific humidity           [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   q_zt     ! specific air humidity at zt             [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   U_zu     ! relative wind module at zu                [m/s]
      LOGICAL , INTENT(in   )                     ::   l_use_cs ! use the cool-skin parameterization
      LOGICAL , INTENT(in   )                     ::   l_use_wl ! use the warm-layer parameterization
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Cd       ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ch       ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ce       ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   t_zu     ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   q_zu     ! spec. humidity adjusted at zu           [kg/kg]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ub    ! bulk wind speed at zu                     [m/s]
      !
      REAL(wp), INTENT(in   ), OPTIONAL, DIMENSION(jpi,jpj) ::   Qsw      !             [W/m^2]
      REAL(wp), INTENT(in   ), OPTIONAL, DIMENSION(jpi,jpj) ::   rad_lw   !             [W/m^2]
      REAL(wp), INTENT(in   ), OPTIONAL, DIMENSION(jpi,jpj) ::   slp      !             [Pa]
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   pdT_cs
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   pdT_wl   !             [K]
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   pHz_wl   !             [m]
      !
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   CdN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   ChN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   CeN
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
         &  Linv,            & !: 1/L (inverse of Obukhov length...
         &  z0, z0t, z0q
      !
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zsst  ! to back up the initial bulk SST
      !
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   func_m, func_h
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   ztmp0, ztmp1, ztmp2
      !
      LOGICAL ::  lreturn_cdn=.FALSE., lreturn_chn=.FALSE., lreturn_cen=.FALSE., &
         &        lreturn_z0=.FALSE., lreturn_ustar=.FALSE., lreturn_L=.FALSE., lreturn_UN10=.FALSE.
      CHARACTER(len=40), PARAMETER :: crtnm = 'turb_ecmwf@mod_blk_ecmwf.f90'
      !!----------------------------------------------------------------------------------

      ALLOCATE ( u_star(jpi,jpj), t_star(jpi,jpj), q_star(jpi,jpj), &
         &     func_m(jpi,jpj), func_h(jpi,jpj),  &
         &     dt_zu(jpi,jpj), dq_zu(jpi,jpj),    &
         &     znu_a(jpi,jpj), Linv(jpi,jpj),   &
         &     z0(jpi,jpj), z0t(jpi,jpj), z0q(jpi,jpj), &
         &     ztmp0(jpi,jpj), ztmp1(jpi,jpj), ztmp2(jpi,jpj) )

      IF( kt == nit000 ) CALL ECMWF_INIT(l_use_cs, l_use_wl)

      IF( PRESENT(CdN) )     lreturn_cdn   = .TRUE.
      IF( PRESENT(ChN) )     lreturn_chn   = .TRUE.
      IF( PRESENT(CeN) )     lreturn_cen   = .TRUE.
      IF( PRESENT(xz0) )     lreturn_z0    = .TRUE.
      IF( PRESENT(xu_star) ) lreturn_ustar = .TRUE.
      IF( PRESENT(xL) )      lreturn_L     = .TRUE.
      IF( PRESENT(xUN10) )   lreturn_UN10  = .TRUE.

      l_zt_equal_zu = ( ABS(zu - zt) < 0.01_wp )

      !! Initializations for cool skin and warm layer:
      IF( l_use_cs .AND. (.NOT.(PRESENT(Qsw) .AND. PRESENT(rad_lw) .AND. PRESENT(slp))) ) &
         &   CALL ctl_stop( '['//TRIM(crtnm)//'] => ' , 'you need to provide Qsw, rad_lw & slp to use cool-skin param!' )

      IF( l_use_wl .AND. (.NOT.(PRESENT(Qsw) .AND. PRESENT(rad_lw) .AND. PRESENT(slp))) ) &
         &   CALL ctl_stop( '['//TRIM(crtnm)//'] => ' , 'you need to provide Qsw, rad_lw & slp to use warm-layer param!' )

      IF( l_use_cs .OR. l_use_wl ) THEN
         ALLOCATE ( zsst(jpi,jpj) )
         zsst = T_s ! backing up the bulk SST
         IF( l_use_cs ) T_s = T_s - 0.25_wp   ! First guess of correction
         q_s    = rdct_qsat_salt*q_sat(MAX(T_s, 200._wp), slp) ! First guess of q_s
      ENDIF


      ! Identical first gess as in COARE, with IFS parameter values though...
      !
      !! First guess of temperature and humidity at height zu:
      t_zu = MAX( t_zt ,  180._wp )   ! who knows what's given on masked-continental regions...
      q_zu = MAX( q_zt , 1.e-6_wp )   !               "

      !! Pot. temp. difference (and we don't want it to be 0!)
      dt_zu = t_zu - T_s ;   dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6_wp), dt_zu )
      dq_zu = q_zu - q_s ;   dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9_wp), dq_zu )

      znu_a = visc_air(t_zu) ! Air viscosity (m^2/s) at zt given from temperature in (K)

      Ub = SQRT(U_zu*U_zu + 0.5_wp*0.5_wp) ! initial guess for wind gustiness contribution

      ztmp0   = LOG(    zu*10000._wp) ! optimization: 10000. == 1/z0 (with z0 first guess == 0.0001)
      ztmp1   = LOG(10._wp*10000._wp) !       "                    "               "
      u_star = 0.035_wp*Ub*ztmp1/ztmp0       ! (u* = 0.035*Un10)

      z0     = charn0_ecmwf*u_star*u_star/grav + 0.11_wp*znu_a/u_star
      z0     = MIN( MAX(ABS(z0), 1.E-9) , 1._wp )                      ! (prevents FPE from stupid values from masked region later on)

      z0t    = 1._wp / ( 0.1_wp*EXP(vkarmn/(0.00115/(vkarmn/ztmp1))) )
      z0t    = MIN( MAX(ABS(z0t), 1.E-9) , 1._wp )                      ! (prevents FPE from stupid values from masked region later on)

      Cd     = (vkarmn/ztmp0)**2    ! first guess of Cd

      ztmp0 = vkarmn*vkarmn/LOG(zt/z0t)/Cd

      ztmp2 = Ri_bulk( zu, T_s, t_zu, q_s, q_zu, Ub ) ! Bulk Richardson Number (BRN)

      !! First estimate of zeta_u, depending on the stability, ie sign of BRN (ztmp2):
      ztmp1 = 0.5 + SIGN( 0.5_wp , ztmp2 )
      func_h = (1._wp - ztmp1) *   ztmp0*ztmp2 / (1._wp - ztmp2*zi0*0.004_wp*Beta0**3/zu) & !  BRN < 0
         &  +       ztmp1      * ( ztmp0*ztmp2 + 27._wp/9._wp*ztmp2*ztmp2 )                 !  BRN > 0
      
      !! First guess M-O stability dependent scaling params.(u*,t*,q*) to estimate z0 and z/L
      ztmp0  = vkarmn/(LOG(zu/z0t) - psi_h_ecmwf(func_h))

      u_star = MAX ( Ub*vkarmn/(LOG(zu) - LOG(z0)  - psi_m_ecmwf(func_h)) , 1.E-9 )  !  (MAX => prevents FPE from stupid values from masked region later on)
      t_star = dt_zu*ztmp0
      q_star = dq_zu*ztmp0

      ! What needs to be done if zt /= zu:
      IF( .NOT. l_zt_equal_zu ) THEN
         !! First update of values at zu (or zt for wind)
         ztmp0 = psi_h_ecmwf(func_h) - psi_h_ecmwf(zt*func_h/zu)    ! zt*func_h/zu == zeta_t
         ztmp1 = LOG(zt/zu) + ztmp0
         t_zu = t_zt - t_star/vkarmn*ztmp1
         q_zu = q_zt - q_star/vkarmn*ztmp1
         q_zu = (0.5_wp + SIGN(0.5_wp,q_zu))*q_zu !Makes it impossible to have negative humidity :
         !
         dt_zu = t_zu - T_s  ; dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6_wp), dt_zu )
         dq_zu = q_zu - q_s  ; dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9_wp), dq_zu )
      ENDIF


      !! => that was same first guess as in COARE...


      !! First guess of inverse of Obukov length (1/L) :
      Linv = One_on_L( t_zu, q_zu, u_star, t_star, q_star )

      !! Functions such as  u* = Ub*vkarmn/func_m
      ztmp0 = zu*Linv
      func_m = LOG(zu) - LOG(z0)  - psi_m_ecmwf(ztmp0) + psi_m_ecmwf( z0*Linv)
      func_h = LOG(zu) - LOG(z0t) - psi_h_ecmwf(ztmp0) + psi_h_ecmwf(z0t*Linv)

      !! ITERATION BLOCK
      DO j_itt = 1, nb_itt

         !! Bulk Richardson Number at z=zu (Eq. 3.25)
         ztmp0 = Ri_bulk( zu, T_s, t_zu, q_s, q_zu, Ub ) ! Bulk Richardson Number (BRN)

         !! New estimate of the inverse of the Obukhon length (Linv == zeta/zu) :
         Linv = ztmp0*func_m*func_m/func_h / zu     ! From Eq. 3.23, Chap.3.2.3, IFS doc - Cy40r1
         !! Note: it is slightly different that the L we would get with the usual
         Linv = SIGN( MIN(ABS(Linv),200._wp), Linv ) ! (prevent FPE from stupid values from masked region later on...) !#LOLO

         !! Update func_m with new Linv:
         func_m = LOG(zu) -LOG(z0) - psi_m_ecmwf(zu*Linv) + psi_m_ecmwf(z0*Linv) ! LB: should be "zu+z0" rather than "zu" alone, but z0 is tiny wrt zu!

         !! Need to update roughness lengthes:
         u_star = Ub*vkarmn/func_m
         ztmp2  = u_star*u_star
         ztmp1  = znu_a/u_star
         z0     = MIN( ABS( alpha_M*ztmp1 + charn0_ecmwf*ztmp2/grav ) , 0.001_wp)
         z0t    = MIN( ABS( alpha_H*ztmp1                           ) , 0.001_wp)   ! eq.3.26, Chap.3, p.34, IFS doc - Cy31r1
         z0q    = MIN( ABS( alpha_Q*ztmp1                           ) , 0.001_wp)

         !! Update wind at zu with convection-related wind gustiness in unstable conditions (Chap. 3.2, IFS doc - Cy40r1, Eq.3.17 and Eq.3.18 + Eq.3.8)
         ztmp2 = Beta0*Beta0*ztmp2*(MAX(-zi0*Linv/vkarmn,0._wp))**(2._wp/3._wp) ! square of wind gustiness contribution  (combining Eq. 3.8 and 3.18, hap.3, IFS doc - Cy31r1)
         !!   ! Only true when unstable (L<0) => when ztmp0 < 0 => explains "-" before zi0
         Ub = MAX(SQRT(U_zu*U_zu + ztmp2), 0.2_wp)        ! include gustiness in bulk wind speed
         ! => 0.2 prevents Ub to be 0 in stable case when U_zu=0.


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
         ENDIF

         !! Updating because of updated z0 and z0t and new Linv...
         ztmp0 = zu*Linv
         func_m = log(zu) - LOG(z0 ) - psi_m_ecmwf(ztmp0) + psi_m_ecmwf(z0 *Linv)
         func_h = log(zu) - LOG(z0t) - psi_h_ecmwf(ztmp0) + psi_h_ecmwf(z0t*Linv)


         IF( l_use_cs ) THEN
            !! Cool-skin contribution

            CALL UPDATE_QNSOL_TAU( zu, T_s, q_s, t_zu, q_zu, u_star, t_star, q_star, U_zu, Ub, slp, rad_lw, &
               &                   ztmp1, ztmp0,  Qlat=ztmp2)  ! Qnsol -> ztmp1 / Tau -> ztmp0

            CALL CS_ECMWF( Qsw, ztmp1, u_star, zsst )  ! Qnsol -> ztmp1

            T_s(:,:) = zsst(:,:) + dT_cs(:,:)
            IF( l_use_wl ) T_s(:,:) = T_s(:,:) + dT_wl(:,:)
            q_s(:,:) = rdct_qsat_salt*q_sat(MAX(T_s(:,:), 200._wp), slp(:,:))

         ENDIF

         IF( l_use_wl ) THEN
            !! Warm-layer contribution
            CALL UPDATE_QNSOL_TAU( zu, T_s, q_s, t_zu, q_zu, u_star, t_star, q_star, U_zu, Ub, slp, rad_lw, &
               &                   ztmp1, ztmp2)  ! Qnsol -> ztmp1 / Tau -> ztmp2
            CALL WL_ECMWF( Qsw, ztmp1, u_star, zsst )
            !! Updating T_s and q_s !!!
            T_s(:,:) = zsst(:,:) + dT_wl(:,:)
            IF( l_use_cs ) T_s(:,:) = T_s(:,:) + dT_cs(:,:)
            q_s(:,:) = rdct_qsat_salt*q_sat(MAX(T_s(:,:), 200._wp), slp(:,:))
         ENDIF

         IF( l_use_cs .OR. l_use_wl .OR. (.NOT. l_zt_equal_zu) ) THEN
            dt_zu = t_zu - T_s ;  dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6_wp), dt_zu )
            dq_zu = q_zu - q_s ;  dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9_wp), dq_zu )
         ENDIF

      END DO !DO j_itt = 1, nb_itt

      Cd = vkarmn*vkarmn/(func_m*func_m)
      Ch = vkarmn*vkarmn/(func_m*func_h)
      ztmp2 = log(zu/z0q) - psi_h_ecmwf(zu*Linv) + psi_h_ecmwf(z0q*Linv)   ! func_q
      Ce = vkarmn*vkarmn/(func_m*ztmp2)


      IF( lreturn_cdn .OR. lreturn_chn .OR. lreturn_cen ) ztmp0 = 1._wp/LOG(zu/z0)
      IF( lreturn_cdn )   CdN = vkarmn2*ztmp0*ztmp0
      IF( lreturn_chn )   ChN = vkarmn2*ztmp0/LOG(zu/z0t)
      IF( lreturn_cen )   CeN = vkarmn2*ztmp0/LOG(zu/z0q)

      IF( lreturn_z0 )    xz0     = z0
      IF( lreturn_ustar ) xu_star = u_star
      IF( lreturn_L )     xL      = 1./Linv
      IF( lreturn_UN10 )  xUN10   = u_star/vkarmn*LOG(10./z0)

      DEALLOCATE ( u_star, t_star, q_star, func_m, func_h, &
         &       dt_zu, dq_zu, z0, z0t, z0q, znu_a, Linv, ztmp0, ztmp1, ztmp2 )

      IF( l_use_cs .AND. PRESENT(pdT_cs) ) pdT_cs = dT_cs
      IF( l_use_wl .AND. PRESENT(pdT_wl) ) pdT_wl = dT_wl
      IF( l_use_wl .AND. PRESENT(pHz_wl) ) pHz_wl = Hz_wl

      IF( l_use_cs .OR. l_use_wl ) DEALLOCATE ( zsst )

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
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: psi_m_ecmwf
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) :: zzeta, zx2, zx, ztmp, zpsi_unst, zpsi_stab, zstab, zc
      !!----------------------------------------------------------------------------------
      zc = 5._wp/0.35_wp
      !
      DO jj = 1, jpj
         DO ji = 1, jpi

            zzeta = MIN( pzeta(ji,jj) , 5._wp ) !! Very stable conditions (L positif and big!):

            ! *** Unstable (Paulson 1970)    [eq.3.20, Chap.3, p.33, IFS doc - Cy31r1] :
            zx2 = SQRT( ABS(1._wp - 16._wp*zzeta) )  ! (1 - 16z)^0.5
            zx  = SQRT(zx2)                          ! (1 - 16z)^0.25
            ztmp = 1._wp + zx
            zpsi_unst = LOG( 0.125_wp*ztmp*ztmp*(1._wp + zx2) ) - 2._wp*ATAN( zx ) + 0.5_wp*rpi

            ! *** Stable                   [eq.3.22, Chap.3, p.33, IFS doc - Cy31r1] :
            zpsi_stab = -2._wp/3._wp*(zzeta - zc)*EXP(-0.35_wp*zzeta) &
               &       - zzeta - 2._wp/3._wp*zc
            !
            zstab = 0.5_wp + SIGN(0.5_wp, zzeta) ! zzeta > 0 => zstab = 1
            !
            psi_m_ecmwf(ji,jj) =         zstab  * zpsi_stab &  ! (zzeta > 0) Stable
               &              + (1._wp - zstab) * zpsi_unst    ! (zzeta < 0) Unstable
            !
         END DO
      END DO
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
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: psi_h_ecmwf
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) ::  zzeta, zx2, zpsi_unst, zpsi_stab, zstab, zc
      !!----------------------------------------------------------------------------------
      zc = 5._wp/0.35_wp
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zzeta = MIN(pzeta(ji,jj) , 5._wp)   ! Very stable conditions (L positif and big!):
            !
            ! *** Unstable (Paulson 1970)   [eq.3.20, Chap.3, p.33, IFS doc - Cy31r1] :
            zx2 = SQRT( ABS(1._wp - 16._wp*zzeta) )  ! (1 -16z)^0.5
            zpsi_unst = 2._wp*LOG( 0.5_wp*(1._wp + zx2) )
            !
            ! *** Stable [eq.3.22, Chap.3, p.33, IFS doc - Cy31r1] :
            zpsi_stab = -2._wp/3._wp*(zzeta - zc)*EXP(-0.35_wp*zzeta) &
               &       - ABS(1._wp + 2._wp/3._wp*zzeta)**1.5_wp - 2._wp/3._wp*zc + 1._wp
            !
            ! LB: added ABS() to avoid NaN values when unstable, which contaminates the unstable solution...
            !
            zstab = 0.5_wp + SIGN(0.5_wp, zzeta) ! zzeta > 0 => zstab = 1
            !
            psi_h_ecmwf(ji,jj) =         zstab  * zpsi_stab &  ! (zzeta > 0) Stable
               &              + (1._wp - zstab) * zpsi_unst    ! (zzeta < 0) Unstable            
            !
         END DO
      END DO
   END FUNCTION psi_h_ecmwf


   !!======================================================================
END MODULE mod_blk_ecmwf
