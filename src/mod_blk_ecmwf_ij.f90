! AeroBulk / 2019 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_blk_ecmwf_ij
   !!====================================================================================
   !!       Computes turbulent components of surface fluxes
   !!         according to IFS of the ECMWF
   !!
   !!       With Cool-Skin and Warm-Layer correction of SST (if needed)
   !!
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (usually 2m) to zu (usually 10m) if needed
   !!   * the "effective" bulk wind speed at zu: Ubzu
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!    Using the bulk formulation/param. of ECMWF
   !!      + consideration of cool-skin warm layer parametrization (CS: Fairall et al. 1996; WL: Zeng & Beljaars, 2005 )
   !!
   !!       Routine turb_ecmwf_ij maintained and developed in AeroBulk
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

   PUBLIC :: ECMWF_INIT, TURB_ECMWF_IJ, psi_m_ecmwf_ij, psi_h_ecmwf_ij

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



   SUBROUTINE turb_ecmwf_ij( kt, zt, zu, T_s, t_zt, q_s, q_zt, U_zu, l_use_cs, l_use_wl,    &
      &                      Cd, Ch, Ce, t_zu, q_zu, Ubzu,                                 &
      &                      Qsw, rad_lw, slp, pdT_cs,                                   & ! optionals for cool-skin (and warm-layer)
      &                      pdT_wl, pHz_wl,                                             & ! optionals for warm-layer only
      &                      CdN, ChN, CeN, xz0, xu_star, xL, xUN10 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_ecmwf_ij  ***
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
      !!    *  Ubzu  : bulk wind speed at zu                                 [m/s]
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
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ubzu    ! bulk wind speed at zu                     [m/s]
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
      INTEGER :: ji, jj, jit
      LOGICAL :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zsst  ! to back up the initial bulk SST
      !
      REAL(wp) :: zdt, zdq, zus, zts, zqs, zNu_a, zRib, z1oL, zpsi_m, zpsi_h, zUn10, z1oL
      REAL(wp) :: zz0, zz0t, zz0q, zprof_m, zprof_h
      !
      LOGICAL ::  lreturn_cdn=.FALSE., lreturn_chn=.FALSE., lreturn_cen=.FALSE., &
         &        lreturn_z0=.FALSE., lreturn_ustar=.FALSE., lreturn_L=.FALSE., lreturn_UN10=.FALSE.
      CHARACTER(len=40), PARAMETER :: crtnm = 'turb_ecmwf_ij@mod_blk_ecmwf_ij.f90'
      !!----------------------------------------------------------------------------------

      IF( kt == nit000 ) CALL ECMWF_INIT(l_use_cs, l_use_wl)

      lreturn_cdn   =  PRESENT(CdN)
      lreturn_chn   =  PRESENT(ChN)
      lreturn_cen   =  PRESENT(CeN)
      lreturn_z0    =  PRESENT(xz0)
      lreturn_ustar =  PRESENT(xu_star)
      lreturn_L     =  PRESENT(xL)
      lreturn_UN10  =  PRESENT(xUN10)

      l_zt_equal_zu = ( ABS(zu - zt) < 0.01_wp )

      !! Initializations for cool skin and warm layer:
      IF( l_use_cs .AND. (.NOT.(PRESENT(Qsw) .AND. PRESENT(rad_lw) .AND. PRESENT(slp))) ) &
         &   CALL ctl_stop( '['//TRIM(crtnm)//'] => ' , 'you need to provide Qsw, rad_lw & slp to use cool-skin param!' )

      IF( l_use_wl .AND. (.NOT.(PRESENT(Qsw) .AND. PRESENT(rad_lw) .AND. PRESENT(slp))) ) &
         &   CALL ctl_stop( '['//TRIM(crtnm)//'] => ' , 'you need to provide Qsw, rad_lw & slp to use warm-layer param!' )

      IF( l_use_cs .OR. l_use_wl ) THEN
         ALLOCATE ( zsst(jpi,jpj) )
         zsst = T_s ! backing up the bulk SST
         IF( l_use_cs ) T_s = T_s - 0.25_wp   ! Crude first guess for skin correction
         q_s    = rdct_qsat_salt*q_sat(MAX(T_s, 200._wp), slp) ! First guess of q_s
      ENDIF


      !! First guess of temperature and humidity at height zu:
      t_zu = MAX( t_zt ,  180._wp )   ! who knows what's given on masked-continental regions...
      q_zu = MAX( q_zt , 1.e-6_wp )   !               "


      DO jj = 1, jpj
         DO ji = 1, jpi

            ! Identical first gess as in COARE, with IFS parameter values though...

            !! Pot. temp. difference (and we don't want it to be 0!)
            zdt = t_zu(ji,jj) - T_s(ji,jj) ;   zdt = SIGN( MAX(ABS(zdt),1.E-6_wp), zdt )
            zdq = q_zu(ji,jj) - q_s(ji,jj) ;   zdq = SIGN( MAX(ABS(zdq),1.E-9_wp), zdq )

            zNu_a = visc_air(t_zu(ji,jj)) ! Air viscosity (m^2/s) at zt given from temperature in (K)

            Ubzu(ji,jj) = SQRT(U_zu(ji,jj)*U_zu(ji,jj) + 0.5_wp*0.5_wp) ! initial guess for wind gustiness contribution

            ztmp0   = LOG(    zu*10000._wp) ! optimization: 10000. == 1/z0 (with z0 first guess == 0.0001)
            ztmp1   = LOG(10._wp*10000._wp) !       "                    "               "
            zus = 0.035_wp*Ubzu(ji,jj)*ztmp1/ztmp0       ! (u* = 0.035*Un10)

            zz0     = charn0_ecmwf*zus*zus/grav + 0.11_wp*zNu_a/zus
            zz0     = MIN( MAX(ABS(zz0), 1.E-9) , 1._wp )                      ! (prevents FPE from stupid values from masked region later on)

            zz0t    = 1._wp / ( 0.1_wp*EXP(vkarmn/(0.00115/(vkarmn/ztmp1))) )
            zz0t    = MIN( MAX(ABS(zz0t), 1.E-9) , 1._wp )                      ! (prevents FPE from stupid values from masked region later on)

            Cd(ji,jj)     = MAX( (vkarmn/ztmp0)**2 , Cx_min )   ! first guess of Cd

            ztmp0 = vkarmn2/LOG(zt/zz0t)/Cd(ji,jj)

            zRib = Ri_bulk( zu, T_s(ji,jj), t_zu(ji,jj), q_s(ji,jj), q_zu(ji,jj), Ubzu(ji,jj) ) ! Bulk Richardson Number (BRN)

            !! First estimate of zeta_u, depending on the stability, ie sign of BRN (zRi):
            ztmp1 = 0.5 + SIGN( 0.5_wp , zRi )
            zprof_h = (1._wp - ztmp1) *   ztmp0*zRi / (1._wp - zRi*zi0*0.004_wp*Beta0**3/zu) & !  BRN < 0
               &  +       ztmp1      * ( ztmp0*zRi + 27._wp/9._wp*zRi*zRi )                 !  BRN > 0

            !! First guess M-O stability dependent scaling params.(u*,t*,q*) to estimate z0 and z/L
            ztmp0  = vkarmn/(LOG(zu/zz0t) - psi_h_ecmwf_ij(zprof_h))

            zus = MAX ( Ubzu(ji,jj)*vkarmn/(LOG(zu) - LOG(zz0)  - psi_m_ecmwf_ij(zprof_h)) , 1.E-9 )  !  (MAX => prevents FPE from stupid values from masked region later on)
            zts = zdt*ztmp0
            zqs = zdq*ztmp0

            ! What needs to be done if zt /= zu:
            IF( .NOT. l_zt_equal_zu ) THEN
               !! First update of values at zu (or zt for wind)
               ztmp0 = psi_h_ecmwf_ij(zprof_h) - psi_h_ecmwf_ij(zt*zprof_h/zu)    ! zt*zprof_h/zu == zeta_t
               ztmp1 = LOG(zt/zu) + ztmp0
               t_zu(ji,jj) = t_zt - zts/vkarmn*ztmp1
               q_zu(ji,jj) = q_zt - zqs/vkarmn*ztmp1
               q_zu(ji,jj) = (0.5_wp + SIGN(0.5_wp,q_zu(ji,jj)))*q_zu(ji,jj) !Makes it impossible to have negative humidity :
               !
               zdt = t_zu(ji,jj) - T_s(ji,jj)  ; zdt = SIGN( MAX(ABS(zdt),1.E-6_wp), zdt )
               zdq = q_zu(ji,jj) - q_s(ji,jj)  ; zdq = SIGN( MAX(ABS(zdq),1.E-9_wp), zdq )
            ENDIF


            !! => that was same first guess as in COARE...


            !! First guess of inverse of Obukov length (1/L) :
            z1oL = One_on_L( t_zu(ji,jj), q_zu(ji,jj), zus, zts, zqs )

            !! Functions such as  u* = Ubzu*vkarmn/zprof_m
            ztmp0 = zu*z1oL
            zprof_m = LOG(zu) - LOG(zz0)  - psi_m_ecmwf_ij(ztmp0) + psi_m_ecmwf_ij( zz0*z1oL)
            zprof_h = LOG(zu) - LOG(zz0t) - psi_h_ecmwf_ij(ztmp0) + psi_h_ecmwf_ij(zz0t*z1oL)

            !! ITERATION BLOCK
            DO jit = 1, nb_iter

               !! Bulk Richardson Number at z=zu (Eq. 3.25)
               ztmp0 = Ri_bulk( zu, T_s(ji,jj), t_zu(ji,jj), q_s(ji,jj), q_zu(ji,jj), Ubzu(ji,jj) ) ! Bulk Richardson Number (BRN)

               !! New estimate of the inverse of the Obukhon length (z1oL == zeta/zu) :
               z1oL = ztmp0*zprof_m*zprof_m/zprof_h / zu     ! From Eq. 3.23, Chap.3.2.3, IFS doc - Cy40r1
               !! Note: it is slightly different that the L we would get with the usual
               z1oL = SIGN( MIN(ABS(z1oL),200._wp), z1oL ) ! (prevent FPE from stupid values from masked region later on...)

               !! Update zprof_m with new z1oL:
               zprof_m = LOG(zu) -LOG(zz0) - psi_m_ecmwf_ij(zu*z1oL) + psi_m_ecmwf_ij(zz0*z1oL) ! LB: should be "zu+z0" rather than "zu" alone, but z0 is tiny wrt zu!

               !! Need to update roughness lengthes:
               zus = Ubzu(ji,jj)*vkarmn/zprof_m
               ztmp2  = zus*zus !
               ztmp1  = zNu_a/zus
               zz0     = MIN( ABS( alpha_M*ztmp1 + charn0_ecmwf*ztmp2/grav ) , 0.001_wp)
               zz0t    = MIN( ABS( alpha_H*ztmp1                           ) , 0.001_wp)   ! eq.3.26, Chap.3, p.34, IFS doc - Cy31r1
               zz0q    = MIN( ABS( alpha_Q*ztmp1                           ) , 0.001_wp)

               !! Update wind at zu with convection-related wind gustiness in unstable conditions (Chap. 3.2, IFS doc - Cy40r1, Eq.3.17 and Eq.3.18 + Eq.3.8)
               ztmp2 = Beta0*Beta0*ztmp2*(MAX(-zi0*z1oL/vkarmn,0._wp))**(2._wp/3._wp) ! square of wind gustiness contribution  (combining Eq. 3.8 and 3.18, hap.3, IFS doc - Cy31r1)
               !!   ! Only true when unstable (L<0) => when ztmp0 < 0 => explains "-" before zi0
               Ubzu(ji,jj) = MAX(SQRT(U_zu(ji,jj)*U_zu(ji,jj) + ztmp2), 0.2_wp)        ! include gustiness in bulk wind speed
               ! => 0.2 prevents Ubzu to be 0 in stable case when U_zu=0.


               !! Need to update "theta" and "q" at zu in case they are given at different heights
               !! as well the air-sea differences:
               IF( .NOT. l_zt_equal_zu ) THEN
                  !! Arrays zprof_m and zprof_h are free for a while so using them as temporary arrays...
                  zprof_h = psi_h_ecmwf_ij(zu*z1oL) ! temporary array !!!
                  zprof_m = psi_h_ecmwf_ij(zt*z1oL) ! temporary array !!!

                  ztmp2  = psi_h_ecmwf_ij(zz0t*z1oL)
                  ztmp0  = zprof_h - ztmp2
                  ztmp1  = vkarmn/(LOG(zu) - LOG(zz0t) - ztmp0)
                  zts = zdt*ztmp1
                  ztmp2  = ztmp0 - zprof_m + ztmp2
                  ztmp1  = LOG(zt/zu) + ztmp2
                  t_zu(ji,jj)   = t_zt - zts/vkarmn*ztmp1

                  ztmp2  = psi_h_ecmwf_ij(zz0q*z1oL)
                  ztmp0  = zprof_h - ztmp2
                  ztmp1  = vkarmn/(LOG(zu) - LOG(zz0q) - ztmp0)
                  zqs = zdq*ztmp1
                  ztmp2  = ztmp0 - zprof_m + ztmp2
                  ztmp1  = LOG(zt/zu) + ztmp2
                  q_zu(ji,jj)   = q_zt - zqs/vkarmn*ztmp1
               ENDIF

               !! Updating because of updated z0 and z0t and new z1oL...
               ztmp0 = zu*z1oL
               zprof_m = log(zu) - LOG(zz0 ) - psi_m_ecmwf_ij(ztmp0) + psi_m_ecmwf_ij(zz0 *z1oL)
               zprof_h = log(zu) - LOG(zz0t) - psi_h_ecmwf_ij(ztmp0) + psi_h_ecmwf_ij(zz0t*z1oL)


               IF( l_use_cs ) THEN
                  !! Cool-skin contribution

                  CALL UPDATE_QNSOL_TAU( zu, T_s(ji,jj), q_s(ji,jj), t_zu(ji,jj), q_zu(ji,jj), zus, zts, zqs, U_zu(ji,jj), Ubzu(ji,jj), slp, rad_lw, &
                     &                   ztmp1, ztmp0,  Qlat=ztmp2)  ! Qnsol -> ztmp1 / Tau -> ztmp0

                  CALL CS_ECMWF( Qsw, ztmp1, zus, zsst )  ! Qnsol -> ztmp1

                  T_s(ji,jj) = zsst(:,:) + dT_cs(:,:)
                  IF( l_use_wl ) T_s(ji,jj) = T_s(ji,jj) + dT_wl(:,:)
                  q_s(ji,jj) = rdct_qsat_salt*q_sat(MAX(T_s(ji,jj), 200._wp), slp(:,:))

               ENDIF

               IF( l_use_wl ) THEN
                  !! Warm-layer contribution
                  CALL UPDATE_QNSOL_TAU( zu, T_s(ji,jj), q_s(ji,jj), t_zu(ji,jj), q_zu(ji,jj), zus, zts, zqs, U_zu(ji,jj), Ubzu(ji,jj), slp, rad_lw, &
                     &                   ztmp1, ztmp2)  ! Qnsol -> ztmp1 / Tau -> ztmp2
                  CALL WL_ECMWF( Qsw, ztmp1, zus, zsst )
                  !! Updating T_s(ji,jj) and q_s !!!
                  T_s(ji,jj) = zsst(:,:) + dT_wl(:,:) !
                  IF( l_use_cs ) T_s(ji,jj) = T_s(ji,jj) + dT_cs(:,:)
                  q_s(ji,jj) = rdct_qsat_salt*q_sat(MAX(T_s(ji,jj), 200._wp), slp(:,:))
               ENDIF

               IF( l_use_cs .OR. l_use_wl .OR. (.NOT. l_zt_equal_zu) ) THEN
                  zdt = t_zu(ji,jj) - T_s(ji,jj) ;  zdt = SIGN( MAX(ABS(zdt),1.E-6_wp), zdt )
                  zdq = q_zu(ji,jj) - q_s(ji,jj) ;  zdq = SIGN( MAX(ABS(zdq),1.E-9_wp), zdq )
               ENDIF

            END DO !DO jit = 1, nb_iter

            Cd(ji,jj) = MAX( vkarmn2/(zprof_m*zprof_m) , Cx_min )
            Ch(ji,jj) = MAX( vkarmn2/(zprof_m*zprof_h) , Cx_min )
            ztmp2     = LOG(zu/zz0q) - psi_h_ecmwf_ij(zu*z1oL) + psi_h_ecmwf_ij(zz0q*z1oL)   ! zprof_q
            Ce(ji,jj) = MAX( vkarmn2/(zprof_m*ztmp2)  , Cx_min )

            IF( lreturn_cdn .OR. lreturn_chn .OR. lreturn_cen ) ztmp0 = 1._wp/LOG(zu/zz0)
            IF( lreturn_cdn )   CdN(ji,jj) = MAX( vkarmn2*ztmp0*ztmp0       , Cx_min )
            IF( lreturn_chn )   ChN(ji,jj) = MAX( vkarmn2*ztmp0/LOG(zu/zz0t) , Cx_min )
            IF( lreturn_cen )   CeN(ji,jj) = MAX( vkarmn2*ztmp0/LOG(zu/zz0q) , Cx_min )

            IF( lreturn_z0 )    xz0(ji,jj)     = zz0
            IF( lreturn_ustar ) xu_star(ji,jj) = zus
            IF( lreturn_L )     xL(ji,jj)      = 1./z1oL
            IF( lreturn_UN10 )  xUN10(ji,jj)   = zus/vkarmn*LOG(10./zz0)
            
            IF( l_use_cs .AND. PRESENT(pdT_cs) ) pdT_cs = dT_cs
            IF( l_use_wl .AND. PRESENT(pdT_wl) ) pdT_wl = dT_wl
            IF( l_use_wl .AND. PRESENT(pHz_wl) ) pHz_wl = Hz_wl

         END DO
      END DO

      IF( l_use_cs .OR. l_use_wl ) DEALLOCATE ( zsst )

   END SUBROUTINE turb_ecmwf_ij


   FUNCTION psi_m_ecmwf_ij( pzeta )
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
      REAL(wp), DIMENSION(jpi,jpj) :: psi_m_ecmwf_ij
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) :: zta, zx2, zx, ztmp, zpsi_unst, zpsi_stab, zstab, zc
      !!----------------------------------------------------------------------------------
      zc = 5._wp/0.35_wp
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zta = MIN( pzeta(ji,jj) , 5._wp ) !! Very stable conditions (L positif and big!):

            ! *** Unstable (Paulson 1970)    [eq.3.20, Chap.3, p.33, IFS doc - Cy31r1] :
            zx2 = SQRT( ABS(1._wp - 16._wp*zta) )  ! (1 - 16z)^0.5
            zx  = SQRT(zx2)                          ! (1 - 16z)^0.25
            ztmp = 1._wp + zx
            zpsi_unst = LOG( 0.125_wp*ztmp*ztmp*(1._wp + zx2) ) - 2._wp*ATAN( zx ) + 0.5_wp*rpi !

            ! *** Stable                   [eq.3.22, Chap.3, p.33, IFS doc - Cy31r1] :
            zpsi_stab = -2._wp/3._wp*(zta - zc)*EXP(-0.35_wp*zta) &
               &       - zta - 2._wp/3._wp*zc
            !
            zstab = 0.5_wp + SIGN(0.5_wp, zta) ! zta > 0 => zstab = 1
            !
            psi_m_ecmwf_ij(ji,jj) =         zstab  * zpsi_stab &  ! (zta > 0) Stable
               &              + (1._wp - zstab) * zpsi_unst    ! (zta < 0) Unstable
            !
         END DO
      END DO
   END FUNCTION psi_m_ecmwf_ij


   FUNCTION psi_h_ecmwf_ij( pzeta )
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
      REAL(wp), DIMENSION(jpi,jpj) :: psi_h_ecmwf_ij
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) ::  zta, zx2, zpsi_unst, zpsi_stab, zstab, zc
      !!----------------------------------------------------------------------------------
      zc = 5._wp/0.35_wp
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zta = MIN(pzeta(ji,jj) , 5._wp)   ! Very stable conditions (L positif and big!):
            !
            ! *** Unstable (Paulson 1970)   [eq.3.20, Chap.3, p.33, IFS doc - Cy31r1] :
            zx2 = SQRT( ABS(1._wp - 16._wp*zta) )  ! (1 -16z)^0.5
            zpsi_unst = 2._wp*LOG( 0.5_wp*(1._wp + zx2) )
            !
            ! *** Stable [eq.3.22, Chap.3, p.33, IFS doc - Cy31r1] :
            zpsi_stab = -2._wp/3._wp*(zta - zc)*EXP(-0.35_wp*zta) &
               &       - ABS(1._wp + 2._wp/3._wp*zta)**1.5_wp - 2._wp/3._wp*zc + 1._wp
            !
            ! LB: added ABS() to avoid NaN values when unstable, which contaminates the unstable solution...
            !
            zstab = 0.5_wp + SIGN(0.5_wp, zta) ! zta > 0 => zstab = 1
            !
            psi_h_ecmwf_ij(ji,jj) =         zstab  * zpsi_stab &  ! (zta > 0) Stable
               &              + (1._wp - zstab) * zpsi_unst    ! (zta < 0) Unstable
            !
         END DO
      END DO
   END FUNCTION psi_h_ecmwf_ij


   !!======================================================================
END MODULE mod_blk_ecmwf_ij
