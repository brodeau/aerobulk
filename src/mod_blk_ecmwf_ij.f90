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
   USE mod_skin_ecmwf_ij  !: cool-skin & warm-layer parameterizations of
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
   END SUBROUTINE ecmwf_init


   
   SUBROUTINE turb_ecmwf_ij( kt, zt, zu, T_s, t_zt, q_s, q_zt, U_zu, l_use_cs, l_use_wl, &
      &                      Cd, Ch, Ce, t_zu, q_zu, Ubzu,                               &
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
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xL     ! Obukhov length [m]
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xUN10  ! Neutral wind at zu
      !
      INTEGER :: ji, jj, jit
      LOGICAL :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zSST     ! to back up the initial bulk SST
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zu_star, zt_star, zq_star
      !
      REAL(wp) :: zdt, zdq, zus, zts, zqs, zNu_a, zRib, z1oL, zpsi_m_u, zpsi_h_u, zpsi_h_t, zdT_cs, zgust2
      REAL(wp) :: zz0, zz0t, zz0q, zprof_m, zprof_h, zpsi_h_z0t, zpsi_h_z0q, zzeta_u, zzeta_t
      REAL(wp) :: ztmp0, ztmp1, ztmp2
      REAL(wp) :: zlog_z0, zlog_zu, zlog_zt_o_zu
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
         ALLOCATE ( zSST(jpi,jpj) )
         zSST = T_s ! backing up the bulk SST
         IF( l_use_cs ) T_s = T_s - 0.25_wp   ! Crude first guess for skin correction
         q_s    = rdct_qsat_salt*q_sat(MAX(T_s, 200._wp), slp) ! First guess of q_s
      ENDIF

      ALLOCATE ( zu_star(jpi,jpj), zt_star(jpi,jpj), zq_star(jpi,jpj) )


      !! Constants:
      zlog_zu      = LOG(zu)
      zlog_zt_o_zu = LOG(zt/zu)


      CALL FIRST_GUESS_COARE( zt, zu, T_s, t_zt, q_s, q_zt, U_zu, U_zu*0.+charn0_ecmwf,  &
         &                    zu_star, zt_star, zq_star, t_zu, q_zu, Ubzu )
      PRINT *, ''
      PRINT *, 'LOLO: apres FIRST_GUESS_COARE turb_ecmwf_ij: zu_star =', zu_star
      !+ check t_zu and q_zu !!!


      DO jj = 1, jpj
         DO ji = 1, jpi

            !! Pot. temp. difference (and we don't want it to be 0!)
            zdt = t_zu(ji,jj) - T_s(ji,jj) ;   zdt = SIGN( MAX(ABS(zdt),1.E-6_wp), zdt )
            zdq = q_zu(ji,jj) - q_s(ji,jj) ;   zdq = SIGN( MAX(ABS(zdq),1.E-9_wp), zdq )

            !! first guess turb. scales:
            zus = zu_star(ji,jj)
            zts = zt_star(ji,jj)
            zqs = zq_star(ji,jj)

            !! First guess of inverse of Obukov length (1/L) :
            z1oL = One_on_L( t_zu(ji,jj), q_zu(ji,jj), zus, zts, zqs )

            zNu_a = visc_air(t_zu(ji,jj)) ! Air viscosity (m^2/s) at zt given from temperature in (K)

            CALL update_z0_ecmwf( zus, zNu_a,  zz0, zz0t, zz0q )
            zlog_z0 = LOG(zz0)
            
            !! Functions such as  u* = Ubzu*vkarmn/zprof_m
            zzeta_u = zu*z1oL
            zprof_m = zlog_zu - zlog_z0   - psi_m_ecmwf_ij(zzeta_u) + psi_m_ecmwf_ij(zz0 *z1oL)
            zprof_h = zlog_zu - LOG(zz0t) - psi_h_ecmwf_ij(zzeta_u) + psi_h_ecmwf_ij(zz0t*z1oL)


            !! ITERATION BLOCK
            DO jit = 1, nb_iter

               !! Bulk Richardson Number at z=zu (Eq. 3.25)
               zRib = Ri_bulk( zu, T_s(ji,jj), t_zu(ji,jj), q_s(ji,jj), q_zu(ji,jj), Ubzu(ji,jj) ) ! Bulk Richardson Number (BRN)

               !! New estimate of the inverse of the Obukhon length (z1oL == zeta/zu) :
               z1oL = zRib*zprof_m*zprof_m/zprof_h / zu     ! From Eq. 3.23, Chap.3.2.3, IFS doc - Cy40r1
               !! Note: it is slightly different that the L we would get with the usual
               z1oL = SIGN( MIN(ABS(z1oL),200._wp), z1oL ) ! (prevent FPE from stupid values from masked region later on...)
               zzeta_u = zu*z1oL
               zzeta_t = zt*z1oL

               !! Update stability corrections:
               zpsi_m_u = psi_m_ecmwf_ij(zzeta_u)
               zpsi_h_u = psi_h_ecmwf_ij(zzeta_u)
               zpsi_h_t = psi_h_ecmwf_ij(zzeta_t)

               !! Update zprof_m with new z1oL:
               zprof_m = zlog_zu - zlog_z0 - zpsi_m_u + psi_m_ecmwf_ij(zz0*z1oL) ! LB: should be "zu+z0" rather than "zu" alone, but z0 is tiny wrt zu!

               !! Need to update roughness lengthes:
               zus = Ubzu(ji,jj)*vkarmn/zprof_m
               CALL update_z0_ecmwf( zus, zNu_a,  zz0, zz0t, zz0q )
               zlog_z0 = LOG(zz0)

               !! Update wind at zu with convection-related wind gustiness in unstable conditions (Chap. 3.2, IFS doc - Cy40r1, Eq.3.17 and Eq.3.18 + Eq.3.8)
               zgust2 = Beta0*Beta0*zus*zus*(MAX(-zi0*z1oL/vkarmn,0._wp))**(2._wp/3._wp) ! square of wind gustiness contribution  (combining Eq. 3.8 and 3.18, hap.3, IFS doc - Cy31r1)
               !!   ! Only true when unstable (L<0) => when zRib < 0 => explains "-" before zi0
               Ubzu(ji,jj) = MAX(SQRT(U_zu(ji,jj)*U_zu(ji,jj) + zgust2), 0.2_wp)        ! include gustiness in bulk wind speed
               ! => 0.2 prevents Ubzu to be 0 in stable case when U_zu=0.


               !! Need to update "theta" and "q" at zu in case they are given at different heights
               !! as well the air-sea differences:
               IF( .NOT. l_zt_equal_zu ) THEN
                  !! Arrays zprof_m and zprof_h are free for a while so using them as temporary arrays...

                  zpsi_h_z0t  = psi_h_ecmwf_ij(zz0t*z1oL)
                  ztmp0  = zpsi_h_u - zpsi_h_z0t
                  ztmp1  = vkarmn/(zlog_zu - LOG(zz0t) - ztmp0)
                  zts = zdt*ztmp1
                  ztmp2  = ztmp0 - zpsi_h_t + zpsi_h_z0t
                  ztmp1  = zlog_zt_o_zu + ztmp2
                  t_zu(ji,jj)   = t_zt(ji,jj) - zts/vkarmn*ztmp1

                  zpsi_h_z0q  = psi_h_ecmwf_ij(zz0q*z1oL)
                  ztmp0  = zpsi_h_u - zpsi_h_z0q
                  ztmp1  = vkarmn/(zlog_zu - LOG(zz0q) - ztmp0)
                  zqs = zdq*ztmp1
                  ztmp2  = ztmp0 - zpsi_h_t + zpsi_h_z0q
                  ztmp1  = zlog_zt_o_zu + ztmp2
                  q_zu(ji,jj)   = q_zt(ji,jj) - zqs/vkarmn*ztmp1
               ENDIF

               !! Updating because of updated z0 and z0t and new zzeta_u...
               zprof_m = zlog_zu - zlog_z0   - zpsi_m_u + psi_m_ecmwf_ij(zz0 *z1oL)
               zprof_h = zlog_zu - LOG(zz0t) - zpsi_h_u + psi_h_ecmwf_ij(zz0t*z1oL)


               IF( l_use_cs ) THEN
                  !! Cool-skin contribution

                  CALL UPDATE_QNSOL_TAU( zu, T_s(ji,jj), q_s(ji,jj), t_zu(ji,jj), q_zu(ji,jj), &
                     &                   zus, zts, zqs, U_zu(ji,jj), Ubzu(ji,jj), slp(ji,jj), rad_lw(ji,jj), &
                     &                   ztmp1, ztmp0,  Qlat=ztmp2)  ! Qnsol -> ztmp1 / Tau -> ztmp0    !LOLO: fix ztmp0 is not what it is supposed to be !!!

                  CALL CS_ECMWF_IJ( Qsw(ji,jj), ztmp1, zus, zSST(ji,jj), zdT_cs )  ! Qnsol -> ztmp1

                  T_s(ji,jj) = zSST(ji,jj) + zdT_cs
                  IF( l_use_wl ) T_s(ji,jj) = T_s(ji,jj) + dT_wl(ji,jj)
                  q_s(ji,jj) = rdct_qsat_salt*q_sat(MAX(T_s(ji,jj), 200._wp), slp(ji,jj))

               ENDIF

               IF( l_use_wl ) THEN
                  !! Warm-layer contribution
                  CALL UPDATE_QNSOL_TAU( zu, T_s(ji,jj), q_s(ji,jj), t_zu(ji,jj), q_zu(ji,jj), &
                     &                   zus, zts, zqs, U_zu(ji,jj), Ubzu(ji,jj), slp(ji,jj), rad_lw(ji,jj), &
                     &                   ztmp1, ztmp2)  ! Qnsol -> ztmp1 / Tau -> ztmp2
                  CALL WL_ECMWF_IJ( ji, jj, Qsw(ji,jj), ztmp1, zus, zSST(ji,jj) )
                  !! Updating T_s(ji,jj) and q_s !!!
                  T_s(ji,jj) = zSST(ji,jj) + dT_wl(ji,jj) !
                  IF( l_use_cs ) T_s(ji,jj) = T_s(ji,jj) + zdT_cs
                  q_s(ji,jj) = rdct_qsat_salt*q_sat(MAX(T_s(ji,jj), 200._wp), slp(ji,jj))
               ENDIF

               IF( l_use_cs .OR. l_use_wl .OR. (.NOT. l_zt_equal_zu) ) THEN
                  zdt = t_zu(ji,jj) - T_s(ji,jj) ;  zdt = SIGN( MAX(ABS(zdt),1.E-6_wp), zdt )
                  zdq = q_zu(ji,jj) - q_s(ji,jj) ;  zdq = SIGN( MAX(ABS(zdq),1.E-9_wp), zdq )
               ENDIF

               
               PRINT *, 'LOLO:                 zus =', zus
               
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

            IF( l_use_cs .AND. PRESENT(pdT_cs) ) pdT_cs(ji,jj) = zdT_cs

            !! Remove:
            zu_star(ji,jj) = zus !LOLO
            
         END DO
      END DO

      IF( l_use_wl .AND. PRESENT(pdT_wl) ) pdT_wl(:,:) = dT_wl(:,:)
      IF( l_use_wl .AND. PRESENT(pHz_wl) ) pHz_wl(:,:) = Hz_wl(:,:)

      IF( l_use_cs .OR. l_use_wl ) DEALLOCATE ( zSST )



      PRINT *, 'LOLO: fin turb_ecmwf_ij: zu_star =', zu_star, '(', nb_iter, 'iterations)'
      PRINT *, ''
      
      DEALLOCATE ( zu_star, zt_star, zq_star )

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
      REAL(wp) :: psi_m_ecmwf_ij
      REAL(wp), INTENT(in) :: pzeta
      !
      REAL(wp) :: zta, zx2, zx, ztmp, zpsi_unst, zpsi_stab, zstab, zc
      !!----------------------------------------------------------------------------------
      zc = 5._wp/0.35_wp
      !
      zta = MIN( pzeta , 5._wp ) !! Very stable conditions (L positif and big!):

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
      psi_m_ecmwf_ij =         zstab  * zpsi_stab &  ! (zta > 0) Stable
         &              + (1._wp - zstab) * zpsi_unst    ! (zta < 0) Unstable
      !
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
      REAL(wp) :: psi_h_ecmwf_ij
      REAL(wp), INTENT(in) :: pzeta
      !
      REAL(wp) ::  zta, zx2, zpsi_unst, zpsi_stab, zstab, zc
      !!----------------------------------------------------------------------------------
      zc = 5._wp/0.35_wp
      !
      zta = MIN(pzeta , 5._wp)   ! Very stable conditions (L positif and big!):
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
      psi_h_ecmwf_ij =         zstab  * zpsi_stab &  ! (zta > 0) Stable
         &              + (1._wp - zstab) * zpsi_unst    ! (zta < 0) Unstable
      !
   END FUNCTION psi_h_ecmwf_ij



   SUBROUTINE update_z0_ecmwf( pus, pnu,  pz0, pz0t, pz0q )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  update_z0_ecmwf  ***
      !!
      !! ** Purpose :
      !!
      !! INPUT :
      !! -------
      !!    *  pus : u* aka friction velocity                            [m/s]
      !!    *  pnu : kinematic viscosity of air                          [lolo]
      !!
      !! OUTPUT :
      !! --------
      !!    * pz0, pz0t, pz0q: roughness lengthes for momentum, theta and q    [m]
      !!
      !! ** Author: L. Brodeau, May 2021 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(in)  :: pus
      REAL(wp), INTENT(in)  :: pnu
      !!
      REAL(wp), INTENT(out) :: pz0
      REAL(wp), INTENT(out) :: pz0t
      REAL(wp), INTENT(out) :: pz0q

      !
      REAL(wp) :: ztmp

      ztmp  = pnu/pus

      pz0   = MIN( ABS( alpha_M*ztmp + charn0_ecmwf*pus*pus/grav ) , 0.001_wp)
      pz0t  = MIN( ABS( alpha_H*ztmp                             ) , 0.001_wp)   ! eq.3.26, Chap.3, p.34, IFS doc - Cy31r1
      pz0q  = MIN( ABS( alpha_Q*ztmp                             ) , 0.001_wp)

   END SUBROUTINE update_z0_ecmwf


END MODULE mod_blk_ecmwf_ij
