! AeroBulk / 2019 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_blk_coare3p0
   !!====================================================================================
   !!       Computes turbulent components of surface fluxes
   !!         according to Fairall et al. 2018 (COARE v3.6)
   !!         "THE TOGA-COARE BULK AIR-SEA FLUX ALGORITHM"
   !!
   !!       With Cool-Skin and Warm-Layer correction of SST (if needed)
   !!
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (usually 2m) to zu (usually 10m) if needed
   !!   * the "effective" bulk wind speed at zu: Ubzu (including gustiness contribution in unstable conditions)
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!       Routine turb_coare3p0 maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, July 2019
   !!
   !!====================================================================================
   USE mod_const        !: physical and othe constants
   USE mod_phymbl       !: thermodynamics
   USE mod_common_coare !: common stuff across COARE versions...
   USE mod_skin_coare   !: cool-skin parameterization

   IMPLICIT NONE

   INTERFACE charn_coare3p0
      MODULE PROCEDURE charn_coare3p0_vctr, charn_coare3p0_sclr
   END INTERFACE charn_coare3p0

   PRIVATE

   PUBLIC :: TURB_COARE3P0, charn_coare3p0

   !! COARE own values for given constants:
   REAL(wp), PARAMETER :: zi0   = 600._wp     ! scale height of the atmospheric boundary layer...
   REAL(wp), PARAMETER :: Beta0 =  1.25_wp    ! gustiness parameter
   REAL(wp), PARAMETER :: zeta_abs_max = 50._wp
   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE COARE3P0_INIT(nx, ny, l_use_wl)
      !!---------------------------------------------------------------------
      !!                  ***  FUNCTION coare3p0_init  ***
      !!
      !! INPUT :
      !! -------
      !!    * l_use_wl : use the warm-layer parameterization
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: nx, ny   ! shape of the domain
      LOGICAL, INTENT(in) :: l_use_wl ! use the warm-layer parameterization
      INTEGER :: ierr
      !!---------------------------------------------------------------------
      IF( l_use_wl ) THEN
         ierr = 0
         ALLOCATE ( Tau_ac(nx,ny) , Qnt_ac(nx,ny), dT_wl(nx,ny), Hz_wl(nx,ny), STAT=ierr )
         IF( ierr > 0 ) CALL ctl_stop( ' COARE3P0_INIT => allocation of Tau_ac, Qnt_ac, dT_wl & Hz_wl failed!' )
         Tau_ac(:,:) = 0._wp
         Qnt_ac(:,:) = 0._wp
         dT_wl(:,:)  = 0._wp
         Hz_wl(:,:)  = Hwl_max
      ENDIF
      !IF( l_use_cs ) THEN
      !   ierr = 0
      !   ALLOCATE ( dT_cs(nx,ny), STAT=ierr )
      !   IF( ierr > 0 ) CALL ctl_stop( ' COARE3P0_INIT => allocation of dT_cs failed!' )
      !   dT_cs(:,:) = -0.25_wp  ! First guess of skin correction
      !ENDIF
   END SUBROUTINE COARE3P0_INIT

   SUBROUTINE COARE3P0_EXIT(l_use_wl)
      !!---------------------------------------------------------------------
      !!                  ***  FUNCTION coare3p0_exit  ***
      !!
      !! INPUT :
      !! -------
      !!    * l_use_wl : use the warm-layer parameterization
      !!---------------------------------------------------------------------
      LOGICAL , INTENT(in) ::   l_use_wl ! use the warm-layer parameterization
      INTEGER :: ierr
      !!---------------------------------------------------------------------
      IF( l_use_wl ) THEN
         ierr = 0
         DEALLOCATE ( Tau_ac , Qnt_ac, dT_wl, Hz_wl, STAT=ierr )
         IF( ierr > 0 ) CALL ctl_stop( ' COARE3P0_EXIT => deallocation of Tau_ac, Qnt_ac, dT_wl & Hz_wl failed!' )
      ENDIF
      !IF( l_use_cs ) THEN
      !   ierr = 0
      !   DEALLOCATE ( dT_cs, STAT=ierr )
      !   IF( ierr > 0 ) CALL ctl_stop( ' COARE3P0_EXIT => deallocation of dT_cs failed!' )
      !ENDIF
   END SUBROUTINE COARE3P0_EXIT




   SUBROUTINE TURB_COARE3P0( kt, zt, zu, T_s, t_zt, q_s, q_zt, U_zu, l_use_cs, l_use_wl, &
      &                      Cd, Ch, Ce, t_zu, q_zu, Ubzu,                               &
      &                      Qsw, rad_lw, slp, pdT_cs,                                   & ! optionals for cool-skin (and warm-layer)
      &                      isecday_utc, plong, pdT_wl, pHz_wl,                         & ! optionals for warm-layer only
      &                      CdN, ChN, CeN, xz0, xu_star, xL, xUN10 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_coare3p0  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to Fairall et al. (2003)
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
      !!    * isecday_utc:
      !!    *  plong  : longitude array                                       [deg.E]
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
      INTEGER,  INTENT(in   )                 ::   kt       ! current time step
      REAL(wp), INTENT(in   )                 ::   zt       ! height for t_zt and q_zt                    [m]
      REAL(wp), INTENT(in   )                 ::   zu       ! height for U_zu                             [m]
      REAL(wp), INTENT(inout), DIMENSION(:,:) ::   T_s      ! sea surface temperature                [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   t_zt     ! potential air temperature              [Kelvin]
      REAL(wp), INTENT(inout), DIMENSION(:,:) ::   q_s      ! sea surface specific humidity           [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   q_zt     ! specific air humidity at zt             [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   U_zu     ! relative wind module at zu                [m/s]
      LOGICAL , INTENT(in   )                 ::   l_use_cs ! use the cool-skin parameterization
      LOGICAL , INTENT(in   )                 ::   l_use_wl ! use the warm-layer parameterization
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   Cd       ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   Ch       ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   Ce       ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   t_zu     ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   q_zu     ! spec. humidity adjusted at zu           [kg/kg]
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   Ubzu    ! bulk wind speed at zu                     [m/s]
      !
      REAL(wp), INTENT(in   ), OPTIONAL, DIMENSION(:,:) ::   Qsw      !             [W/m^2]
      REAL(wp), INTENT(in   ), OPTIONAL, DIMENSION(:,:) ::   rad_lw   !             [W/m^2]
      REAL(wp), INTENT(in   ), OPTIONAL, DIMENSION(:,:) ::   slp      !             [Pa]
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   pdT_cs
      !
      INTEGER,  INTENT(in   ), OPTIONAL                     ::   isecday_utc ! current UTC time, counted in second since 00h of the current day
      REAL(wp), INTENT(in   ), OPTIONAL, DIMENSION(:,:) ::   plong    !             [deg.E]
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   pdT_wl   !             [K]
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   pHz_wl   !             [m]
      !
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   CdN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   ChN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   CeN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   xz0  ! Aerodynamic roughness length   [m]
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   xu_star  ! u*, friction velocity
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   xL  ! zeta (zu/L)
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   xUN10  ! Neutral wind at zu
      !
      INTEGER :: Ni, Nj, ji, jj, jit
      LOGICAL :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: xSST     ! to back up the initial bulk SST
      !
      REAL(wp) :: zdt, zdq, zus, zus2, zUzu, zts, zqs, zNu_a, z1oL, zdT_cs
      REAL(wp) :: zUn10, zgust2, zz0, zz0t, zzta_u, zzta_t
      REAL(wp) :: zlog_10, zlog_zu, zlog_zt, zlog_ztu, zlog_z0, zlog_z0t
      REAL(wp) :: zQns, zQlat, zTau
      REAL(wp) :: ztmp0, ztmp1
      !
      LOGICAL ::  lreturn_cdn=.FALSE., lreturn_chn=.FALSE., lreturn_cen=.FALSE., &
         &        lreturn_z0=.FALSE., lreturn_ustar=.FALSE., lreturn_L=.FALSE., lreturn_UN10=.FALSE.
      CHARACTER(len=40), PARAMETER :: crtnm = 'turb_coare3p0@mod_blk_coare3p0.f90'
      !!----------------------------------------------------------------------------------
      Ni = SIZE(T_s,1)
      Nj = SIZE(T_s,2)

      IF( kt == nit000 ) CALL COARE3P0_INIT( Ni, Nj,  l_use_wl )

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

      IF( l_use_wl .AND. (.NOT.(PRESENT(Qsw) .AND. PRESENT(rad_lw) .AND. PRESENT(slp)     &
         &  .AND. PRESENT(isecday_utc) .AND. PRESENT(plong)) ) ) &
         &   CALL ctl_stop( '['//TRIM(crtnm)//'] => ' , 'you need to provide Qsw, rad_lw, slp, isecday_utc', &
         &   ' & plong to use warm-layer param!'  )

      IF( l_use_cs .OR. l_use_wl ) THEN
         ALLOCATE ( xSST(Ni,Nj) )
         xSST = T_s ! backing up the bulk SST
         IF( l_use_cs ) T_s = T_s - 0.25_wp   ! First guess of correction
         q_s    = rdct_qsat_salt*q_sat(MAX(T_s, 200._wp), slp) ! First guess of q_s
      ENDIF

      !! Constants:
      zlog_10  = LOG(10._wp)
      zlog_zt  = LOG(zt)
      zlog_zu  = LOG(zu)
      zlog_ztu = LOG(zt/zu)

      DO jj = 1, Nj
         DO ji = 1, Ni

            zUzu = U_zu(ji,jj)

            CALL FIRST_GUESS_COARE( zt, zu, T_s(ji,jj), t_zt(ji,jj), q_s(ji,jj), q_zt(ji,jj), zUzu, &
               &                    charn_coare3p0(zUzu),  zus, zts, zqs, &
               &                    t_zu(ji,jj), q_zu(ji,jj), Ubzu(ji,jj),  pz0=zz0 )

            zlog_z0 = LOG(zz0)
            znu_a   = visc_air(t_zu(ji,jj)) ! Air viscosity (m^2/s) at zt given from temperature in (K)

            !! Pot. temp. difference (and we don't want it to be 0!)
            zdt = t_zu(ji,jj) - T_s(ji,jj) ;   zdt = SIGN( MAX(ABS(zdt),1.E-6_wp), zdt )
            zdq = q_zu(ji,jj) - q_s(ji,jj) ;   zdq = SIGN( MAX(ABS(zdq),1.E-9_wp), zdq )


            !! ITERATION BLOCK
            DO jit = 1, nb_iter

               zus2    = zus*zus   ! u*^2

               !!Inverse of Obukov length (1/L) :
               z1oL = One_on_L(t_zu(ji,jj), q_zu(ji,jj), zus, zts, zqs)  ! 1/L == 1/[Obukhov length]
               z1oL = SIGN( MIN(ABS(z1oL),200._wp), z1oL ) ! 1/L (prevents FPE from stupid values from masked region later on...)

               !! Update wind at zu with convection-related wind gustiness in unstable conditions (Fairall et al. 2003, Eq.8):
               zgust2 = Beta0*Beta0*zus2*(MAX(-zi0*z1oL/vkarmn,0._wp))**(2._wp/3._wp) ! square of wind gustiness contribution, zgust2 == Ug^2
               !!   ! Only true when unstable (L<0) => when z1oL < 0 => explains "-" before zi0
               Ubzu(ji,jj) = MAX(SQRT(zUzu*zUzu + zgust2), 0.2_wp)        ! include gustiness in bulk wind speed
               ! => 0.2 prevents Ubzu to be 0 in stable case when U_zu=0.

               !! Stability parameters:
               zzta_u = zu*z1oL
               zzta_u = SIGN( MIN(ABS(zzta_u),zeta_abs_max), zzta_u )
               IF( .NOT. l_zt_equal_zu ) THEN
                  zzta_t = zt*z1oL
                  zzta_t = SIGN( MIN(ABS(zzta_t),zeta_abs_max), zzta_t )
               ENDIF

               !! Adjustment the wind at 10m (not needed in the current algo form):
               !IF( zu \= 10._wp ) U10 = U_zu + zus/vkarmn*(LOG(10._wp/zu) - psi_m_coare(10._wp*z1oL) + psi_m_coare(zzta_u))

               !! Roughness lengthes z0, z0t (z0q = z0t) :
               zUn10 = zus/vkarmn*(zlog_10 - zlog_z0)       ! Neutral wind speed at 10m
               zz0    = charn_coare3p0(zUn10)*zus2/grav + 0.11_wp*znu_a/zus ! Roughness length (eq.6)
               zz0     = MIN( MAX(ABS(zz0), 1.E-9) , 1._wp )  ! (prevents FPE from stupid values from masked region later on)
               zlog_z0 = LOG(zz0)

               ztmp1 = ( znu_a / (zz0*zus) )**0.6_wp         ! (1./Re_r)^0.6 (Re_r: roughness Reynolds number) COARE 3.0 - specific!
               zz0t   = MIN( 1.1E-4_wp , 5.5E-5_wp*ztmp1 )   ! Scalar roughness for temp. and q (eq.28) #LB: some use 1.15 not 1.1 !!!
               zz0t   = MIN( MAX(ABS(zz0t), 1.E-9) , 1._wp ) ! (prevents FPE from stupid values from masked region later on)
               zlog_z0t = LOG(zz0t)

               !! Turbulent scales at zu :
               ztmp0   = psi_h_coare(zzta_u)
               ztmp1   = vkarmn/(zlog_zu - zlog_z0t - ztmp0) ! #LB: in ztmp0, some use psi_h_coare(zzta_t) rather than psi_h_coare(zzta_t) ???

               zts = zdt*ztmp1
               zqs = zdq*ztmp1
               zus = MAX( Ubzu(ji,jj)*vkarmn/(zlog_zu - zlog_z0 - psi_m_coare(zzta_u)) , 1.E-9 )

               IF( .NOT. l_zt_equal_zu ) THEN
                  !! Re-updating temperature and humidity at zu if zt /= zu :
                  ztmp1 = zlog_zt - zlog_zu + ztmp0 - psi_h_coare(zzta_t)
                  t_zu(ji,jj) = t_zt(ji,jj) - zts/vkarmn*ztmp1
                  q_zu(ji,jj) = q_zt(ji,jj) - zqs/vkarmn*ztmp1
               ENDIF

               IF( l_use_cs ) THEN
                  !! Cool-skin contribution
                  CALL UPDATE_QNSOL_TAU( zu, T_s(ji,jj), q_s(ji,jj), t_zu(ji,jj), q_zu(ji,jj), zus, zts, zqs, &
                     &                   zUzu, Ubzu(ji,jj), slp(ji,jj), rad_lw(ji,jj), zQns, zTau, Qlat=zQlat )

                  CALL CS_COARE( Qsw(ji,jj), zQns, zus, xSST(ji,jj), zQlat, zdT_cs )
                  IF( PRESENT(pdT_cs) ) pdT_cs(ji,jj) = zdT_cs
                  T_s(ji,jj) = xSST(ji,jj) + zdT_cs
                  IF( l_use_wl ) T_s(ji,jj) = T_s(ji,jj) + dT_wl(ji,jj)
                  q_s(ji,jj) = rdct_qsat_salt*q_sat(MAX(T_s(ji,jj), 200._wp), slp(ji,jj))
               ENDIF

               IF( l_use_wl ) THEN
                  !! Warm-layer contribution
                  CALL UPDATE_QNSOL_TAU( zu, T_s(ji,jj), q_s(ji,jj), t_zu(ji,jj), q_zu(ji,jj), zus, zts, zqs, &
                     &                   zUzu, Ubzu(ji,jj), slp(ji,jj), rad_lw(ji,jj),  zQns, zTau)
                  !! In WL_COARE or , Tau_ac and Qnt_ac must be updated at the final itteration step => add a flag to do this!
                  CALL WL_COARE( ji, jj, Qsw(ji,jj), zQns, zTau, xSST(ji,jj), plong(ji,jj), isecday_utc, MOD(nb_iter,jit) )

                  !! Updating T_s and q_s !!!
                  T_s(ji,jj) = xSST(ji,jj) + dT_wl(ji,jj)
                  IF( l_use_cs ) T_s(ji,jj) = T_s(ji,jj) + zdT_cs
                  q_s(ji,jj) = rdct_qsat_salt*q_sat(MAX(T_s(ji,jj), 200._wp), slp(ji,jj))
               ENDIF

               IF( l_use_cs .OR. l_use_wl .OR. (.NOT. l_zt_equal_zu) ) THEN
                  zdt = t_zu(ji,jj) - T_s(ji,jj) ;  zdt = SIGN( MAX(ABS(zdt),1.E-6_wp), zdt )
                  zdq = q_zu(ji,jj) - q_s(ji,jj) ;  zdq = SIGN( MAX(ABS(zdq),1.E-9_wp), zdq )
               ENDIF

            END DO !DO jit = 1, nb_iter

            ! compute transfer coefficients at zu :
            ztmp0 = zus/Ubzu(ji,jj)
            Cd(ji,jj) = MAX( ztmp0*ztmp0   , Cx_min )
            Ch(ji,jj) = MAX( ztmp0*zts/zdt , Cx_min )
            Ce(ji,jj) = MAX( ztmp0*zqs/zdq , Cx_min )

            !! Optional output
            IF( lreturn_cdn .OR. lreturn_chn .OR. lreturn_cen ) ztmp0 = 1._wp/(zlog_zu - zlog_z0)
            IF( lreturn_cdn )   CdN(ji,jj) = MAX( vkarmn2*ztmp0*ztmp0 , Cx_min )
            IF( lreturn_chn .OR. lreturn_cen ) ztmp1 = vkarmn2*ztmp0/(zlog_zu - zlog_z0t)
            IF( lreturn_chn )   ChN(ji,jj) = MAX( ztmp1 , Cx_min )
            IF( lreturn_cen )   CeN(ji,jj) = MAX( ztmp1 , Cx_min )

            IF( lreturn_z0 )        xz0(ji,jj) = zz0
            IF( lreturn_ustar ) xu_star(ji,jj) = zus
            IF( lreturn_L )          xL(ji,jj) = 1._wp / z1oL
            IF( lreturn_UN10 )    xUN10(ji,jj) = zus/vkarmn*(zlog_10 - zlog_z0)

         END DO
      END DO

      IF( l_use_wl .AND. PRESENT(pdT_wl) ) pdT_wl = dT_wl
      IF( l_use_wl .AND. PRESENT(pHz_wl) ) pHz_wl = Hz_wl

      IF( l_use_cs .OR. l_use_wl ) DEALLOCATE ( xSST )

      IF( kt == nitend ) CALL COARE3P0_EXIT( l_use_wl )

   END SUBROUTINE TURB_COARE3P0


   !!===============================================================================================
   FUNCTION charn_coare3p0_sclr( pwnd )
      !!-------------------------------------------------------------------
      !! Compute the Charnock parameter as a function of the wind speed
      !!
      !! (Fairall et al., 2003 p.577-578)
      !!
      !! Wind below 10 m/s :  alfa = 0.011
      !! Wind between 10 and 18 m/s : linear increase from 0.011 to 0.018
      !! Wind greater than 18 m/s :  alfa = 0.018
      !!
      !! Author: L. Brodeau, June 2016 / AeroBulk  (https://github.com/brodeau/aerobulk/)
      !!-------------------------------------------------------------------
      REAL(wp), INTENT(in) :: pwnd   ! wind speed
      REAL(wp) :: charn_coare3p0_sclr
      !
      REAL(wp) :: zw, zgt10, zgt18
      !!-------------------------------------------------------------------
      zw = pwnd   ! wind speed
      !!
      !! Charnock's constant, increases with the wind :
      zgt10 = 0.5_wp + SIGN(0.5_wp,(zw - 10._wp))  ! If zw<10. --> 0, else --> 1
      zgt18 = 0.5_wp + SIGN(0.5_wp,(zw - 18._wp))  ! If zw<18. --> 0, else --> 1
      !
      charn_coare3p0_sclr =  (1. - zgt10)*0.011    &    ! wind is lower than 10 m/s
         &                  + zgt10*((1. - zgt18)*(0.011 + (0.018 - 0.011) &
         &                   *(zw - 10.)/(18. - 10.)) + zgt18*( 0.018 ) )    ! Hare et al. (1999)
      !!
   END FUNCTION charn_coare3p0_sclr

   FUNCTION charn_coare3p0_vctr( pwnd )
      REAL(wp), DIMENSION(:,:), INTENT(in)           :: pwnd   ! wind speed
      REAL(wp), DIMENSION(SIZE(pwnd,1),SIZE(pwnd,2)) :: charn_coare3p0_vctr
      INTEGER  :: ji, jj
      DO jj = 1, SIZE(pwnd,2)
         DO ji = 1, SIZE(pwnd,1)
            charn_coare3p0_vctr(ji,jj) = charn_coare3p0_sclr( pwnd(ji,jj) )
         END DO
      END DO
   END FUNCTION charn_coare3p0_vctr

   !!======================================================================
END MODULE mod_blk_coare3p0
