! AeroBulk / 2019 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!

MODULE mod_blk_coare3p0
   !!======================================================================================================================
   !!           COARE 3.0
   !!           +++++++++
   !!       Computes turbulent bulk transfer coefficients according to:
   !!           * Fairall et al., 2003
   !!
   !!======================================================================================================================
   !!
   !!
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
   USE mod_skin_coare   !: cool-skin & warm-layer parameterizations

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: TURB_COARE3P0, charn_coare3p0

   CHARACTER(len=8), PARAMETER :: clbl = 'COARE3P0'

   !! COARE own values for given constants:
   REAL(wp), PARAMETER :: zi0   = 600._wp     ! scale height of the atmospheric boundary layer...
   REAL(wp), PARAMETER :: Beta0 =  1.25_wp    ! gustiness parameter
   REAL(wp), PARAMETER :: zeta_abs_max = 50._wp
   !!----------------------------------------------------------------------
CONTAINS



   SUBROUTINE turb_coare3p0( kt, zt, zu, pT_s, pt_zt, pq_s, pq_zt, pU_zu, l_use_cs, l_use_wl, &
      &                       pCd, pCh, pCe, pt_zu, pq_zu, pUbzu,                             &
      &                       pQsw, prad_lw, pslp, pdT_cs,                                    & ! optionals for cool-skin (and warm-layer)
      &                       isecday_utc, plong,                                             & ! optionals for warm-layer only
      &                       pdT_wl, pHz_wl,                                                 & ! optionals for warm-layer only
      &                       pCdN, pChN, pCeN, pz0, pu_star, pL, pUN10 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_coare3p0  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to Fairall et al. (2003)
      !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!                Returns the effective bulk wind speed at zu to be used in the bulk formulas
      !!
      !!                Applies the cool-skin warm-layer correction of the SST to pT_s
      !!                if the net shortwave flux at the surface (pQsw), the downwelling longwave
      !!                radiative fluxes at the surface (prad_lw), and the sea-leve pressure (pslp)
      !!                are provided as (optional) arguments!
      !!
      !! INPUT :
      !! -------
      !!    *  kt    : current time step (starts at 1)
      !!    *  zt    : height for temperature and spec. hum. of air            [m]
      !!    *  zu    : height for wind speed (usually 10m)                     [m]
      !!    *  pt_zt : potential air temperature at zt                         [K]
      !!    *  pq_zt : specific humidity of air at zt                          [kg/kg]
      !!    *  pU_zu : scalar wind speed at zu                                 [m/s]
      !!    * l_use_cs : use the cool-skin parameterization
      !!    * l_use_wl : use the warm-layer parameterization
      !!
      !! INPUT/OUTPUT:
      !! -------------
      !!    *  pT_s  : always "bulk SST" as input                              [K]
      !!              -> unchanged "bulk SST" as output if CSWL not used      [K]
      !!              -> skin temperature as output if CSWL used              [K]
      !!
      !!    *  pq_s  : SSQ aka saturation specific humidity at temp. pT_s       [kg/kg]
      !!              -> doesn't need to be given a value if skin temp computed (in case l_use_cs=True or l_use_wl=True)
      !!              -> MUST be given the correct value if not computing skint temp. (in case l_use_cs=False or l_use_wl=False)
      !!
      !! OPTIONAL INPUT:
      !! ---------------
      !!    *  pQsw    : net solar flux (after albedo) at the surface (>0)     [W/m^2]
      !!    *  prad_lw : downwelling longwave radiation at the surface  (>0)   [W/m^2]
      !!    *  pslp    : sea-level pressure                                    [Pa]
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
      !!    *  pCd     : drag coefficient
      !!    *  pCh     : sensible heat coefficient
      !!    *  pCe     : evaporation coefficient
      !!    *  pt_zu   : pot. air temperature adjusted at wind height zu       [K]
      !!    *  pq_zu   : specific humidity of air        //                    [kg/kg]
      !!    *  pUbzu   : bulk wind speed at zu                                 [m/s]
      !!
      !! OPTIONAL OUTPUT:
      !! ----------------
      !!    * pCdN      : neutral-stability drag coefficient
      !!    * pChN      : neutral-stability sensible heat coefficient
      !!    * pCeN      : neutral-stability evaporation coefficient
      !!    * pz0      : return the aerodynamic roughness length (integration constant for wind stress) [m]
      !!    * pu_star  : return u* the friction velocity                    [m/s]
      !!    * pL       : return the Obukhov length                          [m]
      !!    * pUN10    : neutral wind speed at 10m                          [m/s]
      !!
      !! ** Author: L. Brodeau, June 2019 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      INTEGER,  INTENT(in   )                     ::   kt       ! current time step
      REAL(wp), INTENT(in   )                     ::   zt       ! height for pt_zt and pq_zt                [m]
      REAL(wp), INTENT(in   )                     ::   zu       ! height for pU_zu                          [m]
      REAL(wp), INTENT(inout), DIMENSION(:,:) ::   pT_s     ! sea surface temperature              [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   pt_zt    ! potential air temperature            [Kelvin]
      REAL(wp), INTENT(inout), DIMENSION(:,:) ::   pq_s     ! sea surface specific humidity         [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   pq_zt    ! specific air humidity at zt           [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   pU_zu    ! relative wind module at zu              [m/s]
      LOGICAL , INTENT(in   )                     ::   l_use_cs ! use the cool-skin parameterization
      LOGICAL , INTENT(in   )                     ::   l_use_wl ! use the warm-layer parameterization
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   pCd      ! transfer coefficient for momentum         [-]
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   pCh      ! transfer coefficient for sensible heat    [-]
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   pCe      ! transfert coefficient for evaporation     [-]
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   pt_zu    ! pot. air temp. adjusted at zu             [K]
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   pq_zu    ! spec. humidity adjusted at zu         [kg/kg]
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   pUbzu    ! bulk wind speed at zu                   [m/s]
      !
      REAL(wp), INTENT(in   ), OPTIONAL, DIMENSION(:,:) ::   pQsw      !             [W/m^2]
      REAL(wp), INTENT(in   ), OPTIONAL, DIMENSION(:,:) ::   prad_lw   !             [W/m^2]
      REAL(wp), INTENT(in   ), OPTIONAL, DIMENSION(:,:) ::   pslp      !             [Pa]
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   pdT_cs
      !
      INTEGER,  INTENT(in   ), OPTIONAL                     ::   isecday_utc ! current UTC time, counted in second since 00h of the current day
      REAL(wp), INTENT(in   ), OPTIONAL, DIMENSION(:,:) ::   plong    !             [deg.E]
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   pdT_wl   !             [K]
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   pHz_wl   !             [m]
      !
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   pCdN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   pChN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   pCeN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   pz0  ! Aerodynamic roughness length   [m]
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   pu_star  ! u*, friction velocity
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   pL  ! zeta (zu/L)
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   pUN10  ! Neutral wind at zu
      !!----------------------------------------------------------------------------------
      LOGICAL  :: l_skin
      INTEGER  :: Ni, Nj
      INTEGER  :: ji, jj, jit
      REAL(wp) :: zdum, zm_ztzu                ! => `1.` if `zu /= zt`, `0.` otherwize
      !
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: xSST     ! to back up the initial bulk SST
      !
      REAL(wp) :: zdt, zdq, zus, zus2, zUzu, zts, zqs, zNu_a, z1oL, zdT_cs
      REAL(wp) :: zUn10, zgust2, zz0, zz0t, zzta_u, zzta_t
      REAL(wp) :: zlog_10, zlog_zu, zlog_zt, zlog_ztu, zlog_z0, zlog_z0t
      REAL(wp) :: zQns, zQlat, zTau
      REAL(wp) :: zSST, zT_s, zq_s, zubzu, zt_zt, zq_zt, zt_zu, zq_zu
      REAL(wp) :: ztmp0, ztmp1
      !
      LOGICAL ::  lreturn_cdn=.FALSE., lreturn_chn=.FALSE., lreturn_cen=.FALSE., &
         &        lreturn_z0=.FALSE., lreturn_ustar=.FALSE., lreturn_L=.FALSE., lreturn_UN10=.FALSE.
      CHARACTER(len=40), PARAMETER :: crtnm = 'turb_coare3p0@mod_blk_coare3p0.f90'
      !!----------------------------------------------------------------------------------
      Ni = SIZE(pT_s,1)
      Nj = SIZE(pT_s,2)

      IF( kt == nit000 ) CALL COARE3P0_INIT( Ni, Nj,  l_use_wl )

      lreturn_cdn   =  PRESENT(pCdN)
      lreturn_chn   =  PRESENT(pChN)
      lreturn_cen   =  PRESENT(pCeN)
      lreturn_z0    =  PRESENT(pz0)
      lreturn_ustar =  PRESENT(pu_star)
      lreturn_L     =  PRESENT(pL)
      lreturn_UN10  =  PRESENT(pUN10)

      zm_ztzu = MERGE( 0._wp, 1._wp,  ABS(zu - zt) < 0.01_wp )  ! => `1.` if `zu /= zt`, `0.` otherwize

      !! Initializations for cool skin and warm layer:
      IF( l_use_cs .AND. (.NOT.(PRESENT(pQsw) .AND. PRESENT(prad_lw) .AND. PRESENT(pslp))) ) &
         &   CALL ctl_stop( '['//TRIM(crtnm)//'] => ' , 'you need to provide pQsw, prad_lw & pslp to use cool-skin param!' )

      IF( l_use_wl .AND. (.NOT.(PRESENT(pQsw) .AND. PRESENT(prad_lw) .AND. PRESENT(pslp)     &
         &  .AND. PRESENT(isecday_utc) .AND. PRESENT(plong)) ) ) &
         &   CALL ctl_stop( '['//TRIM(crtnm)//'] => ' , 'you need to provide pQsw, prad_lw, pslp, isecday_utc', &
         &   ' & plong to use warm-layer param!'  )


      l_skin = l_use_cs .OR. l_use_wl

      ALLOCATE ( xSST(Ni,Nj) )
      xSST = pT_s
      IF( l_skin ) THEN
         IF( l_use_cs ) pT_s = pT_s - 0.25_wp   ! First guess of correction
         pq_s    = rdct_qsat_salt*q_sat(MAX(pT_s, 200._wp), pslp) ! First guess of pq_s
      ENDIF

      !! Constants:
      zlog_10  = LOG(10._wp)
      zlog_zt  = LOG(zt)
      zlog_zu  = LOG(zu)
      zlog_ztu = LOG(zt/zu)

      DO jj = 1, Nj
         DO ji = 1, Ni

            zSST  = xSST(ji,jj)
            zT_s  = pT_s(ji,jj)
            zt_zt = pt_zt(ji,jj)
            zq_s  = pq_s(ji,jj)
            zq_zt = pq_zt(ji,jj)
            zUzu  = pU_zu(ji,jj)

            CALL first_guess_coare( zt, zu, zT_s, zt_zt, zq_s, zq_zt, zUzu, &
               &                    charn_coare3p0(zUzu),  zus, zts, zqs, &
               &                    zt_zu, zq_zu, zUbzu,  pz0=zz0 )

            zlog_z0 = LOG(zz0)
            znu_a   = visc_air(zt_zt) ! Air viscosity (m^2/s) at zt given from temperature in (K)

            !! Pot. temp. difference (and we don't want it to be 0!)
            zdt = zt_zu - zT_s ;   zdt = SIGN( MAX(ABS(zdt),1.E-09_wp), zdt )
            zdq = zq_zu - zq_s ;   zdq = SIGN( MAX(ABS(zdq),1.E-12_wp), zdq )


            !! ITERATION BLOCK
            DO jit = 1, nb_iter

               zus2    = zus*zus   ! u*^2

               !!Inverse of Obukov length (1/L) :
               z1oL = One_on_L(zt_zu, zq_zu, zus, zts, zqs)  ! 1/L == 1/[Obukhov length]
               z1oL = SIGN( MIN(ABS(z1oL),200._wp), z1oL ) ! 1/L (prevents FPE from stupid values from masked region later on...)

               !! Update wind at zu with convection-related wind gustiness in unstable conditions (Fairall et al. 2003, Eq.8):
               zgust2 = Beta0*Beta0*zus2*(MAX(-zi0*z1oL/vkarmn,0._wp))**(2._wp/3._wp) ! square of wind gustiness contribution, zgust2 == Ug^2
               !!   ! Only true when unstable (L<0) => when z1oL < 0 => explains "-" before zi0
               zUbzu = MAX(SQRT(zUzu*zUzu + zgust2), 0.2_wp)        ! include gustiness in bulk wind speed
               ! => 0.2 prevents pUbzu to be 0 in stable case when pU_zu=0.

               !! Stability parameters:
               zzta_u = zu*z1oL
               zzta_u = SIGN( MIN(ABS(zzta_u),zeta_abs_max), zzta_u )
               zzta_t = zt*z1oL
               zzta_t = SIGN( MIN(ABS(zzta_t),zeta_abs_max), zzta_t )


               !! Adjustment the wind at 10m (not needed in the current algo form):
               !IF( zu \= 10._wp ) U10 = pU_zu + zus/vkarmn*(LOG(10._wp/zu) - psi_m_coare(10._wp*z1oL) + psi_m_coare(zzta_u))

               !! Roughness lengthes z0, z0t (z0q = z0t) :
               zUn10 = zus/vkarmn*(zlog_10 - zlog_z0)       ! Neutral wind speed at 10m
               zz0    = charn_coare3p0(zUn10)*zus2/grav + 0.11_wp*znu_a/zus ! Roughness length (eq.6)
               zz0     = MIN( MAX(ABS(zz0), 1.E-9) , 1._wp )  ! (prevents FPE from stupid values from masked region later on)
               zlog_z0 = LOG(zz0)

               ztmp1 = ( znu_a / (zz0*zus) )**0.6_wp         ! (1./Re_r)^0.6 (Re_r: roughness Reynolds number) COARE 3.0 - specific!
               zz0t   = MIN( 1.1E-4_wp , 5.5E-5_wp*ztmp1 )   ! Scalar roughness for temp. and q (eq.28) #LB: some use 1.15 not 1.1 !!!
               zz0t   = MIN( MAX(ABS(zz0t), 1.E-9) , 1._wp ) ! (prevents FPE from stupid values from masked region later on)
               zlog_z0t = LOG(zz0t)

               !! Turbulent scales at zu:
               ztmp0   = psi_h_coare(zzta_u)
               ztmp1   = vkarmn/(zlog_zu - zlog_z0t - ztmp0) ! #LB: in ztmp0, some use psi_h_coare(zzta_t) rather than psi_h_coare(zzta_t) ???

               zts = zdt*ztmp1
               zqs = zdq*ztmp1
               zus = MAX( zUbzu*vkarmn/(zlog_zu - zlog_z0 - psi_m_coare(zzta_u)) , 1.E-9 )

               !! Adjusting temperature and humidity at zu if required by `zm_ztzu`:
               ztmp1 = zlog_zt - zlog_zu + ztmp0 - psi_h_coare(zzta_t)
               zt_zu = zt_zt - zm_ztzu*zts/vkarmn*ztmp1
               zq_zu = zq_zt - zm_ztzu*zqs/vkarmn*ztmp1

               IF( l_use_cs ) THEN
                  !! Cool-skin contribution
                  CALL UPDATE_QNSOL_TAU( zu, zT_s, zq_s, zt_zu, zq_zu, zus, zts, zqs, &
                     &                   zUzu, zUbzu, pslp(ji,jj), prad_lw(ji,jj), zQns, zTau, Qlat=zQlat )

                  CALL CS_COARE( pQsw(ji,jj), zQns, zus, zSST, zQlat, zdT_cs )
                  IF( PRESENT(pdT_cs) ) pdT_cs(ji,jj) = zdT_cs
                  zT_s = zSST + zdT_cs
                  IF( l_use_wl ) zT_s = zT_s + dT_wl(ji,jj)
                  zq_s = rdct_qsat_salt*q_sat(MAX(zT_s, 200._wp), pslp(ji,jj))
               ENDIF

               IF( l_use_wl ) THEN
                  !! Warm-layer contribution
                  CALL UPDATE_QNSOL_TAU( zu, zT_s, zq_s, zt_zu, zq_zu, zus, zts, zqs, &
                     &                   zUzu, zUbzu, pslp(ji,jj), prad_lw(ji,jj),  zQns, zTau)
                  !! In WL_COARE or , Tau_ac and Qnt_ac must be updated at the final itteration step => add a flag to do this!
                  CALL WL_COARE( ji, jj, pQsw(ji,jj), zQns, zTau, zSST, plong(ji,jj), isecday_utc, MOD(nb_iter,jit) )
                  !! Updating pT_s and pq_s !!!
                  zT_s = zSST + dT_wl(ji,jj)
                  IF( l_use_cs ) zT_s = zT_s + zdT_cs
                  zq_s = rdct_qsat_salt*q_sat(MAX(zT_s, 200._wp), pslp(ji,jj))
               ENDIF

               zdt = zt_zu - zT_s ;  zdt = SIGN( MAX(ABS(zdt),1.E-09_wp), zdt )
               zdq = zq_zu - zq_s ;  zdq = SIGN( MAX(ABS(zdq),1.E-12_wp), zdq )

            END DO !DO jit = 1, nb_iter

            !! Update arrays that are returned by the routine:
            pT_s(ji,jj)  = zT_s
            pq_s(ji,jj)  = zq_s
            pt_zu(ji,jj) = zt_zu
            pq_zu(ji,jj) = zq_zu
            pUbzu(ji,jj) = zUbzu


            ! compute transfer coefficients at zu :
            ztmp0 = zus/zUbzu
            pCd(ji,jj) = MAX( ztmp0*ztmp0   , Cx_min )
            pCh(ji,jj) = MAX( ztmp0*zts/zdt , Cx_min )
            pCe(ji,jj) = MAX( ztmp0*zqs/zdq , Cx_min )

            !! Optional output
            IF( lreturn_cdn .OR. lreturn_chn .OR. lreturn_cen ) ztmp0 = 1._wp/(zlog_zu - zlog_z0)
            IF( lreturn_cdn )   pCdN(ji,jj) = MAX( vkarmn2*ztmp0*ztmp0 , Cx_min )
            IF( lreturn_chn .OR. lreturn_cen ) ztmp1 = vkarmn2*ztmp0/(zlog_zu - zlog_z0t)
            IF( lreturn_chn )   pChN(ji,jj) = MAX( ztmp1 , Cx_min )
            IF( lreturn_cen )   pCeN(ji,jj) = MAX( ztmp1 , Cx_min )

            IF( lreturn_z0 )        pz0(ji,jj) = zz0
            IF( lreturn_ustar ) pu_star(ji,jj) = zus
            IF( lreturn_L )          pL(ji,jj) = 1._wp / z1oL
            IF( lreturn_UN10 )    pUN10(ji,jj) = zus/vkarmn*(zlog_10 - zlog_z0)

         END DO
      END DO

      IF( l_use_wl .AND. PRESENT(pdT_wl) ) pdT_wl = dT_wl
      IF( l_use_wl .AND. PRESENT(pHz_wl) ) pHz_wl = Hz_wl

      DEALLOCATE ( xSST )

      IF( kt == nitend ) CALL COARE3P0_EXIT( l_use_wl )

   END SUBROUTINE turb_coare3p0



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
         IF( ierr > 0 ) CALL ctl_stop( ' '//clbl//'_INIT => allocation of Tau_ac, Qnt_ac, dT_wl & Hz_wl failed!' )
         Tau_ac(:,:) = 0._wp
         Qnt_ac(:,:) = 0._wp
         dT_wl(:,:)  = 0._wp
         Hz_wl(:,:)  = Hwl_max
      ENDIF
      !IF( l_use_cs ) THEN
      !   ierr = 0
      !   ALLOCATE ( dT_cs(nx,ny), STAT=ierr )
      !   IF( ierr > 0 ) CALL ctl_stop( ' '//clbl//'_INIT => allocation of dT_cs failed!' )
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
         IF( ierr > 0 ) CALL ctl_stop( ' '//clbl//'_EXIT => deallocation of Tau_ac, Qnt_ac, dT_wl & Hz_wl failed!' )
      ENDIF
      !IF( l_use_cs ) THEN
      !   ierr = 0
      !   DEALLOCATE ( dT_cs, STAT=ierr )
      !   IF( ierr > 0 ) CALL ctl_stop( ' '//clbl//'_EXIT => deallocation of dT_cs failed!' )
      !ENDIF
   END SUBROUTINE COARE3P0_EXIT






   !!===============================================================================================
   FUNCTION charn_coare3p0( pwnd )
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
      REAL(wp) :: charn_coare3p0
      !
      REAL(wp) :: zw, zgt10, zgt18
      !!-------------------------------------------------------------------
      zw = pwnd   ! wind speed
      !!
      !! Charnock's constant, increases with the wind :
      zgt10 = 0.5_wp + SIGN(0.5_wp,(zw - 10._wp))  ! If zw<10. --> 0, else --> 1
      zgt18 = 0.5_wp + SIGN(0.5_wp,(zw - 18._wp))  ! If zw<18. --> 0, else --> 1
      !
      charn_coare3p0 =  (1. - zgt10)*0.011    &    ! wind is lower than 10 m/s
         &              + zgt10*((1. - zgt18)*(0.011 + (0.018 - 0.011) &
         &              *(zw - 10.)/(18. - 10.)) + zgt18*( 0.018 ) )    ! Hare et al. (1999)
      !!
   END FUNCTION charn_coare3p0

   !!======================================================================
END MODULE mod_blk_coare3p0
