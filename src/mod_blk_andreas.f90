!!! TO DO: consistent psi_m and psi_h needed!!! For now is those of NCAR !!!
!!
! AeroBulk / 2020 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_blk_andreas
   !!====================================================================================
   !!       Computes turbulent components of surface fluxes
   !!         according to Andreas et al. (2015)
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !!       Andreas, E.L., Mahrt, L. and Vickers, D. (2015),
   !!       An improved bulk air–sea surface flux algorithm,
   !!       including spray‐mediated transfer.
   !!       Q.J.R. Meteorol. Soc., 141: 642-654. doi:10.1002/qj.2424
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !!
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (2m) to zu (10m) if needed
   !!   * the effective bulk wind speed at 10m Ub
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!    Using the bulk formulation/param. of Large & Yeager 2008
   !!
   !!       Routine turb_andreas maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, 2020
   !!
   !!====================================================================================
   USE mod_const       !: physical and othe constants
   USE mod_phymbl      !: thermodynamics

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: TURB_ANDREAS

   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE turb_andreas( zt, zu, sst, t_zt, ssq, q_zt, U_zu,   &
      &                     Cd, Ch, Ce, t_zu, q_zu, Ub,           &
      &                    CdN, ChN, CeN, xz0, xu_star, xL, xUN10 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_andreas  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to Large & Yeager (2004) and Large & Yeager (2008)
      !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!                Returns the effective bulk wind speed at zu to be used in the bulk formulas
      !!
      !! INPUT :
      !! -------
      !!    *  zt   : height for temperature and spec. hum. of air            [m]
      !!    *  zu   : height for wind speed (usually 10m)                     [m]
      !!    *  sst  : bulk SST                                                [K]
      !!    *  t_zt : potential air temperature at zt                         [K]
      !!    *  ssq  : specific humidity at saturation at SST                  [kg/kg]
      !!    *  q_zt : specific humidity of air at zt                          [kg/kg]
      !!    *  U_zu : scalar wind speed at zu                                 [m/s]
      !!    *  SLP  : sea level pressure (needed if zt /= zu)                 [Pa]
      !!    *  gamma: adiabatic lapse-rate of moist air (needed if zt /= zu)  [K/m]
      !!
      !! OUTPUT :
      !! --------
      !!    *  Cd     : drag coefficient
      !!    *  Ch     : sensible heat coefficient
      !!    *  Ce     : evaporation coefficient
      !!    *  t_zu   : pot. air temperature adjusted at wind height zu       [K]
      !!    *  q_zu   : specific humidity of air        //                    [kg/kg]
      !!    *  Ub  : bulk wind speed at zu                                 [m/s]
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
      REAL(wp), INTENT(in   )                     ::   zt       ! height for t_zt and q_zt                    [m]
      REAL(wp), INTENT(in   )                     ::   zu       ! height for U_zu                             [m]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   sst      ! sea surface temperature                [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   t_zt     ! potential air temperature              [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   ssq      ! sea surface specific humidity           [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   q_zt     ! specific air humidity at zt             [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   U_zu     ! relative wind module at zu                [m/s]
      REAL(wp), INTENT( out), DIMENSION(jpi,jpj) ::   Cd       ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ch       ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ce       ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   t_zu     ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   q_zu     ! spec. humidity adjusted at zu           [kg/kg]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ub    ! bulk wind speed at zu                     [m/s]
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
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   u_star, t_star, q_star
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   z0       ! roughness length (momentum) [m]
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   zeta_u        ! stability parameter at height zu
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   ztmp0, ztmp1, ztmp2
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   sqrtCd       ! square root of Cd
      !
      LOGICAL ::  lreturn_cdn=.FALSE., lreturn_chn=.FALSE., lreturn_cen=.FALSE., &
         &        lreturn_z0=.FALSE., lreturn_ustar=.FALSE., lreturn_L=.FALSE., lreturn_UN10=.FALSE.
      CHARACTER(len=40), PARAMETER :: crtnm = 'turb_andreas@mod_blk_andreas.f90'
      !!----------------------------------------------------------------------------------

      ALLOCATE( u_star(jpi,jpj), t_star(jpi,jpj), q_star(jpi,jpj), &
         &          z0(jpi,jpj), zeta_u(jpi,jpj), sqrtCd(jpi,jpj), &
         &       ztmp0(jpi,jpj),  ztmp1(jpi,jpj),  ztmp2(jpi,jpj)  )

      IF( PRESENT(CdN) )     lreturn_cdn   = .TRUE.
      IF( PRESENT(ChN) )     lreturn_chn   = .TRUE.
      IF( PRESENT(CeN) )     lreturn_cen   = .TRUE.
      IF( PRESENT(xz0) )     lreturn_z0    = .TRUE.
      IF( PRESENT(xu_star) ) lreturn_ustar = .TRUE.
      IF( PRESENT(xL) )      lreturn_L     = .TRUE.
      IF( PRESENT(xUN10) )   lreturn_UN10  = .TRUE.

      l_zt_equal_zu = .FALSE.
      IF( ABS(zu - zt) < 0.01_wp )   l_zt_equal_zu = .TRUE.    ! testing "zu == zt" is risky with double precision

      Ub = MAX( 0.25_wp , U_zu ) !  relative bulk wind speed at zu

      !! First guess:
      Cd = 1.2E-3_wp
      Ch = 1.2E-3_wp
      Ce = 1.2E-3_wp

      !! Initializing values at z_u with z_t values:
      t_zu = t_zt
      q_zu = q_zt

      !! First guess of turbulent scales:
      ztmp0  = SQRT(Cd)
      u_star = ztmp0*Ub       ! u*
      t_star = Ch/ztmp0*(t_zu - sst) ! theta*
      q_star = Ce/ztmp0*(q_zu - ssq) ! q*



      !! ITERATION BLOCK
      DO j_itt = 1, nb_itt

         !! Stability parameter :
         zeta_u   = zu*One_on_L( t_zu, q_zu, u_star, t_star, q_star )
         zeta_u   = sign( min(abs(zeta_u),10._wp), zeta_u )

         ztmp0 = UN10_from_ustar( zu, Ub, u_star, psi_m(zeta_u) ) ! UN10

         u_star = U_STAR_ANDREAS( ztmp0 )

         !! Drag coefficient:
         ztmp0 = u_star/Ub
         Cd    = ztmp0*ztmp0

         !! Roughness length:
         z0    = z0_from_Cd( zu, Cd,  ppsi=psi_m(zeta_u) )

         !! z0t and z0q, based on LKB, just like into COARE 2.5:
         ztmp0 = z0 * u_star / visc_air(t_zu) ! Re_r
         ztmp1 = z0tq_LKB( 1, ztmp0, z0 )    ! z0t
         ztmp2 = z0tq_LKB( 2, ztmp0, z0 )    ! z0q

         !! Turbulent scales at zu :
         ztmp0   = psi_h(zeta_u)   ! lolo: zeta_u for scalars???
         t_star  = (t_zu - sst)*vkarmn/(LOG(zu) - LOG(ztmp1) - ztmp0)  ! theta* (ztmp1 == z0t in rhs term)
         q_star  = (q_zu - ssq)*vkarmn/(LOG(zu) - LOG(ztmp2) - ztmp0)  !   q*   (ztmp2 == z0q in rhs term)

         IF( .NOT. l_zt_equal_zu ) THEN
            !! Re-updating temperature and humidity at zu if zt /= zu:
            ztmp0 = zt*One_on_L( t_zu, q_zu, u_star, t_star, q_star ) ! zeta_t
            ztmp0 = LOG(zt/zu) + psi_h(zeta_u) - psi_h(ztmp0)
            t_zu = t_zt - t_star/vkarmn*ztmp0
            q_zu = q_zt - q_star/vkarmn*ztmp0
         ENDIF

      END DO !DO j_itt = 1, nb_itt

      
      ! compute transfer coefficients at zu:
      ztmp0 = u_star/Ub
      Cd   = ztmp0*ztmp0

      ztmp1 = t_zu - sst ;  ztmp1 = SIGN( MAX(ABS(ztmp1),1.E-6_wp), ztmp1 )  ! dt_zu
      ztmp2 = q_zu - ssq ;  ztmp2 = SIGN( MAX(ABS(ztmp2),1.E-9_wp), ztmp2 )  ! dq_zu
      Ch   = ztmp0*t_star/ztmp1
      Ce   = ztmp0*q_star/ztmp2
      

      IF( lreturn_cdn )   CdN     = -999.
      IF( lreturn_chn )   ChN     = -999.
      IF( lreturn_cen )   CeN     = -999.
      !IF( lreturn_z0 )    xz0     = z0_from_Cd( zu, Cd,  ppsi=psi_m(zeta_u) )
      IF( lreturn_z0 )    xz0     = z0
      IF( lreturn_ustar ) xu_star = u_star
      IF( lreturn_L )     xL      = zu/zeta_u
      IF( lreturn_UN10 )  xUN10   =  UN10_from_ustar( zu, Ub, u_star, psi_m(zeta_u) )

      DEALLOCATE( u_star, t_star, q_star, z0, zeta_u, sqrtCd, ztmp0, ztmp1, ztmp2 ) !

   END SUBROUTINE turb_andreas


   FUNCTION U_STAR_ANDREAS( pun10 )
      !!----------------------------------------------------------------------------------
      !! Estimate of the friction velocity as a function of the neutral-stability wind
      !! speed at at 10m
      !!
      !! Origin: Eq.(2.2) of Andreas et al. (2015)
      !!
      !! ** Author: L. Brodeau, April 2020 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pun10          !: neutral-stability scalar wind speed at 10m (m/s)
      REAL(wp), DIMENSION(jpi,jpj)             :: u_star_andreas !: friction velocity    [m/s]
      !
      INTEGER  ::     ji, jj ! dummy loop indices
      REAL(wp) :: za, zt, zw ! local scalars
      !!----------------------------------------------------------------------------------
      !
      DO jj = 1, jpj
         DO ji = 1, jpi

            zw  = pun10(ji,jj)

            za = zw - 8.271_wp

            zt = za + SQRT( MAX( 0.12_wp*za*za + 0.181_wp , 0._wp ) )

            u_star_andreas(ji,jj) =   0.239_wp + 0.0433_wp * zt

         END DO
      END DO
      !
   END FUNCTION U_STAR_ANDREAS


   FUNCTION psi_m( pzeta )
      !!----------------------------------------------------------------------------------
      !! Universal profile stability function for momentum
      !!    !! Psis, L&Y 2004, Eq. (8c), (8d), (8e)
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: psi_m
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) :: zx2, zx, zstab   ! local scalars
      !!----------------------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
            zx2 = SQRT( ABS( 1._wp - 16._wp*pzeta(ji,jj) ) )
            zx2 = MAX( zx2 , 1._wp )
            zx  = SQRT( zx2 )
            zstab = 0.5_wp + SIGN( 0.5_wp , pzeta(ji,jj) )
            !
            psi_m(ji,jj) =        zstab  * (-5._wp*pzeta(ji,jj))       &          ! Stable
               &          + (1._wp - zstab) * (2._wp*LOG((1._wp + zx)*0.5_wp)   &          ! Unstable
               &               + LOG((1._wp + zx2)*0.5_wp) - 2._wp*ATAN(zx) + rpi*0.5_wp)  !    "
            !
         END DO
      END DO
   END FUNCTION psi_m


   FUNCTION psi_h( pzeta )
      !!----------------------------------------------------------------------------------
      !! Universal profile stability function for temperature and humidity
      !!    !! Psis, L&Y 2004, Eq. (8c), (8d), (8e)
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: psi_h
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) :: zx2, zstab  ! local scalars
      !!----------------------------------------------------------------------------------
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            zx2 = SQRT( ABS( 1._wp - 16._wp*pzeta(ji,jj) ) )
            zx2 = MAX( zx2 , 1._wp )
            zstab = 0.5_wp + SIGN( 0.5_wp , pzeta(ji,jj) )
            !
            psi_h(ji,jj) =         zstab  * (-5._wp*pzeta(ji,jj))        &  ! Stable
               &           + (1._wp - zstab) * (2._wp*LOG( (1._wp + zx2)*0.5_wp ))   ! Unstable
            !
         END DO
      END DO
   END FUNCTION psi_h

   !!======================================================================
END MODULE mod_blk_andreas
