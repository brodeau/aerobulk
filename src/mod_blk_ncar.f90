! AeroBulk / 2019 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_blk_ncar
   !!====================================================================================
   !!       Computes turbulent components of surface fluxes
   !!         according to Large & Yeager (2004,2008)
   !!
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (2m) to zu (10m) if needed
   !!   * the effective bulk wind speed at 10m U_blk
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!    Using the bulk formulation/param. of Large & Yeager 2008
   !!
   !!       Routine turb_ncar maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, 2016
   !!
   !!====================================================================================
   USE mod_const       !: physical and othe constants
   USE mod_phymbl      !: thermodynamics

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: TURB_NCAR, CD_N10_NCAR, CH_N10_NCAR, CE_N10_NCAR

   !!----------------------------------------------------------------------
CONTAINS

   !   SUBROUTINE turb_ncar( zt, zu, sst, t_zt, ssq, q_zt, U_zu, SLP, gamma, &
   SUBROUTINE turb_ncar( zt, zu, sst, t_zt, ssq, q_zt, U_zu, &
      &                  Cd, Ch, Ce, t_zu, q_zu, U_blk,                  &
      &                      CdN, ChN, CeN, xz0, xu_star, xL, xUN10 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_ncar  ***
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
      !!    *  U_blk  : bulk wind speed at zu                                 [m/s]
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
      !REAL(wp), INTENT(in  ), DIMENSION(jpi,jpj) ::   SLP      ! sea level pressure                         [Pa]
      !REAL(wp), INTENT(in  ), DIMENSION(jpi,jpj) ::   gamma    ! adiabatic lapse-rate of moist air         [K/m]
      REAL(wp), INTENT( out), DIMENSION(jpi,jpj) ::   Cd       ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ch       ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ce       ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   t_zu     ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   q_zu     ! spec. humidity adjusted at zu           [kg/kg]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   U_blk    ! bulk wind speed at zu                     [m/s]
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
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   Cx_n10        ! 10m neutral latent/sensible coefficient
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   sqrtCdn10   ! square root of Cd_n10
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   zeta_u        ! stability parameter at height zu
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   ztmp0, ztmp1, ztmp2
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   sqrtCd       ! square root of Cd
      !
      LOGICAL ::  lreturn_cdn=.FALSE., lreturn_chn=.FALSE., lreturn_cen=.FALSE., &
         &        lreturn_z0=.FALSE., lreturn_ustar=.FALSE., lreturn_L=.FALSE., lreturn_UN10=.FALSE.
      CHARACTER(len=40), PARAMETER :: crtnm = 'turb_ncar@mod_blk_ncar.f90'
      !!----------------------------------------------------------------------------------

      ALLOCATE( Cx_n10(jpi,jpj), sqrtCdn10(jpi,jpj), &
         &    zeta_u(jpi,jpj), sqrtCd(jpi,jpj),      &
         &    ztmp0(jpi,jpj),  ztmp1(jpi,jpj), ztmp2(jpi,jpj) )

      IF( PRESENT(CdN) )     lreturn_cdn   = .TRUE.
      IF( PRESENT(ChN) )     lreturn_chn   = .TRUE.
      IF( PRESENT(CeN) )     lreturn_cen   = .TRUE.
      IF( PRESENT(xz0) )     lreturn_z0    = .TRUE.
      IF( PRESENT(xu_star) ) lreturn_ustar = .TRUE.
      IF( PRESENT(xL) )      lreturn_L     = .TRUE.
      IF( PRESENT(xUN10) )   lreturn_UN10  = .TRUE.

      l_zt_equal_zu = .FALSE.
      IF( ABS(zu - zt) < 0.01_wp )   l_zt_equal_zu = .TRUE.    ! testing "zu == zt" is risky with double precision

      U_blk = MAX( 0.5_wp , U_zu )   !  relative wind speed at zu (normally 10m), we don't want to fall under 0.5 m/s

      ztmp0 = cd_n10_ncar( U_blk )
      sqrtCdn10 = SQRT( ztmp0 )

      !! Initializing transf. coeff. with their first guess neutral equivalents :
      Cd = ztmp0

      Ce = CE_N10_NCAR( sqrtCdn10 )

      ztmp0 = 0.5_wp + SIGN(0.5_wp, virt_temp(t_zt, q_zt) - virt_temp(sst, ssq)) ! we guess stability based on delta of virt. pot. temp.
      Ch = CH_N10_NCAR( sqrtCdn10 , ztmp0 )

      sqrtCd = sqrtCdn10

      !! Initializing values at z_u with z_t values:
      t_zu = t_zt
      q_zu = q_zt

      !! ITERATION BLOCK
      DO j_itt = 1, nb_itt
         !
         ztmp1 = t_zu - sst   ! Updating air/sea differences
         ztmp2 = q_zu - ssq

         ! Updating turbulent scales :   (L&Y 2004 Eq. (7))
         ztmp0 = sqrtCd*U_blk       ! u*
         ztmp1 = Ch/sqrtCd*ztmp1    ! theta*
         ztmp2 = Ce/sqrtCd*ztmp2    ! q*

         ! Estimate the inverse of Obukov length (1/L) at height zu:
         ztmp0 = One_on_L( t_zu, q_zu, ztmp0, ztmp1, ztmp2 )
         !ztmp0 = One_on_L( 0.5*(t_zu + sst), 0.5*(q_zu + ssq), ztmp0, ztmp1, ztmp2 ) ! using an approximation of mean
         !                                                                ! theta & q in surface layer rather than values at zu...

         !! Stability parameters :
         zeta_u   = zu*ztmp0
         zeta_u   = sign( min(abs(zeta_u),10._wp), zeta_u )

         !! Shifting temperature and humidity at zu (L&Y 2004 Eq. (9b-9c))
         IF( .NOT. l_zt_equal_zu ) THEN
            ztmp0 = zt*ztmp0 ! zeta_t !
            ztmp0 = SIGN( MIN(ABS(ztmp0),10._wp), ztmp0 )  ! Temporaty array ztmp0 == zeta_t !!!
            ztmp0 = LOG(zt/zu) + psi_h(zeta_u) - psi_h(ztmp0)                   ! ztmp0 just used as temp array again!
            t_zu = t_zt - ztmp1/vkarmn*ztmp0    ! ztmp1 is still theta*  L&Y 2004 Eq. (9b)
            !!
            q_zu = q_zt - ztmp2/vkarmn*ztmp0    ! ztmp2 is still q*      L&Y 2004 Eq. (9c)
            q_zu = MAX(0._wp, q_zu)
            !!
            !! Prevent q_zu to reach beyond saturation:
            !ztmp0 = t_zu - zu*gamma  ! ztmp0 = absolute temp. at zu (slightly colder that pot. temp. at zu)
            !q_zu = MIN( q_sat( ztmp0, SLP ), q_zu )
            !PRINT *, 'LOLO: mod_blk_ncar.f90 => SLP =', SLP
            !PRINT *, 'LOLO: mod_blk_ncar.f90 => gamma_moist =', gamma
            !PRINT *, 'LOLO: mod_blk_ncar.f90 => t_zu =', REAL(ztmp0-rt0, 4)
            !PRINT *, 'LOLO: mod_blk_ncar.f90 => q_zu_sane =', q_sat( ztmp0, SLP ) ;            PRINT *, ''
         END IF

         ! Update neutral wind speed at 10m and neutral Cd at 10m (L&Y 2004 Eq. 9a)...
         !   In very rare low-wind conditions, the old way of estimating the
         !   neutral wind speed at 10m leads to a negative value that causes the code
         !   to crash. To prevent this a threshold of 0.25m/s is imposed.
         ztmp2 = psi_m(zeta_u)
         ztmp0 = MAX( 0.25_wp , UN10_from_CD(zu, U_blk, Cd, ppsi=ztmp2) ) ! U_n10 (ztmp2 == psi_m(zeta_u))
         ztmp0 = CD_N10_NCAR(ztmp0)                                       ! Cd_n10
         sqrtCdn10 = sqrt(ztmp0)

         !! Update of transfer coefficients:
         ztmp1  = 1._wp + sqrtCdn10/vkarmn*(LOG(zu/10._wp) - ztmp2)   ! L&Y 2004 Eq. (10a) (ztmp2 == psi_m(zeta_u))
         Cd     = ztmp0 / ( ztmp1*ztmp1 )
         sqrtCd = SQRT( Cd )

         ztmp0  = ( LOG(zu/10._wp) - psi_h(zeta_u) ) / vkarmn / sqrtCdn10
         ztmp2  = sqrtCd / sqrtCdn10

         ztmp1  = 0.5_wp + sign(0.5_wp,zeta_u)       ! stability flag
         Cx_n10 = CH_N10_NCAR( sqrtCdn10 , ztmp1 )
         ztmp1  = 1._wp + Cx_n10*ztmp0
         Ch     = Cx_n10*ztmp2 / ztmp1   ! L&Y 2004 Eq. (10b)

         Cx_n10 = CE_N10_NCAR( sqrtCdn10 )
         ztmp1  = 1._wp + Cx_n10*ztmp0
         Ce     = Cx_n10*ztmp2 / ztmp1  ! L&Y 2004 Eq. (10c)

      END DO !DO j_itt = 1, nb_itt

      IF( lreturn_cdn )   CdN     = sqrtCdn10*sqrtCdn10
      IF( lreturn_chn )   ChN     = CH_N10_NCAR( sqrtCdn10 , 0.5_wp+sign(0.5_wp,zeta_u) )
      IF( lreturn_cen )   CeN     = CE_N10_NCAR( sqrtCdn10 )
      IF( lreturn_z0 )    xz0     = z0_from_Cd( zu, sqrtCdn10*sqrtCdn10 )
      IF( lreturn_ustar ) xu_star = SQRT( Cd )*U_blk
      IF( lreturn_L )     xL      = zu/zeta_u
      IF( lreturn_UN10 )  xUN10   = UN10_from_CD( zu, U_blk, Cd, ppsi=psi_m(zeta_u) )

      DEALLOCATE( Cx_n10, sqrtCdn10, zeta_u, sqrtCd, ztmp0, ztmp1, ztmp2 ) !

   END SUBROUTINE turb_ncar


   FUNCTION CD_N10_NCAR( pw10 )
      !!----------------------------------------------------------------------------------
      !! Estimate of the neutral drag coefficient at 10m as a function
      !! of neutral wind  speed at 10m
      !!
      !! Origin: Large & Yeager 2008, Eq. (11)
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pw10           ! scalar wind speed at 10m (m/s)
      REAL(wp), DIMENSION(jpi,jpj)             :: cd_n10_ncar
      !
      INTEGER  ::     ji, jj     ! dummy loop indices
      REAL(wp) :: zgt33, zw, zw6 ! local scalars
      !!----------------------------------------------------------------------------------
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zw  = pw10(ji,jj)
            zw6 = zw*zw*zw
            zw6 = zw6*zw6
            !
            ! When wind speed > 33 m/s => Cyclone conditions => special treatment
            zgt33 = 0.5_wp + SIGN( 0.5_wp, (zw - 33._wp) )   ! If pw10 < 33. => 0, else => 1
            !
            cd_n10_ncar(ji,jj) = 1.e-3_wp * ( &
               &       (1._wp - zgt33)*( 2.7_wp/zw + 0.142_wp + zw/13.09_wp - 3.14807E-10_wp*zw6) & ! wind <  33 m/s
               &      +    zgt33   *      2.34_wp )                                                 ! wind >= 33 m/s
            !
            cd_n10_ncar(ji,jj) = MAX(cd_n10_ncar(ji,jj), 1.E-6_wp)
            !
         END DO
      END DO
      !
   END FUNCTION CD_N10_NCAR



   FUNCTION CH_N10_NCAR( psqrtcdn10 , pstab )
      !!----------------------------------------------------------------------------------
      !! Estimate of the neutral heat transfer coefficient at 10m      !!
      !! Origin: Large & Yeager 2008, Eq. (9) and (12)

      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: ch_n10_ncar
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: psqrtcdn10 ! sqrt( CdN10 )
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pstab      ! stable ABL => 1 / unstable ABL => 0
      !!----------------------------------------------------------------------------------
      IF ( ANY(pstab < -0.00001) .OR. ANY(pstab >  1.00001) ) THEN
         PRINT *, 'ERROR: CH_N10_NCAR@mod_blk_ncar.f90: pstab ='
         PRINT *, pstab
         STOP
      END IF
      !
      ch_n10_ncar = 1.e-3_wp * psqrtcdn10*( 18._wp*pstab + 32.7_wp*(1._wp - pstab) )   ! Eq. (9) & (12) Large & Yeager, 2008
      !
   END FUNCTION CH_N10_NCAR

   FUNCTION CE_N10_NCAR( psqrtcdn10 )
      !!----------------------------------------------------------------------------------
      !! Estimate of the neutral heat transfer coefficient at 10m      !!
      !! Origin: Large & Yeager 2008, Eq. (9) and (13)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: ce_n10_ncar
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: psqrtcdn10 ! sqrt( CdN10 )
      !!----------------------------------------------------------------------------------
      ce_n10_ncar = 1.e-3_wp * ( 34.6_wp * psqrtcdn10 )
      !
   END FUNCTION CE_N10_NCAR







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
END MODULE mod_blk_ncar
