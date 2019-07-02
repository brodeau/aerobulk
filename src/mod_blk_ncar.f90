! AeroBulk / 2016 / L. Brodeau
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
   !!    Using the bulk formulation/param. of "CORE2" aka "NCAR", Large & Yeager (2004,2008)
   !!
   !!       Routine turb_ncar maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, 2016
   !!
   !!====================================================================================
   USE mod_const   !: physical and othe constants
   USE mod_phymbl  !: thermodynamics

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: TURB_NCAR

CONTAINS

   SUBROUTINE turb_ncar( zt, zu, sst, t_zt, ssq, q_zt, U_zu, &
      &                  Cd, Ch, Ce, t_zu, q_zu, U_blk,      &
      &                   xz0, xu_star, xL, xUN10 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_ncar  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to Large & Yeager (2004) and Large & Yeager (2008)
      !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!                Returns the effective bulk wind speed at 10m to be used in the bulk formulas
      !!
      !! ** Method : Monin Obukhov Similarity Theory
      !!             + Large & Yeager (2004,2008) closure: CD_n10 = f(U_n10)
      !!
      !! ** References :   Large & Yeager, 2004 / Large & Yeager, 2008
      !!
      !! ** Last update: Laurent Brodeau, June 2014:
      !!    - handles both cases zt=zu and zt/=zu
      !!    - optimized: less 2D arrays allocated and less operations
      !!    - better first guess of stability by checking air-sea difference of virtual temperature
      !!       rather than temperature difference only...
      !!    - added function "cd_neutral_10m" that uses the improved parametrization of
      !!      Large & Yeager 2008. Drag-coefficient reduction for Cyclone conditions!
      !!    - using code-wide physical constants defined into "phycst.mod" rather than redifining them
      !!      => 'vkarmn' and 'grav'
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!
      !! INPUT :
      !! -------
      !!    *  zt   : height for temperature and spec. hum. of air            [m]
      !!    *  zu   : height for wind speed (generally 10m)                   [m]
      !!    *  U_zu : scalar wind speed at 10m                                [m/s]
      !!    *  sst  : SST                                                     [K]
      !!    *  t_zt : potential air temperature at zt                         [K]
      !!    *  ssq  : specific humidity at saturation at SST                  [kg/kg]
      !!    *  q_zt : specific humidity of air at zt                          [kg/kg]
      !!
      !!
      !! OUTPUT :
      !! --------
      !!    *  Cd     : drag coefficient
      !!    *  Ch     : sensible heat coefficient
      !!    *  Ce     : evaporation coefficient
      !!    *  t_zu   : pot. air temperature adjusted at wind height zu       [K]
      !!    *  q_zu   : specific humidity of air        //                    [kg/kg]
      !!    *  U_blk  : bulk wind speed at 10m                                [m/s]
      !!
      !! OPTIONAL OUTPUT:
      !! ----------------
      !!    * xz0         : return the aerodynamic roughness length (integration constant for wind stress) [m]
      !!    * xu_star     : return u* the friction velocity                    [m/s]
      !!    * xL          : return the Monin-Obukhov length                    [m]
      !!    * xUN10       : return the Monin-Obukhov length                    [m/s]
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(in   )                     ::   zt       ! height for t_zt and q_zt                    [m]
      REAL(wp), INTENT(in   )                     ::   zu       ! height for U_zu                             [m]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   sst      ! sea surface temperature                [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   t_zt     ! potential air temperature              [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   ssq      ! sea surface specific humidity           [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   q_zt     ! specific air humidity                   [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   U_zu     ! relative wind module at zu                [m/s]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Cd       ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ch       ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ce       ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   t_zu     ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   q_zu     ! spec. humidity adjusted at zu           [kg/kg]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   U_blk    ! bulk wind at 10m                          [m/s]
      !
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xz0  ! Aerodynamic roughness length   [m]
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xu_star  ! u*, friction velocity
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xL  ! zeta (zu/L)
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xUN10  ! Neutral wind at zu
      !
      INTEGER :: j_itt
      LOGICAL :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   Cx_n10        ! 10m neutral latent/sensible coefficient
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   sqrt_Cd_n10   ! root square of Cd_n10
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   zeta_u        ! stability parameter at height zu
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   zpsi_h_u
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   ztmp0, ztmp1, ztmp2
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   stab          ! stability test integer
      !
      LOGICAL :: lreturn_z0=.FALSE., lreturn_ustar=.FALSE., lreturn_L=.FALSE., lreturn_UN10=.FALSE.
      !!----------------------------------------------------------------------------------

      ALLOCATE( Cx_n10(jpi,jpj), sqrt_Cd_n10(jpi,jpj), &
         &    zeta_u(jpi,jpj), stab(jpi,jpj), zpsi_h_u(jpi,jpj),  &
         &    ztmp0(jpi,jpj),  ztmp1(jpi,jpj), ztmp2(jpi,jpj) )

      IF( PRESENT(xz0) )     lreturn_z0    = .TRUE.
      IF( PRESENT(xu_star) ) lreturn_ustar = .TRUE.
      IF( PRESENT(xL) )      lreturn_L     = .TRUE.
      IF( PRESENT(xUN10) )   lreturn_UN10  = .TRUE.


      l_zt_equal_zu = .FALSE.
      IF( ABS(zu - zt) < 0.01 ) l_zt_equal_zu = .TRUE.    ! testing "zu == zt" is risky with double precision

      U_blk = MAX( 0.5_wp , U_zu )   !  relative wind speed at zu (normally 10m), we don't want to fall under 0.5 m/s

      !! First guess of stability:
      ztmp0 = virt_temp(t_zt, q_zt) - virt_temp(sst, ssq) ! air-sea difference of virtual pot. temp. at zt
      stab  = 0.5_wp + sign(0.5_wp,ztmp0)                           ! stab = 1 if dTv > 0  => STABLE, 0 if unstable

      !! As a first guess: initializing transf. coeff. with the Neutral coefficients at 10m:
      Cd = cd_neutral_10m( U_blk )
      sqrt_Cd_n10 = SQRT( Cd )
      Ce = 1.e-3*( 34.6 * sqrt_Cd_n10 )
      Ch = 1.e-3*sqrt_Cd_n10*(18.*stab + 32.7*(1. - stab))
      stab = sqrt_Cd_n10   ! Temporaty array !!! stab == SQRT(Cd)

      !! Initializing values at z_u with z_t values:
      t_zu = t_zt   ;   q_zu = q_zt

      !!  * Now starting iteration loop
      DO j_itt=1, nb_itt
         !
         ztmp1 = t_zu - sst   ! Updating air/sea differences
         ztmp2 = q_zu - ssq

         ! Updating turbulent scales :   (L&Y 2004 eq. (7))
         ztmp0 = stab*U_blk       ! u*       (stab == SQRT(Cd))
         ztmp1 = Ch/stab*ztmp1    ! theta*   (stab == SQRT(Cd))
         ztmp2 = Ce/stab*ztmp2    ! q*       (stab == SQRT(Cd))

         ! Estimate the inverse of Monin-Obukov length (1/L) at height zu:
         ztmp0 = One_on_L( t_zu, q_zu, ztmp0, ztmp1, ztmp2 )
         
         !! Stability parameters :
         zeta_u   = zu*ztmp0
         zeta_u = sign( min(abs(zeta_u),10.0_wp), zeta_u )
         zpsi_h_u = psi_h( zeta_u )

         !! Shifting temperature and humidity at zu (L&Y 2004 eq. (9b-9c))
         IF( .NOT. l_zt_equal_zu ) THEN
            !! Array 'stab' is free for the moment so using it to store 'zeta_t'
            stab = zt*ztmp0
            stab = SIGN( MIN(ABS(stab),10.0_wp), stab )  ! Temporaty array stab == zeta_t !!!
            stab = LOG(zt/zu) + zpsi_h_u - psi_h(stab)                   ! stab just used as temp array again!
            t_zu = t_zt - ztmp1/vkarmn*stab    ! ztmp1 is still theta*  L&Y 2004 eq.(9b)
            q_zu = q_zt - ztmp2/vkarmn*stab    ! ztmp2 is still q*      L&Y 2004 eq.(9c)
            q_zu = max(0._wp, q_zu)
         END IF

         ! Update neutral wind speed at 10m and neutral Cd at 10m (L&Y 2004 eq. 9a)...
         !   In very rare low-wind conditions, the old way of estimating the
         !   neutral wind speed at 10m leads to a negative value that causes the code
         !   to crash. To prevent this a threshold of 0.25m/s is imposed.
         ztmp2 = psi_m(zeta_u)
         ztmp0 = MAX( 0.25_wp , U_blk/(1._wp + sqrt_Cd_n10/vkarmn*(LOG(zu/10._wp) - ztmp2)) ) ! U_n10 (ztmp2 == psi_m(zeta_u))
         ztmp0 = cd_neutral_10m(ztmp0)                                               ! Cd_n10
         sqrt_Cd_n10 = sqrt(ztmp0)

         stab    = 0.5_wp + sign(0.5_wp,zeta_u)                        ! update stability
         Cx_n10  = 1.e-3*sqrt_Cd_n10*(18.*stab + 32.7*(1. - stab))  ! L&Y 2004 eq. (6c-6d)    (Cx_n10 == Ch_n10)

         !! Update of transfer coefficients:
         ztmp1 = 1. + sqrt_Cd_n10/vkarmn*(LOG(zu/10.) - ztmp2)   ! L&Y 2004 eq. (10a) (ztmp2 == psi_m(zeta_u))
         Cd      = ztmp0 / ( ztmp1*ztmp1 )
         stab = SQRT( Cd ) ! Temporary array !!! (stab == SQRT(Cd))

         ztmp0 = (LOG(zu/10.) - zpsi_h_u) / vkarmn / sqrt_Cd_n10
         ztmp2 = stab / sqrt_Cd_n10   ! (stab == SQRT(Cd))
         ztmp1 = 1. + Cx_n10*ztmp0    ! (Cx_n10 == Ch_n10)
         Ch  = Cx_n10*ztmp2 / ztmp1   ! L&Y 2004 eq. (10b)

         Cx_n10  = 1.e-3 * (34.6 * sqrt_Cd_n10)  ! L&Y 2004 eq. (6b)    ! Cx_n10 == Ce_n10
         ztmp1 = 1. + Cx_n10*ztmp0
         Ce  = Cx_n10*ztmp2 / ztmp1  ! L&Y 2004 eq. (10c)

      END DO

      IF( lreturn_z0 )    xz0     = zu*EXP( -(vkarmn/sqrt_Cd_n10) )
      IF( lreturn_ustar ) xu_star = stab*U_blk
      IF( lreturn_L )     xL      = zu/zeta_u
      IF( lreturn_UN10 )  xUN10   = U_blk/(1. + sqrt_Cd_n10/vkarmn*(LOG(zu/10.) - psi_m(zeta_u)))

      
      DEALLOCATE( Cx_n10, sqrt_Cd_n10, &
         &    zeta_u, stab, zpsi_h_u, ztmp0, &
         &    ztmp1, ztmp2 )

   END SUBROUTINE turb_ncar


   FUNCTION cd_neutral_10m( pw10 )
      !!----------------------------------------------------------------------------------      
      !! Estimate of the neutral drag coefficient at 10m as a function
      !! of neutral wind  speed at 10m
      !!
      !! Origin: Large & Yeager 2008 eq.(11a) and eq.(11b)
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pw10           ! scalar wind speed at 10m (m/s)
      REAL(wp), DIMENSION(jpi,jpj)             :: cd_neutral_10m
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
            cd_neutral_10m(ji,jj) = 1.e-3 * ( &
               &       (1. - zgt33)*( 2.7/zw + 0.142 + zw/13.09 - 3.14807E-10*zw6) & ! wind <  33 m/s
               &      +    zgt33   *      2.34 )                                     ! wind >= 33 m/s
            !
            cd_neutral_10m(ji,jj) = MAX(cd_neutral_10m(ji,jj), 1.E-6_wp)
            !
         END DO
      END DO
      !
   END FUNCTION cd_neutral_10m


   FUNCTION psi_m( pzeta )
      !!----------------------------------------------------------------------------------
      !! Universal profile stability function for momentum
      !!    !! Psis, L&Y 2004 eq. (8c), (8d), (8e)
      !!     
      !! pzet0 : stability paramenter, z/L where z is altitude measurement                                          
      !!         and L is M-O length
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pzeta
      REAL(wp), DIMENSION(jpi,jpj)             ::   psi_m
      !
      INTEGER  ::   ji, jj         ! dummy loop indices
      REAL(wp) :: zx2, zx, zstab   ! local scalars
      !!----------------------------------------------------------------------------------
      !
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
      !
   END FUNCTION psi_m


   FUNCTION psi_h( pzeta )
      !!----------------------------------------------------------------------------------
      !! Universal profile stability function for temperature and humidity
      !!    !! Psis, L&Y 2004 eq. (8c), (8d), (8e)
      !!
      !! pzet0 : stability paramenter, z/L where z is altitude measurement                                          
      !!         and L is M-O length
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      REAL(wp), DIMENSION(jpi,jpj)             :: psi_h
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
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
      !
   END FUNCTION psi_h

   !!======================================================================
END MODULE mod_blk_ncar
