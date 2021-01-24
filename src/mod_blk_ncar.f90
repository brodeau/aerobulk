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
   !!                     (http://aerobulk.sourceforge.net/)
   !!
   !!            Author: Laurent Brodeau, 2016
   !!
   !!====================================================================================
   USE mod_const   !: physical and othe constants
   !USE mod_thermo  !: thermodynamics

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: TURB_NCAR

CONTAINS

   SUBROUTINE turb_ncar( zt, zu, sst, T_zt, q_sat, q_zt, dU,    &
      &                      Cd, Ch, Ce , T_zu, q_zu )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_core  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to Large & Yeager (2004) and Large & Yeager (2008)
      !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
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
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   )                     ::   zt       ! height for T_zt and q_zt                   [m]
      REAL(wp), INTENT(in   )                     ::   zu       ! height for dU                              [m]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   sst      ! sea surface temperature              [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   T_zt     ! potential air temperature            [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   q_sat    ! sea surface specific humidity         [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   q_zt     ! specific air humidity                 [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   dU       ! relative wind module at zu            [m/s]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Cd       ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ch       ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ce       ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   T_zu     ! air temp. shifted at zu                     [K]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   q_zu     ! spec. hum.  shifted at zu               [kg/kg]
      !
      INTEGER ::   j_itt
      LOGICAL ::   l_zt_equal_zu = .FALSE.      ! if q and t are given at different height than U
      !
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   U_zu          ! relative wind at zu                            [m/s]
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   Ce_n10        ! 10m neutral latent coefficient
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   Ch_n10        ! 10m neutral sensible coefficient
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   sqrt_Cd_n10   ! root square of Cd_n10
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   sqrt_Cd       ! root square of Cd
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   zeta_u        ! stability parameter at height zu
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   zeta_t        ! stability parameter at height zt
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   zpsi_h_u, zpsi_m_u
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   ztmp0, ztmp1, ztmp2
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   stab          ! 1st stability test integer
      !!----------------------------------------------------------------------

      ALLOCATE( U_zu(jpi,jpj), Ce_n10(jpi,jpj), Ch_n10(jpi,jpj), sqrt_Cd_n10(jpi,jpj), sqrt_Cd(jpi,jpj) )
      ALLOCATE( zeta_u(jpi,jpj), stab(jpi,jpj) )
      ALLOCATE( zpsi_h_u(jpi,jpj), zpsi_m_u(jpi,jpj), ztmp0(jpi,jpj), ztmp1(jpi,jpj), ztmp2(jpi,jpj) )

      l_zt_equal_zu = .FALSE.
      IF( ABS(zu - zt) < 0.01 ) l_zt_equal_zu = .TRUE.    ! testing "zu == zt" is risky with double precision

      IF( .NOT. l_zt_equal_zu )   ALLOCATE( zeta_t(jpi,jpj) )

      U_zu = MAX( 0.5 , dU )   !  relative wind speed at zu (normally 10m), we don't want to fall under 0.5 m/s

      !! First guess of stability:
      ztmp0 = T_zt*(1. + 0.608*q_zt) - sst*(1. + 0.608*q_sat) ! air-sea difference of virtual pot. temp. at zt
      stab  = 0.5 + sign(0.5,ztmp0)                           ! stab = 1 if dTv > 0  => STABLE, 0 if unstable

      !! Neutral coefficients at 10m:
      ztmp0 = cd_neutral_10m( U_zu )

      sqrt_Cd_n10 = SQRT( ztmp0 )
      Ce_n10  = 1.e-3*( 34.6 * sqrt_Cd_n10 )
      Ch_n10  = 1.e-3*sqrt_Cd_n10*(18.*stab + 32.7*(1. - stab))

      !! Initializing transf. coeff. with their first guess neutral equivalents :
      Cd = ztmp0   ;   Ce = Ce_n10   ;   Ch = Ch_n10   ;   sqrt_Cd = sqrt_Cd_n10

      !! Initializing values at z_u with z_t values:
      T_zu = T_zt   ;   q_zu = q_zt

      !!  * Now starting iteration loop
      DO j_itt=1, 5     !lolo crude
         !
         ztmp1 = T_zu - sst   ! Updating air/sea differences
         ztmp2 = q_zu - q_sat

         ! Updating turbulent scales :   (L&Y 2004 eq. (7))
         ztmp1  = Ch/sqrt_Cd*ztmp1    ! theta*
         ztmp2  = Ce/sqrt_Cd*ztmp2    ! q*

         ztmp0 = T_zu*(1. + 0.608*q_zu) ! virtual potential temperature at zu

         ! Estimate the inverse of Monin-Obukov length (1/L) at height zu:
         ztmp0 =  (vkarmn*grav/ztmp0*(ztmp1*(1.+0.608*q_zu) + 0.608*T_zu*ztmp2)) / (Cd*U_zu*U_zu)
         !                                                                     ( Cd*U_zu*U_zu is U*^2 at zu)

         !! Stability parameters :
         zeta_u   = zu*ztmp0   ;  zeta_u = sign( min(abs(zeta_u),10.0), zeta_u )
         zpsi_h_u = psi_h( zeta_u )
         zpsi_m_u = psi_m( zeta_u )

         !! Shifting temperature and humidity at zu (L&Y 2004 eq. (9b-9c))
         IF ( .NOT. l_zt_equal_zu ) THEN
            zeta_t = zt*ztmp0 ;  zeta_t = sign( min(abs(zeta_t),10.0), zeta_t )
            stab = LOG(zu/zt) - zpsi_h_u + psi_h(zeta_t)  ! stab just used as temp array!!!
            T_zu = T_zt + ztmp1/vkarmn*stab    ! ztmp1 is still theta*
            q_zu = q_zt + ztmp2/vkarmn*stab    ! ztmp2 is still q*
            q_zu = max(0., q_zu)
         END IF

         ! Update neutral wind speed at 10m and neutral Cd at 10m (L&Y 2004 eq. 9a)...
         !   In very rare low-wind conditions, the old way of estimating the
         !   neutral wind speed at 10m leads to a negative value that causes the code
         !   to crash. To prevent this a threshold of 0.25m/s is imposed.
         ztmp0 = MAX( 0.25 , U_zu/(1. + sqrt_Cd_n10/vkarmn*(LOG(zu/10.) - zpsi_m_u)) ) !  U_n10
         ztmp0 = cd_neutral_10m(ztmp0)                                                 ! Cd_n10
         sqrt_Cd_n10 = sqrt(ztmp0)

         Ce_n10  = 1.e-3 * (34.6 * sqrt_Cd_n10)                     ! L&Y 2004 eq. (6b)
         stab    = 0.5 + sign(0.5,zeta_u)                           ! update stability
         Ch_n10  = 1.e-3*sqrt_Cd_n10*(18.*stab + 32.7*(1. - stab))  ! L&Y 2004 eq. (6c-6d)

         !! Update of transfer coefficients:
         ztmp1 = 1. + sqrt_Cd_n10/vkarmn*(LOG(zu/10.) - zpsi_m_u)   ! L&Y 2004 eq. (10a)
         Cd      = ztmp0 / ( ztmp1*ztmp1 )
         sqrt_Cd = SQRT( Cd )
         !
         ztmp0 = (LOG(zu/10.) - zpsi_h_u) / vkarmn / sqrt_Cd_n10
         ztmp2 = sqrt_Cd / sqrt_Cd_n10
         ztmp1 = 1. + Ch_n10*ztmp0
         Ch  = Ch_n10*ztmp2 / ztmp1  ! L&Y 2004 eq. (10b)
         !
         ztmp1 = 1. + Ce_n10*ztmp0
         Ce  = Ce_n10*ztmp2 / ztmp1  ! L&Y 2004 eq. (10c)
         !
      END DO

      DEALLOCATE( U_zu, Ce_n10, Ch_n10, sqrt_Cd_n10, sqrt_Cd )
      DEALLOCATE( zeta_u, stab )
      DEALLOCATE( zpsi_h_u, zpsi_m_u, ztmp0, ztmp1, ztmp2 )

      IF( .NOT. l_zt_equal_zu ) DEALLOCATE( zeta_t )

   END SUBROUTINE turb_ncar


   FUNCTION cd_neutral_10m( zw10 )
      !!----------------------------------------------------------------------
      !! Estimate of the neutral drag coefficient at 10m as a function
      !! of neutral wind  speed at 10m
      !!
      !! Origin: Large & Yeager 2008 eq.(11a) and eq.(11b)
      !!
      !! Author: L. Brodeau, june 2014
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   zw10           ! scalar wind speed at 10m (m/s)
      REAL(wp), DIMENSION(jpi,jpj)             ::   cd_neutral_10m
      !
      REAL(wp), DIMENSION(:,:), POINTER ::   rgt33
      !!----------------------------------------------------------------------
      !
      ALLOCATE( rgt33(jpi,jpj) )
      !
      !! When wind speed > 33 m/s => Cyclone conditions => special treatment
      rgt33 = 0.5_wp + SIGN( 0.5_wp, (zw10 - 33._wp) )   ! If zw10 < 33. => 0, else => 1
      cd_neutral_10m = 1.e-3 * ( &
         &       (1._wp - rgt33)*( 2.7_wp/zw10 + 0.142_wp + zw10/13.09_wp - 3.14807E-10*zw10**6) & ! zw10< 33.
         &      + rgt33         *      2.34   )                                                    ! zw10 >= 33.
      !
      DEALLOCATE( rgt33)
      !
   END FUNCTION cd_neutral_10m


   FUNCTION psi_m(pta)   !! Psis, L&Y 2004 eq. (8c), (8d), (8e)
      !-------------------------------------------------------------------------------
      ! universal profile stability function for momentum
      !-------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pta
      !
      REAL(wp), DIMENSION(jpi,jpj)             :: psi_m
      REAL(wp), DIMENSION(:,:), POINTER        :: X2, X, stabit
      !-------------------------------------------------------------------------------
      !
      ALLOCATE( X2(jpi,jpj), X(jpi,jpj), stabit(jpi,jpj) )
      !
      X2 = SQRT( ABS( 1. - 16.*pta ) )  ;  X2 = MAX( X2 , 1. )   ;   X = SQRT( X2 )
      stabit = 0.5 + SIGN( 0.5 , pta )
      psi_m = -5.*pta*stabit  &                                                          ! Stable
         &    + (1. - stabit)*(2.*LOG((1. + X)*0.5) + LOG((1. + X2)*0.5) - 2.*ATAN(X) + rpi*0.5)  ! Unstable
      !
      DEALLOCATE( X2, X, stabit )
      !
   END FUNCTION psi_m


   FUNCTION psi_h( pta )    !! Psis, L&Y 2004 eq. (8c), (8d), (8e)
      !-------------------------------------------------------------------------------
      ! universal profile stability function for temperature and humidity
      !-------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pta
      !
      REAL(wp), DIMENSION(jpi,jpj)             ::   psi_h
      REAL(wp), DIMENSION(:,:), POINTER        ::   X2, X, stabit
      !-------------------------------------------------------------------------------
      !
      ALLOCATE( X2(jpi,jpj), X(jpi,jpj), stabit(jpi,jpj) )
      !
      X2 = SQRT( ABS( 1. - 16.*pta ) )   ;   X2 = MAX( X2 , 1. )   ;   X = SQRT( X2 )
      stabit = 0.5 + SIGN( 0.5 , pta )
      psi_h = -5.*pta*stabit   &                                       ! Stable
         &    + (1. - stabit)*(2.*LOG( (1. + X2)*0.5 ))                ! Unstable
      !
      DEALLOCATE( X2, X, stabit )
      !
   END FUNCTION psi_h

   !!======================================================================
END MODULE mod_blk_ncar
