MODULE sbcblk_algo_ncar
   !!======================================================================
   !!                   ***  MODULE  sbcblk_algo_ncar  ***
   !! Computes:
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (2m) to zu (10m) if needed
   !!   * the effective bulk wind speed at 10m U_blk
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!    Using the bulk formulation/param. of Large & Yeager 2008
   !!
   !!       Routine turb_ncar maintained and developed in AeroBulk
   !!                     (http://aerobulk.sourceforge.net/)
   !!
   !!                         L. Brodeau, 2015
   !!=====================================================================
   !! History :  3.6  !  2016-02  (L.Brodeau) successor of old turb_ncar of former sbcblk_core.F90
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   turb_ncar  : computes the bulk turbulent transfer coefficients
   !!                   adjusts t_air and q_air from zt to zu m
   !!                   returns the effective bulk wind speed at 10m
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE sbcwave, ONLY   :  cdn_wave ! wave module
#if defined key_si3 || defined key_cice
   USE sbc_ice         ! Surface boundary condition: ice fields
#endif
   !
   USE iom             ! I/O manager library
   USE lib_mpp         ! distribued memory computing library
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control
   USE lib_fortran     ! to use key_nosignedzero


   IMPLICIT NONE
   PRIVATE

   PUBLIC ::   TURB_NCAR   ! called by sbcblk.F90

   !                              ! NCAR own values for given constants:
   REAL(wp), PARAMETER ::   rctv0 = 0.608   ! constant to obtain virtual temperature...
   
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE turb_ncar( zt, zu, sst, t_zt, ssq, q_zt, U_zu, &
      &                  Cd, Ch, Ce, t_zu, q_zu, U_blk,      &
      &                  Cdn, Chn, Cen                       )
      !!----------------------------------------------------------------------------------
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
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
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
      !!    *  U_blk  : bulk wind at 10m                                      [m/s]
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
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Cdn, Chn, Cen ! neutral transfer coefficients
      !
      INTEGER ::   j_itt
      LOGICAL ::   l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      INTEGER , PARAMETER ::   nb_itt = 4       ! number of itterations
      !
      REAL(wp), DIMENSION(jpi,jpj) ::   Cx_n10        ! 10m neutral latent/sensible coefficient
      REAL(wp), DIMENSION(jpi,jpj) ::   sqrt_Cd_n10   ! root square of Cd_n10
      REAL(wp), DIMENSION(jpi,jpj) ::   zeta_u        ! stability parameter at height zu
      REAL(wp), DIMENSION(jpi,jpj) ::   zpsi_h_u
      REAL(wp), DIMENSION(jpi,jpj) ::   ztmp0, ztmp1, ztmp2
      REAL(wp), DIMENSION(jpi,jpj) ::   stab          ! stability test integer
      !!----------------------------------------------------------------------------------
      !
      l_zt_equal_zu = .FALSE.
      IF( ABS(zu - zt) < 0.01 ) l_zt_equal_zu = .TRUE.    ! testing "zu == zt" is risky with double precision

      U_blk = MAX( 0.5 , U_zu )   !  relative wind speed at zu (normally 10m), we don't want to fall under 0.5 m/s

      !! First guess of stability:
      ztmp0 = t_zt*(1. + rctv0*q_zt) - sst*(1. + rctv0*ssq) ! air-sea difference of virtual pot. temp. at zt
      stab  = 0.5 + sign(0.5,ztmp0)                           ! stab = 1 if dTv > 0  => STABLE, 0 if unstable

      !! Neutral coefficients at 10m:
      IF( ln_cdgw ) THEN      ! wave drag case
         cdn_wave(:,:) = cdn_wave(:,:) + rsmall * ( 1._wp - tmask(:,:,1) )
         ztmp0   (:,:) = cdn_wave(:,:)
      ELSE
         ztmp0 = cd_neutral_10m( U_blk )
      ENDIF

      sqrt_Cd_n10 = SQRT( ztmp0 )

      !! Initializing transf. coeff. with their first guess neutral equivalents :
      Cd = ztmp0
      Ce = 1.e-3*( 34.6 * sqrt_Cd_n10 )
      Ch = 1.e-3*sqrt_Cd_n10*(18.*stab + 32.7*(1. - stab))
      stab = sqrt_Cd_n10   ! Temporaty array !!! stab == SQRT(Cd)
 
      IF( ln_cdgw )   Cen = Ce  ; Chn = Ch

      !! Initializing values at z_u with z_t values:
      t_zu = t_zt   ;   q_zu = q_zt

      !!  * Now starting iteration loop
      DO j_itt=1, nb_itt
         !
         ztmp1 = t_zu - sst   ! Updating air/sea differences
         ztmp2 = q_zu - ssq

         ! Updating turbulent scales :   (L&Y 2004 eq. (7))
         ztmp1  = Ch/stab*ztmp1    ! theta*   (stab == SQRT(Cd))
         ztmp2  = Ce/stab*ztmp2    ! q*       (stab == SQRT(Cd))

         ztmp0 = 1. + rctv0*q_zu      ! multiply this with t and you have the virtual temperature

         ! Estimate the inverse of Monin-Obukov length (1/L) at height zu:
         ztmp0 =  (grav*vkarmn/(t_zu*ztmp0)*(ztmp1*ztmp0 + rctv0*t_zu*ztmp2)) / (Cd*U_blk*U_blk)
         !                                                      ( Cd*U_blk*U_blk is U*^2 at zu )

         !! Stability parameters :
         zeta_u   = zu*ztmp0   ;  zeta_u = sign( min(abs(zeta_u),10.0), zeta_u )
         zpsi_h_u = psi_h( zeta_u )

         !! Shifting temperature and humidity at zu (L&Y 2004 eq. (9b-9c))
         IF( .NOT. l_zt_equal_zu ) THEN
            !! Array 'stab' is free for the moment so using it to store 'zeta_t'
            stab = zt*ztmp0 ;  stab = SIGN( MIN(ABS(stab),10.0), stab )  ! Temporaty array stab == zeta_t !!!
            stab = LOG(zt/zu) + zpsi_h_u - psi_h(stab)                   ! stab just used as temp array again!
            t_zu = t_zt - ztmp1/vkarmn*stab    ! ztmp1 is still theta*  L&Y 2004 eq.(9b)
            q_zu = q_zt - ztmp2/vkarmn*stab    ! ztmp2 is still q*      L&Y 2004 eq.(9c)
            q_zu = max(0., q_zu)
         END IF

         ztmp2 = psi_m(zeta_u)
         IF( ln_cdgw ) THEN      ! surface wave case
            stab = vkarmn / ( vkarmn / sqrt_Cd_n10 - ztmp2 )  ! (stab == SQRT(Cd))
            Cd   = stab * stab
            ztmp0 = (LOG(zu/10.) - zpsi_h_u) / vkarmn / sqrt_Cd_n10
            ztmp2 = stab / sqrt_Cd_n10   ! (stab == SQRT(Cd))
            ztmp1 = 1. + Chn * ztmp0     
            Ch    = Chn * ztmp2 / ztmp1  ! L&Y 2004 eq. (10b)
            ztmp1 = 1. + Cen * ztmp0
            Ce    = Cen * ztmp2 / ztmp1  ! L&Y 2004 eq. (10c)

         ELSE
            ! Update neutral wind speed at 10m and neutral Cd at 10m (L&Y 2004 eq. 9a)...
            !   In very rare low-wind conditions, the old way of estimating the
            !   neutral wind speed at 10m leads to a negative value that causes the code
            !   to crash. To prevent this a threshold of 0.25m/s is imposed.
            ztmp0 = MAX( 0.25 , U_blk/(1. + sqrt_Cd_n10/vkarmn*(LOG(zu/10.) - ztmp2)) ) ! U_n10 (ztmp2 == psi_m(zeta_u))
            ztmp0 = cd_neutral_10m(ztmp0)                                               ! Cd_n10
            Cdn(:,:) = ztmp0
            sqrt_Cd_n10 = sqrt(ztmp0)

            stab    = 0.5 + sign(0.5,zeta_u)                           ! update stability
            Cx_n10  = 1.e-3*sqrt_Cd_n10*(18.*stab + 32.7*(1. - stab))  ! L&Y 2004 eq. (6c-6d)    (Cx_n10 == Ch_n10)
            Chn(:,:) = Cx_n10

            !! Update of transfer coefficients:
            ztmp1 = 1. + sqrt_Cd_n10/vkarmn*(LOG(zu/10.) - ztmp2)   ! L&Y 2004 eq. (10a) (ztmp2 == psi_m(zeta_u))
            Cd      = ztmp0 / ( ztmp1*ztmp1 )
            stab = SQRT( Cd ) ! Temporary array !!! (stab == SQRT(Cd))

            ztmp0 = (LOG(zu/10.) - zpsi_h_u) / vkarmn / sqrt_Cd_n10
            ztmp2 = stab / sqrt_Cd_n10   ! (stab == SQRT(Cd))
            ztmp1 = 1. + Cx_n10*ztmp0    ! (Cx_n10 == Ch_n10)
            Ch  = Cx_n10*ztmp2 / ztmp1   ! L&Y 2004 eq. (10b)

            Cx_n10  = 1.e-3 * (34.6 * sqrt_Cd_n10)  ! L&Y 2004 eq. (6b)    ! Cx_n10 == Ce_n10
            Cen(:,:) = Cx_n10
            ztmp1 = 1. + Cx_n10*ztmp0
            Ce  = Cx_n10*ztmp2 / ztmp1  ! L&Y 2004 eq. (10c)
            ENDIF
         !
      END DO
      !
   END SUBROUTINE turb_ncar


   FUNCTION cd_neutral_10m( pw10 )
      !!----------------------------------------------------------------------------------      
      !! Estimate of the neutral drag coefficient at 10m as a function
      !! of neutral wind  speed at 10m
      !!
      !! Origin: Large & Yeager 2008 eq.(11a) and eq.(11b)
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
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
            zgt33 = 0.5 + SIGN( 0.5, (zw - 33.) )   ! If pw10 < 33. => 0, else => 1
            !
            cd_neutral_10m(ji,jj) = 1.e-3 * ( &
               &       (1. - zgt33)*( 2.7/zw + 0.142 + zw/13.09 - 3.14807E-10*zw6) & ! wind <  33 m/s
               &      +    zgt33   *      2.34 )                                     ! wind >= 33 m/s
            !
            cd_neutral_10m(ji,jj) = MAX(cd_neutral_10m(ji,jj), 1.E-6)
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
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
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
            zx2 = SQRT( ABS( 1. - 16.*pzeta(ji,jj) ) )
            zx2 = MAX ( zx2 , 1. )
            zx  = SQRT( zx2 )
            zstab = 0.5 + SIGN( 0.5 , pzeta(ji,jj) )
            !
            psi_m(ji,jj) =        zstab  * (-5.*pzeta(ji,jj))       &          ! Stable
               &          + (1. - zstab) * (2.*LOG((1. + zx)*0.5)   &          ! Unstable
               &               + LOG((1. + zx2)*0.5) - 2.*ATAN(zx) + rpi*0.5)  !    "
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
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
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
            zx2 = SQRT( ABS( 1. - 16.*pzeta(ji,jj) ) )
            zx2 = MAX ( zx2 , 1. )
            zstab = 0.5 + SIGN( 0.5 , pzeta(ji,jj) )
            !
            psi_h(ji,jj) =         zstab  * (-5.*pzeta(ji,jj))        &  ! Stable
               &           + (1. - zstab) * (2.*LOG( (1. + zx2)*0.5 ))   ! Unstable
            !
         END DO
      END DO
      !
   END FUNCTION psi_h

   !!======================================================================
END MODULE sbcblk_algo_ncar
