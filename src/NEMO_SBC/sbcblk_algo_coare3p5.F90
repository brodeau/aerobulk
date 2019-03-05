MODULE sbcblk_algo_coare3p5
   !!======================================================================
   !!                       ***  MODULE  sbcblk_algo_coare3p5  ***
   !! Computes turbulent components of surface fluxes
   !!         according to Edson et al. 2013 (COARE v3.5) /JPO
   !!
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (2m) to zu (10m) if needed
   !!   * the effective bulk wind speed at 10m U_blk
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!    Using the bulk formulation/param. of COARE v3.5, Edson et al. 2013
   !!
   !!
   !!       Routine turb_coare3p5 maintained and developed in AeroBulk
   !!                     (http://aerobulk.sourceforge.net/)
   !!
   !!            Author: Laurent Brodeau, 2016, brodeau@gmail.com
   !!
   !!======================================================================
   !! History :  3.6  !  2016-02  (L.Brodeau)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   turb_coare3p5  : computes the bulk turbulent transfer coefficients
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

   PUBLIC ::   TURB_COARE3P5   ! called by sbcblk.F90

   !                                   ! COARE own values for given constants:
   REAL(wp), PARAMETER ::   charn0_max =   0.028   ! value above which the Charnock paramter levels off for winds > 18
   REAL(wp), PARAMETER ::   zi0        = 600.      ! scale height of the atmospheric boundary layer...1
   REAL(wp), PARAMETER ::   Beta0      =   1.25    ! gustiness parameter
   REAL(wp), PARAMETER ::   rctv0      =   0.608   ! constant to obtain virtual temperature...

   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE turb_coare3p5( zt, zu, sst, t_zt, ssq, q_zt, U_zu,  &
      &                      Cd, Ch, Ce, t_zu, q_zu, U_blk,       &
      &                      Cdn, Chn, Cen                        )
      !!----------------------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_coare3p5  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to Fairall et al. (2003)
      !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!
      !! ** Method : Monin Obukhov Similarity Theory
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
      !!
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(in   )                     ::   zt       ! height for t_zt and q_zt                   [m]
      REAL(wp), INTENT(in   )                     ::   zu       ! height for U_zu                              [m]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   sst      ! sea surface temperature              [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   t_zt     ! potential air temperature            [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   ssq      ! sea surface specific humidity         [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   q_zt     ! specific air humidity                 [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   U_zu     ! relative wind module at zu            [m/s]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Cd       ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ch       ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ce       ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   t_zu     ! pot. air temp. adjusted at zu             [K]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   q_zu     ! spec. humidity adjusted at zu             [kg/kg]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   U_blk    ! bulk wind at 10m                          [m/s]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Cdn, Chn, Cen ! neutral transfer coefficients
      !
      INTEGER :: j_itt
      LOGICAL ::   l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      INTEGER , PARAMETER ::   nb_itt = 4       ! number of itterations
      !
      REAL(wp), DIMENSION(jpi,jpj) ::  &
         &  u_star, t_star, q_star, &
         &  dt_zu, dq_zu,    &
         &  znu_a,           & !: Nu_air, Viscosity of air
         &  z0, z0t
      REAL(wp), DIMENSION(jpi,jpj) ::   zeta_u        ! stability parameter at height zu
      REAL(wp), DIMENSION(jpi,jpj) ::   ztmp0, ztmp1, ztmp2
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   zeta_t        ! stability parameter at height zt
      !!----------------------------------------------------------------------------------
      !
      l_zt_equal_zu = .FALSE.
      IF( ABS(zu - zt) < 0.01 ) l_zt_equal_zu = .TRUE.    ! testing "zu == zt" is risky with double precision

      IF( .NOT. l_zt_equal_zu )   ALLOCATE( zeta_t(jpi,jpj) )

      !! First guess of temperature and humidity at height zu:
      t_zu = MAX(t_zt , 0.0)    ! who knows what's given on masked-continental regions...
      q_zu = MAX(q_zt , 1.E-6)  !               "

      !! Pot. temp. difference (and we don't want it to be 0!)
      dt_zu = t_zu - sst ;  dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6), dt_zu )
      dq_zu = q_zu - ssq ;  dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9), dq_zu )

      znu_a = visc_air(t_zt) ! Air viscosity (m^2/s) at zt given from temperature in (K)

      ztmp2 = 0.5*0.5  ! initial guess for wind gustiness contribution
      U_blk = SQRT(U_zu*U_zu + ztmp2)

      ztmp2   = 10000.     ! optimization: ztmp2 == 1/z0 (with z0 first guess == 0.0001)
      ztmp0   = LOG(zu*ztmp2)
      ztmp1   = LOG(10.*ztmp2)
      u_star = 0.035*U_blk*ztmp1/ztmp0       ! (u* = 0.035*Un10)

      !! COARE 3.5 first guess of UN10 is U_zu
      ztmp2 = MIN( 0.0017*U_zu - 0.005 , charn0_max)  ! alpha Charnock parameter (Eq. 13 Edson al. 2013)
      ztmp2 = MAX( ztmp2 , 0. )                       ! alpha Charnock parameter (Eq. 13 Edson al. 2013)
      z0     = ztmp2*u_star*u_star/grav + 0.11*znu_a/u_star
      z0t    = 0.1*EXP(vkarmn/(0.00115/(vkarmn/ztmp1)))   !  WARNING: 1/z0t !

      ztmp2  = vkarmn/ztmp0
      Cd     = ztmp2*ztmp2    ! first guess of Cd

      ztmp0 = vkarmn*vkarmn/LOG(zt*z0t)/Cd

      !Ribcu = -zu/(zi0*0.004*Beta0**3) !! Saturation Rib, zi0 = tropicalbound. layer depth
      ztmp2  = grav*zu*(dt_zu + rctv0*t_zu*dq_zu)/(t_zu*U_blk*U_blk)  !! Ribu Bulk Richardson number
      ztmp1 = 0.5 + sign(0.5 , ztmp2)
      ztmp0 = ztmp0*ztmp2
      !!             Ribu < 0                                 Ribu > 0   Beta = 1.25
      zeta_u = (1.-ztmp1) * (ztmp0/(1.+ztmp2/(-zu/(zi0*0.004*Beta0**3)))) &
         &  +     ztmp1   * (ztmp0*(1. + 27./9.*ztmp2/ztmp0))

      !! First guess M-O stability dependent scaling params.(u*,t*,q*) to estimate z0 and z/L
      ztmp0  =  vkarmn/(LOG(zu*z0t) - psi_h_coare(zeta_u))

      u_star = U_blk*vkarmn/(LOG(zu) - LOG(z0)  - psi_m_coare(zeta_u))
      t_star = dt_zu*ztmp0
      q_star = dq_zu*ztmp0

      ! What's need to be done if zt /= zu:
      IF( .NOT. l_zt_equal_zu ) THEN

         zeta_t = zt*zeta_u/zu

         !! First update of values at zu (or zt for wind)
         ztmp0 = psi_h_coare(zeta_u) - psi_h_coare(zeta_t)
         ztmp1 = log(zt/zu) + ztmp0
         t_zu = t_zt - t_star/vkarmn*ztmp1
         q_zu = q_zt - q_star/vkarmn*ztmp1
         q_zu = (0.5 + sign(0.5,q_zu))*q_zu !Makes it impossible to have negative humidity :

         dt_zu = t_zu - sst  ; dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6), dt_zu )
         dq_zu = q_zu - ssq  ; dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9), dq_zu )

      END IF

      !! ITERATION BLOCK
      DO j_itt = 1, nb_itt

         !!Inverse of Monin-Obukov length (1/L) :
         ztmp0 = One_on_L(t_zu, q_zu, u_star, t_star, q_star)  ! 1/L == 1/[Monin-Obukhov length]

         ztmp1 = u_star*u_star   ! u*^2

         !! Update wind at 10m taking into acount convection-related wind gustiness:
         ! Ug = Beta*w*  (Beta = 1.25, Fairall et al. 2003, Eq.8):
         ztmp2 = Beta0*Beta0*ztmp1*(MAX(-zi0*ztmp0/vkarmn,0.))**(2./3.)   ! => ztmp2 == Ug^2
         !!   ! Only true when unstable (L<0) => when ztmp0 < 0 => explains "-" before 600.
         U_blk = MAX(sqrt(U_zu*U_zu + ztmp2), 0.2)        ! include gustiness in bulk wind speed
         ! => 0.2 prevents U_blk to be 0 in stable case when U_zu=0.

         !! COARE 3.5: Charnock parameter is computed from the neutral wind speed at 10m: Eq. 13 (Edson al. 2013)
         ztmp2 = u_star/vkarmn*LOG(10./z0)   ! UN10 Neutral wind at 10m!
         ztmp2 = MIN( 0.0017*ztmp2 - 0.005 , charn0_max)  ! alpha Charnock parameter (Eq. 13 Edson al. 2013)
         ztmp2 = MAX( ztmp2 , 0. )

         !! Roughness lengthes z0, z0t (z0q = z0t) :
         z0   = ztmp2*ztmp1/grav + 0.11*znu_a/u_star ! Roughness length (eq.6)
         ztmp1 = z0*u_star/znu_a                             ! Re_r: roughness Reynolds number
         !z0t  = MIN( 1.1E-4 , 5.5E-5*ztmp1**(-0.6) )   ! COARE 3.0
         !! Chris Fairall and Jim Edsson, private communication, March 2016 / COARE 3.5 :
         z0t   = MIN( 1.6e-4 , 5.8E-5*ztmp1**(-0.72)) ! These thermal roughness lengths give Stanton and
         !z0q = z0t                                   ! Dalton numbers that closely approximate COARE3.0

         !! Stability parameters:
         zeta_u = zu*ztmp0 ; zeta_u = sign( min(abs(zeta_u),50.0), zeta_u )
         IF( .NOT. l_zt_equal_zu ) THEN
            zeta_t = zt*ztmp0 ;  zeta_t = sign( min(abs(zeta_t),50.0), zeta_t )
         END IF

         !! Turbulent scales at zu=10m :
         ztmp0   = psi_h_coare(zeta_u)
         ztmp1   = vkarmn/(LOG(zu) -LOG(z0t) - ztmp0)

         t_star = dt_zu*ztmp1
         q_star = dq_zu*ztmp1
         u_star = U_blk*vkarmn/(LOG(zu) -LOG(z0) - psi_m_coare(zeta_u))

         IF( .NOT. l_zt_equal_zu ) THEN
            ! What's need to be done if zt /= zu
            !! Re-updating temperature and humidity at zu :
            ztmp2 = ztmp0 - psi_h_coare(zeta_t)
            ztmp1 = log(zt/zu) + ztmp2
            t_zu = t_zt - t_star/vkarmn*ztmp1
            q_zu = q_zt - q_star/vkarmn*ztmp1
            dt_zu = t_zu - sst ;  dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6), dt_zu )
            dq_zu = q_zu - ssq ;  dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9), dq_zu )
         END IF

      END DO
      !
      ! compute transfer coefficients at zu :
      ztmp0 = u_star/U_blk
      Cd   = ztmp0*ztmp0
      Ch   = ztmp0*t_star/dt_zu
      Ce   = ztmp0*q_star/dq_zu
      !
      ztmp1 = zu + z0
      Cdn = vkarmn*vkarmn / (log(ztmp1/z0 )*log(ztmp1/z0 ))
      Chn = vkarmn*vkarmn / (log(ztmp1/z0t)*log(ztmp1/z0t))
      Cen = Chn
      !
      IF( .NOT. l_zt_equal_zu ) DEALLOCATE( zeta_t )
      !
   END SUBROUTINE turb_coare3p5



   FUNCTION One_on_L( ptha, pqa, pus, pts, pqs )
      !!------------------------------------------------------------------------
      !!
      !! Evaluates the 1./(Monin Obukhov length) from air temperature and
      !!  specific humidity, and frictional scales u*, t* and q*
      !!
      !! Author: L. Brodeau, june 2016 / AeroBulk
      !!         (https://sourceforge.net/p/aerobulk)
      !!------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: One_on_L         !: 1./(Monin Obukhov length) [m^-1]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: ptha,  &  !: average potetntial air temperature [K]
         &                                        pqa,   &  !: average specific humidity of air   [kg/kg]
         &                                      pus, pts, pqs   !: frictional velocity, temperature and humidity
      !
      INTEGER  ::   ji, jj         ! dummy loop indices
      REAL(wp) ::     zqa          ! local scalar
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zqa = (1. + rctv0*pqa(ji,jj))
            !
            One_on_L(ji,jj) =  grav*vkarmn*(pts(ji,jj)*zqa + rctv0*ptha(ji,jj)*pqs(ji,jj)) &
               &                      / ( pus(ji,jj)*pus(ji,jj) * ptha(ji,jj)*zqa )
            !
         END DO
      END DO
      !
   END FUNCTION One_on_L


   FUNCTION psi_m_coare( pzeta )
      !!----------------------------------------------------------------------------------
      !! ** Purpose: compute the universal profile stability function for momentum
      !!             COARE 3.0, Fairall et al. 2003
      !!             pzeta : stability paramenter, z/L where z is altitude
      !!                     measurement and L is M-O length
      !!       Stability function for wind speed and scalars matching Kansas and free
      !!       convection forms with weighting f convective form, follows Fairall et
      !!       al (1996) with profile constants from Grachev et al (2000) BLM stable
      !!       form from Beljaars and Holtslag (1991)
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: psi_m_coare
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) :: zta, zphi_m, zphi_c, zpsi_k, zpsi_c, zf, zc, zstab
      !!----------------------------------------------------------------------------------
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zta = pzeta(ji,jj)
            !
            zphi_m = ABS(1. - 15.*zta)**.25    !!Kansas unstable
            !
            zpsi_k = 2.*LOG((1. + zphi_m)/2.) + LOG((1. + zphi_m*zphi_m)/2.)   &
               & - 2.*ATAN(zphi_m) + 0.5*rpi
            !
            zphi_c = ABS(1. - 10.15*zta)**.3333                   !!Convective
            !
            zpsi_c = 1.5*LOG((1. + zphi_c + zphi_c*zphi_c)/3.) &
               &     - 1.7320508*ATAN((1. + 2.*zphi_c)/1.7320508) + 1.813799447
            !
            zf = zta*zta
            zf = zf/(1. + zf)
            zc = MIN(50., 0.35*zta)
            zstab = 0.5 + SIGN(0.5, zta)
            !
            psi_m_coare(ji,jj) = (1. - zstab) * ( (1. - zf)*zpsi_k + zf*zpsi_c ) & ! (zta < 0)
               &                -   zstab     * ( 1. + 1.*zta     &                ! (zta > 0)
               &                         + 0.6667*(zta - 14.28)/EXP(zc) + 8.525 )   !     "
            !
         END DO
      END DO
      !
   END FUNCTION psi_m_coare


   FUNCTION psi_h_coare( pzeta )
      !!---------------------------------------------------------------------
      !! Universal profile stability function for temperature and humidity
      !! COARE 3.0, Fairall et al. 2003
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! Stability function for wind speed and scalars matching Kansas and free
      !! convection forms with weighting f convective form, follows Fairall et
      !! al (1996) with profile constants from Grachev et al (2000) BLM stable
      !! form from Beljaars and Holtslag (1991)
      !!
      !! Author: L. Brodeau, june 2016 / AeroBulk
      !!         (https://sourceforge.net/p/aerobulk)
      !!----------------------------------------------------------------
      !!
      REAL(wp), DIMENSION(jpi,jpj) :: psi_h_coare
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) :: zta, zphi_h, zphi_c, zpsi_k, zpsi_c, zf, zc, zstab
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zta = pzeta(ji,jj)
            !
            zphi_h = (ABS(1. - 15.*zta))**.5  !! Kansas unstable   (zphi_h = zphi_m**2 when unstable, zphi_m when stable)
            !
            zpsi_k = 2.*LOG((1. + zphi_h)/2.)
            !
            zphi_c = (ABS(1. - 34.15*zta))**.3333   !! Convective
            !
            zpsi_c = 1.5*LOG((1. + zphi_c + zphi_c*zphi_c)/3.) &
               &    -1.7320508*ATAN((1. + 2.*zphi_c)/1.7320508) + 1.813799447
            !
            zf = zta*zta
            zf = zf/(1. + zf)
            zc = MIN(50.,0.35*zta)
            zstab = 0.5 + SIGN(0.5, zta)
            !
            psi_h_coare(ji,jj) = (1. - zstab) * ( (1. - zf)*zpsi_k + zf*zpsi_c ) &
               &                -   zstab     * ( (ABS(1. + 2.*zta/3.))**1.5     &
               &                           + .6667*(zta - 14.28)/EXP(zc) + 8.525 )
            !
         END DO
      END DO
      !
   END FUNCTION psi_h_coare


   FUNCTION visc_air( ptak )
      !!---------------------------------------------------------------------
      !! Air kinetic viscosity (m^2/s) given from temperature in degrees...
      !!
      !! Author: L. Brodeau, june 2016 / AeroBulk
      !!         (https://sourceforge.net/p/aerobulk)
      !!---------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             ::   visc_air
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   ptak       ! air temperature   [K]
      !
      INTEGER  ::   ji, jj         ! dummy loop indices
      REAL(wp) ::   ztc, ztc2      ! local scalar
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            ztc  = ptak(ji,jj) - rt0   ! air temp, in deg. C
            ztc2 = ztc*ztc
            visc_air(ji,jj) = 1.326E-5*(1. + 6.542E-3*ztc + 8.301E-6*ztc2 - 4.84E-9*ztc2*ztc)
         END DO
      END DO
      !
   END FUNCTION visc_air

   !!======================================================================
END MODULE sbcblk_algo_coare3p5
