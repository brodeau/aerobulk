MODULE sbcblk_algo_ecmwf
   !!======================================================================
   !!                       ***  MODULE  sbcblk_algo_ecmwf  ***
   !! Computes turbulent components of surface fluxes
   !!         according to the method in IFS of the ECMWF model
   !!
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (2m) to zu (10m) if needed
   !!   * the effective bulk wind speed at 10m U_blk
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!    Using the bulk formulation/param. of IFS of ECMWF (cycle 31r2)
   !!         based on IFS doc (avaible online on the ECMWF's website)
   !!
   !!
   !!       Routine turb_ecmwf maintained and developed in AeroBulk
   !!                     (http://aerobulk.sourceforge.net/)
   !!
   !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
   !!----------------------------------------------------------------------
   !! History :  4.0  !  2016-02  (L.Brodeau)   Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   turb_ecmwf  : computes the bulk turbulent transfer coefficients
   !!                   adjusts t_air and q_air from zt to zu m
   !!                   returns the effective bulk wind speed at 10m
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE iom             ! I/O manager library
   USE lib_mpp         ! distribued memory computing library
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control
   USE sbcwave, ONLY   :  cdn_wave ! wave module
#if defined key_si3 || defined key_cice
   USE sbc_ice         ! Surface boundary condition: ice fields
#endif
   USE lib_fortran     ! to use key_nosignedzero

   USE sbc_oce         ! Surface boundary condition: ocean fields

   IMPLICIT NONE
   PRIVATE

   PUBLIC ::   TURB_ECMWF   ! called by sbcblk.F90

   !                   !! ECMWF own values for given constants, taken form IFS documentation...
   REAL(wp), PARAMETER ::   charn0 = 0.018    ! Charnock constant (pretty high value here !!!
   !                                          !    =>  Usually 0.011 for moderate winds)
   REAL(wp), PARAMETER ::   zi0     = 1000.   ! scale height of the atmospheric boundary layer...1
   REAL(wp), PARAMETER ::   Beta0    = 1.     ! gustiness parameter ( = 1.25 in COAREv3)
   REAL(wp), PARAMETER ::   rctv0    = 0.608  ! constant to obtain virtual temperature...
   REAL(wp), PARAMETER ::   Cp_dry = 1005.0   ! Specic heat of dry air, constant pressure      [J/K/kg]
   REAL(wp), PARAMETER ::   Cp_vap = 1860.0   ! Specic heat of water vapor, constant pressure  [J/K/kg]
   REAL(wp), PARAMETER ::   alpha_M = 0.11    ! For roughness length (smooth surface term)
   REAL(wp), PARAMETER ::   alpha_H = 0.40    ! (Chapter 3, p.34, IFS doc Cy31r1)
   REAL(wp), PARAMETER ::   alpha_Q = 0.62    !
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE TURB_ECMWF( zt, zu, sst, t_zt, ssq , q_zt , U_zu,   &
      &                   Cd, Ch, Ce , t_zu, q_zu, U_blk,         &
      &                   Cdn, Chn, Cen                           )
      !!----------------------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_ecmwf  ***
      !!
      !!            2015: L. Brodeau (brodeau@gmail.com)
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to IFS doc. (cycle 31)
      !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!
      !! ** Method : Monin Obukhov Similarity Theory
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
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
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
      INTEGER :: j_itt
      LOGICAL ::   l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      INTEGER , PARAMETER ::   nb_itt = 4       ! number of itterations
      !
      REAL(wp), DIMENSION(jpi,jpj) ::   u_star, t_star, q_star,   &
         &  dt_zu, dq_zu,    &
         &  znu_a,           & !: Nu_air, Viscosity of air
         &  Linv,            & !: 1/L (inverse of Monin Obukhov length...
         &  z0, z0t, z0q
      REAL(wp), DIMENSION(jpi,jpj) ::   func_m, func_h
      REAL(wp), DIMENSION(jpi,jpj) ::   ztmp0, ztmp1, ztmp2
      !!----------------------------------------------------------------------------------
      !
      ! Identical first gess as in COARE, with IFS parameter values though
      !
      l_zt_equal_zu = .FALSE.
      IF( ABS(zu - zt) < 0.01 )   l_zt_equal_zu = .TRUE.    ! testing "zu == zt" is risky with double precision


      !! First guess of temperature and humidity at height zu:
      t_zu = MAX( t_zt , 0.0  )   ! who knows what's given on masked-continental regions...
      q_zu = MAX( q_zt , 1.e-6)   !               "

      !! Pot. temp. difference (and we don't want it to be 0!)
      dt_zu = t_zu - sst   ;   dt_zu = SIGN( MAX(ABS(dt_zu),1.e-6), dt_zu )
      dq_zu = q_zu - ssq   ;   dq_zu = SIGN( MAX(ABS(dq_zu),1.e-9), dq_zu )

      znu_a = visc_air(t_zt) ! Air viscosity (m^2/s) at zt given from temperature in (K)

      ztmp2 = 0.5 * 0.5  ! initial guess for wind gustiness contribution
      U_blk = SQRT(U_zu*U_zu + ztmp2)

      ! z0     = 0.0001
      ztmp2   = 10000.     ! optimization: ztmp2 == 1/z0
      ztmp0   = LOG(zu*ztmp2)
      ztmp1   = LOG(10.*ztmp2)
      u_star = 0.035*U_blk*ztmp1/ztmp0       ! (u* = 0.035*Un10)

      z0     = charn0*u_star*u_star/grav + 0.11*znu_a/u_star
      z0t    = 0.1*EXP(vkarmn/(0.00115/(vkarmn/ztmp1)))   !  WARNING: 1/z0t !

      Cd     = (vkarmn/ztmp0)**2    ! first guess of Cd

      ztmp0 = vkarmn*vkarmn/LOG(zt*z0t)/Cd

      ztmp2 = Ri_bulk( zu, t_zu, dt_zu, q_zu, dq_zu, U_blk )   ! Ribu = Bulk Richardson number

      !! First estimate of zeta_u, depending on the stability, ie sign of Ribu (ztmp2):
      ztmp1 = 0.5 + SIGN( 0.5 , ztmp2 )
      func_m = ztmp0*ztmp2 ! temporary array !!
      !!             Ribu < 0                                 Ribu > 0   Beta = 1.25
      func_h = (1.-ztmp1)*(func_m/(1.+ztmp2/(-zu/(zi0*0.004*Beta0**3)))) &  ! temporary array !!! func_h == zeta_u
         &  +     ztmp1*(func_m*(1. + 27./9.*ztmp2/ztmp0))

      !! First guess M-O stability dependent scaling params.(u*,t*,q*) to estimate z0 and z/L
      ztmp0   =        vkarmn/(LOG(zu*z0t) - psi_h_ecmwf(func_h))

      u_star = U_blk*vkarmn/(LOG(zu) - LOG(z0)  - psi_m_ecmwf(func_h))
      t_star = dt_zu*ztmp0
      q_star = dq_zu*ztmp0

      ! What's need to be done if zt /= zu:
      IF( .NOT. l_zt_equal_zu ) THEN
         !
         !! First update of values at zu (or zt for wind)
         ztmp0 = psi_h_ecmwf(func_h) - psi_h_ecmwf(zt*func_h/zu)    ! zt*func_h/zu == zeta_t
         ztmp1 = log(zt/zu) + ztmp0
         t_zu = t_zt - t_star/vkarmn*ztmp1
         q_zu = q_zt - q_star/vkarmn*ztmp1
         q_zu = (0.5 + sign(0.5,q_zu))*q_zu !Makes it impossible to have negative humidity :

         dt_zu = t_zu - sst  ; dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6), dt_zu )
         dq_zu = q_zu - ssq  ; dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9), dq_zu )
         !
      ENDIF


      !! => that was same first guess as in COARE...


      !! First guess of inverse of Monin-Obukov length (1/L) :
      ztmp0 = (1. + rctv0*q_zu)  ! the factor to apply to temp. to get virt. temp...
      Linv  =  grav*vkarmn*(t_star*ztmp0 + rctv0*t_zu*q_star) / ( u_star*u_star * t_zu*ztmp0 )

      !! Functions such as  u* = U_blk*vkarmn/func_m
      ztmp1 = zu + z0
      ztmp0 = ztmp1*Linv
      func_m = LOG(ztmp1) -LOG(z0) - psi_m_ecmwf(ztmp0) + psi_m_ecmwf(z0*Linv)
      func_h = LOG(ztmp1*z0t) - psi_h_ecmwf(ztmp0) + psi_h_ecmwf(1./z0t*Linv)


      !! ITERATION BLOCK
      !! ***************

      DO j_itt = 1, nb_itt

         !! Bulk Richardson Number at z=zu (Eq. 3.25)
         ztmp0 = Ri_bulk(zu, t_zu, dt_zu, q_zu, dq_zu, U_blk)

         !! New estimate of the inverse of the Monin-Obukhon length (Linv == zeta/zu) :
         Linv = ztmp0*func_m*func_m/func_h / zu     ! From Eq. 3.23, Chap.3, p.33, IFS doc - Cy31r1

         !! Update func_m with new Linv:
         ztmp1 = zu + z0
         func_m = LOG(ztmp1) -LOG(z0) - psi_m_ecmwf(ztmp1*Linv) + psi_m_ecmwf(z0*Linv)

         !! Need to update roughness lengthes:
         u_star = U_blk*vkarmn/func_m
         ztmp2  = u_star*u_star
         ztmp1  = znu_a/u_star
         z0    = alpha_M*ztmp1 + charn0*ztmp2/grav
         z0t    = alpha_H*ztmp1                              ! eq.3.26, Chap.3, p.34, IFS doc - Cy31r1
         z0q    = alpha_Q*ztmp1

         !! Update wind at 10m taking into acount convection-related wind gustiness:
         ! Only true when unstable (L<0) => when ztmp0 < 0 => - !!!
         ztmp2 = ztmp2 * (MAX(-zi0*Linv/vkarmn,0.))**(2./3.) ! => w*^2  (combining Eq. 3.8 and 3.18, hap.3, IFS doc - Cy31r1)
         !! => equivalent using Beta=1 (gustiness parameter, 1.25 for COARE, also zi0=600 in COARE..)
         U_blk = MAX(sqrt(U_zu*U_zu + ztmp2), 0.2)              ! eq.3.17, Chap.3, p.32, IFS doc - Cy31r1
         ! => 0.2 prevents U_blk to be 0 in stable case when U_zu=0.


         !! Need to update "theta" and "q" at zu in case they are given at different heights
         !! as well the air-sea differences:
         IF( .NOT. l_zt_equal_zu ) THEN

            !! Arrays func_m and func_h are free for a while so using them as temporary arrays...
            func_h = psi_h_ecmwf((zu+z0)*Linv) ! temporary array !!!
            func_m = psi_h_ecmwf((zt+z0)*Linv) ! temporary array !!!

            ztmp2  = psi_h_ecmwf(z0t*Linv)
            ztmp0  = func_h - ztmp2
            ztmp1  = vkarmn/(LOG(zu+z0) - LOG(z0t) - ztmp0)
            t_star = dt_zu*ztmp1
            ztmp2  = ztmp0 - func_m + ztmp2
            ztmp1  = LOG(zt/zu) + ztmp2
            t_zu   = t_zt - t_star/vkarmn*ztmp1

            ztmp2  = psi_h_ecmwf(z0q*Linv)
            ztmp0  = func_h - ztmp2
            ztmp1  = vkarmn/(LOG(zu+z0) - LOG(z0q) - ztmp0)
            q_star = dq_zu*ztmp1
            ztmp2  = ztmp0 - func_m + ztmp2
            ztmp1  = log(zt/zu) + ztmp2
            q_zu   = q_zt - q_star/vkarmn*ztmp1

            dt_zu = t_zu - sst ;  dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6), dt_zu )
            dq_zu = q_zu - ssq ;  dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9), dq_zu )

         END IF

         !! Updating because of updated z0 and z0t and new Linv...
         ztmp1 = zu + z0
         ztmp0 = ztmp1*Linv
         func_m = log(ztmp1) - LOG(z0 ) - psi_m_ecmwf(ztmp0) + psi_m_ecmwf(z0 *Linv)
         func_h = log(ztmp1) - LOG(z0t) - psi_h_ecmwf(ztmp0) + psi_h_ecmwf(z0t*Linv)

      END DO

      Cd = vkarmn*vkarmn/(func_m*func_m)
      Ch = vkarmn*vkarmn/(func_m*func_h)
      ztmp1 = log((zu + z0)/z0q) - psi_h_ecmwf((zu + z0)*Linv) + psi_h_ecmwf(z0q*Linv)   ! func_q
      Ce = vkarmn*vkarmn/(func_m*ztmp1)

      ztmp1 = zu + z0
      Cdn = vkarmn*vkarmn / (log(ztmp1/z0 )*log(ztmp1/z0 ))
      Chn = vkarmn*vkarmn / (log(ztmp1/z0t)*log(ztmp1/z0t))
      Cen = vkarmn*vkarmn / (log(ztmp1/z0q)*log(ztmp1/z0q))

   END SUBROUTINE TURB_ECMWF


   FUNCTION psi_m_ecmwf( pzeta )
      !!----------------------------------------------------------------------------------
      !! Universal profile stability function for momentum
      !!     ECMWF / as in IFS cy31r1 documentation, available online
      !!     at ecmwf.int
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: psi_m_ecmwf
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) :: zzeta, zx, ztmp, psi_unst, psi_stab, stab
      !!----------------------------------------------------------------------------------
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zzeta = MIN( pzeta(ji,jj) , 5. ) !! Very stable conditions (L positif and big!):
            !
            ! Unstable (Paulson 1970):
            !   eq.3.20, Chap.3, p.33, IFS doc - Cy31r1
            zx = SQRT(ABS(1. - 16.*zzeta))
            ztmp = 1. + SQRT(zx)
            ztmp = ztmp*ztmp
            psi_unst = LOG( 0.125*ztmp*(1. + zx) )   &
               &       -2.*ATAN( SQRT(zx) ) + 0.5*rpi
            !
            ! Unstable:
            ! eq.3.22, Chap.3, p.33, IFS doc - Cy31r1
            psi_stab = -2./3.*(zzeta - 5./0.35)*EXP(-0.35*zzeta) &
               &       - zzeta - 2./3.*5./0.35
            !
            ! Combining:
            stab = 0.5 + SIGN(0.5, zzeta) ! zzeta > 0 => stab = 1
            !
            psi_m_ecmwf(ji,jj) = (1. - stab) * psi_unst & ! (zzeta < 0) Unstable
               &                +      stab  * psi_stab   ! (zzeta > 0) Stable
            !
         END DO
      END DO
      !
   END FUNCTION psi_m_ecmwf

   
   FUNCTION psi_h_ecmwf( pzeta )
      !!----------------------------------------------------------------------------------
      !! Universal profile stability function for temperature and humidity
      !!     ECMWF / as in IFS cy31r1 documentation, available online
      !!     at ecmwf.int
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: psi_h_ecmwf
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) ::  zzeta, zx, psi_unst, psi_stab, stab
      !!----------------------------------------------------------------------------------
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zzeta = MIN(pzeta(ji,jj) , 5.)   ! Very stable conditions (L positif and big!):
            !
            zx  = ABS(1. - 16.*zzeta)**.25        ! this is actually (1/phi_m)**2  !!!
            !                                     ! eq.3.19, Chap.3, p.33, IFS doc - Cy31r1
            ! Unstable (Paulson 1970) :
            psi_unst = 2.*LOG(0.5*(1. + zx*zx))   ! eq.3.20, Chap.3, p.33, IFS doc - Cy31r1
            !
            ! Stable:
            psi_stab = -2./3.*(zzeta - 5./0.35)*EXP(-0.35*zzeta) & ! eq.3.22, Chap.3, p.33, IFS doc - Cy31r1
               &       - ABS(1. + 2./3.*zzeta)**1.5 - 2./3.*5./0.35 + 1. 
            ! LB: added ABS() to avoid NaN values when unstable, which contaminates the unstable solution...
            !
            stab = 0.5 + SIGN(0.5, zzeta) ! zzeta > 0 => stab = 1
            !
            !
            psi_h_ecmwf(ji,jj) = (1. - stab) * psi_unst &   ! (zzeta < 0) Unstable
               &                +    stab    * psi_stab     ! (zzeta > 0) Stable
            !
         END DO
      END DO
      !
   END FUNCTION psi_h_ecmwf


   FUNCTION Ri_bulk( pz, ptz, pdt, pqz, pdq, pub )
      !!----------------------------------------------------------------------------------
      !! Bulk Richardson number (Eq. 3.25 IFS doc)
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) ::   Ri_bulk   !
      !
      REAL(wp)                    , INTENT(in) ::   pz    ! height above the sea        [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   ptz   ! air temperature at pz m     [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pdt   ! ptz - sst                   [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pqz   ! air temperature at pz m [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pdq   ! pqz - ssq               [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pub   ! bulk wind speed           [m/s]
      !!----------------------------------------------------------------------------------
      !
      Ri_bulk =   grav*pz/(pub*pub)                                          &
         &      * ( pdt/(ptz - 0.5_wp*(pdt + grav*pz/(Cp_dry+Cp_vap*pqz)))   &
         &          + rctv0*pdq )
      !
   END FUNCTION Ri_bulk


   FUNCTION visc_air(ptak)
      !!----------------------------------------------------------------------------------
      !! Air kinetic viscosity (m^2/s) given from temperature in degrees...
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             ::   visc_air   !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   ptak       ! air temperature in (K)
      !
      INTEGER  ::   ji, jj      ! dummy loop indices
      REAL(wp) ::   ztc, ztc2   ! local scalar
      !!----------------------------------------------------------------------------------
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            ztc  = ptak(ji,jj) - rt0   ! air temp, in deg. C
            ztc2 = ztc*ztc
            visc_air(ji,jj) = 1.326e-5*(1. + 6.542E-3*ztc + 8.301e-6*ztc2 - 4.84e-9*ztc2*ztc)
         END DO
      END DO
      !
   END FUNCTION visc_air

   !!======================================================================
END MODULE sbcblk_algo_ecmwf
