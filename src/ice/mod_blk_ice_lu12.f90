! AeroBulk / 2020 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
! NEMO actually only needs the flux over sea-ice: F_ice
!  * F_ice is then used by SI3 which returns the basal ice-ocean flux F_b.
!  * As its SBC "Fs", OPA then uses Fs = (1.-ifr)*F_oce + ifr*F_b
!
! So we do not need the present routine to return a flux over open water
! At least when ifr is rather small...
!
!
!
! In "ice.F90" of SI3 (NEMO), the following pond information is available (also the same per ice cat.):
!  at_ip      !: total melt pond concentration
!  hm_ip      !: mean melt pond depth                     [m]
!  vt_ip      !: total melt pond volume per gridcell area [m]
! 
!
!
MODULE mod_blk_ice_lu12
   !!====================================================================================
   !!       Computes turbulent components of surface fluxes over sea-ice
   !!       Routine turb_ice_lu12 maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, January 2020
   !!
   !!====================================================================================
   USE mod_const       !: physical and other constants
   USE mod_phymbl      !: misc. physical functions

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: turb_ice_lu12, CdN10_form_LU12

   REAL(wp), PARAMETER :: rz0_s_0  = 0.69e-3_wp  ! Eq.(43) of Lupkes & Gryanik (2015) [m] => to estimate CdN10 for skin drag!
   REAL(wp), PARAMETER :: rCe_0    = 2.23E-3_wp !LOLO: this one can be more accurate when sea-ice data => Lupkes et al (2013), Eq.(1)
   REAL(wp), PARAMETER :: rNu_0    = 1._wp
   REAL(wp), PARAMETER :: rMu_0    = 1._wp
   REAL(wp), PARAMETER :: rBeta_0  = 1._wp
   
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE turb_ice_lu12( kt, zt, zu, Ts_i, t_zt, qs_i, q_zt, U_zu, frice, &
      &                      Cd, Ch, Ce, t_zu, q_zu, Ub,                      &
      &                      CdN, ChN, CeN, xz0, xu_star, xL, xUN10 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_ice_lu12  ***
      !!
      !! ** Purpose :   Computestransfert coefficients of turbulent surface
      !!                fluxes according
      !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!                Returns the effective bulk wind speed at zu to be used in the bulk formulas
      !!
      !! INPUT :
      !! -------
      !!    *  kt   : current time step (starts at 1)
      !!    *  zt   : height for temperature and spec. hum. of air            [m]
      !!    *  zu   : height for wind speed (usually 10m)                     [m]
      !!    *  Ts_i  : surface temperature of sea-ice                         [K]
      !!    *  t_zt : potential air temperature at zt                         [K]
      !!    *  qs_i  : saturation specific humidity at temp. Ts_i over ice    [kg/kg]
      !!    *  q_zt : specific humidity of air at zt                          [kg/kg]
      !!    *  U_zu : scalar wind speed at zu                                 [m/s]
      !!    * frice : sea-ice concentration        (fraction)
      !!
      !! OUTPUT :
      !! --------
      !!    *  Cd     : drag coefficient
      !!    *  Ch     : sensible heat coefficient
      !!    *  Ce     : evaporation coefficient
      !!    *  t_zu   : pot. air temperature adjusted at wind height zu       [K]
      !!    *  q_zu   : specific humidity of air        //                    [kg/kg]
      !!    *  Ub     : bulk wind speed at zu that we used                    [m/s]
      !!
      !! OPTIONAL OUTPUT:
      !! ----------------
      !!    * CdN      : neutral-stability drag coefficient
      !!    * ChN      : neutral-stability sensible heat coefficient
      !!    * CeN      : neutral-stability evaporation coefficient
      !!    * xz0     : return the aerodynamic roughness length (integration constant for wind stress) [m]
      !!    * xu_star : return u* the friction velocity                    [m/s]
      !!    * xL      : return the Obukhov length                          [m]
      !!    * xUN10   : neutral wind speed at 10m                          [m/s]
      !!
      !! ** Author: L. Brodeau, January 2020 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      INTEGER,  INTENT(in )                     :: kt    ! current time step
      REAL(wp), INTENT(in )                     :: zt    ! height for t_zt and q_zt                    [m]
      REAL(wp), INTENT(in )                     :: zu    ! height for U_zu                             [m]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: Ts_i  ! ice surface temperature                [Kelvin]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: t_zt  ! potential air temperature              [Kelvin]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: qs_i  ! sat. spec. hum. at ice/air interface    [kg/kg]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: q_zt  ! spec. air humidity at zt               [kg/kg]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: U_zu  ! relative wind module at zu                [m/s]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: frice ! sea-ice concentration        (fraction)
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Cd    ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Ch    ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Ce    ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: t_zu  ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: q_zu  ! spec. humidity adjusted at zu           [kg/kg]
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Ub ! bulk wind speed at zu                     [m/s]
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   CdN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   ChN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   CeN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xz0  ! Aerodynamic roughness length   [m]
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xu_star  ! u*, friction velocity
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xL  ! zeta (zu/L)
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xUN10  ! Neutral wind at zu
      !!
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: dt_zu, dq_zu
      !!
      LOGICAL ::  lreturn_cdn=.FALSE., lreturn_chn=.FALSE., lreturn_cen=.FALSE., &
         &        lreturn_z0=.FALSE., lreturn_ustar=.FALSE., lreturn_L=.FALSE., lreturn_UN10=.FALSE.
      !!
      CHARACTER(len=40), PARAMETER :: crtnm = 'turb_ice_lu12@mod_blk_ice_lu12.f90'
      !!----------------------------------------------------------------------------------
      !u_star(jpi,jpj), t_star(jpi,jpj), q_star(jpi,jpj),  &
      ALLOCATE ( dt_zu(jpi,jpj),  dq_zu(jpi,jpj) )
      
      IF( PRESENT(CdN) )     lreturn_cdn   = .TRUE.
      IF( PRESENT(ChN) )     lreturn_chn   = .TRUE.
      IF( PRESENT(CeN) )     lreturn_cen   = .TRUE.
      IF( PRESENT(xz0) )     lreturn_z0    = .TRUE.
      IF( PRESENT(xu_star) ) lreturn_ustar = .TRUE.
      IF( PRESENT(xL) )      lreturn_L     = .TRUE.
      IF( PRESENT(xUN10) )   lreturn_UN10  = .TRUE.
      
      !! Scalar wind speed cannot be below 0.2 m/s
      Ub = MAX( U_zu, wspd_thrshld_ice )

      !! First guess of temperature and humidity at height zu:
      t_zu = MAX( t_zt ,   100._wp )   ! who knows what's given on masked-continental regions...
      q_zu = MAX( q_zt , 0.1e-6_wp )   !               "

      !! Air-Ice differences (and we don't want it to be 0!)
      dt_zu = t_zu - Ts_i ;   dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6_wp), dt_zu )
      dq_zu = q_zu - qs_i ;   dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9_wp), dq_zu )



      !! To estimate CDN10_skin:
      !!  we use the method that comes in LG15, i.e. by starting from a default roughness length z0 for skin drag:

      Ce(:,:) = rz0_s_0 !! temporary array to contain roughness length for skin drag !
      !!
      Cd(:,:) = Cd_from_z0( zu, Ce(:,:) )  + CdN10_form_LU12( frice(:,:) )
      !!          N10 skin drag                     N10 form drag

      Ch(:,:) = Cd(:,:)
      Ce(:,:) = Cd(:,:)
      

      IF( lreturn_cdn )   CdN = Cd(:,:)
      IF( lreturn_chn )   ChN = Ch(:,:)
      IF( lreturn_cen )   CeN = Ce(:,:)

      IF( lreturn_z0 )    xz0     = z0_from_Cd( zu, Cd )
      IF( lreturn_ustar ) xu_star = SQRT(Cd)*Ub
      IF( lreturn_L )     xL      = 1./One_on_L(t_zu, q_zu, SQRT(Cd)*Ub, &
         &                          Cd/SQRT(Cd)*dt_zu, Cd/SQRT(Cd)*dq_zu)
      IF( lreturn_UN10 )  xUN10   = SQRT(Cd)*Ub/vkarmn * LOG( 10._wp / z0_from_Cd( zu, Cd ) )
      
      DEALLOCATE ( dt_zu, dq_zu )

   END SUBROUTINE turb_ice_lu12


   
   FUNCTION CdN10_form_LU12( pfrice )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  CdN10_form_LU12  ***
      !!
      !! ** Purpose :    Computes the "form" contribution of the neutral air-ice
      !!                 drag referenced at 10m to make it dependent on edges at
      !!                 leads, melt ponds and flows (to be added to the "skin"
      !!                 contribution. After some
      !!                 approximations, this can be resumed to a dependency on
      !!                 ice concentration.
      !!
      !! ** Method :     The parameterization is taken from Lupkes et al. (2012) eq.(50)
      !!                 with the highest level of approximation: level4, eq.(59)
      !!                 The generic drag over a cell partly covered by ice can be re-written as follows:
      !!
      !!                 Cd = Cdw * (1-A) + Cdi * A + Ce * (1-A)**(nu+1/(10*beta)) * A**mu
      !!
      !!                    Ce = 2.23e-3       , as suggested by Lupkes (eq. 59)
      !!                    nu = mu = beta = 1 , as suggested by Lupkes (eq. 59)
      !!                    A is the concentration of ice minus melt ponds (if any)
      !!
      !!                 This new drag has a parabolic shape (as a function of A) starting at
      !!                 Cdw(say 1.5e-3) for A=0, reaching 1.97e-3 for A~0.5
      !!                 and going down to Cdi(say 1.4e-3) for A=1
      !!
      !!                 It is theoretically applicable to all ice conditions (not only MIZ)
      !!                 => see Lupkes et al (2013)
      !!
      !! ** References : Lupkes et al. JGR 2012 (theory)
      !!                 Lupkes et al. GRL 2013 (application to GCM)
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)              :: CdN10_form_LU12  ! neutral FORM drag coefficient contribution over sea-ice
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pfrice           ! ice concentration [fraction]  => at_i_b

      !!----------------------------------------------------------------------
      REAL(wp)            ::   zcoef
      !!----------------------------------------------------------------------
      zcoef = rNu_0 + 1._wp / ( 10._wp * rBeta_0 )
      
      !! We are not an AGCM, we are an OGCM!!! => we drop term "(1 - A)*Cd_w"
      !!  => so we keep only the last rhs terms of Eq.(1) of Lupkes et al, 2013 that we divide by "A":
      
      CdN10_form_LU12(:,:) = rCe_0 * pfrice(:,:)**(rMu_0 - 1._wp) * (1._wp - pfrice(:,:))**zcoef
      
      !! => seems okay for winter 100% sea-ice as second rhs term vanishes as pfrice == 1....
      
   END FUNCTION CdN10_form_LU12


   !!======================================================================
END MODULE mod_blk_ice_lu12
