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
! # "h_f" in Lupkes is mean floe freeboard
!
! "ice freeboard" is saved as:
!    ( zrho1 * hm_i(:,:) - zrho2 * hm_s(:,:) )
!    WHERE( z2d < 0._wp )   z2d = 0._wp
! into "icewri.F90"
!
!   possibly "h_f" in Lupkes*
!
MODULE mod_blk_ice_lu12
   !!====================================================================================
   !!       Computes turbulent components of surface fluxes over sea-ice
   !!       Routine turb_ice_lu12 maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!  Lüpkes, C., Gryanik, V. M., Hartmann, J., and Andreas, E. L. ( 2012), A parametrization, based on sea ice morphology,
   !!  of the neutral atmospheric drag coefficients for weather prediction and climate models, J. Geophys. Res., 117, D13112,
   !!  doi:10.1029/2012JD017630.
   !!
   !!       => Sespite the fact that the sea-ice concentration (frice) must be provided,
   !!          only transfer coefficients, and air temp. + hum. height adjustement
   !!          over ice are returned/performed.
   !!        ==> 'frice' is only here to estimate the form drag caused by sea-ice...
   !!
   !!       Routine turb_ice_lu12 maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, January 2020
   !!
   !!====================================================================================
   USE mod_const       !: physical and other constants
   USE mod_phymbl      !: misc. physical functions
   !USE mod_blk_ncar, ONLY: CD_N10_NCAR  !: in order to have a decent estimate of z0 over water
   USE mod_cdn_form_ice

   IMPLICIT NONE
   PRIVATE

   REAL(wp), PARAMETER :: rz0_i_s_0  = 0.69e-3_wp  ! Eq.(43) of Lupkes & Gryanik (2015) [m] => to estimate CdN10 for skin drag!
   REAL(wp), PARAMETER :: rz0_i_f_0  = 4.54e-4_wp  ! bottom p.562 MIZ [m] (LG15)   
   REAL(wp), PARAMETER :: rz0_w_s_0 = 3.27E-4      ! fixed roughness length over water (paragraph below Eq.36)

   PUBLIC :: turb_ice_lu12

   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE turb_ice_lu12( kt, zt, zu, Ts_i, t_zt, qs_i, q_zt, U_zu, frice, &
                                !&                      hf, Di,                                          &
      &                      Cd_i, Ch_i, Ce_i, t_zu_i, q_zu_i, Ubzu,            &
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
      !! OPTIONAL INPUT:
      !! ----------------
      !!    * hf    : mean freeboard of floes    [m]
      !!    * Di    : cross wind dimension of the floe (aka effective edge length for form drag)   [m]
      !!
      !! OUTPUT :
      !! --------
      !!    *  Cd_i   : drag coefficient over sea-ice
      !!    *  Ch_i   : sensible heat coefficient over sea-ice
      !!    *  Ce_i   : sublimation coefficient over sea-ice
      !!    *  t_zu_i   : pot. air temperature adjusted at wind height zu       [K]
      !!    *  q_zu_i   : specific humidity of air        //                    [kg/kg]
      !!    *  Ubzu     : bulk wind speed at zu that we used                    [m/s]
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
      !!
      !REAL(wp), INTENT(in ), DIMENSION(jpi,jpj), OPTIONAL :: hf    ! mean freeboard of floes    [m]
      !REAL(wp), INTENT(in ), DIMENSION(jpi,jpj), OPTIONAL :: Di    ! cross wind dimension of the floe (aka effective edge length for form drag)   [m]
      !!
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Cd_i  ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Ch_i  ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Ce_i  ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: t_zu_i  ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: q_zu_i  ! spec. humidity adjusted at zu           [kg/kg]
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Ubzu ! bulk wind speed at zu                     [m/s]
      !!
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj), OPTIONAL :: CdN
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj), OPTIONAL :: ChN
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj), OPTIONAL :: CeN
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj), OPTIONAL :: xz0  ! Aerodynamic roughness length   [m]
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj), OPTIONAL :: xu_star  ! u*, friction velocity
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj), OPTIONAL :: xL  ! zeta (zu/L)
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj), OPTIONAL :: xUN10  ! Neutral wind at zu
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: dt_zu, dq_zu, z0_w, z0_i
      !!
      LOGICAL :: l_known_hf=.FALSE., l_known_Di=.FALSE.
      LOGICAL :: lreturn_cdfrm=.FALSE., lreturn_cdn=.FALSE., lreturn_chn=.FALSE., lreturn_cen=.FALSE., &
         &       lreturn_z0=.FALSE., lreturn_ustar=.FALSE., lreturn_L=.FALSE., lreturn_UN10=.FALSE.
      !!
      CHARACTER(len=40), PARAMETER :: crtnm = 'turb_ice_lu12@mod_blk_ice_lu12.f90'
      !!----------------------------------------------------------------------------------
      !l_known_hf    = PRESENT(hf)
      !l_known_Di    = PRESENT(Di)
      !lreturn_cdfrm = PRESENT(CdN_frm)
      lreturn_cdn   = PRESENT(CdN)
      lreturn_chn   = PRESENT(ChN)
      lreturn_cen   = PRESENT(CeN)
      lreturn_z0    = PRESENT(xz0)
      lreturn_ustar = PRESENT(xu_star)
      lreturn_L     = PRESENT(xL)
      lreturn_UN10  = PRESENT(xUN10)


      !u_star(jpi,jpj), t_star(jpi,jpj), q_star(jpi,jpj),  &
      ALLOCATE ( dt_zu(jpi,jpj), dq_zu(jpi,jpj) , z0_w(jpi,jpj), z0_i(jpi,jpj) )

      !! Scalar wind speed cannot be below 0.2 m/s
      Ubzu = MAX( U_zu, wspd_thrshld_ice )

      !! First guess of temperature and humidity at height zu:
      t_zu_i = MAX( t_zt ,   100._wp )   ! who knows what's given on masked-continental regions...
      q_zu_i = MAX( q_zt , 0.1e-6_wp )   !               "

      !! Air-Ice differences (and we don't want it to be 0!)
      dt_zu = t_zu_i - Ts_i ;   dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6_wp), dt_zu )
      dq_zu = q_zu_i - qs_i ;   dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9_wp), dq_zu )



      !! To estimate CDN10_skin:
      !!  we use the method that comes in LG15, i.e. by starting from a default roughness length z0 for skin drag:

      Ce_i(:,:) = rz0_i_s_0 !! temporary array to contain roughness length for skin drag !


      !! Method #1:
      !Cd_i(:,:) = Cd_from_z0( zu, Ce_i(:,:) )  + CdN10_f_LU13( frice(:,:) )
      !IF( lreturn_cdfrm ) CdN_frm = CdN10_f_LU13( frice(:,:) )
      !PRINT *, 'LOLO: estimate of Cd_f_i method #1 =>', CdN10_f_LU13( frice(:,:) ); PRINT *, ''

      !! Method #2:
      !! We need an estimate of z0 over water:
      !z0_w(:,:) = z0_from_Cd( zu, CD_N10_NCAR(Ubzu) )
      !!PRINT *, 'LOLO: estimate of z0_w =>', z0_w
      !Cd_i(:,:)   = Cd_from_z0( zu, Ce_i(:,:) )  + CdN10_f_LU12( frice(:,:), z0_w(:,:) )
      !IF( lreturn_cdfrm ) CdN_frm =  CdN10_f_LU12( frice(:,:), z0_w(:,:) )
      !!          N10 skin drag                     N10 form drag

      !! Method #3:
      !Cd_i(:,:)   = Cd_from_z0( zu, Ce_i(:,:) ) + CdN10_f_LU12_eq36( frice(:,:) )
      !IF( lreturn_cdfrm ) CdN_frm = CdN10_f_LU12_eq36( frice(:,:) )
      !PRINT *, 'LOLO: estimate of Cd_f_i method #2 =>', CdN10_f_LU12( frice(:,:), z0_w(:,:) )

      !! Method #4:
      !! using eq.21 of LG15 instead:
      z0_i(:,:) = rz0_i_f_0
      !Cd_i(:,:)   = Cd_from_z0( zu, Ce_i(:,:) )  + CdN_f_LG15( zu, frice(:,:), z0_i(:,:) ) / frice(:,:)
      Cd_i(:,:)   = Cd_from_z0( zu, Ce_i(:,:) )  + CdN_f_LG15( zu, frice(:,:), z0_i(:,:) ) !/ frice(:,:)
      !IF( lreturn_cdfrm ) CdN_frm = CdN_f_LG15( zu, frice(:,:), z0_i(:,:) )


      Ch_i(:,:) = Cd_i(:,:)
      Ce_i(:,:) = Cd_i(:,:)

      IF( lreturn_cdn )   CdN = Cd_i(:,:)
      IF( lreturn_chn )   ChN = Ch_i(:,:)
      IF( lreturn_cen )   CeN = Ce_i(:,:)

      IF( lreturn_z0 )    xz0     = z0_from_Cd( zu, Cd_i )
      IF( lreturn_ustar ) xu_star = SQRT(Cd_i)*Ubzu
      IF( lreturn_L )     xL      = 1./One_on_L(t_zu_i, q_zu_i, SQRT(Cd_i)*Ubzu, &
         &                          Cd_i/SQRT(Cd_i)*dt_zu, Cd_i/SQRT(Cd_i)*dq_zu)
      IF( lreturn_UN10 )  xUN10   = SQRT(Cd_i)*Ubzu/vkarmn * LOG( 10._wp / z0_from_Cd( zu, Cd_i ) )

      DEALLOCATE ( dt_zu, dq_zu, z0_w )

   END SUBROUTINE turb_ice_lu12

   !!======================================================================
END MODULE mod_blk_ice_lu12
