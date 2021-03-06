! AeroBulk / 2019 / L. Brodeau
!
!***********************************************************************************
! MODULE that gathers a collection of usefull functions related to the physics /
! thermodynamics of air within the Marine Boundaty Layer
!***********************************************************************************
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following paper:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.


MODULE mod_phymbl

   !!   pot_temp      : potential temperature
   !!   virt_temp     : virtual (aka sensible) temperature (potential or absolute)
   !!   rho_air       : density of (moist) air (depends on T_air, q_air and SLP
   !!   visc_air      : kinematic viscosity (aka Nu_air) of air from temperature
   !!   L_vap         : latent heat of vaporization of water as a function of temperature
   !!   cp_air        : specific heat of (moist) air (depends spec. hum. q_air)
   !!   gamma_moist   : adiabatic lapse-rate of moist air
   !!   One_on_L      : 1. / ( Obukhov length )
   !!   Ri_bulk       : bulk Richardson number aka BRN
   !!   q_sat         : saturation humidity as a function of SLP and temperature
   !!   q_air_rh      : specific humidity as a function of RH (fraction, not %), t_air and SLP

   USE mod_const

   IMPLICIT NONE

   INTERFACE pot_temp
      MODULE PROCEDURE pot_temp_vctr, pot_temp_sclr
   END INTERFACE pot_temp

   INTERFACE abs_temp
      MODULE PROCEDURE abs_temp_vctr, abs_temp_sclr
   END INTERFACE abs_temp

   INTERFACE virt_temp
      MODULE PROCEDURE virt_temp_vctr, virt_temp_sclr
   END INTERFACE virt_temp

   INTERFACE Pz_from_P0_tz_qz
      MODULE PROCEDURE Pz_from_P0_tz_qz_vctr, Pz_from_P0_tz_qz_sclr
   END INTERFACE Pz_from_P0_tz_qz

   INTERFACE Theta_from_z_P0_T_q
      MODULE PROCEDURE Theta_from_z_P0_T_q_vctr, Theta_from_z_P0_T_q_sclr
   END INTERFACE Theta_from_z_P0_T_q

   INTERFACE T_from_z_P0_Theta_q
      MODULE PROCEDURE T_from_z_P0_Theta_q_vctr, T_from_z_P0_Theta_q_sclr
   END INTERFACE T_from_z_P0_Theta_q

   INTERFACE visc_air
      MODULE PROCEDURE visc_air_vctr, visc_air_sclr
   END INTERFACE visc_air

   INTERFACE gamma_moist
      MODULE PROCEDURE gamma_moist_vctr, gamma_moist_sclr
   END INTERFACE gamma_moist

   INTERFACE e_sat
      MODULE PROCEDURE e_sat_vctr, e_sat_sclr
   END INTERFACE e_sat

   INTERFACE e_sat_ice
      MODULE PROCEDURE e_sat_ice_vctr, e_sat_ice_sclr
   END INTERFACE e_sat_ice

   INTERFACE de_sat_dt_ice
      MODULE PROCEDURE de_sat_dt_ice_vctr, de_sat_dt_ice_sclr
   END INTERFACE de_sat_dt_ice

   INTERFACE One_on_L
      MODULE PROCEDURE One_on_L_vctr, One_on_L_sclr
   END INTERFACE One_on_L

   INTERFACE Ri_bulk
      MODULE PROCEDURE Ri_bulk_vctr, Ri_bulk_sclr
   END INTERFACE Ri_bulk

   INTERFACE q_sat
      MODULE PROCEDURE q_sat_vctr, q_sat_sclr
   END INTERFACE q_sat

   INTERFACE dq_sat_dt_ice
      MODULE PROCEDURE dq_sat_dt_ice_vctr, dq_sat_dt_ice_sclr
   END INTERFACE dq_sat_dt_ice

   INTERFACE L_vap
      MODULE PROCEDURE L_vap_vctr, L_vap_sclr
   END INTERFACE L_vap

   INTERFACE rho_air
      MODULE PROCEDURE rho_air_vctr, rho_air_sclr
   END INTERFACE rho_air

   INTERFACE cp_air
      MODULE PROCEDURE cp_air_vctr, cp_air_sclr
   END INTERFACE cp_air

   INTERFACE alpha_sw
      MODULE PROCEDURE alpha_sw_vctr, alpha_sw_sclr
   END INTERFACE alpha_sw

   INTERFACE update_qnsol_tau
      MODULE PROCEDURE update_qnsol_tau_vctr, update_qnsol_tau_sclr
   END INTERFACE update_qnsol_tau

   INTERFACE bulk_formula
      MODULE PROCEDURE bulk_formula_vctr, bulk_formula_sclr
   END INTERFACE bulk_formula

   INTERFACE qlw_net
      MODULE PROCEDURE qlw_net_vctr, qlw_net_sclr
   END INTERFACE qlw_net

   INTERFACE z0_from_Cd
      MODULE PROCEDURE z0_from_Cd_vctr, z0_from_Cd_sclr
   END INTERFACE z0_from_Cd

   INTERFACE z0_from_ustar
      MODULE PROCEDURE z0_from_ustar_vctr, z0_from_ustar_sclr
   END INTERFACE z0_from_ustar

   INTERFACE UN10_from_CD
      MODULE PROCEDURE UN10_from_CD_vctr, UN10_from_CD_sclr
   END INTERFACE UN10_from_CD

   INTERFACE f_m_louis
      MODULE PROCEDURE f_m_louis_vctr, f_m_louis_sclr
   END INTERFACE f_m_louis

   INTERFACE f_h_louis
      MODULE PROCEDURE f_h_louis_vctr, f_h_louis_sclr
   END INTERFACE f_h_louis


   !PUBLIC virt_temp
   !PUBLIC rho_air
   !PUBLIC visc_air
   !PUBLIC L_vap
   !PUBLIC cp_air
   !!PUBLIC gamma_moist
   !PUBLIC One_on_L
   !PUBLIC Ri_bulk
   !PUBLIC q_sat
   !PUBLIC dq_sat_dt_ice
   !PUBLIC e_sat
   !PUBLIC e_sat_ice
   !PUBLIC de_sat_dt_ice
   !PUBLIC e_air
   !PUBLIC rh_air
   !PUBLIC rho_air_adv
   !PUBLIC dry_static_energy
   !PUBLIC q_air_rh
   !PUBLIC q_air_dp
   !PUBLIC q_sat_crude
   !PUBLIC update_qnsol_tau
   !PUBLIC alpha_sw
   !PUBLIC bulk_formula
   !PUBLIC qlw_net
   !PUBLIC z0_from_Cd, z0_from_ustar
   !PUBLIC Cd_from_z0
   !PUBLIC f_m_louis, f_h_louis
   !PUBLIC UN10_from_ustar
   !PUBLIC UN10_from_CDN
   !PUBLIC UN10_from_CD
   !PUBLIC Re_rough_tq_LKB
   !PUBLIC z0tq_LKB
   !PUBLIC vmean, variance
   !PUBLIC to_kelvin_3d


   !! P R I V A T E :

   REAL(wp), PARAMETER, PRIVATE :: &
                                !! Constants for Goff formula in the presence of ice:
      &      rAg_i = -9.09718_wp, &
      &      rBg_i = -3.56654_wp, &
      &      rCg_i = 0.876793_wp, &
      &      rDg_i = LOG10(6.1071_wp)

   REAL(wp), PARAMETER, PRIVATE :: rc_louis  = 5._wp
   REAL(wp), PARAMETER, PRIVATE :: rc2_louis = rc_louis * rc_louis
   REAL(wp), PARAMETER, PRIVATE :: ram_louis = 2. * rc_louis
   REAL(wp), PARAMETER, PRIVATE :: rah_louis = 3. * rc_louis

   REAL(wp), PARAMETER, PRIVATE :: repsilon = 1.e-6_wp




CONTAINS

   !===============================================================================================
   FUNCTION pot_temp_sclr( pTa, pPz,  pPref )
      !!------------------------------------------------------------------------
      !!
      !! Poisson equation to obtain potential temperature from absolute temperature, pressure, and
      !! the reference (sea-level) pressure.
      !!
      !! Air parcel is at height `z` m above sea-level where the pressure is `ppz`,
      !! its absolute temperature is `pTa`.
      !! `pPref` is the reference pressure at sea level aka P0.
      !!
      !! Author: L. Brodeau, June 2021 / AeroBulk
      !!         (https://github.com/brodeau/aerobulk/)
      !!------------------------------------------------------------------------
      REAL(wp), INTENT(in)           :: pTa            !: absolute air temperature at `z` m above sea level  [K]
      REAL(wp), INTENT(in)           :: pPz            !: pressure at `z` m above sea level                  [Pa]
      REAL(wp), INTENT(in), OPTIONAL :: pPref          !: reference pressure (sea-level)                     [Pa]
      REAL(wp)                       :: pot_temp_sclr  !: potential air temperature at `z` m above sea level [K]
      !!
      REAL(wp) :: zPref = Patm
      !!-------------------------------------------------------------------
      IF( PRESENT(pPref) ) zPref = pPref
      pot_temp_sclr = pTa * ( zPref / pPz )**rpoiss_dry
      !!
   END FUNCTION pot_temp_sclr

   FUNCTION pot_temp_vctr( pTa, pPz,  pPref )
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)           :: pTa
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)           :: pPz
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in), OPTIONAL :: pPref
      REAL(wp), DIMENSION(jpi,jpj)                       :: pot_temp_vctr
      !!-------------------------------------------------------------------
      IF( PRESENT(pPref) ) THEN
         pot_temp_vctr = pTa * ( pPref / pPz )**rpoiss_dry
      ELSE
         pot_temp_vctr = pTa * ( Patm  / pPz )**rpoiss_dry
      END IF
   END FUNCTION pot_temp_vctr
   !===============================================================================================


   !===============================================================================================
   FUNCTION abs_temp_sclr( pThta, pPz,  pPref )
      !!------------------------------------------------------------------------
      !!
      !! Poisson equation to obtain absolute temperature from potential temperature, pressure, and
      !! the reference (sea-level) pressure.
      !!
      !! Air parcel is at height `z` m above sea-level where the pressure is `ppz`,
      !! its potential temperature is `pThta`.
      !! `pPref` is the reference pressure at sea level aka P0.
      !!
      !! Author: L. Brodeau, June 2021 / AeroBulk
      !!         (https://github.com/brodeau/aerobulk/)
      !!------------------------------------------------------------------------
      REAL(wp), INTENT(in)           :: pThta          !: potential air temperature at `z` m above sea level  [K]
      REAL(wp), INTENT(in)           :: pPz            !: pressure at `z` m above sea level                  [Pa]
      REAL(wp), INTENT(in), OPTIONAL :: pPref          !: reference pressure (sea-level)                     [Pa]
      REAL(wp)                       :: abs_temp_sclr  !: potential air temperature at `z` m above sea level [K]
      !!
      REAL(wp) :: zPref = Patm
      !!-------------------------------------------------------------------
      IF( PRESENT(pPref) ) zPref = pPref
      abs_temp_sclr = pThta / MAX ( ( zPref / pPz )**rpoiss_dry, 1.E-9_wp )
      !!
   END FUNCTION abs_temp_sclr

   FUNCTION abs_temp_vctr( pThta, pPz,  pPref )
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)           :: pThta
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)           :: pPz
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in), OPTIONAL :: pPref
      REAL(wp), DIMENSION(jpi,jpj)                       :: abs_temp_vctr
      !!-------------------------------------------------------------------
      IF( PRESENT(pPref) ) THEN
         abs_temp_vctr = pThta / MAX ( ( pPref / pPz )**rpoiss_dry, 1.E-9_wp )
      ELSE
         abs_temp_vctr = pThta / MAX ( ( Patm  / pPz )**rpoiss_dry, 1.E-9_wp )
      END IF
   END FUNCTION abs_temp_vctr
   !===============================================================================================





   !===============================================================================================
   FUNCTION virt_temp_sclr( pTa, pqa )
      !!------------------------------------------------------------------------
      !!
      !! Compute the (absolute/potential) VIRTUAL temperature, based on the
      !! (absolute/potential) temperature and specific humidity
      !!
      !! If input temperature is absolute then output virtual temperature is absolute
      !! If input temperature is potential then output virtual temperature is potential
      !!
      !! Author: L. Brodeau, June 2019 / AeroBulk
      !!         (https://github.com/brodeau/aerobulk/)
      !!------------------------------------------------------------------------
      REAL(wp)             :: virt_temp_sclr !: virtual temperature [K]
      REAL(wp), INTENT(in) :: pTa       !: absolute or potential air temperature [K]
      REAL(wp), INTENT(in) :: pqa       !: specific humidity of air   [kg/kg]
      !!-------------------------------------------------------------------
      virt_temp_sclr = pTa * (1._wp + rctv0*pqa)
      !!
      !! This is exactly the same thing as:
      !! virt_temp_sclr = pTa * ( pwa + reps0) / (reps0*(1.+pwa))
      !! with wpa (mixing ration) defined as : pwa = pqa/(1.-pqa)
      !!
   END FUNCTION virt_temp_sclr
   !!
   FUNCTION virt_temp_vctr( pTa, pqa )
      REAL(wp), DIMENSION(jpi,jpj)             :: virt_temp_vctr !: virtual temperature [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pTa !: absolute or potential air temperature [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pqa !: specific humidity of air   [kg/kg]
      virt_temp_vctr(:,:) = pTa(:,:) * (1._wp + rctv0*pqa(:,:))
   END FUNCTION virt_temp_vctr
   !===============================================================================================




   !===============================================================================================
   FUNCTION Pz_from_P0_tz_qz_sclr( pz, pslp, pTa, pqa,  l_ice )
      !!-------------------------------------------------------------------------------
      !!                           ***  FUNCTION Pz_from_P0_tz_qz  ***
      !!
      !! ** Purpose : compute air pressure at height `pz` m above sea level,
      !!              based barometric equation, from absolute air temperature
      !!              and specific humidity at height `pz` and sea-level pressure
      !!
      !! ** Author: G. Samson,  Feb  2021
      !!            L. Brodeau, June 2021
      !!-------------------------------------------------------------------------------
      REAL(wp), INTENT(in)           :: pz               ! height above sea-level             [m]
      REAL(wp), INTENT(in)           :: pslp             ! pressure at sea-level (z=0)        [Pa]
      REAL(wp), INTENT(in)           :: pTa              ! air abs. temperature at z=`pz`     [K]
      REAL(wp), INTENT(in)           :: pqa              ! air specific humidity at z=`pz`    [kg/kg]
      LOGICAL , INTENT(in), OPTIONAL :: l_ice            ! sea-ice presence
      REAL(wp)                       :: Pz_from_P0_tz_qz_sclr  ! pressure at `pz` m above sea level [Pa]
      !!
      REAL(wp)           :: ztpa, zpa, zxm, zmask, zqsat, zf
      INTEGER            :: it
      LOGICAL            :: lice = .FALSE. ! sea-ice presence
      INTEGER, PARAMETER :: niter = 3      ! iteration indice and number
      !!-------------------------------------------------------------------------------
      IF( PRESENT(l_ice) ) lice = l_ice
      !!
      zpa = pslp              ! first guess of air pressure at zt   [Pa]
      DO it = 1, niter
         zqsat = q_sat( pTa, zpa, l_ice=lice )                               ! saturation specific humidity [kg/kg]
         zf    = pqa/zqsat
         zxm   = (1._wp - zf) * rmm_dryair + zf * rmm_water    ! moist air molar mass [kg/mol]
         zpa   = pslp * EXP( -grav * zxm * pz / ( R_gas * pTa ) )
      END DO
      !!
      Pz_from_P0_tz_qz_sclr = zpa
      !!
   END FUNCTION Pz_from_P0_tz_qz_sclr

   FUNCTION Pz_from_P0_tz_qz_vctr( pz, pslp, pTa, pqa,  l_ice )
      REAL(wp),                     INTENT(in) :: pz               ! height above sea-level             [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pslp             ! pressure at sea-level (z=0)        [Pa]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pTa              ! air abs. temperature at z=`pz`     [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pqa              ! air specific humidity at z=`pz`    [kg/kg]
      LOGICAL , OPTIONAL          , INTENT(in) :: l_ice            ! sea-ice presence
      REAL(wp), DIMENSION(jpi,jpj)             :: Pz_from_P0_tz_qz_vctr  ! pressure at `pz` m above sea level [Pa]
      !!
      INTEGER :: ji, jj         ! loop indices
      LOGICAL :: lice = .FALSE. ! presence of ice ?
      !!
      IF( PRESENT(l_ice) ) lice = l_ice
      DO jj = 1, jpj
         DO ji = 1, jpi
            Pz_from_P0_tz_qz_vctr(ji,jj) = Pz_from_P0_tz_qz_sclr( pz, pslp(ji,jj), pTa(ji,jj), pqa(ji,jj),  l_ice=lice )
         END DO
      END DO
   END FUNCTION Pz_from_P0_tz_qz_vctr

   !===============================================================================================


   !===============================================================================================
   FUNCTION Theta_from_z_P0_T_q_sclr( pz, pslp, pTa, pqa )
      !!-------------------------------------------------------------------------------
      !!                           ***  FUNCTION Theta_from_z_P0_T_q  ***
      !!
      !! ** Purpose : Converts absolute temperature at height `pz` to potential temperature,
      !!              using sea-level pressure and specific humidity at heaight `pz`
      !!
      !! ** Author: G. Samson,  Feb  2021
      !!            L. Brodeau, June 2021
      !!-------------------------------------------------------------------------------
      REAL(wp), INTENT(in) :: pz               ! height above sea-level             [m]
      REAL(wp), INTENT(in) :: pslp             ! pressure at sea-level (z=0)        [Pa]
      REAL(wp), INTENT(in) :: pTa              ! air abs. temperature at z=`pz`     [K]
      REAL(wp), INTENT(in) :: pqa              ! air specific humidity at z=`pz`    [kg/kg]
      REAL(wp)             :: Theta_from_z_P0_T_q_sclr ! air pot. temperature at z=`pz` [K]
      !!
      REAL(wp)           :: zPz
      !!-------------------------------------------------------------------------------
      zPz = Pz_from_P0_tz_qz_sclr( pz, pslp, pTa, pqa ) ! pressure at z=`pz`
      !!
      Theta_from_z_P0_T_q_sclr = pot_temp_sclr( pTa, zPz,  pPref=pslp )
      !!
   END FUNCTION Theta_from_z_P0_T_q_sclr

   FUNCTION Theta_from_z_P0_T_q_vctr( pz, pslp, pTa, pqa )
      REAL(wp),                     INTENT(in) :: pz
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pslp
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pTa
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pqa
      REAL(wp), DIMENSION(jpi,jpj)             :: Theta_from_z_P0_T_q_vctr
      !!-------------------------------------------------------------------------------
      Theta_from_z_P0_T_q_vctr = pot_temp_vctr( pTa, Pz_from_P0_tz_qz_vctr( pz, pslp, pTa, pqa ),  pPref=pslp )
   END FUNCTION Theta_from_z_P0_T_q_vctr
   !===============================================================================================


   !===============================================================================================
   FUNCTION T_from_z_P0_Theta_q_sclr( pz, pslp, pThta, pqa )
      !!-------------------------------------------------------------------------------
      !!                           ***  FUNCTION T_from_z_P0_Theta_q  ***
      !!
      !! ** Purpose : Converts potential temperature at height `pz` to absolute temperature,
      !!              using sea-level pressure and specific humidity at heaight `pz`
      !!
      !! ** Author: G. Samson,  Feb  2021
      !!            L. Brodeau, June 2021
      !!-------------------------------------------------------------------------------
      REAL(wp), INTENT(in) :: pz               ! height above sea-level             [m]
      REAL(wp), INTENT(in) :: pslp             ! pressure at sea-level (z=0)        [Pa]
      REAL(wp), INTENT(in) :: pThta            ! air potential temperature at z=`pz`     [K]
      REAL(wp), INTENT(in) :: pqa              ! air specific humidity at z=`pz`    [kg/kg]
      REAL(wp)             :: T_from_z_P0_Theta_q_sclr ! air pot. temperature at z=`pz` [K]
      !!
      REAL(wp) :: zPz, zTa
      INTEGER  :: it
      !!-------------------------------------------------------------------------------
      zTa = pThta - rgamma_dry*pz   ! 1st guess of abs. temp. based on pot. temp
      !!
      DO it = 1, 4
         zPz = Pz_from_P0_tz_qz_sclr( pz, pslp, zTa, pqa ) ! pressure at z=`pz`
         zTa = abs_temp_sclr( pThta, zPz,  pPref=pslp )    ! update of absolute temperature
      END DO
      T_from_z_P0_Theta_q_sclr = zTa
      !!
   END FUNCTION T_from_z_P0_Theta_q_sclr

   FUNCTION T_from_z_P0_Theta_q_vctr( pz, pslp, pThta, pqa )
      REAL(wp),                     INTENT(in) :: pz
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pslp
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pThta
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pqa
      REAL(wp), DIMENSION(jpi,jpj)             :: T_from_z_P0_Theta_q_vctr
      INTEGER :: ji, jj
      DO jj = 1, jpj
         DO ji = 1, jpi
            T_from_z_P0_Theta_q_vctr(ji,jj) = T_from_z_P0_Theta_q_sclr( pz, pslp(ji,jj), pThta(ji,jj), pqa(ji,jj) )
         END DO
      END DO
   END FUNCTION T_from_z_P0_Theta_q_vctr
   !===============================================================================================





   !   !===============================================================================================
   !   FUNCTION Pz_from_P0_sclr( pz, pslp, pqa, pThta, pTa, l_ice )
   !      !!-------------------------------------------------------------------------------
   !      !!                           ***  FUNCTION Pz_from_P0  ***
   !      !!
   !      !! ** Purpose : compute air pressure at height `pz` m above sea level,
   !      !!              based barometric equation, from potential or absolute air temperature
   !      !!              and specific humidity at height `pz` and pressure at sea level
   !      !!
   !      !! ** Author: G. Samson,  Feb  2021
   !      !!            L. Brodeau, June 2021
   !      !!-------------------------------------------------------------------------------
   !      REAL(wp), INTENT(in )             :: pz               ! height above sea-level             [m]
   !      REAL(wp), INTENT(in )             :: pslp             ! pressure at sea-level (z=0)        [Pa]
   !      REAL(wp), INTENT(in )             :: pqa              ! air specific humidity at z=`pz`    [kg/kg]
   !      REAL(wp), INTENT(in )  , OPTIONAL :: pThta            ! air potential temperature          [K]
   !      REAL(wp), INTENT(inout), OPTIONAL :: pTa              ! air absolute temperature           [K]
   !      LOGICAL , INTENT(in)   , OPTIONAL :: l_ice            ! sea-ice presence
   !      REAL(wp)                          :: Pz_from_P0_sclr  ! pressure at `pz` m above sea level [Pa]
   !      !!
   !      REAL(wp)           :: ztpa, zta, zpa, zxm, zmask, zqsat, zf
   !      INTEGER            :: it
   !      LOGICAL            :: lice             ! sea-ice presence
   !      INTEGER, PARAMETER :: niter = 3    ! iteration indice and number
   !
   !      IF( PRESENT(pThta) ) THEN
   !         zmask = 1._wp
   !         ztpa = pThta
   !      ELSE
   !         zmask = 0._wp
   !         zta   = pTa
   !      ENDIF
   !
   !      lice = .FALSE.
   !      IF( PRESENT(l_ice) ) lice = l_ice
   !
   !      zpa = pslp              ! first guess of air pressure at zt   [Pa]
   !      DO it = 1, niter
   !         zta =            zmask  * ztpa*(zpa/Patm )**rpoiss_dry & ! if provided temp is potential => converstion to absolute temp !
   !            &  + (1._wp - zmask) * zta                            !  "           "      absolute  => keep unchanged
   !         !!
   !         zqsat = q_sat( zta, zpa, l_ice=lice )                               ! saturation specific humidity [kg/kg]
   !         zf    = pqa/zqsat
   !         zxm   = (1._wp - zf) * rmm_dryair + zf * rmm_water    ! moist air molar mass [kg/mol]
   !         zpa   = pslp * EXP( -grav * zxm * pz / ( R_gas * zta ) )
   !      END DO
   !
   !      Pz_from_P0_sclr = zpa
   !      IF(( PRESENT(pTa) ).AND.( PRESENT(pThta) )) pTa = zta   !! Absolute temperature is returned (it was not used as an input)
   !
   !   END FUNCTION Pz_from_P0_sclr
   !
   !   FUNCTION Pz_from_P0_vctr( pqa, pslp, pz, pThta, pTa, l_ice )
   !
   !      !!-------------------------------------------------------------------------------
   !      !!                           ***  FUNCTION Pz_from_P0  ***
   !      !!
   !      !! ** Purpose : compute air pressure using barometric equation
   !      !!              from either potential or absolute air temperature
   !      !! ** Author: G. Samson, Feb 2021
   !      !!-------------------------------------------------------------------------------
   !
   !      REAL(wp), DIMENSION(jpi,jpj)                          :: Pz_from_P0_vctr ! air pressure              [Pa]
   !      REAL(wp), DIMENSION(jpi,jpj), INTENT(in )             :: pqa            ! air specific humidity     [kg/kg]
   !      REAL(wp), DIMENSION(jpi,jpj), INTENT(in )             :: pslp           ! sea-level pressure        [Pa]
   !      REAL(wp),                     INTENT(in )             :: pz             ! height above surface      [m]
   !      REAL(wp), DIMENSION(jpi,jpj), INTENT(in )  , OPTIONAL :: pThta           ! air potential temperature [K]
   !      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout), OPTIONAL :: pTa            ! air absolute temperature  [K]
   !      LOGICAL                     , INTENT(in)   , OPTIONAL :: l_ice          ! sea-ice presence
   !      !!
   !      INTEGER                                               :: ji, jj         ! loop indices
   !      LOGICAL                                               :: lice           ! sea-ice presence
   !      lice = .FALSE.
   !      IF( PRESENT(l_ice) ) lice = l_ice
   !      IF( PRESENT(pThta) ) THEN
   !         DO jj = 1, jpj
   !            DO ji = 1, jpi
   !               Pz_from_P0_vctr(ji,jj) = Pz_from_P0_sclr( pqa(ji,jj), pslp(ji,jj), pz, pThta=pThta(ji,jj), pTa=pTa(ji,jj), l_ice=lice )
   !            END DO
   !         END DO
   !      ELSE
   !         DO jj = 1, jpj
   !            DO ji = 1, jpi
   !               Pz_from_P0_vctr(ji,jj) = Pz_from_P0_sclr( pqa(ji,jj), pslp(ji,jj), pz,                     pTa=pTa(ji,jj), l_ice=lice )
   !            END DO
   !         END DO
   !      ENDIF
   !   END FUNCTION Pz_from_P0_vctr
   !   !===============================================================================================




   !===============================================================================================
   FUNCTION rho_air_sclr( pTa, pqa, pslp )
      !!-------------------------------------------------------------------------------
      !!                           ***  FUNCTION rho_air_sclr  ***
      !!
      !! ** Purpose : compute density of (moist) air using the eq. of state of the atmosphere
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!-------------------------------------------------------------------------------
      REAL(wp), INTENT(in) :: pTa           ! absolute air temperature             [K]
      REAL(wp), INTENT(in) :: pqa            ! air specific humidity   [kg/kg]
      REAL(wp), INTENT(in) :: pslp           ! pressure in                [Pa]
      REAL(wp)             :: rho_air_sclr   ! density of moist air   [kg/m^3]
      !!-------------------------------------------------------------------------------
      rho_air_sclr = MAX( pslp / (R_dry*pTa * ( 1._wp + rctv0*pqa )) , 0.8_wp )
      !!
   END FUNCTION rho_air_sclr
   !!
   FUNCTION rho_air_vctr( pTa, pqa, pslp )
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pTa         ! absolute air temperature             [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pqa          ! air specific humidity   [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pslp          ! pressure in                [Pa]
      REAL(wp), DIMENSION(jpi,jpj)             :: rho_air_vctr ! density of moist air   [kg/m^3]
      !!-------------------------------------------------------------------------------
      rho_air_vctr = MAX( pslp / (R_dry*pTa * ( 1._wp + rctv0*pqa )) , 0.8_wp )
   END FUNCTION rho_air_vctr

   !===============================================================================================
   FUNCTION visc_air_sclr(pTa)
      !!----------------------------------------------------------------------------------
      !! Air kinetic viscosity (m^2/s) given from air temperature in Kelvin
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp)             :: visc_air_sclr   ! kinetic viscosity (m^2/s)
      REAL(wp), INTENT(in) :: pTa       ! absolute air temperature in [K]
      REAL(wp) ::   ztc, ztc2   ! local scalar
      !!----------------------------------------------------------------------------------
      ztc  = pTa - rt0   ! absolute air temp, in deg. C
      ztc2 = ztc*ztc
      visc_air_sclr = 1.326e-5*(1. + 6.542E-3*ztc + 8.301e-6*ztc2 - 4.84e-9*ztc2*ztc)
      !!
   END FUNCTION visc_air_sclr
   !!
   FUNCTION visc_air_vctr(pTa)
      REAL(wp), DIMENSION(jpi,jpj)             ::   visc_air_vctr   ! kinetic viscosity (m^2/s)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pTa       ! absolute air temperature in [K]
      INTEGER  ::   ji, jj      ! dummy loop indices
      DO jj = 1, jpj
         DO ji = 1, jpi
            visc_air_vctr(ji,jj) = visc_air_sclr( pTa(ji,jj) )
         END DO
      END DO
   END FUNCTION visc_air_vctr
   !===============================================================================================


   !===============================================================================================
   FUNCTION L_vap_sclr( psst )
      !!---------------------------------------------------------------------------------
      !!                           ***  FUNCTION L_vap_sclr  ***
      !!
      !! ** Purpose : Compute the latent heat of vaporization of water from temperature
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp)             :: L_vap_sclr  ! latent heat of vaporization   [J/kg]
      REAL(wp), INTENT(in) :: psst        ! water temperature               [K]
      !!----------------------------------------------------------------------------------
      L_vap_sclr = (  2.501_wp - 0.00237_wp * ( psst - rt0)  ) * 1.e6_wp
      !!
   END FUNCTION L_vap_sclr

   FUNCTION L_vap_vctr( psst )
      REAL(wp), DIMENSION(jpi,jpj)             :: L_vap_vctr  ! latent heat of vaporization [J/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: psst        ! water temperature             [K]
      L_vap_vctr = (  2.501_wp - 0.00237_wp * ( psst(:,:) - rt0)  ) * 1.e6_wp
   END FUNCTION L_vap_vctr
   !===============================================================================================


   !===============================================================================================
   FUNCTION cp_air_sclr( pqa )
      !!-------------------------------------------------------------------------------
      !!                           ***  FUNCTION cp_air_sclr  ***
      !!
      !! ** Purpose : Compute specific heat (Cp) of moist air
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!-------------------------------------------------------------------------------
      REAL(wp), INTENT(in) :: pqa           ! air specific humidity         [kg/kg]
      REAL(wp)             :: cp_air_sclr   ! specific heat of moist air   [J/K/kg]
      !!-------------------------------------------------------------------------------
      cp_air_sclr = rCp_dry + rCp_vap * pqa
      !!
   END FUNCTION cp_air_sclr
   !!
   FUNCTION cp_air_vctr( pqa )
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pqa      ! air specific humidity         [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj)             :: cp_air_vctr   ! specific heat of moist air   [J/K/kg]
      cp_air_vctr = rCp_dry + rCp_vap * pqa
   END FUNCTION cp_air_vctr
   !===============================================================================================


   !===============================================================================================
   FUNCTION gamma_moist_sclr( pTa, pqa )
      !!----------------------------------------------------------------------------------
      !! ** Purpose : Compute the moist adiabatic lapse-rate.
      !!     => http://glossary.ametsoc.org/wiki/Moist-adiabatic_lapse_rate
      !!     => http://www.geog.ucsb.edu/~joel/g266_s10/lecture_notes/chapt03/oh10_3_01/oh10_3_01.html
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp)             :: gamma_moist_sclr !                           [K/m]
      REAL(wp), INTENT(in) ::   pTa           ! absolute air temperature  [K] !LOLO: double check it's absolute !!!
      REAL(wp), INTENT(in) ::   pqa            ! specific humidity     [kg/kg]
      REAL(wp) :: zta, zqa, zwa, ziRT, zLvap        ! local scalars
      !!----------------------------------------------------------------------------------
      zta = MAX( pTa,  180._wp) ! prevents screw-up over masked regions where field == 0.
      zqa = MAX( pqa,  1.E-6_wp) !    "                   "                     "
      !!
      zwa = zqa / (1._wp - zqa)   ! w is mixing ratio w = q/(1-q) | q = w/(1+w)
      ziRT = 1._wp / (R_dry*zta)    ! 1/RT
      zLvap = L_vap_sclr( pTa )
      !!
      gamma_moist_sclr = grav * ( 1._wp + zLvap*zwa*ziRT ) / ( rCp_dry + zLvap*zLvap*zwa*reps0*ziRT/zta )
      !!
   END FUNCTION gamma_moist_sclr
   !!
   FUNCTION gamma_moist_vctr( pTa, pqa )
      REAL(wp), DIMENSION(jpi,jpj)             ::   gamma_moist_vctr
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pTa
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pqa
      INTEGER  :: ji, jj
      DO jj = 1, jpj
         DO ji = 1, jpi
            gamma_moist_vctr(ji,jj) = gamma_moist_sclr( pTa(ji,jj), pqa(ji,jj) )
         END DO
      END DO
   END FUNCTION gamma_moist_vctr
   !===============================================================================================


   !===============================================================================================
   FUNCTION One_on_L_sclr( pThta, pqa, pus, pts, pqs )
      !!------------------------------------------------------------------------
      !!
      !! Evaluates the 1./(Obukhov length) from air temperature,
      !! air specific humidity, and frictional scales u*, t* and q*
      !!
      !! Author: L. Brodeau, June 2019 / AeroBulk
      !!         (https://github.com/brodeau/aerobulk/)
      !!------------------------------------------------------------------------
      REAL(wp)             :: One_on_L_sclr !: 1./(Obukhov length) [m^-1]
      REAL(wp), INTENT(in) :: pThta     !: reference potential temperature of air [K]
      REAL(wp), INTENT(in) :: pqa      !: reference specific humidity of air   [kg/kg]
      REAL(wp), INTENT(in) :: pus      !: u*: friction velocity [m/s]
      REAL(wp), INTENT(in) :: pts, pqs !: \theta* and q* friction aka turb. scales for temp. and spec. hum.
      REAL(wp) ::     zqa          ! local scalar
      !!-------------------------------------------------------------------
      zqa = (1._wp + rctv0*pqa)
      !
      ! The main concern is to know whether, the vertical turbulent flux of virtual temperature, < u' theta_v' > is estimated with:
      !  a/  -u* [ theta* (1 + 0.61 q) + 0.61 theta q* ] => this is the one that seems correct! chose this one!
      !                      or
      !  b/  -u* [ theta*              + 0.61 theta q* ]
      !
      One_on_L_sclr = grav*vkarmn*( pts*zqa + rctv0*pThta*pqs ) / MAX( pus*pus*pThta*zqa , 1.E-9_wp )
      !
      One_on_L_sclr = SIGN( MIN(ABS(One_on_L_sclr),200._wp), One_on_L_sclr ) ! (prevent FPE from stupid values over masked regions...)
      !!
   END FUNCTION One_on_L_sclr
   !!
   FUNCTION One_on_L_vctr( pThta, pqa, pus, pts, pqs )
      REAL(wp), DIMENSION(jpi,jpj)             :: One_on_L_vctr !: 1./(Obukhov length) [m^-1]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pThta     !: reference potential temperature of air [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pqa      !: reference specific humidity of air   [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pus      !: u*: friction velocity [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pts, pqs !: \theta* and q* friction aka turb. scales for temp. and spec. hum.
      INTEGER  ::   ji, jj         ! dummy loop indices
      !!-------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
            One_on_L_vctr(ji,jj) = One_on_L_sclr( pThta(ji,jj), pqa(ji,jj), pus(ji,jj), pts(ji,jj), pqs(ji,jj) )
         END DO
      END DO
   END FUNCTION One_on_L_vctr
   !===============================================================================================

   !===============================================================================================
   FUNCTION Ri_bulk_sclr( pz, psst, pThta, pssq, pqa, pub,  pTa_layer, pqa_layer )
      !!----------------------------------------------------------------------------------
      !! Bulk Richardson number according to "wide-spread equation"...
      !!
      !!    Reminder: the Richardson number is the ratio "buoyancy" / "shear"
      !!
      !! ** Author: L. Brodeau, June 2019 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp)             :: Ri_bulk_sclr
      REAL(wp), INTENT(in) :: pz    ! height above the sea (aka "delta z")  [m]
      REAL(wp), INTENT(in) :: psst  ! SST                                   [K]
      REAL(wp), INTENT(in) :: pThta  ! pot. air temp. at height "pz"         [K]
      REAL(wp), INTENT(in) :: pssq  ! 0.98*q_sat(SST)                   [kg/kg]
      REAL(wp), INTENT(in) :: pqa   ! air spec. hum. at height "pz"     [kg/kg]
      REAL(wp), INTENT(in) :: pub   ! bulk wind speed                     [m/s]
      REAL(wp), INTENT(in), OPTIONAL :: pTa_layer ! when possible, a better guess of absolute temperature WITHIN the layer [K]
      REAL(wp), INTENT(in), OPTIONAL :: pqa_layer ! when possible, a better guess of specific humidity    WITHIN the layer [kg/kg]
      LOGICAL  :: l_ptqa_l_prvd = .FALSE.
      REAL(wp) :: zqa, zdthv, ztv, zsstv  ! local scalars
      !!-------------------------------------------------------------------
      IF( PRESENT(pTa_layer) .AND. PRESENT(pqa_layer) ) l_ptqa_l_prvd=.TRUE.
      !
      zsstv = virt_temp_sclr( psst, pssq )          ! virtual SST (absolute==potential because z=0!)
      !
      zdthv = virt_temp_sclr( pThta, pqa  ) - zsstv  ! air-sea delta of "virtual potential temperature"
      !
      !! ztv: estimate of the ABSOLUTE virtual temp. within the layer
      IF( l_ptqa_l_prvd ) THEN
         ztv = virt_temp_sclr( pTa_layer, pqa_layer )
      ELSE
         ztv = 0.5_wp*( zsstv + virt_temp_sclr( pThta-rgamma_dry*pz, pqa ) )
      END IF
      !
      Ri_bulk_sclr = grav*zdthv*pz / ( ztv*pub*pub )      ! the usual definition of Ri_bulk_sclr
      !!
   END FUNCTION Ri_bulk_sclr
   !!
   FUNCTION Ri_bulk_vctr( pz, psst, pThta, pssq, pqa, pub,  pTa_layer, pqa_layer )
      REAL(wp), DIMENSION(jpi,jpj)             :: Ri_bulk_vctr
      REAL(wp)                    , INTENT(in) :: pz
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: psst
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pThta
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pssq
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pqa
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pub
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in), OPTIONAL :: pTa_layer
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in), OPTIONAL :: pqa_layer
      LOGICAL  :: l_ptqa_l_prvd = .FALSE.
      INTEGER  ::   ji, jj
      IF( PRESENT(pTa_layer) .AND. PRESENT(pqa_layer) ) l_ptqa_l_prvd=.TRUE.
      DO jj = 1, jpj
         DO ji = 1, jpi
            IF( l_ptqa_l_prvd ) THEN
               Ri_bulk_vctr(ji,jj) = Ri_bulk_sclr( pz, psst(ji,jj), pThta(ji,jj), pssq(ji,jj), pqa(ji,jj), pub(ji,jj), &
                  &                                pTa_layer=pTa_layer(ji,jj ),  pqa_layer=pqa_layer(ji,jj ) )
            ELSE
               Ri_bulk_vctr(ji,jj) = Ri_bulk_sclr( pz, psst(ji,jj), pThta(ji,jj), pssq(ji,jj), pqa(ji,jj), pub(ji,jj) )
            END IF
         END DO
      END DO
   END FUNCTION Ri_bulk_vctr
   !===============================================================================================


   !===============================================================================================
   FUNCTION e_sat_sclr( pTa )
      !!----------------------------------------------------------------------------------
      !!                   ***  FUNCTION e_sat_sclr  ***
      !!                  < SCALAR argument version >
      !! ** Purpose : water vapor at saturation in [Pa]
      !!              Based on accurate estimate by Goff, 1957
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!
      !!    Note: what rt0 should be here, is 273.16 (triple point of water) and not 273.15 like here
      !!----------------------------------------------------------------------------------
      REAL(wp)             ::   e_sat_sclr   ! water vapor at saturation   [kg/kg]
      REAL(wp), INTENT(in) ::   pTa    ! absolute air temperature                  [K]
      REAL(wp) ::   zta, ztmp   ! local scalar
      !!----------------------------------------------------------------------------------
      zta = MAX( pTa , 180._wp )   ! absolute air temp., prevents fpe0 errors dute to unrealistically low values over masked regions...
      ztmp = rt0 / zta   !LOLO: rt0 or rtt0 ???? (273.15 vs 273.16 )
      !
      ! Vapour pressure at saturation [Pa] : WMO, (Goff, 1957)
      e_sat_sclr = 100.*( 10.**( 10.79574*(1. - ztmp) - 5.028*LOG10(zta/rt0)        &
         &    + 1.50475*10.**(-4)*(1. - 10.**(-8.2969*(zta/rt0 - 1.)) )  &
         &    + 0.42873*10.**(-3)*(10.**(4.76955*(1. - ztmp)) - 1.) + 0.78614) )
      !
   END FUNCTION e_sat_sclr
   !!
   FUNCTION e_sat_vctr(pTa)
      REAL(wp), DIMENSION(jpi,jpj)             :: e_sat_vctr !: vapour pressure at saturation  [Pa]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pTa    !: temperature [K]
      INTEGER  ::   ji, jj         ! dummy loop indices
      DO jj = 1, jpj
         DO ji = 1, jpi
            e_sat_vctr(ji,jj) = e_sat_sclr(pTa(ji,jj))
         END DO
      END DO
   END FUNCTION e_sat_vctr
   !===============================================================================================

   !===============================================================================================
   FUNCTION e_sat_ice_sclr(pTa)
      !!---------------------------------------------------------------------------------
      !! Same as "e_sat" but over ice rather than water!
      !!---------------------------------------------------------------------------------
      REAL(wp)             :: e_sat_ice_sclr !: vapour pressure at saturation in presence of ice [Pa]
      REAL(wp), INTENT(in) :: pTa
      !!
      REAL(wp) :: zta, zle, ztmp
      !!---------------------------------------------------------------------------------
      zta = MAX( pTa , 180._wp )   ! absolute air temp., prevents fpe0 errors dute to unrealistically low values over masked regions...
      ztmp = rtt0/zta
      !!
      zle  = rAg_i*(ztmp - 1._wp) + rBg_i*LOG10(ztmp) + rCg_i*(1._wp - zta/rtt0) + rDg_i
      !!
      e_sat_ice_sclr = 100._wp * 10._wp**zle
   END FUNCTION e_sat_ice_sclr
   !!
   FUNCTION e_sat_ice_vctr(pTa)
      !! Same as "e_sat" but over ice rather than water!
      REAL(wp), DIMENSION(jpi,jpj) :: e_sat_ice_vctr !: vapour pressure at saturation in presence of ice [Pa]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pTa
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
            e_sat_ice_vctr(ji,jj) = e_sat_ice_sclr( pTa(ji,jj) )
         END DO
      END DO
   END FUNCTION e_sat_ice_vctr
   !!
   FUNCTION de_sat_dt_ice_sclr(pTa)
      !!---------------------------------------------------------------------------------
      !! d [ e_sat_ice ] / dT   (derivative / temperature)
      !! Analytical exact formulation: double checked!!!
      !!  => DOUBLE-check possible / finite-difference version with "./bin/test_phymbl.x"
      !!---------------------------------------------------------------------------------
      REAL(wp)             :: de_sat_dt_ice_sclr !:  [Pa/K]
      REAL(wp), INTENT(in) :: pTa
      !!
      REAL(wp) :: zta, zde
      !!---------------------------------------------------------------------------------
      zta = MAX( pTa , 180._wp )   ! absolute air temp., prevents fpe0 errors dute to unrealistically low values over masked regions...
      !!
      zde = -(rAg_i*rtt0)/(zta*zta) - rBg_i/(zta*LOG(10._wp)) - rCg_i/rtt0
      !!
      de_sat_dt_ice_sclr = LOG(10._wp) * zde * e_sat_ice_sclr(zta)
   END FUNCTION de_sat_dt_ice_sclr
   !!
   FUNCTION de_sat_dt_ice_vctr(pTa)
      !! Same as "e_sat" but over ice rather than water!
      REAL(wp), DIMENSION(jpi,jpj) :: de_sat_dt_ice_vctr !: vapour pressure at saturation in presence of ice [Pa]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pTa
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
            de_sat_dt_ice_vctr(ji,jj) = de_sat_dt_ice_sclr( pTa(ji,jj) )
         END DO
      END DO
   END FUNCTION de_sat_dt_ice_vctr



   !===============================================================================================

   !===============================================================================================
   FUNCTION q_sat_sclr( pTa, pslp,  l_ice )
      !!---------------------------------------------------------------------------------
      !!                           ***  FUNCTION q_sat_sclr  ***
      !!
      !! ** Purpose : Conputes specific humidity of air at saturation
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp) :: q_sat_sclr
      REAL(wp), INTENT(in) :: pTa  !: absolute temperature of air [K]
      REAL(wp), INTENT(in) :: pslp  !: atmospheric pressure        [Pa]
      LOGICAL,  INTENT(in), OPTIONAL :: l_ice  !: we are above ice
      REAL(wp) :: ze_s
      LOGICAL  :: lice
      !!----------------------------------------------------------------------------------
      lice = .FALSE.
      IF( PRESENT(l_ice) ) lice = l_ice
      IF( lice ) THEN
         ze_s = e_sat_ice( pTa )
      ELSE
         ze_s = e_sat( pTa ) ! Vapour pressure at saturation (Goff) :
      END IF
      q_sat_sclr = reps0*ze_s/(pslp - (1._wp - reps0)*ze_s)
   END FUNCTION q_sat_sclr
   !!
   FUNCTION q_sat_vctr( pTa, pslp,  l_ice )
      REAL(wp), DIMENSION(jpi,jpj) :: q_sat_vctr
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pTa  !: absolute temperature of air [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pslp  !: atmospheric pressure        [Pa]
      LOGICAL,  INTENT(in), OPTIONAL :: l_ice  !: we are above ice
      LOGICAL  :: lice
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------------------
      lice = .FALSE.
      IF( PRESENT(l_ice) ) lice = l_ice
      DO jj = 1, jpj
         DO ji = 1, jpi
            q_sat_vctr(ji,jj) = q_sat_sclr( pTa(ji,jj) , pslp(ji,jj), l_ice=lice )
         END DO
      END DO
   END FUNCTION q_sat_vctr
   !===============================================================================================


   !===============================================================================================
   FUNCTION dq_sat_dt_ice_sclr( pTa, pslp )
      !!---------------------------------------------------------------------------------
      !!     ***  FUNCTION dq_sat_dt_ice_sclr  ***
      !!    => d [ q_sat_ice(T) ] / dT
      !! Analytical exact formulation: double checked!!!
      !!  => DOUBLE-check possible / finite-difference version with "./bin/test_phymbl.x"
      !!----------------------------------------------------------------------------------
      REAL(wp) :: dq_sat_dt_ice_sclr
      REAL(wp), INTENT(in) :: pTa  !: absolute temperature of air [K]
      REAL(wp), INTENT(in) :: pslp  !: atmospheric pressure        [Pa]
      REAL(wp) :: ze_s, zde_s_dt, ztmp
      !!----------------------------------------------------------------------------------
      ze_s     =  e_sat_ice_sclr( pTa ) ! Vapour pressure at saturation  in presence of ice (Goff)
      zde_s_dt = de_sat_dt_ice(   pTa )
      !
      ztmp = (reps0 - 1._wp)*ze_s + pslp
      !
      dq_sat_dt_ice_sclr = reps0*pslp*zde_s_dt / ( ztmp*ztmp )
      !
   END FUNCTION dq_sat_dt_ice_sclr
   !!
   FUNCTION dq_sat_dt_ice_vctr( pTa, pslp )
      REAL(wp), DIMENSION(jpi,jpj) :: dq_sat_dt_ice_vctr
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pTa  !: absolute temperature of air [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pslp  !: atmospheric pressure        [Pa]
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
            dq_sat_dt_ice_vctr(ji,jj) = dq_sat_dt_ice_sclr( pTa(ji,jj) , pslp(ji,jj) )
         END DO
      END DO
   END FUNCTION dq_sat_dt_ice_vctr
   !===============================================================================================


   !===============================================================================================
   FUNCTION q_air_rh(prha, pTa, pslp)
      !!----------------------------------------------------------------------------------
      !! Specific humidity of air out of Relative Humidity
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: q_air_rh
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: prha        !: relative humidity       [%]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pTa         !: absolute air temperature [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pslp        !: atmospheric pressure    [Pa]
      !
      INTEGER  ::   ji, jj      ! dummy loop indices
      REAL(wp) ::   ze      ! local scalar
      !!----------------------------------------------------------------------------------
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            ze = 0.01_wp*prha(ji,jj)*e_sat_sclr(pTa(ji,jj))
            q_air_rh(ji,jj) = ze*reps0/MAX( pslp(ji,jj) - (1. - reps0)*ze , 1._wp )
         END DO
      END DO
      !
   END FUNCTION q_air_rh




   FUNCTION q_air_dp(da, slp)
      !!
      !! Air specific humidity from dew point temperature
      !!
      REAL(wp), DIMENSION(jpi,jpj) :: q_air_dp  !: kg/kg
      !!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: &
         &     da,     &    !: dew-point temperature   [K]
         &     slp         !: atmospheric pressure    [Pa]
      !!
      q_air_dp = MAX( e_sat(da) , 0._wp ) ! q_air_dp is e_sat !!!
      q_air_dp = q_air_dp*reps0/MAX( slp - (1. - reps0)*q_air_dp, 1._wp )  ! MAX() are solely here to prevent NaN
      !!                                                                   ! over masked regions with silly values
   END FUNCTION q_air_dp







   FUNCTION rho_air_adv(pTa, pqa, pslp)
      !!
      !! Advanced version, using TRUE virtual temperature
      !!
      REAL(wp), DIMENSION(jpi,jpj) :: rho_air_adv      !: density of air [kg/m^3]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pTa  !: absolute air temperature in [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pqa  !: air spec. hum. [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pslp !: atm. pressure in       [Pa]
      !!
      ! " pTa/(1. - e_air(pqa, pslp)/pslp*(1. - reps0)) " is virtual absolute temperature !
      rho_air_adv = pslp / ( R_dry * pTa/(1. - e_air(pqa, pslp)/pslp*(1. - reps0)) )
      !!
   END FUNCTION rho_air_adv








   FUNCTION q_sat_crude(temp, zrho)

      REAL(wp), DIMENSION(jpi,jpj)             :: q_sat_crude
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: &
         &       temp,   &  !: sea surface temperature  [K]
         &       zrho       !: air density         [kg/m^3]

      q_sat_crude = 640380._wp/zrho * exp(-5107.4_wp/temp)

   END FUNCTION q_sat_crude







   FUNCTION dry_static_energy( pz, pTa, pqa )
      !!----------------------------------------------------------------------------------
      !! Dry static energy "s" (Eq. 3.5 IFS doc)
      !!
      !! ** Author: L. Brodeau, June 2019 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             ::   dry_static_energy
      REAL(wp)                    , INTENT(in) ::   pz    ! height above the sea         [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pTa   ! absolute air temp. at pz m   [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pqa   ! air spec. hum. at pz m   [kg/kg]
      !!----------------------------------------------------------------------------------
      dry_static_energy = grav*pz + cp_air(pqa)*pTa
   END FUNCTION dry_static_energy


   !===============================================================================================
   SUBROUTINE UPDATE_QNSOL_TAU_SCLR( pzu, pts, pqs, pThta, pqa, pust, ptst, pqst, pwnd, pUb, pslp, prlw, &
      &                              pQns, pTau,    Qlat )
      !!----------------------------------------------------------------------------------
      !! Purpose: returns the non-solar heat flux to the ocean aka "Qlat + Qsen + Qlw"
      !!          and the module of the wind stress => pTau = Tau
      !! ** Author: L. Brodeau, Sept. 2019 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(in)  :: pzu  ! height above the sea-level where all this takes place (normally 10m)
      REAL(wp), INTENT(in)  :: pts  ! water temperature at the air-sea interface [K]
      REAL(wp), INTENT(in)  :: pqs  ! satur. spec. hum. at T=pts   [kg/kg]
      REAL(wp), INTENT(in)  :: pThta ! potential air temperature at z=pzu [K]
      REAL(wp), INTENT(in)  :: pqa  ! specific humidity at z=pzu [kg/kg]
      REAL(wp), INTENT(in)  :: pust ! u*
      REAL(wp), INTENT(in)  :: ptst ! t*
      REAL(wp), INTENT(in)  :: pqst ! q*
      REAL(wp), INTENT(in)  :: pwnd ! wind speed module at z=pzu [m/s]
      REAL(wp), INTENT(in)  :: pUb  ! bulk wind speed at z=pzu (inc. pot. effect of gustiness etc) [m/s]
      REAL(wp), INTENT(in)  :: pslp ! sea-level atmospheric pressure [Pa]
      REAL(wp), INTENT(in)  :: prlw ! downwelling longwave radiative flux [W/m^2]
      !
      REAL(wp), INTENT(out) :: pQns ! non-solar heat flux to the ocean aka "Qlat + Qsen + Qlw" [W/m^2]]
      REAL(wp), INTENT(out) :: pTau ! module of the wind stress [N/m^2]
      !
      REAL(wp), OPTIONAL, INTENT(out) :: Qlat
      !
      REAL(wp) :: zdt, zdq, zCd, zCh, zCe, zz0, zQlat, zQsen, zQlw
      !!----------------------------------------------------------------------------------
      zdt = pThta - pts ;  zdt = SIGN( MAX(ABS(zdt),1.E-6_wp), zdt )
      zdq = pqa - pqs ;  zdq = SIGN( MAX(ABS(zdq),1.E-9_wp), zdq )
      zz0 = pust/pUb
      zCd = zz0*zz0
      zCh = zz0*ptst/zdt
      zCe = zz0*pqst/zdq

      CALL BULK_FORMULA_SCLR( pzu, pts, pqs, pThta, pqa, zCd, zCh, zCe, &
         &                    pwnd, pUb, pslp,                         &
         &                    pTau, zQsen, zQlat )

      zQlw = qlw_net_sclr( prlw, pts ) ! Net longwave flux

      pQns = zQlat + zQsen + zQlw

      IF( PRESENT(Qlat) ) Qlat = zQlat

   END SUBROUTINE UPDATE_QNSOL_TAU_SCLR
   !!
   SUBROUTINE UPDATE_QNSOL_TAU_VCTR( pzu, pts, pqs, pThta, pqa, pust, ptst, pqst, pwnd, pUb, pslp, prlw, &
      &                              pQns, pTau,    Qlat)
      !!----------------------------------------------------------------------------------
      !! Purpose: returns the non-solar heat flux to the ocean aka "Qlat + Qsen + Qlw"
      !!          and the module of the wind stress => pTau = Tau
      !! ** Author: L. Brodeau, Sept. 2019 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp),                     INTENT(in)  :: pzu  ! height above the sea-level where all this takes place (normally 10m)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pts  ! water temperature at the air-sea interface [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pqs  ! satur. spec. hum. at T=pts   [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pThta  ! potential air temperature at z=pzu [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pqa  ! specific humidity at z=pzu [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pust ! u*
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: ptst ! t*
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pqst ! q*
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pwnd ! wind speed module at z=pzu [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pUb  ! bulk wind speed at z=pzu (inc. pot. effect of gustiness etc) [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pslp ! sea-level atmospheric pressure [Pa]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: prlw ! downwelling longwave radiative flux [W/m^2]
      !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pQns ! non-solar heat flux to the ocean aka "Qlat + Qsen + Qlw" [W/m^2]]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pTau ! module of the wind stress [N/m^2]
      !
      REAL(wp), DIMENSION(jpi,jpj), OPTIONAL, INTENT(out) :: Qlat
      !
      REAL(wp) :: zQlat
      INTEGER  ::   ji, jj     ! dummy loop indices
      LOGICAL  :: lrQlat=.false.
      !!----------------------------------------------------------------------------------
      lrQlat = PRESENT(Qlat)
      DO jj = 1, jpj
         DO ji = 1, jpi
            CALL UPDATE_QNSOL_TAU_SCLR( pzu, pts(ji,jj),  pqs(ji,jj),  pThta(ji,jj), pqa(ji,jj),  &
               &                            pust(ji,jj), ptst(ji,jj), pqst(ji,jj), pwnd(ji,jj), &
               &                             pUb(ji,jj),  pslp(ji,jj), prlw(ji,jj),              &
               &                            pQns(ji,jj), pTau(ji,jj),              Qlat=zQlat   )
            IF( lrQlat ) Qlat(ji,jj) = zQlat
         END DO
      END DO
   END SUBROUTINE UPDATE_QNSOL_TAU_VCTR
   !===============================================================================================


   !===============================================================================================
   SUBROUTINE BULK_FORMULA_SCLR( pzu, pts, pqs, pThta, pqa, &
      &                          pCd, pCh, pCe,            &
      &                          pwnd, pUb, pslp,          &
      &                          pTau, pQsen, pQlat,       &
      &                          pEvap, prhoa, l_ice      )
      !!----------------------------------------------------------------------------------
      REAL(wp),                     INTENT(in)  :: pzu  ! height above the sea-level where all this takes place (normally 10m)
      REAL(wp), INTENT(in)  :: pts  ! water temperature at the air-sea interface [K]
      REAL(wp), INTENT(in)  :: pqs  ! satur. spec. hum. at T=pts   [kg/kg]
      REAL(wp), INTENT(in)  :: pThta  ! potential air temperature at z=pzu [K]
      REAL(wp), INTENT(in)  :: pqa  ! specific humidity at z=pzu [kg/kg]
      REAL(wp), INTENT(in)  :: pCd
      REAL(wp), INTENT(in)  :: pCh
      REAL(wp), INTENT(in)  :: pCe
      REAL(wp), INTENT(in)  :: pwnd ! wind speed module at z=pzu [m/s]
      REAL(wp), INTENT(in)  :: pUb  ! bulk wind speed at z=pzu (inc. pot. effect of gustiness etc) [m/s]
      REAL(wp), INTENT(in)  :: pslp ! sea-level atmospheric pressure [Pa]
      !!
      REAL(wp), INTENT(out) :: pTau  ! module of the wind stress [N/m^2]
      REAL(wp), INTENT(out) :: pQsen !  [W/m^2]
      REAL(wp), INTENT(out) :: pQlat !  [W/m^2]
      !!
      REAL(wp), INTENT(out), OPTIONAL :: pEvap ! Evaporation [kg/m^2/s] (usually <0, as ocean loses water!)
      REAL(wp), INTENT(out), OPTIONAL :: prhoa ! Air density at z=pzu [kg/m^3]
      LOGICAL,  INTENT(in),  OPTIONAL :: l_ice  !: we are above ice
      !!
      REAL(wp) :: zta, zrho, zUrho, zevap
      LOGICAL  :: lice
      !!----------------------------------------------------------------------------------
      lice = .FALSE.
      IF( PRESENT(l_ice) ) lice = l_ice

      !! Need zta, absolute temperature at pzu (formula to estimate rho_air needs absolute temperature, not the potential temperature "pThta")
      zta  = pThta - rgamma_dry*pzu   ! Absolute temp. is slightly colder...
      zrho = rho_air(zta, pqa, pslp)
      zrho = rho_air(zta, pqa, pslp-zrho*grav*pzu) ! taking into account that we are pzu m above the sea level where SLP is given!

      zUrho = pUb*MAX(zrho, 1._wp)     ! rho*U10

      pTau = zUrho * pCd * pwnd ! Wind stress module

      zevap = zUrho * pCe * (pqa - pqs)
      pQsen = zUrho * pCh * (pThta - pts) * cp_air(pqa)

      IF( lice) THEN
         pQlat =      rLsub * zevap
         IF( PRESENT(pEvap) ) pEvap = MIN( zevap , 0._wp )
      ELSE
         pQlat = L_vap(pts) * zevap
         IF( PRESENT(pEvap) ) pEvap = zevap
      END IF

      IF( PRESENT(prhoa) ) prhoa = zrho

   END SUBROUTINE BULK_FORMULA_SCLR
   !!
   SUBROUTINE BULK_FORMULA_VCTR( pzu, pts, pqs, pThta, pqa, &
      &                          pCd, pCh, pCe,           &
      &                          pwnd, pUb, pslp,         &
      &                          pTau, pQsen, pQlat,      &
      &                          pEvap, prhoa, l_ice )
      !!----------------------------------------------------------------------------------
      REAL(wp),                     INTENT(in)  :: pzu  ! height above the sea-level where all this takes place (normally 10m)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pts  ! water temperature at the air-sea interface [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pqs  ! satur. spec. hum. at T=pts   [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pThta  ! potential air temperature at z=pzu [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pqa  ! specific humidity at z=pzu [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pCd
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pCh
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pCe
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pwnd ! wind speed module at z=pzu [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pUb  ! bulk wind speed at z=pzu (inc. pot. effect of gustiness etc) [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pslp ! sea-level atmospheric pressure [Pa]
      !!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pTau  ! module of the wind stress [N/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pQsen !  [W/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pQlat !  [W/m^2]
      !!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out), OPTIONAL :: pEvap ! Evaporation [kg/m^2/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out), OPTIONAL :: prhoa ! Air density at z=pzu [kg/m^3]
      LOGICAL,  INTENT(in),  OPTIONAL :: l_ice  !: we are above ice
      !!
      REAL(wp) :: zevap, zrho
      LOGICAL  :: lice=.FALSE., lrE=.FALSE., lrR=.FALSE.
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------------------
      lice = PRESENT(l_ice)
      lrE  = PRESENT(pEvap)
      lrR  = PRESENT(prhoa)
      DO jj = 1, jpj
         DO ji = 1, jpi
            CALL BULK_FORMULA_SCLR( pzu, pts(ji,jj), pqs(ji,jj), pThta(ji,jj), pqa(ji,jj), &
               &                    pCd(ji,jj), pCh(ji,jj), pCe(ji,jj),                  &
               &                    pwnd(ji,jj), pUb(ji,jj), pslp(ji,jj),                &
               &                    pTau(ji,jj), pQsen(ji,jj), pQlat(ji,jj),             &
               &                    pEvap=zevap, prhoa=zrho, l_ice=lice )
            IF( lrE ) pEvap(ji,jj) = zevap
            IF( lrR ) prhoa(ji,jj) = zrho
         END DO
      END DO
   END SUBROUTINE BULK_FORMULA_VCTR
   !===============================================================================================



   !===============================================================================================
   FUNCTION alpha_sw_sclr( psst )
      !!---------------------------------------------------------------------------------
      !!                           ***  FUNCTION alpha_sw_sclr  ***
      !!
      !! ** Purpose : ROUGH estimate of the thermal expansion coefficient of sea-water at the surface (P =~ 1010 hpa)
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp)             ::   alpha_sw_sclr   ! thermal expansion coefficient of sea-water [1/K]
      REAL(wp), INTENT(in) ::   psst   ! sea-water temperature                   [K]
      !!----------------------------------------------------------------------------------
      alpha_sw_sclr = 2.1e-5_wp * MAX( psst - rt0 + 3.2_wp , 0._wp )**0.79_wp
      !!
   END FUNCTION alpha_sw_sclr
   !!
   FUNCTION alpha_sw_vctr( psst )
      REAL(wp), DIMENSION(jpi,jpj)             ::   alpha_sw_vctr   ! thermal expansion coefficient of sea-water [1/K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   psst   ! water temperature                [K]
      alpha_sw_vctr = 2.1e-5_wp * MAX( psst(:,:) - rt0 + 3.2_wp , 0._wp )**0.79_wp
   END FUNCTION alpha_sw_vctr
   !===============================================================================================


   !===============================================================================================
   FUNCTION qlw_net_sclr( pdwlw, pts,  l_ice )
      !!---------------------------------------------------------------------------------
      !!                           ***  FUNCTION qlw_net_sclr  ***
      !!
      !! ** Purpose : Estimate of the net longwave flux at the surface
      !!----------------------------------------------------------------------------------
      REAL(wp) :: qlw_net_sclr
      REAL(wp), INTENT(in) :: pdwlw !: downwelling longwave (aka infrared, aka thermal) radiation [W/m^2]
      REAL(wp), INTENT(in) :: pts   !: surface temperature [K]
      LOGICAL,  INTENT(in), OPTIONAL :: l_ice  !: we are above ice
      REAL(wp) :: zemiss, zt2
      LOGICAL  :: lice
      !!----------------------------------------------------------------------------------
      lice = .FALSE.
      IF( PRESENT(l_ice) ) lice = l_ice
      IF( lice ) THEN
         zemiss = emiss_i
      ELSE
         zemiss = emiss_w
      END IF
      zt2 = pts*pts
      qlw_net_sclr = zemiss*( pdwlw - stefan*zt2*zt2)  ! zemiss used both as the IR albedo and IR emissivity...
      !!
   END FUNCTION qlw_net_sclr
   !!
   FUNCTION qlw_net_vctr( pdwlw, pts,  l_ice )
      REAL(wp), DIMENSION(jpi,jpj) :: qlw_net_vctr
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pdwlw !: downwelling longwave (aka infrared, aka thermal) radiation [W/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pts   !: surface temperature [K]
      LOGICAL,  INTENT(in), OPTIONAL :: l_ice  !: we are above ice
      LOGICAL  :: lice
      INTEGER  :: ji, jj
      lice = .FALSE.
      IF( PRESENT(l_ice) ) lice = l_ice
      DO jj = 1, jpj
         DO ji = 1, jpi
            qlw_net_vctr(ji,jj) = qlw_net_sclr( pdwlw(ji,jj) , pts(ji,jj), l_ice=lice )
         END DO
      END DO
   END FUNCTION qlw_net_vctr
   !===============================================================================================


   !===============================================================================================
   FUNCTION z0_from_Cd_sclr( pzu, pCd,  ppsi )
      REAL(wp)                       :: z0_from_Cd_sclr        !: roughness length [m]
      REAL(wp), INTENT(in)           :: pzu   !: reference height zu [m]
      REAL(wp), INTENT(in)           :: pCd   !: (neutral or non-neutral) drag coefficient []
      REAL(wp), INTENT(in), OPTIONAL :: ppsi !: "Psi_m(pzu/L)" stability correction profile for momentum []
      !!
      !! If pCd is the NEUTRAL-STABILITY drag coefficient then ppsi must be 0 or not given
      !! If pCd is the drag coefficient (in stable or unstable conditions) then pssi must be provided
      !!----------------------------------------------------------------------------------
      IF( PRESENT(ppsi) ) THEN
         !! Cd provided is the actual Cd (not the neutral-stability CdN) :
         z0_from_Cd_sclr = pzu * EXP( - ( vkarmn/SQRT(pCd) + ppsi ) ) !LB: ok, double-checked!
      ELSE
         !! Cd provided is the neutral-stability Cd, aka CdN :
         z0_from_Cd_sclr = pzu * EXP( - vkarmn/SQRT(pCd) )            !LB: ok, double-checked!
      END IF
      !!
   END FUNCTION z0_from_Cd_sclr

   FUNCTION z0_from_Cd_vctr( pzu, pCd,  ppsi )
      REAL(wp), DIMENSION(jpi,jpj) :: z0_from_Cd_vctr        !: roughness length [m]
      REAL(wp)                    , INTENT(in) :: pzu   !: reference height zu [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pCd   !: (neutral or non-neutral) drag coefficient []
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in), OPTIONAL :: ppsi !: "Psi_m(pzu/L)" stability correction profile for momentum []
      !!----------------------------------------------------------------------------------
      IF( PRESENT(ppsi) ) THEN
         z0_from_Cd_vctr = pzu * EXP( - ( vkarmn/SQRT(pCd(:,:)) + ppsi(:,:) ) ) !LB: ok, double-checked!
      ELSE
         z0_from_Cd_vctr = pzu * EXP( - vkarmn/SQRT(pCd(:,:)) )            !LB: ok, double-checked!
      END IF
      !!
   END FUNCTION z0_from_Cd_vctr
   !===============================================================================================


   !===============================================================================================
   FUNCTION z0_from_ustar_sclr( pzu, pus, puzu )
      !!----------------------------------------------------------------------------------
      REAL(wp)             :: z0_from_ustar_sclr    !: roughness length    [m]
      REAL(wp), INTENT(in) :: pzu  !: reference height zu [m]
      REAL(wp), INTENT(in) :: pus  !: friction velocity   [m/s]
      REAL(wp), INTENT(in) :: puzu !: wind speed at z=pzu [m/s]
      !!----------------------------------------------------------------------------------
      z0_from_ustar_sclr = pzu * EXP( -  vkarmn*puzu/pus )
      !!
   END FUNCTION z0_from_ustar_sclr
   !!
   FUNCTION z0_from_ustar_vctr( pzu, pus, puzu )
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: z0_from_ustar_vctr    !: roughness length    [m]
      REAL(wp)                    , INTENT(in) :: pzu  !: reference height zu [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pus  !: friction velocity   [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: puzu !: wind speed at z=pzu [m/s]
      !!----------------------------------------------------------------------------------
      z0_from_ustar_vctr = pzu * EXP( -  vkarmn*puzu(:,:)/pus(:,:) )
      !!
   END FUNCTION z0_from_ustar_vctr
   !===============================================================================================


   !===============================================================================================
   FUNCTION Cd_from_z0( pzu, pz0,  ppsi )
      REAL(wp), DIMENSION(jpi,jpj) :: Cd_from_z0        !: (neutral or non-neutral) drag coefficient []
      REAL(wp)                    , INTENT(in) :: pzu   !: reference height zu [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pz0   !: roughness length [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in), OPTIONAL :: ppsi !: "Psi_m(pzu/L)" stability correction profile for momentum []
      !!
      !! If we want to return the NEUTRAL-STABILITY drag coefficient then ppsi must be 0 or not given
      !! If we want to return the stability-corrected Cd (i.e. in stable or unstable conditions) then pssi must be provided
      !!----------------------------------------------------------------------------------
      IF( PRESENT(ppsi) ) THEN
         !! The Cd we return is the actual Cd (not the neutral-stability CdN) :
         Cd_from_z0 = 1._wp / ( LOG( pzu / pz0(:,:) ) - ppsi(:,:) )
      ELSE
         !! The Cd we return is the neutral-stability Cd, aka CdN :
         Cd_from_z0 = 1._wp /   LOG( pzu / pz0(:,:) )
      END IF
      Cd_from_z0 = vkarmn2 * Cd_from_z0 * Cd_from_z0
      !!
   END FUNCTION Cd_from_z0
   !===============================================================================================


   !===============================================================================================
   FUNCTION f_m_louis_sclr( pzu, pRib, pCdn, pz0 )
      !!----------------------------------------------------------------------------------
      !!  Stability correction function for MOMENTUM
      !!                 Louis (1979)
      !!----------------------------------------------------------------------------------
      REAL(wp)             :: f_m_louis_sclr ! term "f_m" in Eq.(6) when option "Louis" rather than "Psi(zeta) is chosen, Lupkes & Gryanik (2015),
      REAL(wp), INTENT(in) :: pzu     ! reference height (height for pwnd)  [m]
      REAL(wp), INTENT(in) :: pRib    ! Bulk Richardson number
      REAL(wp), INTENT(in) :: pCdn    ! neutral drag coefficient
      REAL(wp), INTENT(in) :: pz0     ! roughness length                      [m]
      !!----------------------------------------------------------------------------------
      REAL(wp) :: ztu, zts, zstab
      !!----------------------------------------------------------------------------------
      zstab = 0.5 + SIGN(0.5_wp, pRib) ; ! Unstable (Ri<0) => zstab = 0 | Stable (Ri>0) => zstab = 1
      !
      ztu = pRib / ( 1._wp + 3._wp * rc2_louis * pCdn * SQRT( ABS( -pRib * ( pzu / pz0 + 1._wp) ) ) ) ! ABS is just here for when it's stable conditions and ztu is not used anyways
      zts = pRib / SQRT( ABS( 1._wp + pRib ) ) ! ABS is just here for when it's UNstable conditions and zts is not used anyways
      !
      f_m_louis_sclr = (1._wp - zstab) *         ( 1._wp - ram_louis * ztu )  &  ! Unstable Eq.(A6)
         &               +      zstab  * 1._wp / ( 1._wp + ram_louis * zts )     ! Stable   Eq.(A7)
      !!
   END FUNCTION f_m_louis_sclr
   !!
   FUNCTION f_m_louis_vctr( pzu, pRib, pCdn, pz0 )
      REAL(wp), DIMENSION(jpi,jpj)             :: f_m_louis_vctr
      REAL(wp),                     INTENT(in) :: pzu
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pRib
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pCdn
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pz0
      INTEGER  :: ji, jj
      DO jj = 1, jpj
         DO ji = 1, jpi
            f_m_louis_vctr(ji,jj) = f_m_louis_sclr( pzu, pRib(ji,jj), pCdn(ji,jj), pz0(ji,jj) )
         END DO
      END DO
   END FUNCTION f_m_louis_vctr
   !===============================================================================================

   !===============================================================================================
   FUNCTION f_h_louis_sclr( pzu, pRib, pChn, pz0 )
      !!----------------------------------------------------------------------------------
      !!  Stability correction function for HEAT
      !!                 Louis (1979)
      !!----------------------------------------------------------------------------------
      REAL(wp)             :: f_h_louis_sclr ! term "f_h" in Eq.(6) when option "Louis" rather than "Psi(zeta) is chosen, Lupkes & Gryanik (2015),
      REAL(wp), INTENT(in) :: pzu     ! reference height (height for pwnd)  [m]
      REAL(wp), INTENT(in) :: pRib    ! Bulk Richardson number
      REAL(wp), INTENT(in) :: pChn    ! neutral heat transfer coefficient
      REAL(wp), INTENT(in) :: pz0     ! roughness length                      [m]
      !!----------------------------------------------------------------------------------
      REAL(wp) :: ztu, zts, zstab
      !!----------------------------------------------------------------------------------
      zstab = 0.5 + SIGN(0.5_wp, pRib) ; ! Unstable (Ri<0) => zstab = 0 | Stable (Ri>0) => zstab = 1
      !
      ztu = pRib / ( 1._wp + 3._wp * rc2_louis * pChn * SQRT( ABS(-pRib * ( pzu / pz0 + 1._wp) ) ) )
      zts = pRib / SQRT( ABS( 1._wp + pRib ) )
      !
      f_h_louis_sclr = (1._wp - zstab) *         ( 1._wp - rah_louis * ztu )  &  ! Unstable Eq.(A6)
         &              +       zstab  * 1._wp / ( 1._wp + rah_louis * zts )     ! Stable   Eq.(A7)  !LOLO: in paper it's "ram_louis" and not "rah_louis" typo or what????
      !!
   END FUNCTION f_h_louis_sclr
   !!
   FUNCTION f_h_louis_vctr( pzu, pRib, pChn, pz0 )
      REAL(wp), DIMENSION(jpi,jpj)             :: f_h_louis_vctr
      REAL(wp),                     INTENT(in) :: pzu
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pRib
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pChn
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pz0
      INTEGER  :: ji, jj
      DO jj = 1, jpj
         DO ji = 1, jpi
            f_h_louis_vctr(ji,jj) = f_h_louis_sclr( pzu, pRib(ji,jj), pChn(ji,jj), pz0(ji,jj) )
         END DO
      END DO
   END FUNCTION f_h_louis_vctr
   !===============================================================================================


   !===============================================================================================
   FUNCTION UN10_from_ustar( pzu, pUzu, pus, ppsi )
      !!----------------------------------------------------------------------------------
      !!  Provides the neutral-stability wind speed at 10 m
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: UN10_from_ustar  !: neutral stability wind speed at 10m [m/s]
      REAL(wp),                     INTENT(in) :: pzu   !: measurement heigh of wind speed   [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pUzu  !: bulk wind speed at height pzu m   [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pus   !: friction velocity                 [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: ppsi !: "Psi_m(pzu/L)" stability correction profile for momentum []
      !!----------------------------------------------------------------------------------
      UN10_from_ustar(:,:) = pUzu(:,:) - pus(:,:)/vkarmn * ( LOG(pzu/10._wp) - ppsi(:,:) )
      !!
   END FUNCTION UN10_from_ustar
   !===============================================================================================


   !===============================================================================================
   FUNCTION UN10_from_CDN( pzu, pUb, pCdn, ppsi )
      !!----------------------------------------------------------------------------------
      !!  Provides the neutral-stability wind speed at 10 m
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: UN10_from_CDN  !: [m/s]
      REAL(wp),                     INTENT(in) :: pzu  !: measurement heigh of bulk wind speed
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pUb  !: bulk wind speed at height pzu m   [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pCdn !: neutral-stability drag coefficient
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: ppsi !: "Psi_m(pzu/L)" stability correction profile for momentum []
      !!----------------------------------------------------------------------------------
      UN10_from_CDN(:,:) = pUb / ( 1._wp + SQRT(pCdn(:,:))/vkarmn * (LOG(pzu/10._wp) - ppsi(:,:)) )
      !!
   END FUNCTION UN10_from_CDN
   !===============================================================================================


   !===============================================================================================
   FUNCTION UN10_from_CD_sclr( pzu, pUb, pCd, ppsi )
      !!----------------------------------------------------------------------------------
      !!  Provides the neutral-stability wind speed at 10 m
      !!----------------------------------------------------------------------------------
      REAL(wp)             :: UN10_from_CD_sclr  !: [m/s]
      REAL(wp),                     INTENT(in) :: pzu  !: measurement heigh of bulk wind speed
      REAL(wp), INTENT(in) :: pUb  !: bulk wind speed at height pzu m   [m/s]
      REAL(wp), INTENT(in) :: pCd  !: drag coefficient
      REAL(wp), INTENT(in) :: ppsi !: "Psi_m(pzu/L)" stability correction profile for momentum []
      !!----------------------------------------------------------------------------------
      !! Reminder: UN10 = u*/vkarmn * log(10/z0)
      !!     and: u* = sqrt(Cd) * Ub
      !!                                  u*/vkarmn * log(   10   /       z0    )
      UN10_from_CD_sclr = SQRT(pCd)*pUb/vkarmn * LOG( 10._wp / z0_from_Cd_sclr( pzu, pCd, ppsi=ppsi ) )
      !!
   END FUNCTION UN10_from_CD_sclr
   !!
   FUNCTION UN10_from_CD_vctr( pzu, pUb, pCd, ppsi )
      REAL(wp), DIMENSION(jpi,jpj)             :: UN10_from_CD_vctr  !: [m/s]
      REAL(wp),                     INTENT(in) :: pzu  !: measurement heigh of bulk wind speed
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pUb  !: bulk wind speed at height pzu m   [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pCd  !: drag coefficient
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: ppsi !: "Psi_m(pzu/L)" stability correction profile for momentum []
      !!----------------------------------------------------------------------------------
      UN10_from_CD_vctr(:,:) = SQRT(pCd(:,:))*pUb/vkarmn * LOG( 10._wp / z0_from_Cd_vctr( pzu, pCd(:,:), ppsi=ppsi(:,:) ) )
      !!
   END FUNCTION UN10_from_CD_vctr
   !===============================================================================================


   !===============================================================================================
   !FUNCTION Re_rough_tq_LKB( iflag, pRer )
   !   !!---------------------------------------------------------------------------------
   !   !!       ***  FUNCTION Re_rough_tq_LKB  ***
   !   !!
   !   !! ** Purpose : returns the "temperature/humidity roughness Reynolds number"
   !   !!              * iflag==1 => temperature => returns: [z_{0t} u*]/Nu_{air}
   !   !!              * iflag==2 => humidity    => returns: [z_{0q} u*]/Nu_{air}
   !   !!              from roughness reynold number "pRer" (i.e. [z_0 u*]/Nu_{air})
   !   !!              between 0 and 1000. Out of range "pRer" indicated by prt=-999.
   !   !!
   !   !!              Based on Liu et al. (1979) JAS 36 1722-1723s
   !   !!
   !   !!              Note: this is what is used into COARE 2.5 to estimate z_{0t} and z_{0q}
   !   !!
   !   !! ** Author: L. Brodeau, April 2020 / AeroBulk (https://github.com/brodeau/aerobulk/)
   !   !!----------------------------------------------------------------------------------
   !   REAL(wp), DIMENSION(jpi,jpj)             :: Re_rough_tq_LKB
   !   INTEGER,                      INTENT(in) :: iflag     !: 1 => dealing with temperature; 2 => dealing with humidity
   !   REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pRer      !: roughness Reynolds number  [z_0 u*]/Nu_{air}
   !   !-------------------------------------------------------------------
   !   ! Scalar Re_r relation from Liu et al.
   !   REAL(wp), DIMENSION(8,2), PARAMETER :: &
   !      & XA = (/ 0.177, 1.376, 1.026, 1.625, 4.661, 34.904, 1667.19, 5.88e5, &
   !      &         0.292, 1.808, 1.393, 1.956, 4.994, 30.709, 1448.68, 2.98e5  /)
   !   !!
   !   REAL(wp), DIMENSION(8,2), PARAMETER :: &
   !      & XB = (/ 0., 0.929, -0.599, -1.018, -1.475, -2.067, -2.907, -3.935, &
   !      &         0., 0.826, -0.528, -0.870, -1.297, -1.845, -2.682, -3.616  /)
   !   !!
   !   REAL(wp), DIMENSION(0:8),   PARAMETER :: &
   !      & XRAN = (/ 0., 0.11, 0.825, 3.0, 10.0, 30.0, 100., 300., 1000. /)
   !   !-------------------------------------------------------------------
   !   ! Scalar Re_r relation from Moana Wave data.
   !   !
   !   !      real*8 A(9,2),B(9,2),RAN(9),pRer,prt
   !   !      integer iflag
   !   !      DATA A/0.177,2.7e3,1.03,1.026,1.625,4.661,34.904,1667.19,5.88E5,
   !   !     &       0.292,3.7e3,1.4,1.393,1.956,4.994,30.709,1448.68,2.98E5/
   !   !      DATA B/0.,4.28,0,-0.599,-1.018,-1.475,-2.067,-2.907,-3.935,
   !   !     &       0.,4.28,0,-0.528,-0.870,-1.297,-1.845,-2.682,-3.616/
   !   !      DATA RAN/0.11,.16,1.00,3.0,10.0,30.0,100.,300.,1000./
   !   !-------------------------------------------------------------------
   !   LOGICAL  :: lfound=.FALSE.
   !   REAL(wp) :: zrr
   !   INTEGER  :: ji, jj, jm
   !   !!----------------------------------------------------------------------------------
   !   Re_rough_tq_LKB(:,:) = -999._wp
   !   !
   !   DO jj = 1, jpj
   !      DO ji = 1, jpi
   !         zrr    = pRer(ji,jj)
   !         lfound = .FALSE.
   !         IF( (zrr > 0.).AND.(zrr < 1000.) ) THEN
   !            jm = 0
   !            DO WHILE ( .NOT. lfound )
   !               jm = jm + 1
   !               lfound = ( (zrr > XRAN(jm-1)) .AND. (zrr <= XRAN(jm)) )
   !            END DO
   !            Re_rough_tq_LKB(ji,jj) = XA(jm,iflag) * zrr**XB(jm,iflag)
   !         END IF
   !      END DO
   !   END DO
   !   !!
   !END FUNCTION Re_rough_tq_LKB
   !===============================================================================================






   !===============================================================================================
   FUNCTION z0tq_LKB( iflag, pRer, pz0 )
      !!---------------------------------------------------------------------------------
      !!       ***  FUNCTION z0tq_LKB  ***
      !!
      !! ** Purpose : returns the "temperature/humidity roughness lengths"
      !!              * iflag==1 => temperature => returns: z_{0t}
      !!              * iflag==2 => humidity    => returns: z_{0q}
      !!              from roughness reynold number "pRer" (i.e. [z_0 u*]/Nu_{air})
      !!              between 0 and 1000. Out of range "pRer" indicated by prt=-999.
      !!              and roughness length (for momentum)
      !!
      !!              Based on Liu et al. (1979) JAS 36 1722-1723s
      !!
      !!              Note: this is what is used into COARE 2.5 to estimate z_{0t} and z_{0q}
      !!
      !! ** Author: L. Brodeau, April 2020 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: z0tq_LKB
      INTEGER,                      INTENT(in) :: iflag     !: 1 => dealing with temperature; 2 => dealing with humidity
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pRer      !: roughness Reynolds number  [z_0 u*]/Nu_{air}
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pz0       !: roughness length (for momentum) [m]
      !-------------------------------------------------------------------
      ! Scalar Re_r relation from Liu et al.
      REAL(wp), DIMENSION(8,2), PARAMETER :: &
         & XA = RESHAPE ( (/ 0.177, 1.376, 1.026, 1.625, 4.661, 34.904, 1667.19, 5.88e5,  &
         &                   0.292, 1.808, 1.393, 1.956, 4.994, 30.709, 1448.68, 2.98e5 /), (/8,2/) )
      !!
      REAL(wp), DIMENSION(8,2), PARAMETER :: &
         & XB = RESHAPE ( (/ 0., 0.929, -0.599, -1.018, -1.475, -2.067, -2.907, -3.935,  &
         &                   0., 0.826, -0.528, -0.870, -1.297, -1.845, -2.682, -3.616 /), (/8,2/) )
      !!
      REAL(wp), DIMENSION(0:8),   PARAMETER :: &
         & XRAN = (/ 0., 0.11, 0.825, 3.0, 10.0, 30.0, 100., 300., 1000. /)
      !-------------------------------------------------------------------
      ! Scalar Re_r relation from Moana Wave data.
      !
      !      real*8 A(9,2),B(9,2),RAN(9),pRer,prt
      !      integer iflag
      !      DATA A/0.177,2.7e3,1.03,1.026,1.625,4.661,34.904,1667.19,5.88E5,
      !     &       0.292,3.7e3,1.4,1.393,1.956,4.994,30.709,1448.68,2.98E5/
      !      DATA B/0.,4.28,0,-0.599,-1.018,-1.475,-2.067,-2.907,-3.935,
      !     &       0.,4.28,0,-0.528,-0.870,-1.297,-1.845,-2.682,-3.616/
      !      DATA RAN/0.11,.16,1.00,3.0,10.0,30.0,100.,300.,1000./
      !-------------------------------------------------------------------
      LOGICAL  :: lfound=.FALSE.
      REAL(wp) :: zrr
      INTEGER  :: ji, jj, jm
      !!---------------------------------------------------------------------------------
      z0tq_LKB(:,:) = -999._wp
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            zrr    = pRer(ji,jj)
            lfound = .FALSE.
            IF( (zrr > 0.).AND.(zrr < 1000.) ) THEN
               jm = 0
               DO WHILE ( .NOT. lfound )
                  jm = jm + 1
                  lfound = ( (zrr > XRAN(jm-1)) .AND. (zrr <= XRAN(jm)) )
               END DO
               z0tq_LKB(ji,jj) = XA(jm,iflag)*zrr**XB(jm,iflag) * pz0(ji,jj)/zrr
            END IF
         END DO
      END DO
      z0tq_LKB(:,:) = MIN( MAX(ABS(z0tq_LKB(:,:)), 1.E-9) , 0.05_wp )
      !!
   END FUNCTION z0tq_LKB
   !===============================================================================================


   !===============================================================================================
   FUNCTION e_air(pqa, pslp)
      !!--------------------------------------------------------------------
      !!                  **** Function e_air ****
      !!
      !! Gives vapour pressure of air from pressure and specific humidity
      !!
      !!--------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: e_air    !: vapour pressure at saturation  [Pa]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pqa      ! specific humidity of air       [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pslp     ! atmospheric pressure            [Pa]
      !!
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: ee, e_old
      REAL(wp) :: zdiff
      ALLOCATE ( ee(jpi,jpj), e_old(jpi,jpj) )
      !
      zdiff  = 1._wp
      e_old = pqa*pslp/reps0
      !
      DO WHILE ( zdiff > repsilon )
         ee = pqa/reps0*(pslp - (1. - reps0)*e_old)
         zdiff  = SUM( ABS( ee - e_old) )
         e_old = ee
      END DO
      e_air = ee
      !
      DEALLOCATE ( ee, e_old )
      !!
   END FUNCTION e_air
   !===============================================================================================


   !===============================================================================================
   FUNCTION rh_air(pqa, pTa, pslp)
      !!----------------------------------------------------------------------------------
      !!                  **** Function e_air ****
      !!
      !! Gives relative humidity of air of air from spec. hum., temperature and pressure
      !!
      !!---------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: rh_air  !: relative humidity [%]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pqa     !: specific humidity of air      [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pTa     !: absolute air temperature               [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pslp    !: atmospheric pressure          [Pa]
      !!----------------------------------------------------------------------------------
      rh_air = e_sat(pTa)
      rh_air = 100. * ( e_air(pqa, pslp) / rh_air )
      !!
   END FUNCTION rh_air
   !===============================================================================================




   !   FUNCTION f_m_louis_sclr( pzu, pt_zu, pq_zu, pwnd, pts, pqs, pCdn, pz0 )
   !      !!----------------------------------------------------------------------------------
   !      !!  Stability correction function for MOMENTUM
   !      !!                 Louis (1979)
   !      !!----------------------------------------------------------------------------------
   !      REAL(wp)             :: f_m_louis_sclr ! term "f_m" in Eq.(6) when option "Louis" rather than "Psi(zeta) is chosen, Lupkes & Gryanik (2015),
   !      REAL(wp),                     INTENT(in) :: pzu     ! reference height (height for pwnd)  [m]
   !      REAL(wp), INTENT(in) :: pt_zu   ! potential air temperature            [Kelvin]
   !      REAL(wp), INTENT(in) :: pq_zu   ! specific air humidity at zt          [kg/kg]
   !      REAL(wp), INTENT(in) :: pwnd    ! wind speed                           [m/s]
   !      REAL(wp), INTENT(in) :: pts     ! sea or ice surface temperature       [K]
   !      REAL(wp), INTENT(in) :: pqs     ! humidity at saturation over ice at T=Ts_i [kg/kg]
   !      REAL(wp), INTENT(in) :: pCdn    ! neutral drag coefficient
   !      REAL(wp), INTENT(in) :: pz0     ! roughness length                      [m]
   !      !!----------------------------------------------------------------------------------
   !      REAL(wp) :: zrib, ztu, zts, zstab
   !      !!----------------------------------------------------------------------------------
   !      zrib = Ri_bulk( pzu, pts, pt_zu, pqs, pq_zu, pwnd ) ! Bulk Richardson Number:
   !      !
   !      zstab = 0.5 + SIGN(0.5_wp, zrib) ; ! Unstable (Ri<0) => zstab = 0 | Stable (Ri>0) => zstab = 1
   !      !
   !      ztu = zrib / ( 1._wp + 3._wp * rc2_louis * pCdn * SQRT( ABS( -zrib * ( pzu / pz0 + 1._wp) ) ) ) ! ABS is just here for when it's stable conditions and ztu is not used anyways
   !      zts = zrib / SQRT( ABS( 1._wp + zrib ) ) ! ABS is just here for when it's UNstable conditions and zts is not used anyways
   !      !
   !      f_m_louis_sclr = (1._wp - zstab) *         ( 1._wp - ram_louis * ztu )  &  ! Unstable Eq.(A6)
   !         &               +      zstab  * 1._wp / ( 1._wp + ram_louis * zts )     ! Stable   Eq.(A7)
   !      !
   !   END FUNCTION f_m_louis_sclr



   !===============================================================================================
   FUNCTION VARIANCE( pvc )
      REAL(4)                            :: VARIANCE
      REAL(wp), DIMENSION(:), INTENT(in) :: pvc
      !!
      INTEGER :: Nn
      REAL(wp) :: zmean
      !!
      Nn = SIZE(pvc)
      !!
      zmean = SUM(pvc)/Nn
      !!
      VARIANCE = REAL( SQRT( SUM( (pvc(:) - zmean) * (pvc(:) - zmean) ) / Nn ) , 4)
      !!
   END FUNCTION VARIANCE
   !===============================================================================================

   !===============================================================================================
   FUNCTION VMEAN( pvc )
      REAL(4)                            :: VMEAN
      REAL(wp), DIMENSION(:), INTENT(in) :: pvc
      !!
      INTEGER :: Nn
      !!
      Nn = SIZE(pvc)
      !!
      VMEAN = SUM(pvc)/Nn
      !!
   END FUNCTION VMEAN
   !===============================================================================================


   !===============================================================================================
   SUBROUTINE TO_KELVIN_3D( pt, cname )
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) :: pt
      CHARACTER(len=*), OPTIONAL, INTENT(in)    :: cname
      INTEGER :: nt
      REAL(wp) :: zm
      CHARACTER(len=32) :: cvar='...'
      !
      IF(PRESENT(cname)) cvar=trim(cname)
      !
      nt = SIZE( pt )
      !
      zm = SUM(pt)/REAL(nt)
      IF( (zm < 50._wp).AND.(zm > -80._wp) ) THEN
         PRINT *, ' *** Variable ', TRIM(cvar), ' is in [deg.C] => converting to [K] !!!'
         pt  = pt + rt0
      ELSEIF ( (zm > 200._wp).AND.(zm < 320._wp) ) THEN
         PRINT *, ' *** Variable ', TRIM(cvar), ' is already in [K], doing nothing...'
      ELSE
         PRINT *, ' *** PROBLEM: cannot figure out unit of variable ', TRIM(cvar), ' !!!'
         STOP
      END IF
   END SUBROUTINE TO_KELVIN_3D
   !===============================================================================================







   SUBROUTINE FIRST_GUESS_COARE( zt, zu, psst, t_zt, pssq, q_zt, U_zu, pcharn, &
      &                          pus, pts, pqs, t_zu, q_zu, Ubzu,  pz0 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  FIRST_GUESS_COARE  ***
      !!
      !! ** Purpose :  Computes fairly accurate first guess of u*, theta* and q*
      !!               by means of the method developed by Fairall et al in the
      !!               COARE family of algorithms
      !!               Purpose is to limit the number of itteration needed...
      !!
      !! INPUT :
      !! -------
      !!    *  zt   : height for temperature and spec. hum. of air            [m]
      !!    *  zu   : height for wind speed (usually 10m)                     [m]
      !!    *  psst  : bulk SST                                                [K]
      !!    *  t_zt : potential air temperature at zt                         [K]
      !!    *  pssq  : SSQ aka saturation specific humidity at temp. psst       [kg/kg]
      !!    *  q_zt : specific humidity of air at zt                          [kg/kg]
      !!    *  U_zu : scalar wind speed at zu                                 [m/s]
      !!    *  pcharn: Charnock parameter
      !!
      !! OUTPUT :
      !! --------
      !!    *  pus    : FIRST GUESS of u* aka friction velocity                            [m/s]
      !!    *  pts    : FIRST GUESS of theta*                                               [K]
      !!    *  pqs    : FIRST GUESS of q* aka friction velocity                            [kg/kg]
      !!    *  t_zu   : FIRST GUESS of pot. air temperature adjusted at wind height zu      [K]
      !!    *  q_zu   : FIRST GUESS of specific humidity of air        //                   [kg/kg]
      !!    *  Ubzu   : FIRST GUESS of bulk wind speed at zu                                [m/s]
      !!
      !! ** Author: L. Brodeau, May 2021 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(in)                     ::   zt
      REAL(wp), INTENT(in)                     ::   zu
      REAL(wp), INTENT(in),  DIMENSION(jpi,jpj) ::   psst
      REAL(wp), INTENT(in),  DIMENSION(jpi,jpj) ::   t_zt
      REAL(wp), INTENT(in),  DIMENSION(jpi,jpj) ::   pssq
      REAL(wp), INTENT(in),  DIMENSION(jpi,jpj) ::   q_zt
      REAL(wp), INTENT(in),  DIMENSION(jpi,jpj) ::   U_zu
      REAL(wp), INTENT(in),  DIMENSION(jpi,jpj) ::   pcharn
      !!
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) ::   pus
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) ::   pts
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) ::   pqs
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) ::   t_zu
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) ::   q_zu
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) ::   Ubzu
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj), OPTIONAL :: pz0    ! roughness length [m]
      !
      INTEGER :: ji, jj
      LOGICAL :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !
      REAL(wp) :: zdt, zdq, zUb, zCd, zus, zts, zqs, zNu_a, zRib
      REAL(wp) :: zlog_zt_o_zu, zlog_10, zlog_zt, zlog_zu, zlog_z0, zlog_z0t, z1_o_sqrt_Cd10, zcc, zcc_ri, z1_o_Ribcu
      REAL(wp) :: zc_a, zc_b
      REAL(wp) :: zz0, zz0t, zprof, zstab, zzeta_u
      REAL(wp) :: ztmp
      !
      REAL(wp), PARAMETER :: zzi0=600._wp, zBeta0=1.2_wp ! COARE values, after all it's a coare method...
      !
      CHARACTER(len=40), PARAMETER :: crtnm = 'FIRST_GUESS_COARE@mod_phymbl.f90'
      !!----------------------------------------------------------------------------------

      l_zt_equal_zu = ( ABS(zu - zt) < 0.01_wp )

      !! First guess of temperature and humidity at height zu:
      t_zu = MAX( t_zt ,  180._wp )   ! who knows what's given on masked-continental regions...
      q_zu = MAX( q_zt , 1.e-6_wp )   !               "

      zz0 = 0.0001

      !! Constants:
      zlog_10 = LOG(10._wp)
      zlog_zt = LOG(zt)
      zlog_zu = LOG(zu)
      zlog_zt_o_zu = LOG(zt/zu)
      zc_a = 0.035_wp*LOG(10._wp/zz0)/LOG(zu/zz0)   !       "                    "               "
      zc_b = 0.004_wp*zzi0*zBeta0*zBeta0*zBeta0

      DO jj = 1, jpj
         DO ji = 1, jpi

            !! Air-sea differences (and we don't want them to be 0...)
            zdt = t_zu(ji,jj) - psst(ji,jj) ;   zdt = SIGN( MAX(ABS(zdt),1.E-6_wp), zdt )
            zdq = q_zu(ji,jj) - pssq(ji,jj) ;   zdq = SIGN( MAX(ABS(zdq),1.E-9_wp), zdq )

            zNu_a = visc_air(t_zu(ji,jj)) ! Air viscosity (m^2/s) at zt given from temperature in (K)

            zUb = SQRT(U_zu(ji,jj)*U_zu(ji,jj) + 0.5_wp*0.5_wp) ! initial guess for wind gustiness contribution

            zus = zc_a*zUb

            !! Roughness length:
            zz0     = pcharn(ji,jj)*zus*zus/grav + 0.11_wp*zNu_a/zus
            zz0     = MIN( MAX(ABS(zz0), 1.E-8) , 1._wp )      ! (prevents FPE from stupid values from masked region later on)
            zlog_z0 = LOG(zz0)

            zCd          = (vkarmn/(zlog_zu - zlog_z0))**2     ! first guess of Cd
            z1_o_sqrt_Cd10 =       (zlog_10 - zlog_z0)/vkarmn  ! sum less costly than product: log(a/b) == log(a) - log(b)

            !! theta,q roughness length:
            zz0t  = 10._wp / EXP( vkarmn/( 0.00115*z1_o_sqrt_Cd10 ) )
            zz0t    = MIN( MAX(ABS(zz0t), 1.E-8) , 1._wp )         ! (prevents FPE from stupid values from masked region later on)
            zlog_z0t = LOG(zz0t)

            zRib = Ri_bulk( zu, psst(ji,jj), t_zu(ji,jj), pssq(ji,jj), q_zu(ji,jj), zUb ) ! Bulk Richardson Number (BRN)

            !! First estimate of zeta_u, depending on the stability, ie sign of BRN (zRib):
            zcc = vkarmn2/(zCd*(zlog_zt - zlog_z0t))
            zcc_ri  = zcc*zRib
            z1_o_Ribcu = -zc_b/zu
            zstab   = 0.5 + SIGN( 0.5_wp , zRib )
            zzeta_u = (1._wp - zstab) *   zcc_ri / (1._wp + zRib*z1_o_Ribcu) & !  Ri_bulk < 0, Unstable
               &  +        zstab      * ( zcc_ri + 27._wp/9._wp*zRib*zRib )    !  Ri_bulk > 0, Stable

            !! First guess M-O stability dependent scaling params.(u*,t*,q*) to estimate z0 and z/L

            !PRINT *, 'LOLO STOP: need to have generic psi functions available into mod_phymbl.f90 !'
            !PRINT *, ' => preferably those of COARE, because it s a COARE first gues...'
            !STOP

            zus  = MAX ( zUb*vkarmn/(zlog_zu - zlog_z0  - psi_m_coare_sclr(zzeta_u)) , 1.E-9 ) ! (MAX => prevents FPE from stupid values from masked region later on)
            ztmp = vkarmn/(zlog_zu - zlog_z0t - psi_h_coare_sclr(zzeta_u))
            zts  = zdt*ztmp
            zqs  = zdq*ztmp

            ! What needs to be done if zt /= zu:
            IF( .NOT. l_zt_equal_zu ) THEN
               !! First update of values at zu (or zt for wind)
               !lolo               zprof = zlog_zt_o_zu + psi_h_ij(zzeta_u) - psi_h_ij(zt*zzeta_u/zu)   ! zt*zzeta_u/zu == zeta_t
               t_zu(ji,jj) = t_zt(ji,jj) - zts/vkarmn*zprof
               q_zu(ji,jj) = q_zt(ji,jj) - zqs/vkarmn*zprof
               q_zu(ji,jj) = (0.5_wp + SIGN(0.5_wp,q_zu(ji,jj)))*q_zu(ji,jj) !Makes it impossible to have negative humidity :
               !
            ENDIF

            pus(ji,jj)  = zus
            pts(ji,jj)  = zts
            pqs(ji,jj)  = zqs
            Ubzu(ji,jj) = zub

            zz0        = pcharn(ji,jj)*zus*zus/grav + 0.11_wp*zNu_a/zus ! LOLO rm !
            PRINT *, 'LOLO: mod_phymbl.f90 end of "FIRST_GUESS_COARE" => z0 =', REAL(zz0,4)

            IF( PRESENT(pz0) ) THEN
               !! Again, because new zus:
               zz0        = pcharn(ji,jj)*zus*zus/grav + 0.11_wp*zNu_a/zus
               pz0(ji,jj) = MIN( MAX(ABS(zz0), 1.E-8) , 1._wp )      ! (prevents FPE from stupid values from masked region later on)
            END IF

         END DO
      END DO

   END SUBROUTINE FIRST_GUESS_COARE



   FUNCTION psi_m_coare_sclr( pzeta )
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
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp) :: psi_m_coare_sclr
      REAL(wp), INTENT(in) :: pzeta
      !!
      REAL(wp) :: zphi_m, zphi_c, zpsi_k, zpsi_c, zf, zc, zstb
      !!----------------------------------------------------------------------------------
      zphi_m = ABS(1. - 15.*pzeta)**.25    !!Kansas unstable
      !
      zpsi_k = 2.*LOG((1. + zphi_m)/2.) + LOG((1. + zphi_m*zphi_m)/2.)   &
         & - 2.*ATAN(zphi_m) + 0.5*rpi
      !
      zphi_c = ABS(1. - 10.15*pzeta)**.3333                   !!Convective
      !
      zpsi_c = 1.5*LOG((1. + zphi_c + zphi_c*zphi_c)/3.) &
         &     - 1.7320508*ATAN((1. + 2.*zphi_c)/1.7320508) + 1.813799447
      !
      zf = pzeta*pzeta
      zf = zf/(1. + zf)
      zc = MIN(50._wp, 0.35_wp*pzeta)
      zstb = 0.5 + SIGN(0.5_wp, pzeta)
      !
      psi_m_coare_sclr = (1. - zstb) * ( (1. - zf)*zpsi_k + zf*zpsi_c ) & ! (pzeta < 0)
         &           -   zstb  * ( 1. + 1.*pzeta     &                ! (pzeta > 0)
         &                          + 0.6667*(pzeta - 14.28)/EXP(zc) + 8.525 )  !     "
      !!
   END FUNCTION psi_m_coare_sclr

   FUNCTION psi_h_coare_sclr( pzeta )
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
      !! Author: L. Brodeau, June 2016 / AeroBulk
      !!         (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------
      REAL(wp) :: psi_h_coare_sclr
      REAL(wp), INTENT(in) :: pzeta
      !!
      REAL(wp) :: zphi_h, zphi_c, zpsi_k, zpsi_c, zf, zc, zstb
      !!----------------------------------------------------------------
      zphi_h = (ABS(1. - 15.*pzeta))**.5  !! Kansas unstable   (zphi_h = zphi_m**2 when unstable, zphi_m when stable)
      !
      zpsi_k = 2.*LOG((1. + zphi_h)/2.)
      !
      zphi_c = (ABS(1. - 34.15*pzeta))**.3333   !! Convective
      !
      zpsi_c = 1.5*LOG((1. + zphi_c + zphi_c*zphi_c)/3.) &
         &    -1.7320508*ATAN((1. + 2.*zphi_c)/1.7320508) + 1.813799447
      !
      zf = pzeta*pzeta
      zf = zf/(1. + zf)
      zc = MIN(50._wp,0.35_wp*pzeta)
      zstb = 0.5 + SIGN(0.5_wp, pzeta)
      !
      psi_h_coare_sclr = (1.-zstb) * ( (1. - zf)*zpsi_k + zf*zpsi_c ) &
         &                  -zstb  * ( (ABS(1. + 2.*pzeta/3.))**1.5     &
         &                            + .6667*(pzeta - 14.28)/EXP(zc) + 8.525 )
      !!
   END FUNCTION psi_h_coare_sclr



   !!
   
   SUBROUTINE check_unit_consitency( cfield, Xval, mask )

      !! Ignore values where mask==0

      CHARACTER(len=*),         INTENT(in) :: cfield
      REAL(wp),   DIMENSION(:,:), INTENT(in) :: Xval
      INTEGER(1), DIMENSION(:,:), INTENT(in) :: mask

      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: lmask
      INTEGER :: nx, ny
      CHARACTER(len=64) :: cunit
      REAL(wp) :: zmean, zmin, zmax
      LOGICAL  :: l_too_large=.FALSE., l_too_small=.FALSE., l_mean_outside=.FALSE.

      nx = SIZE(mask,1)
      ny = SIZE(mask,2)
      ALLOCATE( lmask(nx,ny) )
      lmask(:,:) = .FALSE.
      WHERE( mask==1 ) lmask = .TRUE.

      zmean = SUM( Xval * REAL(mask,wp) ) / SUM( REAL(mask,wp) )

      !PRINT *, 'LOLO, zmean of '//TRIM(cfield)//' =>', zmean

      SELECT CASE (TRIM(cfield))

      CASE('sst')
         zmax = 313._wp
         zmin = 270._wp
         cunit = 'K'

      CASE('t_air')
         zmax = 323._wp
         zmin = 200._wp
         cunit = 'K'

      CASE('q_air')
         zmax = 0.08_wp
         zmin = 0._wp
         cunit = 'kg/kg'

      CASE('slp')
         zmax = 110000.
         zmin =  80000.
         cunit = 'Pa'

      CASE('u10')
         zmax = 50._wp
         zmin =  0._wp
         cunit = 'm/s'

      CASE('v10')
         zmax = 50._wp
         zmin =  0._wp
         cunit = 'm/s'

      CASE('rad_sw')
         zmax = 1500.0_wp
         zmin =  0._wp
         cunit = 'W/m^2'

      CASE('rad_lw')
         zmax = 700.0_wp
         zmin =  0._wp
         cunit = 'W/m^2'

      CASE DEFAULT
         WRITE(*,'(" *** ERROR (mod_aerobulk_compute.f90): we do not know field ",a," !")') TRIM(cfield)
         STOP
      END SELECT

      IF ( MAXVAL(Xval, MASK=lmask) > zmax )                   l_too_large    = .TRUE.
      IF ( MINVAL(Xval, MASK=lmask) < zmin )                   l_too_small    = .TRUE.
      IF ( (zmean < zmin) .OR. (zmean > zmax) ) l_mean_outside = .TRUE.

      IF ( l_too_large .OR. l_too_small .OR. l_mean_outside ) THEN
         WRITE(*,'(" *** ERROR (mod_aerobulk_compute.f90): field ",a," does not seem to be in ",a," !")') TRIM(cfield), TRIM(cunit)
         WRITE(*,'(" min value = ", es10.3," max value = ", es10.3," mean value = ", es10.3)') MINVAL(Xval), MAXVAL(Xval), zmean
         STOP
      END IF

      DEALLOCATE( lmask )
   END SUBROUTINE check_unit_consitency


   FUNCTION type_of_humidity( Xval, mask )
      !!
      !! Guess the type of humidity contained into array Xval
      !! based on the mean, min and max values of the field!
      !! Values where mask==0 are ignored !
      !!
      REAL(wp),   DIMENSION(:,:), INTENT(in) :: Xval
      INTEGER(1), DIMENSION(:,:), INTENT(in) :: mask
      CHARACTER(len=2)                       :: type_of_humidity !: => 'sh', 'rh' or 'dp'

      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: lmask
      INTEGER :: nx, ny
      CHARACTER(len=64) :: cunit
      REAL(wp) :: zmean, zmin, zmax
      LOGICAL  :: l_too_large=.FALSE., l_too_small=.FALSE., l_mean_outside=.FALSE.

      nx = SIZE(mask,1)
      ny = SIZE(mask,2)
      ALLOCATE( lmask(nx,ny) )
      lmask(:,:) = .FALSE.
      WHERE( mask==1 ) lmask = .TRUE.

      zmean =    SUM( Xval * REAL(mask,wp) ) / SUM( REAL(mask,wp) )
      zmin  = MINVAL( Xval,  MASK=lmask )
      zmax  = MAXVAL( Xval,  MASK=lmask )

      !PRINT *, 'LOLO, zmean of '//TRIM(cfield)//' =>', zmean
      type_of_humidity = '00'

      IF( (zmean >= 0._wp).AND.(zmean < 0.08_wp).AND.(zmin >= 0._wp).AND.(zmax < 0.08_wp) ) THEN
         !! Specific humidity! [kg/kg]
         type_of_humidity = 'sh'

      ELSEIF( (zmean >= 150._wp).AND.(zmean < 325._wp).AND.(zmin >= 150._wp).AND.(zmax < 325._wp) ) THEN
         !! Dew-point temperature [K]
         type_of_humidity = 'dp'

      ELSEIF( (zmean >= 0._wp).AND.(zmean <= 100._wp).AND.(zmin >= 0._wp).AND.(zmax <= 100._wp) ) THEN
         !! Relative humidity [%]
         type_of_humidity = 'rh'

      ELSE
         WRITE(6,*) 'ERROR: type_of_humidity()@mod_aerobulk_compute => un-identified humidity type!'
         WRITE(6,*) '   ==> we could not identify the humidity type based on the mean, min & max of the field:'
         WRITE(6,*) '     * mean =', REAL(zmean,4)
         WRITE(6,*) '     * min  =', REAL(zmin, 4)
         WRITE(6,*) '     * max  =', REAL(zmax, 4)
         STOP
      END IF

      DEALLOCATE( lmask )

   END FUNCTION type_of_humidity




END MODULE mod_phymbl
