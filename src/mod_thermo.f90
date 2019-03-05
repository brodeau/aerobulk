! AeroBulk / 2016 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_thermo

   USE mod_const

   IMPLICIT none

   PRIVATE

   PUBLIC :: visc_air, Lvap, e_sat, e_sat_buck, e_air, cp_air, rh_air, &
      &      rho_air, rho_air_adv, q_sat, q_air_rh, q_air_dp, q_sat_simple, &
      &      One_on_L_MO, gamma_moist

   REAL(wp), PARAMETER  :: &
      &      repsilon = 1.e-6

CONTAINS


   FUNCTION visc_air(Ta)

      !! Air viscosity (m^2/s) given from temperature in degrees...

      REAL(wp), DIMENSION(jpi,jpj) :: visc_air
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: Ta  ! air temperature in (K)

      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: TC, TC2

      ALLOCATE ( TC(jpi,jpj), TC2(jpi,jpj) )

      TC  = Ta - rt0   ! air temp, in deg. C
      TC2 = TC*TC

      visc_air = 1.326E-5*(1. + 6.542E-3*TC + 8.301E-6*TC2 - 4.84E-9*TC2*TC)

      DEALLOCATE ( TC, TC2 )

   END FUNCTION visc_air



   FUNCTION Lvap(zsst)

      !: latent heat of vaporization of water from temperature in (K)

      REAL(wp), DIMENSION(jpi,jpj)             :: Lvap   !: [J/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: zsst   !: water temperature [K]

      Lvap = (2.501 - 0.00237*(zsst - rt0))*1.E6

   END FUNCTION Lvap



   FUNCTION e_sat(rT)

      !!**************************************************
      !! rT:     air temperature [K]
      !! e_sat:  water vapor at saturation [Pa]
      !!
      !! Recommended by WMO
      !!
      !! Goff, J. A., 1957: Saturation pressure of water on the new kelvin
      !! temperature scale. Transactions of the American society of heating
      !! and ventilating engineers, 347–354.
      !!
      !!**************************************************

      REAL(wp), DIMENSION(jpi,jpj)             :: e_sat !: vapour pressure at saturation  [Pa]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: rT    !: temperature (K)

      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: ztmp

      ALLOCATE ( ztmp(jpi,jpj) )

      ztmp(:,:) = rtt0/rT(:,:)

      e_sat = 100.*( 10.**(10.79574*(1. - ztmp) - 5.028*LOG10(rT/rtt0)         &
         &       + 1.50475*10.**(-4)*(1. - 10.**(-8.2969*(rT/rtt0 - 1.)) )   &
         &       + 0.42873*10.**(-3)*(10.**(4.76955*(1. - ztmp)) - 1.) + 0.78614) )

      DEALLOCATE ( ztmp )

   END FUNCTION e_sat



   FUNCTION e_sat_buck(rT, slp)

      !!**************************************************
      !!  rT:     air temperature          [K]
      !! slp:     atmospheric pressure     [Pa]
      !! e_sat:  water vapor at saturation [Pa]
      !!
      !! Based on Buck' formula for saturation vapor pressure
      !! from Buck (1981), J. App. Meteor., 1527-1532.
      !!
      !! This version follows the saturation specific humidity computation in
      !! the COARE Fortran code v2.5b.  This results in an increase of ~5% in
      !! latent heat flux compared to the calculation with Teten's
      !!
      !!**************************************************

      REAL(wp), DIMENSION(jpi,jpj)             :: e_sat_buck !: vapour pressure at saturation [Pa]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: rT, &        !: temperature                   [K]
         &                                        slp          !: atmospheric pressure          [Pa]

      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: ztmp

      ALLOCATE ( ztmp(jpi,jpj) )

      !! Achtung: originaly given with temperature in deg.C and pressure in
      !!          millibars! 1 mb = 100 Pa !!!

      ztmp(:,:) = rT(:,:) - rt0

      !! Buck 1981:
      !! Buck, A. L., New equations for computing vapor pressure and enhancement
      !! factor, J. Appl. Meteorol., 20, 1527-1532, 1981
      !e_sat_buck = 611.21 * EXP( 17.502*ztmp/(ztmp + 240.97) )


      !! Official COARE 3.0 code:
      e_sat_buck = 611.2*(1.0007 + 3.46e-8*slp) * EXP( 17.502*ztmp/(ztmp + 240.97) )

      !! Kara et al. 2000:
      !!e_sat_buck = 611.21*(1. + 3.46E-8*slp)*EXP( (17.5*ztmp)/(240.97 + ztmp) )

      DEALLOCATE ( ztmp )

   END FUNCTION e_sat_buck









   FUNCTION e_air(q_air, slp)

      !!--------------------------------------------------------------------
      !!                  **** Function e_air ****
      !!
      !! Gives vapour pressure of air from pressure and specific humidity
      !!
      !!--------------------------------------------------------------------

      REAL(wp), DIMENSION(jpi,jpj)             ::    e_air      !: vapour pressure at saturation  [Pa]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::          &
         &                 q_air,   &  ! specific humidity of air      [kg/kg]
         &                 slp        ! atmospheric pressure          [Pa]

      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: ee, e_old
      REAL :: diff

      ALLOCATE ( ee(jpi,jpj), e_old(jpi,jpj) )

      diff  = 1.
      e_old = q_air*slp/reps0

      DO WHILE ( diff > repsilon )
         ee = q_air/reps0*(slp - (1. - reps0)*e_old)
         diff  = SUM( abs( ee - e_old) )
         e_old = ee
      END DO

      e_air = ee

      DEALLOCATE ( ee, e_old )

   END FUNCTION e_air


   FUNCTION cp_air( pqa )
      !!-------------------------------------------------------------------------------
      !! ** Purpose : Compute specific heat (Cp) of MOIST air
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
      !!-------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pqa     !: air spec. hum. [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj)             :: cp_air  !: [J/K/kg]
      !!-------------------------------------------------------------------------------
      !
      cp_air = Cp_dry + Cp_vap*pqa
      !
   END FUNCTION cp_air


   FUNCTION rh_air(q_air, t_air, slp)

      !! Relative humidity of air

      REAL(wp), DIMENSION(jpi,jpj)             :: rh_air  !: relative humidity [] (fraction!!!, not percent!)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: &
         &                 q_air,   &     !: specific humidity of air      [kg/kg]
         &                 t_air,   &     !: air temperature               [K]
         &                 slp           !: atmospheric pressure          [Pa]

      rh_air = e_sat(t_air)
      rh_air = e_air(q_air,slp) / rh_air

   END FUNCTION rh_air



   FUNCTION q_air_rh(rha, ta, slp)

      !! Specific humidity of air from Relative humidity

      REAL(wp), DIMENSION(jpi,jpj) :: q_air_rh

      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: &
         &     rha,     &   !: relative humidity      [fraction, not %!!!]
         &     ta,      &   !: air temperature        [K]
         &     slp         !: atmospheric pressure          [Pa]

      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: ztmp

      ALLOCATE ( ztmp(jpi,jpj) )

      ztmp       = rha*e_sat(ta)
      q_air_rh = ztmp*reps0/(slp - (1. - reps0)*ztmp)

      DEALLOCATE ( ztmp )

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
      !q_air_dp = e_sat(da)*reps0/(slp - (1. - reps0)*e_sat(da))
      q_air_dp = e_sat(da)*reps0/(slp - (1. - reps0)*e_sat(da))
      !!
   END FUNCTION q_air_dp



   FUNCTION rho_air( ptak, pqa, pslp )
      !!-------------------------------------------------------------------------------
      !! ** Purpose : compute density of (moist) air with eq. of state
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
      !!-------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: ptak, &  !: air temperature   [K]
         &                                        pqa, &  !: air spec. hum.    [kg/kg]
         &                                        pslp    !: pressure in       [Pa]
      REAL(wp), DIMENSION(jpi,jpj)             ::   rho_air   !:              [kg/m^3]
      !!-------------------------------------------------------------------------------
      !
      rho_air = pslp/(R_dry*ptak*(1._wp + rctv0*pqa))
      !
   END FUNCTION rho_air




   FUNCTION rho_air_adv(zt, zq, zP)
      !!
      !! Advanced version, using TRUE virtual temperature
      !!
      REAL(wp), DIMENSION(jpi,jpj) :: rho_air_adv      !: density of air [kg/m^3]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::  &
         &      zt,       &     !: air temperature in (K)
         &      zq,       &     !: air spec. hum. (kg/kg)
         &      zP              !: pressure in       (Pa)
      !!
      REAL(wp), DIMENSION(jpi,jpj) :: ztv !: virtual temperature
      !!
      ztv = zt/(1. - e_air(zq, zP)/zP*(1. - reps0))
      !!
      rho_air_adv = zP/(R_dry*ztv)
      !!
   END FUNCTION rho_air_adv


   FUNCTION q_sat(temp, slp,  cform)

      !! Specific humidity at saturation

      REAL(wp), DIMENSION(jpi,jpj) :: q_sat
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::  &
         &                  temp,  &   !: sea surface temperature         [K]
         &                  slp       !: sea level atmospheric pressure  [Pa]

      CHARACTER(len=*), OPTIONAL, INTENT(in) :: cform


      !! Local :
      LOGICAL :: lbuck  !: we use Buck formula to compute e_sat instead of Goff 1957
      REAL(wp), DIMENSION(jpi,jpj) ::  &
         &    e_s

      lbuck = .FALSE.
      IF ( PRESENT(cform) ) THEN
         IF ( (TRIM(cform) == 'buck').OR.(TRIM(cform) == 'Buck').OR.(TRIM(cform) == 'BUCK') ) THEN
            lbuck = .TRUE.
         END IF
      END IF

      !! Vapour pressure at saturation :
      IF ( lbuck ) THEN
         e_s = e_sat_buck(temp, slp)
      ELSE
         e_s = e_sat(temp)  ! using Goff !
      END IF

      q_sat = reps0*e_s/(slp - (1. - reps0)*e_s)

   END FUNCTION q_sat




   FUNCTION e_sat_ice(rt)

      REAL(wp), DIMENSION(jpi,jpj) :: e_sat_ice !: vapour pressure at saturation in presence of ice [Pa]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: rt

      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: ztmp

      ALLOCATE ( ztmp(jpi,jpj) )

      ztmp(:,:) = 273.16/rt(:,:)

      e_sat_ice = 100.*(10**( -9.09718*(ztmp - 1.) - 3.56654*LOG10(ztmp) &
         &                + 0.876793*(1. - rt/273.16) + LOG10(6.1071) ) )

      DEALLOCATE ( ztmp )

   END FUNCTION e_sat_ice



   FUNCTION q_sat_simple(temp, zrho)

      REAL(wp), DIMENSION(jpi,jpj)             :: q_sat_simple
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: &
         &       temp,   &  !: sea surface temperature  [K]
         &       zrho       !: air density         [kg/m^3]

      q_sat_simple = 640380./zrho * exp(-5107.4/temp)

   END FUNCTION q_sat_simple



   FUNCTION One_on_L_MO(theta_a, q_a, us, ts, qs)

      !! ********************************************************************************
      !! Evaluates the 1./(Monin Obukhov length) from average temperature, specific humidity
      !! and frictional u, t and q
      !! 2015: L. Brodeau
      !! ********************************************************************************

      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: theta_a,  &  !: average potetntial air temperature [K]
         &                                q_a,      &  !: average specific humidity of air   [kg/kg]
         &                                us, ts, qs   !: frictional velocity, temperature and humidity

      REAL(wp), DIMENSION(jpi,jpj)             :: One_on_L_MO         !: 1./(Monin Obukhov length) [m^-1]

      REAL(wp), DIMENSION(:,:), ALLOCATABLE  :: rqa

      ALLOCATE ( rqa(jpi,jpj) )

      rqa = (1. + rctv0*q_a)

      One_on_L_MO =  grav*vkarmn*(ts*rqa + rctv0*theta_a*qs) / ( us*us * theta_a*rqa )

      DEALLOCATE ( rqa )

   END FUNCTION One_on_L_MO



   FUNCTION gamma_moist( ptak, pqa )
      !!----------------------------------------------------------------------------------
      !! ** Purpose : Compute the moist adiabatic lapse-rate.
      !!     => http://glossary.ametsoc.org/wiki/Moist-adiabatic_lapse_rate
      !!     => http://www.geog.ucsb.edu/~joel/g266_s10/lecture_notes/chapt03/oh10_3_01/oh10_3_01.html
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: gamma_moist
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: ptak, pqa ! air temperature (K) and specific humidity (kg/kg)
      !
      INTEGER  ::   ji, jj         ! dummy loop indices
      REAL(wp) :: zrv, ziRT        ! local scalar
      !!----------------------------------------------------------------------------------
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            zrv = pqa(ji,jj) / (1. - pqa(ji,jj))
            ziRT = 1./(R_dry*ptak(ji,jj))    ! 1/RT
            gamma_moist(ji,jj) = grav * ( 1. + L0vap*zrv*ziRT ) / ( Cp_dry + L0vap*L0vap*zrv*reps0*ziRT/ptak(ji,jj) )
         END DO
      END DO
      !
   END FUNCTION gamma_moist

END MODULE mod_thermo
