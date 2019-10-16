! AeroBulk / 2015 / L. Brodeau

MODULE mod_const

   IMPLICIT NONE

   PUBLIC

   !! Following NEMO:
   INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND( 6, 37)   !: single precision (real 4)
   INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)   !: double precision (real 8)
   INTEGER, PARAMETER :: wp = dp                           !: working precision

   ! THINGS THAT NEED TO BE GIVEN A VALUE, anh have the same name and type as in NEMO...
   ! Stupid values here to prevent to ommit to give them a value:
   INTEGER,                  SAVE :: jpi, jpj    !: 2D dimensions of array to be used in AeroBulk
   INTEGER,  PARAMETER            :: jpk = 1
   REAL(wp), DIMENSION(jpk), SAVE :: gdept_1d = (/ 1._wp /) !: depth at which SST is measured [m]
   REAL(wp),                 SAVE :: rdt = 3600. !: time step for the cool-skin/warm-layer parameterization  [s]
   INTEGER,                  SAVE :: nb_itt=5  !: number of itteration in the bulk algorithm
   
   LOGICAL, SAVE :: l_first_call=.true. , l_last_call=.false.

   LOGICAL, PARAMETER :: ldebug_blk_algos=.false.



   
   REAL(wp), PARAMETER, PUBLIC :: &
      &   rpi  = 3.141592653589793, &
      &   rt0  = 273.15,    &      !: freezing point of fresh water [K]
      &   rtt0 = 273.16,    &      !: riple point of temperature    [K]
      &  grav  = 9.8,       &      !: acceleration of gravity    [m.s^-2]
      &  Patm  = 101000.,   &
      !!
      & rho0_a = 1.2  ,     &  !: Approx. of density of air                     [kg/m^3]
      & rCp0_a  = 1015.0 ,  &  !: Specic heat of moist air                      [J/K/kg]
      & rCp_dry = 1005.0 ,  &  !: Specic heat of dry air, constant pressure      [J/K/kg]
      & rCp_vap = 1860.0 ,  &  !: Specic heat of water vapor, constant pressure  [J/K/kg]
      !!
      &  R_dry = 287.05,    &  !: Specific gas constant for dry air              [J/K/kg]
      &  R_vap = 461.495,   &  !: Specific gas constant for water vapor          [J/K/kg]
      !!
      & rCp0_w  = 4190. ,    &  !: Specific heat capacity of seawater (ECMWF 4190) [J/K/kg]
      & rho0_w = 1025.  ,    &  !: Density of sea-water  (ECMWF->1025)             [kg/m^3]
      &  rnu0_w = 1.e-6,      &  !: kinetic viscosity of water                      [m^2/s]
      &  rk0_w  = 0.6,        &  !: thermal conductivity of water (at 20C)          [W/m/K]
      !!
      &  reps0 = R_dry/R_vap,  &   !: ratio of gas constant for dry air and water vapor => ~ 0.622
      !!
      &  rctv0 = R_vap/R_dry - 1. , &   !: for virtual temperature (== (1-eps)/eps) => ~ 0.608
      !!
      &  rnu0_air  = 1.5E-5,   &   !: kinematic viscosity of air    [m^2/s]
      !!
      &  rLevap = 2.46E6,     &   !: Latent heat of vaporization for sea-water in J/kg
      &  vkarmn = 0.4,        &   !: Von Karman's constant
      &  Pi    = 3.141592654, &
      &  twoPi = 2.*Pi,       &
      &  emiss_w = 0.97,      &   !: emissivity of sea water
      &  stefan = 5.67E-8,    &   !: Stefan Boltzman constant
      !!
      &  oce_alb0  = 0.066,  &   !: Default sea surface albedo over ocean when nothing better is available      
      !!                         !: NEMO: 0.066 / ECMWF: 0.055
      &  Tswf  = 273.,  &        !: BAD!!! because sea-ice not used yet!!!
      !&  Tswf  = 271.4          !: freezing point of sea-water (K)
      &  to_rad = Pi/180., &
      &  R_earth = 6.37E6,        & ! Earth radius (m)
      &  rtilt_earth = 23.5, &
      &  Sol0 = 1366.  , &       !: Solar constant W/m^2
      &  rdct_qsat_salt = 0.98   !: factor to apply to q_sat(SST) to account for salt

   REAL(wp), PARAMETER, PUBLIC :: radrw    = rho0_a/rho0_w !: Density ratio
   REAL(wp), PARAMETER, PUBLIC :: sq_radrw = SQRT(rho0_a/rho0_w)


   
   INTEGER, DIMENSION(12), PARAMETER, PUBLIC :: &
      &   tdmn = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /), &
      &   tdml = (/ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)



   !! IFS:

   !REAL(wp), PARAMETER :: &
   !&    RKBOL = 1.380658E-23, &
   !&    RNAVO = 6.0221367E+23, &
   !&    R = RNAVO*RKBOL, &
   !&    RMD = 28.9644, &
   !&    RMV = 18.0153, &
   !&    RD = 1000.*R/RMD, &
   !&    RV = 1000.*R/RMV, &
   !&    rCp_dry = 3.5*RD, &
   !&    rLevap=2.5008E+6, &
    !&    rctv0 = RV/RD-1.0







   

   !! Max and min values for variable (for netcdf files)

   !REAL(wp), PARAMETER :: &
   !&                qlat_min = -1200., qlat_max = 200., &
   !     &                qsen_min = -900. , qsen_max = 300., &
   !     &                taum_min =    0., taum_max =   4.,  &
   !     !!
   !     &                msl_min  = 85000., msl_max = 105000., &
   !     &                tair_min =  223., tair_max = 323.,  &
   !     &                qair_min =    0., qair_max = 0.03,  &
   !     &                w10_min  =    0., w10_max  = 45.,   &
   !     &                sst_min  =  270., sst_max  = 310.,  &
   !     &                ice_min  =    0., ice_max  =   1.,  &
   !     &                cx_min   =    0., cx_max   = 0.01,  &
   !     &                rho_min  =   0.8, rho_max  =  1.5       ! density of air
   !CCSM/ccsm2/models/ice/dice5/flux_ai.F90



END MODULE mod_const




!! In CCSM:
!! from CCSM/ccsm2/models/ice/dice5/flux_ai.F90 :


!      REAL, PARAMETER :: cpair  =  1.005e3         !: Specific heat of dry air
!      REAL, PARAMETER :: cpwv   =  1.81e3          !: Specific heat of water vapor
!      REAL, PARAMETER :: cpvir  =  cpwv/cpair - 1. !: Defined as cpwv/cpair - 1.
!      REAL, PARAMETER :: gravit =  9.80616         !: Acceleration of gravity
!      REAL, PARAMETER :: stebol = 56.7e-9          !: Stefan-Boltzmann's constant
!      REAL, PARAMETER :: xkar   =  0.4             !: Von Karman constant
!      REAL, PARAMETER :: zvir   =  0.606           !: rh2o/rair - 1.0
!      REAL, PARAMETER :: zzsice =  0.0005          !: ice surface roughness
!      REAL, PARAMETER :: latvap =  2.5e6           !: latent heat of evaporation
!      REAL, PARAMETER :: latice =   .334e6         !: latent heat of fusion
