! AeroBulk / 2015 / L. Brodeau (brodeau@gmail.com)
! https://sourceforge.net/p/aerobulk

PROGRAM TEST_COEF

   USE mod_const
   USE mod_thermo

   USE mod_blk_coare

   IMPLICIT NONE

   INTEGER, PARAMETER :: nb_algos = 1

   !CHARACTER(len=10), DIMENSION(nb_algos), PARAMETER :: vca = (/ 'coare', 'ncar ', 'ecmwf' /)
   CHARACTER(len=10), DIMENSION(nb_algos), PARAMETER :: vca = (/ 'coare' /)

   REAL(4), DIMENSION(nb_algos) :: vCd, vCe, vCh, vTheta_u, vT_u, vQu, vz0, vus, vRho_u, vUg

   REAL(wp), PARAMETER ::   &
      &   zu  = 10.

   INTEGER, PARAMETER ::   &
      &   n_dt   = 21,   &
      &   n_u    = 30,   &
      &   n_q    = 10

   CHARACTER(len=2) :: car

   CHARACTER(len=100) :: &
      &   calgob


   INTEGER :: jarg, ialgo

   INTEGER, PARAMETER :: lx=1, ly=1
   REAL(wp),    DIMENSION(lx,ly) :: Ublk, Ts, ssq_s, zz0, zus

   REAL(wp), DIMENSION(lx,ly) :: sst, qsat_zt, SLP, &
      &  W10, t_zt, q_zt, RH_zt, t_zu, q_zu, qsat_sst, rho_zu, tmp, rad_lw, rad_sw


   REAL(wp), DIMENSION(lx,ly) :: Cd, Ce, Ch, Cp_ma, rgamma

   REAL(wp) :: zt, nu_air

   CHARACTER(len=3) :: czt, czu

   LOGICAL :: l_ask_for_slp = .FALSE. , &  !: ask for SLP, otherwize assume SLP = 1010 hPa
      &     l_use_rh      = .FALSE.      !: ask for RH rather than q for humidity

   jpi = lx ; jpj = ly


   nb_itt = 20  ! 20 itterations in bulk algorithm...


   jarg = 0

   DO WHILE ( jarg < command_argument_count() )

      jarg = jarg + 1
      CALL get_command_argument(jarg,car)

      SELECT CASE (trim(car))

      CASE('-h')
         call usage_test()

      CASE('-p')
         !PRINT *, 'Will ask for SLP'
         l_ask_for_slp = .TRUE.

      CASE('-r')
         !PRINT *, "Using RH instead of q!!!"
         l_use_rh = .TRUE.

      CASE DEFAULT
         PRINT *, 'Unknown option: ', trim(car) ; PRINT *, ''
         CALL usage_test()

      END SELECT

   END DO



   PRINT *, ''
   PRINT *, '  *** Epsilon            = Rd/Rv       (~0.622) =>', reps0
   PRINT *, '  *** Virt. temp. const. = (1-eps)/eps (~0.608) =>', rctv0
   PRINT *, ''





   PRINT *, ''


   PRINT *, 'Give "zt", altitude for air temp. and humidity in meters (generally 2 or 10):'
   READ(*,*) zt
   PRINT *, ''


   IF ( (zt > 99.).OR.(zu > 99.) ) THEN
      PRINT *, 'Be reasonable in your choice of zt or zu, they should not exceed a few tenths of meters!' ; STOP
   END IF
   IF ( zt < 10. ) THEN
      WRITE(czt,'(i1,"m")') INT(zt)
   ELSE
      WRITE(czt,'(i2,"m")') INT(zt)
   END IF
   WRITE(czu,'(i2,"m")') INT(zu)


   IF ( l_ask_for_slp ) THEN
      PRINT *, 'Give sea-level pressure (hPa):'
      READ(*,*) SLP
      SLP = SLP*100.
   ELSE
      SLP = Patm
      PRINT *, 'Using a sea-level pressure of ', Patm
   END IF
   PRINT *, ''

   PRINT *, 'Give SST (deg. C):'
   READ(*,*) sst
   sst = sst + rt0
   PRINT *, 'For this sst the latent heat of vaporization is Lvap =', Lvap(sst), ' [J/kg]'
   PRINT *, ''

   PRINT *, 'Give temperature at ',trim(czt),' (deg. C):'
   READ(*,*) t_zt
   t_zt = t_zt + rt0
   PRINT *, ''


   !! Asking for humidity:
   qsat_zt = q_sat(t_zt, SLP)  ! spec. hum. at saturation [kg/kg]

   IF ( l_use_rh ) THEN
      PRINT *, 'Give relative humidity at ',trim(czt),' (%):'
      READ(*,*) RH_zt
      RH_zt = 1.E-2*RH_zt
      q_zt = q_air_rh(RH_zt, t_zt, SLP)
      PRINT *, 'q_',TRIM(czt),' from RH_',TRIM(czt),' =>', 1000*q_zt, ' [g/kg]'
      !PRINT *, 'Inverse => RH from q_zt:', 100*rh_air(q_zt, t_zt, SLP)
   ELSE
      WRITE(*, '("Give specific humidity at ",a," (g/kg) (saturation is at ",f6.3," g/kg):")') &
         &      TRIM(czt), 1000.*qsat_zt
      READ(*,*) q_zt
      q_zt = 1.E-3*q_zt
      RH_zt   = rh_air(q_zt, t_zt, SLP)
      WRITE(*,'("  => Relative humidity at ",a," = ",f4.1,"%")') TRIM(czt), 100*RH_zt
      !PRINT *, 'Inverse => q_zt from RH :', 1000*q_air_rh(RH_zt, t_zt, SLP)
   END IF
   PRINT *, ''

   IF ( q_zt(1,1) > qsat_zt(1,1) ) THEN
      PRINT *, ' ERROR: you can not go belong saturation!!!' ; STOP
   END IF



   PRINT *, ''
   PRINT *, '===================================================================='
   PRINT *, ' *** density of air at ',TRIM(czt),' => ',  rho_air(t_zt, q_zt, SLP), '[kg/m^3]'

   Cp_ma = cp_air(q_zt)
   PRINT *, ' *** Cp of (moist) air at ',TRIM(czt),' => ', Cp_ma, '[J/K/kg]'
   PRINT *, ''
   rgamma = gamma_moist(t_zt, q_zt)
   PRINT *, ' *** Adiabatic lapse-rate of (moist) air at ',TRIM(czt),' => ', REAL(1000.*rgamma ,4), '[K/1000m]'
   PRINT *, '===================================================================='
   PRINT *, ''
   PRINT *, ''

   qsat_sst = 0.98*q_sat(sst, SLP)

   PRINT *, '' ;  PRINT *, 'BEFORE: SSQ = 0.98*q_sat(sst) =', 1000.*qsat_sst, '[g/kg]' ; PRINT *, ''


   !! Must give something more like a potential temperature at zt:
   !! UPDATING t_zt !!!
   t_zt = t_zt + rgamma*zt

   PRINT *, ''
   PRINT *, 'Pot. temp. at ',TRIM(czt),' (using gamma)  =', t_zt - rt0, ' [deg.C]'



   !! Checking the difference of virtual potential temperature between air at zt and sea surface:
   tmp = t_zt*(1. + rctv0*q_zt)
   PRINT *, 'Virtual pot. temp. at ',TRIM(czt),'   =', REAL(tmp - rt0 , 4), ' [deg.C]'
   PRINT *, 'Pot. temp. diff. air/sea at ',TRIM(czt),' =', REAL(t_zt - sst , 4), ' [deg.C]'
   PRINT *, 'Virt. pot. temp. diff. air/sea at ',TRIM(czt),' =', REAL(tmp - sst*(1. + rctv0*qsat_sst) , 4), ' [deg.C]'
   PRINT *, ''; PRINT *, ''


   PRINT *, 'Give wind speed at zu=10m (m/s):'
   READ(*,*) W10
   PRINT *, ''

   PRINT *, 'Give downwelling shortwave radiation (W/m^2):'
   READ(*,*) rad_sw
   PRINT *, ''

   PRINT *, 'Give downwelling longwave radiation (W/m^2):'
   READ(*,*) rad_lw
   PRINT *, ''




   !DO ialgo = 1, nb_algos
   ialgo = 1

   calgob = trim(vca(ialgo))

   zz0 = 0.
   zus = 0.
   
   Ts = sst
   
   CALL turb_coare( '3.0', zt, zu, Ts, t_zt, ssq_s, q_zt, W10, &
      &             Cd, Ch, Ce, t_zu, q_zu, Ublk, &
      &             rad_sw=rad_sw, rad_lw=rad_lw, slp=SLP, &
      &             xz0=zz0, xu_star=zus )

   

   PRINT *, ''
   PRINT *, ' Skin correction on SST =>', Ts - sst
   PRINT *, ' Skin correction on SSQ =>', 1000.*(ssq_s - qsat_sst)
   PRINT *, ''
   STOP


   !CASE DEFAULT
   !   write(6,*) 'Bulk algorithm #', ialgo, ' is unknown!!!' ; STOP
   !END SELECT



   vTheta_u(ialgo) = REAL(   t_zu(1,1) -rt0 , 4)   ! Potential temp.

   t_zu = t_zu - rgamma*zu   ! Real temp.

   vCd(ialgo) = REAL(1000.*Cd(1,1) ,4)
   vCh(ialgo) = REAL(1000.*Ch(1,1) ,4)
   vCe(ialgo) = REAL(1000.*Ce(1,1) ,4)

   vT_u(ialgo) = REAL( t_zu(1,1) -rt0 , 4)    ! Real temp.
   vQu(ialgo) = REAL(  q_zu(1,1) , 4)

   !! Air density at zu (10m)
   rho_zu = rho_air(t_zu, q_zu, SLP)
   tmp = SLP - rho_zu*grav*zu
   rho_zu = rho_air(t_zu, q_zu, tmp)
   vRho_u(ialgo) = REAL(rho_zu(1,1) ,4)

   !! Gustiness contribution:
   vUg(ialgo) = REAL(Ublk(1,1)-W10(1,1) , 4)

   !! z0 et u*:
   vz0(ialgo) = REAL(zz0(1,1) ,4)
   vus(ialgo) = REAL(zus(1,1) ,4)

   !END DO


   PRINT *, ''; PRINT *, ''


   IF ( zt < zu ) THEN
      PRINT *, ''; PRINT *, 'Temperature and humidity at z = ',trim(czt),' :'
      PRINT *, 't_',TRIM(czt),'  =', t_zt-rt0 ;  PRINT *, 'q_',TRIM(czt),'  =', q_zt
   END IF

   PRINT *, ''; PRINT *, 'Temperatures and humidity at z = ',trim(czu),' :'
   PRINT *, '========================================================================='
   PRINT *, '  Algorithm:           ',trim(vca(1)),'    |    '!,trim(vca(2)),'       |    ',trim(vca(3))
   PRINT *, '========================================================================='
   PRINT *, '    theta_',TRIM(czu),' =   ', vTheta_u       , '[deg.C]'
   PRINT *, '    t_',TRIM(czu),'     =   ', vT_u      , '[deg.C]'
   PRINT *, '    q_',TRIM(czu),'     =   ', REAL(1000.*vQu ,4)  , '[g/kg]'
   PRINT *, ''
   PRINT *, '      SSQ     =   ', REAL(1000.*qsat_sst(1,1), 4)  , '[g/kg]'
   PRINT *, '    Delta t   =   ', REAL(vT_u  - (sst(1,1)-rt0) , 4)      , '[deg.C]'
   PRINT *, '    Delta q   =   ', REAL(1000.*(vQu - qsat_sst(1,1)), 4)  , '[g/kg]'
   PRINT *, ''
   PRINT *, '    Ug (gust) =   ', vUg , '[m/s]'
   PRINT *, ''



   tmp = visc_air(t_zu)
   nu_air = tmp(1,1)


   PRINT *, ''
   PRINT *, 'Kinematic viscosity of air =', nu_air
   PRINT *, ''
   PRINT *, 'With a pressure of', int(SLP)
   PRINT *, 'Density of air at ',TRIM(czu),' =', REAL(vRho_u,4), ' [kg/m^3]'
   IF ( zt < zu )  PRINT *, ' (density of air at ',TRIM(czt),' was ', REAL(rho_air(t_zt, q_zt, SLP),4),')'
   PRINT *, ''
   Cp_ma = cp_air(q_zu)
   PRINT *, ' Cp of moist air at ',TRIM(czu),' => ', REAL(Cp_ma,4), ' [J/K/kg]'
   PRINT *, ''


   PRINT *, ''
   PRINT *, '   *** Bulk Transfer Coefficients:'
   PRINT *, '========================================================================='
   PRINT *, '  Algorithm:           ',trim(vca(1)),'    |    '!,trim(vca(2)),'       |    ',trim(vca(3))
   PRINT *, '========================================================================='
   PRINT *, '      C_D     =   ', vCd        , '[10^-3]'
   PRINT *, '      C_E     =   ', vCe        , '[10^-3]'
   PRINT *, '      C_H     =   ', vCh        , '[10^-3]'
   PRINT *, ''
   PRINT *, '      z_0     =   ', vz0        , '[m]'
   PRINT *, '      u*      =   ', vus        , '[m/s]'
   PRINT *, 'Equ. Charn p. =   ', REAL(vz0*grav/(vus*vus) -0.11*grav*nu_air/(vus*vus*vus),4)
   PRINT *, ''


   PRINT *, ''

END PROGRAM TEST_COEF











SUBROUTINE usage_test()
   !!
   PRINT *,''
   PRINT *,'   List of command line options:'
   PRINT *,''
   PRINT *,' -p   => ask for sea-level pressure, otherwize assume 1010 hPa'
   PRINT *,''
   PRINT *,' -r   => Ask for relative humidity rather than specific humidity'
   PRINT *,''
   PRINT *,' -h   => Show this message'
   PRINT *,''
   STOP
   !!
END SUBROUTINE usage_test
!!
