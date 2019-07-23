! AeroBulk / 2015 / L. Brodeau

PROGRAM TEST_AEROBULK

   USE mod_const
   USE mod_phymbl

   USE mod_blk_coare3p0
   USE mod_blk_coare3p6
   USE mod_blk_ncar
   USE mod_blk_ecmwf

   IMPLICIT NONE

   INTEGER, PARAMETER :: nb_algos = 4

   CHARACTER(len=8), DIMENSION(nb_algos), PARAMETER :: &
      &      vca = (/ 'coare3p0', 'coare3p6', 'ncar    ', 'ecmwf   ' /)

   REAL(4), DIMENSION(nb_algos) ::  &
      &           vCd, vCe, vCh, vTheta_u, vT_u, vQu, vz0, vus, vRho_u, vUg, vL, vBRN, &
      &           vUN10, vQL, vTau, vQH, vEvap, vTs, vqs

   REAL(wp), PARAMETER ::   &
      & to_mm_p_day = 24.*3600.  !: freshwater flux: from kg/s/m^2 == mm/s to mm/day

   INTEGER, PARAMETER ::   &
      &   n_dt   = 21,   &
      &   n_u    = 30,   &
      &   n_q    = 10

   CHARACTER(len=2) :: car

   CHARACTER(len=100) :: &
      &   calgob

   INTEGER :: jarg, ialgo, jq

   INTEGER, PARAMETER :: lx=1, ly=1
   REAL(wp),    DIMENSION(lx,ly) :: Ublk, zz0, zus, zts, zqs, zL, zUN10

   REAL(wp), DIMENSION(lx,ly) :: sst, Ts, qsat_zt, SLP, &
      &  W10, t_zt, theta_zt, q_zt, RH_zt, d_zt, t_zu, theta_zu, q_zu, ssq, qs, rho_zu, rad_sw, rad_lw, &
      &  tmp

   REAL(wp), DIMENSION(lx,ly) :: Cd, Ce, Ch, Cp_ma, rgamma

   REAL(wp) :: zt, zu, nu_air

   CHARACTER(len=3) :: czt, czu

   LOGICAL :: l_ask_for_slp = .FALSE. , &  !: ask for SLP, otherwize assume SLP = 1010 hPa
      &     l_use_rh      = .FALSE. ,   &  !: ask for RH rather than q for humidity
      &     l_use_dp      = .FALSE. ,   &  !: ask for dew-point temperature rather than q for humidity
      &     l_use_cswl    = .FALSE.        !: compute and use the skin temperature
   !                                       !: (Cool Skin Warm Layer parameterization)
   !                                       !:  => only in COARE and ECMWF

   jpi = lx ; jpj = ly

   nb_itt = 20  ! 20 itterations in bulk algorithm...

   OPEN(6, FORM='formatted', RECL=512)

   jarg = 0

   DO WHILE ( jarg < command_argument_count() )

      jarg = jarg + 1
      CALL get_command_argument(jarg,car)

      SELECT CASE (TRIM(car))

      CASE('-h')
         call usage_test()

      CASE('-p')
         l_ask_for_slp = .TRUE.

      CASE('-r')
         l_use_rh = .TRUE.

      CASE('-d')
         l_use_dp = .TRUE.

      CASE('-S')
         l_use_cswl = .TRUE.

      CASE DEFAULT
         WRITE(6,*) 'Unknown option: ', TRIM(car) ; WRITE(6,*) ''
         CALL usage_test()

      END SELECT

   END DO



   WRITE(6,*) ''
   WRITE(6,*) '  *** Epsilon aka reps0  = Rd/Rv       (~0.622) =>', reps0
   WRITE(6,*) '  *** Virt. temp. const. = (1-eps)/eps (~0.608) =>', rctv0
   WRITE(6,*) ''





   WRITE(6,*) ''

   WRITE(6,*) 'Give "zu", height of wind speed measurement in meters (generally 10m):'
   READ(*,*) zu
   WRITE(6,*) ''


   WRITE(6,*) 'Give "zt", height of air temp. and humidity measurement in meters (generally 2 or 10m):'
   READ(*,*) zt
   WRITE(6,*) ''


   IF ( (zt > 99.).OR.(zu > 99.) ) THEN
      WRITE(6,*) 'Be reasonable in your choice of zt or zu, they should not exceed a few tenths of meters!' ; STOP
   END IF
   IF ( zt < 10. ) THEN
      WRITE(czt,'(i1,"m")') INT(zt)
   ELSE
      WRITE(czt,'(i2,"m")') INT(zt)
   END IF
   WRITE(czu,'(i2,"m")') INT(zu)


   IF ( l_ask_for_slp ) THEN
      WRITE(6,*) 'Give sea-level pressure (hPa):'
      READ(*,*) SLP
      SLP = SLP*100.
   ELSE
      SLP = Patm
      WRITE(6,*) 'Using a sea-level pressure of ', Patm
   END IF
   WRITE(6,*) ''

   WRITE(6,*) 'Give SST (deg. C):'
   READ(*,*) sst
   sst = sst + rt0
   WRITE(6,*) 'For this sst the latent heat of vaporization is L_vap =', L_vap(sst), ' [J/kg]'
   WRITE(6,*) ''

   WRITE(6,*) 'Give temperature at ',TRIM(czt),' (deg. C):'
   READ(*,*) t_zt
   t_zt = t_zt + rt0
   WRITE(6,*) ''


   !! Asking for humidity:
   qsat_zt = q_sat(t_zt, SLP)  ! spec. hum. at saturation [kg/kg]

   IF ( l_use_rh ) THEN
      WRITE(6,*) 'Give relative humidity at ',TRIM(czt),' [%]:'
      READ(*,*) RH_zt
      RH_zt = 1.E-2*RH_zt
      q_zt = q_air_rh(RH_zt, t_zt, SLP)
      WRITE(6,*) 'q_',TRIM(czt),' from RH_',TRIM(czt),' =>', 1000*q_zt, ' [g/kg]'
      !WRITE(6,*) 'Inverse => RH from q_zt:', 100*rh_air(q_zt, t_zt, SLP)
   ELSEIF ( l_use_dp ) THEN
      WRITE(6,*) 'Give dew-point temperature at ',TRIM(czt),' (deg. C):'
      READ(*,*) d_zt
      d_zt = d_zt + rt0
      q_zt = q_air_dp(d_zt, SLP)
      WRITE(6,*) 'q_',TRIM(czt),' from d_',TRIM(czt),' =>', 1000*q_zt, ' [g/kg]'
      !WRITE(6,*) 'Inverse => RH from q_zt:', 100*rh_air(q_zt, t_zt, SLP)
   ELSE
      WRITE(*, '("Give specific humidity at ",a," (g/kg) (saturation is at ",f6.3," g/kg):")') &
         &      TRIM(czt), 1000.*qsat_zt
      READ(*,*) q_zt
      q_zt = 1.E-3*q_zt
      RH_zt   = rh_air(q_zt, t_zt, SLP)
      WRITE(*,'("  => Relative humidity at ",a," = ",f4.1,"%")') TRIM(czt), 100*RH_zt
      !WRITE(6,*) 'Inverse => q_zt from RH :', 1000*q_air_rh(RH_zt, t_zt, SLP)
   END IF
   WRITE(6,*) ''

   IF ( q_zt(1,1) > qsat_zt(1,1) ) THEN
      WRITE(6,*) ' ERROR: you can not go belong saturation!!!' ; STOP
   END IF



   WRITE(6,*) ''
   WRITE(6,*) '==========================================================================='
   WRITE(6,*) ' *** density of air at ',TRIM(czt),' => ',  rho_air(t_zt, q_zt, SLP), '[kg/m^3]'

   Cp_ma = cp_air(q_zt)
   WRITE(6,*) ' *** Cp of (moist) air at ',TRIM(czt),' => ', Cp_ma, '[J/K/kg]'
   WRITE(6,*) ''
   rgamma = gamma_moist(t_zt, q_zt)
   WRITE(6,*) ' *** Adiabatic lapse-rate of (moist) air at ',TRIM(czt),' => ', REAL(1000.*rgamma ,4), '[K/1000m]'
   WRITE(6,*) '============================================================================'
   WRITE(6,*) ''
   WRITE(6,*) ''

   ssq = 0.98*q_sat(sst, SLP)

   WRITE(6,*) ''
   WRITE(6,*) ' *** q_',TRIM(czt),'                  =', REAL(1000.*q_zt,4), '[g/kg]'
   WRITE(6,*) ' *** SSQ = 0.98*q_sat(sst) =',            REAL(1000.*ssq ,4), '[g/kg]'
   WRITE(6,*) ''


   !! Must give something more like a potential temperature at zt:
   theta_zt = t_zt + rgamma*zt

   WRITE(6,*) ''
   WRITE(6,*) 'Pot. temp. at ',TRIM(czt),' (using gamma)  =', theta_zt - rt0, ' [deg.C]'



   !! Checking the difference of virtual potential temperature between air at zt and sea surface:
   tmp = virt_temp(theta_zt, q_zt)
   WRITE(6,*) 'Virtual pot. temp. at ',TRIM(czt),'   =', REAL(tmp - rt0 , 4), ' [deg.C]'
   WRITE(6,*) 'Pot. temp. diff. air/sea at ',TRIM(czt),' =', REAL(theta_zt - sst , 4), ' [deg.C]'
   WRITE(6,*) 'Virt. pot. temp. diff. air/sea at ',TRIM(czt),' =', REAL(tmp - virt_temp(sst, ssq), 4), ' [deg.C]'
   WRITE(6,*) ''; WRITE(6,*) ''

   

   WRITE(6,*) 'Give wind speed at zu (m/s):'
   READ(*,*) W10
   WRITE(6,*) ''

   
   !! We have enough to calculate the bulk Richardson number:
   !tmp = Ri_bulk_ecmwf( zt, theta_zt, theta_zt-sst, q_zt, q_zt-ssq, W10 )
   !WRITE(6,*) ' *** Bulk Richardson number "a la ECMWF":', REAL(tmp, 4)
   !tmp = Ri_bulk_ecmwf2( zt, sst, theta_zt, ssq, q_zt, W10 )
   !WRITE(6,*) ' *** Bulk Richardson number "a la ECMWF#2":', REAL(tmp, 4)
   !tmp = Ri_bulk_coare( zt, theta_zt, theta_zt-sst, q_zt, q_zt-ssq, W10 )
   !WRITE(6,*) ' *** Bulk Richardson number "a la COARE":', REAL(tmp, 4)
   tmp = Ri_bulk( zt, sst, theta_zt, ssq, q_zt, W10 )
   WRITE(6,*) ' *** Initial Bulk Richardson number:', REAL(tmp, 4)
   WRITE(6,*) ''
   

   
   IF ( l_use_cswl ) THEN

      WRITE(6,*) ''
      WRITE(6,*) '----------------------------------------------------------'
      WRITE(6,*) '          Will consider the skin temperature!'
      WRITE(6,*) ' => using cool-skin warm-layer param. in COARE and ECMWF'
      WRITE(6,*) ' => need the downwelling radiative fluxes at the surface'
      WRITE(6,*) '----------------------------------------------------------'
      WRITE(6,*) ''
      WRITE(6,*) 'Give downwelling shortwave (solar) radiation at the surface:'
      READ(*,*) rad_sw
      WRITE(6,*)
      WRITE(6,*) 'Give downwelling longwave (infrared) radiation at the surface:'
      READ(*,*) rad_lw
      WRITE(6,*)

   END IF



   DO ialgo = 1, nb_algos

      calgob = TRIM(vca(ialgo))

      zz0 = 0.
      zus = 0. ; zts = 0. ; zqs = 0. ; zL = 0. ; zUN10 = 0.


      !! Mind that TURB_COARE and TURB_ECMWF will modify SST and SSQ if their
      !! respective Cool Skin Warm Layer parameterization is used
      !!  => if optional input arrays "rad_sw, rad_lw, slp" are given a value !

      Ts = sst ; qs = ssq


      SELECT CASE(ialgo)

      CASE(1)

         IF ( l_use_cswl ) THEN

            CALL TURB_COARE3P0( zt, zu, Ts, theta_zt, qs, q_zt, W10, &
               &             Cd, Ch, Ce, theta_zu, q_zu, Ublk,           &
               &             Qsw=(1._wp - oce_alb0)*rad_sw, rad_lw=rad_lw, slp=SLP,      &
               &             xz0=zz0, xu_star=zus, xL=zL, xUN10=zUN10 )
            !! => Ts and qs are updated wrt to skin temperature !

         ELSE
            CALL TURB_COARE3P0( zt, zu, Ts, theta_zt, qs, q_zt, W10, &
               &             Cd, Ch, Ce, theta_zu, q_zu, Ublk,           &
               &             xz0=zz0, xu_star=zus, xL=zL, xUN10=zUN10 )

            !! => Ts and qs are not updated: Ts=sst and qs=ssq


         END IF


      CASE(2)

         IF ( l_use_cswl ) THEN

            CALL TURB_COARE3P6( zt, zu, Ts, theta_zt, qs, q_zt, W10, &
               &             Cd, Ch, Ce, theta_zu, q_zu, Ublk,             &
               &             Qsw=(1._wp - oce_alb0)*rad_sw, rad_lw=rad_lw, slp=SLP,          &
               &             xz0=zz0, xu_star=zus, xL=zL, xUN10=zUN10 )
            !! => Ts and qs are updated wrt to skin temperature !

         ELSE

            CALL TURB_COARE3P6( zt, zu, Ts, theta_zt, qs, q_zt, W10, &
               &             Cd, Ch, Ce, theta_zu, q_zu, Ublk,             &
               &             xz0=zz0, xu_star=zus, xL=zL, xUN10=zUN10 )
            !! => Ts and qs are not updated: Ts=sst and qs=ssq
         END IF


      CASE(3)
         CALL TURB_NCAR( zt, zu, sst, theta_zt, ssq, q_zt, W10, &
            &            Cd, Ch, Ce, theta_zu, q_zu, Ublk,      &
            &            xz0=zz0, xu_star=zus, xL=zL, xUN10=zUN10 )


      CASE(4)

         IF ( l_use_cswl ) THEN

            CALL TURB_ECMWF( zt, zu, Ts, theta_zt, qs, q_zt, W10, &
               &             Cd, Ch, Ce, theta_zu, q_zu, Ublk,        &
               &             Qsw=(1._wp - oce_alb0)*rad_sw, rad_lw=rad_lw, slp=SLP,   &
               &             xz0=zz0, xu_star=zus, xL=zL, xUN10=zUN10)
            !! => Ts and qs are updated wrt to skin temperature !

         ELSE

            CALL TURB_ECMWF( zt, zu, Ts, theta_zt, qs, q_zt, W10, &
               &             Cd, Ch, Ce, theta_zu, q_zu, Ublk,        &
               &             xz0=zz0, xu_star=zus, xL=zL, xUN10=zUN10)
            !! => Ts and qs are not updated: Ts=sst and qs=ssq

         END IF


      CASE DEFAULT
         WRITE(6,*) 'Bulk algorithm #', ialgo, ' is unknown!!!' ; STOP

      END SELECT



      !! Bulk Richardson Number for layer "sea-level -- zu":
      tmp = Ri_bulk(zu, Ts, theta_zu, qs, q_zu, Ublk )
      vBRN(ialgo) = REAL(tmp(1,1),4)


      vTheta_u(ialgo) = REAL(   theta_zu(1,1) -rt0 , 4)   ! Potential temperature at zu

      !! Real temperature at zu
      t_zu = theta_zu ! first guess...
      DO jq = 1, 4
         rgamma = gamma_moist(t_zu, q_zu)
         t_zu = theta_zu - rgamma*zu   ! Real temp.
      END DO

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

      zts = Ch*(theta_zu - Ts)*Ublk/zus
      zqs = Ce*(q_zu     - qs)*Ublk/zus

      vL(ialgo) = zL(1,1)

      vUN10(ialgo) = zUN10(1,1)
      
      !! Turbulent fluxes:
      vTau(ialgo)  = ( rho_zu(1,1) * Cd(1,1) *           W10(1,1)            * Ublk(1,1) )*1000. ! mN/m^2
      tmp = cp_air(q_zu)
      vQH(ialgo)   = rho_zu(1,1)*tmp(1,1)*Ch(1,1) * ( theta_zu(1,1) - Ts(1,1)  ) * Ublk(1,1)
      vEvap(ialgo) = rho_zu(1,1)*Ce(1,1)          * ( qs(1,1)      - q_zu(1,1) ) * Ublk(1,1)  ! mm/s
      tmp = L_vap(Ts)
      vQL(ialgo)   = -1.* ( tmp(1,1)*vEvap(ialgo) )

      vEvap(ialgo) = to_mm_p_day * vEvap(ialgo)  ! mm/day

      vTs(ialgo) = Ts(1,1)
      vqs(ialgo) = qs(1,1)


   END DO


   WRITE(6,*) ''; WRITE(6,*) ''


   IF ( zt < zu ) THEN
      WRITE(6,*) ''; WRITE(6,*) 'Potential temperature and humidity at z = ',TRIM(czt),' :'
      WRITE(6,*) 't_',TRIM(czt),'  =', theta_zt-rt0 ;  WRITE(6,*) 'q_',TRIM(czt),'  =', q_zt
   END IF

   WRITE(6,*) ''; WRITE(6,*) 'Temperatures and humidity at z = ',TRIM(czu),' :'
   WRITE(6,*) '===================================================================================================='
   WRITE(6,*) '  Algorithm:         ',TRIM(vca(1)),'   |   ',TRIM(vca(2)),'    |    ',TRIM(vca(3)),'     |    ',TRIM(vca(4))
   WRITE(6,*) '===================================================================================================='
   WRITE(6,*) '    theta_',TRIM(czu),' =   ', vTheta_u       , '[deg.C]'
   WRITE(6,*) '    t_',TRIM(czu),'     =   ', vT_u      , '[deg.C]'
   WRITE(6,*) '    q_',TRIM(czu),'     =   ', REAL(1000.*vQu ,4)  , '[g/kg]'
   WRITE(6,*) ''
   WRITE(6,*) '      SSQ     =   ', REAL(1000.*qs(1,1), 4)  , '[g/kg]'
   WRITE(6,*) '    Delta t   =   ', REAL(vT_u  - (Ts(1,1)-rt0) , 4)      , '[deg.C]'
   WRITE(6,*) '    Delta q   =   ', REAL(1000.*(vQu - qs(1,1)), 4)  , '[g/kg]'
   WRITE(6,*) ''
   WRITE(6,*) '    Ug (gust) =   ', vUg , '[m/s]'
   WRITE(6,*) ''



   tmp = visc_air(t_zu)
   nu_air = tmp(1,1)


   WRITE(6,*) ''
   WRITE(6,*) 'Kinematic viscosity of air =', nu_air
   WRITE(6,*) ''
   WRITE(6,*) 'With a pressure of', int(SLP)
   WRITE(6,*) 'Density of air at ',TRIM(czu),' =', REAL(vRho_u,4), ' [kg/m^3]'
   IF ( zt < zu )  WRITE(6,*) ' (density of air at ',TRIM(czt),' was ', REAL(rho_air(t_zt, q_zt, SLP),4),')'
   WRITE(6,*) ''
   Cp_ma = cp_air(q_zu)
   WRITE(6,*) ' Cp of moist air at ',TRIM(czu),' => ', REAL(Cp_ma,4), ' [J/K/kg]'
   WRITE(6,*) ''


   WRITE(6,*) ''
   WRITE(6,*) '=============================================================================================='
   WRITE(6,*) '  Algorithm:         ',TRIM(vca(1)),'   |   ',TRIM(vca(2)),'    |    ',TRIM(vca(3)),'     |    ',TRIM(vca(4))
   WRITE(6,*) '=============================================================================================='
   WRITE(6,*) '      C_D     =   ', vCd        , '[10^-3]'
   WRITE(6,*) '      C_E     =   ', vCe        , '[10^-3]'
   WRITE(6,*) '      C_H     =   ', vCh        , '[10^-3]'
   WRITE(6,*) ''
   WRITE(6,*) '      z_0     =   ', vz0        , '[m]'
   WRITE(6,*) '      u*      =   ', vus        , '[m/s]'
   WRITE(6,*) '      L       =   ', vL         , '[m]'
   WRITE(6,*) '      Ri_bulk =   ', vBRN       , '[-]'
   WRITE(6,*) '      UN10    =   ', vUN10      , '[m/s]'
   WRITE(6,*) 'Equ. Charn p. =   ', REAL( grav/(vus*vus)*(vz0 - 0.11*nu_air/vus) , 4)
   WRITE(6,*) ''
   IF ( l_use_cswl ) THEN
      WRITE(6,*) '      Ts      =   ', REAL( vTs-rt0  ,4), '[deg.C]'
      WRITE(6,*) '      qs      =   ', REAL( 1000.*vqs,4), '[g/kg]'
      WRITE(6,*) ''
   END IF
   WRITE(6,*) ' Wind stress  =   ', vTau       , '[mN/m^2]'
   WRITE(6,*) ' Evaporation  =   ', vEvap      , '[mm/day]'
   WRITE(6,*) '    QL        =   ', vQL        , '[W/m^2]'
   WRITE(6,*) '    QH        =   ', vQH        , '[W/m^2]'
   WRITE(6,*) ''
   WRITE(6,*) ''
   CLOSE(6)

END PROGRAM TEST_AEROBULK



SUBROUTINE usage_test()
   !!
   PRINT *,''
   PRINT *,'   List of command line options:'
   PRINT *,''
   PRINT *,' -p   => ask for sea-level pressure, otherwize assume 1010 hPa'
   PRINT *,''
   PRINT *,' -r   => Ask for relative humidity rather than specific humidity'
   PRINT *,''
   PRINT *,' -S   => Use the Cool Skin Warm Layer parameterization to compute'
   PRINT *,'         and use the skin temperature instead of the bulk SST'
   PRINT *,'         only in COARE and ECMWF'
   PRINT *,''
   PRINT *,' -h   => Show this message'
   PRINT *,''
   STOP
   !!
END SUBROUTINE usage_test
!!
