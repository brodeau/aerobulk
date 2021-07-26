! AeroBulk / 2020 / L. Brodeau

PROGRAM TEST_AEROBULK_ICE

   USE mod_const
   USE mod_phymbl
   USE mod_blk_ecmwf    ! we need a bulk algo over water for bulk-ice algos that do not deal with leads...
   USE mod_blk_ice_nemo
   USE mod_blk_ice_an05
   USE mod_blk_ice_lg15_io

   IMPLICIT NONE

   INTEGER, PARAMETER :: nb_algos = 3

   CHARACTER(len=8), DIMENSION(nb_algos), PARAMETER :: &
      &      vca = (/ ' NEMdef ' , ' Andreas' , 'L&G 2015' /)

   REAL(wp), DIMENSION(nb_algos) ::  &
      &           vCd, vCe, vCh, vTheta_u, vT_u, vQu, vz0, vus, vRho_u, vUg, vL, vBRN, &
      &           vUN10, vQL, vTau, vQH, vEvap


   !! Over water only:
   !vTau(ialgo), vQH(ialgo), vQL(ialgo),  pEvap=vEvap

   REAL(wp) :: Tau_w, QH_w, QL_w, Evap_w

   INTEGER, PARAMETER ::   &
      &   n_dt   = 21,   &
      &   n_u    = 30,   &
      &   n_q    = 10

   CHARACTER(len=2) :: car

   CHARACTER(len=100) :: calgob
   CHARACTER(len=512) :: ctitle

   INTEGER :: jarg, ialgo, jq, icpt

   INTEGER, PARAMETER :: lx=1, ly=1

   !! For both sea-ice and water:
   REAL(wp), DIMENSION(lx,ly) :: W10, t_zt, theta_zt, q_zt, RH_zt, SLP, d_zt, frci, qsat_zt, Cp_ma, rgamma, tmp

   !! Over sea-ice only :
   REAL(wp), DIMENSION(lx,ly) :: sit, siq, &
      &                          Cd, Ce, Ch,    &
      &                          t_zu, theta_zu, q_zu, rho_zu, &
      &                          Ublk, zz0, zus, zL, zUN10

   !! Over water only :
   REAL(wp), DIMENSION(lx,ly) :: sst, ssq, &
      &                          Cd_w, Ce_w, Ch_w,    &
      &                          theta_zu_w, q_zu_w, &
      &                          Ublk_w, zz0_w, zus_w, zL_w, zUN10_w

   REAL(wp) :: zt, zu, nu_air

   CHARACTER(len=3) :: czt, czu

   LOGICAL :: l_ask_for_slp = .FALSE. , &  !: ask for SLP, otherwize assume SLP = 1010 hPa
      &     l_use_rh      = .FALSE. ,   &  !: ask for RH rather than q for humidity
      &     l_use_dp      = .FALSE. ,   &  !: ask for dew-point temperature rather than q for humidity
      &     l_force_neutral = .FALSE.      !: force air temp and humi at zt to yield perfectly neutral surface atmospheric layer

   nb_iter = 20  ! 20 itterations in bulk algorithm...

   OPEN(6, FORM='formatted', RECL=512)



   CALL usage_test(1)


   jarg = 0

   DO WHILE ( jarg < command_argument_count() )

      jarg = jarg + 1
      CALL get_command_argument(jarg,car)

      SELECT CASE (TRIM(car))

      CASE('-h')
         call usage_test(0)

      CASE('-p')
         l_ask_for_slp = .TRUE.

      CASE('-r')
         l_use_rh = .TRUE.

      CASE('-d')
         l_use_dp = .TRUE.

      CASE('-N')
         l_force_neutral = .TRUE.

      CASE DEFAULT
         WRITE(6,*) 'Unknown option: ', TRIM(car) ; WRITE(6,*) ''
         CALL usage_test(0)

      END SELECT

   END DO



   WRITE(6,*) ''
   WRITE(6,*) '  *** Epsilon aka reps0  = Rd/Rv       (~0.622) =>', reps0
   WRITE(6,*) '  *** Virt. temp. const. = (1-eps)/eps (~0.608) =>', rctv0
   WRITE(6,*) ''





   WRITE(6,*) ''

   WRITE(6,*) 'Give "zu", height of wind speed measurement in meters (generally 10):'
   READ(*,*) zu
   WRITE(6,*) ''


   WRITE(6,*) 'Give "zt", height of air temp. and humidity measurement in meters (generally 2 or 10):'
   READ(*,*) zt
   IF ( (zt > 99.).OR.(zu > 99.) ) THEN
      WRITE(6,*) 'Be reasonable in your choice of zt or zu, they should not exceed a few tenths of meters!' ; STOP
   END IF
   WRITE(6,*) ''
   CALL prtcol( 6, ' zu', zu, 'm' )
   CALL prtcol( 6, ' zt', zt, 'm' )
   WRITE(6,*) ''

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

   WRITE(6,*) 'Give SIT (deg. C):'
   READ(*,*) sit
   sit = sit + rt0
   WRITE(6,*) ''

   WRITE(6,*) 'Give sea-ice concentration > 0 and < 100 (%):'
   READ(*,*) frci
   frci = frci/100._wp ! back to fraction
   IF( (frci(1,1)<0.01_wp).OR.(frci(1,1)>0.99) ) THEN
      WRITE(6,*) ' Come on! It is supposed to be a mix of water and ice !!!' ; STOP
   END IF
   WRITE(6,*) '  => percentage of lead =', REAL(100.*(1.-frci),4), '%'
   WRITE(6,*) ''

   WRITE(6,*) 'Give SST (in the leads) (deg. C):'
   READ(*,*) sst
   sst = sst + rt0
   WRITE(6,*) 'For this sst the latent heat of vaporization is L_vap =', L_vap(sst), ' [J/kg]'
   WRITE(6,*) ''


   IF ( .NOT. l_force_neutral ) THEN
      WRITE(6,*) 'Give temperature at ',TRIM(czt),' (deg. C):'
      READ(*,*) t_zt
      t_zt = t_zt + rt0
      qsat_zt = q_sat(t_zt, SLP)  ! spec. hum. at saturation [kg/kg]
   ELSE
      WRITE(6,*) 'NOT ASKING FOR temperature at ',TRIM(czt),', because you use the "-N" flag => NEUTRAL!'
   END IF
   WRITE(6,*) ''


   !! Asking for humidity:


   IF ( l_use_rh .OR. l_force_neutral ) THEN
      WRITE(6,*) 'Give relative humidity at ',TRIM(czt),' [%]:'
      READ(*,*) RH_zt
      IF ( .NOT. l_force_neutral ) THEN
         q_zt = q_air_rh(RH_zt, t_zt, SLP)
         WRITE(6,*) 'q_',TRIM(czt),' from RH_',TRIM(czt),' =>', 1000*q_zt, ' [g/kg]'
         !WRITE(6,*) 'Inverse => RH from q_zt:', rh_air(q_zt, t_zt, SLP)
      END IF
   ELSEIF ( l_use_dp ) THEN
      WRITE(6,*) 'Give dew-point temperature at ',TRIM(czt),' (deg. C):'
      READ(*,*) d_zt
      d_zt = d_zt + rt0
      q_zt = q_air_dp(d_zt, SLP)
      WRITE(6,*) 'q_',TRIM(czt),' from d_',TRIM(czt),' =>', 1000*q_zt, ' [g/kg]'
      !WRITE(6,*) 'Inverse => RH from q_zt:', rh_air(q_zt, t_zt, SLP)
   ELSE
      WRITE(*, '("Give specific humidity at ",a," (g/kg) (saturation is at ",f6.3," g/kg):")') &
         &      TRIM(czt), 1000.*qsat_zt
      READ(*,*) q_zt
      q_zt = 1.E-3*q_zt
      RH_zt   = rh_air(q_zt, t_zt, SLP)
      WRITE(*,'("  => Relative humidity at ",a," = ",f4.1,"%")') TRIM(czt), RH_zt
      !WRITE(6,*) 'Inverse => q_zt from RH :', 1000*q_air_rh(RH_zt, t_zt, SLP)
   END IF
   WRITE(6,*) ''


   !! Spec. hum. at saturation at temperature == SIT, over sea-ice:
   siq =                q_sat( sit, SLP, l_ice=.TRUE.  )
   ssq = rdct_qsat_salt*q_sat( sit, SLP, l_ice=.FALSE. )


   IF ( l_force_neutral ) THEN
      !! At this stage we know the relative humidity at zt BUT not the temperature
      !! We need to find the air temperature that yields neutral-stability, i.e. vertical gradient of
      !! virtual potential temperature must be 0 !

      !! Find absolute temp. at zt that yields neutral stability (with a humidity of RH_zt)
      t_zt = sit ! first guess
      DO icpt=1, 10
         q_zt = q_air_rh(RH_zt, t_zt, SLP)
         t_zt = virt_temp(sit, siq) / (1._wp + rctv0*q_air_rh(RH_zt, t_zt, SLP)) - gamma_moist(t_zt, q_zt)*zt ! Eq: theta_v_0 = theta_v_zt
      END DO

      qsat_zt = q_sat(t_zt, SLP)  ! spec. hum. at saturation [kg/kg]
      WRITE(6,*) 'We force t_',TRIM(czt),' to =>', REAL(  t_zt-rt0,4), ' [deg.C]'
      WRITE(6,*) 'We force q_',TRIM(czt),' to =>', REAL(1000.*q_zt,4), ' [g/kg] ', ' (sat:', REAL(1000.*qsat_zt,4),')'

   END IF


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



   WRITE(6,*) ''
   WRITE(6,*) ' *** q_',TRIM(czt),'                     =', REAL(1000.*q_zt,4), '[g/kg]'
   WRITE(6,*) ' *** SSQ over ice   =', REAL(1000.*siq ,4), '[g/kg]'
   WRITE(6,*) ' *** SSQ over water =', REAL(1000.*ssq ,4), '[g/kg]'
   WRITE(6,*) ''

   !! Must give something more like a potential temperature at zt:
   theta_zt = t_zt + rgamma*zt

   WRITE(6,*) ''
   WRITE(6,*) 'Pot. temp. at ',TRIM(czt),' (using gamma)  =', theta_zt - rt0, ' [deg.C]'



   !! Checking the difference of virtual potential temperature between air at zt and sea surface:
   tmp = virt_temp(theta_zt, q_zt)
   WRITE(6,*) 'Virtual pot. temp. at ',TRIM(czt),'   =', REAL(tmp - rt0 , 4), ' [deg.C]'
   WRITE(6,*) 'Pot. temp. diff. air/sea at ',TRIM(czt),' =', REAL(theta_zt - sit , 4), ' [deg.C]'
   WRITE(6,*) 'Virt. pot. temp. diff. air/sea at ',TRIM(czt),' =', REAL(tmp - virt_temp(sit, siq), 4), ' [deg.C]'
   WRITE(6,*) ''; WRITE(6,*) ''

   WRITE(6,*) 'Give wind speed at zu (m/s):'
   READ(*,*) W10 ; CALL prtcol( 6, ' wind speed at zu', W10(1,1), 'm/s' )
   WRITE(6,*) ''


   !! We have enough to calculate the bulk Richardson number based on theta and q at z=zt
   tmp = Ri_bulk( zt, sit, theta_zt, siq, q_zt, W10 )
   WRITE(6,*) ' *** Preliminary bulk Richardson (',TRIM(czt),') over ice:',   REAL(tmp, 4)
   tmp = Ri_bulk( zt, sst, theta_zt, ssq, q_zt, W10 )
   WRITE(6,*) ' *** Preliminary bulk Richardson (',TRIM(czt),') over water:', REAL(tmp, 4)
   WRITE(6,*) ''


   !! ###### Computing fluxes over water (leads):

   CALL TURB_ECMWF( 1, zt, zu, sst, theta_zt, ssq, q_zt, W10, .FALSE., .FALSE., &
      &             Cd_w, Ch_w, Ce_w, theta_zu_w, q_zu_w, Ublk_w,               &
      &             xz0=zz0_w, xu_star=zus_w, xL=zL_w, xUN10=zUN10_w )

   WRITE(6,*) ' With ECMWF bulk algorithm we get the following transfer coeff. over water (leads):'
   WRITE(6,*) '  * Cd =', REAL(1000.*Cd_w,4)
   WRITE(6,*) '  * Ch =', REAL(1000.*Ch_w,4)
   WRITE(6,*) '  * Ce =', REAL(1000.*Ce_w,4)


   !! Turbulent fluxes over water

   CALL BULK_FORMULA( zu, sst(1,1), ssq(1,1), theta_zu_w(1,1), q_zu_w(1,1), &
      &              Cd_w(1,1), Ch_w(1,1), Ce_w(1,1), W10(1,1), Ublk_w(1,1), SLP(1,1), &
      &              Tau_w, QH_w, QL_w,  pEvap=Evap_w )

   WRITE(6,*) ' With ECMWF bulk algorithm we get the following fluxes over water (leads):'
   WRITE(6,*) '  * Tau =', REAL(1000.*Tau_w,4), '[mN/m^2]'
   WRITE(6,*) '  * QH  =', REAL(QH_w,4),          '[W/m^2]'
   WRITE(6,*) '  * QL  =', REAL(QL_w,4),          '[W/m^2]'


   !! ############################################################


   !STOP'LOLO'

   DO ialgo = 1, nb_algos

      calgob = TRIM(vca(ialgo))

      zz0 = 0.
      zus = 0. ; zL = 0. ; zUN10 = 0.


      SELECT CASE(ialgo)

      CASE(1)
         CALL TURB_ICE_NEMO( zt, zu, sit, theta_zt, siq, q_zt, W10,   &
            &                Cd, Ch, Ce, theta_zu, q_zu, Ublk,         &
            &                xz0=zz0, xu_star=zus, xL=zL, xUN10=zUN10 )

      CASE(2)
         CALL TURB_ICE_AN05( zt, zu, sit, theta_zt, siq, q_zt, W10,   &
            &                Cd, Ch, Ce, theta_zu, q_zu, Ublk,         &
            &                xz0=zz0, xu_star=zus, xL=zL, xUN10=zUN10 )

      CASE(3)
         CALL TURB_ICE_LG15_IO( zt, zu, sit, theta_zt, siq, q_zt, W10, frci,  &
            &                   Cd, Ch, Ce, theta_zu, q_zu, Ublk,             &
            &                   xz0=zz0, xu_star=zus, xL=zL, xUN10=zUN10 )

      CASE DEFAULT
         WRITE(6,*) 'Sea-ice bulk algorithm #', ialgo, ' is unknown!!!' ; STOP

      END SELECT



      !! Bulk Richardson Number for layer "sea-level -- zu":
      tmp = Ri_bulk(zu, sit, theta_zu, siq, q_zu, Ublk )
      vBRN(ialgo) = tmp(1,1)


      vTheta_u(ialgo) = theta_zu(1,1) -rt0   ! Potential temperature at zu

      !! Absolute temperature at zu
      t_zu = theta_zu ! first guess...
      DO jq = 1, 4
         rgamma = gamma_moist(0.5*(t_zu+sit), q_zu)
         t_zu = theta_zu - rgamma*zu   ! Absolute temp.
      END DO

      vCd(ialgo) = 1000.*Cd(1,1)
      vCh(ialgo) = 1000.*Ch(1,1)
      vCe(ialgo) = 1000.*Ce(1,1)

      vT_u(ialgo) =  t_zu(1,1) -rt0     ! Absolute temp.
      vQu(ialgo)  =   q_zu(1,1)


      !! Gustiness contribution:
      vUg(ialgo) = Ublk(1,1) - W10(1,1)

      !! z0 et u*:
      vz0(ialgo) = zz0(1,1)
      vus(ialgo) = zus(1,1)

      !zts = Ch*(theta_zu - sit)*Ublk/zus
      !zqs = Ce*(q_zu     - siq)*Ublk/zus

      vL(ialgo) = zL(1,1)

      vUN10(ialgo) = zUN10(1,1)

      !! Air density at zu (10m)
      rho_zu = rho_air(t_zu, q_zu, SLP)
      tmp = SLP - rho_zu*grav*zu
      rho_zu = rho_air(t_zu, q_zu, tmp)
      vRho_u(ialgo) = rho_zu(1,1)



      !! Turbulent fluxes:

      CALL BULK_FORMULA( zu, sit(1,1), siq(1,1), theta_zu(1,1), q_zu(1,1), Cd(1,1), Ch(1,1), Ce(1,1), W10(1,1), Ublk(1,1), SLP(1,1), &
         &              vTau(ialgo), vQH(ialgo), vQL(ialgo),  pEvap=vEvap(ialgo), l_ice=.TRUE. )

      vTau(ialgo)  =       1000. *  vTau(ialgo)  ! mN/m^2
      vEvap(ialgo) = to_mm_p_day * vEvap(ialgo)  ! mm/day


   END DO


   WRITE(6,*) ''; WRITE(6,*) ''

   WRITE(ctitle,*) '  Algorithm:      ',TRIM(vca(1)),'   |   ',TRIM(vca(2)),'    |  ',TRIM(vca(3))!,'   |   ',TRIM(vca(4))

   IF ( zt < zu ) THEN
      WRITE(6,*) ''; WRITE(6,*) 'Potential temperature and humidity at z = ',TRIM(czt),' :'
      WRITE(6,*) 't_',TRIM(czt),'  =', theta_zt-rt0 ;  WRITE(6,*) 'q_',TRIM(czt),'  =', q_zt
   END IF

   WRITE(6,*) ''; WRITE(6,*) 'Temperatures and humidity at z = ',TRIM(czu),' :'
   WRITE(6,*) '===================================================================================================='
   WRITE(6,*) TRIM(ctitle)
   WRITE(6,*) '===================================================================================================='
   WRITE(6,*) '    theta_',TRIM(czu),' =   ', REAL(vTheta_u,  4)       , '[deg.C]'
   WRITE(6,*) '    t_',TRIM(czu),'     =   ', REAL(vT_u    ,  4)       , '[deg.C]'
   WRITE(6,*) '    q_',TRIM(czu),'     =   ', REAL(1000.*vQu, 4)       , '[g/kg]'
   WRITE(6,*) ''
   WRITE(6,*) '      SIQ     =   ', REAL(1000.*siq(1,1), 4)             , '[g/kg]'
   WRITE(6,*) '    Delta t   =   ', REAL(vT_u  - (sit(1,1)-rt0) , 4)    , '[deg.C]'
   WRITE(6,*) '    Delta q   =   ', REAL(1000.*(vQu - siq(1,1)), 4)     , '[g/kg]'
   WRITE(6,*) ''
   WRITE(6,*) '    Ug (gust) =   ', REAL(vUg, 4)                       , '[m/s]'
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
   WRITE(6,*) TRIM(ctitle)
   WRITE(6,*) '=============================================================================================='
   WRITE(6,*) '      C_D     =   ', REAL(vCd  ,4) , '[10^-3]'
   WRITE(6,*) '      C_E     =   ', REAL(vCe  ,4) , '[10^-3]'
   WRITE(6,*) '      C_H     =   ', REAL(vCh  ,4) , '[10^-3]'
   WRITE(6,*) ''
   WRITE(6,*) '      z_0     =   ', REAL(vz0  ,4) , '[m]'
   WRITE(6,*) '      u*      =   ', REAL(vus  ,4) , '[m/s]'
   WRITE(6,*) '      L       =   ', REAL(vL   ,4) , '[m]'
   WRITE(6,*) '      Ri_bulk =   ', REAL(vBRN ,4) , '[-]'
   WRITE(6,*) '      UN10    =   ', REAL(vUN10,4) , '[m/s]'
   WRITE(6,*) 'Equ. Charn p. =   ', REAL( grav/(vus*vus)*(vz0 - 0.11*nu_air/vus) , 4)
   WRITE(6,*) ''
   WRITE(6,*) ' Wind stress  =   ', REAL(vTau ,4) , '[mN/m^2]'
   WRITE(6,*) ' Evaporation  =   ', REAL(vEvap,4) , '[mm/day]'
   WRITE(6,*) '    QL        =   ', REAL(vQL  ,4) , '[W/m^2]'
   WRITE(6,*) '    QH        =   ', REAL(vQH  ,4) , '[W/m^2]'
   WRITE(6,*) ''

   WRITE(6,*) ''
   CLOSE(6)


CONTAINS

   SUBROUTINE prtcol( id, cs1, pval, cs2 )
      INTEGER,          INTENT(in) :: id
      CHARACTER(len=*), INTENT(in) :: cs1
      REAL(wp),         INTENT(in) :: pval
      CHARACTER(len=*), INTENT(in) :: cs2
      WRITE(id,*) ACHAR(27)//'[95m ***'//TRIM(cs1)//' =',REAL(pval,4),TRIM(cs2)//ACHAR(27)//'[0m'
   END SUBROUTINE prtcol


   SUBROUTINE usage_test( icontinue )
      INTEGER, INTENT(in) :: icontinue
      !!
      PRINT *,''
      IF (icontinue>=1) THEN
         PRINT *,''
         PRINT *,'                 ======   A e r o B u l k   ====='
         PRINT *,''
         PRINT *,'                      (L. Brodeau, 2015-2019)'
         PRINT *,''
         PRINT *,'                   Simple interactive test-case'
      END IF
      PRINT *,'##########################################################################'
      PRINT *,'   List of command line options:'
      PRINT *,''
      PRINT *,'   -p   => Ask for sea-level pressure, otherwize assume SLP = 1010 hPa'
      PRINT *,''
      PRINT *,'   -r   => Ask for relative humidity rather than specific humidity'
      PRINT *,''
      PRINT *,'   -N   => Force neutral stability in surface atmospheric layer'
      PRINT *,'           -> will not ask for air temp. at zt, instead will only'
      PRINT *,'           ask for relative humidity at zt and will force the air'
      PRINT *,'           temperature at zt that yields a neutral-stability!'
      PRINT *,''
      PRINT *,'   -h   => Show this message'
      PRINT *,'##########################################################################'
      PRINT *,''
      IF (icontinue==0) STOP
      !!
   END SUBROUTINE usage_test

END PROGRAM TEST_AEROBULK_ICE
