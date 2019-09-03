! AeroBulk / 2015 / L. Brodeau

PROGRAM TEST_SKIN_CORR

   USE mod_const
   USE mod_phymbl

   USE mod_blk_coare3p0
   USE mod_blk_ecmwf
   !USE mod_blk_coare3p6
   USE mod_wl_coare3p6

   IMPLICIT NONE

   CHARACTER(len=8) :: calgo = 'coare3p0'

   REAL(wp), PARAMETER ::   &
      &   zt  =  2. ,  &
      &   zu  = 10. ,  &
      &   wind_max = 50.,  &
      & to_mm_p_day = 24.*3600.  !: freshwater flux: from kg/s/m^2 == mm/s to mm/day

   INTEGER, PARAMETER :: Nw10 = 1201

   REAL, DIMENSION(Nw10) :: &
      &   t_w10, t_dT

   CHARACTER(len=2) :: car
   CHARACTER(len=256) :: cf_out

   INTEGER :: jarg, ians, jw

   INTEGER, PARAMETER :: lx=1, ly=1
   REAL(wp),    DIMENSION(lx,ly) :: Ublk, sst_s, ssq_s
   REAL(wp),    DIMENSION(lx,ly) :: tmp1, tmp2, tmp3 !LOLO remove!!!

   REAL(wp), DIMENSION(lx,ly) :: sst, qsat_zt, rlon, rad_sw, rad_lw, SLP, &
      &  W10, t_zt, theta_zt, q_zt, RH_zt, theta_zu, q_zu, ssq, tmp, dtheta_v

   REAL(wp) :: rt_day_h_utc ! Current UTC time in number of hours since midnight
   INTEGER  :: it_day_s_utc ! Current UTC time in number of hours since midnight

   REAL(wp), DIMENSION(lx,ly) :: Cd, Ce, Ch, Cp_ma, rgamma

   REAL(wp) :: dw, dw0

   CHARACTER(len=3) :: czt, czu

   LOGICAL :: l_ask_for_slp = .FALSE. , &  !: ask for SLP, otherwize assume SLP = 1010 hPa
      &       l_use_rh      !: ask for RH rather than q for humidity

   l_use_rh = .TRUE.

   jpi = lx ; jpj = ly


   nb_itt = 5  ! nb. itterations in bulk algorithm...



   !! Table of wind speeds from 0 to wind_max  m/s :
   dw0 = wind_max/(Nw10 - 1)
   t_w10(1) = 0.
   DO jw = 2, Nw10
      dw = 0.25*dw0
      IF ( t_w10(jw-1) >= 5.  ) dw =    dw0
      IF ( t_w10(jw-1) >= 30. ) dw = 2.*dw0
      t_w10(jw) = t_w10(jw-1) + dw
   END DO


   OPEN(6, FORM='formatted', RECL=512)


   jarg = 0

   DO WHILE ( jarg < command_argument_count() )

      jarg = jarg + 1
      CALL get_command_argument(jarg,car)

      SELECT CASE (trim(car))

      CASE('-h')
         call usage_test()

      CASE('-p')
         !WRITE(6,*) 'Will ask for SLP'
         l_ask_for_slp = .TRUE.

         !CASE('-r')
         !   !PRINT *, "Using RH instead of q!!!"
         !   l_use_rh = .TRUE.

      CASE DEFAULT
         WRITE(6,*) 'Unknown option: ', trim(car) ; WRITE(6,*) ''
         CALL usage_test()

      END SELECT

   END DO



   WRITE(6,*) ''
   WRITE(6,*) '  *** Epsilon            = Rd/Rv       (~0.622) =>', reps0
   WRITE(6,*) '  *** Virt. temp. const. = (1-eps)/eps (~0.608) =>', rctv0
   WRITE(6,*) ''





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

   ians=0
   DO WHILE ( (ians<1).OR.(ians>3) )
      WRITE(6,*) 'Which algo to use? "coare3p0" => 1 , "ecmwf" => 2 , "coare3p6" => 3 :'
      READ(*,*) ians
      IF ( ians == 1 ) calgo = 'coare3p0'
      IF ( ians == 2 ) calgo = 'ecmwf   '
      IF ( ians == 3 ) calgo = 'coare3p6'
   END DO
   WRITE(6,*) '  ==> your choice: ', TRIM(calgo)
   WRITE(6,*) ''
   

   IF ( l_ask_for_slp ) THEN
      WRITE(6,*) 'Give sea-level pressure (hPa):'
      READ(*,*) SLP
      SLP = SLP*100.
   ELSE
      SLP = Patm
      WRITE(6,*) 'Using a sea-level pressure of ', Patm
   END IF
   WRITE(6,*) ''

   WRITE(6,*) 'Give the longitude (deg.East) [negative values accepted!]:'
   READ(*,*) rlon
   WRITE(6,*) ''

   WRITE(6,*) 'Give the current UTC time in number of hours since midnight:'
   READ(*,*) rt_day_h_utc
   it_day_s_utc = INT( rt_day_h_utc*3600._wp )
   WRITE(6,*) '  ==> UTC time in number of SECONDS since midnight:', it_day_s_utc
   WRITE(6,*) ''

   
   WRITE(6,*) 'Give RAD_SW (W/m^2):'
   READ(*,*) rad_sw
   WRITE(6,*) ''

   WRITE(6,*) 'Give RAD_LW (W/m^2):'
   READ(*,*) rad_lw
   WRITE(6,*) ''


   WRITE(6,*) 'Give SST (deg. C):'
   READ(*,*) sst
   sst = sst + rt0
   WRITE(6,*) 'For this sst the latent heat of vaporization is L_vap =', L_vap(sst), ' [J/kg]'
   WRITE(6,*) ''

   WRITE(6,*) 'Give temperature at ',trim(czt),' (deg. C):'
   READ(*,*) t_zt
   t_zt = t_zt + rt0
   WRITE(6,*) ''


   !! Asking for humidity:
   qsat_zt = q_sat(t_zt, SLP)  ! spec. hum. at saturation [kg/kg]

   !IF ( l_use_rh ) THEN
   !WRITE(6,*) 'Give relative humidity at ',TRIM(czt),' (%):'
   !READ(*,*) RH_zt
   RH_zt = 80.
   RH_zt = 1.E-2*RH_zt
   q_zt = q_air_rh(RH_zt, t_zt, SLP)
   WRITE(6,*) 'q_',TRIM(czt),' from RH_',TRIM(czt),' =>', 1000*q_zt, ' [g/kg]'
   !WRITE(6,*) 'Inverse => RH from q_zt:', 100*rh_air(q_zt, t_zt, SLP)
   !ELSE
   !   WRITE(*, '("Give specific humidity at ",a," (g/kg) (saturation is at ",f6.3," g/kg):")') &
   !      &      TRIM(czt), 1000.*qsat_zt
   !   READ(*,*) q_zt
   !   q_zt = 1.E-3*q_zt
   !   RH_zt   = rh_air(q_zt, t_zt, SLP)
   !   WRITE(*,'("  => Relative humidity at ",a," = ",f4.1,"%")') TRIM(czt), 100*RH_zt
   !   !WRITE(6,*) 'Inverse => q_zt from RH :', 1000*q_air_rh(RH_zt, t_zt, SLP)
   !END IF
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

   ssq = rdct_qsat_salt*q_sat(sst, SLP)

   WRITE(6,*) '' ;  WRITE(6,*) 'SSQ = 0.98*q_sat(sst) =', 1000.*ssq, '[g/kg]' ; WRITE(6,*) ''


   !! Must give something more like a potential temperature at zt:
   theta_zt = t_zt + rgamma*zt

   WRITE(6,*) ''
   WRITE(6,*) 'Pot. temp. at ',TRIM(czt),' (using gamma)  =', theta_zt - rt0, ' [deg.C]'



   !! Checking the difference of virtual potential temperature between air at zt and sea surface:
   tmp = theta_zt*(1. + rctv0*q_zt)

   dtheta_v = tmp - sst*(1. + rctv0*ssq)

   WRITE(6,*) 'Virtual pot. temp. at ',TRIM(czt),'   =', REAL(tmp - rt0 , 4), ' [deg.C]'
   WRITE(6,*) 'Pot. temp. diff. air/sea at ',TRIM(czt),' =', REAL(theta_zt - sst , 4), ' [deg.C]'
   WRITE(6,*) 'Virt. pot. temp. diff. air/sea at ',TRIM(czt),' =', REAL(dtheta_v , 4), ' [deg.C]'
   WRITE(6,*) ''; WRITE(6,*) ''


   !WRITE(6,*) 'Give wind speed at zu=10m (m/s):'
   !READ(*,*) W10
   !WRITE(6,*) ''



   DO jw = 1, Nw10

      W10 = t_w10(jw)

      sst_s = sst

      SELECT CASE (TRIM(calgo))

      CASE('coare3p0')
         CALL turb_coare3p0( zt, zu, sst_s, theta_zt, ssq_s, q_zt, W10, &
            &                Cd, Ch, Ce, theta_zu, q_zu, Ublk,          &
            &                Qsw=(1._wp - oce_alb0)*rad_sw, rad_lw=rad_lw, slp=slp )

      CASE('ecmwf')
         CALL turb_ecmwf( zt, zu, sst_s, theta_zt, ssq_s, q_zt, W10, &
            &             Cd, Ch, Ce, theta_zu, q_zu, Ublk,          &
            &             Qsw=(1._wp - oce_alb0)*rad_sw, rad_lw=rad_lw, slp=slp )



      CASE('coare3p6')
         !CALL turb_coare3p0( zt, zu, sst_s, theta_zt, ssq_s, q_zt, W10, &
         !   &                Cd, Ch, Ce, theta_zu, q_zu, Ublk,          &
         !   &                Qsw=(1._wp - oce_alb0)*rad_sw, rad_lw=rad_lw, slp=slp )
         PRINT *, 'Booh!!! Not ready yet!'

         tmp1 = 0.01
         tmp2 = 0.
         tmp3 = 0.
         CALL WL_COARE3P6_2( rad_sw, rad_sw*0.-500., rad_sw*0.+0.001, sst_s, tmp1, tmp2, tmp3, rlon, it_day_s_utc, 60. )
         STOP

         
         
      CASE DEFAULT
         PRINT *, 'Unknown algorithm: ', TRIM(calgo) ; PRINT *, ''
         CALL usage_test()

      END SELECT

      !! Surface temperature correction as a function of wind:
      t_dT(jw) = sst_s(1,1) - sst(1,1)




   END DO



   PRINT *, ''
   PRINT *, ''
   PRINT *, ''

   IF ( dtheta_v(1,1) < 0. ) THEN
      WRITE(cf_out,'("dat/dT_skin_vs_wind_",a,"_SST",i2.2,"_SW",i4.4,"_LW",i4.4,"_RH",i2.2,"_unstable.dat")') &
         &  trim(calgo), INT(SST-rt0), INT(RAD_SW), INT(RAD_LW), INT(100.*RH_zt)
   ELSE
      WRITE(cf_out,'("dat/dT_skin_vs_wind_",a,"_SST",i2.2,"_SW",i4.4,"_LW",i4.4,"_RH",i2.2,"_stable.dat")') &
         &  trim(calgo), INT(SST-rt0), INT(RAD_SW), INT(RAD_LW), INT(100.*RH_zt)
   END IF

   PRINT *, trim(cf_out)

   OPEN(unit = 11, file = TRIM(cf_out), status = 'unknown')
   WRITE(11,'("#     Wind (U0)          Ts - SST")')
   DO jw = 1, Nw10
      WRITE(11,'(1f16.8, " " ,1f16.8, " " ,1f16.8)') t_w10(jw), t_dT(jw)
   END DO
   CLOSE(11)

   PRINT *, ''
   PRINT *, '  *** Just created file ', TRIM(cf_out)

END PROGRAM TEST_SKIN_CORR




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
