! AeroBulk / 2015 / L. Brodeau


!! Goal of this program is to interactively test various functions found into `mod_phymbl.f90`

PROGRAM TEST_PHYMBL

   USE mod_const
   USE mod_phymbl

   IMPLICIT NONE

   !! Tests to conduct:
   LOGICAL, PARAMETER :: l_test_humi = .false.            ! everything related to the different forms of humidity
   LOGICAL, PARAMETER :: l_test_theta_pressure = .true.  ! everything related to the calculation of potential temperature and pressure at a given height
   
   REAL :: dtime, slp, dT

   INTEGER, PARAMETER :: ntemp = 101

   CHARACTER(len=6) :: cslp
   CHARACTER(len=2) :: cdT

   REAL(wp), DIMENSION(1,ntemp) :: vsst, vmsl, vqsat1, vqsat2, vqsat3, vqsat4, vqsat_ice, &
      &                            vsst_k, vrho, vt10, vq10, vrh

   REAL(wp), DIMENSION(1,ntemp) :: ves_ice, vdes_dt_ice, vdes_dt_ice_fd    ! e_sat(T) and d[e_sat(T)]/dT over ice...
   REAL(wp), DIMENSION(1,ntemp) :: vqs_ice, vdqs_dt_ice, vdqs_dt_ice_fd    ! q_sat(T) and d[q_sat(T)]/dT over ice...

   REAL(wp), PARAMETER :: &
      &   t_min =  -5. , &
      &   t_max =  35.

   REAL(wp) :: zslp, zzt, zta, zqa, zqsat, zP_zt, ztpa
   
   INTEGER :: jt

   CHARACTER(len=128) :: cfout1 !, cfout2



   IF ( l_test_theta_pressure ) THEN

      zslp = Patm  ! sea-level pressure
      

      PRINT *, ''
      PRINT *, 'Give local height above sea-level [m]:'
      READ(*,*) zzt

      PRINT *, ''
      PRINT *, 'Give ABSOLUTE temperature at local height [deg.C]:'
      READ(*,*) zta
      zta = zta + rt0

      zqsat = q_sat( zta, zslp )
      
      PRINT *, ''
      PRINT *, 'Give specific humidity at local height [g/kg]:'
      PRINT *, '    note: saturation value is ',REAL(zqsat*1000.,4), ' g/kg'
      READ(*,*) zqa
      zqa = zqa/1000.
      
      !!PRINT *, 'zslp, zzt, zta, zqa, zqsat = ', zslp, zzt, zta, zqa, zqsat

      !! Spoiler:
      PRINT *, ''
      WRITE(*,'(" SPOILER: `Theta_from_z_P0_T_q_vctr@mod_phymbl` gives: theta =", f7.3, " deg.C at ",f3.0," m")') &
         &   Theta_from_z_P0_T_q( zzt, zslp, zta, zqa )-rt0, zzt
      PRINT *, ''


      !! Calculating pressure at zzt m:
      zP_zt = Pz_from_P0_tz_qz( zzt, zslp, zta, zqa )

      PRINT *, ''
      WRITE(*,'("   *** Reference sea-level pressure we use: ", f7.2, " hPa !")')      zslp/100.
      WRITE(*,'("     ==> computed pressure at ",f4.1,"m is: ", f7.2, " hPa  ")') zzt, zP_zt/100.
      PRINT *, ''

      !! Now that we have the pressure at height zt above sea-level we can convert absolute temperature to potential
      !! using the genuine potential temperature formula:



      !! We are dealing with moist air (spec. hum. of zqa)!
      !zCp = cp_air( zqa ) ; ! specific heat capacity of air at a constant pressure, taking into account humidity      
      !zRm = (1. - zqa)*R_dry + zqa*R_vap  ! universal gas constant for MOIST air     
      !PRINT *, '  *** based on humidity, Cp_air = ', REAL(zCp,4), 'J/K/kg'
      !PRINT *, '  *** based on humidity,  R_air = ', REAL(zRm,4), 'J/K/kg, (R_dry=', REAL(R_dry,4),')'
      !PRINT *, '  *** R/Cp is = ', REAL(zRm/zCp,4)
      !ztpa = zta * ( zslp / zP_zt )**(zRm/zCp)

      ztpa = zta * ( zslp / zP_zt )**rpoiss_dry
      
      PRINT *, '  ==> potential '
      WRITE(*,'("   ==> potential temperature of air at ",f3.0," m is: ", f7.3, " deg.C !")') zzt, ztpa-rt0
      PRINT *, '    ====> from `pot_temp@mod_phymbl` function:', REAL(pot_temp( zta, zP_zt, pPref=zslp ) - rt0,4)
      PRINT *, ''


      PRINT *, ' *** R/Cp used (dry air!):', REAL(rpoiss_dry)
      
      !! For comparison checking what we would have gotten with the gamma lapse-rate version:
      !PRINT *, '  *** rgamma_dry =', rgamma_dry, gamma_moist(zta,0._wp)
      ztpa = zta + rgamma_dry*zzt
      PRINT *, '        ===> using `gamma_dry`   =>', ztpa-rt0, 'deg.C'

      ztpa = zta + gamma_moist(zta, zqa)*zzt
      PRINT *, '        ===> using `gamma_moist` =>', ztpa-rt0, 'deg.C'
      
   END IF
   


   
   IF ( l_test_humi ) THEN

      dtime = (t_max - t_min)/(ntemp - 1)

      PRINT *, ' Temperature increment = ', dtime

      vsst(1,:) = (/ ( t_min + (jt-1)*dtime , jt=1,ntemp ) /)

      !PRINT *, vsst

      vsst_k = vsst + rt0

      !vmsl(:,:) = Patm
      !PRINT *, 'Give pressure at 10m (hPa):'
      !READ(*,*) slp
      !slp = slp*100.
      slp = 101000.   ; ! SLP at 10m !!!
      vmsl(:,:) = slp
      WRITE(cslp,'(i6.6)') INT(slp)


      PRINT *, ''
      PRINT *, 'Give air-sea diff of temperature at 10m (Tair - SST) (K):'
      READ(*,*) dT
      IF ( dT >= 0.) WRITE(cdT,'("p",i1)') INT(ABS(dT))
      IF ( dT <  0.) WRITE(cdT,'("m",i1)') INT(ABS(dT))



      !PRINT *, cslp


      !! Advanced formula:
      vqsat1(:,:) = rdct_qsat_salt*q_sat( vsst_k, vmsl )

      !! Over ice
      vqsat_ice(:,:) =             q_sat( vsst_k, vmsl, l_ice=.TRUE. )

      !! Crude with constant density
      vrho  (:,:) = 1.2_wp
      vqsat2(:,:) = rdct_qsat_salt*q_sat_crude(vsst_k, vrho)

      !! Crude with better density (density at SST! 0m, not 10m !)
      vrho  (:,:) = rho_air(vsst_k, vqsat1, vmsl) ! we use the good q_sat to get a good rho at z=0
      vqsat3(:,:) = rdct_qsat_salt*q_sat_crude(vsst_k, vrho)

      !! Same but using density of air at 10m => T=SST-2 and not saturated => 80% hum!
      vrh(:,:)  = 80._wp                   ! Relative humidity at 10m
      vt10(:,:) = vsst_k + dT              ! Air temperature at 10m

      !! Don't need to do what follows cause we specified that SLP is at 10m !
      !vrho(:,:) = 1.2 ! first guess...
      !DO jt = 1, 5
      !   vq10(:,:) = q_air_rh(vrh, vt10, vmsl-vrho*grav*10.)   ! Specific humidity at 10m
      !   vrho(:,:) = rho_air(vt10, vq10, vmsl-vrho*grav*10.)    ! air density at 10m
      !END DO

      vq10(:,:) = q_air_rh(vrh, vt10, vmsl)   ! Specific humidity at 10m
      vrho(:,:) = rho_air(vt10, vq10, vmsl)    ! air density at 10m

      vqsat4(:,:) = rdct_qsat_salt*q_sat_crude(vsst_k, vrho)


      cfout1 = 'qsat_test_'//cslp//'Pa_dt_'//cdT//'.dat'

      OPEN(11, FILE=cfout1, FORM='formatted', STATUS='unknown',RECL=512)
      WRITE(11,*) '#    SST   Goff   Goff ice   (crude+1.2) Method 3 (crude+rho_0m)   Method 4 (crude+rho_10m)'
      DO jt = 1, ntemp
         WRITE(11,*) REAL(vsst(1,jt),4),   REAL(vqsat1(1,jt),4), REAL(vqsat_ice(1,jt),4), &
            &        REAL(vqsat2(1,jt),4), REAL(vqsat3(1,jt),4), REAL(vqsat4(1,jt),4)
      END DO
      CLOSE(11)

      !! Plot errors in %
      !cfout2 = 'error_qsat_test'//cslp//'.dat'
      !OPEN(11, FILE=cfout2, FORM='formatted', STATUS='unknown',RECL=512)
      !WRITE(11,*) '  Temperature  Method 2 / Method 1 (%)   Method 3 / Method 1 (%)    Method 4 / Method 1 (%) '
      !DO jt = 1, ntemp
      !   WRITE(11,*) REAL(vsst(1,jt),4), REAL(100.*(vqsat2(1,jt)-vqsat1(1,jt))/vqsat1(1,jt),4), &
      !      &      REAL(100.*(vqsat3(1,jt)-vqsat1(1,jt))/vqsat1(1,jt),4), &
      !      &      REAL(100.*(vqsat4(1,jt)-vqsat1(1,jt))/vqsat1(1,jt) , 4)
      !END DO
      !CLOSE(11)


      PRINT *, ''
      !PRINT *, trim(cfout1), ' and ', trim(cfout2), ' written!'
      PRINT *, TRIM(cfout1), ' written!'
      PRINT *, ''


      !! e_sat_ice and derivative over ice
      !! **********************************
      ves_ice(:,:)      =   e_sat_ice( vsst_k )
      vdes_dt_ice(:,:)  =   de_sat_dt_ice( vsst_k )

      !! Finite difference version to control:
      DO jt = 2, ntemp-1
         vdes_dt_ice_fd(1,jt) = (ves_ice(1,jt+1) - ves_ice(1,jt-1)) / ( vsst_k(1,jt+1) - vsst_k(1,jt-1) )
      END DO

      cfout1 = 'e_sat_ice_and_derivative.dat'
      OPEN(11, FILE=cfout1, FORM='formatted', STATUS='unknown',RECL=512)
      WRITE(11,*) '#    SST    e_sat_ice   d[e_sat_ice]/dT  FD:d[e_sat_ice]/dT'
      DO jt = 1, ntemp
         WRITE(11,*) REAL(vsst(1,jt),4),   REAL(ves_ice(1,jt),4), REAL(vdes_dt_ice(1,jt),4), REAL(vdes_dt_ice_fd(1,jt),4)
      END DO
      CLOSE(11)
      PRINT *, ''
      PRINT *, TRIM(cfout1), ' written!'
      PRINT *, ''


      !! e_sat_ice and derivative over ice
      !! **********************************
      vqs_ice(:,:)      =   q_sat(         vsst_k, vmsl, l_ice=.TRUE. )
      vdqs_dt_ice(:,:)  =   dq_sat_dt_ice( vsst_k, vmsl )

      !! Finite difference version to control:
      DO jt = 2, ntemp-1
         vdqs_dt_ice_fd(1,jt) = (vqs_ice(1,jt+1) - vqs_ice(1,jt-1)) / ( vsst_k(1,jt+1) - vsst_k(1,jt-1) )
      END DO

      cfout1 = 'q_sat_ice_and_derivative_'//cslp//'Pa.dat'
      OPEN(11, FILE=cfout1, FORM='formatted', STATUS='unknown',RECL=512)
      WRITE(11,*) '#    SST    q_sat_ice   d[q_sat_ice]/dT  FD:d[q_sat_ice]/dT'
      DO jt = 1, ntemp
         WRITE(11,*) REAL(vsst(1,jt),4),   REAL(vqs_ice(1,jt),4), REAL(vdqs_dt_ice(1,jt),4), REAL(vdqs_dt_ice_fd(1,jt),4)
      END DO
      CLOSE(11)
      PRINT *, ''
      PRINT *, TRIM(cfout1), ' written!'
      PRINT *, ''


   END IF





END PROGRAM TEST_PHYMBL
