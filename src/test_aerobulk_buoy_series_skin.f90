! AeroBulk / 2019 / L. Brodeau

!! https://www.pmel.noaa.gov/ocs/flux-documentation


PROGRAM TEST_AEROBULK_BUOY_SERIES_SKIN

   USE mod_const
   USE mod_phymbl

   !USE mod_blk_coare3p0
   USE mod_blk_coare3p6
   USE mod_wl_coare3p6
   !USE mod_blk_ncar
   !USE mod_blk_ecmwf

   IMPLICIT NONE

   LOGICAL, PARAMETER :: ldebug=.true.

   INTEGER, PARAMETER :: nb_measurements = 115

   INTEGER, PARAMETER :: nb_algos = 1
   
   INTEGER, PARAMETER :: nb_itt_wl = 2 !!LOLO
   
   CHARACTER(len=800) :: cf_data='0', cblabla, cf_out='output.dat'

   CHARACTER(len=8), DIMENSION(nb_algos), PARAMETER :: &
                                !&      vca = (/ 'coare3p0', 'coare3p6', 'ncar    ', 'ecmwf   ' /)
      &      vca = (/ 'coare3p6' /)

   REAL(wp), PARAMETER ::   &
      & to_mm_p_day = 24.*3600.  !: freshwater flux: from kg/s/m^2 == mm/s to mm/day

   INTEGER, PARAMETER ::   &
      &   n_dt   = 21,   &
      &   n_u    = 30,   &
      &   n_q    = 10

   CHARACTER(len=2) :: car

   CHARACTER(len=100) :: &
      &   calgob

   INTEGER :: jt, jarg, jl, ialgo, jq, jtt

   INTEGER, PARAMETER :: nx=1, ny=1, Nt=nb_measurements

   REAL(wp),    DIMENSION(nx,ny,Nt) :: Ublk, zz0, zus, zts, zqs, zL, zUN10

   CHARACTER(len=19)                :: cdt
   CHARACTER(len=19), DIMENSION(Nt) :: ctime
   CHARACTER(len=8),  DIMENSION(Nt) :: cdate
   CHARACTER(len=4),  DIMENSION(Nt) :: clock
   CHARACTER(len=2),  DIMENSION(Nt) :: chh, cmn ! hours and minutes
   CHARACTER(len=16), DIMENSION(Nt) :: cldate ! human!
   INTEGER(4) :: iclock, ihh, imm
   INTEGER(4)      , DIMENSION(Nt)  :: isecday


   INTEGER(8), DIMENSION(Nt) :: idate


   !! Input (or deduced from input) variables:
   REAL(wp), DIMENSION(nx,ny,Nt) :: sst, qsat_zt, SLP, W10, t_zt, theta_zt, q_zt, RH_zt, d_zt, &
      &                          rad_sw, rad_lw, precip, rlat, rlon


   REAL(wp), DIMENSION(nx,ny,Nt) :: Ts, t_zu, theta_zu, q_zu, qs, rho_zu, dummy, &
      &                             dT, pTau_ac, pQ_ac

   REAL(wp), DIMENSION(nx,ny) :: ssq, rgamma, Cp_ma, tmp
   
   
   REAL(wp), DIMENSION(nx,ny,Nt) :: Cd, Ce, Ch, QH, QL, EVAP, RiB



   REAL(4), DIMENSION(nb_algos,Nt) ::  &
      &           vCd, vCe, vCh, vTheta_u, vT_u, vQu, vz0, vus, vRho_u, vUg, vL, vBRN, &
      &           vUN10, vQL, vTau, vQH, vEvap, vTs, vqs



   REAL(wp) :: zt, zu

   CHARACTER(len=3) :: czt, czu

   LOGICAL :: l_ask_for_slp = .FALSE. , &  !: ask for SLP, otherwize assume SLP = 1010 hPa
      &     l_use_rh      = .FALSE. ,   &  !: ask for RH rather than q for humidity
      &     l_use_dp      = .FALSE. ,   &  !: ask for dew-point temperature rather than q for humidity
      &     l_use_cswl    = .FALSE.        !: compute and use the skin temperature
   !                                       !: (Cool Skin Warm Layer parameterization)
   !                                       !:  => only in COARE and ECMWF

   jpi = nx ; jpj = ny

   nb_itt = 20  ! 20 itterations in bulk algorithm...

   OPEN(6, FORM='formatted', RECL=512)


   !----------------------------------------------

   jarg = 0

   DO WHILE ( jarg < command_argument_COUNT() )

      jarg = jarg + 1
      CALL get_command_argument(jarg,car)

      SELECT CASE (TRIM(car))

      CASE('-h')
         call usage_test()

      CASE('-f')
         jarg = jarg + 1
         CALL get_command_ARGUMENT(jarg,cf_data)

      CASE('-p')
         l_ask_for_slp = .TRUE.

      CASE('-r')
         l_use_rh = .TRUE.

      CASE('-d')
         l_use_dp = .TRUE.

      CASE('-S')
         l_use_cswl = .TRUE.

      CASE DEFAULT
         WRITE(6,*) 'Unknown option: ', trim(car) ; WRITE(6,*) ''
         CALL usage_test()

      END SELECT

   END DO


   IF ( trim(cf_data) == '0' ) CALL usage_test()


   WRITE(6,*) ''
   WRITE(6,*) ' *** Input file is ',TRIM(cf_data)
   WRITE(6,*) ''
   WRITE(6,*) '  *** Epsilon            = Rd/Rv       (~0.622) =>', reps0
   WRITE(6,*) '  *** Virt. temp. const. = (1-eps)/eps (~0.608) =>', rctv0
   WRITE(6,*) ''





   WRITE(6,*) ''


   OPEN(11, FILE=TRIM(cf_data), FORM='formatted', STATUS='old')
   READ(11,*) cblabla ! first commented line...
   DO jt = 1, Nt
      READ(11,*) ctime(jt), W10(1,1,jt), sst(1,1,jt), t_zt(1,1,jt), q_zt(1,1,jt), rad_sw(1,1,jt), rad_lw(1,1,jt), precip(1,1,jt), rlat(1,1,jt), rlon(1,1,jt), dummy(1,1,jt)
   END DO
   CLOSE(11)


   DO jt = 1, Nt
      cdt = ctime(jt)
      cdate(jt) = cdt(1:8)
      clock(jt) = cdt(9:13)
      chh(jt)   = cdt(9:10)
      cmn(jt)   = cdt(11:12)
      READ(clock(jt),'(i4.4)') iclock
      READ(chh(jt),'(i2.2)') ihh
      READ(cmn(jt),'(i2.2)') imm
      isecday(jt) = ihh*3600. + imm*60.
      WRITE(cldate(jt),'(a4,"/",a2,"/",a2,"-",a2,":",a2)') cdt(1:4), cdt(5:6), cdt(7:8), chh(jt), cmn(jt)
      IF (ldebug) PRINT *, ' *** date = ', cldate(jt)      
   END DO
   IF (ldebug) PRINT *, ''




   !   it_n(1,:)  = isecday(:)
   !   it_b(1,1)  = 0
   !   it_b(1,2:ny) = isecday(:ny-1)
   !DO jt = 1, Nt
   !   PRINT *, ' it_b, it_n =', it_b(1,jt), it_n(1,jt)
   !END DO
   !STOP


   !! zu and zt
   !! ~~~~~~~~~
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

   !! SLP
   !! ~~~
   IF ( l_ask_for_slp ) THEN
      WRITE(6,*) 'Give sea-level pressure (hPa):'
      READ(*,*) SLP
      SLP = SLP*100.
   ELSE
      SLP = Patm
      WRITE(6,*) 'Using a sea-level pressure of ', Patm
   END IF
   WRITE(6,*) ''

   !! Back to SI unit...
   sst  = sst + rt0
   t_zt = t_zt + rt0
   q_zt = 1.E-3*q_zt

   IF (ldebug) THEN
      !                   19921125132100   4.700000       302.1500       300.8500      1.7600000E-02  0.0000000E+00   428.0000
      WRITE(6,*) '*       idate     ,   wind    ,       SST    ,     t_zt     ,      q_zt      ,    rad_sw     , rad_lw  :'
      DO jt = 1, Nt
         WRITE(6,*) cldate(jt), REAL(W10(:,:,jt),4), REAL(sst(:,:,jt),4), REAL(t_zt(:,:,jt),4), REAL(q_zt(:,:,jt),4), REAL(rad_sw(:,:,jt),4), REAL(rad_lw(:,:,jt),4)
      END DO
   END IF



   !! Some initializations:

   ialgo = 1
   
   zz0 = 0.
   zus = 0. ; zts = 0. ; zqs = 0. ; zL = 0. ; zUN10 = 0.

   pTau_ac = 0.
   pQ_ac   = 0.



   !! Time loop:
   DO jt = 1, Nt

      WRITE(6,*) ''; WRITE(6,*) ''
      WRITE(6,*) '##############################################'
      WRITE(6,*) '#### Time = ', cldate(jt)
      WRITE(6,*) '##############################################'
      WRITE(6,*) ''
      WRITE(6,*) '           ---- BEFORE BULK ALGO + CSWL ----'
      WRITE(6,*) ''
      WRITE(6,*) ' *** density of air at ',TRIM(czt),' => ',  rho_air(t_zt(:,:,jt), q_zt(:,:,jt), SLP(:,:,jt)), '[kg/m^3]'
      WRITE(6,*) ''
      
      Cp_ma(:,:) = cp_air(q_zt(:,:,jt))
      WRITE(6,*) ' *** Cp of (moist) air at ',TRIM(czt),' => ', REAL(Cp_ma,4), '[J/K/kg]'
      WRITE(6,*) ''
      
      rgamma(:,:) = gamma_moist(t_zt(:,:,jt), q_zt(:,:,jt))
      WRITE(6,*) ' *** Adiabatic lapse-rate of (moist) air at ',TRIM(czt),' => ', REAL(1000.*rgamma ,4), '[K/1000m]'
      WRITE(6,*) ''
      
      ssq = 0.98*q_sat(sst(:,:,jt), SLP(:,:,jt))
      WRITE(6,*) ' *** SSQ = 0.98*q_sat(sst) =',            REAL(1000.*ssq ,4), '[g/kg]'
      
      !! Must give something more like a potential temperature at zt:
      theta_zt(:,:,jt) = t_zt(:,:,jt) + rgamma(:,:)*zt
      WRITE(6,*) ''
      WRITE(6,*) ' *** Pot. temp. at ',TRIM(czt),' (using gamma)  =', theta_zt(:,:,jt) - rt0, ' [deg.C]'
      WRITE(6,*) ''


      !! Checking the difference of virtual potential temperature between air at zt and sea surface:
      tmp = virt_temp(theta_zt(:,:,jt), q_zt(:,:,jt))
      WRITE(6,*) ' *** Virtual pot. temp. at ',TRIM(czt),'   =', REAL(tmp - rt0 , 4), ' [deg.C]'
      WRITE(6,*) ''
      WRITE(6,*) ' *** Pot. temp. diff. air/sea at ',TRIM(czt),' =', REAL(theta_zt(:,:,jt) - sst(:,:,jt) , 4), ' [deg.C]'
      WRITE(6,*) ''
      WRITE(6,*) ' *** Virt. pot. temp. diff. air/sea at ',TRIM(czt),' =', REAL(tmp - virt_temp(sst(:,:,jt), ssq), 4), ' [deg.C]'
      WRITE(6,*) ''



      !WRITE(6,*) 'Give wind speed at zu (m/s):'
      !READ(*,*) W10
      !WRITE(6,*) ''


      !! We have enough to calculate the bulk Richardson number:
      !tmp = Ri_bulk_ecmwf( zt, theta_zt, theta_zt-sst, q_zt, q_zt-ssq, W10 )
      !WRITE(6,*) ' *** Bulk Richardson number "a la ECMWF":', REAL(tmp, 4)
      !tmp = Ri_bulk_ecmwf2( zt, sst, theta_zt, ssq, q_zt, W10 )
      !WRITE(6,*) ' *** Bulk Richardson number "a la ECMWF#2":', REAL(tmp, 4)
      !tmp = Ri_bulk_coare( zt, theta_zt, theta_zt-sst, q_zt, q_zt-ssq, W10 )
      !WRITE(6,*) ' *** Bulk Richardson number "a la COARE":', REAL(tmp, 4)
      tmp = Ri_bulk( zt, sst(:,:,jt), theta_zt(:,:,jt), ssq, q_zt(:,:,jt), W10(:,:,jt) )
      WRITE(6,*) ' *** Initial Bulk Richardson number:', REAL(tmp, 4)
      WRITE(6,*) ''


      IF ( l_use_cswl ) THEN

         STOP 'LOLO not like this...'

         WRITE(6,*) ''
         WRITE(6,*) '----------------------------------------------------------'
         WRITE(6,*) '          Will consider the skin temperature!'
         WRITE(6,*) ' => using cool-skin warm-layer param. in COARE and ECMWF'
         WRITE(6,*) ' => need the downwelling radiative fluxes at the surface'
         WRITE(6,*) '----------------------------------------------------------'
         WRITE(6,*) ''
         !WRITE(6,*) 'Give downwelling shortwave (solar) radiation at the surface:'
         !READ(*,*) rad_sw
         !WRITE(6,*)
         !WRITE(6,*) 'Give downwelling longwave (infrared) radiation at the surface:'
         !READ(*,*) rad_lw
         WRITE(6,*)

      END IF



      !l_wl_c36_never_called = .TRUE.

      dT(:,:,jt) = 0.  ! skin = SST for first time step
      Ts(:,:,jt) = sst(:,:,jt)

      
      DO jtt = 1, nb_itt_wl

         CALL TURB_COARE3P6( zt, zu, Ts(:,:,jt), theta_zt(:,:,jt), qs(:,:,jt), q_zt(:,:,jt), W10(:,:,jt), &
            &             Cd(:,:,jt), Ch(:,:,jt), Ce(:,:,jt), theta_zu(:,:,jt), q_zu(:,:,jt), Ublk(:,:,jt),             &
            &             xz0=zz0(:,:,jt), xu_star=zus(:,:,jt), xL=zL(:,:,jt), xUN10=zUN10(:,:,jt) )
         !! => Ts and qs are not updated: Ts=sst and qs=ssq

         
         !! Absolute temperature at zu: LOLO: Take the mean ??? => 0.5 * (t_zu + Ts) ????
         t_zu(:,:,jt) = theta_zu(:,:,jt) ! first guess...
         DO jq = 1, 4
            rgamma(:,:) = gamma_moist(t_zu(:,:,jt), q_zu(:,:,jt))
            t_zu(:,:,jt) = theta_zu(:,:,jt) - rgamma(:,:)*zu   ! Real temp.
         END DO
         
         !! Bulk Richardson Number for layer "sea-level -- zu":
         RiB(:,:,jt) = Ri_bulk(zu, Ts(:,:,jt), theta_zu(:,:,jt), qs(:,:,jt), q_zu(:,:,jt), Ublk(:,:,jt) )

         !! Air density at zu (10m)
         rho_zu(:,:,jt) = rho_air(t_zu(:,:,jt), q_zu(:,:,jt), SLP(:,:,jt))
         tmp(:,:) = SLP(:,:,jt) - rho_zu(:,:,jt)*grav*zu
         rho_zu(:,:,jt) = rho_air(t_zu, q_zu(:,:,jt), tmp(:,:))

         !! Turbulent heat fluxes
         QH  (:,:,jt) = rho_zu(:,:,jt)*tmp(1,1)*Ch(:,:,jt) * ( theta_zu(:,:,jt) - Ts(:,:,jt)  ) * Ublk(:,:,jt)
         EVAP(:,:,jt) = rho_zu(:,:,jt)*Ce(:,:,jt)          * ( qs(:,:,jt)      - q_zu(:,:,jt) ) * Ublk(:,:,jt)  ! mm/s
         QL  (:,:,jt) =  -1.* ( L_vap(Ts(:,:,jt))*EVAP(:,:,jt) )
         


         IF ( jt > 1 ) THEN ! NOT jtt !???
            
            tmp(:,:) = Ts(:,:,jt)*Ts(:,:,jt)            
            tmp(:,:) = emiss_w*(rad_lw(:,:,jt) - sigma0*tmp(:,:)*tmp(:,:)) + QH(:,:,jt) + QL(:,:,jt)
            
            IF (ldebug) THEN
               PRINT *, ' *** Non solar flux:', tmp ; PRINT *, ''
               PRINT *, ' *** Solar flux:',  (1._wp - oce_alb0)*rad_sw(:,:,jt) ; PRINT *, ''
            END IF

            CALL WL_COARE3P6( (1._wp - oce_alb0)*rad_sw, tmp, REAL(vTau/1000.,wp), sst, dT, pTau_ac, pQ_ac, isecday(jt-1), isecday(jt) )

            !PRINT *, '  => dT =', dT ; STOP

            Ts(:,:,jt) = sst(:,:,jt) + dT(:,:,jt)

         END IF

      END DO

      WRITE(6,*) ''
      WRITE(6,*) '           ---- AFTER BULK ALGO + CSWL ----'
      WRITE(6,*) ' .... to do ...'


      WRITE(6,*) ''
      WRITE(6,*) ' *** density of air at ',TRIM(czu),' => ',  rho_zu(:,:,jt), '[kg/m^3]'
      WRITE(6,*) ''
      WRITE(6,*) ' *** SST and Ts       => ',  REAL(sst(:,:,jt)-rt0,4), REAL(Ts(:,:,jt)-rt0,4), '[deg.C]'
      WRITE(6,*) ''
      WRITE(6,*) ' *** theta_zt and theta_zu => ',  REAL(theta_zt(:,:,jt)-rt0,4), REAL(theta_zu(:,:,jt)-rt0,4), '[deg.C]'
      WRITE(6,*) ''
      !theta_zt(:,:,jt)

      WRITE(6,*) ''
      WRITE(6,*) '##############################################'
   END DO


   STOP 'LULU'

   
   vBRN(ialgo,:) = REAL(RiB(1,1,:),4)
   vTheta_u(ialgo,:) = REAL(   theta_zu(1,1,:) -rt0 , 4)   ! Potential temperature at zu

   vCd(ialgo,:) = REAL(1000.*Cd(1,1,:) ,4)
   vCh(ialgo,:) = REAL(1000.*Ch(1,1,:) ,4)
   vCe(ialgo,:) = REAL(1000.*Ce(1,1,:) ,4)

   vT_u(ialgo,:) = REAL( t_zu(1,1,:) -rt0 , 4)    ! Real temp.
   vQu(ialgo,:) = REAL(  q_zu(1,1,:) , 4)

   vRho_u(ialgo,:) = REAL(rho_zu(1,1,:) ,4)

   !! Gustiness contribution:
   vUg(ialgo,:) = REAL(Ublk(1,1,:)-W10(1,1,:) , 4)

   !! z0 et u*:
   vz0(ialgo,:) = REAL(zz0(1,1,:) ,4)
   vus(ialgo,:) = REAL(zus(1,1,:) ,4)

   !zts = Ch*(theta_zu - Ts)*Ublk/zus
   !zqs = Ce*(q_zu     - qs)*Ublkblk/zus

   vL(ialgo,:) = zL(1,1,:)

   vUN10(ialgo,:) = zUN10(1,1,:)

   !! Turbulent fluxes:
   vTau(ialgo,:)  = ( rho_zu(1,1,:) * Cd(1,1,:) *           W10(1,1,:)            * Ublk(1,1,:) )*1000. ! mN/m^2

   vQH(ialgo,:)   =   QH(1,1,:)
   vEvap(ialgo,:) = EVAP(1,1,:)
   vQL(ialgo,:)   =   QL(1,1,:)

   vEvap(ialgo,:) = to_mm_p_day * vEvap(ialgo,:)  ! mm/day

   vTs(ialgo,:) = Ts(1,1,:)
   vqs(ialgo,:) = qs(1,1,:)


   WRITE(6,*) ''; WRITE(6,*) ''


   DO ialgo = 1, nb_algos

      calgob = TRIM(vca(ialgo))

      WRITE(cf_out,*) 'data_'//TRIM(calgob)//'.out'

      OPEN( UNIT=12, FILE=TRIM(cf_out), FORM='FORMATTED', RECL=1024, STATUS='unknown' )
      WRITE(12,*) '# k       date         Qsens    Qlat     SSST     Tau       WebbF   RainHF dt_skin'
      !             037 19921126231700   -41.12  -172.80   302.11    92.15      NaN      NaN  -0.238
      DO jt = 1, Nt
         WRITE(12,'(" ",i3.3," ",i14.14," ",f8.2," ",f8.2," ",f8.2," ",f8.2," ",f8.2," ",f8.2," ",f7.3)') &
            &  INT(jt,2), idate(jt), -vQH(ialgo,jt), -vQL(ialgo,jt), vTs(ialgo,jt), vTau(ialgo,jt), -999, -999, REAL(vTs(ialgo,jt)-sst(1,1,jt),4)
      END DO
      CLOSE(12)

   END DO



END PROGRAM TEST_AEROBULK_BUOY_SERIES_SKIN



SUBROUTINE usage_test()
   !!
   PRINT *,''
   PRINT *,'   List of command line options:'
   PRINT *,''
   PRINT *,' -f <ascii_file>  => file containing data'
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
