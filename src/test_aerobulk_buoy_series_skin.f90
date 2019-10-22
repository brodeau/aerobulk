! AeroBulk / 2019 / L. Brodeau

!! https://www.pmel.noaa.gov/ocs/flux-documentation


PROGRAM TEST_AEROBULK_BUOY_SERIES_SKIN

   USE mod_const
   USE mod_phymbl

   USE io_ezcdf     !* routines for netcdf input/output (par of SOSIE package)


   USE mod_blk_coare3p0
   USE mod_blk_coare3p6
   USE mod_skin_coare, ONLY: Hz_wl, Qnt_ac, Tau_ac
   !
   USE mod_blk_ecmwf
   !USE mod_skin_ecmwf, ONLY: Hz_wl !: problem: Hz_wl already exists from mod_skin_coare ...


   IMPLICIT NONE

   !INTEGER :: DISP_DEBUG

   LOGICAL, PARAMETER :: ldebug=.TRUE.
   !LOGICAL, PARAMETER :: ldebug=.FALSE.

   INTEGER, PARAMETER :: nb_algos = 3

   !INTEGER, PARAMETER :: nb_itt_wl = 2 !!LOLO

   CHARACTER(len=800) :: cf_data='0', cunit_t, clnm_t

   REAL(wp), PARAMETER ::   &
      & to_mm_p_day = 24.*3600.  !: freshwater flux: from kg/s/m^2 == mm/s to mm/day

   INTEGER, PARAMETER ::   &
      &   n_dt   = 21,   &
      &   n_u    = 30,   &
      &   n_q    = 10

   CHARACTER(len=2) :: car

   INTEGER :: jj, jt, jarg, ialgo, jq, n0, info

   INTEGER :: nx, ny, Nt, ians

   CHARACTER(len=9) :: calgo

   INTEGER(4)        :: ihh, imm, isecday_utc

   CHARACTER(len=19), DIMENSION(:), ALLOCATABLE :: ctime
   CHARACTER(len=8),  DIMENSION(:), ALLOCATABLE :: cdate
   CHARACTER(len=4),  DIMENSION(:), ALLOCATABLE :: clock
   CHARACTER(len=2),  DIMENSION(:), ALLOCATABLE :: chh, cmn ! hours and minutes
   CHARACTER(len=16), DIMENSION(:), ALLOCATABLE :: cldate ! human!
   INTEGER(8),        DIMENSION(:), ALLOCATABLE :: idate
   REAL(8),           DIMENSION(:), ALLOCATABLE :: vtime, vlon

   !! Input (or deduced from input) variables:
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: SST, SLP, W10, t_zt, theta_zt, q_zt, &
      &                                       rad_sw, rad_lw, precip, dummy

   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: Ublk, zz0, zus, zL, zUN10

   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: Ts, t_zu, theta_zu, q_zu, qs, rho_zu, dTcs, dTwl, dT, zHwl, zQac, zTac

   REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: xlon, ssq, rgamma, Cp_ma, tmp


   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: Cd, Ce, Ch, QH, QL, Qsw, QNS, Qlw, EVAP, RiB, TAU

   REAL(wp) :: zt, zu, rlon

   CHARACTER(len=3) :: czt, czu

   LOGICAL :: &
      &     l_use_rh      = .FALSE. ,   &  !: ask for RH rather than q for humidity
      &     l_use_dp      = .FALSE. ,   &  !: ask for dew-point temperature rather than q for humidity
      &     l_use_cswl    = .FALSE.        !: compute and use the skin temperature
   !                                       !: (Cool Skin Warm Layer parameterization)
   !                                       !:  => only in COARE and ECMWF


   TYPE(t_unit_t0) :: tut_time_unit
   TYPE(date)      :: d_idate

   nb_itt = 20 ! 20 itterations in bulk algorithm...

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


   !! Getting dimmension of the case and allocating arrays:
   !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   CALL DIMS(cf_data, 'sst', nx, ny, n0, Nt) ! Getting dimmensions from field sst...
   IF ( SUM((/nx,ny,n0/)-(/1,1,-1/)) /= 0 ) THEN
      WRITE(6,*) 'ERROR: wrong shape for field sst in input file =>', nx, ny, n0
      STOP
   END IF

   jpi = nx ; jpj = ny

   WRITE(6,*) ''
   WRITE(6,*) ' *** Allocating arrays according to nx,ny,Nt =', nx,ny,Nt
   ALLOCATE ( Ublk(nx,ny,Nt), zz0(nx,ny,Nt), zus(nx,ny,Nt), zL(nx,ny,Nt), zUN10(nx,ny,Nt) )
   ALLOCATE ( ctime(Nt), cdate(Nt), clock(Nt), chh(Nt), cmn(Nt), cldate(Nt), idate(Nt), vtime(Nt), vlon(nx) )
   ALLOCATE (  SST(nx,ny,Nt), SLP(nx,ny,Nt), W10(nx,ny,Nt), t_zt(nx,ny,Nt), theta_zt(nx,ny,Nt), q_zt(nx,ny,Nt),  &
      &        rad_sw(nx,ny,Nt), rad_lw(nx,ny,Nt), precip(nx,ny,Nt) )
   ALLOCATE (  Ts(nx,ny,Nt), t_zu(nx,ny,Nt), theta_zu(nx,ny,Nt), q_zu(nx,ny,Nt), qs(nx,ny,Nt), rho_zu(nx,ny,Nt), &
      &        dummy(nx,ny,Nt), dT(nx,ny,Nt), dTcs(nx,ny,Nt), dTwl(nx,ny,Nt), zHwl(nx,ny,Nt), zQac(nx,ny,Nt), zTac(nx,ny,Nt) )
   ALLOCATE (  xlon(nx,ny), ssq(nx,ny), rgamma(nx,ny), Cp_ma(nx,ny), tmp(nx,ny) )
   ALLOCATE (  Cd(nx,ny,Nt), Ce(nx,ny,Nt), Ch(nx,ny,Nt), QH(nx,ny,Nt), QL(nx,ny,Nt), Qsw(nx,ny,Nt), Qlw(nx,ny,Nt), QNS(nx,ny,Nt), &
      &        EVAP(nx,ny,Nt), RiB(nx,ny,Nt), TAU(nx,ny,Nt) )

   WRITE(6,*) ' *** Allocation completed!'
   WRITE(6,*) ''


   !! Reading data time-series into netcdf file:
   !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   CALL GETVAR_1D(cf_data, 'lon',  vlon ) ; ! (longitude for solar time...)
   rlon = vlon(1)

   DO jj = 1, ny
      xlon(:,jj) = vlon(:)
   END DO

   PRINT *, ' *** Longitude of station: ', rlon


   CALL GETVAR_1D(cf_data, 'time',  vtime ) ; ! (hours since ...)
   CALL GET_VAR_INFO(cf_data, 'time', cunit_t, clnm_t)
   PRINT *, 'time unit = "'//TRIM(cunit_t)//'"'
   tut_time_unit = GET_TIME_UNIT_T0( TRIM(cunit_t) ) ; ! origin
   PRINT *, ' *** Digested time unit is: ', tut_time_unit


   CALL GETVAR_1D(cf_data, 'sst',    SST  )
   SST  = SST + rt0

   CALL GETVAR_1D(cf_data, 'slp',    SLP  )
   SLP = SLP * 100.

   CALL GETVAR_1D(cf_data, 'wndspd', W10  )

   CALL GETVAR_1D(cf_data, 't_air',  t_zt )
   t_zt = t_zt + rt0

   CALL GETVAR_1D(cf_data, 'rh_air', dummy)
   dummy = MIN(99.999 , dummy)
   DO jt = 1, Nt
      q_zt(:,:,jt) = q_air_rh(0.01*dummy(:,:,jt), t_zt(:,:,jt), SLP(:,:,jt))
   END DO

   CALL GETVAR_1D(cf_data, 'rad_sw',  rad_sw  )
   CALL GETVAR_1D(cf_data, 'rad_lw',  rad_lw  )


   ians=0
   DO WHILE ( (ians<1).OR.(ians>nb_algos) )
      WRITE(6,*) 'Which algo to use? "coare3p0" => 1 , "ecmwf" => 2 , "coare3p6" => 3 :'
      READ(*,*) ians
      IF ( ians == 1 ) calgo = 'coare3p0 '
      IF ( ians == 2 ) calgo = 'ecmwf    '
      IF ( ians == 3 ) calgo = 'coare3p6 '
   END DO
   WRITE(6,*) '  ==> your choice: ', TRIM(calgo)
   WRITE(6,*) ''




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


   !IF (ldebug) THEN
   !   WRITE(6,*) '*       idate     ,   wind    ,       SST    ,     t_zt     ,      q_zt      ,    rad_sw     , rad_lw  :'
   !   DO jt = 1, Nt
   !      WRITE(6,*) vtime(jt), REAL(W10(:,:,jt),4), REAL(SST(:,:,jt),4), REAL(t_zt(:,:,jt),4), REAL(q_zt(:,:,jt),4), REAL(rad_sw(:,:,jt),4), REAL(rad_lw(:,:,jt),4)
   !   END DO
   !END IF


   !! Some initializations:

   ialgo = 1

   zz0 = 0.
   zus = 0.
   zL  = 0.
   zUN10 = 0.

   isecday_utc = 0

   dT(:,:,:)    = 0.  ! skin = SST for first time step
   dTcs(:,:,:) = 0.
   dTwl(:,:,:) = 0.
   zHwl(:,:,:)  = 0.
   zQac(:,:,:)  = 0.
   zTac(:,:,:)  = 0.



   !! Time loop:
   DO jt = 1, Nt

      d_idate = time_to_date( tut_time_unit, vtime(jt) )

      ihh     = d_idate%hour
      imm     = d_idate%minute
      WRITE(cldate(jt),'(i4.4,"/",i2.2,"/",i2.2,"-",i2.2,":",i2.2)') d_idate%year,  d_idate%month,  d_idate%day,  ihh, imm
      isecday_utc = ihh*3600 + imm*60 ! UTC time in seconds since midnight...

      IF (ldebug) THEN
         WRITE(6,*) ''; WRITE(6,*) ''
         WRITE(6,*) '##########################################################'
      END IF
      WRITE(6,'(" #### Time = ",a," => isecday_utc = ",i6.6)') cldate(jt), isecday_utc
      IF (ldebug) THEN
         WRITE(6,*) '##########################################################'
         WRITE(6,*) ''
         WRITE(6,*) '           ---- BEFORE BULK ALGO + CSWL ----'
         WRITE(6,*) ''
      END IF

      info = DISP_DEBUG(ldebug, 'scalar wind speed at '//TRIM(czu), W10(:,:,jt), '[m/s]' )

      info = DISP_DEBUG(ldebug, 'density of air at '//TRIM(czt), rho_air(t_zt(:,:,jt), q_zt(:,:,jt), SLP(:,:,jt)), '[kg/m^3]' )

      Cp_ma(:,:) = cp_air(q_zt(:,:,jt))
      info = DISP_DEBUG(ldebug, 'Cp of (moist) air at '//TRIM(czt), Cp_ma, '[J/K/kg]')

      rgamma(:,:) = gamma_moist(t_zt(:,:,jt), q_zt(:,:,jt))
      info = DISP_DEBUG(ldebug, 'Adiabatic lapse-rate of (moist) air', 1000.*rgamma, '[K/1000m]')

      info = DISP_DEBUG(ldebug, 'SST', SST(:,:,jt)-rt0, '[degC]')

      ssq = rdct_qsat_salt*q_sat(SST(:,:,jt), SLP(:,:,jt))
      info = DISP_DEBUG(ldebug, 'SSQ = 0.98*q_sat(SST)', 1000.*ssq, '[g/kg]')

      info = DISP_DEBUG(ldebug, 'Absolute   air temp. at '//TRIM(czt),     t_zt(:,:,jt) - rt0, '[deg.C]') ! Air temperatures at zt...

      theta_zt(:,:,jt) = t_zt(:,:,jt) + rgamma(:,:)*zt ! potential temperature at zt
      info = DISP_DEBUG(ldebug, 'Potential  air temp. at '//TRIM(czt), theta_zt(:,:,jt) - rt0, '[deg.C]')

      tmp = virt_temp(theta_zt(:,:,jt), q_zt(:,:,jt))
      info = DISP_DEBUG(ldebug, 'Virt. pot. air temp. at '//TRIM(czt),          tmp     - rt0, '[deg.C]')
      info = DISP_DEBUG(ldebug, 'Pot. temp. diff. air/sea at '//TRIM(czt),    theta_zt(:,:,jt) - SST(:,:,jt), '[deg.C]')
      info = DISP_DEBUG(ldebug, 'Virt. pot. temp. diff. air/sea at '//TRIM(czt),    tmp - virt_temp(SST(:,:,jt), ssq), '[deg.C]')

      !! We know enough to estimate the bulk Richardson number:
      info = DISP_DEBUG(ldebug, 'Initial Bulk Richardson number', Ri_bulk( zt, SST(:,:,jt), theta_zt(:,:,jt), ssq, q_zt(:,:,jt), W10(:,:,jt) ), '[--]')

      !STOP

      !IF ( l_use_cswl ) THEN
      !   STOP 'LOLO not like this...'
      !   WRITE(6,*) ''
      !   WRITE(6,*) '----------------------------------------------------------'
      !   WRITE(6,*) '          Will consider the skin temperature!'
      !   WRITE(6,*) ' => using cool-skin warm-layer param. in COARE and ECMWF'
      !   WRITE(6,*) ' => need the downwelling radiative fluxes at the surface'
      !   WRITE(6,*) '----------------------------------------------------------'
      !   WRITE(6,*) ''
      !   !WRITE(6,*) 'Give downwelling shortwave (solar) radiation at the surface:'
      !   !READ(*,*) rad_sw
      !   !WRITE(6,*)
      !   !WRITE(6,*) 'Give downwelling longwave (infrared) radiation at the surface:'
      !   !READ(*,*) rad_lw
      !   WRITE(6,*)
      !END IF


      !l_wl_c36_never_called = .TRUE.


      Ts(:,:,jt) = SST(:,:,jt)
      qs(:,:,jt) = ssq(:,:)


      Qsw(:,:,jt) = (1._wp - oce_alb0)*rad_sw(:,:,jt) ! Net solar heat flux into the ocean

      IF ( TRIM(calgo) == 'coare3p6' ) THEN
         CALL TURB_COARE3P6( jt, zt, zu, Ts(:,:,jt), theta_zt(:,:,jt), qs(:,:,jt), q_zt(:,:,jt), W10(:,:,jt), .TRUE., .TRUE.,  &
            &             Cd(:,:,jt), Ch(:,:,jt), Ce(:,:,jt), theta_zu(:,:,jt), q_zu(:,:,jt), Ublk(:,:,jt),    &
            &             Qsw=Qsw(:,:,jt), rad_lw=rad_lw(:,:,jt), slp=SLP(:,:,jt), pdt_cs=dTcs(:,:,jt),       & ! for cool-skin !
            &             isecday_utc=isecday_utc, plong=xlon(:,:), pdT_wl=dTwl(:,:,jt), pHz_wl=zHwl(:,:,jt), &
            &             xz0=zz0(:,:,jt), xu_star=zus(:,:,jt), xL=zL(:,:,jt), xUN10=zUN10(:,:,jt) )

         zQac(:,:,jt) = Qnt_ac(:,:) ! known from module "mod_skin_coare"
         zTac(:,:,jt) = Tau_ac(:,:) !                "

      ELSEIF( TRIM(calgo) == 'ecmwf'    ) THEN
         CALL TURB_ECMWF(   jt, zt, zu, Ts(:,:,jt), theta_zt(:,:,jt), qs(:,:,jt), q_zt(:,:,jt), W10(:,:,jt), .TRUE., .TRUE.,  &
            &             Cd(:,:,jt), Ch(:,:,jt), Ce(:,:,jt), theta_zu(:,:,jt), q_zu(:,:,jt), Ublk(:,:,jt),  &
            &             Qsw=Qsw(:,:,jt), rad_lw=rad_lw(:,:,jt), slp=SLP(:,:,jt), pdt_cs=dTcs(:,:,jt),     & ! for cool-skin !
            &             pdT_wl=dTwl(:,:,jt), pHz_wl=zHwl(:,:,jt),        &
            &             xz0=zz0(:,:,jt), xu_star=zus(:,:,jt), xL=zL(:,:,jt), xUN10=zUN10(:,:,jt) )

      ELSE
         PRINT *, 'UNKNOWN algo: '//TRIM(calgo)//' !!!'
         STOP
      END IF

      dT(:,:,jt) = Ts(:,:,jt) - SST(:,:,jt)

      !! => Ts and qs ARE updated, but only for cool-skin !!!!
      !IF (jtt == 1)  dTcs(:,:,jt) = Ts(:,:,jt) - SST(:,:,jt)

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
      rho_zu(:,:,jt) = rho_air(t_zu(:,:,jt), q_zu(:,:,jt), tmp(:,:))

      !! Turbulent heat fluxes:
      tmp(:,:) = cp_air(q_zu(:,:,jt))
      QH  (:,:,jt) =  rho_zu(:,:,jt)*tmp*Ch(:,:,jt) * ( theta_zu(:,:,jt) - Ts(:,:,jt) ) * Ublk(:,:,jt)
      EVAP(:,:,jt) = MAX( -rho_zu(:,:,jt)    *Ce(:,:,jt) * (     q_zu(:,:,jt) - qs(:,:,jt) ) * Ublk(:,:,jt) , 0._wp )  ! mm/s
      QL  (:,:,jt) =  -1.* ( L_vap(Ts(:,:,jt))*EVAP(:,:,jt) )

      TAU(:,:,jt)  = rho_zu(:,:,jt) * Cd(:,:,jt) * Ublk(:,:,jt)*Ublk(:,:,jt)

      !! Longwave radiative heat fluxes:
      tmp(:,:) = Ts(:,:,jt)*Ts(:,:,jt)
      Qlw(:,:,jt) = emiss_w*(rad_lw(:,:,jt) - stefan*tmp(:,:)*tmp(:,:))

      QNS(:,:,jt) = QH(:,:,jt) + QL(:,:,jt) + Qlw(:,:,jt) ! Non-solar component of net heat flux !





      IF (ldebug) THEN
         WRITE(6,*) ''
         WRITE(6,*) '           ---- AFTER BULK ALGO + CSWL ----'
         WRITE(6,*) ''
      END IF
      info = DISP_DEBUG(ldebug, 'density of air at '//TRIM(czu), rho_zu(:,:,jt),     '[kg/m^3]' )
      info = DISP_DEBUG(ldebug, 'theta_zu',                    theta_zu(:,:,jt)-rt0, '[deg.C]'  )


      IF (ldebug) THEN
         WRITE(6,*) ''
         WRITE(6,*) '##############################################'
      END IF


   END DO !DO jt = 1, Nt


   IF     ( TRIM(calgo) == 'coare3p6' ) THEN

      CALL PT_SERIES(vtime(:), REAL(rho_zu(1,1,:),4), 'lolo_'//TRIM(calgo)//'.nc', 'time', &
         &           'rho_a', 'kg/m^3', 'Density of air at '//TRIM(czu), -9999._4, &
         &           ct_unit=TRIM(cunit_t), &
         &           vdt02=REAL(   QL(1,1,:),4), cv_dt02='Qlat',  cun02='W/m^2', cln02='Latent Heat Flux',       &
         &           vdt03=REAL(   QH(1,1,:),4), cv_dt03='Qsen',  cun03='W/m^2', cln03='Sensible Heat Flux',     &
         &           vdt04=REAL(  Qlw(1,1,:),4), cv_dt04='Qlw',   cun04='W/m^2', cln04='Net Longwave Heat Flux', &
         &           vdt05=REAL(  QNS(1,1,:),4), cv_dt05='QNS',   cun05='W/m^2', cln05='Non-solar Heat Flux',    &
         &           vdt06=REAL(  Qsw(1,1,:),4), cv_dt06='Qsw',   cun06='W/m^2', cln06='Net Solar Heat Flux',    &
         &           vdt07=REAL(dTcs(1,1,:),4), cv_dt07='dTcs', cun07='deg.C', cln07='Cool-Skin dT',           &
         &           vdt08=REAL(dTwl(1,1,:),4), cv_dt08='dTwl', cun08='deg.C', cln08='Warm-Layer dT',          &
         &           vdt09=REAL(  W10(1,1,:),4), cv_dt09='Wind',  cun09='m/s',   cln09='Module of Wind Speed',   &
         &           vdt10=REAL(  TAU(1,1,:),4), cv_dt10='Tau',   cun10='N/m^2', cln10='Module of Wind Stress',  &
         &           vdt11=REAL(   dT(1,1,:),4), cv_dt11='dT',    cun11='deg.C', cln11='SST - Ts',               &
         &           vdt12=REAL( zHwl(1,1,:),4), cv_dt12='H_wl',  cun12='m',     cln12='Estimated depth of warm-layer', &
         &           vdt13=REAL( zQac(1,1,:),4), cv_dt13='Qnt_ac',cun13='J/m2',  cln13='Accumulated absorbed heat in WL', &
         &           vdt14=REAL( zTac(1,1,:),4), cv_dt14='Tau_ac',cun14='N.s/m2',cln14='Accumulated absorbed momentum in WL' )

   ELSE
      
      CALL PT_SERIES(vtime(:), REAL(rho_zu(1,1,:),4), 'lolo_'//TRIM(calgo)//'.nc', 'time', &
         &           'rho_a', 'kg/m^3', 'Density of air at '//TRIM(czu), -9999._4, &
         &           ct_unit=TRIM(cunit_t), &
         &           vdt02=REAL(   QL(1,1,:),4), cv_dt02='Qlat',  cun02='W/m^2', cln02='Latent Heat Flux',       &
         &           vdt03=REAL(   QH(1,1,:),4), cv_dt03='Qsen',  cun03='W/m^2', cln03='Sensible Heat Flux',     &
         &           vdt04=REAL(  Qlw(1,1,:),4), cv_dt04='Qlw',   cun04='W/m^2', cln04='Net Longwave Heat Flux', &
         &           vdt05=REAL(  QNS(1,1,:),4), cv_dt05='QNS',   cun05='W/m^2', cln05='Non-solar Heat Flux',    &
         &           vdt06=REAL(  Qsw(1,1,:),4), cv_dt06='Qsw',   cun06='W/m^2', cln06='Net Solar Heat Flux',    &
         &           vdt07=REAL(dTcs(1,1,:),4),  cv_dt07='dTcs',  cun07='deg.C', cln07='Cool-Skin dT',           &
         &           vdt08=REAL(dTwl(1,1,:),4),  cv_dt08='dTwl',  cun08='deg.C', cln08='Warm-Layer dT',          &
         &           vdt09=REAL(  W10(1,1,:),4), cv_dt09='Wind',  cun09='m/s',   cln09='Module of Wind Speed',   &
         &           vdt10=REAL(  TAU(1,1,:),4), cv_dt10='Tau',   cun10='N/m^2', cln10='Module of Wind Stress',  &
         &           vdt11=REAL(   dT(1,1,:),4), cv_dt11='dT',    cun11='deg.C', cln11='SST - Ts',               &
         &           vdt12=REAL( zHwl(1,1,:),4), cv_dt12='H_wl',  cun12='m',     cln12='Estimated depth of warm-layer'  )
      
   END IF
   
   WRITE(6,*) ''; WRITE(6,*) ''



CONTAINS

   FUNCTION DISP_DEBUG( ldbg, cstr, rval, cunit )
      INTEGER :: DISP_DEBUG
      LOGICAL,                  INTENT(in) :: ldbg
      CHARACTER(len=*),         INTENT(in) :: cstr
      REAL(wp), DIMENSION(:,:), INTENT(in) :: rval
      CHARACTER(len=*),         INTENT(in) :: cunit
      !!
      DISP_DEBUG = 0
      IF ( ldbg ) THEN
         !WRITE(6,*) ' *** '//TRIM(cstr)
         !WRITE(6,*) ' *** '//TRIM(cstr), ' => ', REAL(rval(1,1),4), ' '//TRIM(cunit)
         WRITE(6,'(" *** ",a40," => ",f12.4," ",a9)') TRIM(cstr),  REAL(rval(1,1),4), TRIM(cunit)
         !WRITE(6,*) ''
         DISP_DEBUG = 1
      END IF
   END FUNCTION DISP_DEBUG


END PROGRAM TEST_AEROBULK_BUOY_SERIES_SKIN



SUBROUTINE usage_test()
   !!
   PRINT *,''
   PRINT *,'   List of command line options:'
   PRINT *,''
   PRINT *,' -f <netcdf_file>  => file containing data'
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
