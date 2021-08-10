! AeroBulk / 2019 / L. Brodeau

!! https://www.pmel.noaa.gov/ocs/flux-documentation


PROGRAM TEST_AEROBULK_BUOY_SERIES_OCE

   USE mod_const
   USE mod_phymbl

   USE io_ezcdf     !* routines for netcdf input/output (par of SOSIE package)

   USE mod_blk_ncar
   USE mod_blk_coare3p0
   USE mod_blk_coare3p6
   USE mod_skin_coare, ONLY: Qnt_ac, Tau_ac ! Hz_wl

   !
   USE mod_blk_ecmwf
   !USE mod_skin_ecmwf, ONLY: Hz_wl !: problem: Hz_wl already exists from mod_skin_coare ...

   USE mod_blk_andreas

   IMPLICIT NONE

   !INTEGER :: DISP_DEBUG



   LOGICAL, PARAMETER :: lverbose = .TRUE.
   LOGICAL, PARAMETER :: ldebug   = .TRUE.

   INTEGER, PARAMETER :: nb_algos = 4

   !INTEGER, PARAMETER :: nb_iter_wl = 2 !!LOLO

   CHARACTER(len=800) :: cf_data='0', cn_exp='0', cunit_t, clnm_t, clndr_t

   CHARACTER(len=80 ) :: cv_time

   CHARACTER(len=80) :: csep='#################################################################################'

   INTEGER, PARAMETER ::   &
      &   n_dt   = 21,   &
      &   n_u    = 30,   &
      &   n_q    = 10

   CHARACTER(len=2) :: car

   INTEGER :: jj, jt, jarg, ialgo, jq, info, ifi=0, ivi=0

   INTEGER, PARAMETER :: nx = 1, ny = 1

   INTEGER :: Nt, ians, nd0, nd1, nd2

   CHARACTER(len=10) :: calgo

   INTEGER(4)        :: ihh, imm, isecday_utc

   CHARACTER(len=16), DIMENSION(:), ALLOCATABLE :: cldate ! human!
   REAL(8),           DIMENSION(:), ALLOCATABLE :: vtime, vlon

   !! Input (or deduced from input) variables:
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: SST, SLP, W10, t_zt, theta_zt, q_zt, &
      &                                       rad_sw, rad_lw, precip, dummy

   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: Ublk, zz0, zus, zL, zUN10

   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: Ts, t_zu, theta_zu, q_zu, qs, rho_zu, dTcs, dTwl, dT, zHwl, zQac, zTac

   REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: xlon, ssq, rgamma, Cp_ma, xtmp, X3


   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: Cd, Ce, Ch, QH, QL, Qsw, QNS, Qlw, EVAP, RiB, TAU

   REAL(wp) :: zt, zu, rlon

   CHARACTER(len=3) :: czt, czu

   LOGICAL :: &
      &   l_3x3_ts = .FALSE., &   !: fields in input netcdf are 3x3 in space (come from NEMO STATION_ASF!)
      &   l_hum_rh = .FALSE., &   !: humidity in NetCDF file is Relative Humidity [%]
      &   l_wndspd = .FALSE.      !: wind speed read instead of deduced from "u" and "v"


   TYPE(t_unit_t0) :: tut_time_unit
   TYPE(date)      :: d_idate

   nb_iter = 20 ! 20 itterations in bulk algorithm...

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

      CASE('-n')
         jarg = jarg + 1
         CALL get_command_ARGUMENT(jarg,cn_exp)

      CASE('-3')
         l_3x3_ts = .TRUE.

      CASE('-r')
         l_hum_rh = .TRUE.

      CASE('-w')
         l_wndspd = .TRUE.

      CASE DEFAULT
         WRITE(6,*) 'Unknown option: ', trim(car) ; WRITE(6,*) ''
         CALL usage_test()

      END SELECT

   END DO


   IF ( TRIM(cf_data) == '0' ) CALL usage_test()
   IF ( TRIM(cn_exp)  == '0' ) CALL usage_test()


   WRITE(6,*) ''
   WRITE(6,*) ' *** Input file is ',TRIM(cf_data)
   WRITE(6,*) ''
   WRITE(6,*) ' *** Name of experiment is ',TRIM(cn_exp)
   WRITE(6,*) ''
   WRITE(6,*) '  *** Epsilon            = Rd/Rv       (~0.622) =>', reps0
   WRITE(6,*) '  *** Virt. temp. const. = (1-eps)/eps (~0.608) =>', rctv0
   WRITE(6,*) ''

   WRITE(6,*) ''


   !! Getting dimmension of the case and allocating arrays:
   !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   !! Number of time records?
   cv_time = 'time'
   IF( l_3x3_ts ) cv_time = 'time_counter'
   CALL DIMS(cf_data, cv_time, Nt,  nd0, nd1, nd2 ) ! Getting dimmensions from field sst...
   PRINT *, ' *** Number of time records to treat: ', Nt


   WRITE(6,*) ''
   WRITE(6,*) ' *** Allocating arrays according to nx,ny,Nt =', INT( (/nx,ny,Nt/), 2 )
   ALLOCATE ( Ublk(nx,ny,Nt), zz0(nx,ny,Nt), zus(nx,ny,Nt), zL(nx,ny,Nt), zUN10(nx,ny,Nt) )
   ALLOCATE (   cldate(Nt), vtime(Nt) )
   ALLOCATE (  SST(nx,ny,Nt), SLP(nx,ny,Nt), W10(nx,ny,Nt), t_zt(nx,ny,Nt), theta_zt(nx,ny,Nt), q_zt(nx,ny,Nt),  &
      &        rad_sw(nx,ny,Nt), rad_lw(nx,ny,Nt), precip(nx,ny,Nt) )
   ALLOCATE (  Ts(nx,ny,Nt), t_zu(nx,ny,Nt), theta_zu(nx,ny,Nt), q_zu(nx,ny,Nt), qs(nx,ny,Nt), rho_zu(nx,ny,Nt), &
      &        dummy(nx,ny,Nt), dT(nx,ny,Nt), dTcs(nx,ny,Nt), dTwl(nx,ny,Nt), zHwl(nx,ny,Nt), zQac(nx,ny,Nt), zTac(nx,ny,Nt) )
   ALLOCATE (  xlon(nx,ny), ssq(nx,ny), rgamma(nx,ny), Cp_ma(nx,ny), xtmp(nx,ny) )
   ALLOCATE (  Cd(nx,ny,Nt), Ce(nx,ny,Nt), Ch(nx,ny,Nt), QH(nx,ny,Nt), QL(nx,ny,Nt), Qsw(nx,ny,Nt), Qlw(nx,ny,Nt), QNS(nx,ny,Nt), &
      &        EVAP(nx,ny,Nt), RiB(nx,ny,Nt), TAU(nx,ny,Nt) )

   WRITE(6,*) ' *** Allocation completed!'
   WRITE(6,*) ''


   CALL GETVAR_1D(cf_data, cv_time,  vtime ) ; ! (hours since ...)
   CALL GET_VAR_INFO(cf_data, cv_time, cunit_t, clnm_t,  clndr=clndr_t)
   PRINT *, 'time unit = "'//TRIM(cunit_t)//'"'
   tut_time_unit = GET_TIME_UNIT_T0( TRIM(cunit_t) ) ; ! origin
   PRINT *, ' *** Digested time unit is: '
   PRINT *, tut_time_unit
   PRINT *, ''



   CALL set_variable_names_default()

   IF( .NOT. l_3x3_ts ) THEN

      !! "1D + t" AeroBulk convention:
      !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ALLOCATE( vlon(nx) )
      CALL GETVAR_1D(cf_data, 'lon',  vlon ) ; ! (longitude for solar time...)
      rlon = vlon(1)
      DO jj = 1, ny
         xlon(:,jj) = vlon(:)
      END DO
      DEALLOCATE( vlon )

      CALL GETVAR_1D(cf_data, cv_sst,    SST  )

      CALL GETVAR_1D(cf_data, cv_patm,    SLP  ) ; ! must be in [Pa]

      IF ( l_wndspd ) THEN
         CALL GETVAR_1D(cf_data, cv_wndspd, W10  )
      ELSE
         CALL GETVAR_1D(cf_data, cv_u_wnd, W10  )
         CALL GETVAR_1D(cf_data, cv_v_wnd, dummy  )
         W10 = SQRT ( W10*W10 + dummy*dummy )
      END IF

      CALL GETVAR_1D(cf_data, cv_t_air,  t_zt ) ; ! must be in [deg.C]
      CALL TO_KELVIN_3D(t_zt, cname=TRIM(cv_t_air) )

      IF ( l_hum_rh ) THEN
         !! Relative humidity is read:
         CALL GETVAR_1D(cf_data, cv_rh_air, dummy)
         dummy = MIN(99.999 , dummy)
         DO jt = 1, Nt
            !PRINT *, 'LOLO: rh, t_zt, SLP =', dummy(:,:,jt), t_zt(:,:,jt), SLP(:,:,jt)
            q_zt(:,:,jt) = q_air_rh(dummy(:,:,jt), t_zt(:,:,jt), SLP(:,:,jt))
         END DO
      ELSE
         !! Dew-point is read:
         CALL GETVAR_1D(cf_data, cv_dp_air,  dummy )
         DO jt = 1, Nt
            q_zt(:,:,jt) = q_air_dp( dummy(:,:,jt), SLP(:,:,jt) )
         END DO
      END IF

      CALL GETVAR_1D(cf_data, cv_radsw,  rad_sw  )
      CALL GETVAR_1D(cf_data, cv_radlw,  rad_lw  )

   ELSE

      !! "2D (3x3) + t" NEMO convention:
      !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !CALL GETVAR_1D_R8_3x3_to_1x1

      CALL set_variable_names_ecmwf()  ! ECMWF !

      ALLOCATE( X3(3,3) )
      CALL GETVAR_2D( ifi, ivi, cf_data, 'nav_lon', 1, 0, 0, X3 )
      xlon(:,:) = X3(2,2)
      rlon = xlon(1,1)
      DEALLOCATE( X3 )

      CALL GETVAR_1D_R8_3x3_to_1x1(cf_data, cv_sst,    SST  )

      CALL GETVAR_1D_R8_3x3_to_1x1(cf_data, cv_patm,    SLP  ) ; ! must be in [Pa]

      IF ( l_wndspd ) THEN
         CALL GETVAR_1D_R8_3x3_to_1x1(cf_data, cv_wndspd, W10  )
      ELSE
         CALL GETVAR_1D_R8_3x3_to_1x1(cf_data, cv_u_wnd, W10  )
         CALL GETVAR_1D_R8_3x3_to_1x1(cf_data, cv_v_wnd, dummy  )
         W10 = SQRT ( W10*W10 + dummy*dummy )
      END IF

      CALL GETVAR_1D_R8_3x3_to_1x1(cf_data, cv_v_wnd, dummy  )
      W10 = SQRT ( W10*W10 + dummy*dummy )

      CALL GETVAR_1D_R8_3x3_to_1x1(cf_data, cv_t_air,  t_zt )
      CALL TO_KELVIN_3D(t_zt, cname=TRIM(cv_t_air) )

      IF ( l_hum_rh ) THEN
         !! Relative humidity is read:
         CALL GETVAR_1D_R8_3x3_to_1x1(cf_data, cv_rh_air, dummy)
         dummy = MIN(99.999 , dummy)
         DO jt = 1, Nt
            !PRINT *, 'LOLO: rh, t_zt, SLP =', dummy(:,:,jt), t_zt(:,:,jt), SLP(:,:,jt)
            q_zt(:,:,jt) = q_air_rh(dummy(:,:,jt), t_zt(:,:,jt), SLP(:,:,jt))
         END DO
      ELSE
         !! Dew-point is read:
         CALL GETVAR_1D_R8_3x3_to_1x1(cf_data, cv_dp_air,  dummy )
         DO jt = 1, Nt
            q_zt(:,:,jt) = q_air_dp( dummy(:,:,jt), SLP(:,:,jt) )
         END DO
      END IF



      CALL GETVAR_1D_R8_3x3_to_1x1(cf_data, cv_radsw,  rad_sw  )
      CALL GETVAR_1D_R8_3x3_to_1x1(cf_data, cv_radlw,  rad_lw  )

   END IF !IF( .NOT. l_3x3_ts )


   CALL TO_KELVIN_3D(SST, cname=TRIM(cv_sst) )



   PRINT *, ' *** Longitude of station: ', rlon

   WRITE(6,*) ''
   WRITE(6,*) ''
   !! All input time series read!


   ians=-1
   DO WHILE ( (ians<0).OR.(ians>nb_algos) )
      WRITE(6,*) 'Which algo to use? "ncar" => 0 , "coare3p0" => 1 , "ecmwf" => 2 , "coare3p6" => 3, "andreas" => 4 :'
      READ(*,*) ians
      IF ( ians == 0 ) calgo = 'ncar     '
      IF ( ians == 1 ) calgo = 'coare3p0 '
      IF ( ians == 2 ) calgo = 'ecmwf    '
      IF ( ians == 3 ) calgo = 'coare3p6 '
      IF ( ians == 4 ) calgo = 'andreas  '
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



   !! Some initializations:

   ialgo = 0
   zz0 = 0.
   zus = 0.
   zL  = 0.
   zUN10 = 0.

   isecday_utc = 0

   d_idate%year   = 0
   d_idate%month  = 0
   d_idate%day    = 0
   d_idate%hour   = 0
   d_idate%minute = 0
   d_idate%second = 0

   dT(:,:,:)   = 0.  ! skin = SST for first time step
   dTcs(:,:,:) = 0.
   dTwl(:,:,:) = 0.
   zHwl(:,:,:) = 0.
   zQac(:,:,:) = 0.
   zTac(:,:,:) = 0.



   !! Time loop:
   DO jt = 1, Nt

      d_idate = time_to_date( tut_time_unit, vtime(jt),  date_prev=d_idate )

      ihh     = d_idate%hour
      imm     = d_idate%minute
      WRITE(cldate(jt),'(i4.4,"/",i2.2,"/",i2.2,"-",i2.2,":",i2.2)') d_idate%year,  d_idate%month,  d_idate%day,  ihh, imm
      isecday_utc = ihh*3600 + imm*60 ! UTC time in seconds since midnight...

      IF (lverbose) THEN
         WRITE(6,*) ''; WRITE(6,*) ''
         WRITE(6,*) csep
      END IF
      WRITE(6,'(" #### Time = ",a," => isecday_utc = ",i6.6," => jt = ",i4.4)') cldate(jt), isecday_utc, jt
      IF (lverbose) THEN
         WRITE(6,*) csep
         WRITE(6,*) ''
         WRITE(6,*) '           ---- BEFORE BULK ALGO + CSWL ----'
         WRITE(6,*) ''
      END IF

      info = DISP_DEBUG(lverbose, 'SST', SST(:,:,jt)-rt0, '[degC]')

      info = DISP_DEBUG(lverbose, 'atmospheric pressure', SLP(:,:,jt), '[Pa]' )

      info = DISP_DEBUG(lverbose, 'Absolute   air temp. at '//TRIM(czt),     t_zt(:,:,jt) - rt0, '[deg.C]') ! Air temperatures at zt...

      info = DISP_DEBUG(lverbose, 'Specific humidity of air at '//TRIM(czt), q_zt(:,:,jt),       '[kg/kg]')

      rgamma(:,:) = gamma_moist(t_zt(:,:,jt), q_zt(:,:,jt))
      info = DISP_DEBUG(lverbose, 'Adiabatic lapse-rate of (moist) air', 1000.*rgamma, '[K/1000m]')

      theta_zt(:,:,jt) = t_zt(:,:,jt) + rgamma(:,:)*zt ! potential temperature at zt
      info = DISP_DEBUG(lverbose, 'Potential  air temp. at '//TRIM(czt), theta_zt(:,:,jt) - rt0, '[deg.C]')

      xtmp = virt_temp(theta_zt(:,:,jt), q_zt(:,:,jt))
      info = DISP_DEBUG(lverbose, 'Virt. pot. air temp. at '//TRIM(czt),          xtmp     - rt0, '[deg.C]')

      info = DISP_DEBUG(lverbose, 'density of air at '//TRIM(czt), rho_air(t_zt(:,:,jt), q_zt(:,:,jt), SLP(:,:,jt)), '[kg/m^3]' )

      info = DISP_DEBUG(lverbose, 'scalar wind speed at '//TRIM(czu), W10(:,:,jt), '[m/s]' )

      Cp_ma(:,:) = cp_air(q_zt(:,:,jt))
      info = DISP_DEBUG(lverbose, 'Cp of (moist) air at '//TRIM(czt), Cp_ma, '[J/K/kg]')


      ssq = rdct_qsat_salt*q_sat(SST(:,:,jt), SLP(:,:,jt))
      info = DISP_DEBUG(lverbose, 'SSQ = 0.98*q_sat(SST)', 1000.*ssq, '[g/kg]')

      info = DISP_DEBUG(lverbose, 'Pot. temp. diff. air/sea at '//TRIM(czt),    theta_zt(:,:,jt) - SST(:,:,jt), '[deg.C]')
      info = DISP_DEBUG(lverbose, 'Virt. pot. temp. diff. air/sea at '//TRIM(czt),    xtmp - virt_temp(SST(:,:,jt), ssq), '[deg.C]')

      !! We know enough to estimate the bulk Richardson number:
      info = DISP_DEBUG(lverbose, 'Initial Bulk Richardson number', Ri_bulk( zt, SST(:,:,jt), theta_zt(:,:,jt), ssq, q_zt(:,:,jt), W10(:,:,jt) ), '[--]')

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


      Qsw(:,:,jt) = (1._wp - roce_alb0)*rad_sw(:,:,jt) ! Net solar heat flux into the ocean


      SELECT CASE ( TRIM(calgo) )

      CASE ( 'ncar' )
         CALL TURB_NCAR    (     zt, zu, Ts(:,:,jt), theta_zt(:,:,jt), qs(:,:,jt), q_zt(:,:,jt), W10(:,:,jt),  &
            &             Cd(:,:,jt), Ch(:,:,jt), Ce(:,:,jt), theta_zu(:,:,jt), q_zu(:,:,jt), Ublk(:,:,jt),    &
            &             xz0=zz0(:,:,jt), xu_star=zus(:,:,jt), xL=zL(:,:,jt), xUN10=zUN10(:,:,jt) )

      CASE ( 'coare3p0' )
         CALL TURB_COARE3P6( jt, zt, zu, Ts(:,:,jt), theta_zt(:,:,jt), qs(:,:,jt), q_zt(:,:,jt), W10(:,:,jt), .TRUE., .TRUE.,  &
            &             Cd(:,:,jt), Ch(:,:,jt), Ce(:,:,jt), theta_zu(:,:,jt), q_zu(:,:,jt), Ublk(:,:,jt),    &
            &             Qsw=Qsw(:,:,jt), rad_lw=rad_lw(:,:,jt), slp=SLP(:,:,jt), pdt_cs=dTcs(:,:,jt),       & ! for cool-skin !
            &             isecday_utc=isecday_utc, plong=xlon(:,:), pdT_wl=dTwl(:,:,jt), pHz_wl=zHwl(:,:,jt), &
            &             xz0=zz0(:,:,jt), xu_star=zus(:,:,jt), xL=zL(:,:,jt), xUN10=zUN10(:,:,jt) )

         zQac(:,:,jt) = Qnt_ac(:,:) ! known from module "mod_skin_coare"
         zTac(:,:,jt) = Tau_ac(:,:) !                "

      CASE ( 'coare3p6' )
         CALL TURB_COARE3P6( jt, zt, zu, Ts(:,:,jt), theta_zt(:,:,jt), qs(:,:,jt), q_zt(:,:,jt), W10(:,:,jt), .TRUE., .TRUE.,  &
            &             Cd(:,:,jt), Ch(:,:,jt), Ce(:,:,jt), theta_zu(:,:,jt), q_zu(:,:,jt), Ublk(:,:,jt),    &
            &             Qsw=Qsw(:,:,jt), rad_lw=rad_lw(:,:,jt), slp=SLP(:,:,jt), pdt_cs=dTcs(:,:,jt),       & ! for cool-skin !
            &             isecday_utc=isecday_utc, plong=xlon(:,:), pdT_wl=dTwl(:,:,jt), pHz_wl=zHwl(:,:,jt), &
            &             xz0=zz0(:,:,jt), xu_star=zus(:,:,jt), xL=zL(:,:,jt), xUN10=zUN10(:,:,jt) )

         zQac(:,:,jt) = Qnt_ac(:,:) ! known from module "mod_skin_coare"
         zTac(:,:,jt) = Tau_ac(:,:) !                "

      CASE ( 'ecmwf'    )
         CALL TURB_ECMWF(   jt, zt, zu, Ts(:,:,jt), theta_zt(:,:,jt), qs(:,:,jt), q_zt(:,:,jt), W10(:,:,jt), .TRUE., .TRUE.,  &
            &             Cd(:,:,jt), Ch(:,:,jt), Ce(:,:,jt), theta_zu(:,:,jt), q_zu(:,:,jt), Ublk(:,:,jt),  &
            &             Qsw=Qsw(:,:,jt), rad_lw=rad_lw(:,:,jt), slp=SLP(:,:,jt), pdt_cs=dTcs(:,:,jt),     & ! for cool-skin !
            &             pdT_wl=dTwl(:,:,jt), pHz_wl=zHwl(:,:,jt),        &
            &             xz0=zz0(:,:,jt), xu_star=zus(:,:,jt), xL=zL(:,:,jt), xUN10=zUN10(:,:,jt) )

      CASE ( 'andreas' )
         CALL TURB_ANDREAS(    zt, zu, Ts(:,:,jt), theta_zt(:,:,jt), qs(:,:,jt), q_zt(:,:,jt), W10(:,:,jt),  &
            &             Cd(:,:,jt), Ch(:,:,jt), Ce(:,:,jt), theta_zu(:,:,jt), q_zu(:,:,jt), Ublk(:,:,jt),    &
            &             xz0=zz0(:,:,jt), xu_star=zus(:,:,jt), xL=zL(:,:,jt), xUN10=zUN10(:,:,jt) )


      CASE DEFAULT
         PRINT *, 'UNKNOWN algo: '//TRIM(calgo)//' !!!'
         STOP
      END SELECT

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

      !! Turbulent heat fluxes:
      CALL BULK_FORMULA( zu, Ts(:,:,jt), qs(:,:,jt), theta_zu(:,:,jt), q_zu(:,:,jt), &
         &              Cd(:,:,jt), Ch(:,:,jt), Ce(:,:,jt), W10(:,:,jt), Ublk(:,:,jt), SLP(:,:,jt), &
         &              TAU(:,:,jt), QH(:,:,jt), QL(:,:,jt),  &
         &              pEvap=EVAP(:,:,jt), prhoa=rho_zu(:,:,jt) )

      !! Longwave radiative heat fluxes:
      Qlw(:,:,jt) = qlw_net( rad_lw(:,:,jt), Ts(:,:,jt) )

      !! Non-solar heat flux:
      QNS(:,:,jt) = QH(:,:,jt) + QL(:,:,jt) + Qlw(:,:,jt)


      IF (lverbose) THEN
         WRITE(6,*) ''
         WRITE(6,*) '           ---- AFTER BULK ALGO + CSWL ----'
         WRITE(6,*) ''
      END IF
      info = DISP_DEBUG(lverbose, 'density of air at '//TRIM(czu), rho_zu(:,:,jt),     '[kg/m^3]' )
      info = DISP_DEBUG(lverbose, 'theta_zu',                    theta_zu(:,:,jt)-rt0, '[deg.C]'  )


      IF (lverbose) THEN
         WRITE(6,*) ''
         WRITE(6,*) csep
      END IF


   END DO !DO jt = 1, Nt


   IF     ( TRIM(calgo) == 'coare3p6' ) THEN

      CALL PT_SERIES(vtime(:), REAL(rho_zu(1,1,:),4), TRIM(cn_exp)//'_'//TRIM(calgo)//'.nc', cv_time, &
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
         &           vdt14=REAL( zTac(1,1,:),4), cv_dt14='Tau_ac',cun14='N.s/m2',cln14='Accumulated absorbed momentum in WL', &
         &           vdt15=REAL(   Cd(1,1,:),4), cv_dt15='Cd',    cun15='',      cln15='Drag coefficient', &
         &           vdt16=REAL(   Ce(1,1,:),4), cv_dt16='Ce',    cun16='',      cln16='Evap. coefficient', &
         &           vdt17=REAL(   Ch(1,1,:),4), cv_dt17='Ch',    cun17='',      cln17='Sens. heat. coefficient'        )

   ELSE

      CALL PT_SERIES(vtime(:), REAL(rho_zu(1,1,:),4), TRIM(cn_exp)//'_'//TRIM(calgo)//'.nc', cv_time, &
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
         &           vdt12=REAL( zHwl(1,1,:),4), cv_dt12='H_wl',  cun12='m',     cln12='Estimated depth of warm-layer', &
         &           vdt13=REAL(   Cd(1,1,:),4), cv_dt13='Cd',    cun13='',      cln13='Drag coefficient', &
         &           vdt14=REAL(   Ce(1,1,:),4), cv_dt14='Ce',    cun14='',      cln14='Evap. coefficient', &
         &           vdt15=REAL(   Ch(1,1,:),4), cv_dt15='Ch',    cun15='',      cln15='Sens. heat. coefficient'        )


   END IF

   WRITE(6,*) ''; WRITE(6,*) ''
   CLOSE(6)

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
         WRITE(6,'(" *** ",a48,f12.4," ",a)') &
            &       TRIM(cstr)//' => '//ACHAR(27)//'[91m ',  REAL(rval(1,1),4), ACHAR(27)//'[0m'//TRIM(cunit)
         DISP_DEBUG = 1
      END IF
   END FUNCTION DISP_DEBUG

   SUBROUTINE usage_test()
      !!
      PRINT *,''
      PRINT *,'   List of command line options:'
      PRINT *, '  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      PRINT *,''
      PRINT *,' -f <file>  => NetCDF file containing input data'
      PRINT *,''
      PRINT *,' -n <name>  => name/label for the experiment'
      PRINT *,''
      PRINT *,' -3         => fields in input netcdf are 3x3 in space (those of NEMO/tests/STATION_ASF!)'
      PRINT *,''
      PRINT *,' -r         => humidity in NetCDF file is Relative Humidity [%]'
      PRINT *,''
      PRINT *,' -w         => wind speed read in NetCDF file rather than "u" and "v"'
      PRINT *,''
      PRINT *,' -h         => Show this message'
      PRINT *,''
      STOP
      !!
   END SUBROUTINE usage_test


END PROGRAM TEST_AEROBULK_BUOY_SERIES_OCE
