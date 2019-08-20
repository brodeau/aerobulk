! AeroBulk / 2019 / L. Brodeau

!! https://www.pmel.noaa.gov/ocs/flux-documentation


PROGRAM TEST_AEROBULK_BUOY_SERIES_SKIN

   USE mod_const
   USE mod_phymbl

   USE io_ezcdf     !* routines for netcdf input/output (par of SOSIE package)

   !USE mod_blk_coare3p0
   USE mod_blk_coare3p6
   USE mod_wl_coare3p6
   !USE mod_blk_ncar
   !USE mod_blk_ecmwf

   IMPLICIT NONE

   !INTEGER :: DISP_DEBUG

   LOGICAL, PARAMETER :: ldebug=.TRUE.
   !LOGICAL, PARAMETER :: ldebug=.FALSE.

   INTEGER, PARAMETER :: nb_algos = 1

   INTEGER, PARAMETER :: nb_itt_wl = 2 !!LOLO

   CHARACTER(len=800) :: cf_data='0', cblabla, cf_out='output.dat', cunit_t, clnm_t

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

   INTEGER :: jt, jarg, jl, ialgo, jq, jtt, n0, info

   INTEGER :: nx, ny, Nt

   CHARACTER(len=19) :: cdt
   INTEGER(4)        :: iclock, ihh, imm, isecday_n, isecday_b
   CHARACTER(len=19), DIMENSION(:), ALLOCATABLE :: ctime
   CHARACTER(len=8),  DIMENSION(:), ALLOCATABLE :: cdate
   CHARACTER(len=4),  DIMENSION(:), ALLOCATABLE :: clock
   CHARACTER(len=2),  DIMENSION(:), ALLOCATABLE :: chh, cmn ! hours and minutes
   CHARACTER(len=16), DIMENSION(:), ALLOCATABLE :: cldate ! human!
   INTEGER(8),        DIMENSION(:), ALLOCATABLE :: idate
   REAL(8),           DIMENSION(:), ALLOCATABLE :: vtime

   !! Input (or deduced from input) variables:
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: SST, SLP, W10, t_zt, theta_zt, q_zt, &
      &                                       rad_sw, rad_lw, precip, rlat, rlon, dummy

   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: Ublk, zz0, zus, zL, zUN10

   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: Ts, t_zu, theta_zu, q_zu, qs, rho_zu, &
      &                                       dT, pTau_ac, pQ_ac

   REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: ssq, rgamma, Cp_ma, tmp


   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: Cd, Ce, Ch, QH, QL, Qsw, QNS, Qlw, EVAP, RiB, TAU

   !REAL(4), DIMENSION(:,:),   ALLOCATABLE ::  &
   !   &           vCd, vCe, vCh, vTheta_u, vT_u, vQu, vz0, vus, vRho_u, vUg, vL, vBRN, &
   !   &           vUN10, vQL, vTau, vQH, vEvap, vTs, vqs



   REAL(wp) :: zt, zu

   CHARACTER(len=3) :: czt, czu

   LOGICAL :: &
      &     l_use_rh      = .FALSE. ,   &  !: ask for RH rather than q for humidity
      &     l_use_dp      = .FALSE. ,   &  !: ask for dew-point temperature rather than q for humidity
      &     l_use_cswl    = .FALSE.        !: compute and use the skin temperature
   !                                       !: (Cool Skin Warm Layer parameterization)
   !                                       !:  => only in COARE and ECMWF


   TYPE(t_unit_t0) :: tut_time_unit
   TYPE(date)      :: d_idate



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
   ALLOCATE ( ctime(Nt), cdate(Nt), clock(Nt), chh(Nt), cmn(Nt), cldate(Nt), idate(Nt), vtime(Nt) )
   ALLOCATE (  SST(nx,ny,Nt), SLP(nx,ny,Nt), W10(nx,ny,Nt), t_zt(nx,ny,Nt), theta_zt(nx,ny,Nt), q_zt(nx,ny,Nt),  &
      &        rad_sw(nx,ny,Nt), rad_lw(nx,ny,Nt), precip(nx,ny,Nt), rlat(nx,ny,Nt), rlon(nx,ny,Nt) )
   ALLOCATE (  Ts(nx,ny,Nt), t_zu(nx,ny,Nt), theta_zu(nx,ny,Nt), q_zu(nx,ny,Nt), qs(nx,ny,Nt), rho_zu(nx,ny,Nt), dummy(nx,ny,Nt), &
      &        dT(nx,ny,Nt), pTau_ac(nx,ny,Nt), pQ_ac(nx,ny,Nt) )
   ALLOCATE (   ssq(nx,ny), rgamma(nx,ny), Cp_ma(nx,ny), tmp(nx,ny) )
   ALLOCATE (  Cd(nx,ny,Nt), Ce(nx,ny,Nt), Ch(nx,ny,Nt), QH(nx,ny,Nt), QL(nx,ny,Nt), Qsw(nx,ny,Nt), Qlw(nx,ny,Nt), QNS(nx,ny,Nt), &
      &        EVAP(nx,ny,Nt), RiB(nx,ny,Nt), TAU(nx,ny,Nt) )

   !ALLOCATE ( vCd(nb_algos,Nt), vCe(nb_algos,Nt), vCh(nb_algos,Nt), vTheta_u(nb_algos,Nt), vT_u(nb_algos,Nt), vQu(nb_algos,Nt), &
   !   &       vz0(nb_algos,Nt), vus(nb_algos,Nt), vRho_u(nb_algos,Nt), vUg(nb_algos,Nt), vL(nb_algos,Nt), vBRN(nb_algos,Nt),    &
   !   &       vUN10(nb_algos,Nt), vQL(nb_algos,Nt), vTau(nb_algos,Nt), vQH(nb_algos,Nt), vEvap(nb_algos,Nt), vTs(nb_algos,Nt),  &
   !   &       vqs(nb_algos,Nt) )
   WRITE(6,*) ' *** Allocation completed!'
   WRITE(6,*) ''

   !! Reading data time-series into netcdf file:
   !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

   !PRINT *, rad_lw

   !STOP


   !OPEN(11, FILE=TRIM(cf_data), FORM='formatted', STATUS='old')
   !READ(11,*) cblabla ! first commented line...
   !DO jt = 1, Nt
   !   READ(11,*) ctime(jt), W10(1,1,jt), SST(1,1,jt), t_zt(1,1,jt), q_zt(1,1,jt), rad_sw(1,1,jt), rad_lw(1,1,jt), precip(1,1,jt), rlat(1,1,jt), rlon(1,1,jt), dummy(1,1,jt)
   !END DO
   !CLOSE(11)


   !DO jt = 1, Nt
   !   cdt = ctime(jt)
   !   cdate(jt) = cdt(1:8)
   !   clock(jt) = cdt(9:13)
   !   chh(jt)   = cdt(9:10)
   !   cmn(jt)   = cdt(11:12)
   !   READ(clock(jt),'(i4.4)') iclock
   !   READ(chh(jt),'(i2.2)') ihh
   !   READ(cmn(jt),'(i2.2)') imm
   !   WRITE(cldate(jt),'(a4,"/",a2,"/",a2,"-",a2,":",a2)') cdt(1:4), cdt(5:6), cdt(7:8), chh(jt), cmn(jt)
   !   IF (ldebug) PRINT *, ' *** date = ', cldate(jt)
   !END DO
   !IF (ldebug) PRINT *, ''


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



   IF (ldebug) THEN
      !                   19921125132100   4.700000       302.1500       300.8500      1.7600000E-02  0.0000000E+00   428.0000
      WRITE(6,*) '*       idate     ,   wind    ,       SST    ,     t_zt     ,      q_zt      ,    rad_sw     , rad_lw  :'
      DO jt = 1, Nt
         !WRITE(6,*) cldate(jt), REAL(W10(:,:,jt),4), REAL(SST(:,:,jt),4), REAL(t_zt(:,:,jt),4), REAL(q_zt(:,:,jt),4), REAL(rad_sw(:,:,jt),4), REAL(rad_lw(:,:,jt),4)
         WRITE(6,*) vtime(jt), REAL(W10(:,:,jt),4), REAL(SST(:,:,jt),4), REAL(t_zt(:,:,jt),4), REAL(q_zt(:,:,jt),4), REAL(rad_sw(:,:,jt),4), REAL(rad_lw(:,:,jt),4)
      END DO
   END IF


   !! Some initializations:

   ialgo = 1

   zz0 = 0.
   zus = 0.
   zL  = 0.
   zUN10 = 0.

   pTau_ac = 0.
   pQ_ac   = 0.


   isecday_b = 0.
   isecday_n = 0.


   !! Time loop:
   DO jt = 1, Nt

      ihh = d_idate%hour
      imm = d_idate%minute

      d_idate = time_to_date( tut_time_unit, vtime(jt) )
      WRITE(cldate(jt),'(i4.4,"/",i2.2,"/",i2.2,"-",i2.2,":",i2.2)') d_idate%year,  d_idate%month,  d_idate%day,  ihh, imm

      isecday_n = ihh*3600 + imm*60

      IF (ldebug) THEN
         WRITE(6,*) ''; WRITE(6,*) ''
         WRITE(6,*) '##############################################'
      END IF
      WRITE(6,*) '#### Time = ', cldate(jt), ' (seconds since start of day:',isecday_n,')'
      IF (ldebug) THEN
         WRITE(6,*) '##############################################'
         WRITE(6,*) ''
         WRITE(6,*) '           ---- BEFORE BULK ALGO + CSWL ----'
         WRITE(6,*) ''
      END IF

      info = DISP_DEBUG(ldebug, 'density of air at '//TRIM(czt), rho_air(t_zt(:,:,jt), q_zt(:,:,jt), SLP(:,:,jt)), '[kg/m^3]' )

      Cp_ma(:,:) = cp_air(q_zt(:,:,jt))
      info = DISP_DEBUG(ldebug, 'Cp of (moist) air at '//TRIM(czt), Cp_ma, '[J/K/kg]')

      rgamma(:,:) = gamma_moist(t_zt(:,:,jt), q_zt(:,:,jt))
      info = DISP_DEBUG(ldebug, 'Adiabatic lapse-rate of (moist) air', 1000.*rgamma, '[K/1000m]')

      info = DISP_DEBUG(ldebug, 'SST', SST(:,:,jt)-rt0, '[degC]')

      ssq = rdct_qsat_salt*q_sat(SST(:,:,jt), SLP(:,:,jt))
      info = DISP_DEBUG(ldebug, 'SSQ = 0.98*q_sat(SST)', 1000.*ssq, '[g/kg]')

      ! Air temperatures at zt:
      info = DISP_DEBUG(ldebug, 'Absolute   air temp. at '//TRIM(czt),     t_zt(:,:,jt) - rt0, '[deg.C]')

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

      dT(:,:,jt) = 0.  ! skin = SST for first time step
      Ts(:,:,jt) = SST(:,:,jt)
      qs(:,:,jt) = ssq(:,:)


      !LOLO: DO jtt = 1, nb_itt_wl

      CALL TURB_COARE3P6( zt, zu, Ts(:,:,jt), theta_zt(:,:,jt), qs(:,:,jt), q_zt(:,:,jt), W10(:,:,jt), &
         &             Cd(:,:,jt), Ch(:,:,jt), Ce(:,:,jt), theta_zu(:,:,jt), q_zu(:,:,jt), Ublk(:,:,jt),             &
         &             xz0=zz0(:,:,jt), xu_star=zus(:,:,jt), xL=zL(:,:,jt), xUN10=zUN10(:,:,jt) )
      !! => Ts and qs are not updated: Ts=SST and qs=ssq


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

      !! Turbulent heat fluxes:
      tmp(:,:) = cp_air(q_zu(:,:,jt))
      QH  (:,:,jt) =  rho_zu(:,:,jt)*tmp*Ch(:,:,jt) * ( theta_zu(:,:,jt) - Ts(:,:,jt) ) * Ublk(:,:,jt)
      EVAP(:,:,jt) = -rho_zu(:,:,jt)    *Ce(:,:,jt) * (     q_zu(:,:,jt) - qs(:,:,jt) ) * Ublk(:,:,jt)  ! mm/s
      QL  (:,:,jt) =  -1.* ( L_vap(Ts(:,:,jt))*EVAP(:,:,jt) )

      TAU(:,:,jt)  = rho_zu(:,:,jt) * Cd(:,:,jt) * Ublk(:,:,jt)*Ublk(:,:,jt)

      !! Radiative heat fluxes:
      Qsw(:,:,jt) = (1._wp - oce_alb0)*rad_sw(:,:,jt)
      tmp(:,:) = Ts(:,:,jt)*Ts(:,:,jt)
      Qlw(:,:,jt) = emiss_w*(rad_lw(:,:,jt) - sigma0*tmp(:,:)*tmp(:,:))

      QNS(:,:,jt) = QH(:,:,jt) + QL(:,:,jt) + Qlw(:,:,jt) ! Non-solar component of net heat flux !

      !
      !IF ( jt > 1 ) THEN ! NOT jtt !???
      !
      !   tmp(:,:) = Ts(:,:,jt)*Ts(:,:,jt)
      !   tmp(:,:) = emiss_w*(rad_lw(:,:,jt) - sigma0*tmp(:,:)*tmp(:,:)) + QH(:,:,jt) + QL(:,:,jt)
      !
      !   IF (ldebug) THEN
      !      PRINT *, ' *** Non solar flux:', tmp ; PRINT *, ''
      !      PRINT *, ' *** Solar flux:',  (1._wp - oce_alb0)*rad_sw(:,:,jt) ; PRINT *, ''
      !   END IF
      !


      !CALL WL_COARE3P6( Qsw(:,:,jt), QNS(:,:,jt), TAU(:,:,jt), SST(:,:,jt), dT, pTau_ac, pQ_ac, isecday_b, isecday_n )
      !
      !   PRINT *, '  => dT =', dT ; STOP
      !
      !   Ts(:,:,jt) = sst(:,:,jt) + dT(:,:,jt)
      !
      !END IF

      !LOLO: END DO

      IF (ldebug) THEN
         WRITE(6,*) ''
         WRITE(6,*) '           ---- AFTER BULK ALGO + CSWL ----'
         WRITE(6,*) ''
      END IF

      info = DISP_DEBUG(ldebug, 'density of air at '//TRIM(czu), rho_zu(:,:,jt),     '[kg/m^3]' )
      info = DISP_DEBUG(ldebug, 'Ts',                                Ts(:,:,jt)-rt0, '[deg.C]'  )
      info = DISP_DEBUG(ldebug, 'qs',                          1000.*qs(:,:,jt),     '[g/kg]'   )
      info = DISP_DEBUG(ldebug, 'theta_zu',                    theta_zu(:,:,jt)-rt0, '[deg.C]'  )
      info = DISP_DEBUG(ldebug, 'dT_skin',                    Ts(:,:,jt)-sst(:,:,jt), '[deg.C]'  )

      IF (ldebug) THEN
         WRITE(6,*) ''
         WRITE(6,*) '##############################################'
      END IF


      isecday_b = isecday_n

   END DO

   CALL PT_SERIES(vtime(:), REAL(rho_zu(1,1,:),4), 'lolo.nc', 'time', &
      &           'rho_a', 'kg/m^3', 'Density of air at '//TRIM(czu), -9999._4, &
      &           ct_unit=TRIM(cunit_t), &
      &           vdt2=REAL(  QL(1,1,:),4), cv_dt2='Qlat', cun2='W/m^2', cln2='Latent Heat Flux',       &
      &           vdt3=REAL(  QH(1,1,:),4), cv_dt3='Qsen', cun3='W/m^2', cln3='Sensible Heat Flux',     &
      &           vdt4=REAL( Qlw(1,1,:),4), cv_dt4='Qlw',  cun4='W/m^2', cln4='Net Longwave Heat Flux', &
      &           vdt5=REAL( QNS(1,1,:),4), cv_dt5='QNS',  cun5='W/m^2', cln5='Non-solar Heat Flux',    &
      &           vdt6=REAL( Qsw(1,1,:),4), cv_dt6='Qsw',  cun6='W/m^2', cln6='Net Solar Heat Flux',    &
      &           vdt7=REAL(EVAP(1,1,:),4), cv_dt7='E',    cun7='mm/s',  cln7='Evaporation',            &
      &           vdt8=REAL(TAU(1,1,:),4),  cv_dt8='Tau',  cun8='N/m^2', cln8='Module of Wind Stress'   )

   !,             &



   STOP 'LULU'


   WRITE(6,*) ''; WRITE(6,*) ''


   !DO ialgo = 1, nb_algos
   !
   !   calgob = TRIM(vca(ialgo))
   !
   !   WRITE(cf_out,*) 'data_'//TRIM(calgob)//'.out'
   !
   !   OPEN( UNIT=12, FILE=TRIM(cf_out), FORM='FORMATTED', RECL=1024, STATUS='unknown' )
   !   WRITE(12,*) '# k       date         Qsens    Qlat     SSST     Tau       WebbF   RainHF dt_skin'
   !   !             037 19921126231700   -41.12  -172.80   302.11    92.15      NaN      NaN  -0.238
   !   DO jt = 1, Nt
   !      WRITE(12,'(" ",i3.3," ",i14.14," ",f8.2," ",f8.2," ",f8.2," ",f8.2," ",f8.2," ",f8.2," ",f7.3)') &
   !         &  INT(jt,2), idate(jt), -vQH(ialgo,jt), -vQL(ialgo,jt), vTs(ialgo,jt), vTau(ialgo,jt), -999, -999, REAL(vTs(ialgo,jt)-sst(1,1,jt),4)
   !   END DO
   !   CLOSE(12)
   !
   !END DO

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
         WRITE(6,'(" *** ",a40," => ",f10.4," ",a9)') TRIM(cstr),  REAL(rval(1,1),4), TRIM(cunit)
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
