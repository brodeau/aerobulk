! AeroBulk / 2019 / L. Brodeau

!! https://www.pmel.noaa.gov/ocs/flux-documentation


PROGRAM TEST_AEROBULK_BUOY_SERIES_ICE

   USE mod_const
   USE mod_phymbl

   USE io_ezcdf     !* routines for netcdf input/output (par of SOSIE package)

   USE mod_blk_ice_nemo
   USE mod_blk_ice_an05
   USE mod_blk_ice_lg15
   !USE mod_blk_ice_best

   IMPLICIT NONE

   !INTEGER :: DISP_DEBUG

   LOGICAL, PARAMETER :: ldebug=.TRUE.
   !LOGICAL, PARAMETER :: ldebug=.FALSE.

   INTEGER, PARAMETER :: nb_algos = 3

   CHARACTER(len=800) :: cf_data='0', cunit_t, clnm_t, clndr_t

   INTEGER, PARAMETER ::   &
      &   n_dt   = 21,   &
      &   n_u    = 30,   &
      &   n_q    = 10

   CHARACTER(len=2) :: car

   INTEGER :: jt, jarg, ialgo, jq, n0, info

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
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: SIT, SLP, W10, t_zt, theta_zt, q_zt, rad_sw, rad_lw, dummy

   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: Ublk, zz0, zus, zL, zUN10

   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: t_zu, theta_zu, q_zu, rho_zu

   REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: xlon, siq, rgamma, Cp_ma, tmp


   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: Cd, Ce, Ch, QH, QL, Qsw, QNS, Qlw, RiB, TAU, SBLM

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

   nb_itt = 20 ! 20 itterations in bulk algorithm...

   !OPEN(6, FORM='formatted', RECL=512)


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

   CALL DIMS(cf_data, 'istl1', nx, ny, n0, Nt) ! Getting dimmensions from field sst...
   IF ( SUM((/nx,ny,n0/)-(/1,1,-1/)) /= 0 ) THEN
      WRITE(6,*) 'ERROR: wrong shape for field sst in input file =>', nx, ny, n0, Nt
      STOP
   END IF

   jpi = nx ; jpj = ny

   WRITE(6,*) ''
   WRITE(6,*) ' *** Allocating arrays according to nx,ny,Nt =', nx,ny,Nt
   ALLOCATE ( Ublk(nx,ny,Nt), zz0(nx,ny,Nt), zus(nx,ny,Nt), zL(nx,ny,Nt), zUN10(nx,ny,Nt) )
   ALLOCATE ( ctime(Nt), cdate(Nt), clock(Nt), chh(Nt), cmn(Nt), cldate(Nt), idate(Nt), vtime(Nt), vlon(nx) )
   ALLOCATE ( SIT(nx,ny,Nt), SLP(nx,ny,Nt), W10(nx,ny,Nt), t_zt(nx,ny,Nt), theta_zt(nx,ny,Nt), q_zt(nx,ny,Nt),  &
      &       rad_sw(nx,ny,Nt), rad_lw(nx,ny,Nt) )
   ALLOCATE ( t_zu(nx,ny,Nt), theta_zu(nx,ny,Nt), q_zu(nx,ny,Nt), rho_zu(nx,ny,Nt), dummy(nx,ny,Nt) )
   ALLOCATE ( xlon(nx,ny), siq(nx,ny), rgamma(nx,ny), Cp_ma(nx,ny), tmp(nx,ny) )
   ALLOCATE ( Cd(nx,ny,Nt), Ce(nx,ny,Nt), Ch(nx,ny,Nt), QH(nx,ny,Nt), QL(nx,ny,Nt), Qsw(nx,ny,Nt), Qlw(nx,ny,Nt), QNS(nx,ny,Nt), &
      &       RiB(nx,ny,Nt), TAU(nx,ny,Nt), SBLM(nx,ny,Nt) )

   WRITE(6,*) ' *** Allocation completed!'
   WRITE(6,*) ''


   CALL GETVAR_1D(cf_data, 'time',  vtime ) ; ! (hours since ...)
   CALL GET_VAR_INFO(cf_data, 'time', cunit_t, clnm_t,  clndr=clndr_t)
   PRINT *, 'time unit = "'//TRIM(cunit_t)//'"'
   tut_time_unit = GET_TIME_UNIT_T0( TRIM(cunit_t) ) ; ! origin
   PRINT *, ' *** Digested time unit is: ', tut_time_unit

   CALL GETVAR_1D(cf_data, 'istl1',    SIT  )  ; ! istl1 is in K !

   CALL GETVAR_1D(cf_data, 'msl',    SLP  )

   CALL GETVAR_1D(cf_data, 'u10', W10  )
   CALL GETVAR_1D(cf_data, 'v10', dummy  )
   W10 = SQRT ( W10*W10 + dummy*dummy )

   CALL GETVAR_1D(cf_data, 't2m',  t_zt )

   CALL GETVAR_1D(cf_data, 'd2m',  dummy )
   DO jt = 1, Nt
      q_zt(:,:,jt) = q_air_dp( dummy(:,:,jt), SLP(:,:,jt) )
   END DO

   CALL GETVAR_1D(cf_data, 'radsw',  rad_sw  )
   CALL GETVAR_1D(cf_data, 'radlw',  rad_lw  )


   WRITE(6,*) ''
   WRITE(6,*) ''
   !! All input time series read!


   ians=-1
   DO WHILE ( (ians<0).OR.(ians>nb_algos) )
      WRITE(6,*) 'Which algo to use? "NEMO default (v4)" => 1 | "Andreas (2005)" => 2 | "Lupkes & Gryanik (2015)" => 3 :'
      READ(*,*) ians
      IF ( ians == 1 ) calgo = 'nemo'
      IF ( ians == 2 ) calgo = 'an05'
      IF ( ians == 3 ) calgo = 'lg15'
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



   !! Time loop:
   DO jt = 1, Nt

      d_idate = time_to_date( tut_time_unit, vtime(jt),  date_prev=d_idate )

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
         WRITE(6,*) '           ---- BEFORE BULK ALGO ----'
         WRITE(6,*) ''
      END IF

      info = DISP_DEBUG(ldebug, 'scalar wind speed at '//TRIM(czu), W10(:,:,jt), '[m/s]' )

      info = DISP_DEBUG(ldebug, 'density of air at '//TRIM(czt), rho_air(t_zt(:,:,jt), q_zt(:,:,jt), SLP(:,:,jt)), '[kg/m^3]' )

      Cp_ma(:,:) = cp_air(q_zt(:,:,jt))
      info = DISP_DEBUG(ldebug, 'Cp of (moist) air at '//TRIM(czt), Cp_ma, '[J/K/kg]')

      rgamma(:,:) = gamma_moist(t_zt(:,:,jt), q_zt(:,:,jt))
      info = DISP_DEBUG(ldebug, 'Adiabatic lapse-rate of (moist) air', 1000.*rgamma, '[K/1000m]')

      info = DISP_DEBUG(ldebug, 'SIT', SIT(:,:,jt)-rt0, '[degC]')

      siq = q_sat( SIT(:,:,jt), SLP(:,:,jt), l_ice=.TRUE. )
      info = DISP_DEBUG(ldebug, 'SIQ = ', 1000.*siq, '[g/kg] (over ice!!!)')

      info = DISP_DEBUG(ldebug, 'Absolute   air temp. at '//TRIM(czt),     t_zt(:,:,jt) - rt0, '[deg.C]') ! Air temperatures at zt...

      theta_zt(:,:,jt) = t_zt(:,:,jt) + rgamma(:,:)*zt ! potential temperature at zt
      info = DISP_DEBUG(ldebug, 'Potential  air temp. at '//TRIM(czt), theta_zt(:,:,jt) - rt0, '[deg.C]')

      tmp = virt_temp(theta_zt(:,:,jt), q_zt(:,:,jt))
      info = DISP_DEBUG(ldebug, 'Virt. pot. air temp. at '//TRIM(czt),          tmp     - rt0, '[deg.C]')
      info = DISP_DEBUG(ldebug, 'Pot. temp. diff. air/ice at '//TRIM(czt),    theta_zt(:,:,jt) - SIT(:,:,jt), '[deg.C]')
      info = DISP_DEBUG(ldebug, 'Virt. pot. temp. diff. air/ice at '//TRIM(czt),    tmp - virt_temp(SIT(:,:,jt), siq), '[deg.C]')

      !! We know enough to estimate the bulk Richardson number:
      info = DISP_DEBUG(ldebug, 'Initial Bulk Richardson number', Ri_bulk( zt, SIT(:,:,jt), theta_zt(:,:,jt), siq, q_zt(:,:,jt), W10(:,:,jt) ), '[--]')

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



      Qsw(:,:,jt) = (1._wp - rice_alb0)*rad_sw(:,:,jt) ! Net solar heat flux into sea-ice...


      SELECT CASE ( TRIM(calgo) )

      CASE ( 'nemo' )
         CALL turb_ice_nemo( jt, zt, zu, SIT(:,:,jt), theta_zt(:,:,jt), siq(:,:), q_zt(:,:,jt), W10(:,:,jt),   &
            &                Cd(:,:,jt), Ch(:,:,jt), Ce(:,:,jt), theta_zu(:,:,jt), q_zu(:,:,jt), Ublk(:,:,jt),    &
            &                xz0=zz0(:,:,jt), xu_star=zus(:,:,jt), xL=zL(:,:,jt), xUN10=zUN10(:,:,jt) )
         
      CASE ( 'an05' )
         CALL turb_ice_an05( jt, zt, zu, SIT(:,:,jt), theta_zt(:,:,jt), siq(:,:), q_zt(:,:,jt), W10(:,:,jt),   &
            &                Cd(:,:,jt), Ch(:,:,jt), Ce(:,:,jt), theta_zu(:,:,jt), q_zu(:,:,jt), Ublk(:,:,jt),    &
            &                xz0=zz0(:,:,jt), xu_star=zus(:,:,jt), xL=zL(:,:,jt), xUN10=zUN10(:,:,jt) )

      CASE ( 'lg15' )
         CALL turb_ice_lg15( jt, zt, zu, SIT(:,:,jt), theta_zt(:,:,jt), siq(:,:), q_zt(:,:,jt), W10(:,:,jt),   &
            &                Cd(:,:,jt), Ch(:,:,jt), Ce(:,:,jt), theta_zu(:,:,jt), q_zu(:,:,jt), Ublk(:,:,jt),    &
            &                xz0=zz0(:,:,jt), xu_star=zus(:,:,jt), xL=zL(:,:,jt), xUN10=zUN10(:,:,jt) )

      CASE DEFAULT
         PRINT *, 'UNKNOWN algo: '//TRIM(calgo)//' !!!'
         STOP
      END SELECT


      !! Absolute temperature at zu: LOLO: Take the mean ??? => 0.5 * (t_zu + Ts) ????
      t_zu(:,:,jt) = theta_zu(:,:,jt) ! first guess...
      DO jq = 1, 4
         rgamma(:,:) = gamma_moist(t_zu(:,:,jt), q_zu(:,:,jt))
         t_zu(:,:,jt) = theta_zu(:,:,jt) - rgamma(:,:)*zu   ! Real temp.
      END DO

      !! Bulk Richardson Number for layer "sea-level -- zu":
      RiB(:,:,jt) = Ri_bulk(zu, SIT(:,:,jt), theta_zu(:,:,jt), siq(:,:), q_zu(:,:,jt), Ublk(:,:,jt) )

      !! Turbulent heat fluxes:
      CALL BULK_FORMULA( zu, SIT(:,:,jt), siq(:,:), theta_zu(:,:,jt), q_zu(:,:,jt), &
         &              Cd(:,:,jt), Ch(:,:,jt), Ce(:,:,jt), W10(:,:,jt), Ublk(:,:,jt), SLP(:,:,jt), &
         &              TAU(:,:,jt), QH(:,:,jt), QL(:,:,jt),  &
         &              pEvap=SBLM(:,:,jt), prhoa=rho_zu(:,:,jt), l_ice=.TRUE. )

      !! Longwave radiative heat fluxes:
      Qlw(:,:,jt) = qlw_net( rad_lw(:,:,jt), SIT(:,:,jt), l_ice=.TRUE. )

      !! Non-solar heat flux:
      QNS(:,:,jt) = QH(:,:,jt) + QL(:,:,jt) + Qlw(:,:,jt)


      IF (ldebug) THEN
         WRITE(6,*) ''
         WRITE(6,*) '           ---- AFTER BULK ALGO ----'
         WRITE(6,*) ''
      END IF
      info = DISP_DEBUG(ldebug, 'density of air at '//TRIM(czu), rho_zu(:,:,jt),     '[kg/m^3]' )
      info = DISP_DEBUG(ldebug, 'theta_zu',                    theta_zu(:,:,jt)-rt0, '[deg.C]'  )


      IF (ldebug) THEN
         WRITE(6,*) ''
         WRITE(6,*) '##############################################'
      END IF


   END DO !DO jt = 1, Nt






   CALL PT_SERIES(vtime(:), REAL(rho_zu(1,1,:),4), 'lolo_'//TRIM(calgo)//'.nc', 'time', &
      &           'rho_a', 'kg/m^3', 'Density of air at '//TRIM(czu), -9999._4, &
      &           ct_unit=TRIM(cunit_t), ct_clnd=TRIM(clndr_t), &
      &           vdt02=REAL(   QL(1,1,:),4), cv_dt02='Qlat',   cun02='W/m^2', cln02='Latent Heat Flux',       &
      &           vdt03=REAL(   QH(1,1,:),4), cv_dt03='Qsen',   cun03='W/m^2', cln03='Sensible Heat Flux',     &
      &           vdt04=REAL(  Qlw(1,1,:),4), cv_dt04='Qlw',    cun04='W/m^2', cln04='Net Longwave Heat Flux', &
      &           vdt05=REAL(  QNS(1,1,:),4), cv_dt05='QNS',    cun05='W/m^2', cln05='Non-solar Heat Flux',    &
      &           vdt06=REAL(  Qsw(1,1,:),4), cv_dt06='Qsw',    cun06='W/m^2', cln06='Net Solar Heat Flux',    &
      &           vdt07=REAL(  W10(1,1,:),4), cv_dt07='Wind',   cun07='m/s',   cln07='Module of Wind Speed',   &
      &           vdt08=REAL(  TAU(1,1,:),4), cv_dt08='Tau',    cun08='N/m^2', cln08='Module of Wind Stress',  &
      &           vdt09=REAL(to_mm_p_day*SBLM(1,1,:),4), cv_dt09='SBLM',       cun09='mm/day', cln09='Sublimation of ice', &
      &           vdt10=REAL( 1000.*Cd(1,1,:),4), cv_dt10='Cd', cun10='',      cln10='Drag coefficient',       &
      &           vdt11=REAL( 1000.*Ch(1,1,:),4), cv_dt11='Ch', cun11='',      cln11='Sens. Heat coeff.',      &
      &           vdt12=REAL(  zz0(1,1,:),4), cv_dt12='z0',     cun12='m',     cln12='Roughness length',       &
      &           vdt13=REAL(  RiB(1,1,:),4), cv_dt13='Rib',    cun13='m',     cln13='Bulk Richardson number', &
      &           vdt14=REAL(SIT(1,1,:)-rt0,4), cv_dt14='SIT',  cun14='deg.C', cln14='Sea-ice temperature',    &
      &           vdt15=REAL(t_zt(1,1,:)-SIT(1,1,:),4), cv_dt15='t2m-SIT', cun15='deg.C', cln15='2m air-sea temperature difference' )
   
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


END PROGRAM TEST_AEROBULK_BUOY_SERIES_ICE



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
   PRINT *,'         and use the skin temperature instead of the SIT'
   PRINT *,'         only in COARE and ECMWF'
   PRINT *,''
   PRINT *,' -h   => Show this message'
   PRINT *,''
   STOP
   !!
END SUBROUTINE usage_test
!!
