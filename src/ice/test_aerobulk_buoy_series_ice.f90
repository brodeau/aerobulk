! AeroBulk / 2019 / L. Brodeau

!! https://www.pmel.noaa.gov/ocs/flux-documentation


PROGRAM TEST_AEROBULK_BUOY_SERIES_ICE

   USE mod_const
   USE mod_phymbl

   USE io_ezcdf     !* routines for netcdf input/output (par of SOSIE package)

   USE mod_blk_ice_nemo
   USE mod_blk_ice_an05
   USE mod_blk_ice_lg15
   USE mod_blk_ice_lu12

   IMPLICIT NONE

   !INTEGER :: DISP_DEBUG

   LOGICAL, PARAMETER :: lverbose = .TRUE.
   LOGICAL, PARAMETER :: ldebug   = .FALSE.
   INTEGER, PARAMETER :: jtdbg    = 1   ! if (ldebug) that's the time step we'll start from...
   !INTEGER, PARAMETER :: jtdbg    = 1447   ! if (ldebug) that's the time step we'll start from...

   INTEGER, PARAMETER :: nb_algos = 4

   CHARACTER(len=800) :: cf_data='0', cn_exp='0', cunit_t, clnm_t, clndr_t
   CHARACTER(len=80 ) :: cv_time

   CHARACTER(len=80) :: csep='#################################################################################'

   CHARACTER(len=2) :: car

   INTEGER :: jt, jarg, ialgo, jq, info

   INTEGER, PARAMETER :: nx = 1, ny = 1

   INTEGER :: Nt, ians, nd0, nd1, nd2

   CHARACTER(len=10) :: calgo

   INTEGER(4)        :: ihh, imm, isecday_utc

   CHARACTER(len=19), DIMENSION(:), ALLOCATABLE :: ctime
   CHARACTER(len=8),  DIMENSION(:), ALLOCATABLE :: cdate
   CHARACTER(len=4),  DIMENSION(:), ALLOCATABLE :: clock
   CHARACTER(len=2),  DIMENSION(:), ALLOCATABLE :: chh, cmn ! hours and minutes
   CHARACTER(len=16), DIMENSION(:), ALLOCATABLE :: cldate ! human!
   INTEGER(8),        DIMENSION(:), ALLOCATABLE :: idate
   REAL(8),           DIMENSION(:), ALLOCATABLE :: vtime, vlon

   !! Input (or deduced from input) variables:
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: SIT, SST, SLP, W10, SIC, t_zt, theta_zt, q_zt, rad_sw, rad_lw, dummy

   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: Ublk, zz0, zus, zL, zUN10, zCdN

   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: t_zu, theta_zu, q_zu, rho_zt, rho_zu

   REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: xlon, SIQ, SSQ, rgamma, Cp_ma, tmp

   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: Cd_i, Ce_i, Ch_i, QH, QL, Qsw, QNS, Qlw, RiB_zt, RiB_zu, TAU, SBLM

   REAL(wp) :: zt, zu

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
   WRITE(6,*) ' *** Allocating arrays according to nx,ny,Nt =', INT( (/nx,ny,Nt/),2 )
   ALLOCATE ( Ublk(nx,ny,Nt), zz0(nx,ny,Nt), zus(nx,ny,Nt), zL(nx,ny,Nt), zUN10(nx,ny,Nt), zCdN(nx,ny,Nt) )
   ALLOCATE ( ctime(Nt), cdate(Nt), clock(Nt), chh(Nt), cmn(Nt), cldate(Nt), idate(Nt), vtime(Nt), vlon(nx) )
   ALLOCATE ( SIT(nx,ny,Nt), SST(nx,ny,Nt), SLP(nx,ny,Nt), W10(nx,ny,Nt), t_zt(nx,ny,Nt), theta_zt(nx,ny,Nt), q_zt(nx,ny,Nt),  &
      &       rad_sw(nx,ny,Nt), rad_lw(nx,ny,Nt), SIC(nx,ny,Nt) )
   ALLOCATE ( t_zu(nx,ny,Nt), theta_zu(nx,ny,Nt), q_zu(nx,ny,Nt), rho_zt(nx,ny,Nt), rho_zu(nx,ny,Nt), dummy(nx,ny,Nt) )
   ALLOCATE ( xlon(nx,ny), SIQ(nx,ny), SSQ(nx,ny), rgamma(nx,ny), Cp_ma(nx,ny), tmp(nx,ny) )
   ALLOCATE ( Cd_i(nx,ny,Nt), Ce_i(nx,ny,Nt), Ch_i(nx,ny,Nt), QH(nx,ny,Nt), QL(nx,ny,Nt), Qsw(nx,ny,Nt), Qlw(nx,ny,Nt), QNS(nx,ny,Nt), &
      &       RiB_zt(nx,ny,Nt), RiB_zu(nx,ny,Nt), TAU(nx,ny,Nt), SBLM(nx,ny,Nt) )

   t_zu = 0. ; theta_zu = 0. ; q_zu = 0. ; rho_zt = 0. ; rho_zu = 0. ; dummy = 0.
   Cd_i = 0. ; Ce_i = 0. ; Ch_i = 0. ; QH = 0. ; QL = 0. ; Qsw = 0. ; Qlw = 0. ; QNS = 0.
   RiB_zt = 0. ; RiB_zu = 0. ; TAU = 0. ; SBLM = 0.

   WRITE(6,*) ' *** Allocation completed!'
   WRITE(6,*) ''


   CALL GETVAR_1D(cf_data, cv_time,  vtime ) ; ! (hours since ...)
   CALL GET_VAR_INFO(cf_data, cv_time, cunit_t, clnm_t,  clndr=clndr_t)
   PRINT *, 'time unit = "'//TRIM(cunit_t)//'"'
   tut_time_unit = GET_TIME_UNIT_T0( TRIM(cunit_t) ) ; ! origin
   PRINT *, ' *** Digested time unit is: '
   PRINT *, tut_time_unit
   PRINT *, ''



   IF( .NOT. l_3x3_ts ) THEN

      !! "1D + t" AeroBulk convention:
      !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL GETVAR_1D(cf_data, 'siconc',   SIC  )  ; ! sea-ice concentration (0-1)

      CALL GETVAR_1D(cf_data, 'istl1',    SIT  )  ; ! istl1 is in K !

      CALL GETVAR_1D(cf_data, 'sst',    SST  )  ; ! sst is in K !

      CALL GETVAR_1D(cf_data, 'msl',    SLP  )

      CALL GETVAR_1D(cf_data, 'u10', W10  )
      CALL GETVAR_1D(cf_data, 'v10', dummy  )
      W10 = SQRT ( W10*W10 + dummy*dummy )

      CALL GETVAR_1D(cf_data, 't2m',  t_zt )
      CALL TO_KELVIN_3D(t_zt, cname=TRIM('t2m') )

      IF ( l_hum_rh ) THEN
         !! Relative humidity is read:
         CALL GETVAR_1D(cf_data, 'rh_air', dummy)
         dummy = MIN(99.999 , dummy)
         DO jt = 1, Nt
            !PRINT *, 'LOLO: rh, t_zt, SLP =', dummy(:,:,jt), t_zt(:,:,jt), SLP(:,:,jt)
            q_zt(:,:,jt) = q_air_rh(dummy(:,:,jt), t_zt(:,:,jt), SLP(:,:,jt))
         END DO
      ELSE
         !! Dew-point is read:
         CALL GETVAR_1D(cf_data, 'd2m',  dummy )
         DO jt = 1, Nt
            q_zt(:,:,jt) = q_air_dp( dummy(:,:,jt), SLP(:,:,jt) )
         END DO
      END IF

      CALL GETVAR_1D(cf_data, 'ssrd',  rad_sw  )
      CALL GETVAR_1D(cf_data, 'strd',  rad_lw  )


   ELSE

      !! "2D (3x3) + t" NEMO convention:
      !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      CALL GETVAR_1D_R8_3x3_to_1x1(cf_data, 'siconc',   SIC  )  ; ! sea-ice concentration (0-1)

      CALL GETVAR_1D_R8_3x3_to_1x1(cf_data, 'istl1',    SIT  )  ; ! istl1 is in K !

      CALL GETVAR_1D_R8_3x3_to_1x1(cf_data, 'sst',    SST  )  ; ! sst is in K !

      CALL GETVAR_1D_R8_3x3_to_1x1(cf_data, 'msl',    SLP  ) ; ! must be in [Pa]

      CALL GETVAR_1D_R8_3x3_to_1x1(cf_data, 'u10', W10  )
      CALL GETVAR_1D_R8_3x3_to_1x1(cf_data, 'v10', dummy  )
      W10 = SQRT ( W10*W10 + dummy*dummy )

      CALL GETVAR_1D_R8_3x3_to_1x1(cf_data, 't2m',  t_zt )
      CALL TO_KELVIN_3D(t_zt, cname=TRIM('t2m') )

      IF ( l_hum_rh ) THEN
         !! Relative humidity is read:
         CALL GETVAR_1D_R8_3x3_to_1x1(cf_data, 'rh_air', dummy)
         dummy = MIN(99.999 , dummy)
         DO jt = 1, Nt
            !PRINT *, 'LOLO: rh, t_zt, SLP =', dummy(:,:,jt), t_zt(:,:,jt), SLP(:,:,jt)
            q_zt(:,:,jt) = q_air_rh(dummy(:,:,jt), t_zt(:,:,jt), SLP(:,:,jt))
         END DO
      ELSE
         !! Dew-point is read:
         CALL GETVAR_1D_R8_3x3_to_1x1(cf_data, 'd2m',  dummy )
         DO jt = 1, Nt
            q_zt(:,:,jt) = q_air_dp( dummy(:,:,jt), SLP(:,:,jt) )
         END DO
      END IF



      CALL GETVAR_1D_R8_3x3_to_1x1(cf_data, 'ssrd',  rad_sw  )
      CALL GETVAR_1D_R8_3x3_to_1x1(cf_data, 'strd',  rad_lw  )

   END IF !IF( .NOT. l_3x3_ts )

   CALL TO_KELVIN_3D(SST, cname=TRIM('sst')   )
   CALL TO_KELVIN_3D(SIT, cname=TRIM('istl1') )


   WRITE(6,*) ''
   WRITE(6,*) ''
   !! All input time series read!


   ians=-1
   DO WHILE ( (ians<1).OR.(ians>nb_algos) )
      WRITE(6,*) 'Which algo to use?'
      WRITE(6,*) ' * "NEMO default (v4)"       => 1'
      WRITE(6,*) ' * "Andreas (2005)"          => 2'
      WRITE(6,*) ' * "Lupkes et al (2012)"     => 3'
      WRITE(6,*) ' * "Lupkes & Gryanik (2015)" => 4'

      READ(*,*) ians
      IF ( ians == 1 ) calgo = 'nemo'
      IF ( ians == 2 ) calgo = 'an05'
      IF ( ians == 3 ) calgo = 'lu12'
      IF ( ians == 4 ) calgo = 'lg15'
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
   zCdN = 0.

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

      IF (lverbose) THEN
         WRITE(6,*) ''; WRITE(6,*) ''
         WRITE(6,*) csep
      END IF
      WRITE(6,'(" #### Time = ",a," => isecday_utc = ",i6.6," => jt = ",i4.4)') cldate(jt), isecday_utc, jt
      IF (lverbose) THEN
         WRITE(6,*) csep
         WRITE(6,*) ''
         WRITE(6,*) '           ---- BEFORE BULK ALGO ('//TRIM(calgo)//') ----'
         WRITE(6,*) ''
      END IF

      info = DISP_DEBUG(lverbose, 'SIC', SIC(:,:,jt), '[fraction]')

      info = DISP_DEBUG(lverbose, 'Scalar wind speed at '//TRIM(czu), W10(:,:,jt), '[m/s]' )

      !! Air density at zt:
      rho_zt(:,:,jt) = rho_air( t_zt(:,:,jt), q_zt(:,:,jt), SLP(:,:,jt) )
      info = DISP_DEBUG(lverbose, 'Air density at zt', rho_zt(:,:,jt) , '[kg/m^3]')

      Cp_ma(:,:) = cp_air(q_zt(:,:,jt))
      info = DISP_DEBUG(lverbose, 'Cp of (moist) air at '//TRIM(czt), Cp_ma, '[J/K/kg]')

      rgamma(:,:) = gamma_moist(t_zt(:,:,jt), q_zt(:,:,jt))
      info = DISP_DEBUG(lverbose, 'Adiabatic lapse-rate of (moist) air', 1000.*rgamma, '[K/1000m]')

      info = DISP_DEBUG(lverbose, 'SIT', SIT(:,:,jt)-rt0, '[degC]')
      info = DISP_DEBUG(lverbose, 'SST', SST(:,:,jt)-rt0, '[degC]')

      SIQ = q_sat( SIT(:,:,jt), SLP(:,:,jt), l_ice=.TRUE. )
      info = DISP_DEBUG(lverbose, 'SIQ over ice = ', 1000.*SIQ, '[g/kg]')

      SSQ = rdct_qsat_salt*q_sat( SST(:,:,jt), SLP(:,:,jt), l_ice=.FALSE. )
      info = DISP_DEBUG(lverbose, 'SSQ over water = ', 1000.*SSQ, '[g/kg]')

      info = DISP_DEBUG(lverbose, 'Absolute   air temp. at '//TRIM(czt),     t_zt(:,:,jt) - rt0, '[deg.C]') ! Air temperatures at zt...
      info = DISP_DEBUG(lverbose, 'Specific   air humi. at '//TRIM(czt),     1000.*q_zt(:,:,jt), '[g/kg]') ! Air temperatures at zt...

      theta_zt(:,:,jt) = t_zt(:,:,jt) + rgamma(:,:)*zt ! potential temperature at zt
      info = DISP_DEBUG(lverbose, 'Potential  air temp. at '//TRIM(czt), theta_zt(:,:,jt) - rt0, '[deg.C]')

      tmp = virt_temp(theta_zt(:,:,jt), q_zt(:,:,jt))
      info = DISP_DEBUG(lverbose, 'Virt. pot. air temp. at '//TRIM(czt),          tmp     - rt0, '[deg.C]')
      info = DISP_DEBUG(lverbose, 'Pot. temp. diff. air-ice at '//TRIM(czt),    theta_zt(:,:,jt) - SIT(:,:,jt), '[deg.C]')
      info = DISP_DEBUG(lverbose, 'Virt. pot. temp. diff. at '//TRIM(czt),    tmp - virt_temp(SIT(:,:,jt), SIQ), '[deg.C]')

      !! We know enough to estimate the bulk Richardson number at zt:
      RiB_zt(:,:,jt) = Ri_bulk(zt, SIT(:,:,jt), theta_zt(:,:,jt), SIQ(:,:), q_zt(:,:,jt),  MAX(W10(:,:,jt),wspd_thrshld_ice) )
      info = DISP_DEBUG(lverbose, 'Bulk Richardson number at zt', RiB_zt(:,:,jt) , '[--]')

      !! Net solar heat flux into sea-ice:
      Qsw(:,:,jt) = (1._wp - rice_alb0)*rad_sw(:,:,jt)


      !! Only computing if there is sea-ice !!
      IF( SIC(1,1,jt) > 0.01_wp ) THEN

         SELECT CASE ( TRIM(calgo) )

         CASE ( 'nemo' )
            CALL turb_ice_nemo( zt, zu, SIT(:,:,jt), theta_zt(:,:,jt), SIQ(:,:), q_zt(:,:,jt), W10(:,:,jt),            &
               &                Cd_i(:,:,jt), Ch_i(:,:,jt), Ce_i(:,:,jt), theta_zu(:,:,jt), q_zu(:,:,jt), Ublk(:,:,jt),    &
               &                CdN=zCdN(:,:,jt), xz0=zz0(:,:,jt), xu_star=zus(:,:,jt), xL=zL(:,:,jt), xUN10=zUN10(:,:,jt) )

         CASE ( 'an05' )
            CALL turb_ice_an05( zt, zu, SIT(:,:,jt), theta_zt(:,:,jt), SIQ(:,:), q_zt(:,:,jt), W10(:,:,jt),             &
               &                Cd_i(:,:,jt), Ch_i(:,:,jt), Ce_i(:,:,jt), theta_zu(:,:,jt), q_zu(:,:,jt), Ublk(:,:,jt),     &
               &                CdN=zCdN(:,:,jt), xz0=zz0(:,:,jt), xu_star=zus(:,:,jt), xL=zL(:,:,jt), xUN10=zUN10(:,:,jt) )

         CASE ( 'lu12' )
            CALL turb_ice_lu12( zt, zu, SIT(:,:,jt), theta_zt(:,:,jt), SIQ(:,:), q_zt(:,:,jt), W10(:,:,jt), SIC(:,:,jt), &
               &                Cd_i(:,:,jt), Ch_i(:,:,jt), Ce_i(:,:,jt), theta_zu(:,:,jt), q_zu(:,:,jt), Ublk(:,:,jt),      &
               &                CdN=zCdN(:,:,jt), xz0=zz0(:,:,jt), xu_star=zus(:,:,jt), xL=zL(:,:,jt), xUN10=zUN10(:,:,jt) )

         CASE ( 'lg15' )
            CALL turb_ice_lg15( zt, zu, SIT(:,:,jt), theta_zt(:,:,jt), SIQ(:,:), q_zt(:,:,jt), W10(:,:,jt), SIC(:,:,jt), &
               &                Cd_i(:,:,jt), Ch_i(:,:,jt), Ce_i(:,:,jt), theta_zu(:,:,jt), q_zu(:,:,jt), Ublk(:,:,jt),      &
               &                CdN=zCdN(:,:,jt), xz0=zz0(:,:,jt), xu_star=zus(:,:,jt), xL=zL(:,:,jt), xUN10=zUN10(:,:,jt) )

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
         RiB_zu(:,:,jt) = Ri_bulk(zu, SIT(:,:,jt), theta_zu(:,:,jt), SIQ(:,:), q_zu(:,:,jt), Ublk(:,:,jt) )

         !! Turbulent heat fluxes:
         CALL BULK_FORMULA( zu, SIT(:,:,jt), SIQ(:,:), theta_zu(:,:,jt), q_zu(:,:,jt), &
            &              Cd_i(:,:,jt), Ch_i(:,:,jt), Ce_i(:,:,jt), W10(:,:,jt), Ublk(:,:,jt), SLP(:,:,jt), &
            &              TAU(:,:,jt), QH(:,:,jt), QL(:,:,jt),  &
            &              pEvap=SBLM(:,:,jt), prhoa=rho_zu(:,:,jt), l_ice=.TRUE. )

         !! Longwave radiative heat fluxes:
         Qlw(:,:,jt) = qlw_net( rad_lw(:,:,jt), SIT(:,:,jt), l_ice=.TRUE. )

         !! Non-solar heat flux:
         QNS(:,:,jt) = QH(:,:,jt) + QL(:,:,jt) + Qlw(:,:,jt)


         IF (lverbose) THEN
            WRITE(6,*) ''
            WRITE(6,*) '           ---- AFTER BULK ALGO ('//TRIM(calgo)//') ----'
            WRITE(6,*) ''
         END IF
         info = DISP_DEBUG(lverbose, 'density of air at '//TRIM(czu), rho_zu(:,:,jt),     '[kg/m^3]' )
         info = DISP_DEBUG(lverbose, 'theta_zu',                    theta_zu(:,:,jt)-rt0, '[deg.C]'  )
         info = DISP_DEBUG(lverbose, 'q_zu',                        q_zu(:,:,jt)*1000.,    '[g/kg]'  )
         info = DISP_DEBUG(lverbose, 'Ublk',                        Ublk(:,:,jt),          '[m/s]'   )
         info = DISP_DEBUG(lverbose, 'RiB_zu',                      RiB_zu(:,:,jt),       '[]'       )

         !! Sanity check:
         IF( ABS(RiB_zu(1,1,jt)) > 100._wp ) THEN
            PRINT *, 'Fucked up Richardson number!!!, is everything okay???', RiB_zu(1,1,jt)
            STOP
         END IF
         IF( (ABS(rho_zu(1,1,jt)) > 1.7_wp).OR.(ABS(rho_zu(1,1,jt)) < 0.9_wp) ) THEN
            PRINT *, 'We fucked up, density at zu is irrealistic!!!'
            STOP
         END IF

         IF (lverbose) THEN
            WRITE(6,*) ''
            WRITE(6,*) csep
         END IF

      ELSE

         WRITE(6,*) ''
         WRITE(6,*) ' *** DOING NOTHING! NO SEA-ICE present for this time step !!!'
         WRITE(6,*) ''

      END IF

   END DO !DO jt = 1, Nt






   CALL PT_SERIES(vtime(:), REAL(rho_zu(1,1,:),4), 'aerobulk_test_ice_series_'//TRIM(cn_exp)//'_'//TRIM(calgo)//'.nc', 'time', &
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
      &           vdt10=REAL( 1000.*Cd_i(1,1,:),4), cv_dt10='Cd_i', cun10='',  cln10='Drag coefficient',       &
      &           vdt11=REAL( 1000.*Ch_i(1,1,:),4), cv_dt11='Ch_i', cun11='',  cln11='Sens. Heat coeff.',      &
      &           vdt12=REAL(  zz0(1,1,:),4), cv_dt12='z0',     cun12='m',     cln12='Roughness length',       &
      &           vdt13=REAL(  RiB_zt(1,1,:),4), cv_dt13='Rib_zt', cun13='',   cln13='Bulk Richardson number at zt', &
      &           vdt14=REAL(  RiB_zu(1,1,:),4), cv_dt14='Rib_zu', cun14='m',  cln14='Bulk Richardson number at zu', &
      &           vdt15=REAL(SIT(1,1,:)-rt0,4), cv_dt15='SIT',  cun15='deg.C', cln15='Sea-ice temperature',    &
      &           vdt16=REAL(t_zt(1,1,:)-SIT(1,1,:),4), cv_dt16='t2m-SIT',     cun16='deg.C', cln16='2m air-sea temperature difference', &
      &           vdt17=REAL(1000.*zCdN(1,1,:),4), cv_dt17='CdN', cun17='',    cln17='Neutral-stability drag coefficient',  &
      &           vdt18=REAL(SIC(1,1,:),4),        cv_dt18='A',   cun18='',    cln18='Sea-ice concentration'  &
      &           )

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
      PRINT *,' -3         => fields in input netcdf are 3x3 in space (those of NEMO/tests/STATION_ASF!)'
      PRINT *,''
      PRINT *,' -r         => humidity in NetCDF file is Relative Humidity [%]'
      PRINT *,''
      PRINT *,' -h         => Show this message'
      PRINT *,''
      STOP
      !!
   END SUBROUTINE usage_test

END PROGRAM TEST_AEROBULK_BUOY_SERIES_ICE
