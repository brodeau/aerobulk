! AeroBulk / 2019 / L. Brodeau

!! https://www.pmel.noaa.gov/ocs/flux-documentation


PROGRAM TEST_AEROBULK_CDNF_SERIES

   USE mod_const
   USE mod_phymbl

   USE io_ezcdf     !* routines for netcdf input/output (par of SOSIE package)

   USE mod_cdn_form_ice

   IMPLICIT NONE

   REAL(wp), PARAMETER :: rz0_i_s_0  = 0.69e-3_wp  ! default z0 for skin drag over ice [m]   Eq. 43
   REAL(wp), PARAMETER :: rz0_i_f_0  = 4.54e-4_wp  ! default z0 for form drag over ice [m] bottom p.562 MIZ
   REAL(wp), PARAMETER :: rz0_w_s_0 = 3.27E-4      ! fixed roughness length over water (paragraph below Eq.36)



   LOGICAL, PARAMETER :: lverbose = .TRUE.
   LOGICAL, PARAMETER :: ldebug   = .FALSE.
   INTEGER, PARAMETER :: jtdbg    = 1   ! if (ldebug) that's the time step we'll start from...

   INTEGER, PARAMETER :: nb_algos = 5

   CHARACTER(len=800) :: cf_data='0', cunit_t, clnm_t, clndr_t

   CHARACTER(len=80) :: csep='#################################################################################'

   CHARACTER(len=2) :: car

   INTEGER :: jt0, jt, jarg, ialgo, n0, info

   INTEGER :: nx, ny, Nt, ians

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
   !REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: SIT, SST, SKT, SLP, W10, SIC, t_zt, theta_zt, q_zt, rad_sw, rad_lw, dummy
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: dummy, W10, SIC

   !REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zUN10

   !REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: t_zu, theta_zu, q_zu, rho_zt, rho_zu

   REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: z0_ice, z0_oce, tmp

   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: CdNFi !, Ce_i, Ch_i, QH, QL, Qsw, QNS, Qlw, RiB_zt, RiB_zu, TAU, SBLM

   REAL(wp) :: zu

   CHARACTER(len=3) :: czu



   TYPE(t_unit_t0) :: tut_time_unit
   TYPE(date)      :: d_idate

   nb_iter = 8 ! 8 itterations in bulk algorithm...

   !OPEN(6, FORM='formatted', RECL=512)

   !----------------------------------------------

   jarg = 0

   DO WHILE ( jarg < command_argument_COUNT() )

      jarg = jarg + 1
      CALL get_command_ARGUMENT(jarg,car)

      SELECT CASE (TRIM(car))

      CASE('-h')
         call usage_test()

      CASE('-f')
         jarg = jarg + 1
         CALL get_command_ARGUMENT(jarg,cf_data)

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

   WRITE(6,*) ''
   WRITE(6,*) ' *** Allocating arrays according to nx,ny,Nt =', nx,ny,Nt


   ALLOCATE ( ctime(Nt), cdate(Nt), clock(Nt), chh(Nt), cmn(Nt), cldate(Nt), idate(Nt), vtime(Nt), vlon(nx) )
   ALLOCATE ( z0_ice(nx,ny), z0_oce(nx,ny), tmp(nx,ny) )
   !ALLOCATE ( Ublk(nx,ny,Nt), zUN10(nx,ny,Nt) )

   !ALLOCATE ( SIT(nx,ny,Nt), SST(nx,ny,Nt), SKT(nx,ny,Nt), SLP(nx,ny,Nt), , t_zt(nx,ny,Nt), theta_zt(nx,ny,Nt), q_zt(nx,ny,Nt),  &
   !   &       rad_sw(nx,ny,Nt), rad_lw(nx,ny,Nt), SIC(nx,ny,Nt) )
   !ALLOCATE ( t_zu(nx,ny,Nt), theta_zu(nx,ny,Nt), q_zu(nx,ny,Nt), rho_zt(nx,ny,Nt), rho_zu(nx,ny,Nt), dummy(nx,ny,Nt) )
   !ALLOCATE ( CdNFi(nx,ny,Nt), Ce_i(nx,ny,Nt), Ch_i(nx,ny,Nt), QH(nx,ny,Nt), QL(nx,ny,Nt), Qsw(nx,ny,Nt), Qlw(nx,ny,Nt), QNS(nx,ny,Nt), &
   !   &       RiB_zt(nx,ny,Nt), RiB_zu(nx,ny,Nt), TAU(nx,ny,Nt), SBLM(nx,ny,Nt) )

   ALLOCATE ( SIC(nx,ny,Nt), dummy(nx,ny,Nt), W10(nx,ny,Nt), CdNFi(nx,ny,Nt) )

   WRITE(6,*) ' *** Allocation completed!'
   WRITE(6,*) ''


   CALL GETVAR_1D(cf_data, 'time',  vtime ) ; ! (hours since ...)
   CALL GET_VAR_INFO(cf_data, 'time', cunit_t, clnm_t,  clndr=clndr_t)
   PRINT *, 'time unit = "'//TRIM(cunit_t)//'"'
   tut_time_unit = GET_TIME_UNIT_T0( TRIM(cunit_t) ) ; ! origin
   PRINT *, ' *** Digested time unit is: ', tut_time_unit


   PRINT *, 'LOLO1'
   CALL GETVAR_1D(cf_data, 'siconc',   SIC  )  ; ! sea-ice concentration (0-1)

   PRINT *, 'LOLO2'

   !CALL GETVAR_1D(cf_data, 'istl1',    SIT  )  ; ! istl1 is in K !

   !CALL GETVAR_1D(cf_data, 'sst',    SST  )  ; ! sst is in K !

   !CALL GETVAR_1D(cf_data, 'skt',    SKT  )  ; ! skt is in K !

   !CALL GETVAR_1D(cf_data, 'msl',    SLP  )

   CALL GETVAR_1D(cf_data, 'u10', W10  )
   CALL GETVAR_1D(cf_data, 'v10', dummy  )
   W10 = SQRT ( W10*W10 + dummy*dummy )

   !CALL GETVAR_1D(cf_data, 't2m',  t_zt )

   !CALL GETVAR_1D(cf_data, 'd2m',  dummy )
   !DO jt = 1, Nt
   !   q_zt(:,:,jt) = q_air_dp( dummy(:,:,jt), SLP(:,:,jt) )
   !END DO

   !CALL GETVAR_1D(cf_data, 'radsw',  rad_sw  )
   !CALL GETVAR_1D(cf_data, 'radlw',  rad_lw  )


   WRITE(6,*) ''
   WRITE(6,*) ''
   !! All input time series read!


   ians=-1
   DO WHILE ( (ians<1).OR.(ians>nb_algos) )
      WRITE(6,*) 'Which algo to use?'
      WRITE(6,*) ' * "LU12"  => 1'
      WRITE(6,*) ' * "LU12l" => 2'
      WRITE(6,*) ' * "LU13"  => 3'
      WRITE(6,*) ' * "LG15"  => 4'
      WRITE(6,*) ' * "LG15l" => 5'

      READ(*,*) ians
      IF ( ians == 1 ) calgo = 'LU12'
      IF ( ians == 2 ) calgo = 'LU12l'
      IF ( ians == 3 ) calgo = 'LU13'
      IF ( ians == 4 ) calgo = 'LG15'
      IF ( ians == 5 ) calgo = 'LG15l'
   END DO
   WRITE(6,*) '  ==> your choice: ', TRIM(calgo)
   WRITE(6,*) ''



   !! zu
   !! ~~
   WRITE(6,*) 'Give "zu", height of wind speed measurement in meters (generally 10m):'
   READ(*,*) zu
   WRITE(6,*) ''
   WRITE(czu,'(i2,"m")') INT(zu)



   !! Some initializations:

   ialgo = 0

   isecday_utc = 0

   d_idate%year   = 0
   d_idate%month  = 0
   d_idate%day    = 0
   d_idate%hour   = 0
   d_idate%minute = 0
   d_idate%second = 0


   z0_ice(:,:) = rz0_i_s_0
   z0_oce(:,:) = rz0_w_s_0


   jt0 = 1
   IF( ldebug .AND. (jt>0) ) jt0 = jtdbg

   !! Time loop:
   DO jt = jt0, Nt

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
      END IF


      info = DISP_DEBUG(lverbose, 'SIC', SIC(:,:,jt), '[fraction]')

      info = DISP_DEBUG(lverbose, 'Scalar wind speed at '//TRIM(czu), W10(:,:,jt), '[m/s]' )


      SELECT CASE ( TRIM(calgo) )

      CASE ( 'LU12' )
         CdNFi(:,:,jt) = CdN10_f_LU12(          SIC(:,:,jt), z0_oce(:,:) ) !!!,  pSc, phf, pDi  )

      CASE ( 'LU12l' )
         CdNFi(:,:,jt) = CdN_f_LU12_eq36(  zu, SIC(:,:,jt)       )

      CASE ( 'LU13' )
         CdNFi(:,:,jt) = CdN10_f_LU13(          SIC(:,:,jt)       )

      CASE ( 'LG15' )
         CdNFi(:,:,jt) = CdN_f_LG15(       zu, SIC(:,:,jt), z0_ice(:,:) ) !!! ,  pSc, phf, pDi  )

      CASE ( 'LG15l' )
         CdNFi(:,:,jt) = CdN_f_LG15_light( zu, SIC(:,:,jt), z0_oce(:,:) )

      CASE DEFAULT
         PRINT *, 'UNKNOWN algo: '//TRIM(calgo)//' !!!'
         STOP
      END SELECT

   END DO !DO jt = 1, Nt




   CALL PT_SERIES(vtime(:), REAL(SIC(1,1,:),4), 'CDN'//TRIM(czu)//'_form_'//TRIM(calgo)//'.nc', 'time', &
      &           'A', '[0-1]', 'Sea-ice concentration', -9999._4, &
      &           ct_unit=TRIM(cunit_t), ct_clnd=TRIM(clndr_t), &
      &           vdt02=REAL(  W10(1,1,:),4), cv_dt02='Wind',   cun02='m/s',   cln02='Module of Wind Speed',   &
      &           vdt03=REAL( 1000.*CdNFi(1,1,:),4), cv_dt03='CdN_f_i', cun03='',  cln03='Drag coefficient due to form / sea-ice'  &
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
      PRINT *,''
      PRINT *,' -f <netcdf_file>  => file containing data'
      PRINT *,''
      PRINT *,' -h   => Show this message'
      PRINT *,''
      STOP
      !!
   END SUBROUTINE usage_test

END PROGRAM TEST_AEROBULK_CDNF_SERIES
