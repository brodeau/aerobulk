! AeroBulk / 2015 / L. Brodeau

PROGRAM TEST_COEF_N10

   USE mod_const
   USE mod_phymbl

   USE mod_blk_neutral_10m

   IMPLICIT NONE

   INTEGER, PARAMETER :: nb_algos = 4

   CHARACTER(len=10), DIMENSION(nb_algos), PARAMETER :: vca = (/ 'coare3p0' , 'coare3p6', 'ecmwf   ', 'ncar    ' /)

   REAL(wp), PARAMETER ::   &
      &   wind_max = 60.,   &
      &   zu  = 10.

   INTEGER, PARAMETER ::   &
      &   n_dt   = 21,   &
      &   n_u    = 30,   &
      &   n_q    = 10

   CHARACTER(len=2) :: car

   CHARACTER(len=256) :: cf_out

   CHARACTER(len=100) :: calgob

   INTEGER :: jarg, ialgo, jw

   INTEGER, PARAMETER :: lx=1, ly=1

   REAL(wp), DIMENSION(lx,ly) :: SLP

   REAL(wp) :: dw0, dw

   LOGICAL :: l_ask_for_slp = .FALSE. , &  !: ask for SLP, otherwize assume SLP = 1010 hPa
      &     l_use_rh      = .FALSE.      !: ask for RH rather than q for humidity

   INTEGER, PARAMETER :: &
      &   n_w    = 1201

   REAL, DIMENSION(n_w,1) :: &
      &   t_w10, t_cdn10, t_cen10, t_chn10, t_z0

   OPEN(6, FORM='formatted', RECL=512)

   nb_iter = 50  ! 50 itterations in bulk algorithm...

   jarg = 0

   DO WHILE ( jarg < iargc() )

      jarg = jarg + 1
      CALL getarg(jarg,car)

      SELECT CASE (trim(car))

      CASE('-h')
         call usage_test()

      CASE('-p')
         !WRITE(6,*) 'Will ask for SLP'
         l_ask_for_slp = .TRUE.

      CASE('-r')
         !PRINT *, "Using RH instead of q!!!"
         l_use_rh = .TRUE.

      CASE DEFAULT
         WRITE(6,*) 'Unknown option: ', trim(car) ; WRITE(6,*) ''
         CALL usage_test()

      END SELECT

   END DO

   WRITE(6,*) ''


   SLP = Patm



   !!
   !! Table of wind speeds from 0 to wind_max  m/s :
   dw0 = wind_max/REAL(n_w-1, wp)
   !FORALL (jw = 1:n_ew)  t_w10(jw) = (jw-1)*dw
   !WRITE(6,*) 'dw0 =', dw0 ; STOP
   !!
   t_w10(1,1) = 0.0
   DO jw = 2, n_w
      dw = 0.25*dw0
      IF ( t_w10(jw-1,1) >= 5.  ) dw =    dw0
      IF ( t_w10(jw-1,1) >= 30. ) dw = 2.*dw0
      t_w10(jw,1) = t_w10(jw-1,1) + dw
   END DO


   !! For all wind speeds:

   DO ialgo = 1, nb_algos

      calgob = TRIM(vca(ialgo))

      CALL TURB_NEUTRAL_10M(calgob, t_w10, t_CdN10, t_ChN10, t_Cen10, t_z0)

      cf_out = 'dat/Neutral_coeff_U10N_'//TRIM(calgob)//'.dat'

      OPEN(unit = 11, file = TRIM(cf_out), status = 'unknown')
      !                 0.00000000       1.57248293       1.89781282       1.80325771       0.00041612
      WRITE(11,'("#   N10 Wind (m/s)       Cd_N10           Ce_N10           Ch_N10           z0 (mm)")')
      DO jw = 1, n_w
         WRITE(11,'(1f16.8, " " ,1f16.8, " " ,1f16.8, " " ,1f16.8, " " ,1f16.8)') &
            &   t_w10(jw,1), 1000.*t_cdn10(jw,1), 1000.*t_cen10(jw,1), 1000.*t_chn10(jw,1), 1000.*t_z0(jw,1)
      ENDDO
      CLOSE(11)


      WRITE(6,*) '  *** Just created file ', trim(cf_out)
      WRITE(6,*) ''

   END DO

   WRITE(6,*) ''
   CLOSE(6)


CONTAINS

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

END PROGRAM TEST_COEF_N10
