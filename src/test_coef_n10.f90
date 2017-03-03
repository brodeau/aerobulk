! AeroBulk / 2015 / L. Brodeau (brodeau@gmail.com)
! https://sourceforge.net/p/aerobulk

PROGRAM TEST_COEF_N10

   USE mod_const
   USE mod_thermo

   USE mod_blk_neutral_10m

   IMPLICIT NONE

   INTEGER, PARAMETER :: nb_algos = 4

   CHARACTER(len=10), DIMENSION(nb_algos), PARAMETER :: vca = (/ 'coare30' , 'coare35', 'ncar', 'ecmwf' /)

   REAL(4), DIMENSION(nb_algos) :: vCdn10, vCen10, vChn10

   REAL(wp), PARAMETER ::   &
      &   wind_max = 60.,  &
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
   REAL(wp),    DIMENSION(lx,ly) :: zz0, zus

   REAL(wp), DIMENSION(lx,ly) :: SLP, &
      &  U_N10

   REAL(wp) :: dw0, dw

   REAL(wp), DIMENSION(lx,ly) :: Cdn10, Cen10, Chn10


   LOGICAL :: l_ask_for_slp = .FALSE. , &  !: ask for SLP, otherwize assume SLP = 1010 hPa
      &     l_use_rh      = .FALSE.      !: ask for RH rather than q for humidity

   INTEGER, PARAMETER :: &
      &   n_w    = 1201

   REAL, DIMENSION(n_w,1) :: &
      &   t_w10, t_cdn10, t_cen10, t_chn10


   jpi = lx ; jpj = ly


   nb_itt = 50  ! 50 itterations in bulk algorithm...

   jarg = 0

   DO WHILE ( jarg < iargc() )

      jarg = jarg + 1
      CALL getarg(jarg,car)

      SELECT CASE (trim(car))

      CASE('-h')
         call usage_test()

      CASE('-p')
         !PRINT *, 'Will ask for SLP'
         l_ask_for_slp = .TRUE.

      CASE('-r')
         !PRINT *, "Using RH instead of q!!!"
         l_use_rh = .TRUE.

      CASE DEFAULT
         PRINT *, 'Unknown option: ', trim(car) ; PRINT *, ''
         CALL usage_test()

      END SELECT

   END DO

   PRINT *, ''


   SLP = Patm



   !!
   !! Table of wind speeds from 0 to wind_max  m/s :
   dw0 = wind_max/(n_w - 1)
   !FORALL (jw = 1:n_ew)  t_w10(jw) = (jw-1)*dw
   !PRINT *, 'dw0 =', dw0 ; STOP
   !!
   t_w10(1,1) = 0.0
   DO jw = 2, n_w
      dw = 0.25*dw0
      IF ( t_w10(jw-1,1) >= 5.  ) dw =    dw0
      IF ( t_w10(jw-1,1) >= 30. ) dw = 2.*dw0
      t_w10(jw,1) = t_w10(jw-1,1) + dw
   END DO
















   PRINT *, 'Give neutral wind speed at 10m (m/s):'
   READ(*,*) U_N10
   PRINT *, ''


   DO ialgo = 1, nb_algos

      calgob = trim(vca(ialgo))

      zz0 = 0.
      zus = 0.

      CALL TURB_NEUTRAL_10M(calgob, U_N10, CdN10, ChN10, Cen10) !,    xz0=zz0, xu_star=zus)

      vCdn10(ialgo) = REAL(1000.*Cdn10(1,1) ,4)
      vChn10(ialgo) = REAL(1000.*Chn10(1,1) ,4)
      vCen10(ialgo) = REAL(1000.*Cen10(1,1) ,4)


   END DO


   PRINT *, ''; PRINT *, ''

   PRINT *, ''
   PRINT *, '   *** Bulk Transfer Coefficients:'
   PRINT *, '========================================================================='
   PRINT *, '  Algorithm:           ',TRIM(vca(1)) ,'    |    ',TRIM(vca(2)),'       |    ',TRIM(vca(3)),'       |    ',TRIM(vca(4))
   PRINT *, '========================================================================='
   PRINT *, '      C_D_N10   =   ', vCdn10        , '[10^-3]'
   PRINT *, '      C_E_N10   =   ', vCen10        , '[10^-3]'
   PRINT *, '      C_H_N10   =   ', vChn10        , '[10^-3]'
   PRINT *, ''
   PRINT *, ''



   !! For all wind speeds:


   jpi = n_w ; jpj = 1


   DO ialgo = 1, nb_algos

      calgob = trim(vca(ialgo))

      CALL TURB_NEUTRAL_10M(calgob, t_w10, t_CdN10, t_ChN10, t_Cen10)




      cf_out = 'dat/Neutral_coeff_U10N_'//TRIM(calgob)//'.dat'

      OPEN(unit = 11, file = trim(cf_out), status = 'unknown')
      WRITE(11,'("#     Wind (U_N10)       Cd_N10             Ce_N10")')
      DO jw = 1, n_w
         WRITE(11,'(1f16.8, " " ,1f16.8, " " ,1f16.8)') t_w10(jw,1), 1000.*t_cdn10(jw,1), 1000.*t_cen10(jw,1)
      ENDDO
      CLOSE(11)


      PRINT *, ''
      PRINT *, '  *** Just created file ', trim(cf_out)



   END DO





END PROGRAM TEST_COEF_N10











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
