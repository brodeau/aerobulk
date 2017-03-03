! AeroBulk / 2015 / L. Brodeau (brodeau@gmail.com)
! https://sourceforge.net/p/aerobulk

PROGRAM TEST_THERMO

   USE mod_const
   USE mod_thermo

   IMPLICIT NONE

   REAL :: dtime, slp, dT

   INTEGER, PARAMETER :: ntemp = 101

   CHARACTER(len=6) :: cslp
   CHARACTER(len=2) :: cdT

   REAL, DIMENSION(1,ntemp) :: vsst, vmsl, vqsat1, vqsat2, vqsat3, vqsat4, vqsat5, &
      &                        vsst_k, vrho, vt10, vq10, vrh

   REAL, PARAMETER :: &
      &   t_min =  -5. , &
      &   t_max =  35.

   INTEGER :: jt

   CHARACTER(len=128) :: cfout1 !, cfout2

   jpi = 1
   jpj = ntemp

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
   vqsat1(:,:) = 0.98*q_sat(vsst_k, vmsl)

   !! Simple with constant density
   vrho  (:,:) = 1.2
   vqsat2(:,:) = 0.98*q_sat_simple(vsst_k, vrho)

   !! Simple with better density (density at SST! 0m, not 10m !)
   vrho  (:,:) = rho_air(vsst_k, vqsat1, vmsl) ! we use the good q_sat to get a good rho at z=0
   vqsat3(:,:) = 0.98*q_sat_simple(vsst_k, vrho)

   !! Same but using density of air at 10m => T=SST-2 and not saturated => 80% hum!
   vrh(:,:) = 0.8                          ! Relative humidity at 10m
   vt10(:,:) = vsst_k + dT              ! Air temperature at 10m

   !! Don't need to do what follows cause we specified that SLP is at 10m !
   !vrho(:,:) = 1.2 ! first guess...
   !DO jt = 1, 5
   !   vq10(:,:) = q_air_rh(vrh, vt10, vmsl-vrho*grav*10.)   ! Specific humidity at 10m
   !   vrho(:,:) = rho_air(vt10, vq10, vmsl-vrho*grav*10.)    ! air density at 10m
   !END DO

   vq10(:,:) = q_air_rh(vrh, vt10, vmsl)   ! Specific humidity at 10m
   vrho(:,:) = rho_air(vt10, vq10, vmsl)    ! air density at 10m
   
   vqsat4(:,:) = 0.98*q_sat_simple(vsst_k, vrho)
   
   ! Buck Formula !
   vqsat5(:,:) = 0.98*q_sat(vsst_k, vmsl, cform='buck')
   


   
   cfout1 = 'qsat_test_'//cslp//'Pa_dt_'//cdT//'.dat'
   
   OPEN(11, FILE=cfout1, FORM='formatted', STATUS='unknown',RECL=512)
   WRITE(11,*) '**  SST   Method 1 (advanced)  Method 2 (simple+1.2)&
      & Method 3 (simple+rho_0m)   Method 4 (simple+rho_10m)  Bulk fomula'
   DO jt = 1, ntemp
      WRITE(11,*) REAL(vsst(1,jt),4), REAL(vqsat1(1,jt),4), REAL(vqsat2(1,jt),4), REAL(vqsat3(1,jt),4), &
         &        REAL(vqsat4(1,jt),4), REAL(vqsat5(1,jt),4)
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
   


END PROGRAM TEST_THERMO
