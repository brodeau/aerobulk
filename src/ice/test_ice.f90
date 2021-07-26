! AeroBulk / 2020 / L. Brodeau

PROGRAM TEST_ICE

   USE mod_const
   USE mod_phymbl
   USE mod_blk_ice_an05

   IMPLICIT NONE

   REAL :: dus, tair

   INTEGER, PARAMETER :: nustar = 101

   REAL(wp), DIMENSION(1,nustar)   :: vus, vz0
   REAL(wp), DIMENSION(1,nustar,2) :: vz0tq
   REAL(wp), DIMENSION(1,nustar)   :: vtair, vnua
   
   REAL(wp), PARAMETER :: &
      &   us_min =  0.001 , &
      &   us_max =  0.701

   INTEGER :: ju

   CHARACTER(len=128) :: cfout1 !, cfout2

   dus = (us_max - us_min)/(nustar - 1)

   PRINT *, ' Temperature increment = ', dus

   vus(1,:) = (/ ( us_min + (ju-1)*dus , ju=1,nustar ) /)

   PRINT *, vus
   PRINT *, ''


   vtair(:,:) = 270.0
   PRINT *, 'Give air temperature (deg.C):'
   READ(*,*) tair
   tair = tair + rt0
   vtair(:,:) = tair

   !! Computing kinetic viscosity:
   vnua(:,:) = visc_air( vtair(:,:) )

   !PRINT *, vnua

   !! Computing z0 as a function of u*:
   vz0(:,:) = rough_leng_m( vus(:,:) , vnua(:,:) )


   PRINT *, ' vz0 =', vz0
   

   vz0tq(:,:,:) = rough_leng_tq( vz0(:,:), vus(:,:) , vnua(:,:) )

   PRINT *, ' vz0t =', vz0tq(:,:,1)
   PRINT *, ' vz0q =', vz0tq(:,:,2)

   !! Advanced formula:
   !vqsat1(:,:) = rdct_qsat_salt*q_sat(vus_k, vtair)

   !! Simple with constant density
   !vrho  (:,:) = 1.2
   !vqsat2(:,:) = rdct_qsat_salt*q_sat_simple(vus_k, vrho)


   !vq10(:,:) = q_air_rh(vrh, vt10, vtair)   ! Specific humidity at 10m
   !vrho(:,:) = rho_air(vt10, vq10, vtair)    ! air density at 10m
   
   !vqsat4(:,:) = rdct_qsat_salt*q_sat_simple(vus_k, vrho)
   
   ! Buck Formula !
   !vqsat5(:,:) = rdct_qsat_salt*q_sat(vus_k, vtair, cform='buck')
   


   
   cfout1 = 'z0_z0t_z0q__ustar_test.dat'
   
   OPEN(11, FILE=cfout1, FORM='formatted', STATUS='unknown',RECL=512)
   WRITE(11,*) '#    u*         z0           z0t          z0q'
   DO ju = 1, nustar
      WRITE(11,*) REAL(vus(1,ju),4), REAL(vz0(1,ju),4), REAL(vz0tq(1,ju,1),4), REAL(vz0tq(1,ju,2),4)
   END DO
   CLOSE(11)
   
   PRINT *, ''
   PRINT *, TRIM(cfout1), ' written!'
   PRINT *, ''
   
END PROGRAM TEST_ICE
