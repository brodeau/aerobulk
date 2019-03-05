! AeroBulk / 2015 / L. Brodeau

PROGRAM EXAMPLE_CALL_AEROBULK_COMPUTE
   
   USE mod_aerobulk
   USE mod_const
   
   IMPLICIT NONE

   INTEGER, PARAMETER :: nx = 1, &
      &                  ny = 1
   
   REAL(wp), DIMENSION(nx,ny) :: zsst, zt_zt, zq_zt, zU_zu, zV_zu, zslp, &
      &                         zRsw, zRlw,           &
      &                         zQL, zQH, zTau_x, zTau_y, zTs

   PRINT *, ''

   zsst = 22.
   PRINT *, ' *** SST = ', zsst, ' deg.C'
   zsst  = zsst + rt0  ! [K]
   
   zt_zt = rt0 + 20.  ! [K]
   zq_zt = 0.012         ! [kg/kg]
   zU_zu = 4.            ! [m/s]
   zV_zu = 10.           ! [m/s]
   zslp  = 101000.       ! [Pa]
   
   zRsw  = 0. ! (night)  ! [W/m^2]
   zRlw  = 350.          ! [W/m^2]

   PRINT *, ''
   PRINT *, ' *********** COARE 3.0 *****************'   
   CALL aerobulk_model( 'coare', 2._wp, 10._wp, zsst, zt_zt, &
      &                 zq_zt, zU_zu, zV_zu, zslp,   &
      &                 zQL, zQH, zTau_x, zTau_y,    &
      &                 Niter=20, rad_sw=zRsw, rad_lw=zRlw, T_s=zTs )
   PRINT *, ''
   PRINT *, ' Wind speed             =', SQRT(zU_zu*zU_zu + zV_zu*zV_zu), ' m/s' ; PRINT *, ''
   PRINT *, ' Sensible heat flux: QH =', zQH ; PRINT *, ''
   PRINT *, '  Latent  heat flux: QL =', zQL ; PRINT *, ''
   PRINT *, ' Skin temperature: SSST =', zTs-rt0 ; PRINT *, ''
   PRINT *, ' Tau_x =', zTau_x ; PRINT *, ''
   PRINT *, ' Tau_y =', zTau_y ; PRINT *, ''
   PRINT *, ' Tau   =', SQRT(zTau_x*zTau_x + zTau_y*zTau_y) ; PRINT *, ''
   PRINT *, ''

   PRINT *, ''
   PRINT *, ' *********** ECMWF *****************'   
   CALL aerobulk_model( 'ecmwf', 2._wp, 10._wp, zsst, zt_zt, &
      &                 zq_zt, zU_zu, zV_zu, zslp,   &
      &                 zQL, zQH, zTau_x, zTau_y,    &
      &                 Niter=20, rad_sw=zRsw, rad_lw=zRlw, T_s=zTs )
   PRINT *, ''
   PRINT *, ' Wind speed             =', SQRT(zU_zu*zU_zu + zV_zu*zV_zu), ' m/s' ; PRINT *, ''
   PRINT *, ' Sensible heat flux: QH =', zQH ; PRINT *, ''
   PRINT *, '  Latent  heat flux: QL =', zQL ; PRINT *, ''
   PRINT *, ' Skin temperature: SSST =', zTs-rt0 ; PRINT *, ''
   PRINT *, ' Tau_x =', zTau_x ; PRINT *, ''
   PRINT *, ' Tau_y =', zTau_y ; PRINT *, ''
   PRINT *, ' Tau   =', SQRT(zTau_x*zTau_x + zTau_y*zTau_y) ; PRINT *, ''
   PRINT *, ''
   
   PRINT *, ''
   PRINT *, ' *********** NCAR *****************'   
   CALL aerobulk_model( 'ncar', 2._wp, 10._wp, zsst, zt_zt, &
      &                 zq_zt, zU_zu, zV_zu, zslp,   &
      &                 zQL, zQH, zTau_x, zTau_y,    &
      &                 Niter=20 )
   PRINT *, ''
   PRINT *, ' Wind speed             =', SQRT(zU_zu*zU_zu + zV_zu*zV_zu), ' m/s' ; PRINT *, ''
   PRINT *, ' Sensible heat flux: QH =', zQH ; PRINT *, ''
   PRINT *, '  Latent  heat flux: QL =', zQL ; PRINT *, ''
   PRINT *, ' Tau_x =', zTau_x ; PRINT *, ''
   PRINT *, ' Tau_y =', zTau_y ; PRINT *, ''
   PRINT *, ' Tau   =', SQRT(zTau_x*zTau_x + zTau_y*zTau_y) ; PRINT *, ''
   PRINT *, ''
   
   
END PROGRAM EXAMPLE_CALL_AEROBULK_COMPUTE
