! AeroBulk / 2015 / L. Brodeau

PROGRAM EXAMPLE_CALL_AEROBULK_COMPUTE

   USE mod_aerobulk
   USE mod_const

   IMPLICIT NONE

   INTEGER, PARAMETER :: nx = 2, &
      &                  ny = 1

   REAL(wp), DIMENSION(nx,ny) :: zsst, zt_zt, zq_zt, zU_zu, zV_zu, zslp, &
      &                         zRsw, zRlw,           &
      &                         zQL, zQH, zTau_x, zTau_y, zTs

   PRINT *, ''

   zsst = 22.
   PRINT *, ' *** SST = ', zsst, ' deg.C'
   zsst  = zsst + rt0  ! [K]

   zt_zt(1,1) = rt0 + 20.  ! [K]
   zt_zt(2,1) = rt0 + 25.  ! [K] ! second case is stable ABL as t_air > SST (25>22)!


   zq_zt = 0.012         ! [kg/kg]
   zU_zu = 4.            ! [m/s]
   zV_zu = 10.           ! [m/s]
   zslp  = 101000.       ! [Pa]

   zRsw  = 0. ! (night)  ! [W/m^2]
   zRlw  = 350.          ! [W/m^2]


   PRINT *, ''
   PRINT *, ' *********** COARE 3.6 *****************'
   CALL aerobulk_model( 'coare3p6', 2._wp, 10._wp, zsst, zt_zt, &
      &                 zq_zt, zU_zu, zV_zu, zslp,   &
      &                 zQL, zQH, zTau_x, zTau_y,    &
      &                 Niter=20, rad_sw=zRsw, rad_lw=zRlw, T_s=zTs )
   PRINT *, ''
   PRINT *, ' Wind speed at zu       =', REAL(SQRT(zU_zu*zU_zu + zV_zu*zV_zu),4), ' m/s'
   PRINT *, '    SST                 =', REAL(zsst,4), ' K'
   PRINT *, ' Air temperature at zt  =', REAL(zt_zt,4), ' K'  ; PRINT *, ''
   PRINT *, ' Sensible heat flux: QH =', REAL(zQH,4), ' W/m**2'
   PRINT *, '  Latent  heat flux: QL =', REAL(zQL,4), ' W/m**2'
   PRINT *, ' Skin temperature: SSST =', REAL(zTs-rt0,4), ' K'
   PRINT *, ' Tau_x =', REAL(zTau_x,4), ' N/m**2'
   PRINT *, ' Tau_y =', REAL(zTau_y,4), ' N/m**2'
   PRINT *, ' Tau   =', REAL(SQRT(zTau_x*zTau_x + zTau_y*zTau_y),4), ' N/m**2' ; PRINT *, ''
   PRINT *, ''

   PRINT *, ''
   PRINT *, ' *********** ECMWF *****************'
   CALL aerobulk_model( 'ecmwf', 2._wp, 10._wp, zsst, zt_zt, &
      &                 zq_zt, zU_zu, zV_zu, zslp,   &
      &                 zQL, zQH, zTau_x, zTau_y,    &
      &                 Niter=20, rad_sw=zRsw, rad_lw=zRlw, T_s=zTs )
   PRINT *, ''
   PRINT *, ' Wind speed at zu       =', REAL(SQRT(zU_zu*zU_zu + zV_zu*zV_zu),4), ' m/s'
   PRINT *, '    SST                 =', REAL(zsst,4), ' K' 
   PRINT *, ' Air temperature at zt  =', REAL(zt_zt,4), ' K' ; PRINT *, ''
   PRINT *, ' Sensible heat flux: QH =', REAL(zQH,4), ' W/m**2'
   PRINT *, '  Latent  heat flux: QL =', REAL(zQL,4), ' W/m**2'
   PRINT *, ' Skin temperature: SSST =', REAL(zTs-rt0,4), ' K'
   PRINT *, ' Tau_x =', REAL(zTau_x,4), ' N/m**2'
   PRINT *, ' Tau_y =', REAL(zTau_y,4), ' N/m**2'
   PRINT *, ' Tau   =', REAL(SQRT(zTau_x*zTau_x + zTau_y*zTau_y),4), ' N/m**2' ; PRINT *, ''
   PRINT *, ''

   PRINT *, ''
   PRINT *, ' *********** NCAR *****************'
   CALL aerobulk_model( 'ncar', 2._wp, 10._wp, zsst, zt_zt, &
      &                 zq_zt, zU_zu, zV_zu, zslp,   &
      &                 zQL, zQH, zTau_x, zTau_y,    &
      &                 Niter=20 )
   PRINT *, ''
   PRINT *, ' Wind speed at zu       =', REAL(SQRT(zU_zu*zU_zu + zV_zu*zV_zu),4), ' m/s'
   PRINT *, '    SST                 =', REAL(zsst,4), 'K'
   PRINT *, ' Air temperature at zt  =', REAL(zt_zt,4), 'K' ; PRINT *, ''
   PRINT *, ' Sensible heat flux: QH =', REAL(zQH,4), ' W/m**2'
   PRINT *, '  Latent  heat flux: QL =', REAL(zQL,4), ' W/m**2'
   PRINT *, ' Tau_x =', REAL(zTau_x,4), ' N/m**2'
   PRINT *, ' Tau_y =', REAL(zTau_y,4), ' N/m**2'
   PRINT *, ' Tau   =', REAL(SQRT(zTau_x*zTau_x + zTau_y*zTau_y),4), ' N/m**2' ; PRINT *, ''
   PRINT *, ''

END PROGRAM EXAMPLE_CALL_AEROBULK_COMPUTE
