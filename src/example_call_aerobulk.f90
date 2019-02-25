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
   zU_zu = 5.            ! [m/s]
   zV_zu = 0.            ! [m/s]
   zslp  = 101000.       ! [Pa]
   
   zRsw  = 0. ! (night)  ! [W/m^2]
   zRlw  = 350.          ! [W/m^2]
   
   CALL aerobulk_model( 'coare', 2._wp, 10._wp, zsst, zt_zt, &
      &                 zq_zt, zU_zu, zV_zu, zslp,   &
      &                 zQL, zQH, zTau_x, zTau_y,    &
      &                 Niter=10, rad_sw=zRsw, rad_lw=zRlw, T_s=zTs )

   
   PRINT *, ''
   PRINT *, ' Sensible heat flux: QH =', zQH ; PRINT *, ''
   PRINT *, '  Latent  heat flux: QL =', zQL ; PRINT *, ''
   PRINT *, ' Skin temperature: SSST =', zTs-rt0 ; PRINT *, ''
   PRINT *, ' Tau_x=', zTau_x ; PRINT *, ''
   PRINT *, ' Tau_y=', zTau_y ; PRINT *, ''
   
   
END PROGRAM EXAMPLE_CALL_AEROBULK_COMPUTE
