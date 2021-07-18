! AeroBulk / 2015 / L. Brodeau

PROGRAM EXAMPLE_CALL_AEROBULK_COMPUTE

   USE mod_aerobulk
   USE mod_const
   USE mod_phymbl, ONLY: Theta_from_z_P0_T_q

   
   IMPLICIT NONE

   INTEGER, PARAMETER :: nx = 2, &
      &                  ny = 1
   
   REAL(wp), PARAMETER :: zt = 2.  , & ! ref. height for T and q
      &                   zu = 10.     ! ref. height for wind
   
   REAL(wp), DIMENSION(nx,ny) :: zsst, zt_zt, zq_zt, zU_zu, zV_zu, zslp, &
      &                         zRsw, zRlw,           &
      &                         zQL, zQH, zTau_x, zTau_y, zE, zTs, ztmp

   
   PRINT *, ''

   zsst = 22.
   PRINT *, ' *** SST = ', zsst, ' deg.C'
   zsst  = zsst + rt0  ! [K]

   zt_zt(1,1) = rt0 + 20.  ! [K]
   zt_zt(2,1) = rt0 + 25.  ! [K] ! second case is stable ABL as t_air > SST (25>22)!


   zq_zt = 0.012         ! [kg/kg]
   zU_zu = 4.            ! [m/s]
   zV_zu = 9.            ! [m/s]
   zslp  = 101000.       ! [Pa]

   zRsw  = 0. ! (night)  ! [W/m^2]
   zRlw  = 350.          ! [W/m^2]


   PRINT *, ''
   PRINT *, ' *********** COARE 3.6 *****************'
   CALL aerobulk_model( 1, 1, 'coare3p6', zt, zu, zsst, zt_zt,     &
      &                 zq_zt, zU_zu, zV_zu, zslp,                  &
      &                 zQL, zQH, zTau_x, zTau_y, zE,               &
      &                 Niter=20, rad_sw=zRsw, rad_lw=zRlw, T_s=zTs )
   PRINT *, ''
   PRINT *, ' Wind speed at zu       =', REAL(SQRT(zU_zu*zU_zu + zV_zu*zV_zu),4), ' m/s'
   PRINT *, '    SST                 =', REAL(zsst-rt0,4), ' deg.C'
   PRINT *, ' Abs. temperature at zt =', REAL(zt_zt-rt0,4), ' deg.C'
   ztmp = Theta_from_z_P0_T_q( zt, zslp, zt_zt, zq_zt ) - rt0
   PRINT *, ' Pot. temperature at zt =', REAL(ztmp,4), ' deg.C'  ; PRINT *, ''
   PRINT *, ' Sensible heat flux: QH =', REAL(zQH,4), ' W/m**2'
   PRINT *, '  Latent  heat flux: QL =', REAL(zQL,4), ' W/m**2'
   PRINT *, '  Evaporation:     Evap =', REAL(zE*3600.*24.,4), ' mm/day'
   PRINT *, ' Skin temperature: SSST =', REAL(zTs-rt0,4), ' deg.C'
   PRINT *, ' Tau_x =', REAL(zTau_x,4), ' N/m**2'
   PRINT *, ' Tau_y =', REAL(zTau_y,4), ' N/m**2'
   PRINT *, ' Tau   =', REAL(SQRT(zTau_x*zTau_x + zTau_y*zTau_y),4), ' N/m**2' ; PRINT *, ''
   PRINT *, ''
   

   !l_1st_call_ab_init=.TRUE.
   PRINT *, ''
   PRINT *, ' *********** ECMWF *****************'
   CALL aerobulk_model( 1, 1, 'ecmwf', zt, zu, zsst, zt_zt,        &
      &                 zq_zt, zU_zu, zV_zu, zslp,                  & 
      &                 zQL, zQH, zTau_x, zTau_y, zE,               &
      &                 Niter=20, rad_sw=zRsw, rad_lw=zRlw, T_s=zTs )
   PRINT *, ''
   PRINT *, ' Wind speed at zu       =', REAL(SQRT(zU_zu*zU_zu + zV_zu*zV_zu),4), ' m/s'
   PRINT *, '    SST                 =', REAL(zsst-rt0,4), ' deg.C' 
   PRINT *, ' Abs. temperature at zt =', REAL(zt_zt-rt0,4), ' deg.C'
   ztmp = Theta_from_z_P0_T_q( zt, zslp, zt_zt, zq_zt ) - rt0
   PRINT *, ' Pot. temperature at zt =', REAL(ztmp,4), ' deg.C'  ; PRINT *, ''
   PRINT *, ' Sensible heat flux: QH =', REAL(zQH,4), ' W/m**2'
   PRINT *, '  Latent  heat flux: QL =', REAL(zQL,4), ' W/m**2'
   PRINT *, '  Evaporation:     Evap =', REAL(zE*3600.*24.,4), ' mm/day'
   PRINT *, ' Skin temperature: SSST =', REAL(zTs-rt0,4), ' deg.C'
   PRINT *, ' Tau_x =', REAL(zTau_x,4), ' N/m**2'
   PRINT *, ' Tau_y =', REAL(zTau_y,4), ' N/m**2'
   PRINT *, ' Tau   =', REAL(SQRT(zTau_x*zTau_x + zTau_y*zTau_y),4), ' N/m**2' ; PRINT *, ''
   PRINT *, ''

   !l_1st_call_ab_init=.TRUE.
   PRINT *, ''
   PRINT *, ' *********** NCAR *****************'
   CALL aerobulk_model( 1, 1, 'ncar', zt, zu, zsst, zt_zt, &
      &                 zq_zt, zU_zu, zV_zu, zslp,          &
      &                 zQL, zQH, zTau_x, zTau_y, zE,       &
      &                 Niter=20 )
   PRINT *, ''
   PRINT *, ' Wind speed at zu       =', REAL(SQRT(zU_zu*zU_zu + zV_zu*zV_zu),4), ' m/s'
   PRINT *, '    SST                 =', REAL(zsst-rt0,4), ' deg.C'
   PRINT *, ' Abs. temperature at zt =', REAL(zt_zt-rt0,4), ' deg.C'
   ztmp = Theta_from_z_P0_T_q( zt, zslp, zt_zt, zq_zt ) - rt0
   PRINT *, ' Pot. temperature at zt =', REAL(ztmp,4), ' deg.C'  ; PRINT *, ''
   PRINT *, ' Sensible heat flux: QH =', REAL(zQH,4), ' W/m**2'
   PRINT *, '  Latent  heat flux: QL =', REAL(zQL,4), ' W/m**2'
   PRINT *, '  Evaporation:     Evap =', REAL(zE*3600.*24.,4), ' mm/day'
   PRINT *, ' Tau_x =', REAL(zTau_x,4), ' N/m**2'
   PRINT *, ' Tau_y =', REAL(zTau_y,4), ' N/m**2'
   PRINT *, ' Tau   =', REAL(SQRT(zTau_x*zTau_x + zTau_y*zTau_y),4), ' N/m**2' ; PRINT *, ''
   PRINT *, ''


   
   PRINT *, ''
   PRINT *, ' *********** ANDREAS *****************'
   CALL aerobulk_model( 1, 1, 'andreas', zt, zu, zsst, zt_zt, &
      &                 zq_zt, zU_zu, zV_zu, zslp,          &
      &                 zQL, zQH, zTau_x, zTau_y, zE,       &
      &                 Niter=20 )
   PRINT *, ''
   PRINT *, ' Wind speed at zu       =', REAL(SQRT(zU_zu*zU_zu + zV_zu*zV_zu),4), ' m/s'
   PRINT *, '    SST                 =', REAL(zsst-rt0,4), ' deg.C'
   PRINT *, ' Abs. temperature at zt =', REAL(zt_zt-rt0,4), ' deg.C'
   ztmp = Theta_from_z_P0_T_q( zt, zslp, zt_zt, zq_zt ) - rt0
   PRINT *, ' Pot. temperature at zt =', REAL(ztmp,4), ' deg.C'  ; PRINT *, ''
   PRINT *, ' Sensible heat flux: QH =', REAL(zQH,4), ' W/m**2'
   PRINT *, '  Latent  heat flux: QL =', REAL(zQL,4), ' W/m**2'
   PRINT *, '  Evaporation:     Evap =', REAL(zE*3600.*24.,4), ' mm/day'
   PRINT *, ' Tau_x =', REAL(zTau_x,4), ' N/m**2'
   PRINT *, ' Tau_y =', REAL(zTau_y,4), ' N/m**2'
   PRINT *, ' Tau   =', REAL(SQRT(zTau_x*zTau_x + zTau_y*zTau_y),4), ' N/m**2' ; PRINT *, ''
   PRINT *, ''


   
END PROGRAM EXAMPLE_CALL_AEROBULK_COMPUTE
