! AeroBulk / 2015 / L. Brodeau (brodeau@gmail.com)
! https://sourceforge.net/p/aerobulk

PROGRAM EXAMPLE_CALL_AEROBULK_COMPUTE

   USE mod_aerobulk
   USE mod_const
   
   IMPLICIT NONE


   INTEGER, PARAMETER :: nx = 2, &
      &                  ny = 1
   
   REAL(wp), DIMENSION(nx,ny) :: zsst, zt_zt, zq_zt, zU_zu, zV_zu, zslp, &
      &                         zRsw, zRlw,           &
      &                         zQL, zQH, zTau_x, zTau_y

   
   zsst  = 273.15 + 22.
   zt_zt = 273.15 + 20.
   zq_zt = 0.012
   zU_zu = 5.
   zV_zu = 0.
   zslp  = 101000.
   
   zRsw  = 0. ! (night)
   zRlw  = 350. 
   
   CALL aerobulk_model( 'coare', 2._wp, 10._wp, zsst, zt_zt, &
      &                 zq_zt, zU_zu, zV_zu, zslp,   &
      &                 zQL, zQH, zTau_x, zTau_y,    &
      &                 Niter=10, rad_sw=zRsw, rad_lw=zRlw )

   
   PRINT *, ''
   PRINT *, ' QH=', zQH ; PRINT *, ''
   PRINT *, ' QL=', zQL ; PRINT *, ''
   PRINT *, ' Tau_x=', zTau_x ; PRINT *, ''
   PRINT *, ' Tau_y=', zTau_y ; PRINT *, ''
   
   
END PROGRAM EXAMPLE_CALL_AEROBULK_COMPUTE
