! AeroBulk / 2020 / L. Brodeau
!
!   ORIGIN: https://github.com/brodeau/aerobulk/
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
PROGRAM test_psi_stab
   !!====================================================================================
   !!       Computes profile stability correction fonctions of all algorithms
   !!       in order to compare them and debug them....
   !!
   !!
   !!            Author: Laurent Brodeau, 2020
   !!
   !!====================================================================================
   USE mod_const                                         !: physical and othe constants
   USE mod_phymbl                                        !: thermodynamics
   USE io_ezcdf
   
   USE mod_blk_ncar,     ONLY: psi_m_ncar,  psi_h_ncar
   USE mod_blk_coare3p6, ONLY: psi_m_coare, psi_h_coare
   
   IMPLICIT NONE
   
   INTEGER,  PARAMETER :: Nz = 201
   REAL(wp), PARAMETER :: zeta_min = -10._wp, zeta_max = 10._wp
   
   INTEGER  :: jz
   REAL(wp) :: dzeta
   
   REAL(wp), DIMENSION(Nz,1)   :: vzeta
   REAL(4),  DIMENSION(Nz,1,2) :: vpsi_m, vpsi_h

   !! Aerobulk initialization:
   jpi=Nz
   jpj=1
   
   !! Building zeta axis:
   dzeta = (zeta_max - zeta_min)/(Nz - 1)
   FORALL (jz = 1:Nz)  vzeta(jz,1) = zeta_min + REAL(jz-1,wp)*dzeta

   PRINT *, ''
   PRINT *, ' *** dzeta =', dzeta
   PRINT *, ''
   PRINT *, ' *** vzeta =', vzeta



   vpsi_m(:,:,1) = REAL( psi_m_ncar( vzeta ) , 4 )
   vpsi_h(:,:,1) = REAL( psi_h_ncar( vzeta ) , 4 )

   vpsi_m(:,:,2) = REAL( psi_m_coare( vzeta ) , 4 )
   vpsi_h(:,:,2) = REAL( psi_h_coare( vzeta ) , 4 )


   PRINT *, ' *** vpsi_m =', vpsi_m
   PRINT *, ''
   PRINT *, ' *** vpsi_h =', vpsi_h

   CALL PT_SERIES( vzeta(:,1), &
      &            vpsi_m(:,1,1), 'psi.nc', 'zeta', 'Psi_m_NCAR', '[]', 'Stability profile correction for momentum', REAL(-9999.,4), &
      &             '[]', '[]',              &
      &            vdt02=vpsi_h(:,1,1), cv_dt02='Psi_h_NCAR', &
      &            vdt03=vpsi_m(:,1,2), cv_dt03='Psi_m_COARE', &
      &            vdt04=vpsi_h(:,1,2), cv_dt04='Psi_h_COARE' &      
      &  )









   !!======================================================================
END PROGRAM test_psi_stab
