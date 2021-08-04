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

   USE mod_blk_ncar,     ONLY: psi_m_ncar,    psi_h_ncar
   USE mod_common_coare, ONLY: psi_m_coare,   psi_h_coare
   USE mod_blk_ecmwf,    ONLY: psi_m_ecmwf,   psi_h_ecmwf
   USE mod_blk_andreas,  ONLY: psi_m_andreas, psi_h_andreas

   IMPLICIT NONE

   CHARACTER(len=256), PARAMETER :: cf_out = 'psi.nc'

   INTEGER,  PARAMETER :: Nz = 1001, &
      &                   nalgos = 4
   REAL(wp), PARAMETER :: zeta_min = -15._wp, zeta_max = 15._wp

   INTEGER  :: jz
   REAL(wp) :: dzeta

   REAL(wp), DIMENSION(Nz,1)   :: vzeta
   REAL(4),  DIMENSION(Nz,1,nalgos) :: vpsi_m, vpsi_h

   !! Aerobulk initialization:
   !! ...

   !! Building zeta axis:
   dzeta = (zeta_max - zeta_min)/(Nz - 1)
   FORALL (jz = 1:Nz)  vzeta(jz,1) = zeta_min + REAL(jz-1,wp)*dzeta

   !PRINT *, ''
   !PRINT *, ' *** dzeta =', dzeta
   !PRINT *, ''
   !PRINT *, ' *** vzeta =', vzeta



   vpsi_m(:,:,1) = REAL( psi_m_ncar( vzeta ) , 4 )
   vpsi_h(:,:,1) = REAL( psi_h_ncar( vzeta ) , 4 )

   vpsi_m(:,:,2) = REAL( psi_m_coare( vzeta ) , 4 )
   vpsi_h(:,:,2) = REAL( psi_h_coare( vzeta ) , 4 )

   vpsi_m(:,:,3) = REAL( psi_m_ecmwf( vzeta ) , 4 )
   vpsi_h(:,:,3) = REAL( psi_h_ecmwf( vzeta ) , 4 )

   vpsi_m(:,:,4) = REAL( psi_m_andreas( vzeta ) , 4 )
   vpsi_h(:,:,4) = REAL( psi_h_andreas( vzeta ) , 4 )


   !PRINT *, ' *** vpsi_m =', vpsi_m
   !PRINT *, ''
   !PRINT *, ' *** vpsi_h =', vpsi_h

   CALL PT_SERIES( vzeta(:,1), vpsi_m(:,1,1), TRIM(cf_out), 'zeta', 'Psi_m_NCAR', &
      &            '[]', 'Stability profile correction for momentum', REAL(-9999.,4), &
      &            '[]', '[]',                                   &
      &            vdt02=vpsi_h(:,1,1), cv_dt02='Psi_h_NCAR',    &
      &            vdt03=vpsi_m(:,1,2), cv_dt03='Psi_m_COARE',   &
      &            vdt04=vpsi_h(:,1,2), cv_dt04='Psi_h_COARE',   &
      &            vdt05=vpsi_m(:,1,3), cv_dt05='Psi_m_ECMWF',   &
      &            vdt06=vpsi_h(:,1,3), cv_dt06='Psi_h_ECMWF',   &
      &            vdt07=vpsi_m(:,1,4), cv_dt07='Psi_m_ANDREAS', &
      &            vdt08=vpsi_h(:,1,4), cv_dt08='Psi_h_ANDREAS'  &
      &  )


   PRINT *, ''
   PRINT *, 'Voila! Check "'//TRIM(cf_out)//'" out!'
   PRINT *, 'Alternatively, plot the comparison of all Psi functions:'
   PRINT *, ''
   PRINT *, '    python ./misc/python/plot_Psi_profiles.py psi.nc'
   PRINT *, ''
   PRINT *, ' => will generate figures "Comparaison_Psi_m.png" and "Comparaison_Psi_h.png"...'
   PRINT *, ''


   !!======================================================================
END PROGRAM test_psi_stab
