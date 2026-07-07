! AeroBulk / 2019 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_blk_grachev07
   !!====================================================================================
   !!           Grachev et al., 2007
   !!
   !!   For now: only focus on `psi` function for stable SBL !
   !!
   !!    SHEBA flux–profile relationships in the stable atmospheric boundary layer,
   !!      Boundary-Layer Meteorology (2007)  |   DOI 10.1007/s10546-007-9177-6
   !!
   !!
   !!              (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, July 2026
   !!
   !!====================================================================================
   USE mod_const, ONLY: wp, rpi !: physical and othe constants

   IMPLICIT NONE

   PRIVATE

   INTERFACE psi_m_grachev07
      MODULE PROCEDURE psi_m_grachev07_scl, psi_m_grachev07_vct
   END INTERFACE psi_m_grachev07
   INTERFACE psi_h_grachev07
      MODULE PROCEDURE psi_h_grachev07_scl, psi_h_grachev07_vct
   END INTERFACE psi_h_grachev07

   PUBLIC :: psi_m_grachev07, psi_h_grachev07

   CHARACTER(len=8), PARAMETER :: clbl = 'GRACHEV07'

CONTAINS





   !!===============================================================================================
   FUNCTION psi_m_grachev07_scl( pzeta )
      !!--------------------------------------------------------------------------------------------
      !!--------------------------------------------------------------------------------------------
      REAL(wp)             :: psi_m_grachev07_scl
      REAL(wp), INTENT(in) :: pzeta
      !!--------------------------------------------------------------------------------------------
      REAL(wp) :: zta, zx, zpsi_u, zpsi_s
      !!----------------------------------------------------------------------------------
      zta = pzeta
      !
      ! Unstable stratification:
      zx = ABS(1._wp - 16._wp*zta)**0.25_wp              !  (16 here, not 15!)

      zpsi_u = LOG( 0.5_wp*(1._wp + zx*zx) ) + 2._wp*LOG( 0.5_wp*(1._wp + zx) ) - 2._wp*ATAN( zx ) + 0.5_wp*rpi  ! Eq.(30) Jordan et al. 1999:

      ! Stable stratification:
      zpsi_s = 1._wp + 6.5_wp*zta*( 1._wp + zta )**0.3333333_wp / ( 1.3_wp + zta )     ! Eq.(9.a) Grachev et al. 2007

      !! Combine:
      psi_m_grachev07_scl = MERGE( zpsi_u  ,  -1._wp*zpsi_s  ,  zta < 0._wp )   ! Unstable <= zta < 0)  |  Stable <= zta > 0)
      !!
   END FUNCTION psi_m_grachev07_scl

   FUNCTION psi_m_grachev07_vct( pzeta )
      !!--------------------------------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in)             :: pzeta
      REAL(wp), DIMENSION(SIZE(pzeta,1),SIZE(pzeta,2)) :: psi_m_grachev07_vct
      !!--------------------------------------------------------------------------------------------
      INTEGER  ::   ji, jj     ! dummy loop indices
      !!--------------------------------------------------------------------------------------------
      DO jj = 1, SIZE(pzeta,2)
         DO ji = 1, SIZE(pzeta,1)
            psi_m_grachev07_vct(ji,jj) = psi_m_grachev07_scl( pzeta(ji,jj) )
         END DO
      END DO
   END FUNCTION psi_m_grachev07_vct


   !!===============================================================================================


   !!===============================================================================================
   FUNCTION psi_h_grachev07_scl( pzeta )
      !!--------------------------------------------------------------------------------------------
      !!--------------------------------------------------------------------------------------------
      REAL(wp)             :: psi_h_grachev07_scl
      REAL(wp), INTENT(in) :: pzeta
      !!--------------------------------------------------------------------------------------------
      REAL(wp) :: zta, zx, zpsi_u, zpsi_s
      !!----------------------------------------------------------------------------------
      zta = pzeta

      ! Unstable stratification:
      zx = ABS(1._wp - 16._wp*zta)**0.25_wp              !  (16 here, not 15!)

      zpsi_u =   2._wp*LOG( 0.5_wp * (1._wp + zx*zx) )  ! Eq.(31) Jordan et al. 1999

      ! Stable stratification (identical to Psi_m!):
      !zpsi_s = - ( 0.7_wp*zta + 0.75_wp*(zta - 14.3_wp)*EXP( -0.35_wp*zta) + 10.7_wp )  ! Eq.(33) Jordan et al. 1999
      zpsi_s = 1._wp + 5._wp*zta*( 1._wp + zta ) / ( 1._wp + 3._wp*zta + zta*zta )       ! Eq.(9.b) Grachev et al. 2007

      !! Combine:
      psi_h_grachev07_scl = MERGE( zpsi_u  ,  -1._wp*zpsi_s  ,  zta < 0._wp )   ! Unstable <= zta < 0)  |  Stable <= zta > 0)
      !!
   END FUNCTION psi_h_grachev07_scl

   FUNCTION psi_h_grachev07_vct( pzeta )
      !!--------------------------------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in)             :: pzeta
      REAL(wp), DIMENSION(SIZE(pzeta,1),SIZE(pzeta,2)) :: psi_h_grachev07_vct
      !!--------------------------------------------------------------------------------------------
      INTEGER  ::   ji, jj     ! dummy loop indices
      !!--------------------------------------------------------------------------------------------
      DO jj = 1, SIZE(pzeta,2)
         DO ji = 1, SIZE(pzeta,1)
            psi_h_grachev07_vct(ji,jj) = psi_h_grachev07_scl( pzeta(ji,jj) )
         END DO
      END DO
   END FUNCTION psi_h_grachev07_vct

   !!===============================================================================================

   !SUBROUTINE cap_zeta( pzeta )
   !   !!--------------------------------------------------------------------------------------------
   !   REAL(wp), INTENT(inout) :: pzeta
   !   REAL(wp) ::  zta
   !   !!--------------------------------------------------------------------------------------------
   !   !#fixme: Jean Bildot & Sam Hatfield @ ECMWF, complain that
   !   !        `EXP(-0.35_wp*zta)` later blows up in single precision when unstable with big `zta`
   !   !#fixme: LB suggests:
   !   zta = MAX( pzeta , -50._wp ) ! => regions where `zeta<-50.` are given value -50 (still unrealistic but numerically safe?)
   !   !                            !  ==> prevents numerical problems such as overflows...
   !   zta = MIN(  zta ,   5._wp )  !`zeta` plateaus at 5 in very stable conditions (L>0 and small!), inherent to ECMWF algo!
   !   !
   !   pzeta = zta
   !END SUBROUTINE cap_zeta
   !!======================================================================

END MODULE mod_blk_grachev07
