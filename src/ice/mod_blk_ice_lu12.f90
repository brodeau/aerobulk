! AeroBulk / 2020 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_blk_ice_lu12
   !!====================================================================================
   !!       Computes turbulent components of surface fluxes over sea-ice
   !!       Routine turb_ice_lu12 maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, January 2020
   !!
   !!====================================================================================
   USE mod_const       !: physical and othe constants
   !USE mod_phymbl      !: thermodynamics

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: Cdn10_Lupkes2012
   
   REAL(wp), PARAMETER ::   zCe   = 2.23e-03_wp
   REAL(wp), PARAMETER ::   znu   = 1._wp
   REAL(wp), PARAMETER ::   zmu   = 1._wp
   REAL(wp), PARAMETER ::   zbeta = 1._wp
   
   !!----------------------------------------------------------------------
CONTAINS
   
   SUBROUTINE Cdn10_Lupkes2012( pcd_i_s, pfrice, pcdn )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  Cdn10_Lupkes2012  ***
      !!
      !! ** Purpose :    Recompute the neutral air-ice drag referenced at 10m
      !!                 to make it dependent on edges at leads, melt ponds and flows.
      !!                 After some approximations, this can be resumed to a dependency
      !!                 on ice concentration.
      !!
      !! ** Method :     The parameterization is taken from Lupkes et al. (2012) eq.(50)
      !!                 with the highest level of approximation: level4, eq.(59)
      !!                 The generic drag over a cell partly covered by ice can be re-written as follows:
      !!
      !!                 Cd = Cdw * (1-A) + Cdi * A + Ce * (1-A)**(nu+1/(10*beta)) * A**mu
      !!
      !!                    Ce = 2.23e-3       , as suggested by Lupkes (eq. 59)
      !!                    nu = mu = beta = 1 , as suggested by Lupkes (eq. 59)
      !!                    A is the concentration of ice minus melt ponds (if any)
      !!
      !!                 This new drag has a parabolic shape (as a function of A) starting at
      !!                 Cdw(say 1.5e-3) for A=0, reaching 1.97e-3 for A~0.5
      !!                 and going down to Cdi(say 1.4e-3) for A=1
      !!
      !!                 It is theoretically applicable to all ice conditions (not only MIZ)
      !!                 => see Lupkes et al (2013)
      !!
      !! ** References : Lupkes et al. JGR 2012 (theory)
      !!                 Lupkes et al. GRL 2013 (application to GCM)
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pcd_i_s  ! skin drag coefficient for a 100% ice covered region
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pfrice   ! ice concentration [fraction]  => at_i_b
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pcdn     ! neutral drag coefficient over sea-ice
      !!----------------------------------------------------------------------
      REAL(wp)            ::   zcoef
      !!----------------------------------------------------------------------
      zcoef = znu + 1._wp / ( 10._wp * zbeta )
      
      !! We are not an AGCM, we are an OGCM!!! => we drop term "(1 - A)*Cd_w"
      !!  => so we keep only the two last rhs terms of Eq.(1) of Lupkes et al, 2013 divided by "A"...
      
      pcdn(:,:) = pcd_i_s +  zCe * pfrice(:,:)**(zmu - 1._wp) * (1._wp - pfrice(:,:))**zcoef

      !! => seems okay for winter 100% sea-ice as second rhs term vanishes as pfrice == 1....

      
   END SUBROUTINE Cdn10_Lupkes2012


   !!======================================================================
END MODULE mod_blk_ice_lu12
