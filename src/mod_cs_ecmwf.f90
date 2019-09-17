! AeroBulk / 2016 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_cs_ecmwf
   !!====================================================================================
   !!       Cool-Skin correction of SST
   !!
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (2m) to zu (10m) if needed
   !!   * the effective bulk wind speed at 10m U_blk
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!    Using the bulk formulation/param. of ECMWF
   !!      + consideration of cool-skin parametrization (Fairall et al. 1996)
   !!
   !!       Routine "cs_ecmwf" maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, 2019
   !!
   !!====================================================================================
   USE mod_const   !: physical and othe constants
   USE mod_phymbl  !: thermodynamics

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: CS_ECMWF

   !! Cool-Skin / Warm-Layer related parameters:
   REAL(wp), PARAMETER :: rNu0 = 0.5       !: Nu (exponent of temperature profile) Eq.11
   !                                       !: (Zeng & Beljaars 2005) !: set to 0.5 instead of
   !                                       !: 0.3 to respect a warming of +3 K in calm
   !                                       !: condition for the insolation peak of +1000W/m^2
   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE CS_ECMWF( pQsw, pQnsol, pustar, pSST, pdT )
      !!---------------------------------------------------------------------
      !!
      !!  Cool-Skin scheme according to Fairall et al. 1996
      !!  " A prognostic scheme of sea surface skin temperature for modeling and data assimilation "
      !!
      !!    As included in IFS Cy40   /  E.C.M.W.F.
      !!     ------------------------------------------------------------------
      !!
      !!  **   INPUT:
      !!     *pQsw*       surface net solar radiation into the ocean     [W/m^2] => >= 0 !
      !!     *pQnsol*     surface net non-solar heat flux into the ocean [W/m^2] => normally < 0 !
      !!     *pustar*     friction velocity u*                           [m/s]
      !!     *pSST*       bulk SST at depth z_sst                        [K]
      !!
      !!   **  INPUT/OUTPUT:
      !!     *pdT*  : as input =>  previous estimate of dT cool-skin
      !!             as output =>  new estimate of dT cool-skin
      !!
      !!------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pQsw     ! surface net solar radiation into the ocean [W/m^2]     => >= 0 !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pQnsol   ! surface net non-solar heat flux into the ocean [W/m^2] => normally < 0 !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pustar
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pSST     ! bulk SST at depth z_sst [K]
      !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pdT    ! ?
      !
      INTEGER :: ji,jj
      !
      REAL(wp) :: &
         & Ts,       & !: skin temperature ( = SST + dT_coolskin )
         & zalpha_w, & !: thermal expansion coefficient of sea-water
         & zRhoCp_w, &
         & ZCON2, zQnsol ,zQnet, zlamb, zdelta,&
         & ZROADRW, &
         & zfs
      !
      REAL(wp), DIMENSION(jpi,jpj) :: zus_w, zus_w2  !: u* and u*^2 in water

      !
      !!------------------------------------------------------------------
      !
      !     1. Initialize constants for ocean warm layer and cool skin
      !
      !     1.1 General
      !
      ZROADRW = rho0_a/rho0_w          ! Density ratio                      (-)
      zRhoCp_w = rho0_w*rCp0_w
      !
      !     1.3 Cool skin parametrization constants
      ZCON2 = 16._wp*grav*zRhoCp_w*rnu0_w**3/(rk0_w*rk0_w)
      !
      ! Friction velocities
      ! "MAX( pustar(:,:), 1.E-4)" is u* in the air !
      zus_w(:,:)  = MAX( pustar(:,:), 1.E-4_wp)*SQRT(ZROADRW)       ! u* in the water
      zus_w2(:,:) = zus_w(:,:)*zus_w(:,:)
      !
      !pdT(:,:) = 0._wp

      !  3. Cool skin (Fairall et al. 1996)
      !------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi

            Ts = pSST(ji,jj) + pdT(ji,jj) ! Skin temperature

            zalpha_w = MAX( 1.E-5_wp , 1.E-5_wp*(Ts - rt0) ) ! thermal expansion coefficient of water


            ! Non-solar heat loss to the atmosphere:
            zQnsol = MAX( 1._wp , - pQnsol(ji,jj) )

            zlamb = 6._wp*(1._wp + (zQnsol*zalpha_w*ZCON2/(zus_w2(ji,jj)*zus_w2(ji,jj)))**0.75)**(-1._wp/3._wp)

            zdelta = zlamb*rnu0_w/zus_w(ji,jj)

            !    Solar absorption
            zfs   = 0.065_wp + 11._wp*zdelta - (6.6E-5_wp/zdelta)*(1._wp - EXP(-zdelta/8.E-4_wp))  ! Eq. 8.131 / IFS cy40r1, doc, Part IV,
            zfs   = MAX(zfs , 0.01_wp)
            zQnet = MAX( 1._wp , -zfs*pQsw(ji,jj) + zQnsol )
            pdT(ji,jj) = -zdelta*zQnet/rk0_w

         END DO
      END DO


   END SUBROUTINE CS_ECMWF

   !!======================================================================
END MODULE mod_cs_ecmwf
