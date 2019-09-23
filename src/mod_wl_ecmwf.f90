! AeroBulk / 2016 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_wl_ecmwf
   !!====================================================================================
   !!       Warm-Layer correction of SST
   !!
   !!       Routine "wl_ecmwf" maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, 2019
   !!
   !!====================================================================================
   USE mod_const   !: physical and othe constants
   USE mod_phymbl  !: thermodynamics

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: WL_ECMWF

   !!  Warm-Layer related parameters:
   REAL(wp), PARAMETER :: rd0  = 3.        !: Depth scale [m] of warm layer, "d" in Eq.11 (Zeng & Beljaars 2005)

   REAL(wp), PARAMETER :: z_sst  = 1._wp    !: depth at which bulk SST is taken... [m]
   
   REAL(wp), PARAMETER :: zroadrw = rho0_a/rho0_w , &         ! Density ratio
      &                   zRhoCp_w = rho0_w*rCp0_w


   REAL(wp), PARAMETER :: rNu0 = 1.0       !:  be closer to COARE3p6 ???!LOLO
   !REAL(wp), PARAMETER :: rNu0 = 0.5       !: Nu (exponent of temperature profile) Eq.11
   !                                       !: (Zeng & Beljaars 2005) !: set to 0.5 instead of
   !                                       !: 0.3 to respect a warming of +3 K in calm
   !                                       !: condition for the insolation peak of +1000W/m^2
   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE WL_ECMWF( pQsw, pQnsol, pustar, pSST, rdt, pdT )
      !!---------------------------------------------------------------------
      !!
      !!  Warm-Layer scheme according to Zeng & Beljaars, 2005 (GRL)
      !!  " A prognostic scheme of sea surface skin temperature for modeling and data assimilation "
      !!
      !!    As included in IFS Cy45r1   /  E.C.M.W.F.
      !!     ------------------------------------------------------------------
      !!
      !!  **   INPUT:
      !!     *pQsw*       surface net solar radiation into the ocean     [W/m^2] => >= 0 !
      !!     *pQnsol*     surface net non-solar heat flux into the ocean [W/m^2] => normally < 0 !
      !!     *pustar*     friction velocity u*                           [m/s]
      !!     *pSST*       bulk SST at depth z_sst                        [K]
      !!     *rdt*        physical time step between two successive calls to this routine [s]
      !!
      !!   **  INPUT/OUTPUT:
      !!     *pdT*  : as input =>  previous estimate of dT warm-layer
      !!             as output =>  new estimate of dT warm-layer
      !!
      !!------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pQsw     ! surface net solar radiation into the ocean [W/m^2]     => >= 0 !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pQnsol   ! surface net non-solar heat flux into the ocean [W/m^2] => normally < 0 !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pustar   ! friction velocity [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pSST     ! bulk SST at depth z_sst [K]
      REAL(wp),                     INTENT(in)  :: rdt      ! physical time step between two successive call to this routine [s]
      !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pdT    ! difference between "almost surface (right below viscous layer) and depth of bulk SST
      !
      INTEGER :: ji,jj
      !
      REAL(wp) :: &
         & zdz,    & !: thickness of the warm-layer [m]
         & zalpha_w, & !: thermal expansion coefficient of sea-water
         & ZSRD,   &
         & dT_wl,   & ! temp. diff. between "almost surface (right below viscous layer) and bottom of WL
         & zfr,zdL,zdL2, ztmp, &
         & ZPHI, &
         & zus_a, zusw, zusw2, &
         & flg, zQabs, ZL1, ZL2
      !!---------------------------------------------------------------------


      DO jj = 1, jpj
         DO ji = 1, jpi

            zdz = rd0 ! first guess for warm-layer depth (and unique..., less advanced than COARE3p6 !)
            
            ! dT_wl is the difference between "almost surface (right below viscous layer) and bottom of WL (here zdz)
            ! pdT         "                          "                                    and depth of bulk SST (here z_sst)!
            !! => but of course in general the bulk SST is taken shallower than zdz !!! So correction less pronounced!
            !! => so here since pdT is difference between surface and z_sst, need to increase fof dT_wl !
            flg = 0.5_wp + SIGN( 0.5_wp , z_sst-zdz )               ! => 1 when z_sst>zdz (pdT(ji,jj) = dT_wl) | 0 when z_s$
            dT_wl = pdT(ji,jj) / ( flg + (1._wp-flg)*z_sst/zdz )
            !PRINT *, 'LOLO/mod_wl_ecmwf.f90: dT_wl2=', dT_wl
            !PRINT *, ''

            zalpha_w = alpha_sw( pSST(ji,jj) ) ! thermal expansion coefficient of sea-water (SST accurate enough!)
            
            
            ! *** zfr = Fraction of solar radiation absorbed in warm layer (-)
            zfr = 1._wp - 0.28_wp*EXP(-71.5_wp*zdz) - 0.27_wp*EXP(-2.8_wp*zdz) - 0.45_wp*EXP(-0.07_wp*zdz)  !: Eq. 8.157            
            
            zQabs = zfr*pQsw(ji,jj) + pQnsol(ji,jj)       ! tot heat absorbed in warm layer

            zusw  = MAX(pustar(ji,jj), 1.E-4_wp)*SQRT(zroadrw)    ! u* in the water
            zusw2 = zusw*zusw

            
            !! *** 1st rhs term in eq. 8.156 (IFS doc Cy45r1):
            ZL1 = zQabs / ( zdz * zRhoCp_w * rNu0 ) * (rNu0 + 1._wp)


            !! Buoyancy flux and stability parameter (zdl = -z/L) in water
            ZSRD = zQabs/zRhoCp_w
            !
            flg = 0.5_wp + SIGN(0.5_wp, ZSRD)  ! ZSRD > 0. => 1.  / ZSRD < 0. => 0.
            ztmp = MAX(dT_wl,0._wp)
            zdl = (flg+1._wp) * ( zusw2 * SQRT(ztmp/(5._wp*zdz*grav*zalpha_w/rNu0)) ) & ! (dT_wl > 0.0 .AND. ZSRD < 0.0)
               &  +    flg    *  ZSRD                                                                  !   otherwize
            !
            zus_a = MAX( pustar(ji,jj), 1.E-4_wp )
            zdL = zdz*vkarmn*grav/(zroadrw)**1.5_wp*zalpha_w*zdL/(zus_a*zus_a*zus_a)

            !! Stability function Phi_t(-z/L) (zdL is -z/L) :
            flg = 0.5_wp + SIGN(0.5_wp, zdL)  ! zdl > 0. => 1.  / zdl < 0. => 0.
            zdL2 = zdL*zdL            
            ZPHI =    flg      * ( 1._wp + (5._wp*zdL + 4._wp*zdL2)/(1._wp + 3._wp*zdL + 0.25_wp*zdL2) ) &  ! (zdL > 0) Takaya et al.
               & + (flg+1._wp) * ( 1._wp/SQRT(1._wp - 16._wp*(-ABS(zdL))) )        ! (zdl < 0) Eq. 8.136
            !! FOR zdL > 0.0, old relations:
            !         ZPHI = 1.+5._wp*zdL                                ! Eq. 8.136 (Large et al. 1994)
            !         ZPHI = 1.+5.0*(zdL+zdL**2)/(1.0+3.0*zdL+zdL**2) ! SHEBA, Grachev et al. 2007
            
            !! *** 2nd rhs term in eq. 8.156 (IFS doc Cy45r1):
            ZL2 = - (rNu0 + 1._wp) * vkarmn * zusw / ( zdz * ZPHI )
                        
            ! Forward time / explicit solving of eq. 8.156 (IFS doc Cy45r1): (f_n+1 == pdT(ji,jj) ; f_n == dT_wl)
            dT_wl = MAX ( dT_wl + rdt*ZL1 + rdt*ZL2*dT_wl , 0._wp )

            ! dT_wl is the difference between "almost surface (right below viscous layer) and bottom of WL (here zdz)
            !! => but of course in general the bulk SST is taken shallower than zdz !!! So correction less pronounced!

            flg = 0.5_wp + SIGN( 0.5_wp , z_sst-zdz )               ! => 1 when z_sst>zdz (pdT(ji,jj) = dT_wl) | 0 when z_s$
            pdT(ji,jj) = dT_wl * ( flg + (1._wp-flg)*z_sst/zdz )
            

         END DO ! DO ji = 1, jpi
      END DO ! DO jj = 1, jpj

   END SUBROUTINE WL_ECMWF

   !!======================================================================
END MODULE mod_wl_ecmwf
