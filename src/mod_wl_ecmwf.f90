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
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (2m) to zu (10m) if needed
   !!   * the effective bulk wind speed at 10m U_blk
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!    Using the bulk formulation/param. of ECMWF
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
   REAL(wp), PARAMETER :: rd0  = 1.        !: Depth scale [m], "d" in Eq.11 (Zeng & Beljaars 2005)

   REAL(wp), PARAMETER :: dz_max = 20._wp  !: maximum depth of warm layer (adjustable)
   
   REAL(wp), PARAMETER :: rNu0 = 0.5       !: Nu (exponent of temperature profile) Eq.11
   !                                       !: (Zeng & Beljaars 2005) !: set to 0.5 instead of
   !                                       !: 0.3 to respect a warming of +3 K in calm
   !                                       !: condition for the insolation peak of +1000W/m^2

   INTEGER, PARAMETER :: nb_itt_wl = 10    !: number of sub-itterations for solving the differential equation in warm-layer part
   !                                       !:  => use "nb_itt_wl = 1" for No itteration! => way cheaper !!!
   !                                       !:    => assumes balance between the last 2 terms of Eq.11 (lhs of eq.11 = 0)
   !                                       !:    => in that case no need for sub-itterations !
   !                                       !:    => ACCEPTABLE IN MOST CONDITIONS ! (UNLESS: sunny + very calm/low-wind conditions)
   !                                       !:  => Otherwize use "nb_itt_wl = 10"
   !
   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE WL_ECMWF( pQsw, pQnsol, pustar, pSST, rdt, pdT )
      !!---------------------------------------------------------------------
      !!
      !!  Warm-Layer scheme according to Zeng & Beljaars, 2005 (GRL)
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
      !!     *rdt*        physical time step between two successive calls to this routine [s]
      !!
      !!   **  INPUT/OUTPUT:
      !!     *pdT*  : as input =>  previous estimate of dT warm-layer
      !!             as output =>  new estimate of dT warm-layer
      !!
      !!------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pQsw     ! surface net solar radiation into the ocean [W/m^2]     => >= 0 !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pQnsol   ! surface net non-solar heat flux into the ocean [W/m^2] => normally < 0 !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pustar
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pSST     ! bulk SST at depth z_sst [K]
      REAL(wp),                     INTENT(in)  :: rdt      ! physical time step between two successive call to this routine [s]
      !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pdT    ! ?
      !
      INTEGER :: ji,jj
      !
      REAL(wp) :: &
         & Ts,       & !: skin temperature ( = SST + dT_coolskin )
         & dz_wl,    & !: thickness of the warm-layer [m]
         & zalpha_w, & !: thermal expansion coefficient of sea-water
         & zRhoCp_w, &
         & ZCON3,ZCON4,ZCON5, &
         & ZSRD,ZDT,ZZ,ZEPDU2,&
         & zfr0,zdL,zdL2, ztmp, &
         & ZPHI,ZROADRW, &
         & zus_a, &
         & zsgn
      !
      REAL(wp), DIMENSION(jpi,jpj) :: zus_w, zus_w2  !: u* and u*^2 in water

      !
      !!------------------------------------------------------------------
      !
      !     1. Initialize constants for ocean warm layer and cool skin
      !
      !     1.1 General
      !
      ZEPDU2  = 0.01_wp   !    security constant for velocity**2   (m2/s2)
      ZROADRW = rho0_a/rho0_w          ! Density ratio                      (-)
      zRhoCp_w = rho0_w*rCp0_w
      !
      !     1.2C Warm layer parametrization constants
      !
      dz_wl = dz_max ! initial depth set to max value

      !LOLO: => so normally we should use dz_wl instead of rd0, but then dz_wl also needs to be updated!!!
      
      !    zfr0 = Fraction of solar radiation absorbed in warm layer (-)
      zfr0 = 1._wp - 0.28_wp*EXP(-71.5_wp*rd0) - 0.27_wp*EXP(-2.8_wp*rd0) - 0.45_wp*EXP(-0.07_wp*rd0)  !: Eq. 8.157
      
      
      !
      ZCON3 = rd0*vkarmn*grav/(ZROADRW)**1.5_wp
      ZCON4 = (rNu0 + 1._wp)*vkarmn/rd0
      ZCON5 = (rNu0 + 1._wp)/(rNu0*rd0)
      !
      ! Friction velocities
      ! "MAX( pustar(:,:), 1.E-4)" is u* in the air !
      zus_w(:,:)  = MAX( pustar(:,:), 1.E-4_wp)*SQRT(ZROADRW)       ! u* in the water
      zus_w2(:,:) = zus_w(:,:)*zus_w(:,:)
      !

      ! 2.2 Warm layer; formulation C (Xubin Zeng)
      !!--------------------------------------------

      DO jj = 1, jpj
         DO ji = 1, jpi

            Ts = pSST(ji,jj) + pdT(ji,jj) ! Skin temperature

            zalpha_w = alpha_sw( pSST(ji,jj) ) ! thermal expansion coefficient of sea-water (SST accurate enough!)

            ZDT = Ts - pSST(ji,jj)

            !! Buoyancy flux and stability parameter (zdl = -z/L) in water
            !
            !! Qt/(rho_w*Cpw):
            ZSRD = ( pQsw(ji,jj)*zfr0 + pQnsol(ji,jj) )/zRhoCp_w
            !
            zsgn = 0.5_wp + SIGN(0.5_wp, ZSRD)  ! ZSRD > 0. => 1.  / ZSRD < 0. => 0.
            ztmp = MAX(ZDT,0._wp)
            zdl = (zsgn+1._wp) * ( zus_w2(ji,jj) * SQRT(ztmp/(5._wp*rd0*grav*zalpha_w/rNu0)) ) & ! (ZDT > 0.0 .AND. ZSRD < 0.0)
               &  +    zsgn    *  ZSRD                                                                  !   otherwize
            !
            zus_a = MAX( pustar(ji,jj), 1.E-4_wp )
            zdL = ZCON3*zalpha_w*zdL/(zus_a*zus_a*zus_a)

            !! Stability function Phi_t(-z/L) (zdL is -z/L) :
            zsgn = 0.5_wp + SIGN(0.5_wp, zdL)  ! zdl > 0. => 1.  / zdl < 0. => 0.
            zdL2 = zdL*zdL
            ZPHI =    zsgn      * ( 1._wp + (5._wp*zdL + 4._wp*zdL2)/(1._wp + 3._wp*zdL + 0.25_wp*zdL2) ) &  ! (zdL > 0) Takaya et al.
               & + (zsgn+1._wp) * ( 1._wp/SQRT(1._wp - 16._wp*(-ABS(zdL))) )        ! (zdl < 0) Eq. 8.136
            !
            !! FOR zdL > 0.0, old relations:
            !         ZPHI = 1.+5._wp*zdL                                ! Eq. 8.136 (Large et al. 1994)
            !         ZPHI = 1.+5.0*(zdL+zdL**2)/(1.0+3.0*zdL+zdL**2) ! SHEBA, Grachev et al. 2007

            !! Solving 11 by itteration with time step of zdt...
            ZZ = 1._wp + ZCON4*rdt*zus_w(ji,jj)/ZPHI
            ZZ = SIGN( MAX(ABS(ZZ) , 1e-4_wp), ZZ )
            pdT(ji,jj) = MAX( 0._wp , (ZDT + ZCON5*ZSRD*rdt)/ZZ )

         END DO ! DO ji = 1, jpi
      END DO ! DO jj = 1, jpj

   END SUBROUTINE WL_ECMWF

   !!======================================================================
END MODULE mod_wl_ecmwf
