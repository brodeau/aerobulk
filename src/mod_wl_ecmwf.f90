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
   !!       Warm-Layer correction of SST (if needed)
   !!
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (2m) to zu (10m) if needed
   !!   * the effective bulk wind speed at 10m U_blk
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!    Using the bulk formulation/param. of ECMWF
   !!      + consideration of cool-skin warm layer parametrization (CS: Fairall et al. 1996; WL: Zeng & Beljaars, 2005 )
   !!
   !!       Routine "cswl_ecmwf" maintained and developed in AeroBulk
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

   !! Cool-Skin / Warm-Layer related parameters:
   REAL(wp), PARAMETER :: rd0  = 3.        !: Depth scale [m], "d" in Eq.11 (Zeng & Beljaars 2005)
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


   SUBROUTINE WL_ECMWF( pQsw, pQnsol, pustar, pSST, pTs )
      !!---------------------------------------------------------------------
      !!
      !!  Cool-Skin Warm-Layer scheme according to Zeng & Beljaars, 2005 (GRL)
      !!  " A prognostic scheme of sea surface skin temperature for modeling and data assimilation "
      !!
      !!    As included in IFS Cy40   /  E.C.M.W.F.
      !!     ------------------------------------------------------------------
      !!
      !!  **   INPUT:
      !!
      !!     *pQsw*       net solar radiative flux to the ocean
      !!     *pQnsol*     net non-solar heat flux to the ocean
      !!     *pustar*     friction velocity u*
      !!     *pSST*       SST
      !!
      !!   **  INPUT/OUTPUT:
      !!     *pTs*  : as input  =>  previous estimate of skin temperature
      !!             as output =>  new estimate of skin temperature
      !!
      !!------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pQsw     ! net solar radiation into the sea [W/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pQnsol   ! net non-solar heat flux into the sea [W/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pustar
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pSST     ! bulk SST
      !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pTs
      !
      INTEGER :: ji,jj
      !
      REAL(wp) :: &
         & zRhoCp_w, &
         & ZCON3,ZCON4,ZCON5, zQnsol ,zQnet, zlamb, zdelta,&
         & ZSRD,ZDSST,ZZ,ZEPDU2,&
         & ZFI,zdL,zdL2, ztmp, &
         & ZPHI,ZROADRW, &
         & zus_a, &
         & zfs, zsgn
      !
      REAL(wp), DIMENSION(jpi,jpj) :: &
         &  zalpha_w, &       !: thermal expansion coefficient of seawater
         & zus_w, zus_w2, &   !: u* and u*^2 in water
         &       zdT_w        !: warm skin temperature increment
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
      !    ZFI = Fraction of solar radiation absorbed in warm layer (-)
      ZFI = 1._wp - 0.28_wp*EXP(-71.5_wp*rd0) - 0.27_wp*EXP(-2.8_wp*rd0) - 0.45_wp*EXP(-0.07_wp*rd0)  !: Eq. 8.157
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
      ! Ocean buoyancy
      zalpha_w(:,:) = MAX( 1.E-5_wp , 1.E-5_wp*(pTs(:,:) - rt0) ) ! thermal expansion coefficient of water
      !
      zdT_w = 0._wp
      !
      !  3. Cool skin (Fairall et al. 1996)
      !------------------------------------
      !DO jj = 1, jpj
      !   DO ji = 1, jpi
      !      !
      !      ! Non-solar heat loss to the atmosphere:
      !      zQnsol = MAX( 1._wp , - pQnsol(ji,jj) )!
      !
      !      zlamb = 6._wp*(1._wp + (zQnsol*zalpha_w(ji,jj)*ZCON2/(zus_w2(ji,jj)*zus_w2(ji,jj)))**0.75)**(-1._wp/3._wp)
      !
      !      zdelta = zlamb*rnu0_w/zus_w(ji,jj)
      !
      !      !    Solar absorption
      !      zfs   = 0.065_wp + 11._wp*zdelta - (6.6E-5_wp/zdelta)*(1._wp - EXP(-zdelta/8.E-4_wp))  ! Eq. 8.131 / IFS cy40r1, doc, Part IV,
      !      zfs   = MAX(zfs , 0.01_wp)
      !      zQnet = MAX( 1._wp , -zfs*pQsw(ji,jj) + zQnsol )
      !      zdT_c(ji,jj) = -zdelta*zQnet/rk0_w
      !
      !   END DO
      !END DO

      ! 2.2 Warm layer; formulation C (Xubin Zeng)
      !!--------------------------------------------

      DO jj = 1, jpj
         DO ji = 1, jpi

            ZDSST = pTs(ji,jj) - pSST(ji,jj)

            !! Buoyancy flux and stability parameter (zdl = -z/L) in water
            !
            !! Qt/(rho_w*Cpw):
            ZSRD = ( pQsw(ji,jj)*ZFI + pQnsol(ji,jj) )/zRhoCp_w
            !
            zsgn = 0.5_wp + SIGN(0.5_wp, ZSRD)  ! ZSRD > 0. => 1.  / ZSRD < 0. => 0.
            ztmp = MAX(ZDSST,0._wp)
            zdl = (zsgn+1._wp) * ( zus_w2(ji,jj) * SQRT(ztmp/(5._wp*rd0*grav*zalpha_w(ji,jj)/rNu0)) ) & ! (ZDSST > 0.0 .AND. ZSRD < 0.0)
               &  +    zsgn    *  ZSRD                                                                  !   otherwize
            !
            zus_a = MAX( pustar(ji,jj), 1.E-4_wp )
            zdL = ZCON3*zalpha_w(ji,jj)*zdL/(zus_a*zus_a*zus_a)

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
            zdT_w(ji,jj) = MAX( 0._wp , (ZDSST + ZCON5*ZSRD*rdt)/ZZ )

         END DO ! DO ji = 1, jpi
      END DO ! DO jj = 1, jpj

      ! 3. Apply warm layer and cool skin effects
      !------------------------------------------
      pTs(:,:) = pSST(:,:) + zdT_w(:,:)


   END SUBROUTINE WL_ECMWF

   !!======================================================================
END MODULE mod_wl_ecmwf
