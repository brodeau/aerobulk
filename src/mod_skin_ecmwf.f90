! AeroBulk / 2016 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_skin_ecmwf
   !!====================================================================================
   !!       Cool-Skin and warm-layer correction of SST
   !!    Cool-skin and warm-layer parametrization (Fairall et al. 1996)
   !!
   !!       Routine "cs_coare3p6" and "wl_coare3p6" maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, 2019
   !!
   !!====================================================================================
   USE mod_const   !: physical and othe constants
   USE mod_phymbl  !: thermodynamics

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: CS_ECMWF, WL_ECMWF

   !! Cool-skin related parameters:
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC :: &
      &                        dT_cs         !: dT due to cool-skin effect => temperature difference between air-sea interface (z=0) and right below viscous layer (z=delta)
   REAL(wp), PARAMETER :: zcon0 = -16._wp * grav * rho0_w * rCp0_w * rnu0_w*rnu0_w*rnu0_w / ( rk0_w*rk0_w ) ! "-" because ocean convention: Qabs > 0 => gain of heat for ocean!
   !!                             => see eq.(14) in Fairall et al. 1996   (eq.(6) of Zeng aand Beljaars is WRONG! (typo?)

   !! Warm-layer related parameters:
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC :: &
      &                        dT_wl         !: dT due to warm-layer effect => difference between "almost surface (right below viscous layer, z=delta) and depth of bulk SST (z=gdept_1d(1))

   REAL(wp), PARAMETER :: rd0  = 3.        !: Depth scale [m] of warm layer, "d" in Eq.11 (Zeng & Beljaars 2005)
   REAL(wp), PARAMETER :: zRhoCp_w = rho0_w*rCp0_w
   REAL(wp), PARAMETER :: rNu0 = 1.0       !:  be closer to COARE3p6 ???!LOLO
   !REAL(wp), PARAMETER :: rNu0 = 0.5       !: Nu (exponent of temperature profile) Eq.11
   !                                       !: (Zeng & Beljaars 2005) !: set to 0.5 instead of
   !                                       !: 0.3 to respect a warming of +3 K in calm
   !                                       !: condition for the insolation peak of +1000W/m^2
   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE CS_ECMWF( pQsw, pQnsol, pustar, pSST )
      !!---------------------------------------------------------------------
      !!
      !! Cool-skin parameterization, based on Fairall et al., 1996:
      !!
      !! Fairall, C. W., Bradley, E. F., Godfrey, J. S., Wick, G. A.,
      !! Edson, J. B., and Young, G. S. ( 1996), Cool‐skin and warm‐layer
      !! effects on sea surface temperature, J. Geophys. Res., 101( C1), 1295-1308,
      !! doi:10.1029/95JC03190.
      !!
      !!------------------------------------------------------------------
      !!
      !!  **   INPUT:
      !!     *pQsw*       surface net solar radiation into the ocean     [W/m^2] => >= 0 !
      !!     *pQnsol*     surface net non-solar heat flux into the ocean [W/m^2] => normally < 0 !
      !!     *pustar*     friction velocity u*                           [m/s]
      !!     *pSST*       bulk SST (taken at depth gdept_1d(1))          [K]
      !!------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pQsw   ! net solar a.k.a shortwave radiation into the ocean (after albedo) [W/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pQnsol ! non-solar heat flux to the ocean [W/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pustar  ! friction velocity, temperature and humidity (u*,t*,q*)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pSST ! bulk SST [K]
      !!---------------------------------------------------------------------
      INTEGER  :: ji, jj, jc
      REAL(wp) :: zQabs, zdelta, zfr
      !!---------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi

            zQabs  = MIN( -0.1_wp , pQnsol(ji,jj) ) ! first guess, we do not miss a lot assuming 0 solar flux absorbed in the tiny layer of thicknes$
            !                                       ! also, we ONLY consider when the viscous layer is loosing heat to the atmosphere, we only deal with cool-skin! => hence the "MIN( -0$

            zdelta = delta_skin_layer( pSST(ji,jj), zQabs, pustar(ji,jj) )

            DO jc = 1, 4 ! because implicit in terms of zdelta...
               zfr    = MAX( 0.065_wp + 11._wp*zdelta - 6.6E-5_wp/zdelta*(1._wp - EXP(-zdelta/8.E-4_wp)) , 0.01_wp ) ! Solar absorption, Eq.(5) Zeng & Beljaars, 2005
               !              =>  (WARNING: 0.065 rather than 0.137 in Fairal et al. 1996)
               zQabs  = MIN( -0.1_wp , pQnsol(ji,jj) + zfr*pQsw(ji,jj) ) ! Total cooling at the interface
               zdelta = delta_skin_layer( pSST(ji,jj), zQabs, pustar(ji,jj) )
            END DO

            dT_cs(ji,jj) = MIN( zQabs*zdelta/rk0_w , 0._wp )   ! temperature increment

         END DO
      END DO

   END SUBROUTINE CS_ECMWF



   SUBROUTINE WL_ECMWF( pQsw, pQnsol, pustar, pSST )
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
      !!     *pSST*       bulk SST  (taken at depth gdept_1d(1))         [K]
      !!------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pQsw     ! surface net solar radiation into the ocean [W/m^2]     => >= 0 !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pQnsol   ! surface net non-solar heat flux into the ocean [W/m^2] => normally < 0 !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pustar   ! friction velocity [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pSST     ! bulk SST at depth gdept_1d(1) [K]
      !
      INTEGER :: ji,jj
      !
      REAL(wp) :: &
         & zdz,    & !: thickness of the warm-layer [m]
         & zalpha_w, & !: thermal expansion coefficient of sea-water
         & ZSRD,    &
         & zdTwl,   & ! temp. diff. between "almost surface (right below viscous layer) and bottom of WL
         & zfr,zdL, ztmp, &
         & zus_a, zusw, zusw2, &
         & flg, zQabs, ZL1, ZL2
      !!---------------------------------------------------------------------


      DO jj = 1, jpj
         DO ji = 1, jpi

            zdz = rd0 ! first guess for warm-layer depth (and unique..., less advanced than COARE3p6 !)

            ! zdTwl is the difference between "almost surface (right below viscous layer) and bottom of WL (here zdz)
            ! pdT         "                          "                                    and depth of bulk SST (here gdept_1d(1))!
            !! => but of course in general the bulk SST is taken shallower than zdz !!! So correction less pronounced!
            !! => so here since pdT is difference between surface and gdept_1d(1), need to increase fof zdTwl !
            flg = 0.5_wp + SIGN( 0.5_wp , gdept_1d(1)-zdz )               ! => 1 when gdept_1d(1)>zdz (dT_wl(ji,jj) = zdTwl) | 0 when z_s$
            zdTwl = dT_wl(ji,jj) / ( flg + (1._wp-flg)*gdept_1d(1)/zdz )
            !PRINT *, 'LOLO/mod_wl_ecmwf.f90: zdTwl2=', zdTwl
            !PRINT *, ''

            zalpha_w = alpha_sw( pSST(ji,jj) ) ! thermal expansion coefficient of sea-water (SST accurate enough!)


            ! *** zfr = Fraction of solar radiation absorbed in warm layer (-)
            zfr = 1._wp - 0.28_wp*EXP(-71.5_wp*zdz) - 0.27_wp*EXP(-2.8_wp*zdz) - 0.45_wp*EXP(-0.07_wp*zdz)  !: Eq. 8.157

            zQabs = zfr*pQsw(ji,jj) + pQnsol(ji,jj)       ! tot heat absorbed in warm layer

            zusw  = MAX(pustar(ji,jj), 1.E-4_wp)*SQRT(roadrw)    ! u* in the water
            zusw2 = zusw*zusw


            !! *** 1st rhs term in eq. 8.156 (IFS doc Cy45r1):
            ZL1 = zQabs / ( zdz * zRhoCp_w * rNu0 ) * (rNu0 + 1._wp)


            !! Buoyancy flux and stability parameter (zdl = -z/L) in water
            ZSRD = zQabs/zRhoCp_w
            !
            flg = 0.5_wp + SIGN(0.5_wp, ZSRD)  ! ZSRD > 0. => 1.  / ZSRD < 0. => 0.
            ztmp = MAX(zdTwl,0._wp)
            zdl = (1.-flg) * ( zusw2 * SQRT(ztmp/(5._wp*zdz*grav*zalpha_w/rNu0)) ) & ! (zdTwl > 0.0 .AND. ZSRD < 0.0)
               & +  flg    *  ZSRD                                                                  !   otherwize
            !
            zus_a = MAX( pustar(ji,jj), 1.E-4_wp )
            zdL = zdz*vkarmn*grav/(roadrw)**1.5_wp*zalpha_w*zdL/(zus_a*zus_a*zus_a)

            !! *** 2nd rhs term in eq. 8.156 (IFS doc Cy45r1):
            ZL2 = - (rNu0 + 1._wp) * vkarmn * zusw / ( zdz * PHI(zdl) )

            ! Forward time / explicit solving of eq. 8.156 (IFS doc Cy45r1): (f_n+1 == dT_wl(ji,jj) ; f_n == zdTwl)
            zdTwl = MAX ( zdTwl + rdt*ZL1 + rdt*ZL2*zdTwl , 0._wp )

            ! zdTwl is the difference between "almost surface (right below viscous layer) and bottom of WL (here zdz)
            !! => but of course in general the bulk SST is taken shallower than zdz !!! So correction less pronounced!

            flg = 0.5_wp + SIGN( 0.5_wp , gdept_1d(1)-zdz )               ! => 1 when gdept_1d(1)>zdz (dT_wl(ji,jj) = zdTwl) | 0 when z_s$
            dT_wl(ji,jj) = zdTwl * ( flg + (1._wp-flg)*gdept_1d(1)/zdz )

         END DO
      END DO

   END SUBROUTINE WL_ECMWF



   FUNCTION delta_skin_layer( pSST, pQabs, pustar_a )
      !!---------------------------------------------------------------------
      !! Computes the thickness (m) of the viscous skin layer.
      !! Based on Fairall et al., 1996
      !!
      !! Fairall, C. W., Bradley, E. F., Godfrey, J. S., Wick, G. A.,
      !! Edson, J. B., and Young, G. S. ( 1996), Cool‐skin and warm‐layer
      !! effects on sea surface temperature, J. Geophys. Res., 101( C1), 1295-1308,
      !! doi:10.1029/95JC03190.
      !!
      !! L. Brodeau, october 2019
      !!---------------------------------------------------------------------
      REAL(wp)                :: delta_skin_layer
      REAL(wp), INTENT(in)    :: pSST     ! bulk SST [K] => to know the thermal expansion [K]
      REAL(wp), INTENT(in)    :: pQabs    ! < 0 !!! part of the net heat flux actually absorbed in the WL [W/m^2] => term "Q + Rs*fs" in eq.6 of Fairall et al. 1996
      REAL(wp), INTENT(in)    :: pustar_a ! friction velocity in the air (u*) [m/s]
      !!---------------------------------------------------------------------
      REAL(wp) :: zusw, zusw2, zlamb, zalpha_w, zQb !, ztf, zQ
      !!---------------------------------------------------------------------

      zQb = pQabs

      !zQ = MIN( -0.1_wp , pQabs )

      !ztf = 0.5_wp + SIGN(0.5_wp, zQ)  ! Qabs < 0 => cooling of the layer => ztf = 0 (normal case)
      !                                   ! Qabs > 0 => warming of the layer => ztf = 1 (ex: weak evaporation and strong positive sensible heat flux)
      zalpha_w = alpha_sw( pSST ) ! thermal expansion coefficient of sea-water (SST accurate enough!)

      zusw  = MAX(pustar_a, 1.E-4_wp) * sq_radrw    ! u* in the water
      zusw2 = zusw*zusw

      zlamb = 6._wp*( 1._wp + (zalpha_w*zcon0/(zusw2*zusw2)*zQb)**0.75 )**(-1./3.) ! see eq.(14) in Fairall et al., 1996
      !zlamb = 6._wp*( 1._wp + MAX(zalpha_w*zcon0/(zusw2*zusw2)*zQ, 0._wp)**0.75 )**(-1./3.) ! see eq.(14) in Fairall et al., 1996

      delta_skin_layer = zlamb*rnu0_w/zusw

      !delta_skin_layer =  (1._wp - ztf) * zlamb*rnu0_w/zusw    &         ! see eq.(12) in Fairall et al., 1996
      !   &               +     ztf  * MIN(6._wp*rnu0_w/zusw , 0.007_wp)
   END FUNCTION delta_skin_layer



   FUNCTION PHI( pzeta)
      !!---------------------------------------------------------------------
      !!
      !! Zeng & Beljaars
      !!
      !! L. Brodeau, october 2019
      !!---------------------------------------------------------------------
      REAL(wp)                :: PHI
      REAL(wp), INTENT(in)    :: pzeta    ! stability parameter
      !!---------------------------------------------------------------------
      REAL(wp) :: ztf, zzt2
      !!---------------------------------------------------------------------
      !
!!! Stability function PHI_t(-z/L) (zdL is -z/L) :
      !flg = 0.5_wp + SIGN(0.5_wp, zdL)  ! zdl > 0. => 1.  / zdl < 0. => 0.
      !zdL2 = zdL*zdL
      !ZPHI =    flg   * ( 1._wp + (5._wp*zdL + 4._wp*zdL2)/(1._wp + 3._wp*zdL + 0.25_wp*zdL2) ) &  ! (zdL > 0)
      !   & + (1.-flg) * ( 1._wp/SQRT(1._wp - 16._wp*(-ABS(zdL))) )                                 ! (zdl < 0)
!!! FOR zdL > 0.0, old relations:
      !!         ZPHI = 1.+5._wp*zdL                                ! Eq. 8.136 (Large et al. 1994)
      !!         ZPHI = 1.+5.0*(zdL+zdL**2)/(1.0+3.0*zdL+zdL**2) ! SHEBA, Grachev et al. 2007
      !!
      zzt2 = pzeta*pzeta
      !
      ztf = 0.5_wp + SIGN(0.5_wp, pzeta)  ! zeta > 0 => ztf = 1
      !                                   ! zeta < 0 => ztf = 0
      PHI =      ztf     * ( 1. + (5.*pzeta + 4.*zzt2)/(1. + 3.*pzeta + 0.25*zzt2) ) &   ! zeta > 0 Takaya et al.
         &  + (1. - ztf) * 1./SQRT( 1. - 16.*(-ABS(pzeta)) )                             ! zeta < 0 Eq.(8.136)
      !
   END FUNCTION PHI

   !!======================================================================
END MODULE mod_skin_ecmwf
