! AeroBulk / 2016 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_skin_new
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

   PUBLIC :: CS_NEW, WL_NEW

   !! Cool-skin related parameters:
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC :: &
      &                        dT_cs         !: dT due to cool-skin effect => temperature difference between air-sea interface (z=0) and right below viscous layer (z=delta)
   REAL(wp), PARAMETER :: zcon0 = -16._wp * grav * rho0_w * rCp0_w * rnu0_w*rnu0_w*rnu0_w / ( rk0_w*rk0_w ) ! "-" because ocean convention: Qabs > 0 => gain of heat for ocean!
   !!                             => see eq.(14) in Fairall et al. 1996   (eq.(6) of Zeng aand Beljaars is WRONG! (typo?)

   !! Warm-layer related parameters:
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC :: &
      &                        dT_wl,  &   !: dT due to warm-layer effect => difference between "almost surface (right below viscous layer, z=delta) and depth of bulk SST (z=gdept_1d(1))
      &                        Hz_wl       !: depth of warm-layer [m]
   !
   REAL(wp), PARAMETER, PUBLIC :: rd0  = 3.    !: Depth scale [m] of warm layer, "d" in Eq.11 (Zeng & Beljaars 2005)
   REAL(wp), PARAMETER         :: zRhoCp_w = rho0_w*rCp0_w
   !
   REAL(wp), PARAMETER         :: rNuwl0 = 0.5  !: Nu (exponent of temperature profile) Eq.11
   !                                            !: (Zeng & Beljaars 2005) !: set to 0.5 instead of
   !                                            !: 0.3 to respect a warming of +3 K in calm
   !                                            !: condition for the insolation peak of +1000W/m^2
   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE CS_NEW( pQsw, pQnsol, pustar, pSST )
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
            !zQabs  = pQnsol(ji,jj)

            zdelta = delta_skin_layer( pSST(ji,jj), zQabs, pustar(ji,jj) )

            DO jc = 1, 4 ! because implicit in terms of zdelta...
               zfr    = MAX( 0.065_wp + 11._wp*zdelta - 6.6E-5_wp/zdelta*(1._wp - EXP(-zdelta/8.E-4_wp)) , 0.01_wp ) ! Solar absorption, Eq.(5) Zeng & Beljaars, 2005
               !              =>  (WARNING: 0.065 rather than 0.137 in Fairal et al. 1996)
               zQabs  = MIN( -0.1_wp , pQnsol(ji,jj) + zfr*pQsw(ji,jj) ) ! Total cooling at the interface
               !zQabs = pQnsol(ji,jj) + zfr*pQsw(ji,jj)
               zdelta = delta_skin_layer( pSST(ji,jj), zQabs, pustar(ji,jj) )
            END DO

            dT_cs(ji,jj) = MIN( zQabs*zdelta/rk0_w , 0._wp )   ! temperature increment

         END DO
      END DO

   END SUBROUTINE CS_NEW



   SUBROUTINE WL_NEW( pQsw, pQnsol, pustar, pSST,  pustk )
      !!---------------------------------------------------------------------
      !!
      !!  Warm-Layer scheme according to Zeng & Beljaars, 2005 (GRL)
      !!  " A prognostic scheme of sea surface skin temperature for modeling and data assimilation "
      !!
      !!  STIL NO PROGNOSTIC EQUATION FOR THE DEPTH OF THE WARM-LAYER!
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
      !!
      REAL(wp), DIMENSION(jpi,jpj), OPTIONAL, INTENT(in) :: pustk ! surface Stokes velocity [m/s]
      !
      INTEGER :: ji, jj, jc
      !
      REAL(wp) :: &
         & zHwl,    &  !: thickness of the warm-layer [m]
         & ztcorr,  &  !: correction of dT w.r.t measurement depth of bulk SST (first T-point)
         & zalpha_w, & !: thermal expansion coefficient of sea-water [1/K]
         & zdTwl_b, zdTwl_n, & ! temp. diff. between "almost surface (right below viscous layer) and bottom of WL
         & zfr, zeta, ztmp, &
         & zusw, zusw2, &
         & zLa, zfLa, &
         & flg, zwf, zQabs, &
         & ZA, ZB, zL1, zL2, &
         &  zcst0, zcst1, zcst2, zcst3
      !
      LOGICAL :: l_pustk_known
      !!---------------------------------------------------------------------

      l_pustk_known = .FALSE.
      IF ( PRESENT(pustk) ) l_pustk_known = .TRUE.

      DO jj = 1, jpj
         DO ji = 1, jpi

            zHwl = Hz_wl(ji,jj) ! first guess for warm-layer depth (and unique..., less advanced than COARE3p6 !)
            ! it is = rd0 (3m) in default Zeng & Beljaars case...

            !! Previous value of dT / warm-layer, adapted to depth:
            flg = 0.5_wp + SIGN( 0.5_wp , gdept_1d(1)-zHwl )               ! => 1 when gdept_1d(1)>zHwl (dT_wl(ji,jj) = zdTwl) | 0 when z_s$
            ztcorr = flg + (1._wp - flg)*gdept_1d(1)/zHwl
            zdTwl_b = MAX ( dT_wl(ji,jj) / ztcorr , 0._wp )
            ! zdTwl is the difference between "almost surface (right below viscous layer) and bottom of WL (here zHwl)
            ! pdT         "                          "                                    and depth of bulk SST (here gdept_1d(1))!
            !! => but of course in general the bulk SST is taken shallower than zHwl !!! So correction less pronounced!
            !! => so here since pdT is difference between surface and gdept_1d(1), need to increase fof zdTwl !

            zalpha_w = alpha_sw( pSST(ji,jj) ) ! thermal expansion coefficient of sea-water (SST accurate enough!)


            ! *** zfr = Fraction of solar radiation absorbed in warm layer (-)
            zfr = 1._wp - 0.28_wp*EXP(-71.5_wp*zHwl) - 0.27_wp*EXP(-2.8_wp*zHwl) - 0.45_wp*EXP(-0.07_wp*zHwl)  !: Eq. 8.157

            zQabs = zfr*pQsw(ji,jj) + pQnsol(ji,jj)       ! tot heat absorbed in warm layer

            zusw  = MAX( pustar(ji,jj), 1.E-4_wp ) * sq_radrw    ! u* in the water
            zusw2 = zusw*zusw

            ! Langmuir:
            IF ( l_pustk_known ) THEN
               zLa = SQRT(zusw/MAX(pustk(ji,jj),1.E-6))
            ELSE
               zla = 0.3_wp
            END IF
            zfLa = MAX( zla**(-2._wp/3._wp) , 1._wp )   ! Eq.(6)

            zwf = 0.5_wp + SIGN(0.5_wp, zQabs)  ! zQabs > 0. => 1.  / zQabs < 0. => 0.

            zcst1 = vkarmn*grav*zalpha_w

            ! 1/L when zQabs > 0 :
            zL2 = zcst1*zQabs / (zRhoCp_w*zusw2*zusw)
               
            zcst2 = zcst1 / ( 5._wp*zHwl*zusw2 )  !OR: zcst2 = zcst1*rNuwl0 / ( 5._wp*zHwl*zusw2 ) ???

            zcst0 = rdt * (rNuwl0 + 1._wp) / zHwl
            
            ZA = zcst0 * zQabs / ( rNuwl0 * zRhoCp_w )

            zcst3 = -zcst0 * vkarmn * zusw * zfLa

            !! Sorry about all these constants ( constant w.r.t zdTwl), it's for
            !! the sake of optimizations... So all these operations are not done
            !! over and over within the iteration loop...
            
            !! T R U L L Y   I M P L I C I T => needs itteration
            !! => have to itterate just because the 1/(Monin-Obukhov length), zL1, uses zdTwl when zQabs < 0..
            !!    (without this term otherwize the implicit analytical solution is straightforward...)
            zdTwl_n = zdTwl_b
            DO jc = 1, 10
               
               zdTwl_n = 0.5_wp * ( zdTwl_n + zdTwl_b ) ! semi implicit, for faster convergence
               
               ! 1/L when zdTwl > 0 .AND. zQabs < 0 :
               zL1 =         SQRT( zdTwl_n * zcst2 ) ! / zusw !!! Or??? => vkarmn * SQRT( zdTwl_n*grav*zalpha_w/( 5._wp*zHwl ) ) / zusw
               !zL1 = vkarmn*SQRT( zdTwl_n       *grav*zalpha_w        / ( 5._wp*zHwl ) ) / zusw   ! => vkarmn outside, not inside zcst1 (just for this particular line) ???
               
               ! Stability parameter (z/L):
               zeta =  (1._wp - zwf) * zHwl*zL1   +   zwf * zHwl*zL2

               ZB = zcst3 / PHI(zeta)

               zdTwl_n = zdTwl_b + ZA + ZB*zdTwl_n ! Eq.(6)

            END DO
            
            !! Update:
            dT_wl(ji,jj) = MAX ( zdTwl_n , 0._wp ) * ztcorr
            
         END DO
      END DO

   END SUBROUTINE WL_NEW



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
      !! Takaya et al., 2010
      !!  Eq.(5)
      !! L. Brodeau, october 2019
      !!---------------------------------------------------------------------
      REAL(wp)                :: PHI
      REAL(wp), INTENT(in)    :: pzeta    ! stability parameter
      !!---------------------------------------------------------------------
      REAL(wp) :: ztf, zzt2
      !!---------------------------------------------------------------------
      !
      zzt2 = pzeta*pzeta
      !
      ztf = 0.5_wp + SIGN(0.5_wp, pzeta)  ! zeta > 0 => ztf = 1
      !                                   ! zeta < 0 => ztf = 0
      PHI =      ztf     * ( 1. + (5.*pzeta + 4.*zzt2)/(1. + 3.*pzeta + 0.25*zzt2) ) &   ! zeta > 0
         &  + (1. - ztf) * 1./SQRT( 1. - 16.*(-ABS(pzeta)) )                             ! zeta < 0
      !
   END FUNCTION PHI

   !!======================================================================
END MODULE mod_skin_new
