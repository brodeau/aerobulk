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
   !!       Routine "cs_ecmwf" and "wl_ecmwf" maintained and developed in AeroBulk
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

   !! Warm-layer related arrays:
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC :: &
      &                        dT_wl,     &  !: dT due to warm-layer effect => difference between "almost surface (right below viscous layer, z=delta)
                                !!                                     !: and depth of bulk SST (z=gdept_1d(1))
      &                        Hz_wl         !: depth of warm-layer [m]   ! for now constant for ECMWF...

   REAL(wp), PARAMETER, PUBLIC :: rd0  = 3.    !: Depth scale [m] of warm layer, "d" in Eq.11 (Zeng & Beljaars 2005)
   REAL(wp), PARAMETER         :: zRhoCp_w = rho0_w*rCp0_w
   !
   REAL(wp), PARAMETER         :: rNuwl0 = 0.5  !: Nu (exponent of temperature profile) Eq.11
   !                                            !: (Zeng & Beljaars 2005) !: set to 0.5 instead of
   !                                            !: 0.3 to respect a warming of +3 K in calm
   !                                            !: condition for the insolation peak of +1000W/m^2
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE CS_ECMWF( pQsw, pQnsol, pustar, pSST, pdT_cs )
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
      REAL(wp), INTENT(in) :: pQsw   ! net solar a.k.a shortwave radiation into the ocean (after albedo) [W/m^2]
      REAL(wp), INTENT(in) :: pQnsol ! non-solar heat flux to the ocean [W/m^2]
      REAL(wp), INTENT(in) :: pustar ! friction velocity, temperature and humidity (u*,t*,q*)
      REAL(wp), INTENT(in) :: pSST   ! bulk SST [K]
      REAL(wp), INTENT(out):: pdT_cs !: dT due to cool-skin effect => temperature difference between
      !!                             !: air-sea interface (z=0) and right below viscous layer (z=delta)
      !!---------------------------------------------------------------------
      INTEGER  :: jc
      REAL(wp) :: zQabs, zdelta, zfr
      !!---------------------------------------------------------------------
      zQabs = pQnsol ! first guess of heat flux absorbed within the viscous sublayer of thicknes delta,
      !              !   => we DO not miss a lot assuming 0 solar flux absorbed in the tiny layer of thicknes zdelta...

      zdelta = delta_skin_layer_sclr( alpha_sw(pSST), zQabs, pustar )

      DO jc = 1, 4 ! because implicit in terms of zdelta...
         ! Solar absorption, Eq.(5) Zeng & Beljaars, 2005:
         zfr = MAX( 0.065_wp + 11._wp*zdelta - 6.6E-5_wp/zdelta*(1._wp - EXP(-zdelta/8.E-4_wp)) , 0.01_wp )
         !              =>  (WARNING: 0.065 rather than 0.137 in Fairal et al. 1996)
         zQabs = pQnsol + zfr*pQsw
         zdelta = delta_skin_layer_sclr( alpha_sw(pSST), zQabs, pustar )
      END DO

      pdT_cs = zQabs*zdelta/rk0_w   ! temperature increment, yes dT_cs can actually > 0, if Qabs > 0 (rare but possible!)

   END SUBROUTINE CS_ECMWF


   SUBROUTINE WL_ECMWF( ki, kj, pQsw, pQnsol, pustar, pSST,  pustk )
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
      !!---------------------------------------------------------------------
      INTEGER , INTENT(in) :: ki, kj
      REAL(wp), INTENT(in) :: pQsw     ! surface net solar radiation into the ocean [W/m^2]     => >= 0 !
      REAL(wp), INTENT(in) :: pQnsol   ! surface net non-solar heat flux into the ocean [W/m^2] => normally < 0 !
      REAL(wp), INTENT(in) :: pustar   ! friction velocity [m/s]
      REAL(wp), INTENT(in) :: pSST     ! bulk SST at depth gdept_1d(1) [K]
      !!
      REAL(wp), OPTIONAL, INTENT(in) :: pustk ! surface Stokes velocity [m/s]
      !
      INTEGER :: jc
      !
      REAL(wp) :: &
         & zHwl,    &  !: thickness of the warm-layer [m]
         & ztcorr,  &  !: correction of dT w.r.t measurement depth of bulk SST (first T-point)
         & zalpha, & !: thermal expansion coefficient of sea-water [1/K]
         & zdTwl_b, zdTwl_n, & ! temp. diff. between "almost surface (right below viscous layer) and bottom of WL
         & zfr, zeta, &
         & zusw, zusw2, &
         & zLa, zfLa, &
         & flg, zwf, zQabs, &
         & zA, zB, zL1, zL2, &
         &  zcst0, zcst1, zcst2, zcst3
      !!
      LOGICAL :: l_pustk_known
      !!---------------------------------------------------------------------
      l_pustk_known = ( PRESENT(pustk) )

      zHwl = Hz_wl(ki,kj) ! first guess for warm-layer depth, and final choice because `Hz_wl` fixed in present ECMWF algo (as opposed to COARE)
      !!                  ! it is constant! => rd0 == 3m ! Zeng & Beljaars....

      !! Previous value of dT / warm-layer, adapted to depth:
      flg = 0.5_wp + SIGN( 0.5_wp , gdept_1d(1)-zHwl )               ! => 1 when gdept_1d(1)>zHwl (dT_wl(ki,kj) = zdTwl) | 0 when z_s$
      ztcorr = flg + (1._wp - flg)*gdept_1d(1)/zHwl
      zdTwl_b = MAX ( dT_wl(ki,kj) / ztcorr , 0._wp )
      ! zdTwl is the difference between "almost surface (right below viscous layer) and bottom of WL (here zHwl)
      ! pdT         "                          "                                    and depth of bulk SST (here gdept_1d(1))!
      !! => but of course in general the bulk SST is taken shallower than zHwl !!! So correction less pronounced!
      !! => so here since pdT is difference between surface and gdept_1d(1), need to increase fof zdTwl !

      zalpha = alpha_sw( pSST ) ! thermal expansion coefficient of sea-water (SST accurate enough!)


      ! *** zfr = Fraction of solar radiation absorbed in warm layer (-)
      zfr = 1._wp - 0.28_wp*EXP(-71.5_wp*zHwl) - 0.27_wp*EXP(-2.8_wp*zHwl) - 0.45_wp*EXP(-0.07_wp*zHwl)  !: Eq. 8.157

      zQabs = zfr*pQsw + pQnsol       ! tot heat absorbed in warm layer

      zusw  = MAX( pustar, 1.E-4_wp ) * sq_radrw    ! u* in the water
      zusw2 = zusw*zusw

      ! Langmuir:
      IF ( l_pustk_known ) THEN
         zLa = SQRT(zusw/MAX(pustk,1.E-6))
      ELSE
         zla = 0.3_wp
      END IF
      zfLa = MAX( zla**(-2._wp/3._wp) , 1._wp )   ! Eq.(6)

      zwf = 0.5_wp + SIGN(0.5_wp, zQabs)  ! zQabs > 0. => 1.  / zQabs < 0. => 0.

      zcst1 = vkarmn*grav*zalpha

      ! 1/L when zQabs > 0 :
      zL2 = zcst1*zQabs / (zRhoCp_w*zusw2*zusw)

      zcst2 = zcst1 / ( 5._wp*zHwl*zusw2 )  !OR: zcst2 = zcst1*rNuwl0 / ( 5._wp*zHwl*zusw2 ) ???

      zcst0 = rdt * (rNuwl0 + 1._wp) / zHwl

      zA = zcst0 * zQabs / ( rNuwl0 * zRhoCp_w )

      zcst3 = -zcst0 * vkarmn * zusw * zfLa

      !! Sorry about all these constants ( constant w.r.t zdTwl), it's for
      !! the sake of optimizations... So all these operations are not done
      !! over and over within the iteration loop...

      !! T R U L L Y   I M P L I C I T => needs itteration
      !! => have to itterate just because the 1/(Obukhov length), zL1, uses zdTwl when zQabs < 0..
      !!    (without this term otherwize the implicit analytical solution is straightforward...)
      zdTwl_n = zdTwl_b
      DO jc = 1, 10
         
         zdTwl_n = 0.5_wp * ( zdTwl_n + zdTwl_b ) ! semi implicit, for faster convergence
         
         ! 1/L when zdTwl > 0 .AND. zQabs < 0 :
         zL1 =         SQRT( zdTwl_n * zcst2 ) ! / zusw !!! Or??? => vkarmn * SQRT( zdTwl_n*grav*zalpha/( 5._wp*zHwl ) ) / zusw
         !zL1 = vkarmn*SQRT( zdTwl_n       *grav*zalpha        / ( 5._wp*zHwl ) ) / zusw   ! => vkarmn outside, not inside zcst1 (just for this particular line) ???
         
         ! Stability parameter (z/L):
         zeta =  (1._wp - zwf) * zHwl*zL1   +   zwf * zHwl*zL2
         
         zB = zcst3 / PHI(zeta)

         zdTwl_n = MAX ( zdTwl_b + zA + zB*zdTwl_n , 0._wp )  ! Eq.(6)

      END DO

      !! Update:
      dT_wl(ki,kj) = zdTwl_n * ztcorr

   END SUBROUTINE WL_ECMWF


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
END MODULE mod_skin_ecmwf
