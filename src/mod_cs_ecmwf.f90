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
   !!    Cool-skin parametrization (doc IFS@ECMWF, cy45r1)
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

   !! Cool-skin:
   REAL(wp), PARAMETER :: zroadrw = rho0_a/rho0_w , &         ! Density ratio
      &                   zcon2   = 16._wp * grav * rho0_w * rCp0_w * rnu0_w*rnu0_w*rnu0_w / (rk0_w*rk0_w), &
      &                   rNu0     = 0.5                !: Nu (exponent of temperature profile) Eq.11
   !                                                    !: (Zeng & Beljaars 2005) !: set to 0.5 instead of
   !                                                    !: 0.3 to respect a warming of +3 K in calm
   !                                                    !: condition for the insolation peak of +1000W/m^2
   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE CS_ECMWF( pQsw, pQnsol, pustar, pSST,  pdT )
      !!---------------------------------------------------------------------
      !!
      !!  Cool-Skin scheme according to Fairall et al. 1996
      !!  " A prognostic scheme of sea surface skin temperature for modeling and data assimilation "
      !!     as parameterized in IFS Cy45r1 / E.C.M.W.F.
      !!     ------------------------------------------------------------------
      !!
      !!  **   INPUT:
      !!     *pQsw*       surface net solar radiation into the ocean     [W/m^2] => >= 0 !
      !!     *pQnsol*     surface net non-solar heat flux into the ocean [W/m^2] => normally < 0 !
      !!     *pustar*     friction velocity u*                           [m/s]
      !!     *pSST*       bulk SST at depth z_sst                        [K]
      !!
      !!  **  INPUT/OUTPUT:
      !!     *pdT*  : as input =>  previous estimate of dT cool-skin
      !!              as output =>  new estimate of dT cool-skin
      !!
      !!------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pQsw   ! net solar a.k.a shortwave radiation into the ocean (after albedo) [W/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pQnsol ! non-solar heat flux to the ocean [W/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pustar  ! friction velocity, temperature and humidity (u*,t*,q*)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pSST ! bulk SST [K]
      !!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pdT    ! dT due to cool-skin effect
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) :: zQnet, zQnsol, zlamb, zdelta, zalpha_w, zfr, &
         & zusw, zusw2
      !!---------------------------------------------------------------------

      DO jj = 1, jpj
         DO ji = 1, jpi

            zalpha_w = alpha_sw( pSST(ji,jj) ) ! thermal expansion coefficient of sea-water (SST accurate enough!)

            zQnsol = MAX( 1._wp , - pQnsol(ji,jj) ) ! Non-solar heat loss to the atmosphere

            zusw  = MAX(pustar(ji,jj), 1.E-4_wp)*SQRT(zroadrw)    ! u* in the water
            zusw2 = zusw*zusw

            zlamb = 6._wp*( 1._wp + (zQnsol*zalpha_w*zcon2/(zusw2*zusw2 ))**0.75 )**(-1./3.)   ! w.r.t COARE 3p6 => seems to ommit absorbed zfr*Qsw (Qnet i.o. Qnsol) and effect of evap...
            !                                                                                  ! so zlamb not implicit in terms of zdelta (zfr(delta)), so no need to have last guess of delta as in COARE 3.6...
            zdelta = zlamb*rnu0_w/zusw

            zfr   = MAX( 0.065_wp + 11._wp*zdelta - 6.6E-5_wp/zdelta*(1._wp - EXP(-zdelta/8.E-4_wp)) , 0.01_wp ) ! Solar absorption / Eq. 8.131 / IFS cy40r1, doc, Part IV,
            zQnet = MAX( 1._wp , zQnsol - zfr*pQsw(ji,jj) ) ! Total cooling at the interface

            !! Update!
            pdT(ji,jj) =  MIN( - zQnet*zdelta/rk0_w , 0._wp )   ! temperature increment

         END DO
      END DO

   END SUBROUTINE CS_ECMWF

   !!======================================================================
END MODULE mod_cs_ecmwf
