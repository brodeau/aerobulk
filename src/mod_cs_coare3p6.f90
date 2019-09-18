! AeroBulk / 2016 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_cs_coare3p6
   !!====================================================================================
   !!       Cool-Skin correction of SST
   !!    Cool-skin parametrization (Fairall et al. 1996)
   !!
   !!       Routine "cs_coare3p6" maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, 2019
   !!
   !!====================================================================================
   USE mod_const   !: physical and othe constants
   USE mod_phymbl  !: thermodynamics

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: CS_COARE3P6

   !! Cool-skin:
   REAL(wp), PARAMETER :: zroadrw = rho0_a/rho0_w , &         ! Density ratio
      &                   zcon2   = 16._wp * grav * rho0_w * rCp0_w * rnu0_w*rnu0_w*rnu0_w / (rk0_w*rk0_w)
   
   !!----------------------------------------------------------------------
CONTAINS

   
   SUBROUTINE CS_COARE3P6( pQsw, pQnsol, pustar, pSST, pQlat, pdelta,  pdT )
      !!
      !!  **   OUTPUT:
      !!     *pdT*        dT due to warming at depth of pSST such that SST_actual = pSST + pdT
      !!---------------------------------------------------------------------
      !!
      !!  Cool-Skin Warm-Layer scheme according to Fairall et al. 1996
      !!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pQsw   ! net solar a.k.a shortwave radiation into the ocean (after albedo) [W/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pQnsol ! non-solar heat flux to the ocean [W/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pustar  ! friction velocity, temperature and humidity (u*,t*,q*)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pSST ! bulk SST [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pQlat  ! latent heat flux [W/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pdelta ! thickness of the viscous (skin) layer [m]
      !!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pdT    ! dT due to cooling such as SSST = pSST + pdT
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) :: zz1, zz2, zus, zfr, &
         &        zdt, zdq, ztf, zQnsol, &
         &        zdelta, zlamb, zalpha_w, zQnet
      !!---------------------------------------------------------------------

      DO jj = 1, jpj
         DO ji = 1, jpi

            zalpha_w = 2.1e-5*MAX(pSST(ji,jj)-rt0 + 3.2_wp, 0._wp)**0.79  ! thermal expansion coefficient of sea-water (SST accurate enough!) ! alpha = thermal expansion of water (~2.5E-4) LB: remove from loop, sst accurate enough!

            zQnsol = MAX( 1._wp , - pQnsol(ji,jj) ) ! Non-solar heat loss to the atmosphere
            
            zdelta = pdelta(ji,jj)
            
            !! Fraction of the shortwave flux absorbed by the cool-skin sublayer:
            !zQsw_f = 0.065 + 11.*zdelta - 6.6e-5/zdelta*(1. - EXP(-zdelta/8.e-4)) ! Eq.16 (Fairall al. 1996b)
            zfr = 0.137 + 11.*zdelta - 6.6e-5/zdelta*(1. - EXP(-zdelta/8.e-4)) ! Eq.16 (Fairall al. 1996b)
            !LB: why 0.065 and not 0.137 like in the paper??? Beljaars & Zeng use 0.065
            !LB: maybe comes from Wick et al 2005 ...
            
            zQnet = MAX( 1._wp , zQnsol - zfr*pQsw(ji,jj) ) ! Total cooling at the interface

            ztf = 0.5 + SIGN(0.5_wp, zQnet) ! Qt > 0 => cooling of the layer => ztf = 1
            !                               Qt < 0 => warming of the layer => ztf = 0
            

            !! Term alpha*Qb (Qb is the virtual surface cooling inc. buoyancy effect of salinity due to evap):
            zz1 = zalpha_w*zQnet - 0.026*pQlat(ji,jj)*rCp0_w/rLevap  ! alpha*(Eq.8) == alpha*Qb "-" because Qlat < 0
            !! LB: this terms only makes sense if > 0 i.e. in the cooling case
            !! so similar to what's done in ECMWF:
            zz1 = MAX(0._wp , zz1)    ! 1. instead of 0.1 though ZQ = MAX(1.0,-pQlw(ji,jj) - pQsen(ji,jj) - pQlat(ji,jj))

            zus = MAX(pustar(ji,jj), 1.E-4_wp) ! Laurent: too low wind (u*) might cause problem in stable cases:            
            zz2 = zus*zus * zroadrw
            zz2 = zz2*zz2
            zlamb =  6._wp*( 1._wp + (zcon2*zz1/zz2)**0.75 )**(-1./3.) ! Lambda (Eq.14) (Saunders)
            
            ! Updating molecular sublayer thickness (delta):
            zz2    = rnu0_w/(SQRT(zroadrw)*zus)
            zdelta =      ztf    *          zlamb*zz2   &  ! Eq.12 (when alpha*Qb>0 / cooling of layer)
               &    + (1._wp - ztf) * MIN(0.007_wp , 6._wp*zz2 )    ! Eq.12 (when alpha*Qb<0 / warming of layer)
            !LB: changed 0.01 to 0.007

            ! Updating temperature increment:
            pdT(ji,jj) =  MIN( - zQnet*zdelta/rk0_w , 0._wp )   ! temperature increment !  Eq.13 Cool skin !LOLO get rid of warming that comes from I don't know which term...
            
            pdelta(ji,jj) = zdelta

         END DO
      END DO
      
   END SUBROUTINE CS_COARE3P6
   
   !!======================================================================
END MODULE mod_cs_coare3p6
