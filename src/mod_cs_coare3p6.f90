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

   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE CS_COARE3P6( pTzu, pqzu, pSST, pslp, pUzu, pus, pts, pqs, &
      &                    pQnsol, pQsw, pQlat, pdelta,  pdT )
      !!
      !!  **   OUTPUT:
      !!     *pdT*        dT due to warming at depth of pSST such that SST_actual = pSST + pdT
      !!---------------------------------------------------------------------
      !!
      !!  Cool-Skin Warm-Layer scheme according to Fairall et al. 1996
      !!     
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pTzu ! air temperature at height zu above sea surface [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pqzu ! air specific humidity at height zu above sea surface [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pSST ! bulk SST [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pslp ! sea-level atmospheric pressure [Pa]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pUzu ! scalar wind speed at height zu above sea surface [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pus, pts, pqs ! friction velocity, temperature and humidity (u*,t*,q*)
      !!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pQnsol ! non-solar heat flux to the ocean [W/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pQsw   ! net solar a.k.a shortwave radiation into the ocean (after albedo) [W/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pQlat  ! latent heat flux [W/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pdelta ! thickness of the viscous (skin) layer [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pdT    ! dT due to cooling such as SSST = pSST + pdT
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) :: zz0, zz1, zz2, zus, zfr, zrho_a, &
         &        T_s, q_s, zdt, zdq, ztf, &
         &        zdelta, zlamb, zalpha, zQt
      !!---------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi

            zdelta = pdelta(ji,jj)
            
            T_s = pSST(ji,jj) + pdT(ji,jj) ! actual skin temperature
            zz0 = e_sat(MAX(T_s, 200._wp)) !
            q_s = rdct_qsat_salt*reps0*zz0/(pslp(ji,jj) - (1. - reps0)*zz0) ! actual specific hum. at saturation at T=T_s

            zdt = pTzu(ji,jj) - T_s  ; zdt = SIGN( MAX(ABS(zdt),1e-6_wp), zdt )
            zdq = pqzu(ji,jj) - q_s  ; zdq = SIGN( MAX(ABS(zdq),1e-9_wp), zdq )

            
            !! Fraction of the shortwave flux absorbed by the cool-skin sublayer:
            !zQsw_f = 0.065 + 11.*zdelta - 6.6e-5/zdelta*(1. - EXP(-zdelta/8.e-4)) ! Eq.16 (Fairall al. 1996b)
            zfr = 0.137 + 11.*zdelta - 6.6e-5/zdelta*(1. - EXP(-zdelta/8.e-4)) ! Eq.16 (Fairall al. 1996b)
            !LB: why 0.065 and not 0.137 like in the paper??? Beljaars & Zeng use 0.065
            !LB: maybe comes from Wick et al 2005 ...
            
            zQt = -(pQnsol(ji,jj) + zfr*pQsw(ji,jj))  ! Total cooling at the interface

            ztf = 0.5 + SIGN(0.5_wp, zQt) ! Qt > 0 => cooling of the layer => ztf = 1
            !                               Qt < 0 => warming of the layer => ztf = 0
            
            zalpha = 2.1e-5*MAX(T_s-rt0 + 3.2_wp, 0._wp)**0.79  ! alpha = thermal expansion of water (~2.5E-4) LB: remove from loop, sst accurate enough!

            !! Term alpha*Qb (Qb is the virtual surface cooling inc. buoyancy effect of salinity due to evap):
            zz1 = zalpha*zQt - 0.026*pQlat(ji,jj)*rCp0_w/rLevap  ! alpha*(Eq.8) == alpha*Qb "-" because Qlat < 0
            !! LB: this terms only makes sense if > 0 i.e. in the cooling case
            !! so similar to what's done in ECMWF:
            zz1 = MAX(0._wp , zz1)    ! 1. instead of 0.1 though ZQ = MAX(1.0,-pQlw(ji,jj) - pQsen(ji,jj) - pQlat(ji,jj))

            !! Laurent: too low wind (u*) might cause problem in stable cases:
            zus = MAX(pus(ji,jj), 1.E-4_wp)

            ! Lambda (=> zz0, empirical coeff.) (Eq.14):
            zz0 = 16._wp * zz1 * grav * rho0_w * rCp0_w * rnu0_w*rnu0_w*rnu0_w  ! (numerateur) zz1 == alpha*Q
            zrho_a = rho_air( pTzu(ji,jj), pqzu(ji,jj), pslp(ji,jj) )
            zz2 = zus*zus * zrho_a / rho0_w * rk0_w
            zz2 = zz2*zz2                                             ! denominateur
            !LB:  zz0 has the sign of zz1 and therefore of Qb !
            zlamb =  6._wp*( 1._wp + (zz0/zz2)**(3./4.) )**(-1./3.) !  Eq.14   (Saunders)
            
            ! Updating molecular sublayer thickness (delta):
            zz2    = rnu0_w/(SQRT(zrho_a/rho0_w)*zus)
            zdelta =      ztf    *          zlamb*zz2   &  ! Eq.12 (when alpha*Qb>0 / cooling of layer)
               &    + (1._wp - ztf) * MIN(0.007_wp , 6._wp*zz2 )    ! Eq.12 (when alpha*Qb<0 / warming of layer)
            !LB: changed 0.01 to 0.007
            pdelta(ji,jj) = zdelta

            ! Updating temperature increment:
            pdT(ji,jj) =  MIN( - zQt*zdelta/rk0_w , 0._wp )   ! temperature increment !  Eq.13 Cool skin !LOLO get rid of warming that comes from I don't know which term...
            !
         END DO
      END DO
      
   END SUBROUTINE CS_COARE3P6
   
   !!======================================================================
END MODULE mod_cs_coare3p6
