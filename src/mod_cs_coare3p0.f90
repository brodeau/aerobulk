! AeroBulk / 2016 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_cs_coare3p0
   !!====================================================================================
   !!       Cool-Skin correction of SST
   !!    Cool-skin parametrization (Fairall et al. 1996)
   !!
   !!       Routine CS_COARE3P0 maintained and developed in AeroBulk
   !!             (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, 2016
   !!
   !!====================================================================================
   USE mod_const   !: physical and othe constants
   USE mod_phymbl  !: thermodynamics

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: CS_COARE3P0
   
   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE CS_COARE3P0( pTzu, pqzu, psst, pslp, pU10, pus, pts, pqs, &
      &                   prhoa, pRlw, pQsw, pdelta, pT_s, pq_s )
      !!---------------------------------------------------------------------
      !!
      !!  Cool-Skin Warm-Layer scheme according to Fairall et al. 1996
      !!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pTzu, pqzu, psst, pslp, &
         &                                           pU10, pus, pts, pqs
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: prhoa, pRlw, pQsw
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pdelta, pT_s, pq_s
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) :: zz0, zz1, zz2, zCe, zCh, zus, zQsen, zQlat, zQlw, zfr, &
         &        zdt, zdq, ztf, &
         &        zdelta, zlamb, zalpha, zQt
      !!---------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi

            zdelta = pdelta(ji,jj)

            zdt = pTzu(ji,jj) - pT_s(ji,jj)  ; zdt = SIGN( MAX(ABS(zdt),1e-6_wp), zdt )
            zdq = pqzu(ji,jj) - pq_s(ji,jj)  ; zdq = SIGN( MAX(ABS(zdq),1e-9_wp), zdq )

            !! compute transfer coefficients at zu :
            zz0 = pus(ji,jj)/pU10(ji,jj)    !LB! not needed inside the loop !
            zCh = zz0*pts(ji,jj)/zdt
            zCe = zz0*pqs(ji,jj)/zdq

            ! Turbulent heat fluxes:
            zz1 = prhoa(ji,jj)*pU10(ji,jj)
            zQlat = MIN( rLevap*zCe*zz1*(pqzu(ji,jj) - pq_s(ji,jj)) , 0._wp )
            zQsen =     rCp0_a*zCh*zz1*(pTzu(ji,jj) - pT_s(ji,jj))

            ! Net longwave flux:
            zz1  = pT_s(ji,jj)*pT_s(ji,jj)
            zQlw = emiss_w*(pRlw(ji,jj) - stefan*zz1*zz1)

            !! Fraction of the shortwave flux absorbed by the cool-skin sublayer:
            !zQsw_f = 0.065 + 11.*zdelta - 6.6e-5/zdelta*(1. - EXP(-zdelta/8.e-4)) ! Eq.16 (Fairall al. 1996b)
            zfr = 0.137 + 11.*zdelta - 6.6e-5/zdelta*(1. - EXP(-zdelta/8.e-4)) ! Eq.16 (Fairall al. 1996b)
            !LB: why 0.065 and not 0.137 like in the paper??? Beljaars & Zeng use 0.065
            !LB: maybe comes from Wick et al 2005 ...

            zQt = -(zQlw + zQsen + zQlat + zfr*pQsw(ji,jj))  ! Total cooling at the interface

            ztf    = 0.5 + SIGN(0.5_wp, zQt) ! Qt > 0 => cooling of the layer => ztf = 1
            !                               Qt < 0 => warming of the layer => ztf = 0

            zalpha = 2.1e-5*MAX(pT_s(ji,jj)-rt0 + 3.2_wp, 0._wp)**0.79  ! alpha = thermal expansion of water (~2.5E-4) LB: remove from loop, sst accurate enough!

            !! Term alpha*Qb (Qb is the virtual surface cooling inc. buoyancy effect of salinity due to evap):
            zz1 = zalpha*zQt - 0.026*zQlat*rCp0_w/rLevap  ! alpha*(Eq.8) == alpha*Qb "-" because Qlat < 0
            !! LB: this terms only makes sense if > 0 i.e. in the cooling case
            !! so similar to what's donce in ECMWF:
            zz1 = MAX(0._wp , zz1)    ! 1. instead of 0.1 though ZQ = MAX(1.0,-pQlw(ji,jj) - pQsen(ji,jj) - pQlat(ji,jj))

            !! Laurent: too low wind (u*) might cause problem in stable cases:
            zus = MAX(pus(ji,jj), 1.E-4_wp)

            ! Lambda (=> zz0, empirical coeff.) (Eq.14):
            zz0 = 16. * zz1 * grav * rho0_w * rCp0_w * rnu0_w*rnu0_w*rnu0_w  ! (numerateur) zz1 == alpha*Q
            zz2 = zus*zus * prhoa(ji,jj) / rho0_w * rk0_w
            zz2 =  zz2*zz2                                             ! denominateur
            !LB:  zz0 has the sign of zz1 and therefore of Qb !
            zlamb =  6.*( 1. + (zz0/zz2)**(3./4.) )**(-1./3.) !  Eq.14   (Saunders)

            ! Updating molecular sublayer thickness (delta):
            zz2    = rnu0_w/(SQRT(prhoa(ji,jj)/rho0_w)*zus)
            zdelta =      ztf    *          zlamb*zz2   &  ! Eq.12 (when alpha*Qb>0 / cooling of layer)
               &    + (1. - ztf) * MIN(0.007_wp , 6._wp*zz2 )    ! Eq.12 (when alpha*Qb<0 / warming of layer)
            !LB: changed 0.01 to 0.007
            pdelta(ji,jj) = zdelta

            ! Updating pT_s and q_s ...
            zz2 =  - zQt*zdelta/rk0_w   ! temperature increment !  Eq.13 Cool skin
            !
            pT_s(ji,jj) = psst(ji,jj) + zz2
            !
         END DO
      END DO

      pq_s = rdct_qsat_salt*q_sat(MAX(pT_s, 200._wp), pslp)   !skin !LB: just to avoid problem on masked regions

   END SUBROUTINE CS_COARE3P0

   !!======================================================================
END MODULE mod_cs_coare3p0
