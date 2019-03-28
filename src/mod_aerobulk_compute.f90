! AeroBulk / 2015 / L. Brodeau

MODULE mod_aerobulk_compute

   USE mod_const       !: physical constants
   USE mod_thermo      !: thermodynamics functions

   USE mod_blk_coare   !: COAREv3   algorithm
   USE mod_blk_ncar    !: Large & Yeager algorithm
   USE mod_blk_ecmwf   !: following ECMWF doc...

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: aerobulk_compute

CONTAINS

   SUBROUTINE aerobulk_compute( calgo, zt, zu, sst, t_zt, &
      &                         q_zt, U_zu, V_zu, slp,    &
      &                         QL, QH, Tau_x, Tau_y,     &
      &                         mask, rad_sw, rad_lw, T_s )
      !!
      !!******************************
      !! 2015: L. Brodeau
      !!  => all constants taken from mod_thermo and mod_const must
      !!     be done or used from NEMO constant bank...    ... vkarmn ... grav ...
      !!******************************
      !!
      !!======================================================================================
      !!
      !! INPUT :
      !! -------
      !!    *  calgo: what bulk algorithm to use => 'coare'/'coare35'/'ncar'/'ecmwf'
      !!    *  zt   : height for temperature and spec. hum. of air           [m]
      !!    *  zu   : height for wind (10m = traditional anemometric height  [m]
      !!    *  sst  : bulk SST                                               [K]
      !!    *  t_zt : air temperature at zt                                  [K]
      !!    *  q_zt : specific humidity of air at zt                         [kg/kg]
      !!    *  U_zu : zonal wind speed at zu                                 [m/s]
      !!    *  V_zu : meridional wind speed at zu                            [m/s]
      !!    *  slp  : mean sea-level pressure                                [Pa] ~101000 Pa
      !!
      !! OPTIONAL INPUT (will trigger l_use_skin=TRUE if present!):
      !! ---------------
      !!    *  mask:    regions where mask==0 are ignored, mask is 1 elsewhere    !!!
      !!    *  rad_sw : downwelling shortwave radiation at the surface (>0)   [W/m^2]
      !!    *  rad_lw : downwelling longwave radiation at the surface  (>0)   [W/m^2]
      !!
      !! OUTPUT :
      !! --------
      !!    *  QL     : Latent heat flux                                     [W/m^2]
      !!    *  QH     : Sensible heat flux                                   [W/m^2]
      !!    *  Tau_x  : zonal wind stress                                    [N/m^2]
      !!    *  Tau_y  : zonal wind stress                                    [N/m^2]
      !!
      !! OPTIONAL OUTPUT
      !! ---------------
      !!    *  T_s : skin temperature    [K]
      !!             (only when l_use_skin=TRUE)
      !!
      !!============================================================================
      !!
      !! I/O ARGUMENTS:
      CHARACTER(len=*),             INTENT(in)  :: calgo
      REAL(wp),                     INTENT(in)  :: zt, zu
      REAL(wp),   DIMENSION(jpi,jpj), INTENT(in)  :: sst, t_zt, q_zt, U_zu, V_zu, slp
      REAL(wp),   DIMENSION(jpi,jpj), INTENT(out) :: QL, QH, Tau_x, Tau_y
      INTEGER(1), DIMENSION(jpi,jpj), INTENT(in), OPTIONAL :: mask
      REAL(wp),   DIMENSION(jpi,jpj), INTENT(in), OPTIONAL :: rad_sw, rad_lw
      REAL(wp),   DIMENSION(jpi,jpj), INTENT(out),OPTIONAL :: T_s
      !!
      INTEGER(1), DIMENSION(:,:), ALLOCATABLE  :: pmask
      REAL(wp),   DIMENSION(:,:), ALLOCATABLE  ::  &
         &     pWzu,            & !: Scalar wind speed at zu m
         &   pSSQ,              & !: Specific humidiyt at the air-sea interface
         &   pCd, pCh, pCe,     & !: bulk transfer coefficients
         &  pTzt,               & !: potential temperature at zt meters
         &  pTzu, pQzu,         & !: potential temperature and specific humidity at zu meters
         &  pTs, pqs,           & !:
         &   pUblk                !: Bulk scalar wind speed (pWzu corrected for low wind and unstable conditions)

      LOGICAL :: l_use_skin
      !!------------------------------------------------------------------------------


      l_use_skin = .FALSE.

      ALLOCATE ( pmask(jpi,jpj), pWzu(jpi,jpj), pSSQ(jpi,jpj), &
         &     pCd(jpi,jpj), pCh(jpi,jpj), pCe(jpi,jpj),       &
         &     pTzt(jpi,jpj), pTzu(jpi,jpj), pQzu(jpi,jpj), &
         &     pUblk(jpi,jpj), pTs(jpi,jpj), pqs(jpi,jpj)  )

      ! Masked region ?
      IF( PRESENT(mask) ) THEN
         pmask(:,:) = mask(:,:)
      ELSE
         pmask(:,:) = 1
      END IF


      ! Cool skin ?
      IF( PRESENT(rad_sw) .AND. PRESENT(rad_lw) ) THEN
         IF((TRIM(calgo) == 'coare').OR.(TRIM(calgo) == 'coare35').OR.(TRIM(calgo) == 'ecmwf')) THEN
            l_use_skin = .TRUE.
            PRINT *, ''; PRINT *, ' *** Will use the cool-skin warm-layer scheme of ', TRIM(calgo(1:5)), '!'
         END IF
         CALL check_unit_consitency( 'rad_sw', rad_sw, pmask )
         CALL check_unit_consitency( 'rad_lw', rad_lw, pmask )
      END IF

      CALL check_unit_consitency( 'sst',   sst,  pmask )
      CALL check_unit_consitency( 't_air', t_zt, pmask )
      CALL check_unit_consitency( 'q_air', q_zt, pmask )
      CALL check_unit_consitency( 'slp',   slp,  pmask )
      CALL check_unit_consitency( 'u10', ABS(U_zu), pmask )
      CALL check_unit_consitency( 'v10', ABS(V_zu), pmask )


      !! Scalar wind:
      pWzu = sqrt( U_zu*U_zu + V_zu*V_zu )

      !! Computing specific humidity at saturation at sea surface temperature :
      pSSQ (:,:) = 0.98*q_sat(sst, slp) !! lolo/crude / NEMO 3.6  (slp not used!)

      !! Approximate potential temperarure at zt meters above sea surface:
      !pTzt = t_zt + gamma_moist(t_zt, q_zt)*zt
      pTzt = t_zt !! lolo/crude / NEMO 3.6

      !! Mind that TURB_COARE and TURB_ECMWF will modify SST and SSQ if their
      !! respective Cool Skin Warm Layer parameterization is used
      pTs = sst
      pqs = pSSQ


      SELECT CASE(TRIM(calgo))
         !!
      CASE('coare')
         IF( l_use_skin ) THEN
            CALL TURB_COARE ( '3.0', zt, zu, pTs, pTzt, pqs, q_zt, pWzu,  &
               &              pCd, pCh, pCe, pTzu, pQzu, pUblk,            &
               &              rad_sw=rad_sw, rad_lw=rad_lw, slp=slp )
         ELSE
            CALL TURB_COARE ( '3.0', zt, zu, pTs, pTzt, pqs, q_zt, pWzu,  &
               &              pCd, pCh, pCe, pTzu, pQzu, pUblk )
         END IF
         !!
      CASE('coare35')
         IF( l_use_skin ) THEN
            CALL TURB_COARE ( '3.5', zt, zu, pTs, pTzt, pqs, q_zt, pWzu, &
               &              pCd, pCh, pCe, pTzu, pQzu, pUblk,           &
               &              rad_sw=rad_sw, rad_lw=rad_lw, slp=slp )
         ELSE
            CALL TURB_COARE ('3.5', zt, zu, pTs, pTzt, pqs, q_zt, pWzu,  &
               &              pCd, pCh, pCe, pTzu, pQzu, pUblk )
         END IF
         !!
         !!
      CASE('ncar')
         CALL TURB_NCAR( zt, zu, pTs, pTzt, pqs, q_zt, pWzu, &
            &            pCd, pCh, pCe, pTzu, pQzu, pUblk)
         !!
         !!
      CASE('ecmwf')
         IF( l_use_skin ) THEN
            CALL TURB_ECMWF ( zt, zu, pTs, pTzt, pqs, q_zt, pWzu,   &
               &              pCd, pCh, pCe, pTzu, pQzu, pUblk,      &
               &              rad_sw=rad_sw, rad_lw=rad_lw, slp=slp )
         ELSE
            CALL TURB_ECMWF ( zt, zu, pTs, pTzt, pqs, q_zt, pWzu,   &
               &              pCd, pCh, pCe, pTzu, pQzu, pUblk)
         END IF
         !!
         !!
      CASE DEFAULT
         write(6,*) 'ERROR: mod_aerobulk_compute.f90 => bulk algorithm ', trim(calgo), ' is unknown!!!'
         STOP
      END SELECT

      !! Skin temperature:
      !! IF( l_use_skin ), pTs has been updated from SST to skin temperature !
      IF( l_use_skin ) T_s = pTs

      !! Need the air density at zu m, so using t and q corrected at zu m:
      !pRHO = rho_air(pTzu, pQzu, slp)
      !QH   = slp - pRHO*grav*zu      ! QH used as temporary array!
      !pRHO = rho_air(pTzu, pQzu, QH) !! lolo/crude / NEMO 3.6

      !! *** Wind stress ***
      Tau_x = pCd*1.22_wp * U_zu * pUblk !! lolo/crude / NEMO 3.6
      Tau_y = pCd*1.22_wp * V_zu * pUblk !! lolo/crude / NEMO 3.6

      !! *** Latent and Sensible heat fluxes ***
      QL = pCe*1.22_wp * 2.5e6     * (pQzu - pqs) * pUblk  !! lolo/crude / NEMO 3.6
      QH = pCh*1.22_wp * 1000.5_wp * (pTzu - pTs) * pUblk  !! lolo/crude / NEMO 3.6

      DEALLOCATE ( pmask, pWzu, pSSQ, pCd, pCh, pCe, pTzt, pTzu, pQzu, pUblk, pTs, pqs )

   END SUBROUTINE aerobulk_compute




   SUBROUTINE check_unit_consitency( cfield, Xval, mask )

      !! Ignore values where mask==0

      CHARACTER(len=*),         INTENT(in) :: cfield
      REAL(wp),   DIMENSION(:,:), INTENT(in) :: Xval
      INTEGER(1), DIMENSION(:,:), INTENT(in) :: mask

      CHARACTER(len=64) :: cunit
      REAL(wp) :: zmean, vmin, vmax
      LOGICAL  :: l_too_large=.FALSE., l_too_small=.FALSE., l_mean_outside=.FALSE.

      zmean = SUM( Xval * REAL(mask,wp) ) / SUM( REAL(mask,wp) )

      !PRINT *, 'LOLO, zmean of '//TRIM(cfield)//' =>', zmean

      SELECT CASE (TRIM(cfield))

      CASE('sst')
         vmax = 313._wp
         vmin = 270._wp
         cunit = 'K'

      CASE('t_air')
         vmax = 323._wp
         vmin = 200._wp
         cunit = 'K'

      CASE('q_air')
         vmax = 0.08_wp
         vmin = 0._wp
         cunit = 'kg/kg'

      CASE('slp')
         vmax = 108000.
         vmin =  87000.
         cunit = 'Pa'

      CASE('u10')
         vmax = 50._wp
         vmin =  0._wp
         cunit = 'm/s'

      CASE('v10')
         vmax = 50._wp
         vmin =  0._wp
         cunit = 'm/s'

      CASE('rad_sw')
         vmax = 1500.0_wp
         vmin =  0._wp
         cunit = 'W/m^2'

      CASE('rad_lw')
         vmax = 700.0_wp
         vmin =  0._wp
         cunit = 'W/m^2'

      CASE DEFAULT
         WRITE(*,'(" *** ERROR (mod_aerobulk_compute.f90): we do not know field ",a," !")') TRIM(cfield)
         STOP
      END SELECT

      IF ( MAXVAL(Xval) > vmax )                   l_too_large    = .TRUE.
      IF ( MINVAL(Xval) < vmin )                   l_too_small    = .TRUE.
      IF ( (zmean < vmin) .OR. (zmean > vmax) ) l_mean_outside = .TRUE.

      IF ( l_too_large .OR. l_too_small .OR. l_mean_outside ) THEN
         WRITE(*,'(" *** ERROR (mod_aerobulk_compute.f90): field ",a," does not seem to be in ",a," !")') TRIM(cfield), TRIM(cunit)
         WRITE(*,'(" min value = ", es9.3," max value = ", es9.3," mean value = ", es9.3)') MINVAL(Xval), MAXVAL(Xval), zmean
         STOP
      END IF

   END SUBROUTINE check_unit_consitency


END MODULE mod_aerobulk_compute
