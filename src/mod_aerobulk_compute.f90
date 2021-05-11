! AeroBulk / 2015 / L. Brodeau

MODULE mod_aerobulk_compute

   USE mod_const        !: physical constants
   USE mod_phymbl       !: thermodynamics functions
   
   USE mod_blk_coare3p6 !: COARE v3.5 algorithm
   USE mod_blk_ncar     !: Large & Yeager algorithm
   USE mod_blk_ecmwf    !: following ECMWF doc...

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
      !!  => all constants taken from mod_phymbl and mod_const must
      !!     be done or used from NEMO constant bank...    ... vkarmn ... grav ...
      !!******************************
      !!
      !!======================================================================================
      !!
      !! INPUT :
      !! -------
      !!    *  calgo: what bulk algorithm to use => 'coare3p6'/'ncar'/'ecmwf'
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
      INTEGER(1), DIMENSION(:,:), ALLOCATABLE  :: zmask
      REAL(wp),   DIMENSION(:,:), ALLOCATABLE  ::  &
         &     zWzu,            & !: Scalar wind speed at zu m
         &   zSSQ,              & !: Specific humidiyt at the air-sea interface
         &   zCd, zCh, zCe,     & !: bulk transfer coefficients
         &  zTzt,               & !: potential temperature at zt meters
         &  pTzu, zQzu,         & !: potential temperature and specific humidity at zu meters
         &  zTs, zqs,           & !:
         &  zTaum,              & !: wind stress module
         &  zUblk,              & !: Bulk scalar wind speed (zWzu corrected for low wind and unstable conditions)
         &  ztmp                  !: temporary array

      LOGICAL :: l_use_skin
      !!------------------------------------------------------------------------------


      l_use_skin = .FALSE.

      ALLOCATE ( zmask(jpi,jpj), zWzu(jpi,jpj), zSSQ(jpi,jpj), &
         &     zCd(jpi,jpj), zCh(jpi,jpj), zCe(jpi,jpj),       &
         &     zTzt(jpi,jpj), pTzu(jpi,jpj), zQzu(jpi,jpj),    &
         &     zUblk(jpi,jpj), zTs(jpi,jpj), zqs(jpi,jpj),     &
         &     zTaum(jpi,jpj), ztmp(jpi,jpj)  )

      ! Masked region ?
      IF( PRESENT(mask) ) THEN
         zmask(:,:) = mask(:,:)
      ELSE
         zmask(:,:) = 1
      END IF


      ! Cool skin ?
      IF((TRIM(calgo(1:5)) == 'coare').OR.(TRIM(calgo) == 'ecmwf')) THEN
         IF( PRESENT(rad_sw) .AND. PRESENT(rad_lw) ) THEN
            l_use_skin = .TRUE.
            PRINT *, ''; PRINT *, ' *** Will use the cool-skin warm-layer scheme of ', TRIM(calgo(1:5)), '!'; !lilo
         END IF
         CALL check_unit_consitency( 'rad_sw', rad_sw, zmask )
         CALL check_unit_consitency( 'rad_lw', rad_lw, zmask )
      END IF

      CALL check_unit_consitency( 'sst',   sst,  zmask )
      CALL check_unit_consitency( 't_air', t_zt, zmask )
      CALL check_unit_consitency( 'q_air', q_zt, zmask )
      CALL check_unit_consitency( 'slp',   slp,  zmask )
      CALL check_unit_consitency( 'u10', ABS(U_zu), zmask )
      CALL check_unit_consitency( 'v10', ABS(V_zu), zmask )


      !! Scalar wind:
      zWzu = sqrt( U_zu*U_zu + V_zu*V_zu )

      !! Computing specific humidity at saturation at sea surface temperature :
      zSSQ (:,:) = rdct_qsat_salt*q_sat(sst, slp)

      !! Approximate potential temperarure at zt meters above sea surface:
      zTzt = t_zt + gamma_moist(t_zt, q_zt)*zt

      !! Mind that TURB_COARE* and TURB_ECMWF will modify SST and SSQ if their
      !! respective Cool Skin Warm Layer parameterization is used
      zTs = sst
      zqs = zSSQ


      ztmp(:,:) = 0._wp   ! longitude fixed to 0, (for COAREx when using cool-skin/warm-layer)
      

      SELECT CASE(TRIM(calgo))
         !!
      CASE('coare3p6')
         IF( l_use_skin ) THEN
            CALL TURB_COARE3P6 ( 1, zt, zu, zTs, zTzt, zqs, q_zt, zWzu, .TRUE., .TRUE., &
               &              zCd, zCh, zCe, pTzu, zQzu, zUblk,            &
               &              Qsw=(1._wp - roce_alb0)*rad_sw, rad_lw=rad_lw, slp=slp, isecday_utc=12, plong=ztmp  )
         ELSE
            CALL TURB_COARE3P6 ( 1, zt, zu, zTs, zTzt, zqs, q_zt, zWzu, .FALSE., .FALSE.,  &
               &              zCd, zCh, zCe, pTzu, zQzu, zUblk )
         END IF
         !!
      CASE('ncar')
         CALL TURB_NCAR( zt, zu, zTs, zTzt, zqs, q_zt, zWzu, &
            &            zCd, zCh, zCe, pTzu, zQzu, zUblk)
         !!
         !!
      CASE('ecmwf')
         IF( l_use_skin ) THEN
            CALL TURB_ECMWF ( 1, zt, zu, zTs, zTzt, zqs, q_zt, zWzu, .TRUE., .TRUE.,  &
               &              zCd, zCh, zCe, pTzu, zQzu, zUblk,      &
               &              Qsw=(1._wp - roce_alb0)*rad_sw, rad_lw=rad_lw, slp=slp  )
         ELSE
            CALL TURB_ECMWF ( 1, zt, zu, zTs, zTzt, zqs, q_zt, zWzu, .FALSE., .FALSE.,  &
               &              zCd, zCh, zCe, pTzu, zQzu, zUblk)
         END IF
         !!
         !!
      CASE DEFAULT
         write(6,*) 'ERROR: mod_aerobulk_compute.f90 => bulk algorithm ', trim(calgo), ' is unknown!!!'
         STOP
      END SELECT

      !! Skin temperature:
      !! IF( l_use_skin ) => zTs and zqs have been updated from SST to skin temperature !

      CALL BULK_FORMULA( zu, zTs, zqs, pTzu, zQzu, zCd, zCh, zCe, zWzu, zUblk, slp, &
         &                                 zTaum, QH, QL )

      Tau_x = zTaum / zWzu * U_zu
      Tau_y = zTaum / zWzu * V_zu

      !PRINT *, 'LOLO DEBUG INTO mod_aerobulk_compute !!! ', TRIM(calgo)
      !PRINT *, 'zCe =', zCe
      !PRINT *, 'Qlat =', QL
      !PRINT *, 'Ublk =', zUblk
      !PRINT *, 'zCe/Ublk =', zCe/zUblk
      !PRINT *, 't_zu =', pTzu
      !PRINT *, 'q_zu =', zQzu
      !PRINT *, 'ssq =', zSSQ
      !PRINT *, ''

      DEALLOCATE ( zmask, zWzu, zSSQ, zCd, zCh, zCe, zTzt, pTzu, zQzu, zUblk, zTs, zqs, zTaum )

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
         WRITE(*,'(" min value = ", es10.3," max value = ", es10.3," mean value = ", es10.3)') MINVAL(Xval), MAXVAL(Xval), zmean
         STOP
      END IF

   END SUBROUTINE check_unit_consitency


END MODULE mod_aerobulk_compute
