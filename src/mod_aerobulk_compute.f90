! AeroBulk / 2015 / L. Brodeau

MODULE mod_aerobulk_compute

   USE mod_const        !: physical constants
   USE mod_phymbl       !: thermodynamics functions

   USE mod_blk_coare3p0 !: COARE v3.0 algorithm
   USE mod_blk_coare3p6 !: COARE v3.6 algorithm
   USE mod_blk_ncar     !: Large & Yeager algorithm
   USE mod_blk_ecmwf    !: following ECMWF doc...
   USE mod_blk_andreas  !: Andreas et al.

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: aerobulk_compute

CONTAINS

   SUBROUTINE aerobulk_compute( calgo, zt, zu, sst, t_zt, &
      &                         hum_zt, U_zu, V_zu, slp,    &
      &                         QL, QH, Tau_x, Tau_y,     &
      &                         mask, rad_sw, rad_lw,     &
      &                         T_s, Evp )
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
      !!    *  calgo  : what bulk algorithm to use => 'coare3p6'/'ncar'/'ecmwf'
      !!    *  zt     : height for temperature and spec. hum. of air           [m]
      !!    *  zu     : height for wind (10m = traditional anemometric height  [m]
      !!    *  sst    : bulk SST                                               [K]
      !!    *  t_zt   : ABSOLUTE air temperature at zt                         [K]
      !!    *  hum_zt : air humidity at zt, can be given as:
      !!                - specific humidity                                [kg/kg]
      !!                - dew-point temperature                                [K]
      !!                - relative humidity                                    [%]
      !!               => type should normally be recognized based on value range
      !!    *  U_zu   : zonal wind speed at zu                                 [m/s]
      !!    *  V_zu   : meridional wind speed at zu                            [m/s]
      !!    *  slp    : mean sea-level pressure                                [Pa] ~101000 Pa
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
      !!    *  T_s  : skin temperature    [K]    (only when l_use_skin=TRUE)
      !!    *  Evp : evaporation         [mm/s] aka [kg/m^2/s] (usually <0, as ocean loses water!)
      !!
      !!
      !!============================================================================
      !!
      !! I/O ARGUMENTS:
      CHARACTER(len=*),             INTENT(in)  :: calgo
      REAL(wp),                     INTENT(in)  :: zt, zu
      REAL(wp),   DIMENSION(jpi,jpj), INTENT(in)  :: sst, t_zt, hum_zt, U_zu, V_zu, slp
      REAL(wp),   DIMENSION(jpi,jpj), INTENT(out) :: QL, QH, Tau_x, Tau_y
      INTEGER(1), DIMENSION(jpi,jpj), INTENT(in), OPTIONAL :: mask
      REAL(wp),   DIMENSION(jpi,jpj), INTENT(in), OPTIONAL :: rad_sw, rad_lw
      REAL(wp),   DIMENSION(jpi,jpj), INTENT(out),OPTIONAL :: T_s, Evp
      !!
      INTEGER(1), DIMENSION(:,:), ALLOCATABLE  :: zmask
      REAL(wp),   DIMENSION(:,:), ALLOCATABLE  ::  &
         &     zWzu,            & !: Scalar wind speed at zu m
         &   zSSQ,              & !: Specific humidiyt at the air-sea interface
         &   zCd, zCh, zCe,     & !: bulk transfer coefficients
         &  zThtzt, zQzt,       & !: POTENTIAL air temperature and specific humidity at zt meters
         &  zThtzu, zQzu,       & !: POTENTIAL air temperature and specific humidity at zu meters
         &  zTs, zqs, zEvap,    & !:
         &  zTaum,              & !: wind stress module
         &  zUblk,              & !: Bulk scalar wind speed (zWzu corrected for low wind and unstable conditions)
         &  ztmp                  !: temporary array
      !
      CHARACTER(len=2) :: chum
      LOGICAL          :: l_use_skin
      !!------------------------------------------------------------------------------


      l_use_skin = .FALSE.

      ALLOCATE ( zmask(jpi,jpj),  zWzu(jpi,jpj),   zSSQ(jpi,jpj), &
         &       zCd(jpi,jpj),    zCh(jpi,jpj),    zCe(jpi,jpj),  &
         &       zThtzt(jpi,jpj), zQzt(jpi,jpj),                  &
         &       zThtzu(jpi,jpj), zQzu(jpi,jpj),                  &
         &       zUblk(jpi,jpj),  zTs(jpi,jpj),    zqs(jpi,jpj),  &
         &       zEvap(jpi,jpj),  zTaum(jpi,jpj),  ztmp(jpi,jpj)  )

      ! Masked region ?
      IF( PRESENT(mask) ) THEN
         zmask(:,:) = mask(:,:)
      ELSE
         zmask(:,:) = 1
      END IF

      ! Type of humidity provided?
      chum = type_of_humidity(hum_zt, zmask)
      PRINT *, 'LOLO: humidity type is: "'//chum//'" !'

      ! Conversion to specific humidity when needed:
      SELECT CASE(chum)
      CASE('sh')
         zQzt(:,:) =           hum_zt(:,:)                        ! already specific humidity!
      CASE('dp')
         zQzt(:,:) = q_air_dp( hum_zt(:,:),            slp(:,:) ) ! dew-point to specific humidity
      CASE('rh')
         zQzt(:,:) = q_air_rh( hum_zt(:,:), t_zt(:,:), slp(:,:) ) ! relative to specific humidity
      CASE DEFAULT
         WRITE(6,*) 'ERROR: mod_aerobulk_compute.f90 => humidty type "',chum,'" is unknown!!!' ; STOP
      END SELECT


      PRINT *, 'LOLO: as "'//chum//'", "sh" =>', hum_zt(10,10), zQzt(10,10)


      ! Cool skin ?
      IF((TRIM(calgo(1:4)) == 'coar').OR.(TRIM(calgo) == 'ecmwf')) THEN
         IF( PRESENT(rad_sw) .AND. PRESENT(rad_lw) ) THEN
            l_use_skin = .TRUE.
            PRINT *, ' *** we use the cool-skin/warm-layer scheme of ', TRIM(calgo(1:5)), '!'; !lilo
         END IF
         CALL check_unit_consitency( 'rad_sw', rad_sw, zmask )
         CALL check_unit_consitency( 'rad_lw', rad_lw, zmask )
      END IF

      CALL check_unit_consitency( 'sst',   sst,  zmask )
      CALL check_unit_consitency( 't_air', t_zt, zmask )
      CALL check_unit_consitency( 'q_air', zQzt, zmask )
      CALL check_unit_consitency( 'slp',   slp,  zmask )
      CALL check_unit_consitency( 'u10', ABS(U_zu), zmask )
      CALL check_unit_consitency( 'v10', ABS(V_zu), zmask )

      !! Scalar wind:
      zWzu = sqrt( U_zu*U_zu + V_zu*V_zu )

      !! Computing specific humidity at saturation at sea surface temperature :
      zSSQ (:,:) = rdct_qsat_salt*q_sat(sst, slp)

      !! Approximate potential temperarure at zt meters above sea surface:
      !zThtzt = t_zt + gamma_moist(t_zt, zQzt)*zt  !BAD! and should have used `gamma_dry` anyway...
      zThtzt = Theta_from_z_P0_T_q( zt, slp, t_zt, zQzt )

      !! Mind that TURB_COARE* and TURB_ECMWF will modify SST and SSQ if their
      !! respective Cool Skin Warm Layer parameterization is used
      zTs = sst
      zqs = zSSQ


      ztmp(:,:) = 0._wp   ! longitude fixed to 0, (for COAREx when using cool-skin/warm-layer)


      SELECT CASE(TRIM(calgo))
         !!
      CASE('coare3p0')
         IF( l_use_skin ) THEN
            CALL TURB_COARE3P0 ( 1, zt, zu, zTs, zThtzt, zqs, zQzt, zWzu, .TRUE., .TRUE., &
               &                zCd, zCh, zCe, zThtzu, zQzu, zUblk,            &
               &                Qsw=(1._wp - roce_alb0)*rad_sw, rad_lw=rad_lw, slp=slp, isecday_utc=12, plong=ztmp  )
         ELSE
            CALL TURB_COARE3P0 ( 1, zt, zu, zTs, zThtzt, zqs, zQzt, zWzu, .FALSE., .FALSE.,  &
               &                zCd, zCh, zCe, zThtzu, zQzu, zUblk )
         END IF
         !!
      CASE('coare3p6')
         IF( l_use_skin ) THEN
            CALL TURB_COARE3P6 ( 1, zt, zu, zTs, zThtzt, zqs, zQzt, zWzu, .TRUE., .TRUE., &
               &                zCd, zCh, zCe, zThtzu, zQzu, zUblk,            &
               &                Qsw=(1._wp - roce_alb0)*rad_sw, rad_lw=rad_lw, slp=slp, isecday_utc=12, plong=ztmp  )
         ELSE
            CALL TURB_COARE3P6 ( 1, zt, zu, zTs, zThtzt, zqs, zQzt, zWzu, .FALSE., .FALSE.,  &
               &                 zCd, zCh, zCe, zThtzu, zQzu, zUblk )
         END IF
         !!
      CASE('ncar')
         CALL TURB_NCAR( zt, zu, zTs, zThtzt, zqs, zQzt, zWzu, &
            &            zCd, zCh, zCe, zThtzu, zQzu, zUblk)
         !!
         !!
      CASE('ecmwf')
         IF( l_use_skin ) THEN
            CALL TURB_ECMWF ( 1, zt, zu, zTs, zThtzt, zqs, zQzt, zWzu, .TRUE., .TRUE.,  &
               &              zCd, zCh, zCe, zThtzu, zQzu, zUblk,      &
               &              Qsw=(1._wp - roce_alb0)*rad_sw, rad_lw=rad_lw, slp=slp  )
         ELSE
            CALL TURB_ECMWF ( 1, zt, zu, zTs, zThtzt, zqs, zQzt, zWzu, .FALSE., .FALSE.,  &
               &              zCd, zCh, zCe, zThtzu, zQzu, zUblk)
         END IF
         !!
         !!
      CASE('andreas')
         CALL TURB_ANDREAS( zt, zu, zTs, zThtzt, zqs, zQzt, zWzu, &
            &               zCd, zCh, zCe, zThtzu, zQzu, zUblk)
         !!
         !!
      CASE DEFAULT
         WRITE(6,*) 'ERROR: mod_aerobulk_compute.f90 => bulk algorithm ', trim(calgo), ' is unknown!!!'
         STOP
      END SELECT




      !! Skin temperature:
      !! IF( l_use_skin ) => zTs and zqs have been updated from SST to skin temperature !


      !PRINT *, 'LOLO DEBUG INTO mod_aerobulk_compute !!! ', TRIM(calgo)
      !PRINT *, 'LOLO: Ts =', zTs
      !PRINT *, 'LOLO: (sst was) =', sst
      !PRINT *, ''


      CALL BULK_FORMULA( zu, zTs, zqs, zThtzu, zQzu, zCd, zCh, zCe, zWzu, zUblk, slp, &
         &                                 zTaum, QH, QL, pEvap=zEvap )

      Tau_x = zTaum / zWzu * U_zu
      Tau_y = zTaum / zWzu * V_zu

      !PRINT *, 'LOLO DEBUG INTO mod_aerobulk_compute !!! ', TRIM(calgo)
      !PRINT *, 'zCe =', zCe
      !PRINT *, 'Qlat =', QL
      !PRINT *, 'Ublk =', zUblk
      !PRINT *, 'zCe/Ublk =', zCe/zUblk
      !PRINT *, 't_zu =', zThtzu
      !PRINT *, 'q_zu =', zQzu
      !PRINT *, 'ssq =', zSSQ
      !PRINT *, ''

      IF( PRESENT(T_s) ) T_s(:,:)   = zTs(:,:)

      IF( PRESENT(Evp) ) Evp(:,:) = zEvap(:,:)

      DEALLOCATE ( zmask, zWzu, zSSQ, zCd, zCh, zCe, zThtzt, zQzt, zThtzu, zQzu, &
         &         zUblk, zTs, zqs, zEvap, zTaum, ztmp  )


   END SUBROUTINE aerobulk_compute




   SUBROUTINE check_unit_consitency( cfield, Xval, mask )

      !! Ignore values where mask==0

      CHARACTER(len=*),         INTENT(in) :: cfield
      REAL(wp),   DIMENSION(:,:), INTENT(in) :: Xval
      INTEGER(1), DIMENSION(:,:), INTENT(in) :: mask

      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: lmask
      INTEGER :: nx, ny
      CHARACTER(len=64) :: cunit
      REAL(wp) :: zmean, zmin, zmax
      LOGICAL  :: l_too_large=.FALSE., l_too_small=.FALSE., l_mean_outside=.FALSE.

      nx = SIZE(mask,1)
      ny = SIZE(mask,2)
      ALLOCATE( lmask(nx,ny) )
      lmask(:,:) = .FALSE.
      WHERE( mask==1 ) lmask = .TRUE.


      zmean = SUM( Xval * REAL(mask,wp) ) / SUM( REAL(mask,wp) )

      !PRINT *, 'LOLO, zmean of '//TRIM(cfield)//' =>', zmean

      SELECT CASE (TRIM(cfield))

      CASE('sst')
         zmax = 313._wp
         zmin = 270._wp
         cunit = 'K'

      CASE('t_air')
         zmax = 323._wp
         zmin = 200._wp
         cunit = 'K'

      CASE('q_air')
         zmax = 0.08_wp
         zmin = 0._wp
         cunit = 'kg/kg'

      CASE('slp')
         zmax = 108000.
         zmin =  87000.
         cunit = 'Pa'

      CASE('u10')
         zmax = 50._wp
         zmin =  0._wp
         cunit = 'm/s'

      CASE('v10')
         zmax = 50._wp
         zmin =  0._wp
         cunit = 'm/s'

      CASE('rad_sw')
         zmax = 1500.0_wp
         zmin =  0._wp
         cunit = 'W/m^2'

      CASE('rad_lw')
         zmax = 700.0_wp
         zmin =  0._wp
         cunit = 'W/m^2'

      CASE DEFAULT
         WRITE(*,'(" *** ERROR (mod_aerobulk_compute.f90): we do not know field ",a," !")') TRIM(cfield)
         STOP
      END SELECT

      IF ( MAXVAL(Xval, MASK=lmask) > zmax )                   l_too_large    = .TRUE.
      IF ( MINVAL(Xval, MASK=lmask) < zmin )                   l_too_small    = .TRUE.
      IF ( (zmean < zmin) .OR. (zmean > zmax) ) l_mean_outside = .TRUE.

      IF ( l_too_large .OR. l_too_small .OR. l_mean_outside ) THEN
         WRITE(*,'(" *** ERROR (mod_aerobulk_compute.f90): field ",a," does not seem to be in ",a," !")') TRIM(cfield), TRIM(cunit)
         WRITE(*,'(" min value = ", es10.3," max value = ", es10.3," mean value = ", es10.3)') MINVAL(Xval), MAXVAL(Xval), zmean
         STOP
      END IF
      
      DEALLOCATE( lmask )
      
   END SUBROUTINE check_unit_consitency




   FUNCTION type_of_humidity( Xval, mask )

      !! Ignore values where mask==0

      REAL(wp),   DIMENSION(:,:), INTENT(in) :: Xval
      INTEGER(1), DIMENSION(:,:), INTENT(in) :: mask
      CHARACTER(len=2)                       :: type_of_humidity

      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: lmask
      INTEGER :: nx, ny
      CHARACTER(len=64) :: cunit
      REAL(wp) :: zmean, zmin, zmax
      LOGICAL  :: l_too_large=.FALSE., l_too_small=.FALSE., l_mean_outside=.FALSE.
      
      nx = SIZE(mask,1)
      ny = SIZE(mask,2)
      ALLOCATE( lmask(nx,ny) )
      lmask(:,:) = .FALSE.
      WHERE( mask==1 ) lmask = .TRUE.

      zmean =    SUM( Xval * REAL(mask,wp) ) / SUM( REAL(mask,wp) )
      zmin  = MINVAL( Xval,  MASK=lmask )
      zmax  = MAXVAL( Xval,  MASK=lmask )
      
      !PRINT *, 'LOLO, zmean of '//TRIM(cfield)//' =>', zmean
      type_of_humidity = '00'

      IF( (zmean >= 0._wp).AND.(zmean < 0.08_wp).AND.(zmin >= 0._wp).AND.(zmax < 0.08_wp) ) THEN
         !! Specific humidity! [kg/kg]
         type_of_humidity = 'sh'

      ELSEIF( (zmean >= 150._wp).AND.(zmean < 325._wp).AND.(zmin >= 150._wp).AND.(zmax < 325._wp) ) THEN
         !! Dew-point temperature [K]
         type_of_humidity = 'dp'

      ELSEIF( (zmean >= 0._wp).AND.(zmean <= 100._wp).AND.(zmin >= 0._wp).AND.(zmax <= 100._wp) ) THEN
         !! Relative humidity [%]
         type_of_humidity = 'rh'

      ELSE
         WRITE(6,*) 'ERROR: type_of_humidity()@mod_aerobulk_compute => un-identified humidity type!'
         WRITE(6,*) '   ==> we could not identify the humidity type based on the mean, min & max of the field:'
         WRITE(6,*) '     * mean =', REAL(zmean,4)
         WRITE(6,*) '     * min  =', REAL(zmin, 4)
         WRITE(6,*) '     * max  =', REAL(zmax, 4)
         STOP
      END IF
      
      DEALLOCATE( lmask )

   END FUNCTION type_of_humidity


END MODULE mod_aerobulk_compute
