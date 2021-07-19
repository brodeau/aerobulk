! AeroBulk / 2016 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_aerobulk

   USE mod_const
   USE mod_phymbl, ONLY: type_of_humidity, check_unit_consistency
   USE mod_aerobulk_compute

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: AEROBULK_INIT, AEROBULK_MODEL, AEROBULK_BYE

CONTAINS


   SUBROUTINE aerobulk_init( Nt, calgo, psst, pta, pha, pU, pV, pslp,  prsw, prlw )
      !!===================================================================================================================
      !! 1. Set the official 2D shape of the problem based on the `psst` array =>[jpi,jpj] (saved and shared via mod_const)
      !! 2. Check on size agreement between input arrays (must all be [jpi,jpj])
      !! 3. Allocate and fill the `imask` array to disregar "apparent problematic" regions...
      !! 4. Decide the type of humidity in use: specific? relative? dew-point? (saved and shared via mod_const)
      !!
      INTEGER,                  INTENT(in)  :: Nt    !: number of time records to go for...
      CHARACTER(len=*),         INTENT(in)  :: calgo !: what bulk algorithm to use => 'coare'/'ncar'/'ecmwf'/'andreas'
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: psst  !: sea surface temperature             [K]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pta   !: absolute air temperature at z=zt    [K]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pha   !: air humidity at z=zt, can be:
      !!                                                - specific humidity                 [kg/kg]
      !!                                                - dew-point temperature             [K]
      !!                                                - relative humidity                 [%]
      !!                                               => type is normally be recognized based on value range
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pU, pV !: wind vector components at z=zu     [m/s]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pslp   !: sea-level atmospheric pressure     [Pa]
      !!
      !! Output:
      !INTEGER(1), DIMENSION(:,:), INTENT(out) :: imask   !: mask array: masked=>0, elsewhere=>1
      !!
      !! OPTIONAL Output
      REAL(wp), DIMENSION(:,:), INTENT(in), OPTIONAL :: prsw   !: downwelling shortwave radiation  [W/m^2]
      REAL(wp), DIMENSION(:,:), INTENT(in), OPTIONAL :: prlw   !: downwelling  longwave radiation  [W/m^2]
      !LOGICAL,     OPTIONAL   , INTENT(in)  :: l_cswl
      !!===================================================================================================================
      LOGICAL :: lcswl=.FALSE.
      INTEGER :: Ni, Nj, np
      INTEGER(1), DIMENSION(:,:), ALLOCATABLE :: imask   !: mask array: masked=>0, elsewhere=>1
      !REAL(wp), DIMENSION(:,:), ALLOCATABLE  :: ztmp ! kitchen sink array...
      !!===================================================================================================================

      lcswl = ( PRESENT(prsw) .AND. PRESENT(prlw) )
      
      WRITE(6,*)''
      WRITE(6,*)'==================================================================='
      WRITE(6,*)'                   ----- AeroBulk_init -----'
      WRITE(6,*)''
      
      WRITE(6,*)'    *** Bulk parameterization to be used => "', TRIM(calgo), '"'
      ! 1.
      Ni = SIZE(psst,1)
      Nj = SIZE(psst,2)
      
      ! 2.
      IF( ANY(SHAPE(pta) /=(/Ni,Nj/)) ) CALL ctl_stop(' aerobulk_init => SST and t_air arrays do not agree in shape!')
      IF( ANY(SHAPE(pha) /=(/Ni,Nj/)) ) CALL ctl_stop(' aerobulk_init => SST and hum_air arrays do not agree in shape!')
      IF( ANY(SHAPE(pU)  /=(/Ni,Nj/)) ) CALL ctl_stop(' aerobulk_init => SST and U arrays do not agree in shape!')
      IF( ANY(SHAPE(pV)  /=(/Ni,Nj/)) ) CALL ctl_stop(' aerobulk_init => SST and V arrays do not agree in shape!')
      IF( ANY(SHAPE(pslp)/=(/Ni,Nj/)) ) CALL ctl_stop(' aerobulk_init => SST and SLP arrays do not agree in shape!')
      !IF( ANY(SHAPE(imask)/=(/Ni,Nj/)) ) CALL ctl_stop(' aerobulk_init => SST and SLP arrays do not agree in shape!')
      IF( lcswl ) THEN
         IF( ANY(SHAPE(prsw)/=(/Ni,Nj/)) ) CALL ctl_stop(' aerobulk_init => SST and Rad_SW arrays do not agree in shape!')
         IF( ANY(SHAPE(prlw)/=(/Ni,Nj/)) ) CALL ctl_stop(' aerobulk_init => SST and Rad_LW arrays do not agree in shape!')
      END IF
      
      WRITE(6,*)'    *** Computational domain shape: Ni, Nj =', INT(Ni,2), ',', INT(Nj,2)
      
      nitend = Nt
      WRITE(6,*)'    *** Number of time records that will be treated:', nitend      
      
      WRITE(6,*)'    *** Number of iterations in bulk algos: nb_iter  =', INT(nb_iter,1)

      ! 4. Cool-skin/Warm-layer schemes ???
      IF( lcswl ) THEN
         IF( .NOT. ((TRIM(calgo(1:4)) == 'coar').OR.(TRIM(calgo) == 'ecmwf')) ) &
            & CALL ctl_stop('only COARE* and ECMWF algos support Cool-skin/Warm-layer schemes','so do not provide radiative fields')
         WRITE(6,*)'    *** Cool-skin/Warm-layer schemes are going to be used!'
      ELSE
         WRITE(6,*)'    *** Cool-skin/Warm-layer schemes will NOT be used!'
      END IF
      
      ! 3. Allocation and creation of the mask
      !WRITE(6,*)'    *** Allocating the `imask` array!'
      ALLOCATE( imask(Ni,Nj) )
      WRITE(6,*)'    *** Filling the `mask` array...'
      imask(:,:) = 1
      WHERE( (psst < ref_sst_min) .OR. (psst > ref_sst_max) ) imask = 0 ! silly SST
      WHERE( ( pta < ref_taa_min) .OR. (pta  > ref_taa_max) ) imask = 0 ! silly air temperature
      WHERE( (pslp < ref_slp_min) .OR. (pslp > ref_slp_max) ) imask = 0 ! silly atmospheric pressure
      WHERE(     SQRT( pU*pU + pV*pV )       > ref_wnd_max  ) imask = 0 ! silly scalar wind speed
      IF( lcswl ) THEN
         WHERE( ( prsw < ref_rsw_min) .OR. (prsw  > ref_rsw_max) ) imask = 0 ! silly air temperature
         WHERE( ( prlw < ref_rlw_min) .OR. (prlw  > ref_rlw_max) ) imask = 0 ! silly air temperature
      END IF

      
      np = SUM(INT(imask,4))  ! number of valid points (convert to INT(4) before SUM otherwize SUM is INT(1) => overflow!!!)
      IF( np == Ni*Nj ) THEN
         WRITE(6,*)'        ==> no points need to be masked! :)'
      ELSEIF ( np > 0 ) THEN
         WRITE(6,*)'        ==> number of points to mask: ', Ni*Nj-np, ' (out of ',Ni*Nj,')'
      ELSE
         CALL ctl_stop( 'the whole domain is masked!', 'check unit consistency of input fields')
         !STOP
      END IF
      
      !ALLOCATE( ztmp(Ni,Nj) )
      !ztmp(:,:) = SQRT( pU(:,:)*pU(:,:) + pV(:,:)*pV(:,:) )
      
      ! 4. Type of humidity provided?
      ctype_humidity = type_of_humidity( pha, mask=imask )
      WRITE(6,'("     *** Type of air humidity : `",a2,"`")') ctype_humidity

      ! 5. Check unit consistency of input fields:
      CALL check_unit_consistency( 'sst',   psst, mask=imask )
      CALL check_unit_consistency( 't_air', pta,  mask=imask )
      CALL check_unit_consistency( 'slp',   pslp, mask=imask )
      CALL check_unit_consistency( 'u10',   pU,   mask=imask )
      CALL check_unit_consistency( 'v10',   pV,   mask=imask )
      CALL check_unit_consistency( 'wnd',   SQRT(pU*pU + pV*pV), mask=imask )
      CALL check_unit_consistency( ctype_humidity, pha,          mask=imask )
      IF( lcswl ) THEN
         CALL check_unit_consistency( 'rad_sw',   prsw, mask=imask )
         CALL check_unit_consistency( 'rad_lw',   prlw, mask=imask )
      END IF
      
      DEALLOCATE( imask )
      
      WRITE(6,*)'==================================================================='
      !WRITE(6,*)''

      !l_1st_call_ab_init = .FALSE.
      
      !DEALLOCATE( ztmp )

      !STOP'mod_aerobulk.f90'
      
   END SUBROUTINE aerobulk_init



   SUBROUTINE aerobulk_bye()
      WRITE(6,*)'==================================================================='
      WRITE(6,*)'                   ----- AeroBulk_bye -----'
      !DEALLOCATE( imask )
      WRITE(6,*)'==================================================================='
      WRITE(6,*)''
   END SUBROUTINE aerobulk_bye
   




   SUBROUTINE AEROBULK_MODEL( jt, Nt, &
      &                       calgo, zt, zu, sst, t_zt,   &
      &                       hum_zt, U_zu, V_zu, slp,    &
      &                       QL, QH, Tau_x, Tau_y, Evap, &
      &                       Niter, rad_sw, rad_lw, T_s  )
      !!======================================================================================
      !!
      !! INPUT :
      !! -------
      !!    *  jt,Nt  : current time snaphot and total number of time snaphots to go for
      !!    *  calgo  : what bulk algorithm to use => 'coare'/'ncar'/'ecmwf'
      !!    *  zt     : height for temp. & spec. hum. of air (usually 2 or 10) [m]
      !!    *  zu     : height for wind (usually 10)                           [m]
      !!    *  sst    : SST                                                    [K]
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
      !! OUTPUT :
      !! --------
      !!    *  QL     : Latent heat flux                                     [W/m^2]
      !!    *  QH     : Sensible heat flux                                   [W/m^2]
      !!    *  Tau_x  : zonal wind stress                                    [N/m^2]
      !!    *  Tau_y  : zonal wind stress                                    [N/m^2]
      !!    *  Evap    : Evaporation                                          [mm/s] aka [kg/m^2/s] (usually <0, as ocean loses water!)
      !!
      !! OPTIONAL Input:
      !! ---------------
      !!    *  Niter  : number of itterattions in the bulk algorithm (default is 4)
      !!    *  rad_sw : downwelling shortwave radiation at the surface (>0)   [W/m^2]
      !!    *  rad_lw : downwelling longwave radiation at the surface  (>0)   [W/m^2]
      !!
      !! OPTIONAL Output:
      !! ---------------
      !!    *  T_s    : skin temperature                                      [K]
      !!
      !!============================================================================
      INTEGER,                  INTENT(in)  :: jt, Nt
      CHARACTER(len=*),         INTENT(in)  :: calgo
      REAL(wp),                 INTENT(in)  :: zt, zu
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: sst, t_zt, hum_zt, U_zu, V_zu, slp
      REAL(wp), DIMENSION(:,:), INTENT(out) :: QL, QH, Tau_x, Tau_y, Evap
      REAL(wp), DIMENSION(:,:), INTENT(in), OPTIONAL :: rad_sw, rad_lw
      REAL(wp), DIMENSION(:,:), INTENT(out),OPTIONAL :: T_s

      INTEGER, INTENT(in), OPTIONAL :: Niter
      LOGICAL :: l_do_cswl = .FALSE.

      IF( PRESENT(Niter) ) nb_iter = Niter  ! Updating number of itterations (define in mod_const)

      l_do_cswl = ( PRESENT(rad_sw) .AND. PRESENT(rad_lw) ) ! if theser 2 are provided we plan to use the CSWL schemes!

      IF( jt <1 ) CALL ctl_stop('AEROBULK_MODEL => jt < 1 !??', 'we are in a Fortran world here...')

      
      IF( l_do_cswl ) THEN

         IF( jt==1 ) CALL aerobulk_init( Nt, calgo, sst, t_zt, hum_zt, U_zu, V_zu, slp,  prsw=rad_lw, prlw=rad_lw )
         
         CALL aerobulk_compute( jt, calgo, zt, zu, sst, t_zt, &
            &                   hum_zt, U_zu, V_zu, slp,    &
            &                   QL, QH, Tau_x, Tau_y,     &
            &                   rad_sw=rad_sw, rad_lw=rad_lw, T_s=T_s, Evp=Evap )
         
      ELSE
         
         IF( jt==1 ) CALL aerobulk_init( Nt, calgo, sst, t_zt, hum_zt, U_zu, V_zu, slp )
         
         CALL aerobulk_compute( jt, calgo, zt, zu, sst, t_zt, &
            &                   hum_zt, U_zu, V_zu, slp,  &
            &                   QL, QH, Tau_x, Tau_y,     &
            &                   Evp=Evap)
         
      END IF
      

      IF( jt==Nt ) CALL aerobulk_bye()
      
   END SUBROUTINE AEROBULK_MODEL


   !FUNCTION sum_mask( imsk )
   !   INTEGER(1), DIMENSION(:,:), INTENT(in) :: imsk
   !   INTEGER(1), DIMENSION(:,:), INTENT(in :: imsk
   !   
   !   INTEGER
   !
   !END FUNCTION sum_mask
   
END MODULE mod_aerobulk
