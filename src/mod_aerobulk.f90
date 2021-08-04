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

   SUBROUTINE AEROBULK_INIT( Nt, calgo, psst, pta, pha, pU, pV, pslp,  l_use_skin, prsw, prlw )
      !!===========================================================================================
      !! 1. Set the official 2D shape of the problem based on the `psst` array =>[Ni,Nj]
      !! 2. Check on size agreement between input arrays (must all be [Ni,Nj])
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
      !! OPTIONAL Input
      LOGICAL,                  INTENT(in), OPTIONAL :: l_use_skin
      REAL(wp), DIMENSION(:,:), INTENT(in), OPTIONAL :: prsw   !: downwelling shortwave radiation  [W/m^2]
      REAL(wp), DIMENSION(:,:), INTENT(in), OPTIONAL :: prlw   !: downwelling  longwave radiation  [W/m^2]
      !!==================================================================================================
      LOGICAL :: lsrad, lskin
      INTEGER :: Ni, Nj, np
      CHARACTER(len=64) :: chum_ln
      INTEGER(1), DIMENSION(:,:), ALLOCATABLE :: imask   !: mask array: masked=>0, elsewhere=>1
      !!==================================================================================================
      lsrad = .FALSE.
      lsrad = ( PRESENT(prsw) .AND. PRESENT(prlw) )

      lskin = .FALSE.
      IF( PRESENT(l_use_skin) ) lskin = l_use_skin

      WRITE(6,*)''
      WRITE(6,*)'==================================================================='
      WRITE(6,*)'                   ----- AeroBulk_init -----'
      WRITE(6,*)''

      WRITE(6,*)'    *** Bulk parameterization to be used => "', TRIM(calgo), '"'

      ! 1. Check if scheme use compatible with chosen algorithm and presence of radiative fields
      IF( lskin ) THEN

         IF( .NOT.((TRIM(calgo(1:4)) == 'coar').OR.(TRIM(calgo) == 'ecmwf')) ) &
            & CALL ctl_stop(' AEROBULK_INIT => Only `COARE*` and `ECMWF` algorithms support cool-skin & warm/layer schemes')

         IF( .NOT.(lsrad) ) CALL ctl_stop(' AEROBULK_INIT => provide SW and LW rad. input if you want to use skin schemes')

         l_use_skin_schemes = .TRUE. ! if we go here it's fine!
         WRITE(6,*)'       ==> will use the Cool-skin & Warm-ayer scheme of `'//TRIM(calgo)//'` !'

      ELSE
         WRITE(6,*)'    *** Cool-skin & Warm-layer schemes will NOT be used!'
      END IF


      ! 2.
      Ni = SIZE(psst,1)
      Nj = SIZE(psst,2)

      ! 3.
      IF( ANY(SHAPE(pta) /=(/Ni,Nj/)) ) CALL ctl_stop(' AEROBULK_INIT => SST and t_air arrays do not agree in shape!')
      IF( ANY(SHAPE(pha) /=(/Ni,Nj/)) ) CALL ctl_stop(' AEROBULK_INIT => SST and hum_air arrays do not agree in shape!')
      IF( ANY(SHAPE(pU)  /=(/Ni,Nj/)) ) CALL ctl_stop(' AEROBULK_INIT => SST and U arrays do not agree in shape!')
      IF( ANY(SHAPE(pV)  /=(/Ni,Nj/)) ) CALL ctl_stop(' AEROBULK_INIT => SST and V arrays do not agree in shape!')
      IF( ANY(SHAPE(pslp)/=(/Ni,Nj/)) ) CALL ctl_stop(' AEROBULK_INIT => SST and SLP arrays do not agree in shape!')
      IF( lsrad ) THEN
         IF( ANY(SHAPE(prsw)/=(/Ni,Nj/)) ) CALL ctl_stop(' AEROBULK_INIT => SST and Rad_SW arrays do not agree in shape!')
         IF( ANY(SHAPE(prlw)/=(/Ni,Nj/)) ) CALL ctl_stop(' AEROBULK_INIT => SST and Rad_LW arrays do not agree in shape!')
      END IF

      WRITE(6,'("     *** Computational domain shape: Ni x Nj = ",i5.5," x ",i5.5)') Ni, Nj

      nitend = Nt   ! important: nitend is a global variable shared by `mod_const.f90`
      WRITE(6,*)'    *** Number of time records that will be treated:', nitend

      WRITE(6,*)'    *** Number of iterations in bulk algos: nb_iter  =', INT(nb_iter,1)

      ! 4. Allocation and creation of the mask
      ALLOCATE( imask(Ni,Nj) )
      WRITE(6,*)'    *** Filling the `mask` array...'
      imask(:,:) = 1
      WHERE( (psst < ref_sst_min) .OR. (psst > ref_sst_max) ) imask = 0 ! silly SST
      WHERE( ( pta < ref_taa_min) .OR. (pta  > ref_taa_max) ) imask = 0 ! silly air temperature
      WHERE( (pslp < ref_slp_min) .OR. (pslp > ref_slp_max) ) imask = 0 ! silly atmospheric pressure
      WHERE(     SQRT( pU*pU + pV*pV )       > ref_wnd_max  ) imask = 0 ! silly scalar wind speed
      IF( lsrad ) THEN
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

      ! 5. Type of humidity provided?
      ctype_humidity = type_of_humidity( pha, mask=imask )
      SELECT CASE(ctype_humidity)
      CASE('sh')
         chum_ln = 'specific humidity [kg/kg]'
      CASE('rh')
         chum_ln = 'relative humidity [%]'
      CASE('dp')
         chum_ln = 'dew-point temperature [K]'
      CASE DEFAULT
         CALL ctl_stop( ' AEROBULK_INIT => humidty type "',ctype_humidity,'" is unknown!!!' )
      END SELECT


      WRITE(6,'("     *** Type of prescribed air humidity  `",a,"`")') TRIM(chum_ln)

      ! 6. Check unit consistency of input fields:
      CALL check_unit_consistency( 'sst',   psst, mask=imask )
      CALL check_unit_consistency( 't_air', pta,  mask=imask )
      CALL check_unit_consistency( 'slp',   pslp, mask=imask )
      CALL check_unit_consistency( 'u10',   pU,   mask=imask )
      CALL check_unit_consistency( 'v10',   pV,   mask=imask )
      CALL check_unit_consistency( 'wnd',   SQRT(pU*pU + pV*pV), mask=imask )
      CALL check_unit_consistency( ctype_humidity, pha,          mask=imask )
      IF( lsrad ) THEN
         CALL check_unit_consistency( 'rad_sw',   prsw, mask=imask )
         CALL check_unit_consistency( 'rad_lw',   prlw, mask=imask )
      END IF

      DEALLOCATE( imask )

      WRITE(6,*)'==================================================================='
      !WRITE(6,*)''

   END SUBROUTINE AEROBULK_INIT



   SUBROUTINE AEROBULK_BYE()
      WRITE(6,*)'==================================================================='
      WRITE(6,*)'                   ----- AeroBulk_bye -----'
      !DEALLOCATE( imask )
      WRITE(6,*)'==================================================================='
      WRITE(6,*)''
   END SUBROUTINE AEROBULK_BYE





   SUBROUTINE AEROBULK_MODEL( jt, Nt, &
      &                       calgo, zt, zu, sst, t_zt,   &
      &                       hum_zt, U_zu, V_zu, slp,    &
      &                       QL, QH, Tau_x, Tau_y, Evap, &
      &                       Niter, l_use_skin, rad_sw, rad_lw, T_s  )
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
      !!    *  l_use_skin: should we use the cool-skin & warm-layer schemes fot skin temperature ?
      !!    *  rad_sw : downwelling shortwave radiation at the surface (>0)   [W/m^2]
      !!    *  rad_lw : downwelling longwave radiation at the surface  (>0)   [W/m^2]
      !!
      !! OPTIONAL Output:
      !! ---------------
      !!    *  T_s    : skin temperature                                      [K]
      !!
      !!====================================================================================================
      INTEGER,                  INTENT(in)  :: jt, Nt
      CHARACTER(len=*),         INTENT(in)  :: calgo
      REAL(wp),                 INTENT(in)  :: zt, zu
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: sst, t_zt, hum_zt, U_zu, V_zu, slp
      REAL(wp), DIMENSION(:,:), INTENT(out) :: QL, QH, Tau_x, Tau_y, Evap
      !! Optional input:
      INTEGER,                  INTENT(in), OPTIONAL :: Niter
      LOGICAL,                  INTENT(in), OPTIONAL :: l_use_skin
      REAL(wp), DIMENSION(:,:), INTENT(in), OPTIONAL :: rad_sw, rad_lw
      !! Optional output:
      REAL(wp), DIMENSION(:,:), INTENT(out),OPTIONAL :: T_s
      !!
      !! Local:
      LOGICAL :: lskin, lsrad
      !!====================================================================================================

      IF( PRESENT(Niter) )      nb_iter = Niter  ! Updating number of itterations (define in mod_const)

      lskin = .FALSE.
      IF( PRESENT(l_use_skin) ) lskin = l_use_skin

      lsrad = .FALSE.
      lsrad = ( PRESENT(rad_sw) .AND. PRESENT(rad_lw) ) ! if theser 2 are provided we NORMALLY plan to use the CSWL schemes!

      IF( jt <1 ) CALL ctl_stop('AEROBULK_MODEL => jt < 1 !??', 'we are in a Fortran world here...')

      IF( lsrad ) THEN

         IF( jt==1 ) CALL AEROBULK_INIT( Nt, calgo, sst, t_zt, hum_zt, U_zu, V_zu, slp,  l_use_skin=lskin, prsw=rad_lw, prlw=rad_lw )

         CALL AEROBULK_COMPUTE( jt, calgo, zt, zu, sst, t_zt, &
            &                   hum_zt, U_zu, V_zu, slp,      &
            &                   QL, QH, Tau_x, Tau_y,         &
            &                   rad_sw=rad_sw, rad_lw=rad_lw, T_s=T_s, Evp=Evap )

      ELSE

         IF( jt==1 ) CALL AEROBULK_INIT( Nt, calgo, sst, t_zt, hum_zt, U_zu, V_zu, slp,  l_use_skin=lskin )

         CALL AEROBULK_COMPUTE( jt, calgo, zt, zu, sst, t_zt, &
            &                   hum_zt, U_zu, V_zu, slp,  &
            &                   QL, QH, Tau_x, Tau_y,     &
            &                   Evp=Evap)

      END IF


      IF( jt==Nt ) CALL AEROBULK_BYE()

   END SUBROUTINE AEROBULK_MODEL

END MODULE mod_aerobulk
