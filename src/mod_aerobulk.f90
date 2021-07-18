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
   USE mod_phymbl, ONLY: type_of_humidity
   USE mod_aerobulk_compute

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: AEROBULK_INIT, AEROBULK_MODEL, AEROBULK_BYE

CONTAINS



   SUBROUTINE aerobulk_init( Nt, calgo, psst, pta, pha, pU, pV, pslp,  l_cswl )
      !!
      !! 1. Set the official 2D shape of the problem based on the `psst` array =>[jpi,jpj] (saved and shared via mod_const)
      !! 2. Check on size agreement between input arrays (must all be [jpi,jpj])
      !! 3. Allocate and fill the `mask` array to disregar "apparent problematic" regions...
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
      LOGICAL,     OPTIONAL   , INTENT(in)  :: l_cswl

      LOGICAL :: lcswl=.FALSE.
      INTEGER :: ni, nj

      !REAL(wp), DIMENSION(:,:), ALLOCATABLE  :: ztmp ! kitchen sink array...
      IF( PRESENT(l_cswl) ) lcswl=l_cswl
      
      WRITE(6,*)''
      WRITE(6,*)'==================================================================='
      WRITE(6,*)'                   ----- AeroBulk_init -----'
      WRITE(6,*)''

      WRITE(6,*)'    *** Bulk parameterization to be used => "', TRIM(calgo), '"'
      ! 1.
      jpi = SIZE(psst,1) ; jpj = SIZE(psst,2)

      ! 2.
      ni = SIZE(pta,1)  ; nj = SIZE(pta,2)      
      IF( (ni /= jpi).OR.(ni /= jpi) ) CALL ctl_stop(' aerobulk_init => SST and t_air arrays do not agree in shape!')
      ni = SIZE(pha,1)  ; nj = SIZE(pha,2)
      IF( (ni /= jpi).OR.(ni /= jpi) ) CALL ctl_stop(' aerobulk_init => SST and hum_air arrays do not agree in shape!')
      ni = SIZE(pU,1)  ; nj = SIZE(pU,2)
      IF( (ni /= jpi).OR.(ni /= jpi) ) CALL ctl_stop(' aerobulk_init => SST and U arrays do not agree in shape!')
      ni = SIZE(pV,1)  ; nj = SIZE(pV,2)
      IF( (ni /= jpi).OR.(ni /= jpi) ) CALL ctl_stop(' aerobulk_init => SST and V arrays do not agree in shape!')
      ni = SIZE(pslp,1)  ; nj = SIZE(pslp,2)
      IF( (ni /= jpi).OR.(ni /= jpi) ) CALL ctl_stop(' aerobulk_init => SST and SLP arrays do not agree in shape!')

      WRITE(6,*)'    *** Computational domain shape: jpi, jpj =', INT(jpi,2), ',', INT(jpj,2)
      
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
      WRITE(6,*)'    *** Allocating the `mask` array!'
      ALLOCATE( mask(jpi,jpj) )
      WRITE(6,*)'    *** Filling the `mask` array...'
      mask(:,:) = 1
      WHERE( (psst < 270._wp)   .OR. (psst > 320._wp)    ) mask = 0 ! silly SST
      WHERE( ( pta < 180._wp)   .OR. (pta  > 330._wp)    ) mask = 0 ! silly air temperature
      WHERE( (pslp < 80000._wp) .OR. (pslp > 110000._wp) ) mask = 0 ! silly atmospheric pressure
      WHERE(     SQRT( pU*pU + pV*pV ) > 50._wp          ) mask = 0 ! silly scalar wind speed
      
      ni = SUM(mask)  ! number of valid points
      IF( ni == jpi*jpj ) THEN
         WRITE(6,*)'        ==> no points need to be masked! :)'
      ELSEIF( ni > 0 )    THEN
         WRITE(6,*)'        ==> number of points we need to mask: ', jpi*jpj-ni, ' (out of ',jpi*jpj,')'
      ELSE
         WRITE(6,*)'        ==> the whole domain would be masked! :('
         CALL ctl_stop(' one of your input fields must have the wrong unit!', 'check them and come back')
      END IF
      
      !ALLOCATE( ztmp(jpi,jpj) )
      !ztmp(:,:) = SQRT( pU(:,:)*pU(:,:) + pV(:,:)*pV(:,:) )

      
      ! 4. Type of humidity provided?
      ctype_humidity = type_of_humidity( pha, mask )
      WRITE(6,*)'    *** Air humidity type :   ctype_humidity  = ', ctype_humidity
      
      WRITE(6,*)'==================================================================='
      !WRITE(6,*)''

      !l_1st_call_ab_init = .FALSE.
      
      !DEALLOCATE( ztmp )
      
   END SUBROUTINE aerobulk_init



   SUBROUTINE aerobulk_bye()
      WRITE(6,*)'==================================================================='
      WRITE(6,*)'                   ----- AeroBulk_bye -----'
      DEALLOCATE( mask )
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
      !! OPTIONAL :
      !! ----------
      !!    *  Niter  : number of itterattions in the bulk algorithm (default is 4)
      !!    *  rad_sw : downwelling shortwave radiation at the surface (>0)   [W/m^2]
      !!    *  rad_lw : downwelling longwave radiation at the surface  (>0)   [W/m^2]
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
      
      IF( jt==1 ) CALL aerobulk_init( Nt, calgo, sst, t_zt, hum_zt, U_zu, V_zu, slp,  l_cswl=l_do_cswl )

      IF( l_do_cswl ) THEN
         
         CALL aerobulk_compute( jt, calgo, zt, zu, sst, t_zt, &
            &                   hum_zt, U_zu, V_zu, slp,    &
            &                   QL, QH, Tau_x, Tau_y,     &
            &                   rad_sw=rad_sw, rad_lw=rad_lw, T_s=T_s, Evp=Evap )

         !WRITE(6,*)'LOLO DEBUG INTO mod_aerobulk after CALL aerobulk_compute !!! ', TRIM(calgo)
         !WRITE(6,*)'LOLO: Ts =', T_s
         !WRITE(6,*)'LOLO: (sst was) =', sst
         !WRITE(6,*)'LOLO: Evap =', Evap*3600.*24., ' mm/day'
         !WRITE(6,*)''

      ELSE

         CALL aerobulk_compute( jt, calgo, zt, zu, sst, t_zt, &
            &                   hum_zt, U_zu, V_zu, slp,  &
            &                   QL, QH, Tau_x, Tau_y,     &
            &                   Evp=Evap)
         
      END IF


      IF( jt==Nt ) CALL aerobulk_bye()
      
   END SUBROUTINE AEROBULK_MODEL

END MODULE mod_aerobulk
