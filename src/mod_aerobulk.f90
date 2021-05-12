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
   USE mod_aerobulk_compute

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: AEROBULK_MODEL
   
CONTAINS



   SUBROUTINE aerobulk_init(psst, pt, pq, pU, pV, pslp)

      !! 1. Check on correct size for array
      !! 2. Set the official 2D shape of the problem jpi,jpj (saved and shared via mod_const)
      !! 3. Check the values in the array make sense
      !! 4. Give some info

      REAL(wp), DIMENSION(:,:), INTENT(in)  :: psst, pt, pq, pU, pV, pslp

      INTEGER :: ni, nj

      PRINT *, ''
      PRINT *, '****************************************************'
      PRINT *, '             ----- AeroBulk_init -----'
      ! 1.
      jpi = size(psst,1) ; jpj = size(psst,2)

      ! 2.
      ni = size(pt,1)  ; nj = size(pt,2)
      IF ( (ni /= jpi).OR.(ni /= jpi) ) THEN
         PRINT *, 'ERROR: aerobulk_init => SST and t_air arrays do not agree in shape!' ; STOP
      END IF

      ! 3.
      ni = SIZE(pq,1)  ; nj = SIZE(pq,2)
      IF ( (ni /= jpi).OR.(ni /= jpi) ) THEN
         PRINT *, 'ERROR: aerobulk_init => SST and q_air arrays do not agree in shape!' ; STOP
      END IF

      ni = SIZE(pU,1)  ; nj = SIZE(pU,2)
      IF ( (ni /= jpi).OR.(ni /= jpi) ) THEN
         PRINT *, 'ERROR: aerobulk_init => SST and U arrays do not agree in shape!' ; STOP
      END IF

      ni = SIZE(pV,1)  ; nj = SIZE(pV,2)
      IF ( (ni /= jpi).OR.(ni /= jpi) ) THEN
         PRINT *, 'ERROR: aerobulk_init => SST and V arrays do not agree in shape!' ; STOP
      END IF

      ni = SIZE(pslp,1)  ; nj = SIZE(pslp,2)
      IF ( (ni /= jpi).OR.(ni /= jpi) ) THEN
         PRINT *, 'ERROR: aerobulk_init => SST and SLP arrays do not agree in shape!' ; STOP
      END IF

      PRINT *, '    *** Computational domain shape: jpi, jpj =', INT(jpi,2), ',', INT(jpj,2)
      PRINT *, '    *** Number of iterations:       nb_iter  =', INT(nb_iter,1)
      PRINT *, '****************************************************'
      PRINT *, ''

      l_first_call = .FALSE.

   END SUBROUTINE aerobulk_init



   SUBROUTINE aerobulk_bye()
      PRINT *, ''
      PRINT *, 'AeroBulk_bye'
      PRINT *, ''
   END SUBROUTINE aerobulk_bye




   
   SUBROUTINE AEROBULK_MODEL( calgo, zt, zu, sst, t_zt,   &
      &                       q_zt, U_zu, V_zu, slp,      &
      &                       QL, QH, Tau_x, Tau_y, Evap,  &
      &                       Niter, rad_sw, rad_lw, T_s  )
      !!======================================================================================
      !!
      !! INPUT :
      !! -------
      !!    *  calgo: what bulk algorithm to use => 'coare'/'ncar'/'ecmwf'
      !!    *  zt   : height for temp. & spec. hum. of air (usually 2 or 10) [m]
      !!    *  zu   : height for wind (usually 10)                           [m]
      !!    *  sst  : SST                                                    [K]
      !!    *  t_zt : potential air temperature at zt                        [K]
      !!    *  q_zt : specific humidity of air at zt                         [kg/kg]
      !!    *  U_zu : zonal wind speed at zu                                 [m/s]
      !!    *  V_zu : meridional wind speed at zu                            [m/s]
      !!    *  slp  : mean sea-level pressure                                [Pa] ~101000 Pa
      !!
      !! OUTPUT :
      !! --------
      !!    *  QL     : Latent heat flux                                     [W/m^2]
      !!    *  QH     : Sensible heat flux                                   [W/m^2]
      !!    *  Tau_x  : zonal wind stress                                    [N/m^2]
      !!    *  Tau_y  : zonal wind stress                                    [N/m^2]
      !!    *  Evap    : Evaporation                                          [mm/s]
      !!
      !! OPTIONAL :
      !! ----------
      !!    *  Niter  : number of itterattions in the bulk algorithm (default is 4)
      !!    *  rad_sw : downwelling shortwave radiation at the surface (>0)   [W/m^2]
      !!    *  rad_lw : downwelling longwave radiation at the surface  (>0)   [W/m^2]
      !!    *  T_s    : skin temperature                                      [K]
      !!
      !!============================================================================
      CHARACTER(len=*),         INTENT(in)  :: calgo
      REAL(wp),                 INTENT(in)  :: zt, zu
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: sst, t_zt, q_zt, U_zu, V_zu, slp
      REAL(wp), DIMENSION(:,:), INTENT(out) :: QL, QH, Tau_x, Tau_y, Evap
      REAL(wp), DIMENSION(:,:), INTENT(in), OPTIONAL :: rad_sw, rad_lw
      REAL(wp), DIMENSION(:,:), INTENT(out),OPTIONAL :: T_s

      INTEGER, INTENT(in), OPTIONAL :: Niter
      
      IF ( PRESENT(Niter) ) nb_iter = Niter  ! Updating number of itterations (define in mod_const)

      IF ( l_first_call ) CALL aerobulk_init(sst, t_zt, q_zt, U_zu, V_zu, slp)


      IF( PRESENT(rad_sw) .AND. PRESENT(rad_lw) ) THEN

         CALL aerobulk_compute(calgo, zt, zu, sst, t_zt, &
            &                  q_zt, U_zu, V_zu, slp,    &
            &                  QL, QH, Tau_x, Tau_y,     &
            &                  rad_sw=rad_sw, rad_lw=rad_lw, T_s=T_s, Evp=Evap )
         
         !PRINT *, 'LOLO DEBUG INTO mod_aerobulk after CALL aerobulk_compute !!! ', TRIM(calgo)
         !PRINT *, 'LOLO: Ts =', T_s
         !PRINT *, 'LOLO: (sst was) =', sst
         !PRINT *, 'LOLO: Evap =', Evap*3600.*24., ' mm/day'
         !PRINT *, ''

      ELSE

         CALL aerobulk_compute(calgo, zt, zu, sst, t_zt, &
            &                  q_zt, U_zu, V_zu, slp,    &
            &                  QL, QH, Tau_x, Tau_y,     &
            &                  Evp=Evap)

      END IF


      IF ( l_last_call ) CALL aerobulk_bye()

   END SUBROUTINE AEROBULK_MODEL

END MODULE mod_aerobulk
