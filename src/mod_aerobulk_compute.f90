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
      &                         rad_sw, rad_lw, T_s )
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
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: sst, t_zt, q_zt, U_zu, V_zu, slp
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: QL, QH, Tau_x, Tau_y
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in), OPTIONAL :: rad_sw, rad_lw
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out),OPTIONAL :: T_s

      REAL(wp), DIMENSION(:,:), ALLOCATABLE  ::  &
         &     XWzu,            & !: Scalar wind speed at zu m
         &   XSSQ,              & !: Specific humidiyt at the air-sea interface
         &   Cd, Ch, Ce,        & !: bulk transfer coefficients
         &  XTzt,               & !: potential temperature at zt meters
         &  XTzu, XQzu,         & !: potential temperature and specific humidity at zu meters
         &  Ts, qs,             & !:
         &   XUblk,             & !: Bulk scalar wind speed (XWzu corrected for low wind and unstable conditions)
         &   XRHO                 !: density of air

      LOGICAL :: l_use_skin = .FALSE.
      !!------------------------------------------------------------------------------

      ALLOCATE ( XWzu(jpi,jpj), XSSQ(jpi,jpj), &
         &     Cd(jpi,jpj), Ch(jpi,jpj), Ce(jpi,jpj),       &
         &     XTzt(jpi,jpj), XTzu(jpi,jpj), XQzu(jpi,jpj), &
         &     XUblk(jpi,jpj), XRHO(jpi,jpj), Ts(jpi,jpj), qs(jpi,jpj)  )


      ! Cool skin ?
      IF( PRESENT(rad_sw) .AND. PRESENT(rad_lw) ) THEN
         IF((TRIM(calgo) == 'coare').OR.(TRIM(calgo) == 'coare35').OR.(TRIM(calgo) == 'ecmwf')) THEN
            l_use_skin = .TRUE.
            PRINT *, ' *** Will use the cool-skin warm-layer scheme of ', TRIM(calgo(1:5)), '!'
         END IF
      END IF

      !! Scalar wind:
      XWzu = sqrt( U_zu*U_zu + V_zu*V_zu )

      !! Computing specific humidity at saturation at sea surface temperature :
      XSSQ (:,:) = 0.98*q_sat(sst, slp)

      !! Approximate potential temperarure at zt meters high:
      XTzt = t_zt + gamma_moist(t_zt, q_zt)*zt

      !! Mind that TURB_COARE and TURB_ECMWF will modify SST and SSQ if their
      !! respective Cool Skin Warm Layer parameterization is used
      Ts = sst ; qs = XSSQ


      SELECT CASE(TRIM(calgo))
         !!
      CASE('coare')
         IF( l_use_skin ) THEN
            CALL TURB_COARE ( '3.0', zt, zu, Ts, XTzt, qs, q_zt, XWzu,  &
               &              Cd, Ch, Ce, XTzu, XQzu, XUblk,            &
               &              rad_sw=rad_sw, rad_lw=rad_lw, slp=slp )            
         ELSE
            CALL TURB_COARE ( '3.0', zt, zu, Ts, XTzt, qs, q_zt, XWzu,  &
               &              Cd, Ch, Ce, XTzu, XQzu, XUblk )
         END IF
         !!
      CASE('coare35')
         IF( l_use_skin ) THEN
            CALL TURB_COARE ( '3.5', zt, zu, Ts, XTzt, qs, q_zt, XWzu, &
               &              Cd, Ch, Ce, XTzu, XQzu, XUblk,           &
               &              rad_sw=rad_sw, rad_lw=rad_lw, slp=slp )
         ELSE
            CALL TURB_COARE ('3.5', zt, zu, Ts, XTzt, qs, q_zt, XWzu,  &
               &              Cd, Ch, Ce, XTzu, XQzu, XUblk )
         END IF
         !!
         !!
      CASE('ncar')
         CALL TURB_NCAR( zt, zu, Ts, XTzt, qs, q_zt, XWzu, &
            &            Cd, Ch, Ce, XTzu, XQzu, XUblk)
         !!
         !!
      CASE('ecmwf')
         IF( l_use_skin ) THEN
            CALL TURB_ECMWF ( zt, zu, Ts, XTzt, qs, q_zt, XWzu,   &
               &              Cd, Ch, Ce, XTzu, XQzu, XUblk,      &
               &              rad_sw=rad_sw, rad_lw=rad_lw, slp=slp )
         ELSE
            CALL TURB_ECMWF ( zt, zu, Ts, XTzt, qs, q_zt, XWzu,   &
               &              Cd, Ch, Ce, XTzu, XQzu, XUblk)
         END IF
         !!
         !!
      CASE DEFAULT
         write(6,*) 'ERROR: mod_aerobulk_compute.f90 => bulk algorithm ', trim(calgo), ' is unknown!!!'
         STOP
      END SELECT

      !! Skin temperature:
      !! IF( l_use_skin ), Ts has been updated from SST to skin temperature !
      T_s = Ts
      
      !! Need the air density at zu m, so using t and q corrected at zu m:
      XRHO = rho_air(XTzu, XQzu, slp)
      QH   = slp - XRHO*grav*zu      ! QH used as temporary array!
      XRHO = rho_air(XTzu, XQzu, QH)


      !! *** Wind stress ***
      Tau_x = Cd*XRHO * U_zu * XUblk
      Tau_y = Cd*XRHO * V_zu * XUblk


      !! *** Latent and Sensible heat fluxes ***
      QL = Ce*XRHO*Lvap(Ts)     * (XQzu - qs) * XUblk
      QH = Ch*XRHO*cp_air(XQzu) * (XTzu - Ts) * XUblk


      !PRINT *, 'LOLO DEBUG INTO mod_aerobulk_compute !!! ', TRIM(calgo)
      !PRINT *, 'Ce =', Ce
      !PRINT *, 'Qlat =', QL
      !PRINT *, 'Ublk =', XUblk
      !PRINT *, 'Ce/Ublk =', Ce/XUblk
      !PRINT *, 't_zu =', XTzu
      !PRINT *, 'q_zu =', XQzu
      !PRINT *, 'Rho =', XRHO
      !PRINT *, 'ssq =', XSSQ
      !PRINT *, 'Lvap =', Lvap(Ts)
      !PRINT *, ''


      DEALLOCATE ( XWzu, XSSQ, Cd, Ch, Ce, XTzt, XTzu, XQzu, XUblk, XRHO, Ts, qs )

   END SUBROUTINE aerobulk_compute

END MODULE mod_aerobulk_compute
