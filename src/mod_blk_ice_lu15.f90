! AeroBulk / 2020 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_blk_ice_lu15
   !!====================================================================================
   !!       Computes turbulent components of surface fluxes over sea-ice




   !!       Routine turb_ice_lu15 maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, January 2020
   !!
   !!====================================================================================
   USE mod_const       !: physical and othe constants
   USE mod_phymbl      !: thermodynamics

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: Cdn10_Lupkes2015
   !PUBLIC :: rough_leng_m, rough_leng_tq

   ! ECHAM6 constants
   REAL(wp), PARAMETER ::   z0_skin_ice  = 0.69e-3_wp  ! Eq. 43 [m]
   REAL(wp), PARAMETER ::   z0_form_ice  = 0.57e-3_wp  ! Eq. 42 [m]
   REAL(wp), PARAMETER ::   z0_ice       = 1.00e-3_wp  ! Eq. 15 [m]
   REAL(wp), PARAMETER ::   zce10        = 2.80e-3_wp  ! Eq. 41
   REAL(wp), PARAMETER ::   zbeta        = 1.1_wp      ! Eq. 41
   REAL(wp), PARAMETER ::   zc           = 5._wp       ! Eq. 13
   REAL(wp), PARAMETER ::   zc2          = zc * zc
   REAL(wp), PARAMETER ::   zam          = 2. * zc     ! Eq. 14
   REAL(wp), PARAMETER ::   zah          = 3. * zc     ! Eq. 30
   REAL(wp), PARAMETER ::   z1_alpha     = 1._wp / 0.2_wp  ! Eq. 51
   REAL(wp), PARAMETER ::   z1_alphaf    = z1_alpha    ! Eq. 56
   REAL(wp), PARAMETER ::   zbetah       = 1.e-3_wp    ! Eq. 26
   REAL(wp), PARAMETER ::   zgamma       = 1.25_wp     ! Eq. 26
   REAL(wp), PARAMETER ::   z1_gamma     = 1._wp / zgamma
   REAL(wp), PARAMETER ::   r1_3         = 1._wp / 3._wp

   !!----------------------------------------------------------------------
CONTAINS

   
   SUBROUTINE Cdn10_Lupkes2015( zu, To_s, t_zu, q_zu, Uo_zu, Ui_zu, Cdn_o, Chn_o, &
      &                         pic, ptm_su, pslp, pcd, pch )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  Cdn10_Lupkes2015  ***
      !!
      !! ** pUrpose :    Alternative turbulent transfert coefficients formulation
      !!                 between sea-ice and atmosphere with distinct momentum
      !!                 and heat coefficients depending on sea-ice concentration
      !!                 and atmospheric stability (no meltponds effect for now).
      !!
      !! ** Method :     The parameterization is adapted from Lupkes et al. (2015)
      !!                 and ECHAM6 atmospheric model. Compared to Lupkes2012 scheme,
      !!                 it considers specific skin and form drags (Andreas et al. 2010)
      !!                 to compute neutral transfert coefficients for both heat and
      !!                 momemtum fluxes. Atmospheric stability effect on transfert
      !!                 coefficient is also taken into account following Louis (1979).
      !!
      !! ** References : Lupkes et al. JGR 2015 (theory)
      !!                 Lupkes et al. ECHAM6 documentation 2015 (implementation)
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   )                  ::   zu    ! reference height (height for Uo_zu)   [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   To_s  ! sea surface temperature  [Kelvin]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   t_zu  ! potential air temperature              [Kelvin]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   q_zu  ! specific air humidity at zt             [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   Uo_zu ! relative wind module at zu over ocean  [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   Ui_zu ! relative wind module at zu over ice    [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   Cdn_o ! neutral drag coefficient over ocean
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   Chn_o ! neutral heat coefficient over ocean
      REAL(wp), DIMENSION(:,:), INTENT(in)    :: pic      ! ice concentration [fraction]  => at_i_b
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   ptm_su ! sea-ice surface temperature [K]
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   pslp   ! sea-level pressure [Pa]
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   pcd    ! momentum transfert coefficient
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   pch    ! heat transfert coefficient
      REAL(wp), DIMENSION(jpi,jpj)            ::   zqo_sat, zqi_sat
      !
      !
      INTEGER  ::   ji, jj         ! dummy loop indices
      REAL(wp) ::   zthetav_os, zthetav_is, zthetav_zu
      REAL(wp) ::   zrib_o, zrib_i
      REAL(wp) ::   zCdn_skin_ice, zCdn_form_ice, zCdn_ice
      REAL(wp) ::   zChn_skin_ice, zChn_form_ice
      REAL(wp) ::   z0w, z0i, zfmi, zfmw, zfhi, zfhw
      REAL(wp) ::   zCdn_form_tmp
      !!----------------------------------------------------------------------

      ! Momentum Neutral Transfert Coefficients (should be a constant)
      zCdn_form_tmp = zce10 * ( LOG( 10._wp / z0_form_ice + 1._wp ) / LOG( zu / z0_form_ice + 1._wp ) )**2   ! Eq. 40
      zCdn_skin_ice = ( vkarmn                                      / LOG( zu / z0_skin_ice + 1._wp ) )**2   ! Eq. 7
      zCdn_ice      = zCdn_skin_ice   ! Eq. 7
      !zCdn_ice     = 1.89e-3         ! old ECHAM5 value (cf Eq. 32)

      ! Heat Neutral Transfert Coefficients
      zChn_skin_ice = vkarmn**2 / ( LOG( zu / z0_ice + 1._wp ) * LOG( zu * z1_alpha / z0_skin_ice + 1._wp ) )   ! Eq. 50 + Eq. 52

      ! Atmospheric and Surface Variables
      zqo_sat(:,:) = rdct_qsat_salt * q_sat( To_s(:,:)   , pslp(:,:) )   ! saturation humidity over ocean [kg/kg]
      zqi_sat(:,:) =                  q_sat( ptm_su(:,:), pslp(:,:) )   ! saturation humidity over ice   [kg/kg]

      DO jj = 1, jpj
         DO ji = 1, jpi

            ! Virtual potential temperature [K]
            zthetav_os = To_s(ji,jj)    * ( 1._wp + rctv0 * zqo_sat(ji,jj) )   ! over ocean
            zthetav_is = ptm_su(ji,jj) * ( 1._wp + rctv0 * zqi_sat(ji,jj) )   ! ocean ice
            zthetav_zu = t_zu (ji,jj)  * ( 1._wp + rctv0 * q_zu(ji,jj)    )   ! at zu

            ! Bulk Richardson Number (could use Ri_bulk function from aerobulk instead)
            zrib_o = grav / zthetav_os * ( zthetav_zu - zthetav_os ) * zu / MAX( 0.5, Uo_zu(ji,jj) )**2   ! over ocean
            zrib_i = grav / zthetav_is * ( zthetav_zu - zthetav_is ) * zu / MAX( 0.5, Ui_zu(ji,jj) )**2   ! over ice

            ! Momentum and Heat Neutral Transfert Coefficients
            zCdn_form_ice = zCdn_form_tmp * pic(ji,jj) * ( 1._wp - pic(ji,jj) )**zbeta  ! Eq. 40
            zChn_form_ice = zCdn_form_ice / ( 1._wp + ( LOG( z1_alphaf ) / vkarmn ) * SQRT( zCdn_form_ice ) )               ! Eq. 53

            ! Momentum and Heat Stability functions (possibility to use psi_m_ecmwf instead ?)
            z0w = zu * EXP( -1._wp * vkarmn / SQRT( Cdn_o(ji,jj) ) ) ! over water
            z0i = z0_skin_ice                                             ! over ice
            IF( zrib_o <= 0._wp ) THEN
               zfmw = 1._wp - zam * zrib_o / ( 1._wp + 3._wp * zc2 * Cdn_o(ji,jj) * SQRT( -zrib_o * ( zu / z0w + 1._wp ) ) )  ! Eq. 10
               zfhw = ( 1._wp + ( zbetah * ( zthetav_os - zthetav_zu )**r1_3 / ( Chn_o(ji,jj) * MAX(0.01, Uo_zu(ji,jj)) )   &     ! Eq. 26
                  &             )**zgamma )**z1_gamma
            ELSE
               zfmw = 1._wp / ( 1._wp + zam * zrib_o / SQRT( 1._wp + zrib_o ) )   ! Eq. 12
               zfhw = 1._wp / ( 1._wp + zah * zrib_o / SQRT( 1._wp + zrib_o ) )   ! Eq. 28
            ENDIF

            IF( zrib_i <= 0._wp ) THEN
               zfmi = 1._wp - zam * zrib_i / (1._wp + 3._wp * zc2 * zCdn_ice * SQRT( -zrib_i * ( zu / z0i + 1._wp)))   ! Eq.  9
               zfhi = 1._wp - zah * zrib_i / (1._wp + 3._wp * zc2 * zCdn_ice * SQRT( -zrib_i * ( zu / z0i + 1._wp)))   ! Eq. 25
            ELSE
               zfmi = 1._wp / ( 1._wp + zam * zrib_i / SQRT( 1._wp + zrib_i ) )   ! Eq. 11
               zfhi = 1._wp / ( 1._wp + zah * zrib_i / SQRT( 1._wp + zrib_i ) )   ! Eq. 27
            ENDIF

            ! Momentum Transfert Coefficients (Eq. 38)
            pcd(ji,jj) = zCdn_skin_ice *   zfmi +  &
               &        zCdn_form_ice * ( zfmi * pic(ji,jj) + zfmw * ( 1._wp - pic(ji,jj) ) ) / MAX( 1.e-06, pic(ji,jj) )

            ! Heat Transfert Coefficients (Eq. 49)
            pch(ji,jj) = zChn_skin_ice *   zfhi +  &
               &        zChn_form_ice * ( zfhi * pic(ji,jj) + zfhw * ( 1._wp - pic(ji,jj) ) ) / MAX( 1.e-06, pic(ji,jj) )
            !
         END DO
      END DO
      !
      !CALL lbc_lnk_multi( 'sbcblk', pcd, 'T',  1., pch, 'T', 1. )
      !
   END SUBROUTINE Cdn10_Lupkes2015


   !!======================================================================

END MODULE mod_blk_ice_lu15
