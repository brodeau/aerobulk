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
   USE mod_const       !: physical and other constants
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

   
   SUBROUTINE Cdn10_Lupkes2015( zu, Ts_o, t_zu, q_zu, Uo_zu, Ui_zu, Cdn_o, Chn_o, &
      &                         pic, Ts_i, pslp, pcd, pch )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  Cdn10_Lupkes2015  ***
      !!
      !! ** Purpose :    Alternative turbulent transfer coefficients formulation
      !!                 between sea-ice and atmosphere with distinct momentum
      !!                 and heat coefficients depending on sea-ice concentration
      !!                 and atmospheric stability (no meltponds effect for now).
      !!
      !! ** Method :     The parameterization is adapted from Lupkes et al. (2015)
      !!                 and ECHAM6 atmospheric model. Compared to Lupkes2012 scheme,
      !!                 it considers specific skin and form drags (Andreas et al. 2010)
      !!                 to compute neutral transfer coefficients for both heat and
      !!                 momemtum fluxes. Atmospheric stability effect on transfer
      !!                 coefficient is also taken into account following Louis (1979).
      !!
      !! ** References : Lupkes et al. JGR 2015 (theory)
      !!                 Lupkes et al. ECHAM6 documentation 2015 (implementation)
      !!
      !! ** History :
      !!              - G. Samson (2018,2019) original code
      !!              - L. Brodeau (2020) AeroBulk
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   )                   :: zu     ! reference height (height for Uo_zu)   [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: Ts_o   ! sea surface temperature  [Kelvin]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: t_zu   ! potential air temperature              [Kelvin]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: q_zu   ! specific air humidity at zt             [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: Uo_zu  ! relative wind module at zu over ocean  [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: Ui_zu  ! relative wind module at zu over ice    [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: Cdn_o  ! neutral drag coefficient over ocean
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: Chn_o  ! neutral heat coefficient over ocean
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pic    ! ice concentration [fraction]  => at_i_b
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: Ts_i ! sea-ice surface temperature [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pslp   ! sea-level pressure [Pa]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pcd    ! momentum transfer coefficient
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pch    ! heat transfer coefficient
      !
      REAL(wp) ::   zfi, zfo, zqsat_o, zqsat_i
      REAL(wp) ::   zrib_o, zrib_i, ztmp
      REAL(wp) ::   zCdn_skin_ice, zCdn_form_ice, zCdn_ice
      REAL(wp) ::   zChn_skin_ice, zChn_form_ice
      REAL(wp) ::   z0w, z0i, zfmi, zfmw, zfhi, zfhw
      REAL(wp) ::   zCdn_form_tmp, zwndspd_o, zwndspd_i
      INTEGER  ::   ji, jj         ! dummy loop indices
      !!----------------------------------------------------------------------

      ! Momentum Neutral Transfer Coefficients (should be a constant)
      zCdn_form_tmp = zce10 * ( LOG( 10._wp / z0_form_ice + 1._wp ) / LOG( zu / z0_form_ice + 1._wp ) )**2   ! Eq. 40
      zCdn_skin_ice = ( vkarmn                                      / LOG( zu / z0_skin_ice + 1._wp ) )**2   ! Eq. 7
      zCdn_ice      = zCdn_skin_ice   ! Eq. 7
      !zCdn_ice     = 1.89e-3         ! old ECHAM5 value (cf Eq. 32)

      ! Heat Neutral Transfer Coefficients
      zChn_skin_ice = vkarmn**2 / ( LOG( zu / z0_ice + 1._wp ) * LOG( zu * z1_alpha / z0_skin_ice + 1._wp ) )   ! Eq. 50 + Eq. 52
      
      DO jj = 1, jpj
         DO ji = 1, jpi

            zfi = pic(ji,jj)    ! fraction of sea-ice
            zfo = 1._wp - zfi    ! fraction of open ocean
            
            zwndspd_o = MAX( 0.5, Uo_zu(ji,jj) )
            zwndspd_i = MAX( 0.5, Ui_zu(ji,jj) )
            
            ! Specific humidities at saturation at air-sea and air-ice interface [kg/kg]:
            zqsat_o = rdct_qsat_salt*q_sat( Ts_o(ji,jj), pslp(ji,jj) )
            zqsat_i =                q_sat( Ts_i(ji,jj), pslp(ji,jj) ) !LOLO!!! No!!! Must use a special 1 over ice!
            
            ! Bulk Richardson Number:
            zrib_o = Ri_bulk( zu, Ts_o(ji,jj), t_zu(ji,jj), zqsat_o, q_zu(ji,jj), zwndspd_o ) !LOLO: is t_zu potential T? It should!!!!
            zrib_i = Ri_bulk( zu, Ts_i(ji,jj), t_zu(ji,jj), zqsat_i, q_zu(ji,jj), zwndspd_i )
            
            ! Momentum and Heat Neutral Transfer Coefficients
            zCdn_form_ice = zCdn_form_tmp * zfi * zfo**zbeta                          ! Eq. 40
            zChn_form_ice = zCdn_form_ice / ( 1._wp + ( LOG( z1_alphaf ) / vkarmn ) * SQRT( zCdn_form_ice ) )   ! Eq. 53

            ! Momentum and Heat Stability functions (possibility to use psi_m_ecmwf instead ?)
            z0w = zu * EXP( -1._wp * vkarmn / SQRT( Cdn_o(ji,jj) ) ) ! over water
            z0i = z0_skin_ice                                        ! over ice

            IF( zrib_o <= 0._wp ) THEN
               ztmp = virt_temp( Ts_o(ji,jj), zqsat_o ) - virt_temp( t_zu(ji,jj), q_zu(ji,jj) ) ! difference of potential virtual temperature
               zfmw = 1._wp - zam*zrib_o / ( 1._wp + 3._wp*zc2*Cdn_o(ji,jj)*SQRT( -zrib_o*( zu / z0w + 1._wp ) ) )  ! Eq. 10
               zfhw = ( 1._wp + ( zbetah*ztmp**r1_3 / ( Chn_o(ji,jj) * zwndspd_o ) )**zgamma )**z1_gamma      ! Eq. 26
            ELSE
               ztmp = zrib_o / SQRT( 1._wp + zrib_o )
               zfmw = 1._wp / ( 1._wp + zam * ztmp )   ! Eq. 12
               zfhw = 1._wp / ( 1._wp + zah * ztmp )   ! Eq. 28
            ENDIF

            IF( zrib_i <= 0._wp ) THEN
               ztmp = zrib_i / ( 1._wp + 3._wp * zc2 * zCdn_ice * SQRT( -zrib_i * ( zu / z0i + 1._wp) ) )
               zfmi = 1._wp - zam * ztmp   ! Eq.  9
               zfhi = 1._wp - zah * ztmp   ! Eq. 25
            ELSE
               ztmp = zrib_i / SQRT( 1._wp + zrib_i )
               zfmi = 1._wp / ( 1._wp + zam * ztmp )   ! Eq. 11
               zfhi = 1._wp / ( 1._wp + zah * ztmp )   ! Eq. 27
            ENDIF

            
            
            ! Momentum and Heat transfer coefficients (Eq. 38) and (Eq. 49):
            ztmp       = 1._wp / MAX( 1.e-06, zfi )
            pcd(ji,jj) = zCdn_skin_ice * zfmi + zCdn_form_ice * ( zfmi*zfi + zfmw*zfo ) * ztmp
            pch(ji,jj) = zChn_skin_ice * zfhi + zChn_form_ice * ( zfhi*zfi + zfhw*zfo ) * ztmp
            
         END DO
      END DO
      !
      !CALL lbc_lnk_multi( 'sbcblk', pcd, 'T',  1., pch, 'T', 1. )
      !
   END SUBROUTINE Cdn10_Lupkes2015
   
   !!======================================================================

END MODULE mod_blk_ice_lu15
