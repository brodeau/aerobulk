MODULE sbcblk
   !!======================================================================
   !!                       ***  MODULE  sbcblk  ***
   !! Ocean forcing:  momentum, heat and freshwater flux formulation
   !!                         Aerodynamic Bulk Formulas
   !!                        SUCCESSOR OF "sbcblk_core"
   !!=====================================================================
   !! History :  1.0  !  2004-08  (U. Schweckendiek)  Original CORE code
   !!            2.0  !  2005-04  (L. Brodeau, A.M. Treguier)  improved CORE bulk and its user interface
   !!            3.0  !  2006-06  (G. Madec)  sbc rewritting
   !!             -   !  2006-12  (L. Brodeau)  Original code for turb_core
   !!            3.2  !  2009-04  (B. Lemaire)  Introduce iom_put
   !!            3.3  !  2010-10  (S. Masson)  add diurnal cycle
   !!            3.4  !  2011-11  (C. Harris)  Fill arrays required by CICE
   !!            3.7  !  2014-06  (L. Brodeau)  simplification and optimization of CORE bulk
   !!            4.0  !  2016-06  (L. Brodeau)  sbcblk_core becomes sbcblk and is not restricted to the CORE algorithm anymore
   !!                 !                        ==> based on AeroBulk (http://aerobulk.sourceforge.net/)
   !!            4.0  !  2016-10  (G. Madec)  introduce a sbc_blk_init routine
   !!            4.0  !  2016-10  (M. Vancoppenolle)  Introduce Jules emulator (M. Vancoppenolle) 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_blk_init  : initialisation of the chosen bulk formulation as ocean surface boundary condition
   !!   sbc_blk       : bulk formulation as ocean surface boundary condition
   !!   blk_oce       : computes momentum, heat and freshwater fluxes over ocean
   !!   rho_air       : density of (moist) air (depends on T_air, q_air and SLP
   !!   cp_air        : specific heat of (moist) air (depends spec. hum. q_air)
   !!   q_sat         : saturation humidity as a function of SLP and temperature
   !!   L_vap         : latent heat of vaporization of water as a function of temperature
   !!             sea-ice case only : 
   !!   blk_ice_tau   : provide the air-ice stress
   !!   blk_ice_flx   : provide the heat and mass fluxes at air-ice interface
   !!   blk_ice_qcn   : provide ice surface temperature and snow/ice conduction flux (emulating JULES coupler)
   !!   Cdn10_Lupkes2012 : Lupkes et al. (2012) air-ice drag
   !!   Cdn10_Lupkes2015 : Lupkes et al. (2015) air-ice drag 
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE fldread        ! read input fields
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE cyclone        ! Cyclone 10m wind form trac of cyclone centres
   USE sbcdcy         ! surface boundary condition: diurnal cycle
   USE sbcwave , ONLY :   cdn_wave ! wave module
   USE sbc_ice        ! Surface boundary condition: ice fields
   USE lib_fortran    ! to use key_nosignedzero
#if defined key_si3
   USE ice     , ONLY :   u_ice, v_ice, jpl, a_i_b, at_i_b, tm_su, rn_cnd_s, hfx_err_dif
   USE icethd_dh      ! for CALL ice_thd_snwblow
#endif
   USE sbcblk_algo_ncar     ! => turb_ncar     : NCAR - CORE (Large & Yeager, 2009) 
   USE sbcblk_algo_coare    ! => turb_coare    : COAREv3.0 (Fairall et al. 2003) 
   USE sbcblk_algo_coare3p5 ! => turb_coare3p5 : COAREv3.5 (Edson et al. 2013)
   USE sbcblk_algo_ecmwf    ! => turb_ecmwf    : ECMWF (IFS cycle 31) 
   !
   USE iom            ! I/O manager library
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE prtctl         ! Print control

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_blk_init  ! called in sbcmod
   PUBLIC   sbc_blk       ! called in sbcmod
#if defined key_si3
   PUBLIC   blk_ice_tau   ! routine called in iceforcing
   PUBLIC   blk_ice_flx   ! routine called in iceforcing
   PUBLIC   blk_ice_qcn   ! routine called in iceforcing
#endif 

!!Lolo: should ultimately be moved in the module with all physical constants ?
!!gm  : In principle, yes.
   REAL(wp), PARAMETER ::   Cp_dry = 1005.0       !: Specic heat of dry air, constant pressure      [J/K/kg]
   REAL(wp), PARAMETER ::   Cp_vap = 1860.0       !: Specic heat of water vapor, constant pressure  [J/K/kg]
   REAL(wp), PARAMETER ::   R_dry = 287.05_wp     !: Specific gas constant for dry air              [J/K/kg]
   REAL(wp), PARAMETER ::   R_vap = 461.495_wp    !: Specific gas constant for water vapor          [J/K/kg]
   REAL(wp), PARAMETER ::   reps0 = R_dry/R_vap   !: ratio of gas constant for dry air and water vapor => ~ 0.622
   REAL(wp), PARAMETER ::   rctv0 = R_vap/R_dry   !: for virtual temperature (== (1-eps)/eps) => ~ 0.608

   INTEGER , PARAMETER ::   jpfld   =10           ! maximum number of files to read
   INTEGER , PARAMETER ::   jp_wndi = 1           ! index of 10m wind velocity (i-component) (m/s)    at T-point
   INTEGER , PARAMETER ::   jp_wndj = 2           ! index of 10m wind velocity (j-component) (m/s)    at T-point
   INTEGER , PARAMETER ::   jp_tair = 3           ! index of 10m air temperature             (Kelvin)
   INTEGER , PARAMETER ::   jp_humi = 4           ! index of specific humidity               ( % )
   INTEGER , PARAMETER ::   jp_qsr  = 5           ! index of solar heat                      (W/m2)
   INTEGER , PARAMETER ::   jp_qlw  = 6           ! index of Long wave                       (W/m2)
   INTEGER , PARAMETER ::   jp_prec = 7           ! index of total precipitation (rain+snow) (Kg/m2/s)
   INTEGER , PARAMETER ::   jp_snow = 8           ! index of snow (solid prcipitation)       (kg/m2/s)
   INTEGER , PARAMETER ::   jp_slp  = 9           ! index of sea level pressure              (Pa)
   INTEGER , PARAMETER ::   jp_tdif =10           ! index of tau diff associated to HF tau   (N/m2)   at T-point

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf   ! structure of input fields (file informations, fields read)

   !                                             !!! Bulk parameters
   REAL(wp), PARAMETER ::   cpa    = 1000.5         ! specific heat of air (only used for ice fluxes now...)
   REAL(wp), PARAMETER ::   Ls     =    2.839e6     ! latent heat of sublimation
   REAL(wp), PARAMETER ::   Stef   =    5.67e-8     ! Stefan Boltzmann constant
   REAL(wp), PARAMETER ::   Cd_ice =    1.4e-3      ! transfer coefficient over ice
   REAL(wp), PARAMETER ::   albo   =    0.066       ! ocean albedo assumed to be constant
   !
   !                           !!* Namelist namsbc_blk : bulk parameters
   LOGICAL  ::   ln_NCAR        ! "NCAR"      algorithm   (Large and Yeager 2008)
   LOGICAL  ::   ln_COARE_3p0   ! "COARE 3.0" algorithm   (Fairall et al. 2003)
   LOGICAL  ::   ln_COARE_3p5   ! "COARE 3.5" algorithm   (Edson et al. 2013)
   LOGICAL  ::   ln_ECMWF       ! "ECMWF"     algorithm   (IFS cycle 31)
   !
   LOGICAL  ::   ln_taudif      ! logical flag to use the "mean of stress module - module of mean stress" data
   REAL(wp) ::   rn_pfac        ! multiplication factor for precipitation
   REAL(wp) ::   rn_efac        ! multiplication factor for evaporation
   REAL(wp) ::   rn_vfac        ! multiplication factor for ice/ocean velocity in the calculation of wind stress
   REAL(wp) ::   rn_zqt         ! z(q,t) : height of humidity and temperature measurements
   REAL(wp) ::   rn_zu          ! z(u)   : height of wind measurements
!!gm ref namelist initialize it so remove the setting to false below
   LOGICAL  ::   ln_Cd_L12 = .FALSE. !  Modify the drag ice-atm depending on ice concentration (from Lupkes et al. JGR2012)
   LOGICAL  ::   ln_Cd_L15 = .FALSE. !  Modify the drag ice-atm depending on ice concentration (from Lupkes et al. JGR2015)
   !
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   Cd_atm                    ! transfer coefficient for momentum      (tau)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   Ch_atm                    ! transfer coefficient for sensible heat (Q_sens)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   Ce_atm                    ! tansfert coefficient for evaporation   (Q_lat)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   t_zu                      ! air temperature at wind speed height (needed by Lupkes 2015 bulk scheme)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   q_zu                      ! air spec. hum.  at wind speed height (needed by Lupkes 2015 bulk scheme)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   cdn_oce, chn_oce, cen_oce ! needed by Lupkes 2015 bulk scheme

   INTEGER  ::   nblk           ! choice of the bulk algorithm
   !                            ! associated indices:
   INTEGER, PARAMETER ::   np_NCAR      = 1   ! "NCAR" algorithm        (Large and Yeager 2008)
   INTEGER, PARAMETER ::   np_COARE_3p0 = 2   ! "COARE 3.0" algorithm   (Fairall et al. 2003)
   INTEGER, PARAMETER ::   np_COARE_3p5 = 3   ! "COARE 3.5" algorithm   (Edson et al. 2013)
   INTEGER, PARAMETER ::   np_ECMWF     = 4   ! "ECMWF" algorithm       (IFS cycle 31)

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sbcblk.F90 10425 2018-12-19 21:54:16Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sbc_blk_alloc()
      !!-------------------------------------------------------------------
      !!             ***  ROUTINE sbc_blk_alloc ***
      !!-------------------------------------------------------------------
      ALLOCATE( Cd_atm (jpi,jpj), Ch_atm (jpi,jpj), Ce_atm (jpi,jpj), t_zu(jpi,jpj), q_zu(jpi,jpj), &
         &      cdn_oce(jpi,jpj), chn_oce(jpi,jpj), cen_oce(jpi,jpj), STAT=sbc_blk_alloc )
      !
      CALL mpp_sum ( 'sbcblk', sbc_blk_alloc )
      IF( sbc_blk_alloc /= 0 )   CALL ctl_stop( 'STOP', 'sbc_blk_alloc: failed to allocate arrays' )
   END FUNCTION sbc_blk_alloc


   SUBROUTINE sbc_blk_init
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_blk_init  ***
      !!
      !! ** Purpose :   choose and initialize a bulk formulae formulation
      !!
      !! ** Method  : 
      !!
      !!----------------------------------------------------------------------
      INTEGER  ::   ifpr, jfld            ! dummy loop indice and argument
      INTEGER  ::   ios, ierror, ioptio   ! Local integer
      !!
      CHARACTER(len=100)            ::   cn_dir                ! Root directory for location of atmospheric forcing files
      TYPE(FLD_N), DIMENSION(jpfld) ::   slf_i                 ! array of namelist informations on the fields to read
      TYPE(FLD_N) ::   sn_wndi, sn_wndj, sn_humi, sn_qsr       ! informations about the fields to be read
      TYPE(FLD_N) ::   sn_qlw , sn_tair, sn_prec, sn_snow      !       "                        "
      TYPE(FLD_N) ::   sn_slp , sn_tdif                        !       "                        "
      NAMELIST/namsbc_blk/ sn_wndi, sn_wndj, sn_humi, sn_qsr, sn_qlw ,                &   ! input fields
         &                 sn_tair, sn_prec, sn_snow, sn_slp, sn_tdif,                &
         &                 ln_NCAR, ln_COARE_3p0, ln_COARE_3p5, ln_ECMWF,             &   ! bulk algorithm
         &                 cn_dir , ln_taudif, rn_zqt, rn_zu,                         & 
         &                 rn_pfac, rn_efac, rn_vfac, ln_Cd_L12, ln_Cd_L15
      !!---------------------------------------------------------------------
      !
      !                                      ! allocate sbc_blk_core array
      IF( sbc_blk_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'sbc_blk : unable to allocate standard arrays' )
      !
      !                             !** read bulk namelist  
      REWIND( numnam_ref )                !* Namelist namsbc_blk in reference namelist : bulk parameters
      READ  ( numnam_ref, namsbc_blk, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc_blk in reference namelist', lwp )
      !
      REWIND( numnam_cfg )                !* Namelist namsbc_blk in configuration namelist : bulk parameters
      READ  ( numnam_cfg, namsbc_blk, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc_blk in configuration namelist', lwp )
      !
      IF(lwm) WRITE( numond, namsbc_blk )
      !
      !                             !** initialization of the chosen bulk formulae (+ check)
      !                                   !* select the bulk chosen in the namelist and check the choice
                                                               ioptio = 0
      IF( ln_NCAR      ) THEN   ;   nblk =  np_NCAR        ;   ioptio = ioptio + 1   ;   ENDIF
      IF( ln_COARE_3p0 ) THEN   ;   nblk =  np_COARE_3p0   ;   ioptio = ioptio + 1   ;   ENDIF
      IF( ln_COARE_3p5 ) THEN   ;   nblk =  np_COARE_3p5   ;   ioptio = ioptio + 1   ;   ENDIF
      IF( ln_ECMWF     ) THEN   ;   nblk =  np_ECMWF       ;   ioptio = ioptio + 1   ;   ENDIF
      !
      IF( ioptio /= 1 )   CALL ctl_stop( 'sbc_blk_init: Choose one and only one bulk algorithm' )
      !
      IF( ln_dm2dc ) THEN                 !* check: diurnal cycle on Qsr
         IF( sn_qsr%nfreqh /= 24 )   CALL ctl_stop( 'sbc_blk_init: ln_dm2dc=T only with daily short-wave input' )
         IF( sn_qsr%ln_tint ) THEN 
            CALL ctl_warn( 'sbc_blk_init: ln_dm2dc=T daily qsr time interpolation done by sbcdcy module',   &
               &           '              ==> We force time interpolation = .false. for qsr' )
            sn_qsr%ln_tint = .false.
         ENDIF
      ENDIF
      !                                   !* set the bulk structure
      !                                      !- store namelist information in an array
      slf_i(jp_wndi) = sn_wndi   ;   slf_i(jp_wndj) = sn_wndj
      slf_i(jp_qsr ) = sn_qsr    ;   slf_i(jp_qlw ) = sn_qlw
      slf_i(jp_tair) = sn_tair   ;   slf_i(jp_humi) = sn_humi
      slf_i(jp_prec) = sn_prec   ;   slf_i(jp_snow) = sn_snow
      slf_i(jp_slp)  = sn_slp    ;   slf_i(jp_tdif) = sn_tdif
      !
      lhftau = ln_taudif                     !- add an extra field if HF stress is used
      jfld = jpfld - COUNT( (/.NOT.lhftau/) )
      !
      !                                      !- allocate the bulk structure
      ALLOCATE( sf(jfld), STAT=ierror )
      IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_blk_init: unable to allocate sf structure' )
      DO ifpr= 1, jfld
         ALLOCATE( sf(ifpr)%fnow(jpi,jpj,1) )
         IF( slf_i(ifpr)%ln_tint )   ALLOCATE( sf(ifpr)%fdta(jpi,jpj,1,2) )
         IF( slf_i(ifpr)%nfreqh > 0. .AND. MOD( 3600. * slf_i(ifpr)%nfreqh , REAL(nn_fsbc) * rdt) /= 0. )   &
            &  CALL ctl_warn( 'sbc_blk_init: sbcmod timestep rdt*nn_fsbc is NOT a submultiple of atmospheric forcing frequency.',   &
            &                 '               This is not ideal. You should consider changing either rdt or nn_fsbc value...' )

      END DO
      !                                      !- fill the bulk structure with namelist informations
      CALL fld_fill( sf, slf_i, cn_dir, 'sbc_blk_init', 'surface boundary condition -- bulk formulae', 'namsbc_blk' )
      !
      IF ( ln_wave ) THEN
      !Activated wave module but neither drag nor stokes drift activated
         IF ( .NOT.(ln_cdgw .OR. ln_sdw .OR. ln_tauwoc .OR. ln_stcor ) )   THEN
            CALL ctl_stop( 'STOP',  'Ask for wave coupling but ln_cdgw=F, ln_sdw=F, ln_tauwoc=F, ln_stcor=F' )
      !drag coefficient read from wave model definable only with mfs bulk formulae and core 
         ELSEIF (ln_cdgw .AND. .NOT. ln_NCAR )       THEN       
             CALL ctl_stop( 'drag coefficient read from wave model definable only with NCAR and CORE bulk formulae')
         ELSEIF (ln_stcor .AND. .NOT. ln_sdw)                             THEN
             CALL ctl_stop( 'Stokes-Coriolis term calculated only if activated Stokes Drift ln_sdw=T')
         ENDIF
      ELSE
      IF ( ln_cdgw .OR. ln_sdw .OR. ln_tauwoc .OR. ln_stcor )                & 
         &   CALL ctl_stop( 'Not Activated Wave Module (ln_wave=F) but asked coupling ',    &
         &                  'with drag coefficient (ln_cdgw =T) '  ,                        &
         &                  'or Stokes Drift (ln_sdw=T) ' ,                                 &
         &                  'or ocean stress modification due to waves (ln_tauwoc=T) ',      &  
         &                  'or Stokes-Coriolis term (ln_stcori=T)'  )
      ENDIF 
      !
      !           
      IF(lwp) THEN                     !** Control print
         !
         WRITE(numout,*)                  !* namelist 
         WRITE(numout,*) '   Namelist namsbc_blk (other than data information):'
         WRITE(numout,*) '      "NCAR"      algorithm   (Large and Yeager 2008)     ln_NCAR      = ', ln_NCAR
         WRITE(numout,*) '      "COARE 3.0" algorithm   (Fairall et al. 2003)       ln_COARE_3p0 = ', ln_COARE_3p0
         WRITE(numout,*) '      "COARE 3.5" algorithm   (Edson et al. 2013)         ln_COARE_3p5 = ', ln_COARE_3p0
         WRITE(numout,*) '      "ECMWF"     algorithm   (IFS cycle 31)              ln_ECMWF     = ', ln_ECMWF
         WRITE(numout,*) '      add High freq.contribution to the stress module     ln_taudif    = ', ln_taudif
         WRITE(numout,*) '      Air temperature and humidity reference height (m)   rn_zqt       = ', rn_zqt
         WRITE(numout,*) '      Wind vector reference height (m)                    rn_zu        = ', rn_zu
         WRITE(numout,*) '      factor applied on precipitation (total & snow)      rn_pfac      = ', rn_pfac
         WRITE(numout,*) '      factor applied on evaporation                       rn_efac      = ', rn_efac
         WRITE(numout,*) '      factor applied on ocean/ice velocity                rn_vfac      = ', rn_vfac
         WRITE(numout,*) '         (form absolute (=0) to relative winds(=1))'
         WRITE(numout,*) '      use ice-atm drag from Lupkes2012                    ln_Cd_L12    = ', ln_Cd_L12
         WRITE(numout,*) '      use ice-atm drag from Lupkes2015                    ln_Cd_L15    = ', ln_Cd_L15
         !
         WRITE(numout,*)
         SELECT CASE( nblk )              !* Print the choice of bulk algorithm
         CASE( np_NCAR      )   ;   WRITE(numout,*) '   ==>>>   "NCAR" algorithm        (Large and Yeager 2008)'
         CASE( np_COARE_3p0 )   ;   WRITE(numout,*) '   ==>>>   "COARE 3.0" algorithm   (Fairall et al. 2003)'
         CASE( np_COARE_3p5 )   ;   WRITE(numout,*) '   ==>>>   "COARE 3.5" algorithm   (Edson et al. 2013)'
         CASE( np_ECMWF     )   ;   WRITE(numout,*) '   ==>>>   "ECMWF" algorithm       (IFS cycle 31)'
         END SELECT
         !
      ENDIF
      !
   END SUBROUTINE sbc_blk_init


   SUBROUTINE sbc_blk( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_blk  ***
      !!
      !! ** Purpose :   provide at each time step the surface ocean fluxes
      !!              (momentum, heat, freshwater and runoff)
      !!
      !! ** Method  : (1) READ each fluxes in NetCDF files:
      !!      the 10m wind velocity (i-component) (m/s)    at T-point
      !!      the 10m wind velocity (j-component) (m/s)    at T-point
      !!      the 10m or 2m specific humidity     ( % )
      !!      the solar heat                      (W/m2)
      !!      the Long wave                       (W/m2)
      !!      the 10m or 2m air temperature       (Kelvin)
      !!      the total precipitation (rain+snow) (Kg/m2/s)
      !!      the snow (solid prcipitation)       (kg/m2/s)
      !!      the tau diff associated to HF tau   (N/m2)   at T-point   (ln_taudif=T)
      !!              (2) CALL blk_oce
      !!
      !!      C A U T I O N : never mask the surface stress fields
      !!                      the stress is assumed to be in the (i,j) mesh referential
      !!
      !! ** Action  :   defined at each time-step at the air-sea interface
      !!              - utau, vtau  i- and j-component of the wind stress
      !!              - taum        wind stress module at T-point
      !!              - wndm        wind speed  module at T-point over free ocean or leads in presence of sea-ice
      !!              - qns, qsr    non-solar and solar heat fluxes
      !!              - emp         upward mass flux (evapo. - precip.)
      !!              - sfx         salt flux due to freezing/melting (non-zero only if ice is present)
      !!
      !! ** References :   Large & Yeager, 2004 / Large & Yeager, 2008
      !!                   Brodeau et al. Ocean Modelling 2010
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !!---------------------------------------------------------------------
      !
      CALL fld_read( kt, nn_fsbc, sf )             ! input fields provided at the current time-step
      !
      !                                            ! compute the surface ocean fluxes using bulk formulea
      IF( MOD( kt - 1, nn_fsbc ) == 0 )   CALL blk_oce( kt, sf, sst_m, ssu_m, ssv_m )

#if defined key_cice
      IF( MOD( kt - 1, nn_fsbc ) == 0 )   THEN
         qlw_ice(:,:,1)   = sf(jp_qlw )%fnow(:,:,1)
         IF( ln_dm2dc ) THEN ; qsr_ice(:,:,1) = sbc_dcy( sf(jp_qsr)%fnow(:,:,1) )
         ELSE                ; qsr_ice(:,:,1) =          sf(jp_qsr)%fnow(:,:,1) 
         ENDIF 
         tatm_ice(:,:)    = sf(jp_tair)%fnow(:,:,1)
         qatm_ice(:,:)    = sf(jp_humi)%fnow(:,:,1)
         tprecip(:,:)     = sf(jp_prec)%fnow(:,:,1) * rn_pfac
         sprecip(:,:)     = sf(jp_snow)%fnow(:,:,1) * rn_pfac
         wndi_ice(:,:)    = sf(jp_wndi)%fnow(:,:,1)
         wndj_ice(:,:)    = sf(jp_wndj)%fnow(:,:,1)
      ENDIF
#endif
      !
   END SUBROUTINE sbc_blk


   SUBROUTINE blk_oce( kt, sf, pst, pu, pv )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_oce  ***
      !!
      !! ** Purpose :   provide the momentum, heat and freshwater fluxes at
      !!      the ocean surface at each time step
      !!
      !! ** Method  :   bulk formulea for the ocean using atmospheric
      !!      fields read in sbc_read
      !!
      !! ** Outputs : - utau    : i-component of the stress at U-point  (N/m2)
      !!              - vtau    : j-component of the stress at V-point  (N/m2)
      !!              - taum    : Wind stress module at T-point         (N/m2)
      !!              - wndm    : Wind speed module at T-point          (m/s)
      !!              - qsr     : Solar heat flux over the ocean        (W/m2)
      !!              - qns     : Non Solar heat flux over the ocean    (W/m2)
      !!              - emp     : evaporation minus precipitation       (kg/m2/s)
      !!
      !!  ** Nota  :   sf has to be a dummy argument for AGRIF on NEC
      !!---------------------------------------------------------------------
      INTEGER  , INTENT(in   )                 ::   kt    ! time step index
      TYPE(fld), INTENT(inout), DIMENSION(:)   ::   sf    ! input data
      REAL(wp) , INTENT(in)   , DIMENSION(:,:) ::   pst   ! surface temperature                      [Celcius]
      REAL(wp) , INTENT(in)   , DIMENSION(:,:) ::   pu    ! surface current at U-point (i-component) [m/s]
      REAL(wp) , INTENT(in)   , DIMENSION(:,:) ::   pv    ! surface current at V-point (j-component) [m/s]
      !
      INTEGER  ::   ji, jj               ! dummy loop indices
      REAL(wp) ::   zztmp                ! local variable
      REAL(wp), DIMENSION(jpi,jpj) ::   zwnd_i, zwnd_j    ! wind speed components at T-point
      REAL(wp), DIMENSION(jpi,jpj) ::   zsq               ! specific humidity at pst
      REAL(wp), DIMENSION(jpi,jpj) ::   zqlw, zqsb        ! long wave and sensible heat fluxes
      REAL(wp), DIMENSION(jpi,jpj) ::   zqla, zevap       ! latent heat fluxes and evaporation
      REAL(wp), DIMENSION(jpi,jpj) ::   zst               ! surface temperature in Kelvin
      REAL(wp), DIMENSION(jpi,jpj) ::   zU_zu             ! bulk wind speed at height zu  [m/s]
      REAL(wp), DIMENSION(jpi,jpj) ::   ztpot             ! potential temperature of air at z=rn_zqt [K]
      REAL(wp), DIMENSION(jpi,jpj) ::   zrhoa             ! density of air   [kg/m^3]
      !!---------------------------------------------------------------------
      !
      ! local scalars ( place there for vector optimisation purposes)
      zst(:,:) = pst(:,:) + rt0      ! convert SST from Celcius to Kelvin (and set minimum value far above 0 K)

      ! ----------------------------------------------------------------------------- !
      !      0   Wind components and module at T-point relative to the moving ocean   !
      ! ----------------------------------------------------------------------------- !

      ! ... components ( U10m - U_oce ) at T-point (unmasked)
!!gm    move zwnd_i (_j) set to zero  inside the key_cyclone ???
      zwnd_i(:,:) = 0._wp
      zwnd_j(:,:) = 0._wp
#if defined key_cyclone
      CALL wnd_cyc( kt, zwnd_i, zwnd_j )    ! add analytical tropical cyclone (Vincent et al. JGR 2012)
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vect. opt.
            sf(jp_wndi)%fnow(ji,jj,1) = sf(jp_wndi)%fnow(ji,jj,1) + zwnd_i(ji,jj)
            sf(jp_wndj)%fnow(ji,jj,1) = sf(jp_wndj)%fnow(ji,jj,1) + zwnd_j(ji,jj)
         END DO
      END DO
#endif
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vect. opt.
            zwnd_i(ji,jj) = (  sf(jp_wndi)%fnow(ji,jj,1) - rn_vfac * 0.5 * ( pu(ji-1,jj  ) + pu(ji,jj) )  )
            zwnd_j(ji,jj) = (  sf(jp_wndj)%fnow(ji,jj,1) - rn_vfac * 0.5 * ( pv(ji  ,jj-1) + pv(ji,jj) )  )
         END DO
      END DO
      CALL lbc_lnk_multi( 'sbcblk', zwnd_i, 'T', -1., zwnd_j, 'T', -1. )
      ! ... scalar wind ( = | U10m - U_oce | ) at T-point (masked)
      wndm(:,:) = SQRT(  zwnd_i(:,:) * zwnd_i(:,:)   &
         &             + zwnd_j(:,:) * zwnd_j(:,:)  ) * tmask(:,:,1)

      ! ----------------------------------------------------------------------------- !
      !      I   Radiative FLUXES                                                     !
      ! ----------------------------------------------------------------------------- !

      ! ocean albedo assumed to be constant + modify now Qsr to include the diurnal cycle                    ! Short Wave
      zztmp = 1. - albo
      IF( ln_dm2dc ) THEN   ;   qsr(:,:) = zztmp * sbc_dcy( sf(jp_qsr)%fnow(:,:,1) ) * tmask(:,:,1)
      ELSE                  ;   qsr(:,:) = zztmp *          sf(jp_qsr)%fnow(:,:,1)   * tmask(:,:,1)
      ENDIF

      zqlw(:,:) = (  sf(jp_qlw)%fnow(:,:,1) - Stef * zst(:,:)*zst(:,:)*zst(:,:)*zst(:,:)  ) * tmask(:,:,1)   ! Long  Wave

      ! ----------------------------------------------------------------------------- !
      !     II    Turbulent FLUXES                                                    !
      ! ----------------------------------------------------------------------------- !

      ! ... specific humidity at SST and IST tmask(
      zsq(:,:) = 0.98 * q_sat( zst(:,:), sf(jp_slp)%fnow(:,:,1) )
      !!
      !! Estimate of potential temperature at z=rn_zqt, based on adiabatic lapse-rate
      !!    (see Josey, Gulev & Yu, 2013) / doi=10.1016/B978-0-12-391851-2.00005-2
      !!    (since reanalysis products provide T at z, not theta !)
      ztpot = sf(jp_tair)%fnow(:,:,1) + gamma_moist( sf(jp_tair)%fnow(:,:,1), sf(jp_humi)%fnow(:,:,1) ) * rn_zqt

      SELECT CASE( nblk )        !==  transfer coefficients  ==!   Cd, Ch, Ce at T-point
      !
      CASE( np_NCAR      )   ;   CALL turb_ncar    ( rn_zqt, rn_zu, zst, ztpot, zsq, sf(jp_humi)%fnow, wndm,   &  ! NCAR-COREv2
         &                                           Cd_atm, Ch_atm, Ce_atm, t_zu, q_zu, zU_zu, cdn_oce, chn_oce, cen_oce )
      CASE( np_COARE_3p0 )   ;   CALL turb_coare   ( rn_zqt, rn_zu, zst, ztpot, zsq, sf(jp_humi)%fnow, wndm,   &  ! COARE v3.0
         &                                           Cd_atm, Ch_atm, Ce_atm, t_zu, q_zu, zU_zu, cdn_oce, chn_oce, cen_oce )
      CASE( np_COARE_3p5 )   ;   CALL turb_coare3p5( rn_zqt, rn_zu, zst, ztpot, zsq, sf(jp_humi)%fnow, wndm,   &  ! COARE v3.5
         &                                           Cd_atm, Ch_atm, Ce_atm, t_zu, q_zu, zU_zu, cdn_oce, chn_oce, cen_oce )
      CASE( np_ECMWF     )   ;   CALL turb_ecmwf   ( rn_zqt, rn_zu, zst, ztpot, zsq, sf(jp_humi)%fnow, wndm,   &  ! ECMWF
         &                                           Cd_atm, Ch_atm, Ce_atm, t_zu, q_zu, zU_zu, cdn_oce, chn_oce, cen_oce )
      CASE DEFAULT
         CALL ctl_stop( 'STOP', 'sbc_oce: non-existing bulk formula selected' )
      END SELECT

      !                          ! Compute true air density :
      IF( ABS(rn_zu - rn_zqt) > 0.01 ) THEN     ! At zu: (probably useless to remove zrho*grav*rn_zu from SLP...)
         zrhoa(:,:) = rho_air( t_zu(:,:)              , q_zu(:,:)              , sf(jp_slp)%fnow(:,:,1) )
      ELSE                                      ! At zt:
         zrhoa(:,:) = rho_air( sf(jp_tair)%fnow(:,:,1), sf(jp_humi)%fnow(:,:,1), sf(jp_slp)%fnow(:,:,1) )
      END IF

!!      CALL iom_put( "Cd_oce", Cd_atm)  ! output value of pure ocean-atm. transfer coef.
!!      CALL iom_put( "Ch_oce", Ch_atm)  ! output value of pure ocean-atm. transfer coef.

      DO jj = 1, jpj             ! tau module, i and j component
         DO ji = 1, jpi
            zztmp = zrhoa(ji,jj)  * zU_zu(ji,jj) * Cd_atm(ji,jj)   ! using bulk wind speed
            taum  (ji,jj) = zztmp * wndm  (ji,jj)
            zwnd_i(ji,jj) = zztmp * zwnd_i(ji,jj)
            zwnd_j(ji,jj) = zztmp * zwnd_j(ji,jj)
         END DO
      END DO

      !                          ! add the HF tau contribution to the wind stress module
      IF( lhftau )   taum(:,:) = taum(:,:) + sf(jp_tdif)%fnow(:,:,1)

      CALL iom_put( "taum_oce", taum )   ! output wind stress module

      ! ... utau, vtau at U- and V_points, resp.
      !     Note the use of 0.5*(2-umask) in order to unmask the stress along coastlines
      !     Note the use of MAX(tmask(i,j),tmask(i+1,j) is to mask tau over ice shelves
      DO jj = 1, jpjm1
         DO ji = 1, fs_jpim1
            utau(ji,jj) = 0.5 * ( 2. - umask(ji,jj,1) ) * ( zwnd_i(ji,jj) + zwnd_i(ji+1,jj  ) ) &
               &          * MAX(tmask(ji,jj,1),tmask(ji+1,jj,1))
            vtau(ji,jj) = 0.5 * ( 2. - vmask(ji,jj,1) ) * ( zwnd_j(ji,jj) + zwnd_j(ji  ,jj+1) ) &
               &          * MAX(tmask(ji,jj,1),tmask(ji,jj+1,1))
         END DO
      END DO
      CALL lbc_lnk_multi( 'sbcblk', utau, 'U', -1., vtau, 'V', -1. )

      !  Turbulent fluxes over ocean
      ! -----------------------------

      ! zqla used as temporary array, for rho*U (common term of bulk formulae):
      zqla(:,:) = zrhoa(:,:) * zU_zu(:,:) * tmask(:,:,1)

      IF( ABS( rn_zu - rn_zqt) < 0.01_wp ) THEN
         !! q_air and t_air are given at 10m (wind reference height)
         zevap(:,:) = rn_efac*MAX( 0._wp,             zqla(:,:)*Ce_atm(:,:)*(zsq(:,:) - sf(jp_humi)%fnow(:,:,1)) ) ! Evaporation, using bulk wind speed
         zqsb (:,:) = cp_air(sf(jp_humi)%fnow(:,:,1))*zqla(:,:)*Ch_atm(:,:)*(zst(:,:) - ztpot(:,:)             )   ! Sensible Heat, using bulk wind speed
      ELSE
         !! q_air and t_air are not given at 10m (wind reference height)
         ! Values of temp. and hum. adjusted to height of wind during bulk algorithm iteration must be used!!!
         zevap(:,:) = rn_efac*MAX( 0._wp,             zqla(:,:)*Ce_atm(:,:)*(zsq(:,:) - q_zu(:,:) ) ) ! Evaporation, using bulk wind speed
         zqsb (:,:) = cp_air(sf(jp_humi)%fnow(:,:,1))*zqla(:,:)*Ch_atm(:,:)*(zst(:,:) - t_zu(:,:) )   ! Sensible Heat, using bulk wind speed
      ENDIF

      zqla(:,:) = L_vap(zst(:,:)) * zevap(:,:)     ! Latent Heat flux


      IF(ln_ctl) THEN
         CALL prt_ctl( tab2d_1=zqla  , clinfo1=' blk_oce: zqla   : ', tab2d_2=Ce_atm , clinfo2=' Ce_oce  : ' )
         CALL prt_ctl( tab2d_1=zqsb  , clinfo1=' blk_oce: zqsb   : ', tab2d_2=Ch_atm , clinfo2=' Ch_oce  : ' )
         CALL prt_ctl( tab2d_1=zqlw  , clinfo1=' blk_oce: zqlw   : ', tab2d_2=qsr, clinfo2=' qsr : ' )
         CALL prt_ctl( tab2d_1=zsq   , clinfo1=' blk_oce: zsq    : ', tab2d_2=zst, clinfo2=' zst : ' )
         CALL prt_ctl( tab2d_1=utau  , clinfo1=' blk_oce: utau   : ', mask1=umask,   &
            &          tab2d_2=vtau  , clinfo2=           ' vtau : ', mask2=vmask )
         CALL prt_ctl( tab2d_1=wndm  , clinfo1=' blk_oce: wndm   : ')
         CALL prt_ctl( tab2d_1=zst   , clinfo1=' blk_oce: zst    : ')
      ENDIF

      ! ----------------------------------------------------------------------------- !
      !     III    Total FLUXES                                                       !
      ! ----------------------------------------------------------------------------- !
      !
      emp (:,:) = (  zevap(:,:)                                          &   ! mass flux (evap. - precip.)
         &         - sf(jp_prec)%fnow(:,:,1) * rn_pfac  ) * tmask(:,:,1)
      !
      qns(:,:) = zqlw(:,:) - zqsb(:,:) - zqla(:,:)                                &   ! Downward Non Solar
         &     - sf(jp_snow)%fnow(:,:,1) * rn_pfac * rLfus                        &   ! remove latent melting heat for solid precip
         &     - zevap(:,:) * pst(:,:) * rcp                                      &   ! remove evap heat content at SST
         &     + ( sf(jp_prec)%fnow(:,:,1) - sf(jp_snow)%fnow(:,:,1) ) * rn_pfac  &   ! add liquid precip heat content at Tair
         &     * ( sf(jp_tair)%fnow(:,:,1) - rt0 ) * rcp                          &
         &     + sf(jp_snow)%fnow(:,:,1) * rn_pfac                                &   ! add solid  precip heat content at min(Tair,Tsnow)
         &     * ( MIN( sf(jp_tair)%fnow(:,:,1), rt0 ) - rt0 ) * rcpi
      qns(:,:) = qns(:,:) * tmask(:,:,1)
      !
#if defined key_si3
      qns_oce(:,:) = zqlw(:,:) - zqsb(:,:) - zqla(:,:)                                ! non solar without emp (only needed by SI3)
      qsr_oce(:,:) = qsr(:,:)
#endif
      !
      IF ( nn_ice == 0 ) THEN
         CALL iom_put( "qlw_oce" ,   zqlw )                 ! output downward longwave heat over the ocean
         CALL iom_put( "qsb_oce" , - zqsb )                 ! output downward sensible heat over the ocean
         CALL iom_put( "qla_oce" , - zqla )                 ! output downward latent   heat over the ocean
         CALL iom_put( "qemp_oce",   qns-zqlw+zqsb+zqla )   ! output downward heat content of E-P over the ocean
         CALL iom_put( "qns_oce" ,   qns  )                 ! output downward non solar heat over the ocean
         CALL iom_put( "qsr_oce" ,   qsr  )                 ! output downward solar heat over the ocean
         CALL iom_put( "qt_oce"  ,   qns+qsr )              ! output total downward heat over the ocean
         tprecip(:,:) = sf(jp_prec)%fnow(:,:,1) * rn_pfac * tmask(:,:,1) ! output total precipitation [kg/m2/s]
         sprecip(:,:) = sf(jp_snow)%fnow(:,:,1) * rn_pfac * tmask(:,:,1) ! output solid precipitation [kg/m2/s]
         CALL iom_put( 'snowpre', sprecip )                 ! Snow
         CALL iom_put( 'precip' , tprecip )                 ! Total precipitation
      ENDIF
      !
      IF(ln_ctl) THEN
         CALL prt_ctl(tab2d_1=zqsb , clinfo1=' blk_oce: zqsb   : ', tab2d_2=zqlw , clinfo2=' zqlw  : ')
         CALL prt_ctl(tab2d_1=zqla , clinfo1=' blk_oce: zqla   : ', tab2d_2=qsr  , clinfo2=' qsr   : ')
         CALL prt_ctl(tab2d_1=pst  , clinfo1=' blk_oce: pst    : ', tab2d_2=emp  , clinfo2=' emp   : ')
         CALL prt_ctl(tab2d_1=utau , clinfo1=' blk_oce: utau   : ', mask1=umask,   &
            &         tab2d_2=vtau , clinfo2=              ' vtau  : ' , mask2=vmask )
      ENDIF
      !
   END SUBROUTINE blk_oce



   FUNCTION rho_air( ptak, pqa, pslp )
      !!-------------------------------------------------------------------------------
      !!                           ***  FUNCTION rho_air  ***
      !!
      !! ** Purpose : compute density of (moist) air using the eq. of state of the atmosphere
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk) 
      !!-------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   ptak      ! air temperature             [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pqa       ! air specific humidity   [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pslp      ! pressure in                [Pa]
      REAL(wp), DIMENSION(jpi,jpj)             ::   rho_air   ! density of moist air   [kg/m^3]
      !!-------------------------------------------------------------------------------
      !
      rho_air = pslp / (  R_dry*ptak * ( 1._wp + rctv0*pqa )  )
      !
   END FUNCTION rho_air


   FUNCTION cp_air( pqa )
      !!-------------------------------------------------------------------------------
      !!                           ***  FUNCTION cp_air  ***
      !!
      !! ** Purpose : Compute specific heat (Cp) of moist air
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
      !!-------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pqa      ! air specific humidity         [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj)             ::   cp_air   ! specific heat of moist air   [J/K/kg]
      !!-------------------------------------------------------------------------------
      !
      Cp_air = Cp_dry + Cp_vap * pqa
      !
   END FUNCTION cp_air


   FUNCTION q_sat( ptak, pslp )
      !!----------------------------------------------------------------------------------
      !!                           ***  FUNCTION q_sat  ***
      !!
      !! ** Purpose : Specific humidity at saturation in [kg/kg]
      !!              Based on accurate estimate of "e_sat"
      !!              aka saturation water vapor (Goff, 1957)
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   ptak    ! air temperature                       [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pslp    ! sea level atmospheric pressure       [Pa]
      REAL(wp), DIMENSION(jpi,jpj)             ::   q_sat   ! Specific humidity at saturation   [kg/kg]
      !
      INTEGER  ::   ji, jj         ! dummy loop indices
      REAL(wp) ::   ze_sat, ztmp   ! local scalar
      !!----------------------------------------------------------------------------------
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            ztmp = rt0 / ptak(ji,jj)
            !
            ! Vapour pressure at saturation [hPa] : WMO, (Goff, 1957)
            ze_sat = 10.**( 10.79574*(1. - ztmp) - 5.028*LOG10(ptak(ji,jj)/rt0)        &
               &    + 1.50475*10.**(-4)*(1. - 10.**(-8.2969*(ptak(ji,jj)/rt0 - 1.)) )  &
               &    + 0.42873*10.**(-3)*(10.**(4.76955*(1. - ztmp)) - 1.) + 0.78614  )
               !
            q_sat(ji,jj) = reps0 * ze_sat/( 0.01_wp*pslp(ji,jj) - (1._wp - reps0)*ze_sat )   ! 0.01 because SLP is in [Pa]
            !
         END DO
      END DO
      !
   END FUNCTION q_sat


   FUNCTION gamma_moist( ptak, pqa )
      !!----------------------------------------------------------------------------------
      !!                           ***  FUNCTION gamma_moist  ***
      !!
      !! ** Purpose : Compute the moist adiabatic lapse-rate.
      !!     => http://glossary.ametsoc.org/wiki/Moist-adiabatic_lapse_rate
      !!     => http://www.geog.ucsb.edu/~joel/g266_s10/lecture_notes/chapt03/oh10_3_01/oh10_3_01.html
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   ptak          ! air temperature       [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pqa           ! specific humidity [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj)             ::   gamma_moist   ! moist adiabatic lapse-rate
      !
      INTEGER  ::   ji, jj         ! dummy loop indices
      REAL(wp) :: zrv, ziRT        ! local scalar
      !!----------------------------------------------------------------------------------
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            zrv = pqa(ji,jj) / (1. - pqa(ji,jj))
            ziRT = 1. / (R_dry*ptak(ji,jj))    ! 1/RT
            gamma_moist(ji,jj) = grav * ( 1. + rLevap*zrv*ziRT ) / ( Cp_dry + rLevap*rLevap*zrv*reps0*ziRT/ptak(ji,jj) )
         END DO
      END DO
      !
   END FUNCTION gamma_moist


   FUNCTION L_vap( psst )
      !!---------------------------------------------------------------------------------
      !!                           ***  FUNCTION L_vap  ***
      !!
      !! ** Purpose : Compute the latent heat of vaporization of water from temperature
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             ::   L_vap   ! latent heat of vaporization   [J/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   psst   ! water temperature                [K]
      !!----------------------------------------------------------------------------------
      !
      L_vap = (  2.501 - 0.00237 * ( psst(:,:) - rt0)  ) * 1.e6
      !
   END FUNCTION L_vap

#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   blk_ice_tau : provide the air-ice stress
   !!   blk_ice_flx : provide the heat and mass fluxes at air-ice interface
   !!   blk_ice_qcn : provide ice surface temperature and snow/ice conduction flux (emulating JULES coupler)
   !!   Cdn10_Lupkes2012 : Lupkes et al. (2012) air-ice drag
   !!   Cdn10_Lupkes2015 : Lupkes et al. (2015) air-ice drag 
   !!----------------------------------------------------------------------

   SUBROUTINE blk_ice_tau
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_ice_tau  ***
      !!
      !! ** Purpose :   provide the surface boundary condition over sea-ice
      !!
      !! ** Method  :   compute momentum using bulk formulation
      !!                formulea, ice variables and read atmospheric fields.
      !!                NB: ice drag coefficient is assumed to be a constant
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) ::   zwndi_f , zwndj_f, zwnorm_f   ! relative wind module and components at F-point
      REAL(wp) ::   zwndi_t , zwndj_t             ! relative wind components at T-point
      REAL(wp), DIMENSION(jpi,jpj) ::   zrhoa     ! transfer coefficient for momentum      (tau)
      !!---------------------------------------------------------------------
      !
      ! set transfer coefficients to default sea-ice values
      Cd_atm(:,:) = Cd_ice
      Ch_atm(:,:) = Cd_ice
      Ce_atm(:,:) = Cd_ice

      wndm_ice(:,:) = 0._wp      !!gm brutal....

      ! ------------------------------------------------------------ !
      !    Wind module relative to the moving ice ( U10m - U_ice )   !
      ! ------------------------------------------------------------ !
      ! C-grid ice dynamics :   U & V-points (same as ocean)
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vect. opt.
            zwndi_t = (  sf(jp_wndi)%fnow(ji,jj,1) - rn_vfac * 0.5 * ( u_ice(ji-1,jj  ) + u_ice(ji,jj) )  )
            zwndj_t = (  sf(jp_wndj)%fnow(ji,jj,1) - rn_vfac * 0.5 * ( v_ice(ji  ,jj-1) + v_ice(ji,jj) )  )
            wndm_ice(ji,jj) = SQRT( zwndi_t * zwndi_t + zwndj_t * zwndj_t ) * tmask(ji,jj,1)
         END DO
      END DO
      CALL lbc_lnk( 'sbcblk', wndm_ice, 'T',  1. )
      !
      ! Make ice-atm. drag dependent on ice concentration
      IF    ( ln_Cd_L12 ) THEN   ! calculate new drag from Lupkes(2012) equations
         CALL Cdn10_Lupkes2012( Cd_atm )
         Ch_atm(:,:) = Cd_atm(:,:)       ! momentum and heat transfer coef. are considered identical
      ELSEIF( ln_Cd_L15 ) THEN   ! calculate new drag from Lupkes(2015) equations
         CALL Cdn10_Lupkes2015( Cd_atm, Ch_atm ) 
      ENDIF

!!      CALL iom_put( "Cd_ice", Cd_atm)  ! output value of pure ice-atm. transfer coef.
!!      CALL iom_put( "Ch_ice", Ch_atm)  ! output value of pure ice-atm. transfer coef.

      ! local scalars ( place there for vector optimisation purposes)
      ! Computing density of air! Way denser that 1.2 over sea-ice !!!
      zrhoa (:,:) =  rho_air(sf(jp_tair)%fnow(:,:,1), sf(jp_humi)%fnow(:,:,1), sf(jp_slp)%fnow(:,:,1))

      !!gm brutal....
      utau_ice  (:,:) = 0._wp
      vtau_ice  (:,:) = 0._wp
      !!gm end

      ! ------------------------------------------------------------ !
      !    Wind stress relative to the moving ice ( U10m - U_ice )   !
      ! ------------------------------------------------------------ !
      ! C-grid ice dynamics :   U & V-points (same as ocean)
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vect. opt.
            utau_ice(ji,jj) = 0.5 * zrhoa(ji,jj) * Cd_atm(ji,jj) * ( wndm_ice(ji+1,jj  ) + wndm_ice(ji,jj) )            &
               &          * ( 0.5 * (sf(jp_wndi)%fnow(ji+1,jj,1) + sf(jp_wndi)%fnow(ji,jj,1) ) - rn_vfac * u_ice(ji,jj) )
            vtau_ice(ji,jj) = 0.5 * zrhoa(ji,jj) * Cd_atm(ji,jj) * ( wndm_ice(ji,jj+1  ) + wndm_ice(ji,jj) )            &
               &          * ( 0.5 * (sf(jp_wndj)%fnow(ji,jj+1,1) + sf(jp_wndj)%fnow(ji,jj,1) ) - rn_vfac * v_ice(ji,jj) )
         END DO
      END DO
      CALL lbc_lnk_multi( 'sbcblk', utau_ice, 'U', -1., vtau_ice, 'V', -1. )
      !
      !
      IF(ln_ctl) THEN
         CALL prt_ctl(tab2d_1=utau_ice  , clinfo1=' blk_ice: utau_ice : ', tab2d_2=vtau_ice  , clinfo2=' vtau_ice : ')
         CALL prt_ctl(tab2d_1=wndm_ice  , clinfo1=' blk_ice: wndm_ice : ')
      ENDIF
      !
   END SUBROUTINE blk_ice_tau


   SUBROUTINE blk_ice_flx( ptsu, phs, phi, palb )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_ice_flx  ***
      !!
      !! ** Purpose :   provide the heat and mass fluxes at air-ice interface
      !!
      !! ** Method  :   compute heat and freshwater exchanged
      !!                between atmosphere and sea-ice using bulk formulation
      !!                formulea, ice variables and read atmmospheric fields.
      !!
      !! caution : the net upward water flux has with mm/day unit
      !!---------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   ptsu   ! sea ice surface temperature
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   phs    ! snow thickness
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   phi    ! ice thickness
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   palb   ! ice albedo (all skies)
      !!
      INTEGER  ::   ji, jj, jl               ! dummy loop indices
      REAL(wp) ::   zst3                     ! local variable
      REAL(wp) ::   zcoef_dqlw, zcoef_dqla   !   -      -
      REAL(wp) ::   zztmp, z1_rLsub           !   -      -
      REAL(wp) ::   zfr1, zfr2               ! local variables
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   z1_st         ! inverse of surface temperature
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   z_qlw         ! long wave heat flux over ice
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   z_qsb         ! sensible  heat flux over ice
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   z_dqlw        ! long wave heat sensitivity over ice
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   z_dqsb        ! sensible  heat sensitivity over ice
      REAL(wp), DIMENSION(jpi,jpj)     ::   zevap, zsnw   ! evaporation and snw distribution after wind blowing (SI3)
      REAL(wp), DIMENSION(jpi,jpj)     ::   zrhoa
      !!---------------------------------------------------------------------
      !
      zcoef_dqlw = 4.0 * 0.95 * Stef             ! local scalars
      zcoef_dqla = -Ls * 11637800. * (-5897.8)
      !
      zrhoa(:,:) = rho_air( sf(jp_tair)%fnow(:,:,1), sf(jp_humi)%fnow(:,:,1), sf(jp_slp)%fnow(:,:,1) )
      !
      zztmp = 1. / ( 1. - albo )
      WHERE( ptsu(:,:,:) /= 0._wp )   ;   z1_st(:,:,:) = 1._wp / ptsu(:,:,:)
      ELSEWHERE                       ;   z1_st(:,:,:) = 0._wp
      END WHERE
      !                                     ! ========================== !
      DO jl = 1, jpl                        !  Loop over ice categories  !
         !                                  ! ========================== !
         DO jj = 1 , jpj
            DO ji = 1, jpi
               ! ----------------------------!
               !      I   Radiative FLUXES   !
               ! ----------------------------!
               zst3 = ptsu(ji,jj,jl) * ptsu(ji,jj,jl) * ptsu(ji,jj,jl)
               ! Short Wave (sw)
               qsr_ice(ji,jj,jl) = zztmp * ( 1. - palb(ji,jj,jl) ) * qsr(ji,jj)
               ! Long  Wave (lw)
               z_qlw(ji,jj,jl) = 0.95 * ( sf(jp_qlw)%fnow(ji,jj,1) - Stef * ptsu(ji,jj,jl) * zst3 ) * tmask(ji,jj,1)
               ! lw sensitivity
               z_dqlw(ji,jj,jl) = zcoef_dqlw * zst3

               ! ----------------------------!
               !     II    Turbulent FLUXES  !
               ! ----------------------------!

               ! ... turbulent heat fluxes with Ch_atm recalculated in blk_ice_tau
               ! Sensible Heat
               z_qsb(ji,jj,jl) = zrhoa(ji,jj) * cpa * Ch_atm(ji,jj) * wndm_ice(ji,jj) * (ptsu(ji,jj,jl) - sf(jp_tair)%fnow(ji,jj,1))
               ! Latent Heat
               qla_ice(ji,jj,jl) = rn_efac * MAX( 0.e0, zrhoa(ji,jj) * Ls  * Ch_atm(ji,jj) * wndm_ice(ji,jj) *  &
                  &                ( 11637800. * EXP( -5897.8 * z1_st(ji,jj,jl) ) / zrhoa(ji,jj) - sf(jp_humi)%fnow(ji,jj,1) ) )
               ! Latent heat sensitivity for ice (Dqla/Dt)
               IF( qla_ice(ji,jj,jl) > 0._wp ) THEN
                  dqla_ice(ji,jj,jl) = rn_efac * zcoef_dqla * Ch_atm(ji,jj) * wndm_ice(ji,jj) *  &
                     &                 z1_st(ji,jj,jl)*z1_st(ji,jj,jl) * EXP(-5897.8 * z1_st(ji,jj,jl))
               ELSE
                  dqla_ice(ji,jj,jl) = 0._wp
               ENDIF

               ! Sensible heat sensitivity (Dqsb_ice/Dtn_ice)
               z_dqsb(ji,jj,jl) = zrhoa(ji,jj) * cpa * Ch_atm(ji,jj) * wndm_ice(ji,jj)

               ! ----------------------------!
               !     III    Total FLUXES     !
               ! ----------------------------!
               ! Downward Non Solar flux
               qns_ice (ji,jj,jl) =     z_qlw (ji,jj,jl) - z_qsb (ji,jj,jl) - qla_ice (ji,jj,jl)
               ! Total non solar heat flux sensitivity for ice
               dqns_ice(ji,jj,jl) = - ( z_dqlw(ji,jj,jl) + z_dqsb(ji,jj,jl) + dqla_ice(ji,jj,jl) )
            END DO
            !
         END DO
         !
      END DO
      !
      tprecip(:,:) = sf(jp_prec)%fnow(:,:,1) * rn_pfac * tmask(:,:,1)  ! total precipitation [kg/m2/s]
      sprecip(:,:) = sf(jp_snow)%fnow(:,:,1) * rn_pfac * tmask(:,:,1)  ! solid precipitation [kg/m2/s]
      CALL iom_put( 'snowpre', sprecip )                    ! Snow precipitation
      CALL iom_put( 'precip' , tprecip )                    ! Total precipitation

      ! --- evaporation --- !
      z1_rLsub = 1._wp / rLsub
      evap_ice (:,:,:) = rn_efac * qla_ice (:,:,:) * z1_rLsub    ! sublimation
      devap_ice(:,:,:) = rn_efac * dqla_ice(:,:,:) * z1_rLsub    ! d(sublimation)/dT
      zevap    (:,:)   = rn_efac * ( emp(:,:) + tprecip(:,:) )   ! evaporation over ocean

      ! --- evaporation minus precipitation --- !
      zsnw(:,:) = 0._wp
      CALL ice_thd_snwblow( (1.-at_i_b(:,:)), zsnw )  ! snow distribution over ice after wind blowing
      emp_oce(:,:) = ( 1._wp - at_i_b(:,:) ) * zevap(:,:) - ( tprecip(:,:) - sprecip(:,:) ) - sprecip(:,:) * (1._wp - zsnw )
      emp_ice(:,:) = SUM( a_i_b(:,:,:) * evap_ice(:,:,:), dim=3 ) - sprecip(:,:) * zsnw
      emp_tot(:,:) = emp_oce(:,:) + emp_ice(:,:)

      ! --- heat flux associated with emp --- !
      qemp_oce(:,:) = - ( 1._wp - at_i_b(:,:) ) * zevap(:,:) * sst_m(:,:) * rcp                  & ! evap at sst
         &          + ( tprecip(:,:) - sprecip(:,:) ) * ( sf(jp_tair)%fnow(:,:,1) - rt0 ) * rcp  & ! liquid precip at Tair
         &          +   sprecip(:,:) * ( 1._wp - zsnw ) *                                        & ! solid precip at min(Tair,Tsnow)
         &              ( ( MIN( sf(jp_tair)%fnow(:,:,1), rt0 ) - rt0 ) * rcpi * tmask(:,:,1) - rLfus )
      qemp_ice(:,:) =   sprecip(:,:) * zsnw *                                                    & ! solid precip (only)
         &              ( ( MIN( sf(jp_tair)%fnow(:,:,1), rt0 ) - rt0 ) * rcpi * tmask(:,:,1) - rLfus )

      ! --- total solar and non solar fluxes --- !
      qns_tot(:,:) = ( 1._wp - at_i_b(:,:) ) * qns_oce(:,:) + SUM( a_i_b(:,:,:) * qns_ice(:,:,:), dim=3 )  &
         &           + qemp_ice(:,:) + qemp_oce(:,:)
      qsr_tot(:,:) = ( 1._wp - at_i_b(:,:) ) * qsr_oce(:,:) + SUM( a_i_b(:,:,:) * qsr_ice(:,:,:), dim=3 )

      ! --- heat content of precip over ice in J/m3 (to be used in 1D-thermo) --- !
      qprec_ice(:,:) = rhos * ( ( MIN( sf(jp_tair)%fnow(:,:,1), rt0 ) - rt0 ) * rcpi * tmask(:,:,1) - rLfus )

      ! --- heat content of evap over ice in W/m2 (to be used in 1D-thermo) ---
      DO jl = 1, jpl
         qevap_ice(:,:,jl) = 0._wp ! should be -evap_ice(:,:,jl)*( ( Tice - rt0 ) * rcpi * tmask(:,:,1) )
         !                         ! But we do not have Tice => consider it at 0degC => evap=0 
      END DO

      ! --- shortwave radiation transmitted below the surface (W/m2, see Grenfell Maykut 77) --- !
      zfr1 = ( 0.18 * ( 1.0 - cldf_ice ) + 0.35 * cldf_ice )            ! transmission when hi>10cm
      zfr2 = ( 0.82 * ( 1.0 - cldf_ice ) + 0.65 * cldf_ice )            ! zfr2 such that zfr1 + zfr2 to equal 1
      !
      WHERE    ( phs(:,:,:) <= 0._wp .AND. phi(:,:,:) <  0.1_wp )       ! linear decrease from hi=0 to 10cm  
         qtr_ice_top(:,:,:) = qsr_ice(:,:,:) * ( zfr1 + zfr2 * ( 1._wp - phi(:,:,:) * 10._wp ) )
      ELSEWHERE( phs(:,:,:) <= 0._wp .AND. phi(:,:,:) >= 0.1_wp )       ! constant (zfr1) when hi>10cm
         qtr_ice_top(:,:,:) = qsr_ice(:,:,:) * zfr1
      ELSEWHERE                                                         ! zero when hs>0
         qtr_ice_top(:,:,:) = 0._wp 
      END WHERE
      !
      IF(ln_ctl) THEN
         CALL prt_ctl(tab3d_1=qla_ice , clinfo1=' blk_ice: qla_ice  : ', tab3d_2=z_qsb   , clinfo2=' z_qsb    : ', kdim=jpl)
         CALL prt_ctl(tab3d_1=z_qlw   , clinfo1=' blk_ice: z_qlw    : ', tab3d_2=dqla_ice, clinfo2=' dqla_ice : ', kdim=jpl)
         CALL prt_ctl(tab3d_1=z_dqsb  , clinfo1=' blk_ice: z_dqsb   : ', tab3d_2=z_dqlw  , clinfo2=' z_dqlw   : ', kdim=jpl)
         CALL prt_ctl(tab3d_1=dqns_ice, clinfo1=' blk_ice: dqns_ice : ', tab3d_2=qsr_ice , clinfo2=' qsr_ice  : ', kdim=jpl)
         CALL prt_ctl(tab3d_1=ptsu    , clinfo1=' blk_ice: ptsu     : ', tab3d_2=qns_ice , clinfo2=' qns_ice  : ', kdim=jpl)
         CALL prt_ctl(tab2d_1=tprecip , clinfo1=' blk_ice: tprecip  : ', tab2d_2=sprecip , clinfo2=' sprecip  : ')
      ENDIF
      !
   END SUBROUTINE blk_ice_flx
   

   SUBROUTINE blk_ice_qcn( k_virtual_itd, ptsu, ptb, phs, phi )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_ice_qcn  ***
      !!
      !! ** Purpose :   Compute surface temperature and snow/ice conduction flux
      !!                to force sea ice / snow thermodynamics
      !!                in the case JULES coupler is emulated
      !!                
      !! ** Method  :   compute surface energy balance assuming neglecting heat storage
      !!                following the 0-layer Semtner (1976) approach
      !!
      !! ** Outputs : - ptsu    : sea-ice / snow surface temperature (K)
      !!              - qcn_ice : surface inner conduction flux (W/m2)
      !!
      !!---------------------------------------------------------------------
      INTEGER                   , INTENT(in   ) ::   k_virtual_itd   ! single-category option
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   ptsu            ! sea ice / snow surface temperature
      REAL(wp), DIMENSION(:,:)  , INTENT(in   ) ::   ptb             ! sea ice base temperature
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   phs             ! snow thickness
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   phi             ! sea ice thickness
      !
      INTEGER , PARAMETER ::   nit = 10                  ! number of iterations
      REAL(wp), PARAMETER ::   zepsilon = 0.1_wp         ! characteristic thickness for enhanced conduction
      !
      INTEGER  ::   ji, jj, jl           ! dummy loop indices
      INTEGER  ::   iter                 ! local integer
      REAL(wp) ::   zfac, zfac2, zfac3   ! local scalars
      REAL(wp) ::   zkeff_h, ztsu, ztsu0 !
      REAL(wp) ::   zqc, zqnet           !
      REAL(wp) ::   zhe, zqa0            !
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zgfac   ! enhanced conduction factor
      !!---------------------------------------------------------------------
      
      ! -------------------------------------!
      !      I   Enhanced conduction factor  !
      ! -------------------------------------!
      ! Emulates the enhancement of conduction by unresolved thin ice (k_virtual_itd = 1/2)
      ! Fichefet and Morales Maqueda, JGR 1997
      !
      zgfac(:,:,:) = 1._wp
      
      SELECT CASE ( k_virtual_itd )
      !
      CASE ( 1 , 2 )
         !
         zfac  = 1._wp /  ( rn_cnd_s + rcnd_i )
         zfac2 = EXP(1._wp) * 0.5_wp * zepsilon
         zfac3 = 2._wp / zepsilon
         !   
         DO jl = 1, jpl                
            DO jj = 1 , jpj
               DO ji = 1, jpi
                  zhe = ( rn_cnd_s * phi(ji,jj,jl) + rcnd_i * phs(ji,jj,jl) ) * zfac                            ! Effective thickness
                  IF( zhe >=  zfac2 )   zgfac(ji,jj,jl) = MIN( 2._wp, 0.5_wp * ( 1._wp + LOG( zhe * zfac3 ) ) ) ! Enhanced conduction factor
               END DO
            END DO
         END DO
         !      
      END SELECT
      
      ! -------------------------------------------------------------!
      !      II   Surface temperature and conduction flux            !
      ! -------------------------------------------------------------!
      !
      zfac = rcnd_i * rn_cnd_s
      !
      DO jl = 1, jpl
         DO jj = 1 , jpj
            DO ji = 1, jpi
               !                    
               zkeff_h = zfac * zgfac(ji,jj,jl) / &                                    ! Effective conductivity of the snow-ice system divided by thickness
                  &      ( rcnd_i * phs(ji,jj,jl) + rn_cnd_s * MAX( 0.01, phi(ji,jj,jl) ) )
               ztsu    = ptsu(ji,jj,jl)                                                ! Store current iteration temperature
               ztsu0   = ptsu(ji,jj,jl)                                                ! Store initial surface temperature
               zqa0    = qsr_ice(ji,jj,jl) - qtr_ice_top(ji,jj,jl) + qns_ice(ji,jj,jl) ! Net initial atmospheric heat flux
               !
               DO iter = 1, nit     ! --- Iterative loop
                  zqc   = zkeff_h * ( ztsu - ptb(ji,jj) )                              ! Conduction heat flux through snow-ice system (>0 downwards)
                  zqnet = zqa0 + dqns_ice(ji,jj,jl) * ( ztsu - ptsu(ji,jj,jl) ) - zqc  ! Surface energy budget
                  ztsu  = ztsu - zqnet / ( dqns_ice(ji,jj,jl) - zkeff_h )              ! Temperature update
               END DO
               !
               ptsu   (ji,jj,jl) = MIN( rt0, ztsu )
               qcn_ice(ji,jj,jl) = zkeff_h * ( ptsu(ji,jj,jl) - ptb(ji,jj) )
               qns_ice(ji,jj,jl) = qns_ice(ji,jj,jl) + dqns_ice(ji,jj,jl) * ( ptsu(ji,jj,jl) - ztsu0 )
               qml_ice(ji,jj,jl) = ( qsr_ice(ji,jj,jl) - qtr_ice_top(ji,jj,jl) + qns_ice(ji,jj,jl) - qcn_ice(ji,jj,jl) )  &
                             &   * MAX( 0._wp , SIGN( 1._wp, ptsu(ji,jj,jl) - rt0 ) )

               ! --- Diagnose the heat loss due to changing non-solar flux (as in icethd_zdf_bl99) --- !
               hfx_err_dif(ji,jj) = hfx_err_dif(ji,jj) - ( dqns_ice(ji,jj,jl) * ( ptsu(ji,jj,jl) - ztsu0 ) ) * a_i_b(ji,jj,jl) 

            END DO
         END DO
         !
      END DO 
      !      
   END SUBROUTINE blk_ice_qcn
   

   SUBROUTINE Cdn10_Lupkes2012( Cd )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  Cdn10_Lupkes2012  ***
      !!
      !! ** Purpose :    Recompute the neutral air-ice drag referenced at 10m 
      !!                 to make it dependent on edges at leads, melt ponds and flows.
      !!                 After some approximations, this can be resumed to a dependency
      !!                 on ice concentration.
      !!                
      !! ** Method :     The parameterization is taken from Lupkes et al. (2012) eq.(50)
      !!                 with the highest level of approximation: level4, eq.(59)
      !!                 The generic drag over a cell partly covered by ice can be re-written as follows:
      !!
      !!                 Cd = Cdw * (1-A) + Cdi * A + Ce * (1-A)**(nu+1/(10*beta)) * A**mu
      !!
      !!                    Ce = 2.23e-3       , as suggested by Lupkes (eq. 59)
      !!                    nu = mu = beta = 1 , as suggested by Lupkes (eq. 59)
      !!                    A is the concentration of ice minus melt ponds (if any)
      !!
      !!                 This new drag has a parabolic shape (as a function of A) starting at
      !!                 Cdw(say 1.5e-3) for A=0, reaching 1.97e-3 for A~0.5 
      !!                 and going down to Cdi(say 1.4e-3) for A=1
      !!
      !!                 It is theoretically applicable to all ice conditions (not only MIZ)
      !!                 => see Lupkes et al (2013)
      !!
      !! ** References : Lupkes et al. JGR 2012 (theory)
      !!                 Lupkes et al. GRL 2013 (application to GCM)
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   Cd
      REAL(wp), PARAMETER ::   zCe   = 2.23e-03_wp
      REAL(wp), PARAMETER ::   znu   = 1._wp
      REAL(wp), PARAMETER ::   zmu   = 1._wp
      REAL(wp), PARAMETER ::   zbeta = 1._wp
      REAL(wp)            ::   zcoef
      !!----------------------------------------------------------------------
      zcoef = znu + 1._wp / ( 10._wp * zbeta )

      ! generic drag over a cell partly covered by ice
      !!Cd(:,:) = Cd_oce(:,:) * ( 1._wp - at_i_b(:,:) ) +  &                        ! pure ocean drag
      !!   &      Cd_ice      *           at_i_b(:,:)   +  &                        ! pure ice drag
      !!   &      zCe         * ( 1._wp - at_i_b(:,:) )**zcoef * at_i_b(:,:)**zmu   ! change due to sea-ice morphology

      ! ice-atm drag
      Cd(:,:) = Cd_ice +  &                                                         ! pure ice drag
         &      zCe    * ( 1._wp - at_i_b(:,:) )**zcoef * at_i_b(:,:)**(zmu-1._wp)  ! change due to sea-ice morphology
      
   END SUBROUTINE Cdn10_Lupkes2012


   SUBROUTINE Cdn10_Lupkes2015( Cd, Ch )
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
      !
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   Cd
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   Ch
      REAL(wp), DIMENSION(jpi,jpj)            ::   zst, zqo_sat, zqi_sat
      !
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
      zCdn_form_tmp = zce10 * ( LOG( 10._wp / z0_form_ice + 1._wp ) / LOG( rn_zu / z0_form_ice + 1._wp ) )**2   ! Eq. 40
      zCdn_skin_ice = ( vkarmn                                      / LOG( rn_zu / z0_skin_ice + 1._wp ) )**2   ! Eq. 7
      zCdn_ice      = zCdn_skin_ice   ! Eq. 7 (cf Lupkes email for details)
      !zCdn_ice     = 1.89e-3         ! old ECHAM5 value (cf Eq. 32)

      ! Heat Neutral Transfert Coefficients
      zChn_skin_ice = vkarmn**2 / ( LOG( rn_zu / z0_ice + 1._wp ) * LOG( rn_zu * z1_alpha / z0_skin_ice + 1._wp ) )   ! Eq. 50 + Eq. 52 (cf Lupkes email for details)
     
      ! Atmospheric and Surface Variables
      zst(:,:)     = sst_m(:,:) + rt0                                       ! convert SST from Celcius to Kelvin
      zqo_sat(:,:) = 0.98_wp * q_sat( zst(:,:)  , sf(jp_slp)%fnow(:,:,1) )  ! saturation humidity over ocean [kg/kg]
      zqi_sat(:,:) = 0.98_wp * q_sat( tm_su(:,:), sf(jp_slp)%fnow(:,:,1) )  ! saturation humidity over ice   [kg/kg]
      !
      DO jj = 2, jpjm1           ! reduced loop is necessary for reproducibility
         DO ji = fs_2, fs_jpim1
            ! Virtual potential temperature [K]
            zthetav_os = zst(ji,jj)   * ( 1._wp + rctv0 * zqo_sat(ji,jj) )   ! over ocean
            zthetav_is = tm_su(ji,jj) * ( 1._wp + rctv0 * zqi_sat(ji,jj) )   ! ocean ice
            zthetav_zu = t_zu (ji,jj) * ( 1._wp + rctv0 * q_zu(ji,jj)    )   ! at zu
            
            ! Bulk Richardson Number (could use Ri_bulk function from aerobulk instead)
            zrib_o = grav / zthetav_os * ( zthetav_zu - zthetav_os ) * rn_zu / MAX( 0.5, wndm(ji,jj)     )**2   ! over ocean
            zrib_i = grav / zthetav_is * ( zthetav_zu - zthetav_is ) * rn_zu / MAX( 0.5, wndm_ice(ji,jj) )**2   ! over ice
            
            ! Momentum and Heat Neutral Transfert Coefficients
            zCdn_form_ice = zCdn_form_tmp * at_i_b(ji,jj) * ( 1._wp - at_i_b(ji,jj) )**zbeta  ! Eq. 40
            zChn_form_ice = zCdn_form_ice / ( 1._wp + ( LOG( z1_alphaf ) / vkarmn ) * SQRT( zCdn_form_ice ) )               ! Eq. 53 
                       
            ! Momentum and Heat Stability functions (possibility to use psi_m_ecmwf instead)
            z0w = rn_zu * EXP( -1._wp * vkarmn / SQRT( Cdn_oce(ji,jj) ) ) ! over water
            z0i = z0_skin_ice                                             ! over ice (cf Lupkes email for details)
            IF( zrib_o <= 0._wp ) THEN
               zfmw = 1._wp - zam * zrib_o / ( 1._wp + 3._wp * zc2 * Cdn_oce(ji,jj) * SQRT( -zrib_o * ( rn_zu / z0w + 1._wp ) ) )  ! Eq. 10
               zfhw = ( 1._wp + ( zbetah * ( zthetav_os - zthetav_zu )**r1_3 / ( Chn_oce(ji,jj) * MAX(0.01, wndm(ji,jj)) )   &     ! Eq. 26
                  &             )**zgamma )**z1_gamma
            ELSE
               zfmw = 1._wp / ( 1._wp + zam * zrib_o / SQRT( 1._wp + zrib_o ) )   ! Eq. 12
               zfhw = 1._wp / ( 1._wp + zah * zrib_o / SQRT( 1._wp + zrib_o ) )   ! Eq. 28
            ENDIF
            
            IF( zrib_i <= 0._wp ) THEN
               zfmi = 1._wp - zam * zrib_i / (1._wp + 3._wp * zc2 * zCdn_ice * SQRT( -zrib_i * ( rn_zu / z0i + 1._wp)))   ! Eq.  9
               zfhi = 1._wp - zah * zrib_i / (1._wp + 3._wp * zc2 * zCdn_ice * SQRT( -zrib_i * ( rn_zu / z0i + 1._wp)))   ! Eq. 25
            ELSE
               zfmi = 1._wp / ( 1._wp + zam * zrib_i / SQRT( 1._wp + zrib_i ) )   ! Eq. 11
               zfhi = 1._wp / ( 1._wp + zah * zrib_i / SQRT( 1._wp + zrib_i ) )   ! Eq. 27
            ENDIF
            
            ! Momentum Transfert Coefficients (Eq. 38)
            Cd(ji,jj) = zCdn_skin_ice *   zfmi +  &
               &        zCdn_form_ice * ( zfmi * at_i_b(ji,jj) + zfmw * ( 1._wp - at_i_b(ji,jj) ) ) / MAX( 1.e-06, at_i_b(ji,jj) )
            
            ! Heat Transfert Coefficients (Eq. 49)
            Ch(ji,jj) = zChn_skin_ice *   zfhi +  &
               &        zChn_form_ice * ( zfhi * at_i_b(ji,jj) + zfhw * ( 1._wp - at_i_b(ji,jj) ) ) / MAX( 1.e-06, at_i_b(ji,jj) )
            !
         END DO
      END DO
      CALL lbc_lnk_multi( 'sbcblk', Cd, 'T',  1., Ch, 'T', 1. )
      !
   END SUBROUTINE Cdn10_Lupkes2015

#endif

   !!======================================================================
END MODULE sbcblk
