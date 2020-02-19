! AeroBulk / 2020 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_blk_ice_lg15oi
   !!====================================================================================
   !!       Computes turbulent components of surface fluxes over sea-ice
   !!       Following Lupkes & Gryanik, 2015
   !!
   !!              => case when 100 % sea-ice
   !!
   !!       Routine turb_ice_lg15oi maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, January 2020
   !!
   !!====================================================================================
   USE mod_const       !: physical and other constants
   USE mod_phymbl      !: misc. physical functions 

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: turb_ice_lg15oi

   REAL(wp), PARAMETER ::   rce10_i_0 = 3.46e-3_wp ! (Eq.48) MIZ
   REAL(wp), PARAMETER ::   rbeta_0   = 1.4_wp     ! (Eq.47) MIZ
   REAL(wp), PARAMETER ::   ralpha_0  = 0.2_wp     ! (Eq.12) (ECHAM6 value)
   
   !! To be namelist parameters in NEMO:
   REAL(wp), PARAMETER ::   rz0_s_0  = 0.69e-3_wp  ! Eq. 43 [m]
   REAL(wp), PARAMETER ::   rz0_f_0  = 4.54e-4_wp  ! bottom p.562 MIZ [m]   
   LOGICAL,  PARAMETER :: l_use_form_drag = .TRUE.
   LOGICAL,  PARAMETER :: l_use_pond_info = .FALSE.

   
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE turb_ice_lg15oi( kt, zt, zu, Ti_s, Tw_s, t_zt, qi_s, qw_s, q_zt, &
      &                        U_zu, frice,                                    &
      &                        Cd, Ch, Ce, t_zu, q_zu, U_blk,                  &
      &                       xz0, xu_star, xL, xUN10 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_ice_lg15oi  ***
      !!
      !! ** Purpose :   Computestransfert coefficients of turbulent surface
      !!                fluxes according
      !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!                Returns the effective bulk wind speed at zu to be used in the bulk formulas
      !!
      !! INPUT :
      !! -------
      !!    *  kt   : current time step (starts at 1)
      !!    *  zt   : height for temperature and spec. hum. of air            [m]
      !!    *  zu   : height for wind speed (usually 10m)                     [m]
      !!    *  Ti_s  : surface temperature of sea-ice                         [K]
      !!    *  Tw_s  : surface temperature of water (sea)                     [K]
      !!    *  t_zt : potential air temperature at zt                         [K]
      !!    *  qi_s  : saturation specific humidity at temp. Ti_s over ice    [kg/kg]
      !!    *  qw_s  : saturation specific humidity at temp. Tw_s over water  [kg/kg]
      !!    *  q_zt : specific humidity of air at zt                          [kg/kg]
      !!    *  U_zu : scalar wind speed at zu                                 [m/s]
      !!    * frice : sea-ice concentration        (fraction)
      !!
      !! OUTPUT :
      !! --------
      !!    *  Cd     : drag coefficient
      !!    *  Ch     : sensible heat coefficient
      !!    *  Ce     : evaporation coefficient
      !!    *  t_zu   : pot. air temperature adjusted at wind height zu       [K]
      !!    *  q_zu   : specific humidity of air        //                    [kg/kg]
      !!    *  U_blk  : bulk wind speed at zu                                 [m/s]
      !!
      !! OPTIONAL OUTPUT:
      !! ----------------
      !!    * xz0         : return the aerodynamic roughness length (integration constant for wind stress) [m]
      !!    * xu_star     : return u* the friction velocity                    [m/s]
      !!    * xL          : return the Obukhov length                          [m]
      !!    * xUN10       : neutral wind speed at 10m                          [m/s]
      !!
      !! ** Author: L. Brodeau, January 2020 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      INTEGER,  INTENT(in )                     :: kt    ! current time step
      REAL(wp), INTENT(in )                     :: zt    ! height for t_zt and q_zt                    [m]
      REAL(wp), INTENT(in )                     :: zu    ! height for U_zu                             [m]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: Ti_s  ! ice surface temperature                [Kelvin]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: Tw_s  ! water surface temperature              [Kelvin]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: t_zt  ! potential air temperature              [Kelvin]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: qi_s  ! sat. spec. hum. at ice/air interface    [kg/kg]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: qw_s  ! sat. spec. hum. at water/air interface  [kg/kg]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: q_zt  ! spec. air humidity at zt               [kg/kg]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: U_zu  ! relative wind module at zu                [m/s]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: frice ! sea-ice concentration        (fraction)
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Cd    ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Ch    ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Ce    ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: t_zu  ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: q_zu  ! spec. humidity adjusted at zu           [kg/kg]
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: U_blk ! bulk wind speed at zu                     [m/s]
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(jpi,jpj) :: xz0  ! Aerodynamic roughness length   [m]
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(jpi,jpj) :: xu_star  ! u*, friction velocity
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(jpi,jpj) :: xL  ! zeta (zu/L)
      REAL(wp), INTENT(out), OPTIONAL, DIMENSION(jpi,jpj) :: xUN10  ! Neutral wind at zu
      !
      INTEGER :: j_itt
      LOGICAL :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !
      REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: xtmp1, xtmp2      ! temporary stuff
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: dt_zu, dq_zu, zt_zu, zq_zu  ! third dimension
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zz0_s, zz0_f, RiB, zCdN_s, zChN_s, zCdN_f, zChN_f ! third dimensions (size=2): 1 => ice, 2 => water

      LOGICAL :: lreturn_z0=.FALSE., lreturn_ustar=.FALSE., lreturn_L=.FALSE., lreturn_UN10=.FALSE.
      CHARACTER(len=40), PARAMETER :: crtnm = 'turb_ice_lg15oi@mod_blk_ice_lg15oi.f90'
      !!----------------------------------------------------------------------------------
      ALLOCATE ( xtmp1(jpi,jpj), xtmp2(jpi,jpj) )
      ALLOCATE ( dt_zu(jpi,jpj,2), dq_zu(jpi,jpj,2), zt_zu(jpi,jpj,2), zq_zu(jpi,jpj,2) )
      ALLOCATE ( zz0_s(jpi,jpj,2), zz0_f(jpi,jpj,2), RiB(jpi,jpj,2), zCdN_s(jpi,jpj,2), zChN_s(jpi,jpj,2), zCdN_f(jpi,jpj,2), zChN_f(jpi,jpj,2) )

      
      IF( PRESENT(xz0) )     lreturn_z0    = .TRUE.
      IF( PRESENT(xu_star) ) lreturn_ustar = .TRUE.
      IF( PRESENT(xL) )      lreturn_L     = .TRUE.
      IF( PRESENT(xUN10) )   lreturn_UN10  = .TRUE.

      l_zt_equal_zu = .FALSE.
      IF( ABS(zu - zt) < 0.01_wp )   l_zt_equal_zu = .TRUE. ! testing "zu == zt" is risky with double precision

      !! Scalar wind speed cannot be below 0.2 m/s
      U_blk = MAX( U_zu, 0.2_wp )
           
      !! First guess of temperature and humidity at height zu:
      t_zu = MAX( t_zt ,   100._wp )   ! who knows what's given on masked-continental regions...
      q_zu = MAX( q_zt , 0.1e-6_wp )   !               "
      
      !! Air-Ice & Air-Sea differences (and we don't want them to be 0!)
      dt_zu(:,:,1) = t_zu - Ti_s    
      dq_zu(:,:,1) = q_zu - qi_s    
      dt_zu(:,:,2) = t_zu - Tw_s    
      dq_zu(:,:,2) = q_zu - qw_s    
      dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6_wp), dt_zu )
      dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9_wp), dq_zu )
      
      !! Very crude first guess:
      Cd(:,:) = rCd_ice
      Ch(:,:) = rCd_ice
      Ce(:,:) = rCd_ice

      !! For skin drag :
      zz0_s(:,:,1) = rz0_s_0        !LOLO/RFI! ! Room for improvement. We use the same z0_skin everywhere (= rz0_s_0)...
      xtmp1(:,:) = LOG( zu / zz0_s(:,:,1) )
      zCdN_s(:,:,1) = vkarmn2 / ( xtmp1(:,:) * xtmp1(:,:) )                          ! (Eq.7)   [ index 1 is for ice, 2 for water ]
      zChN_s(:,:,1) = vkarmn2 / ( xtmp1(:,:) * LOG( zu / (ralpha_0*zz0_s(:,:,1)) ) )     ! (Eq.11,12)  [ "" ]

      !! For form drag in MIZ:
      zz0_f(:,:,:)  = 0._wp
      zCdN_f(:,:,:) = 0._wp
      zChN_f(:,:,:) = 0._wp
      IF ( l_use_form_drag ) THEN
         zz0_f(:,:,1) = rz0_f_0        !LOLO/RFI! ! Room for improvement. We use the same z0_form everywhere !!!
         xtmp1(:,:) = 1._wp / zz0_f(:,:,1)
         xtmp2(:,:) = rce10_i_0 * ( LOG( 10._wp * xtmp1(:,:) ) / LOG( zu * xtmp1(:,:) ) )**2       ! part of (Eq.46)
         zCdN_f(:,:,1) = xtmp2(:,:) * frice(:,:) * (1._wp - frice(:,:))**rbeta_0           ! (Eq.46)  [ index 1 is for ice, 2 for water ]
         zChN_f(:,:,1) = zCdN_f(:,:,1) / ( 1._wp + LOG(1._wp/ralpha_0)/vkarmn * SQRT(zCdN_f(:,:,1)) ) ! (Eq.60,61)   [ "" ]
      END IF
      
      !! ITERATION BLOCK
      DO j_itt = 1, nb_itt
                  
         ! Bulk Richardson Number:
         RiB(:,:,1) = Ri_bulk( zu, Ti_s(:,:), t_zu(:,:), qi_s(:,:), q_zu(:,:), U_blk(:,:) )  ! over ice (index=1)
         
         ! Momentum and Heat transfer coefficients WITHOUT FORM DRAG / (Eq.6) and (Eq.10):
         Cd(:,:) = zCdN_s(:,:,1) * f_m_louis( zu, RiB(:,:,1), zCdN_s(:,:,1), zz0_s(:,:,1) ) ! (Eq.6)
         Ch(:,:) = zChN_s(:,:,1) * f_h_louis( zu, RiB(:,:,1), zCdN_s(:,:,1), zz0_s(:,:,1) ) ! (Eq.10) / LOLO: why "zCdN_s" (xtmp1) and not "zChn" ???
         PRINT *, 'LOLO: Cd / skin only / ice =', REAL(Cd,4)


         IF ( l_use_form_drag ) THEN
            !! Form-drag-related NEUTRAL momentum and Heat transfer coefficients:
            !!   MIZ:
            Cd(:,:) = Cd(:,:) + zCdN_f(:,:,1) * f_m_louis( zu, RiB(:,:,1), zCdN_f(:,:,1), zz0_f(:,:,1) ) ! (Eq.6)
            Ch(:,:) = Ch(:,:) + zChN_f(:,:,1) * f_h_louis( zu, RiB(:,:,1), zCdN_f(:,:,1), zz0_f(:,:,1) ) ! (Eq.10) / LOLO: why "zCdN_f" (xtmp1) and not "zChn" ???
         END IF
         
         Ce(:,:) = Ch(:,:)
         PRINT *, 'LOLO: Cd / total / ice =', REAL(Cd,4)
         
         !! Adjusting temperature and humidity from zt to zu:
         IF( .NOT. l_zt_equal_zu ) THEN

            !! Over ice:
            xtmp1(:,:) = zCdN_s(:,:,1) + zCdN_f(:,:,1)    ! total neutral drag coeff!
            xtmp2(:,:) = zz0_s(:,:,1) + zz0_f(:,:,1)      ! total roughness length z0
            xtmp1 = LOG(zt/zu) + f_h_louis( zu, RiB(:,:,1), xtmp1(:,:), xtmp2(:,:) ) &
               &               - f_h_louis( zt, RiB(:,:,1), xtmp1(:,:), xtmp2(:,:) )
            xtmp2 = 1._wp/SQRT(Cd)
            zt_zu(:,:,1) = t_zt - (Ch * dt_zu(:,:,1) * xtmp2) / vkarmn * xtmp1   ! t_star = Ch * dt_zu / SQRT(Cd)
            zq_zu(:,:,1) = q_zt - (Ce * dq_zu(:,:,1) * xtmp2) / vkarmn * xtmp1   ! q_star = Ce * dq_zu / SQRT(Cd)
            zq_zu(:,:,1) = MAX(0._wp, q_zu)
            dt_zu(:,:,1) = zt_zu(:,:,1) - Ti_s
            dq_zu(:,:,1) = zq_zu(:,:,1) - qi_s
            
            !! Over water:
            xtmp1(:,:) = zCdN_s(:,:,2) + zCdN_f(:,:,2)    ! total neutral drag coeff!
            xtmp2(:,:) = zz0_s(:,:,2) + zz0_f(:,:,2)      ! total roughness length z0
            xtmp1 = LOG(zt/zu) + f_h_louis( zu, RiB(:,:,2), xtmp1(:,:), xtmp2(:,:) ) &
               &               - f_h_louis( zt, RiB(:,:,2), xtmp1(:,:), xtmp2(:,:) )
            xtmp2 = 1._wp/SQRT(Cd)
            zt_zu(:,:,2) = t_zt - (Ch * dt_zu(:,:,2) * xtmp2) / vkarmn * xtmp1   ! t_star = Ch * dt_zu / SQRT(Cd)
            zq_zu(:,:,2) = q_zt - (Ce * dq_zu(:,:,2) * xtmp2) / vkarmn * xtmp1   ! q_star = Ce * dq_zu / SQRT(Cd)
            zq_zu(:,:,2) = MAX(0._wp, q_zu)
            dt_zu(:,:,2) = zt_zu(:,:,2) - Tw_s
            dq_zu(:,:,2) = zq_zu(:,:,2) - qw_s
            
            dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6_wp), dt_zu )
            dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9_wp), dq_zu )
         END IF

         PRINT *, ''!LOLO         
      END DO !DO j_itt = 1, nb_itt


      PRINT *, 'LOLO: MUST combine theta_zu and q_zu for the mix of water + ice over the mesh !!!'
      t_zu(:,:) = mix_val_msh(zt_zu, frice)
      q_zu(:,:) = mix_val_msh(zq_zu, frice)
      
      
      IF( lreturn_z0 ) xz0   = z0_from_Cd( zu, zCdN_s(:,:,1)+zCdN_f(:,:,1) )
      
      IF( lreturn_ustar ) xu_star = SQRT(Cd) * U_blk
      IF( lreturn_L ) THEN
         xtmp1 = SQRT(Cd)
         xL    = 1./One_on_L(t_zu, q_zu, xtmp1*U_blk, Ch*mix_val_msh(dt_zu,frice)/xtmp1, Ce*mix_val_msh(dq_zu,frice)/xtmp1)
      END IF
      
      IF( lreturn_UN10 ) xUN10   = (SQRT(Cd) * U_blk) / vkarmn*LOG(10./(z0_from_Cd( zu, zCdN_s(:,:,1)+zCdN_f(:,:,1) )) )
      
      DEALLOCATE ( xtmp1, xtmp2 )
      DEALLOCATE ( dt_zu, dq_zu, zt_zu, zq_zu )      
      DEALLOCATE ( zz0_s, zz0_f, RiB, zCdN_s, zChN_s, zCdN_f, zChN_f )
      
   END SUBROUTINE turb_ice_lg15oi
   
   !!======================================================================
   
   

   FUNCTION mix_val_msh( pfld, pfri )
      REAL(wp), DIMENSION(jpi,jpj) :: mix_val_msh
      REAL(wp), DIMENSION(jpi,jpj,2), INTENT(in) :: pfld  ! field array to "water/ice average" on the mesh [ over ice=>(:,:,1), over water=>(:,:,2) ]
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in) :: pfri  ! sea-ice concentration (fraction)
      mix_val_msh(:,:) = pfri(:,:) * pfld(:,:,1) + (1._wp - pfri(:,:)) * pfld(:,:,2)
   END FUNCTION mix_val_msh
   
END MODULE mod_blk_ice_lg15oi
