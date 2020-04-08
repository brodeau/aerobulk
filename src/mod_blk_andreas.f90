!!! TO DO: consistent psi_m and psi_h needed!!! For now is those of NCAR !!!
!!
! AeroBulk / 2020 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_blk_andreas
   !!====================================================================================
   !!       Computes turbulent components of surface fluxes
   !!         according to Andreas et al. (2015)
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !!       Andreas, E.L., Mahrt, L. and Vickers, D. (2015),
   !!       An improved bulk air–sea surface flux algorithm,
   !!       including spray‐mediated transfer.
   !!       Q.J.R. Meteorol. Soc., 141: 642-654. doi:10.1002/qj.2424
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !!
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (2m) to zu (10m) if needed
   !!   * the effective bulk wind speed at 10m Ub
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!    Using the bulk formulation/param. of Large & Yeager 2008
   !!
   !!       Routine turb_andreas maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, 2020
   !!
   !!====================================================================================
   USE mod_const                                         !: physical and othe constants
   USE mod_phymbl                                        !: thermodynamics
   USE mod_blk_coare3p0, ONLY: psi_m_coare, psi_h_coare
   !USE mod_blk_ncar    , ONLY: cd_n10_ncar, ch_n10_ncar, ce_n10_ncar

   
   IMPLICIT NONE
   PRIVATE
   
   !REAL(wp), PARAMETER :: zeta_abs_max = 50._wp
   REAL(wp), PARAMETER :: L_min    = 1._wp  ! Limits L to L_min when ultra stable (stable => L > 0)
   
   PUBLIC :: TURB_ANDREAS

   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE turb_andreas( zt, zu, sst, t_zt, ssq, q_zt, U_zu,   &
      &                     Cd, Ch, Ce, t_zu, q_zu, Ub,           &
      &                    CdN, ChN, CeN, xz0, xu_star, xL, xUN10 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_andreas  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to Large & Yeager (2004) and Large & Yeager (2008)
      !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!                Returns the effective bulk wind speed at zu to be used in the bulk formulas
      !!
      !! INPUT :
      !! -------
      !!    *  zt   : height for temperature and spec. hum. of air            [m]
      !!    *  zu   : height for wind speed (usually 10m)                     [m]
      !!    *  sst  : bulk SST                                                [K]
      !!    *  t_zt : potential air temperature at zt                         [K]
      !!    *  ssq  : specific humidity at saturation at SST                  [kg/kg]
      !!    *  q_zt : specific humidity of air at zt                          [kg/kg]
      !!    *  U_zu : scalar wind speed at zu                                 [m/s]
      !!    *  SLP  : sea level pressure (needed if zt /= zu)                 [Pa]
      !!
      !! OUTPUT :
      !! --------
      !!    *  Cd     : drag coefficient
      !!    *  Ch     : sensible heat coefficient
      !!    *  Ce     : evaporation coefficient
      !!    *  t_zu   : pot. air temperature adjusted at wind height zu       [K]
      !!    *  q_zu   : specific humidity of air        //                    [kg/kg]
      !!    *  Ub  : bulk wind speed at zu                                 [m/s]
      !!
      !! OPTIONAL OUTPUT:
      !! ----------------
      !!    * CdN      : neutral-stability drag coefficient
      !!    * ChN      : neutral-stability sensible heat coefficient
      !!    * CeN      : neutral-stability evaporation coefficient
      !!    * xz0      : return the aerodynamic roughness length (integration constant for wind stress) [m]
      !!    * xu_star  : return u* the friction velocity                    [m/s]
      !!    * xL       : return the Obukhov length                          [m]
      !!    * xUN10    : neutral wind speed at 10m                          [m/s]
      !!
      !! ** Author: L. Brodeau, June 2019 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(in   )                     ::   zt       ! height for t_zt and q_zt                    [m]
      REAL(wp), INTENT(in   )                     ::   zu       ! height for U_zu                             [m]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   sst      ! sea surface temperature                [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   t_zt     ! potential air temperature              [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   ssq      ! sea surface specific humidity           [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   q_zt     ! specific air humidity at zt             [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   U_zu     ! relative wind module at zu                [m/s]
      !!
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Cd       ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ch       ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ce       ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   t_zu     ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   q_zu     ! spec. humidity adjusted at zu           [kg/kg]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ub    ! bulk wind speed at zu                     [m/s]
      !!
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   CdN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   ChN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   CeN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xz0     ! Aerodynamic roughness length   [m]
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xu_star ! u*, friction velocity
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xL      ! zeta (zu/L)
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   xUN10   ! Neutral wind at zu
      !!
      INTEGER :: j_itt
      LOGICAL :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !!
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   u_star, t_star, q_star
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   z0       ! roughness length (momentum) [m]
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   UN10     ! Neutral wind speed at zu [m/s]
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   zeta_u        ! stability parameter at height zu
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   ztmp0, ztmp1, ztmp2
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   sqrtCd       ! square root of Cd
      !!
      LOGICAL ::  lreturn_cdn=.FALSE., lreturn_chn=.FALSE., lreturn_cen=.FALSE., &
         &        lreturn_z0=.FALSE., lreturn_ustar=.FALSE., lreturn_L=.FALSE., lreturn_UN10=.FALSE.
      CHARACTER(len=40), PARAMETER :: crtnm = 'turb_andreas@mod_blk_andreas.f90'
      !!----------------------------------------------------------------------------------
      
      ALLOCATE( u_star(jpi,jpj), t_star(jpi,jpj), q_star(jpi,jpj), &
         &          z0(jpi,jpj),   UN10(jpi,jpj), zeta_u(jpi,jpj), sqrtCd(jpi,jpj), &
         &       ztmp0(jpi,jpj),  ztmp1(jpi,jpj),  ztmp2(jpi,jpj)  )

      IF( PRESENT(CdN) )     lreturn_cdn   = .TRUE.
      IF( PRESENT(ChN) )     lreturn_chn   = .TRUE.
      IF( PRESENT(CeN) )     lreturn_cen   = .TRUE.
      IF( PRESENT(xz0) )     lreturn_z0    = .TRUE.
      IF( PRESENT(xu_star) ) lreturn_ustar = .TRUE.
      IF( PRESENT(xL) )      lreturn_L     = .TRUE.
      IF( PRESENT(xUN10) )   lreturn_UN10  = .TRUE.

      l_zt_equal_zu = ( ABS(zu - zt) < 0.01_wp ) ! testing "zu == zt" is risky with double precision
      
      Ub = MAX( 0.25_wp , U_zu ) !  relative bulk wind speed at zu

      !! First guess:
      UN10 = Ub
      Cd   = 1.1E-3_wp
      Ch   = 1.1E-3_wp
      Ce   = 1.1E-3_wp
      t_zu = t_zt
      q_zu = q_zt
      
      !! First guess of turbulent scales for scalars:
      ztmp0  = SQRT(Cd)
      t_star = Ch/ztmp0*(t_zu - sst) ! theta*
      q_star = Ce/ztmp0*(q_zu - ssq) ! q*


      !! ITERATION BLOCK
      DO j_itt = 1, nb_itt
      !DO j_itt = 1, 5

         PRINT *, 'LOLO'

         u_star = U_STAR_ANDREAS( UN10 )

         PRINT *, 'LOLO *** u* =', u_star, j_itt
         PRINT *, 'LOLO *** t_zu =', t_zu, j_itt
         PRINT *, 'LOLO *** q_zu =', q_zu, j_itt
         PRINT *, 'LOLO *** theta* =', t_star, j_itt
         PRINT *, 'LOLO *** q* =', q_star, j_itt

         !! Stability parameter :
         ztmp0 = One_on_L( t_zu, q_zu, u_star, t_star, q_star )   ! 1/L
         ztmp0 = MIN( ztmp0 , 1._wp/L_min )  ! 1/L LOLO: needed WHY ????
         ztmp0 = SIGN( MIN(ABS(ztmp0),200._wp), ztmp0 ) ! 1/L (prevents FPE from stupid values from masked region later on...) 
         zeta_u = zu*ztmp0
         !zeta_u = SIGN( MIN(ABS(zeta_u),zeta_abs_max), zeta_u )

         PRINT *, 'LOLO *** L =', zu/zeta_u, j_itt
         PRINT *, 'LOLO *** zeta_u =', zeta_u, j_itt
         PRINT *, 'LOLO *** Ub =', Ub, j_itt         
         
         
         !! Drag coefficient:
         ztmp0 = u_star/Ub
         Cd    = ztmp0*ztmp0
         PRINT *, 'LOLO *** CD =', Cd, j_itt
         
         !! Roughness length:
         IF( j_itt > 1 ) THEN
            !ztmp1 = MAX( psi_m_coare(zeta_u) , psi_min )   ! Psi_m
            !z0    = z0_from_Cd( zu, Cd,  ppsi=ztmp1 )
            z0    = z0_from_Cd( zu, Cd,  ppsi=psi_m_coare(zeta_u) )
         ELSE
            z0    = z0_from_Cd( zu, Cd )
         END IF

         PRINT *, 'LOLO *** z0 =', z0, j_itt
         
         !! z0t and z0q, based on LKB, just like into COARE 2.5:
         ztmp0 = z0 * u_star / visc_air(t_zu) ! Re_r
         ztmp1 = z0tq_LKB( 1, ztmp0, z0 )    ! z0t
         ztmp2 = z0tq_LKB( 2, ztmp0, z0 )    ! z0q

         !! Turbulent scales at zu :
         !ztmp0   = MAX( psi_h_coare(zeta_u) , psi_min )  ! lolo: zeta_u for scalars???
         ztmp0 = psi_h_coare(zeta_u)  ! lolo: zeta_u for scalars???
         PRINT *, 'LOLO *** psi_h(zeta_u) =', ztmp0, j_itt
         t_star  = (t_zu - sst)*vkarmn/(LOG(zu) - LOG(ztmp1) - ztmp0)  ! theta* (ztmp1 == z0t in rhs term)
         q_star  = (q_zu - ssq)*vkarmn/(LOG(zu) - LOG(ztmp2) - ztmp0)  !   q*   (ztmp2 == z0q in rhs term)

         IF( (.NOT. l_zt_equal_zu).AND.( j_itt > 1 ) ) THEN
            !! Re-updating temperature and humidity at zu if zt /= zu:
            ztmp0 = zeta_u/zu*zt   ! zeta_t
            !ztmp0 = SIGN( MIN(ABS(ztmp0),zeta_abs_max), ztmp0 )  ! zeta_t
            ztmp0 = LOG(zt/zu) + psi_h_coare(zeta_u) - psi_h_coare(ztmp0)
            t_zu = t_zt - t_star/vkarmn*ztmp0
            q_zu = q_zt - q_star/vkarmn*ztmp0
         ENDIF

         !! Update neutral-stability wind at zu:
         UN10 = MAX( 0.25_wp , UN10_from_ustar( zu, Ub, u_star, psi_m_coare(zeta_u) ) ) ! UN10
         PRINT *, 'LOLO *** UN10 =', UN10, j_itt
         PRINT *, 'LOLO'
         
      END DO !DO j_itt = 1, nb_itt

      
      ! compute transfer coefficients at zu:
      ztmp0 = u_star/Ub
      Cd   = ztmp0*ztmp0

      ztmp1 = t_zu - sst ;  ztmp1 = SIGN( MAX(ABS(ztmp1),1.E-6_wp), ztmp1 )  ! dt_zu
      ztmp2 = q_zu - ssq ;  ztmp2 = SIGN( MAX(ABS(ztmp2),1.E-9_wp), ztmp2 )  ! dq_zu
      Ch   = ztmp0*t_star/ztmp1
      Ce   = ztmp0*q_star/ztmp2
      

      IF( lreturn_cdn .OR. lreturn_chn .OR. lreturn_cen ) ztmp0 = 1._wp/LOG(zu/z0)
      IF( lreturn_cdn )   CdN     = vkarmn2*ztmp0*ztmp0
      
      IF( lreturn_chn .OR. lreturn_cen ) ztmp1 = z0 * u_star / visc_air(t_zu)  ! Re_r
      IF( lreturn_chn )   ChN     = vkarmn2*ztmp0/LOG(zu/z0tq_LKB( 1, ztmp1, z0 ))
      IF( lreturn_cen )   CeN     = vkarmn2*ztmp0/LOG(zu/z0tq_LKB( 2, ztmp1, z0 ))

      !IF( lreturn_z0 )    xz0     = z0_from_Cd( zu, Cd,  ppsi=psi_m_coare(zeta_u) )
      IF( lreturn_z0 )    xz0     = z0
      IF( lreturn_ustar ) xu_star = u_star
      IF( lreturn_L )     xL      = zu/zeta_u
      IF( lreturn_UN10 )  xUN10   =  UN10_from_ustar( zu, Ub, u_star, psi_m_coare(zeta_u) )
      

      DEALLOCATE( u_star, t_star, q_star, z0, UN10, zeta_u, sqrtCd, ztmp0, ztmp1, ztmp2 ) !

   END SUBROUTINE turb_andreas


   FUNCTION U_STAR_ANDREAS( pun10 )
      !!----------------------------------------------------------------------------------
      !! Estimate of the friction velocity as a function of the neutral-stability wind
      !! speed at at 10m
      !!
      !! Origin: Eq.(2.2) of Andreas et al. (2015)
      !!
      !! ** Author: L. Brodeau, April 2020 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pun10          !: neutral-stability scalar wind speed at 10m (m/s)
      REAL(wp), DIMENSION(jpi,jpj)             :: u_star_andreas !: friction velocity    [m/s]
      !
      INTEGER  ::     ji, jj ! dummy loop indices
      REAL(wp) :: za, zt, zw ! local scalars
      !!----------------------------------------------------------------------------------
      !
      DO jj = 1, jpj
         DO ji = 1, jpi

            zw  = pun10(ji,jj)

            za = zw - 8.271_wp

            zt = za + SQRT( MAX( 0.12_wp*za*za + 0.181_wp , 0._wp ) )

            u_star_andreas(ji,jj) =   0.239_wp + 0.0433_wp * zt

         END DO
      END DO
      !
   END FUNCTION U_STAR_ANDREAS



   !!======================================================================
END MODULE mod_blk_andreas
