MODULE mod_blk_ice_nanuq
   !!======================================================================
   !!                   ***  MODULE  mod_blk_ice_nanuq  ***
   !!       Computes turbulent components of surface fluxes over sea-ice
   !!
   !!   What is the "NANUQ" algorithm ?
   !!    => Given a constant value for the neutral coefficients C_D_N, C_E_N and C_H_N over sea-ice
   !!       it will commpute the C_D, C_E and C_H consistent with the near-surface atmospheric stability
   !!    ==> which is already a better approach than using a constant value for C_D, C_E and C_H as in
   !!        the default (and simplest) option in NEMO for instance.
   !!    ==> the only room for improvement is the pick of the stability functions: we use
   !!        Paulson 1970 (as in Jordan et al. 1999) for unstable conditions and Grachev et al. 2007 for stable conditions
   !!
   !! What it computes:
   !!   * bulk transfer coefficients C_D, C_E and C_H over sea-ice
   !!   * air temp. and spec. hum. adjusted from zt (usually 2m) to zu (usually 10m) if needed
   !!   * the "effective" bulk wind speed at zu: Ub (including gustiness contribution in unstable conditions)
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!       Routine turb_ice_nanuq maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, Summer 2026
   !!
   !!----------------------------------------------------------------------
   USE mod_const        !: physical and othe constants
   USE mod_phymbl       !: thermodynamics

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: turb_ice_nanuq

   ! Note: the number of iterations `nb_iter` is the one of module `mod_const` (shared by all algorithms,
   ! set by the user), it is not redefined here (a local `nb_iter` would clash with the one of `mod_const`).
   REAL(wp), PARAMETER ::   epsi10 = 1.E-10_wp
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !! NANUQ 1.0.1, Brodeau (2026)
   !! -- a standalone fork of NEMO/SI3 version 5 ---
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE turb_ice_nanuq( zt, zu, pts_i, pt_zt, pqs_i, pq_zt, pU_zu, &
      &                      pCdN, pChN, pCeN,                          &
      &                      pCd_i, pCh_i, pCe_i, pt_zu_i, pq_zu_i )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_ice_nanuq  ***
      !!
      !! ** Purpose :   Computes the turbulent transfer coefficients over sea-ice C_D, C_H and C_E
      !!                from given constant NEUTRAL coefficients (pCdN, pChN, pCeN), corrected for the
      !!                stability of the surface atmospheric layer (Monin-Obukhov theory, following the
      !!                method of Large & Yeager 2004, Eq.(9) and (10)), with the stability functions:
      !!                  * unstable: Paulson 1970, as in Eq.(30) & (31) of Jordan et al. 1999
      !!                  * stable:   Grachev et al. 2007, Eq.(12) & (13)
      !!
      !!           If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!           (the bulk wind speed used is the wind speed at zu, not below `wspd_thrshld_ice`)
      !!
      !! INPUT :
      !! -------
      !!    *  zt   : height for temperature and spec. hum. of air            [m]
      !!    *  zu   : height for wind speed (usually 10m)                     [m]
      !!    *  pts_i  : surface temperature of sea-ice                         [K]
      !!    *  pt_zt : potential air temperature at zt                         [K]
      !!    *  pqs_i  : saturation specific humidity at temp. pts_i over ice    [kg/kg]
      !!    *  pq_zt : specific humidity of air at zt                          [kg/kg]
      !!    *  pU_zu : scalar wind speed at zu                                 [m/s]
      !!    * pCdN     : neutral-stability drag coefficient
      !!    * pChN     : neutral-stability sensible heat coefficient
      !!    * pCeN     : neutral-stability evaporation coefficient
      !!
      !! OUTPUT :
      !! --------
      !!    *  pCd_i   : drag coefficient over sea-ice
      !!    *  pCh_i   : sensible heat coefficient over sea-ice
      !!    *  pCe_i   : sublimation coefficient over sea-ice
      !!    *  pt_zu_i : pot. air temp. adjusted at zu over sea-ice             [K]
      !!    *  pq_zu_i : spec. hum. of air adjusted at zu over sea-ice          [kg/kg]
      !!
      !! ** Author: L. Brodeau, July 2023 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(in )                     :: zt    ! height for pt_zt and pq_zt                    [m]
      REAL(wp), INTENT(in )                     :: zu    ! height for pU_zu                             [m]
      REAL(wp), INTENT(in ), DIMENSION(:,:)     :: pts_i  ! ice surface temperature                [Kelvin]
      REAL(wp), INTENT(in ), DIMENSION(:,:)     :: pt_zt  ! potential air temperature              [Kelvin]
      REAL(wp), INTENT(in ), DIMENSION(:,:)     :: pqs_i  ! sat. spec. hum. at ice/air interface    [kg/kg]
      REAL(wp), INTENT(in ), DIMENSION(:,:)     :: pq_zt  ! spec. air humidity at zt               [kg/kg]
      REAL(wp), INTENT(in ), DIMENSION(:,:)     :: pU_zu  ! absolute wind module at zu (not corrected for ice drift) [m/s]
      REAL(wp), INTENT(in )                     :: pCdN, pChN, pCeN  ! neutral coefficients             [-]
      REAL(wp), INTENT(out), DIMENSION(:,:)     :: pCd_i  ! drag coefficient over sea-ice
      REAL(wp), INTENT(out), DIMENSION(:,:)     :: pCh_i  ! transfert coefficient for heat over ice
      REAL(wp), INTENT(out), DIMENSION(:,:)     :: pCe_i  ! transfert coefficient for sublimation over ice
      REAL(wp), INTENT(out), DIMENSION(:,:)     :: pt_zu_i ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(out), DIMENSION(:,:)     :: pq_zu_i ! spec. humidity adjusted at zu           [kg/kg]
      !!----------------------------------------------------------------------------------
      INTEGER  :: Ni, Nj, ji, jj, jit
      REAL(wp) :: zUbzu, z1_vkarmn
      REAL(wp) :: zdum, zdvm, z1_L, zsqrtCD
      REAL(wp) :: zdt_zu, zdq_zu
      REAL(wp) :: zus, zts, zqs
      REAL(wp) :: zzeta_u, zzeta_t           ! stability parameter at height zu
      REAL(wp) :: zCD, zCE, zCH, zt_0, zq_0, zt_zu, zq_zu, zt_zt, zq_zt
      REAL(wp) :: zsqrtCDN, z1_sqrtCDN, zlog1, zlog2
      LOGICAL  :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !!----------------------------------------------------------------------------------
      Ni = SIZE(pts_i,1)
      Nj = SIZE(pts_i,2)

      l_zt_equal_zu = ( ABS(zu - zt) < 0.01_wp )

      ! Optimization:
      zlog1    = LOG(zt/zu)
      zlog2    = LOG(zu/10._wp)

      z1_vkarmn  = 1._wp / vkarmn

      zsqrtCDN   = SQRT( pCdN )
      z1_sqrtCDN = 1._wp / zsqrtCDN

      DO jj=1, Nj
         DO ji=1, Ni

            !! Theta & q at the air-sea interface:
            zt_0 = pts_i(ji,jj)
            zq_0 = pqs_i(ji,jj)

            !! Scalar wind speed cannot be below 0.2 m/s
            zUbzu = MAX( pU_zu(ji,jj), wspd_thrshld_ice )

            zt_zt = pt_zt(ji,jj)
            zq_zt = pq_zt(ji,jj)

            !! First guess of temperature and humidity at height zu:
            zt_zu = MAX( zt_zt , 100._wp )   ! who knows what's given on masked-continental regions...
            zq_zu = MAX( zq_zt , epsi10  )   !               "

            !! First guess of transfer coefficients:
            zCD = pCdN
            zCH = pChN
            zCE = pCeN

            zsqrtCD  = SQRT(zCD)             ! == u*/zUbzu

            !! ITERATION BLOCK
            DO jit = 1, nb_iter

               !! Air-Ice differences:
               zdt_zu = zt_zu - zt_0   !  zdt_zu = SIGN( MAX(ABS(zdt_zu),1.E-6_wp), zdt_zu )   !RM 2nd part!!!
               zdq_zu = zq_zu - zq_0   !  zdq_zu = SIGN( MAX(ABS(zdq_zu),1.E-9_wp), zdq_zu )

               !! Now we can get the urbulent scales:
               zus = zsqrtCD * zUbzu
               zdum  = 1._wp / MAX( zsqrtCD, epsi10 )   ! == zUbzu/u*
               zts =      zCH  * zdt_zu * zdum
               zqs =      zCE  * zdq_zu * zdum

               !!Inverse of Obukov length (1/L) :
               z1_L = One_on_L(zt_zu, zq_zu, zus, zts, zqs)  ! 1/L == 1/[Obukhov length]
               z1_L = SIGN( MIN(ABS(z1_L),200._wp), z1_L ) ! (prevents FPE from stupid values from masked region later on...)

               !! Stability parameters "zeta" :
               zzeta_u = zu*z1_L
               zzeta_t = zt*z1_L
               zzeta_u = SIGN( MIN(ABS(zzeta_u),50.0_wp), zzeta_u )
               zzeta_t = SIGN( MIN(ABS(zzeta_t),50.0_wp), zzeta_t )

               !! Update C_D:
               zdum = 1._wp + zsqrtCDN*z1_vkarmn*(zlog2 - psi_m_ice(zzeta_u))
               zCD  = MAX( pCdN / ( zdum*zdum ) , Cx_min )
               !zCD  = MIN( MAX( pCdN / ( zdum*zdum ), Cx_min ) , 1.9E-3_wp ) ! capped version

               zsqrtCD  = SQRT(zCD)             ! == u*/zUbzu

               !! Update C_H and C_E
               zdum = ( zlog2 - psi_h_ice(zzeta_u) ) * z1_vkarmn * z1_sqrtCDN
               zdvm = zsqrtCD * z1_sqrtCDN
               zCH  = MAX( pChN*zdvm / ( 1._wp + pChN*zdum ) , Cx_min )
               zCE  = MAX( pCeN*zdvm / ( 1._wp + pCeN*zdum ) , Cx_min )
               !zCH  = MIN( MAX( pChN*zdvm / ( 1._wp + pChN*zdum ) , Cx_min ) , 1.999E-3_wp ) ! capped version
               !zCE  = MIN( MAX( pCeN*zdvm / ( 1._wp + pCeN*zdum ) , Cx_min ) , 1.999E-3_wp ) ! capped version

               !! Re-updating temperature and humidity at zu if zt /= zu :
               zdum  = z1_vkarmn * ( psi_h_ice(zzeta_u) - psi_h_ice(zzeta_t) + zlog1 )
               zt_zu = MERGE( zt_zu  ,             zt_zt - zts*zdum   ,  l_zt_equal_zu )
               zq_zu = MERGE( zq_zu  ,  MAX(0._wp, zq_zt - zqs*zdum)  ,  l_zt_equal_zu )

            END DO !DO jit = 1, nb_iter

            pCd_i(ji,jj) = zCD
            pCh_i(ji,jj) = zCH
            pCe_i(ji,jj) = zCE

            pt_zu_i(ji,jj) = zt_zu
            pq_zu_i(ji,jj) = zq_zu

         END DO
      END DO

   END SUBROUTINE turb_ice_nanuq




   FUNCTION psi_m_ice( pzeta )
      !!----------------------------------------------------------------------------------
      !! ** Purpose: compute the universal profile stability function for momentum
      !!
      !!     Psi:
      !!     Unstable => Paulson 1970, as in Eq.(30) of Jordan et al. 1999
      !!     Stable   => Grachev et al. 2007, Eq.(12) (integrated form of their phi_m of Eq.(9a))
      !!
      !!             pzeta : stability paramenter, z/L where z is altitude
      !!                     measurement and L is M-O length
      !!
      !! ** Author: L. Brodeau, 2020 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      !$acc routine seq
      !!----------------------------------------------------------------------------------
      REAL(wp)             :: psi_m_ice
      REAL(wp), INTENT(in) :: pzeta
      !!----------------------------------------------------------------------------------
      REAL(wp), PARAMETER :: zam  = 5._wp                               ! a_m (just below Eq.(9b), Grachev et al. 2007)
      REAL(wp), PARAMETER :: zbm  = zam/6.5_wp                          ! b_m (just below Eq.(9b))
      REAL(wp), PARAMETER :: zbbm = ((1._wp - zbm)/zbm)**(1._wp/3._wp)  ! B_m (just below Eq.(12))
      REAL(wp), PARAMETER :: zsr3 = SQRT(3._wp)
      REAL(wp) :: zta, zts, zx, zpsi_u, zpsi_s
      !!----------------------------------------------------------------------------------
      zta = pzeta
      !
      ! Unstable stratification:
      zx = ABS(1._wp - 16._wp*zta)**0.25_wp              !  (16 here, not 15!)

      zpsi_u = LOG( 0.5_wp*(1._wp + zx*zx) ) + 2._wp*LOG( 0.5_wp*(1._wp + zx) ) - 2._wp*ATAN( zx ) + 0.5_wp*rpi  ! Eq.(30) Jordan et al. 1999:

      ! Stable stratification, Eq.(12) Grachev et al. 2007 (integrated form of phi_m of their Eq.(9a);
      ! the former code used phi_m itself, i.e. 1 + 6.5*zta*(1+zta)**(1/3)/(1.3+zta), which is wrong for psi).
      ! `zts` is forced to be >= 0 so that this branch never gives NaN when zeta < -1
      ! (MERGE below evaluates both branches, even if only one is used).
      !
      ! Alternative stable formula (not used): Holtslag & De Bruin 1988 as given in Jordan et al. 1999
      ! (this formula gives -psi, so it would need a minus sign in the MERGE below):
      !zpsi_s = 0.7_wp*zta + 0.75_wp*(zta - 14.3_wp)*EXP( -0.35_wp*zta) + 10.7_wp  ! Eq.(33) Jordan et al. 1999
      zts = MAX( zta , 0._wp )
      zx  = (1._wp + zts)**(1._wp/3._wp)
      zpsi_s = -3._wp*zam/zbm*(zx - 1._wp) + zam*zbbm/(2._wp*zbm) * (                                     &
         &        2._wp*LOG( (zx + zbbm)/(1._wp + zbbm) )                                              &
         &      -       LOG( (zx*zx - zx*zbbm + zbbm*zbbm)/(1._wp - zbbm + zbbm*zbbm) )                 &
         &      + 2._wp*zsr3*( ATAN( (2._wp*zx - zbbm)/(zsr3*zbbm) ) - ATAN( (2._wp - zbbm)/(zsr3*zbbm) ) ) )

      !! Combine: unstable when zeta < 0, stable when zeta >= 0 (zpsi_s is already psi, no minus sign needed)
      psi_m_ice = MERGE( zpsi_u , zpsi_s , zta < 0._wp )
      !
   END FUNCTION psi_m_ice


   FUNCTION psi_h_ice( pzeta )
      !!----------------------------------------------------------------------------------
      !! ** Purpose: compute the universal profile stability function for
      !!             temperature and humidity
      !!
      !!     Psi:
      !!     Unstable => Paulson 1970, as in Eq.(31) of Jordan et al. 1999
      !!     Stable   => Grachev et al. 2007, Eq.(13) (integrated form of their phi_h of Eq.(9b))
      !!
      !!             pzeta : stability paramenter, z/L where z is altitude
      !!                     measurement and L is M-O length
      !!
      !! ** Author: L. Brodeau, 2020 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      !$acc routine seq
      !!----------------------------------------------------------------------------------
      REAL(wp)             :: psi_h_ice
      REAL(wp), INTENT(in) :: pzeta
      !!----------------------------------------------------------------------------------
      REAL(wp), PARAMETER :: zah  = 5._wp        ! a_h (just below Eq.(9b), Grachev et al. 2007)
      REAL(wp), PARAMETER :: zbh  = 5._wp        ! b_h (just below Eq.(9b))
      REAL(wp), PARAMETER :: zch  = 3._wp        ! c_h (just below Eq.(9b))
      REAL(wp), PARAMETER :: zbbh = SQRT(5._wp)  ! B_h = sqrt(c_h^2 - 4) (just below Eq.(13))
      REAL(wp) :: zta, zts, zx, zz, zpsi_u, zpsi_s
      !!----------------------------------------------------------------------------------
      zta = pzeta

      ! Unstable stratification:
      zx = ABS(1._wp - 16._wp*zta)**0.25_wp              !  (16 here, not 15!)

      zpsi_u =   2._wp*LOG( 0.5_wp * (1._wp + zx*zx) )  ! Eq.(31) Jordan et al. 1999

      ! Stable stratification, Eq.(13) Grachev et al. 2007 (integrated form of phi_h of their Eq.(9b);
      ! the former code used phi_h itself, i.e. 1 + 5*zta*(1+zta)/(1+3*zta+zta**2), which is wrong for psi).
      ! `zts` is forced to be >= 0 so that this branch stays finite when zeta < 0
      ! (MERGE below evaluates both branches, even if only one is used).
      !
      ! Alternative stable formula (not used): Holtslag & De Bruin 1988 as given in Jordan et al. 1999
      ! (identical to Psi_m; this formula gives -psi, so it would need a minus sign in the MERGE below):
      !zpsi_s = 0.7_wp*zta + 0.75_wp*(zta - 14.3_wp)*EXP( -0.35_wp*zta) + 10.7_wp       ! Eq.(33) Jordan et al. 1999
      zts = MAX( zta , 0._wp )
      zz  = 2._wp*zts + zch
      zpsi_s = - 0.5_wp*zbh*LOG( 1._wp + zch*zts + zts*zts )                   &
         &     + ( -zah/zbbh + 0.5_wp*zbh*zch/zbbh )                           &
         &       * ( LOG( (zz  - zbbh)/(zz  + zbbh) ) - LOG( (zch - zbbh)/(zch + zbbh) ) )

      !! Combine: unstable when zeta < 0, stable when zeta >= 0 (zpsi_s is already psi, no minus sign needed)
      psi_h_ice = MERGE( zpsi_u , zpsi_s , zta < 0._wp )

   END FUNCTION psi_h_ice

   !!======================================================================
END MODULE mod_blk_ice_nanuq
