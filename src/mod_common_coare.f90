! AeroBulk / 2021 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following paper:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.


MODULE mod_common_coare

   USE mod_const,   ONLY: wp, rpi, grav, vkarmn, vkarmn2
   USE  mod_phymbl, ONLY: visc_air, Ri_bulk

   IMPLICIT NONE

   INTERFACE first_guess_coare
      MODULE PROCEDURE first_guess_coare_vctr, first_guess_coare_sclr
   END INTERFACE first_guess_coare

   INTERFACE psi_m_coare
      MODULE PROCEDURE psi_m_coare_vctr, psi_m_coare_sclr
   END INTERFACE psi_m_coare

   INTERFACE psi_h_coare
      MODULE PROCEDURE psi_h_coare_vctr, psi_h_coare_sclr
   END INTERFACE psi_h_coare

CONTAINS


   !===============================================================================================
   SUBROUTINE FIRST_GUESS_COARE_SCLR( zt, zu, psst, t_zt, pssq, q_zt, U_zu, pcharn, &
      &                               pus, pts, pqs, t_zu, q_zu, Ubzu,  pz0 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  FIRST_GUESS_COARE_SCLR  ***
      !!
      !! ** Purpose :  Computes fairly accurate first guess of u*, theta* and q*
      !!               by means of the method developed by Fairall et al in the
      !!               COARE family of algorithms
      !!               Purpose is to limit the number of itteration needed...
      !!
      !! INPUT :
      !! -------
      !!    *  zt   : height for temperature and spec. hum. of air            [m]
      !!    *  zu   : height for wind speed (usually 10m)                     [m]
      !!    *  psst  : bulk SST                                                [K]
      !!    *  t_zt : potential air temperature at zt                         [K]
      !!    *  pssq  : SSQ aka saturation specific humidity at temp. psst       [kg/kg]
      !!    *  q_zt : specific humidity of air at zt                          [kg/kg]
      !!    *  U_zu : scalar wind speed at zu                                 [m/s]
      !!    *  pcharn: Charnock parameter
      !!
      !! OUTPUT :
      !! --------
      !!    *  pus    : FIRST GUESS of u* aka friction velocity                            [m/s]
      !!    *  pts    : FIRST GUESS of theta*                                               [K]
      !!    *  pqs    : FIRST GUESS of q* aka friction velocity                            [kg/kg]
      !!    *  t_zu   : FIRST GUESS of pot. air temperature adjusted at wind height zu      [K]
      !!    *  q_zu   : FIRST GUESS of specific humidity of air        //                   [kg/kg]
      !!    *  Ubzu   : FIRST GUESS of bulk wind speed at zu                                [m/s]
      !!
      !! ** Author: L. Brodeau, May 2021 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(in)  ::   zt
      REAL(wp), INTENT(in)  ::   zu
      REAL(wp), INTENT(in)  ::   psst
      REAL(wp), INTENT(in)  ::   t_zt
      REAL(wp), INTENT(in)  ::   pssq
      REAL(wp), INTENT(in)  ::   q_zt
      REAL(wp), INTENT(in)  ::   U_zu
      REAL(wp), INTENT(in)  ::   pcharn
      !!
      REAL(wp), INTENT(out) ::   pus
      REAL(wp), INTENT(out) ::   pts
      REAL(wp), INTENT(out) ::   pqs
      REAL(wp), INTENT(out) ::   t_zu
      REAL(wp), INTENT(out) ::   q_zu
      REAL(wp), INTENT(out) ::   Ubzu
      REAL(wp), INTENT(out), OPTIONAL :: pz0    ! roughness length [m]
      !
      LOGICAL :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !
      REAL(wp) :: zdt, zdq, zUb, zCd, zus, zts, zqs, zNu_a, zRib
      REAL(wp) :: zlog_zt_o_zu, zlog_10, zlog_zt, zlog_zu, zlog_z0, zlog_z0t, z1_o_sqrt_Cd10, zcc, zcc_ri, z1_o_Ribcu
      REAL(wp) :: zprf, zc_a, zc_b
      REAL(wp) :: zz0, zz0t, zstab, zzeta_u, zzeta_t
      REAL(wp) :: ztmp
      !
      REAL(wp), PARAMETER :: zzi0=600._wp, zBeta0=1.2_wp ! COARE values, after all it's a coare method...
      !
      CHARACTER(len=40), PARAMETER :: crtnm = 'FIRST_GUESS_COARE_SCLR@mod_common_coare.f90'
      !!----------------------------------------------------------------------------------
      l_zt_equal_zu = ( ABS(zu - zt) < 0.01_wp )

      !! First guess of temperature and humidity at height zu:
      t_zu = MAX( t_zt ,  180._wp )   ! who knows what's given on masked-continental regions...
      q_zu = MAX( q_zt , 1.e-6_wp )   !               "

      zz0 = 0.0001 ! "rough" first guess of roughness length of sea surface...

      !! Constants:
      zlog_10 = LOG(10._wp)
      zlog_zt = LOG(zt)
      zlog_zu = LOG(zu)
      zlog_zt_o_zu = LOG(zt/zu)
      zc_a = 0.035_wp*LOG(10._wp/zz0)/LOG(zu/zz0)   !       "                    "               "
      zc_b = 0.004_wp*zzi0*zBeta0*zBeta0*zBeta0

      !! Air-sea differences (and we don't want them to be 0...)
      zdt = t_zu - psst ;   zdt = SIGN( MAX(ABS(zdt),1.E-6_wp), zdt )
      zdq = q_zu - pssq ;   zdq = SIGN( MAX(ABS(zdq),1.E-9_wp), zdq )

      zNu_a = visc_air(t_zu) ! Air viscosity (m^2/s) at zt given from temperature in (K)

      zUb = SQRT(U_zu*U_zu + 0.5_wp*0.5_wp) ! initial guess for wind gustiness contribution

      zus = zc_a*zUb

      !! Update roughness length:
      zz0     = pcharn*zus*zus/grav + 0.11_wp*zNu_a/zus
      zz0     = MIN( MAX(ABS(zz0), 1.E-8) , 1._wp )      ! (prevents FPE from stupid values from masked region later on)
      zlog_z0 = LOG(zz0)

      !! First guess of Cd
      zCd          = (vkarmn/(zlog_zu - zlog_z0))**2
      z1_o_sqrt_Cd10 =       (zlog_10 - zlog_z0)/vkarmn  ! sum less costly than product: log(a/b) == log(a) - log(b)

      !! Temperature/humidity roughness lengthes:
      zz0t  = 10._wp / EXP( vkarmn/( 0.00115*z1_o_sqrt_Cd10 ) )
      zz0t    = MIN( MAX(ABS(zz0t), 1.E-8) , 1._wp )         ! (prevents FPE from stupid values from masked region later on)
      zlog_z0t = LOG(zz0t)

      !! Bulk Richardson Number (BRN)
      zRib = Ri_bulk( zu, psst, t_zu, pssq, q_zu, zUb )

      !! First estimate of zeta_u, depending on the stability, ie sign of BRN (zRib):
      zcc = vkarmn2/(zCd*(zlog_zt - zlog_z0t))
      zcc_ri  = zcc*zRib
      z1_o_Ribcu = -zc_b/zu
      zstab   = 0.5 + SIGN( 0.5_wp , zRib )
      zzeta_u = (1._wp - zstab) *   zcc_ri / (1._wp + zRib*z1_o_Ribcu) & !  Ri_bulk < 0, Unstable
         &  +        zstab      * ( zcc_ri + 27._wp/9._wp*zRib*zRib )    !  Ri_bulk > 0, Stable


      !! u*, theta*, q* :
      zus  = MAX ( zUb*vkarmn/(zlog_zu - zlog_z0  - psi_m_coare(zzeta_u)) , 1.E-9 ) ! (MAX => prevents FPE from stupid values from masked region later on)
      ztmp = vkarmn/(zlog_zu - zlog_z0t - psi_h_coare(zzeta_u))
      zts  = zdt*ztmp
      zqs  = zdq*ztmp
      
      !! Adjustment of theta and q from zt to zu if relevant:
      IF( .NOT. l_zt_equal_zu ) THEN
         zzeta_t = zt*zzeta_u/zu
         zprf = LOG(zt/zu) + psi_h_coare(zzeta_u) - psi_h_coare(zzeta_t)
         t_zu = t_zt - zts/vkarmn*zprf
         q_zu = q_zt - zqs/vkarmn*zprf
         q_zu = (0.5_wp + SIGN(0.5_wp,q_zu))*q_zu ! prevents negative humidity...
         !!
         !! Update of theta and q air-sea differences and theta*, q* :
         zdt = t_zu - psst  ; zdt = SIGN( MAX(ABS(zdt),1.E-6_wp), zdt )
         zdq = q_zu - pssq  ; zdq = SIGN( MAX(ABS(zdq),1.E-9_wp), zdq )
         zts = zdt*ztmp
         zqs = zdq*ztmp
      ENDIF
      
      !! Output result:
      pus  = zus
      pts  = zts
      pqs  = zqs
      Ubzu = zUb

      IF( PRESENT(pz0) ) THEN
         !! Again, because new zus:
         zz0 = pcharn*zus*zus/grav + 0.11_wp*zNu_a/zus
         pz0 = MIN( MAX(ABS(zz0), 1.E-8) , 1._wp )    ! (prevents FPE from stupid values from masked region later on)
      END IF

   END SUBROUTINE FIRST_GUESS_COARE_SCLR

   SUBROUTINE FIRST_GUESS_COARE_VCTR( zt, zu, psst, t_zt, pssq, q_zt, U_zu, pcharn, &
      &                               pus, pts, pqs, t_zu, q_zu, Ubzu,  qz0 )
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(in)                  ::   zt
      REAL(wp), INTENT(in)                  ::   zu
      REAL(wp), INTENT(in),  DIMENSION(:,:) ::   psst
      REAL(wp), INTENT(in),  DIMENSION(:,:) ::   t_zt
      REAL(wp), INTENT(in),  DIMENSION(:,:) ::   pssq
      REAL(wp), INTENT(in),  DIMENSION(:,:) ::   q_zt
      REAL(wp), INTENT(in),  DIMENSION(:,:) ::   U_zu
      REAL(wp), INTENT(in),  DIMENSION(:,:) ::   pcharn
      !!
      REAL(wp), INTENT(out), DIMENSION(:,:) ::   pus
      REAL(wp), INTENT(out), DIMENSION(:,:) ::   pts
      REAL(wp), INTENT(out), DIMENSION(:,:) ::   pqs
      REAL(wp), INTENT(out), DIMENSION(:,:) ::   t_zu
      REAL(wp), INTENT(out), DIMENSION(:,:) ::   q_zu
      REAL(wp), INTENT(out), DIMENSION(:,:) ::   Ubzu
      REAL(wp), INTENT(out), DIMENSION(:,:), OPTIONAL :: qz0    ! roughness length [m]
      !
      INTEGER  :: ji, jj
      REAL(wp) :: zz0
      !!----------------------------------------------------------------------------------
      DO jj = 1, SIZE(psst,2)
         DO ji = 1, SIZE(psst,1)
            CALL FIRST_GUESS_COARE_SCLR( zt, zu, psst(ji,jj), t_zt(ji,jj), pssq(ji,jj), q_zt(ji,jj), U_zu(ji,jj), pcharn(ji,jj), &
               &                         pus(ji,jj), pts(ji,jj), pqs(ji,jj), t_zu(ji,jj), q_zu(ji,jj), Ubzu(ji,jj),  pz0=zz0 )
            IF( PRESENT(qz0) ) qz0(ji,jj) = zz0
         END DO
      END DO
   END SUBROUTINE FIRST_GUESS_COARE_VCTR
   !!==============================================================================================



   !!==============================================================================================
   FUNCTION psi_m_coare_sclr( pzeta )
      !!----------------------------------------------------------------------------------
      !! ** Purpose: compute the universal profile stability function for momentum
      !!             COARE 3.0, Fairall et al. 2003
      !!             pzeta : stability paramenter, z/L where z is altitude
      !!                     measurement and L is M-O length
      !!       Stability function for wind speed and scalars matching Kansas and free
      !!       convection forms with weighting f convective form, follows Fairall et
      !!       al (1996) with profile constants from Grachev et al (2000) BLM stable
      !!       form from Beljaars and Holtslag (1991)
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp) :: psi_m_coare_sclr
      REAL(wp), INTENT(in) :: pzeta
      !!
      REAL(wp) :: zphi_m, zphi_c, zpsi_k, zpsi_c, zf, zc, zstb
      !!----------------------------------------------------------------------------------
      zphi_m = ABS(1. - 15.*pzeta)**.25    !!Kansas unstable
      !
      zpsi_k = 2.*LOG((1. + zphi_m)/2.) + LOG((1. + zphi_m*zphi_m)/2.)   &
         & - 2.*ATAN(zphi_m) + 0.5*rpi
      !
      zphi_c = ABS(1. - 10.15*pzeta)**.3333                   !!Convective
      !
      zpsi_c = 1.5*LOG((1. + zphi_c + zphi_c*zphi_c)/3.) &
         &     - 1.7320508*ATAN((1. + 2.*zphi_c)/1.7320508) + 1.813799447
      !
      zf = pzeta*pzeta
      zf = zf/(1. + zf)
      zc = MIN(50._wp, 0.35_wp*pzeta)
      zstb = 0.5 + SIGN(0.5_wp, pzeta)
      !
      psi_m_coare_sclr = (1. - zstb) * ( (1. - zf)*zpsi_k + zf*zpsi_c ) & ! (pzeta < 0)
         &           -   zstb  * ( 1. + 1.*pzeta     &                ! (pzeta > 0)
         &                          + 0.6667*(pzeta - 14.28)/EXP(zc) + 8.525 )  !     "
      !!
   END FUNCTION psi_m_coare_sclr

   FUNCTION psi_m_coare_vctr( pzeta )
      !!----------------------------------------------------------------------------------
      !! ** Purpose: compute the universal profile stability function for momentum
      !!             COARE 3.0, Fairall et al. 2003
      !!             pzeta : stability paramenter, z/L where z is altitude
      !!                     measurement and L is M-O length
      !!       Stability function for wind speed and scalars matching Kansas and free
      !!       convection forms with weighting f convective form, follows Fairall et
      !!       al (1996) with profile constants from Grachev et al (2000) BLM stable
      !!       form from Beljaars and Holtslag (1991)
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in)             :: pzeta
      REAL(wp), DIMENSION(SIZE(pzeta,1),SIZE(pzeta,2)) :: psi_m_coare_vctr
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) :: zta, zphi_m, zphi_c, zpsi_k, zpsi_c, zf, zc, zstab
      !!----------------------------------------------------------------------------------
      DO jj = 1, SIZE(pzeta,2)
         DO ji = 1, SIZE(pzeta,1)
            !
            zta = pzeta(ji,jj)
            !
            zphi_m = ABS(1. - 15.*zta)**.25    !!Kansas unstable
            !
            zpsi_k = 2.*LOG((1. + zphi_m)/2.) + LOG((1. + zphi_m*zphi_m)/2.)   &
               & - 2.*ATAN(zphi_m) + 0.5*rpi
            !
            zphi_c = ABS(1. - 10.15*zta)**.3333                   !!Convective
            !
            zpsi_c = 1.5*LOG((1. + zphi_c + zphi_c*zphi_c)/3.) &
               &     - 1.7320508*ATAN((1. + 2.*zphi_c)/1.7320508) + 1.813799447
            !
            zf = zta*zta
            zf = zf/(1. + zf)
            zc = MIN(50._wp, 0.35_wp*zta)
            zstab = 0.5 + SIGN(0.5_wp, zta)
            !
            psi_m_coare_vctr(ji,jj) = (1. - zstab) * ( (1. - zf)*zpsi_k + zf*zpsi_c ) & ! (zta < 0)
               &                -   zstab     * ( 1. + 1.*zta     &                ! (zta > 0)
               &                         + 0.6667*(zta - 14.28)/EXP(zc) + 8.525 )  !     "
         END DO
      END DO
   END FUNCTION psi_m_coare_vctr
   !!==============================================================================================


   !!==============================================================================================
   FUNCTION psi_h_coare_sclr( pzeta )
      !!---------------------------------------------------------------------
      !! Universal profile stability function for temperature and humidity
      !! COARE 3.0, Fairall et al. 2003
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! Stability function for wind speed and scalars matching Kansas and free
      !! convection forms with weighting f convective form, follows Fairall et
      !! al (1996) with profile constants from Grachev et al (2000) BLM stable
      !! form from Beljaars and Holtslag (1991)
      !!
      !! Author: L. Brodeau, June 2016 / AeroBulk
      !!         (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------
      REAL(wp) :: psi_h_coare_sclr
      REAL(wp), INTENT(in) :: pzeta
      !!
      REAL(wp) :: zphi_h, zphi_c, zpsi_k, zpsi_c, zf, zc, zstb
      !!----------------------------------------------------------------
      zphi_h = (ABS(1. - 15.*pzeta))**.5  !! Kansas unstable   (zphi_h = zphi_m**2 when unstable, zphi_m when stable)
      !
      zpsi_k = 2.*LOG((1. + zphi_h)/2.)
      !
      zphi_c = (ABS(1. - 34.15*pzeta))**.3333   !! Convective
      !
      zpsi_c = 1.5*LOG((1. + zphi_c + zphi_c*zphi_c)/3.) &
         &    -1.7320508*ATAN((1. + 2.*zphi_c)/1.7320508) + 1.813799447
      !
      zf = pzeta*pzeta
      zf = zf/(1. + zf)
      zc = MIN(50._wp,0.35_wp*pzeta)
      zstb = 0.5 + SIGN(0.5_wp, pzeta)
      !
      psi_h_coare_sclr = (1.-zstb) * ( (1. - zf)*zpsi_k + zf*zpsi_c ) &
         &                  -zstb  * ( (ABS(1. + 2.*pzeta/3.))**1.5     &
         &                            + .6667*(pzeta - 14.28)/EXP(zc) + 8.525 )
      !!
   END FUNCTION psi_h_coare_sclr

   FUNCTION psi_h_coare_vctr( pzeta )
      !!---------------------------------------------------------------------
      !! Universal profile stability function for temperature and humidity
      !! COARE 3.0, Fairall et al. 2003
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! Stability function for wind speed and scalars matching Kansas and free
      !! convection forms with weighting f convective form, follows Fairall et
      !! al (1996) with profile constants from Grachev et al (2000) BLM stable
      !! form from Beljaars and Holtslag (1991)
      !!
      !! Author: L. Brodeau, June 2016 / AeroBulk
      !!         (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in)             :: pzeta
      REAL(wp), DIMENSION(SIZE(pzeta,1),SIZE(pzeta,2)) :: psi_h_coare_vctr
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) :: zta, zphi_h, zphi_c, zpsi_k, zpsi_c, zf, zc, zstab
      !!----------------------------------------------------------------
      DO jj = 1, SIZE(pzeta,2)
         DO ji = 1, SIZE(pzeta,1)
            !
            zta = pzeta(ji,jj)
            !
            zphi_h = (ABS(1. - 15.*zta))**.5  !! Kansas unstable   (zphi_h = zphi_m**2 when unstable, zphi_m when stable)
            !
            zpsi_k = 2.*LOG((1. + zphi_h)/2.)
            !
            zphi_c = (ABS(1. - 34.15*zta))**.3333   !! Convective
            !
            zpsi_c = 1.5*LOG((1. + zphi_c + zphi_c*zphi_c)/3.) &
               &    -1.7320508*ATAN((1. + 2.*zphi_c)/1.7320508) + 1.813799447
            !
            zf = zta*zta
            zf = zf/(1. + zf)
            zc = MIN(50._wp,0.35_wp*zta)
            zstab = 0.5 + SIGN(0.5_wp, zta)
            !
            psi_h_coare_vctr(ji,jj) = (1. - zstab) * ( (1. - zf)*zpsi_k + zf*zpsi_c ) &
               &                -   zstab     * ( (ABS(1. + 2.*zta/3.))**1.5     &
               &                           + .6667*(zta - 14.28)/EXP(zc) + 8.525 )
         END DO
      END DO
   END FUNCTION psi_h_coare_vctr
   !!==============================================================================================

END MODULE mod_common_coare
