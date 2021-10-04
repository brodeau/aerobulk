!!! TO DO: consistent psi_m and psi_h needed!!! For now is those of NCAR !!!
!!
!! DAns les psi_x_ecmwf: zeta est corrigé tel qu'il ne peut pas dépasser 5 en conditions très stable!
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
   !!                   ***  MODULE  mod_blk_andreas  ***
   !! Computes:
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (2m) to zu (10m) if needed
   !!   * the effective bulk wind speed at 10m Ubzu
   !!  according to Andreas et al. (2015)
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !!       Andreas, E.L., Mahrt, L. and Vickers, D. (2015),
   !!       An improved bulk air–sea surface flux algorithm,
   !!       including spray‐mediated transfer.
   !!       Q.J.R. Meteorol. Soc., 141: 642-654. doi:10.1002/qj.2424
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !!
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (2m) to zu (10m) if needed
   !!   * the effective bulk wind speed at z=zu: Ubzu
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
   USE mod_const       !: physical and other constants
   USE mod_phymbl      !: misc. physical functions

   IMPLICIT NONE
   PRIVATE

   INTERFACE u_star_andreas
      MODULE PROCEDURE u_star_andreas_vctr, u_star_andreas_sclr
   END INTERFACE u_star_andreas

   
   !! Important (Brodeau fix):
   REAL(wp), PARAMETER :: rRi_max = 0.15_wp   ! Bulk Ri above which the algorithm fucks up!
   !                                          ! (increasing (>0) Ri means that surface layer increasingly stable and/or wind increasingly weak)
   REAL(wp), PARAMETER :: rCs_min = 0.35E-3_wp ! minimum value to tolarate for CE and CH ! Must be larger than "Cx_min" !!!

   !INTEGER, PARAMETER :: iverbose = 1
   INTEGER, PARAMETER :: iverbose = 0

   PUBLIC :: TURB_ANDREAS, u_star_andreas, psi_m_andreas, psi_h_andreas

   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE turb_andreas( zt, zu, sst, t_zt, ssq, q_zt, U_zu, &
      &                     Cd, Ch, Ce, t_zu, q_zu, Ubzu,       &
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
      !!
      !! OUTPUT :
      !! --------
      !!    *  Cd     : drag coefficient
      !!    *  Ch     : sensible heat coefficient
      !!    *  Ce     : evaporation coefficient
      !!    *  t_zu   : pot. air temperature adjusted at wind height zu       [K]
      !!    *  q_zu   : specific humidity of air        //                    [kg/kg]
      !!    *  Ubzu   : bulk wind speed at zu                                 [m/s]
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
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   sst      ! sea surface temperature                [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   t_zt     ! potential air temperature              [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   ssq      ! sea surface specific humidity           [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   q_zt     ! specific air humidity at zt             [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   U_zu     ! relative wind module at zu                [m/s]
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   Cd       ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   Ch       ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   Ce       ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   t_zu     ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   q_zu     ! spec. humidity adjusted at zu           [kg/kg]
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   Ubzu    ! bulk wind speed at zu                     [m/s]
      !!
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   CdN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   ChN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   CeN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   xz0     ! Aerodynamic roughness length   [m]
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   xu_star ! u*, friction velocity
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   xL      ! zeta (zu/L)
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(:,:) ::   xUN10   ! Neutral wind at zu
      !!
      INTEGER :: Ni, Nj, jit
      LOGICAL :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !!
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   u_star, t_star, q_star
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   z0       ! roughness length (momentum) [m]
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   UN10, UN10_old     ! Neutral wind speed at zu [m/s]
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   zeta_u        ! stability parameter at height zu
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   ztmp0, ztmp1, ztmp2
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   RiB       ! square root of Cd
      !!
      LOGICAL ::  lreturn_cdn=.FALSE., lreturn_chn=.FALSE., lreturn_cen=.FALSE., &
         &        lreturn_z0=.FALSE., lreturn_ustar=.FALSE., lreturn_L=.FALSE., lreturn_UN10=.FALSE.
      CHARACTER(len=40), PARAMETER :: crtnm = 'turb_andreas@mod_blk_andreas.f90'
      !!----------------------------------------------------------------------------------
      Ni = SIZE(sst,1)
      Nj = SIZE(sst,2)

      ALLOCATE( u_star(Ni,Nj), t_star(Ni,Nj), q_star(Ni,Nj), z0(Ni,Nj), &
         &        UN10(Ni,Nj), UN10_old(Ni,Nj), zeta_u(Ni,Nj), RiB(Ni,Nj), &
         &       ztmp0(Ni,Nj),  ztmp1(Ni,Nj),  ztmp2(Ni,Nj) )

      lreturn_cdn   = PRESENT(CdN)
      lreturn_chn   = PRESENT(ChN)
      lreturn_cen   = PRESENT(CeN)
      lreturn_z0    = PRESENT(xz0)
      lreturn_ustar = PRESENT(xu_star)
      lreturn_L     = PRESENT(xL)
      lreturn_UN10  = PRESENT(xUN10)

      l_zt_equal_zu = ( ABS(zu - zt) < 0.01_wp ) ! testing "zu == zt" is risky with double precision

      Ubzu = MAX( 0.25_wp , U_zu )   !  relative wind speed at zu (normally 10m), we don't want to fall under 0.5 m/s

      !! First guess:
      UN10 = Ubzu
      Cd   = 1.1E-3_wp
      Ch   = 1.1E-3_wp
      Ce   = 1.1E-3_wp
      t_zu = t_zt
      q_zu = q_zt

      !! First guess of turbulent scales for scalars:
      ztmp0  = SQRT(Cd)
      t_star = Ch/ztmp0*(t_zu - sst) ! theta*
      q_star = Ce/ztmp0*(q_zu - ssq) ! q*

      ! Bulk Richardson number:
      RiB(:,:) = Ri_bulk( zu, sst, t_zu, ssq, q_zu, Ubzu )


      !! ITERATION BLOCK
      DO jit = 1, nb_iter

         IF(iverbose==1) PRINT *, 'LOLO'

         IF(iverbose==1) PRINT *, 'LOLO *** RiB =', RiB, jit

         WHERE ( RiB < rRi_max )
            !! Normal condition case:
            u_star = u_star_andreas( UN10 )
         ELSEWHERE
            !! Extremely stable + weak wind !!!
            !!  => for we force u* to be consistent with minimum value for CD:
            !!  (otherwize algorithm becomes nonsense...)
            u_star = SQRT(Cx_min) * Ubzu     ! Cd does not go below Cx_min !
         ENDWHERE

         IF(iverbose==1) PRINT *, 'LOLO *** u* =', u_star, jit
         IF(iverbose==2) PRINT *, 'LOLO *** t_zu =', t_zu, jit
         IF(iverbose==2) PRINT *, 'LOLO *** q_zu =', q_zu, jit
         IF(iverbose==1) PRINT *, 'LOLO *** theta* =', t_star, jit
         IF(iverbose==1) PRINT *, 'LOLO *** q* =', q_star, jit

         !! Stability parameter :
         zeta_u = zu*One_on_L( t_zu, q_zu, u_star, t_star, q_star )   ! zu * 1/L

         IF(iverbose==1) PRINT *, 'LOLO *** L =', zu/zeta_u, jit
         IF(iverbose==1) PRINT *, 'LOLO *** zeta_u =', zeta_u, jit
         IF(iverbose==1) PRINT *, 'LOLO *** Ubzu =', Ubzu, jit

         !! Drag coefficient:
         ztmp0 = u_star/Ubzu

         Cd = MAX( ztmp0*ztmp0 , Cx_min )

         IF(iverbose==1) PRINT *, 'LOLO *** CD=', Cd, jit

         !! Roughness length:
         z0 = MIN( z0_from_Cd( zu, Cd,  ppsi=psi_m_andreas(zeta_u) ) , z0_sea_max )
         IF(iverbose==1) PRINT *, 'LOLO *** z0 =', z0, jit
         IF(iverbose==1) PRINT *, 'LOLO'

         !! z0t and z0q, based on LKB, just like into COARE 2.5:
         ztmp0 = z0 * u_star / visc_air(t_zu) ! Re_r
         ztmp1 = z0tq_LKB( 1, ztmp0, z0 )     ! z0t
         ztmp2 = z0tq_LKB( 2, ztmp0, z0 )     ! z0q

         !! Turbulent scales at zu :
         ztmp0 = psi_h_andreas(zeta_u)  ! lolo: zeta_u for scalars???
         IF(iverbose==1) PRINT *, 'LOLO *** psi_h(zeta_u) =', ztmp0, jit
         t_star  = (t_zu - sst)*vkarmn/(LOG(zu) - LOG(ztmp1) - ztmp0)  ! theta* (ztmp1 == z0t in rhs term)
         q_star  = (q_zu - ssq)*vkarmn/(LOG(zu) - LOG(ztmp2) - ztmp0)  !   q*   (ztmp2 == z0q in rhs term)

         IF( (.NOT. l_zt_equal_zu).AND.( jit > 1 ) ) THEN
            !! Re-updating temperature and humidity at zu if zt /= zu:
            ztmp0 = zeta_u/zu*zt   ! zeta_t
            ztmp0 = LOG(zt/zu) + psi_h_andreas(zeta_u) - psi_h_andreas(ztmp0)
            t_zu = t_zt - t_star/vkarmn*ztmp0
            q_zu = q_zt - q_star/vkarmn*ztmp0
            RiB  = Ri_bulk( zu, sst, t_zu, ssq, q_zu, Ubzu ) !LOLO
         ENDIF

         !! Update neutral-stability wind at zu:
         UN10 = MAX( 0.1_wp , UN10_from_ustar( zu, Ubzu, u_star, psi_m_andreas(zeta_u) ) ) ! UN10

         IF(iverbose==1) PRINT *, 'LOLO *** UN10 =', UN10, jit
         IF(iverbose==1) PRINT *, 'LOLO'

      END DO !DO jit = 1, nb_iter

      ! Compute transfer coefficients at zu:
      ztmp0 = u_star/Ubzu

      Cd = MAX( ztmp0*ztmp0 , Cx_min )   ! the earlier use of Cx_min on u* should make use of Cx_min here unnecessary!

      ztmp1 = t_zu - sst ;  ztmp1 = SIGN( MAX(ABS(ztmp1),1.E-6_wp), ztmp1 )  ! dt_zu
      ztmp2 = q_zu - ssq ;  ztmp2 = SIGN( MAX(ABS(ztmp2),1.E-9_wp), ztmp2 )  ! dq_zu
      Ch   = MAX( ztmp0*t_star/ztmp1 , rCs_min )
      Ce   = MAX( ztmp0*q_star/ztmp2 , rCs_min )

      IF( lreturn_cdn .OR. lreturn_chn .OR. lreturn_cen ) ztmp0 = 1._wp/LOG(zu/z0)
      IF( lreturn_cdn )   CdN     = MAX( vkarmn2*ztmp0*ztmp0 , Cx_min )

      IF( lreturn_chn .OR. lreturn_cen ) ztmp1 = z0 * u_star / visc_air(t_zu)  ! Re_r
      IF( lreturn_chn )   ChN     = vkarmn2*ztmp0/LOG(zu/z0tq_LKB( 1, ztmp1, z0 ))
      IF( lreturn_cen )   CeN     = vkarmn2*ztmp0/LOG(zu/z0tq_LKB( 2, ztmp1, z0 ))

      !IF( lreturn_z0 )    xz0     = z0_from_Cd( zu, Cd,  ppsi=psi_m_andreas(zeta_u) )
      IF( lreturn_z0 )    xz0     = z0
      IF( lreturn_ustar ) xu_star = u_star
      IF( lreturn_L )     xL      = zu/zeta_u
      IF( lreturn_UN10 )  xUN10   =  UN10_from_ustar( zu, Ubzu, u_star, psi_m_andreas(zeta_u) )


      DEALLOCATE( u_star, t_star, q_star, z0, UN10, zeta_u, RiB, ztmp0, ztmp1, ztmp2 ) !

   END SUBROUTINE turb_andreas

   
   FUNCTION u_star_andreas_sclr( pun10 )
      !!----------------------------------------------------------------------------------
      !! Estimate of the friction velocity as a function of the neutral-stability wind
      !! speed at at 10m
      !!
      !! Origin: Eq.(2.2) of Andreas et al. (2015)
      !!
      !! ** Author: L. Brodeau, April 2020 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(in) :: pun10          !: neutral-stability scalar wind speed at 10m (m/s)
      REAL(wp)             :: u_star_andreas_sclr !: friction velocity    [m/s]
      !
      REAL(wp) :: za, zt ! local scalars
      !!----------------------------------------------------------------------------------
      za = pun10 - 8.271_wp
      zt = za + SQRT( 0.12_wp*za*za + 0.181_wp )
      u_star_andreas_sclr = 0.239_wp + 0.0433_wp * zt
      !!
   END FUNCTION u_star_andreas_sclr
   !!
   FUNCTION u_star_andreas_vctr( pun10 )
      REAL(wp), DIMENSION(:,:), INTENT(in)             :: pun10               !: neutral-stability scalar wind speed at 10m (m/s)
      REAL(wp), DIMENSION(SIZE(pun10,1),SIZE(pun10,2)) :: u_star_andreas_vctr !: friction velocity    [m/s]
      INTEGER  ::     ji, jj ! dummy loop indices
      DO jj = 1, SIZE(pun10,2)
         DO ji = 1, SIZE(pun10,1)
            u_star_andreas_vctr(ji,jj) = u_star_andreas_sclr( pun10(ji,jj) )
         END DO
      END DO
   END FUNCTION u_star_andreas_vctr


   FUNCTION psi_m_andreas( pzeta )
      !!----------------------------------------------------------------------------------
      !!      Universal profile stability function for momentum
      !!  TO DO !!!!!!!!!!!!!!!!!!!!!
      !! LOLO: paper says Paulson 1970 when unstable and Grachev et al 2007 for STABLE
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! ** Author: L. Brodeau, April 2020 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in) :: pzeta
      REAL(wp), DIMENSION(SIZE(pzeta,1),SIZE(pzeta,2)) :: psi_m_andreas
      !
      REAL(wp), PARAMETER :: zam  = 5._wp      ! a_m (just below Eq.(9b)
      REAL(wp), PARAMETER :: zbm = zam/6.5_wp  ! b_m (just below Eq.(9b)
      !
      REAL(wp), PARAMETER :: z1o3 = 1._wp/3._wp
      REAL(wp), PARAMETER :: zsr3 = SQRT(3._wp)
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) :: zta, zx2, zx, zpsi_unst, zbbm, zpsi_stab,  zstab   ! local scalars
      !!----------------------------------------------------------------------------------
      DO jj = 1, SIZE(pzeta,2)
         DO ji = 1, SIZE(pzeta,1)
            !
            zta = MIN( pzeta(ji,jj) , 15._wp ) !! Very stable conditions (L positif and big!)
            !
            !! *** Unstable: Paulson (1970): #LOLO: DOUBLE CHECK IT IS PAULSON!!!!!
            zx2 = SQRT( ABS(1._wp - 16._wp*zta) )  ! (1 - 16z)^0.5
            zx2 = MAX( zx2 , 1._wp )
            zx  = SQRT(zx2)                          ! (1 - 16z)^0.25
            zpsi_unst = 2._wp*LOG(ABS( (1._wp + zx )*0.5_wp ))   &
               &            + LOG(ABS( (1._wp + zx2)*0.5_wp ))   &
               &          - 2._wp*ATAN(zx) + rpi*0.5_wp
            !
            !! *** Stable: Grachev et al 2007 (SHEBA) [Eq.(12) Grachev et al 2007]:
            zx   = ABS(1._wp + zta)**z1o3
            zbbm = ABS( (1._wp - zbm)/zbm )**z1o3 ! B_m
            !
            zpsi_stab = -3.*zam/zbm*(zx - 1._wp) + zam*zbbm/(2.*zbm) * ( &
               &        2.*LOG(ABS( (   zx     +   zbbm         )/(1._wp        +   zbbm   ) )) &
               &         - LOG(ABS( (zx*zx - zx*zbbm + zbbm*zbbm)/(1._wp - zbbm + zbbm*zbbm) )) &
               & + 2.*zsr3*( ATAN( (2.*zx - zbbm)/(zsr3*zbbm) ) - ATAN( (2._wp - zbbm)/(zsr3*zbbm) ) ) )
            !
            !
            zstab = 0.5_wp + SIGN(0.5_wp, zta) ! zta > 0 => zstab = 1
            !
            psi_m_andreas(ji,jj) =       zstab  * zpsi_stab &  ! (zta > 0) Stable
               &              + (1._wp - zstab) * zpsi_unst    ! (zta < 0) Unstable
            !
         END DO
      END DO
   END FUNCTION psi_m_andreas


   FUNCTION psi_h_andreas( pzeta )
      !!----------------------------------------------------------------------------------
      !! Universal profile stability function for temperature and humidity
      !!
      !!    TO DO
      !!       !! LOLO: paper says Paulson 1970 when unstable and Grachev et al 2007 for STABLE
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in) :: pzeta
      REAL(wp), DIMENSION(SIZE(pzeta,1),SIZE(pzeta,2)) :: psi_h_andreas
      !
      REAL(wp), PARAMETER ::  zah = 5._wp       ! a_h (just below Eq.(9b)
      REAL(wp), PARAMETER ::  zbh = 5._wp       ! b_h (just below Eq.(9b)
      REAL(wp), PARAMETER ::  zch = 3._wp       ! c_h (just below Eq.(9b)
      REAL(wp), PARAMETER :: zbbh = SQRT(5._wp) ! B_h (just below Eq.(13)
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) :: zta, zz, zx2, zpsi_unst, zpsi_stab, zstab  ! local scalars
      !!----------------------------------------------------------------------------------
      DO jj = 1, SIZE(pzeta,2)
         DO ji = 1, SIZE(pzeta,1)
            !
            zta = MIN( pzeta(ji,jj) , 15._wp ) !! Very stable conditions (L positif and large!)
            !
            !! *** Unstable: Paulson (1970): #LOLO: DOUBLE CHECK IT IS PAULSON!!!!!
            zx2 = SQRT( ABS(1._wp - 16._wp*zta) )  ! (1 -16z)^0.5
            zx2 = MAX( zx2 , 1._wp )
            zpsi_unst = 2._wp*LOG( 0.5_wp*(1._wp + zx2) )
            !
            !! *** Stable: Grachev et al 2007 (SHEBA) [Eq.(13) Grachev et al 2007]:
            zz = 2.*zta + zch
            zpsi_stab = - 0.5*zbh*LOG(ABS(1._wp + zch*zta + zta*zta)) &
               &        +  (-zah/zbbh + 0.5*zbh*zch/zbbh)  &
               &          *( LOG(ABS((zz  - zbbh)/(zz  + zbbh))) &
               &           - LOG(ABS((zch - zbbh)/(zch + zbbh)))    )
            !
            zstab = 0.5_wp + SIGN(0.5_wp, zta) ! zta > 0 => zstab = 1
            !
            psi_h_andreas(ji,jj) =            zstab  * zpsi_stab &  ! (zta > 0) Stable
               &                   + (1._wp - zstab) * zpsi_unst    ! (zta < 0) Unstable
            !
         END DO
      END DO
   END FUNCTION psi_h_andreas

   !!======================================================================
END MODULE mod_blk_andreas
