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
   
   IMPLICIT NONE
   PRIVATE

   !! Important (Brodeau fix):
   REAL(wp), PARAMETER :: rCd_min = 0.2E-3_wp ! minimum value to tolarate for CD !
   REAL(wp), PARAMETER :: rRi_max = 0.15_wp   ! Bulk Ri above which the algorithm fucks up! ... and CD (hence u*) forced to stick to rCd_min (sqrt(rCd_min)*U)
   !                                          ! (increasing (>0) Ri means that surface layer increasingly stable and/or wind increasingly weak)
   REAL(wp), PARAMETER :: rCs_min = 0.35E-3_wp ! minimum value to tolarate for CE and CH !
   
   !INTEGER, PARAMETER :: iverbose = 1
   INTEGER, PARAMETER :: iverbose = 0
   
   PUBLIC :: TURB_ANDREAS, psi_m_andreas, psi_h_andreas
   
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
      INTEGER :: j_itt, ja
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
      
      ALLOCATE( u_star(jpi,jpj), t_star(jpi,jpj), q_star(jpi,jpj), z0(jpi,jpj), &
         &        UN10(jpi,jpj), UN10_old(jpi,jpj), zeta_u(jpi,jpj), RiB(jpi,jpj), &
         &       ztmp0(jpi,jpj),  ztmp1(jpi,jpj),  ztmp2(jpi,jpj) )
      
      lreturn_cdn   = PRESENT(CdN)   
      lreturn_chn   = PRESENT(ChN)   
      lreturn_cen   = PRESENT(CeN)   
      lreturn_z0    = PRESENT(xz0)   
      lreturn_ustar = PRESENT(xu_star)
      lreturn_L     = PRESENT(xL)    
      lreturn_UN10  = PRESENT(xUN10) 
      
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

      ! Bulk Richardson number:
      RiB(:,:) = Ri_bulk( zu, sst, t_zu, ssq, q_zu, Ub )
      

      !! ITERATION BLOCK
      DO j_itt = 1, nb_itt
      !DO j_itt = 1, 4

         IF(iverbose==1) PRINT *, 'LOLO'

         IF(iverbose==1) PRINT *, 'LOLO *** RiB =', RiB, j_itt
         
         WHERE ( RiB < rRi_max )
            !! Normal condition case:
            u_star = U_STAR_ANDREAS( UN10 )
         ELSEWHERE
            !! Extremely stable + weak wind !!!
            !!  => for we force u* to be consistent with minimum value for CD:
            !!  (otherwize algorithm becomes nonsense...)
            u_star = SQRT(rCd_min) * Ub     ! Cd does not go below rCd_min !
         ENDWHERE
         
         !u_star = MAX( U_STAR_ANDREAS(UN10) , SQRT(rCd_min)*Ub )


         
         IF(iverbose==1) PRINT *, 'LOLO *** u* =', u_star, j_itt
         IF(iverbose==2) PRINT *, 'LOLO *** t_zu =', t_zu, j_itt
         IF(iverbose==2) PRINT *, 'LOLO *** q_zu =', q_zu, j_itt
         IF(iverbose==1) PRINT *, 'LOLO *** theta* =', t_star, j_itt
         IF(iverbose==1) PRINT *, 'LOLO *** q* =', q_star, j_itt

         !! Stability parameter :
         zeta_u = zu*One_on_L( t_zu, q_zu, u_star, t_star, q_star )   ! zu * 1/L
         
         IF(iverbose==1) PRINT *, 'LOLO *** L =', zu/zeta_u, j_itt
         IF(iverbose==1) PRINT *, 'LOLO *** zeta_u =', zeta_u, j_itt
         IF(iverbose==1) PRINT *, 'LOLO *** Ub =', Ub, j_itt         

         !! Drag coefficient:
         ztmp0 = u_star/Ub
         !Cd    = MAX( ztmp0*ztmp0, rCd_min )
         Cd = ztmp0*ztmp0
         
         IF(iverbose==1) PRINT *, 'LOLO *** CD=', Cd, j_itt

         !! Roughness length:
         z0 = MIN( z0_from_Cd( zu, Cd,  ppsi=psi_m_andreas(zeta_u) ) , z0_sea_max )
         IF(iverbose==1) PRINT *, 'LOLO *** z0 =', z0, j_itt
         IF(iverbose==1) PRINT *, 'LOLO'
         
         !! z0t and z0q, based on LKB, just like into COARE 2.5:
         ztmp0 = z0 * u_star / visc_air(t_zu) ! Re_r
         ztmp1 = z0tq_LKB( 1, ztmp0, z0 )     ! z0t
         ztmp2 = z0tq_LKB( 2, ztmp0, z0 )     ! z0q

         !! Turbulent scales at zu :
         ztmp0 = psi_h_andreas(zeta_u)  ! lolo: zeta_u for scalars???
         IF(iverbose==1) PRINT *, 'LOLO *** psi_h(zeta_u) =', ztmp0, j_itt
         t_star  = (t_zu - sst)*vkarmn/(LOG(zu) - LOG(ztmp1) - ztmp0)  ! theta* (ztmp1 == z0t in rhs term)
         q_star  = (q_zu - ssq)*vkarmn/(LOG(zu) - LOG(ztmp2) - ztmp0)  !   q*   (ztmp2 == z0q in rhs term)
         
         IF( (.NOT. l_zt_equal_zu).AND.( j_itt > 1 ) ) THEN
            !! Re-updating temperature and humidity at zu if zt /= zu:
            ztmp0 = zeta_u/zu*zt   ! zeta_t
            ztmp0 = LOG(zt/zu) + psi_h_andreas(zeta_u) - psi_h_andreas(ztmp0)
            t_zu = t_zt - t_star/vkarmn*ztmp0
            q_zu = q_zt - q_star/vkarmn*ztmp0
            RiB  = Ri_bulk( zu, sst, t_zu, ssq, q_zu, Ub ) !LOLO
         ENDIF

         !! Update neutral-stability wind at zu:
         UN10 = MAX( 0.1_wp , UN10_from_ustar( zu, Ub, u_star, psi_m_andreas(zeta_u) ) ) ! UN10

         IF(iverbose==1) PRINT *, 'LOLO *** UN10 =', UN10, j_itt
         IF(iverbose==1) PRINT *, 'LOLO'

      END DO !DO j_itt = 1, nb_itt
      
      ! Compute transfer coefficients at zu:      
      ztmp0 = u_star/Ub
      
      !Cd    = MAX( ztmp0*ztmp0, rCd_min )
      Cd = ztmp0*ztmp0   ! the earlier use of rCd_min on u* should make use of rCd_min here unnecessary!
      
      ztmp1 = t_zu - sst ;  ztmp1 = SIGN( MAX(ABS(ztmp1),1.E-6_wp), ztmp1 )  ! dt_zu
      ztmp2 = q_zu - ssq ;  ztmp2 = SIGN( MAX(ABS(ztmp2),1.E-9_wp), ztmp2 )  ! dq_zu
      Ch   = MAX( ztmp0*t_star/ztmp1 , rCs_min )
      Ce   = MAX( ztmp0*q_star/ztmp2 , rCs_min )
      
      IF( lreturn_cdn .OR. lreturn_chn .OR. lreturn_cen ) ztmp0 = 1._wp/LOG(zu/z0)
      IF( lreturn_cdn )   CdN     = vkarmn2*ztmp0*ztmp0

      IF( lreturn_chn .OR. lreturn_cen ) ztmp1 = z0 * u_star / visc_air(t_zu)  ! Re_r
      IF( lreturn_chn )   ChN     = vkarmn2*ztmp0/LOG(zu/z0tq_LKB( 1, ztmp1, z0 ))
      IF( lreturn_cen )   CeN     = vkarmn2*ztmp0/LOG(zu/z0tq_LKB( 2, ztmp1, z0 ))

      !IF( lreturn_z0 )    xz0     = z0_from_Cd( zu, Cd,  ppsi=psi_m_andreas(zeta_u) )
      IF( lreturn_z0 )    xz0     = z0
      IF( lreturn_ustar ) xu_star = u_star
      IF( lreturn_L )     xL      = zu/zeta_u
      IF( lreturn_UN10 )  xUN10   =  UN10_from_ustar( zu, Ub, u_star, psi_m_andreas(zeta_u) )


      DEALLOCATE( u_star, t_star, q_star, z0, UN10, zeta_u, RiB, ztmp0, ztmp1, ztmp2 ) !

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
      REAL(wp) :: za, zt, zw, zc ! local scalars
      !!----------------------------------------------------------------------------------
      !zc = SQRT(0.0008_wp)
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            zw  = pun10(ji,jj)
            za = zw - 8.271_wp
            zt = za + SQRT( 0.12_wp*za*za + 0.181_wp )
            u_star_andreas(ji,jj) =   0.239_wp + 0.0433_wp * zt
         END DO
      END DO
      !
   END FUNCTION U_STAR_ANDREAS
   

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
      REAL(wp), DIMENSION(jpi,jpj) :: psi_m_andreas
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      REAL(wp), PARAMETER :: zam  = 5._wp      ! a_m (just below Eq.(9b)
      REAL(wp), PARAMETER :: zbm = zam/6.5_wp  ! b_m (just below Eq.(9b)
      !
      REAL(wp), PARAMETER :: z1o3 = 1._wp/3._wp
      REAL(wp), PARAMETER :: zsr3 = SQRT(3._wp)
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) :: zzeta, zx2, zx, zpsi_unst, zbbm, zpsi_stab,  zstab   ! local scalars
      !!----------------------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi

            zzeta = MIN( pzeta(ji,jj) , 15._wp ) !! Very stable conditions (L positif and big!)
            !
            !! *** Unstable: Paulson (1970): #LOLO: DOUBLE CHECK IT IS PAULSON!!!!!
            zx2 = SQRT( ABS(1._wp - 16._wp*zzeta) )  ! (1 - 16z)^0.5
            zx2 = MAX( zx2 , 1._wp )
            zx  = SQRT(zx2)                          ! (1 - 16z)^0.25
            zpsi_unst = 2._wp*LOG(ABS( (1._wp + zx )*0.5_wp ))   &
               &            + LOG(ABS( (1._wp + zx2)*0.5_wp ))   &
               &          - 2._wp*ATAN(zx) + rpi*0.5_wp
            !
            !! *** Stable: Grachev et al 2007 (SHEBA) [Eq.(12) Grachev et al 2007]:
            zx   = ABS(1._wp + zzeta)**z1o3
            zbbm = ABS( (1._wp - zbm)/zbm )**z1o3 ! B_m
            !
            zpsi_stab = -3.*zam/zbm*(zx - 1._wp) + zam*zbbm/(2.*zbm) * ( &
               &        2.*LOG(ABS( (   zx     +   zbbm         )/(1._wp        +   zbbm   ) )) &
               &         - LOG(ABS( (zx*zx - zx*zbbm + zbbm*zbbm)/(1._wp - zbbm + zbbm*zbbm) )) &
               & + 2.*zsr3*( ATAN( (2.*zx - zbbm)/(zsr3*zbbm) ) - ATAN( (2._wp - zbbm)/(zsr3*zbbm) ) ) )
            !
            !
            zstab = 0.5_wp + SIGN(0.5_wp, zzeta) ! zzeta > 0 => zstab = 1
            !
            psi_m_andreas(ji,jj) =       zstab  * zpsi_stab &  ! (zzeta > 0) Stable
               &              + (1._wp - zstab) * zpsi_unst    ! (zzeta < 0) Unstable
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
      REAL(wp), DIMENSION(jpi,jpj) :: psi_h_andreas
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      REAL(wp), PARAMETER ::  zah = 5._wp       ! a_h (just below Eq.(9b)
      REAL(wp), PARAMETER ::  zbh = 5._wp       ! b_h (just below Eq.(9b)
      REAL(wp), PARAMETER ::  zch = 3._wp       ! c_h (just below Eq.(9b)
      REAL(wp), PARAMETER :: zbbh = SQRT(5._wp) ! B_h (just below Eq.(13)
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) :: zzeta, zz, zx2, zpsi_unst, zpsi_stab, zstab  ! local scalars
      !!----------------------------------------------------------------------------------
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zzeta = MIN( pzeta(ji,jj) , 15._wp ) !! Very stable conditions (L positif and large!)
            !
            !! *** Unstable: Paulson (1970): #LOLO: DOUBLE CHECK IT IS PAULSON!!!!!
            zx2 = SQRT( ABS(1._wp - 16._wp*zzeta) )  ! (1 -16z)^0.5
            zx2 = MAX( zx2 , 1._wp )
            zpsi_unst = 2._wp*LOG( 0.5_wp*(1._wp + zx2) )
            !            
            !! *** Stable: Grachev et al 2007 (SHEBA) [Eq.(13) Grachev et al 2007]:
            zz = 2.*zzeta + zch
            zpsi_stab = - 0.5*zbh*LOG(ABS(1._wp + zch*zzeta + zzeta*zzeta)) &
               &        +  (-zah/zbbh + 0.5*zbh*zch/zbbh)  &
               &          *( LOG(ABS((zz  - zbbh)/(zz  + zbbh))) &
               &           - LOG(ABS((zch - zbbh)/(zch + zbbh)))    )
            !
            zstab = 0.5_wp + SIGN(0.5_wp, zzeta) ! zzeta > 0 => zstab = 1
            !
            psi_h_andreas(ji,jj) =            zstab  * zpsi_stab &  ! (zzeta > 0) Stable
               &                   + (1._wp - zstab) * zpsi_unst    ! (zzeta < 0) Unstable
            !
         END DO
      END DO
   END FUNCTION psi_h_andreas

   !!======================================================================
END MODULE mod_blk_andreas
