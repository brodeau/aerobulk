! AeroBulk / 2020 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_cdn_form_ice
   !!====================================================================================
   !!            Author: Laurent Brodeau, January 2020
   !!====================================================================================
   USE mod_const       !: physical and other constants
   USE mod_phymbl      !: misc. physical functions

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: CdN10_f_LU12, CdN_f_LU12_eq36, CdN10_f_LU13, CdN_f_LG15, CdN_f_LG15_light

   REAL(wp), PARAMETER :: rCe_0    = 2.23E-3_wp !LOLO: this one can be more accurate when sea-ice data => Lupkes et al (2013), Eq.(1)
   REAL(wp), PARAMETER :: rNu_0    = 1._wp
   REAL(wp), PARAMETER :: rMu_0    = 1._wp
   REAL(wp), PARAMETER :: rbeta_0  = 1.4_wp     ! (Eq.47) MIZ

   REAL(wp), PARAMETER :: rhmin_0 = 0.286_wp  ! Eq.(25)
   REAL(wp), PARAMETER :: rhmax_0 = 0.534_wp  ! Eq.(25)
   REAL(wp), PARAMETER :: rDmin_0 =   8._wp      ! Eq.(27)
   REAL(wp), PARAMETER :: rDmax_0 = 300._wp      ! Eq.(27)
   REAL(wp), PARAMETER :: rz0_w_0 = 3.27E-4   ! fixed roughness length over water (paragraph below Eq.36)

   !!============================================================
   REAL(wp), PARAMETER ::   rce10_i_0 = 3.46e-3_wp ! (Eq.48) MIZ

   REAL(wp), PARAMETER ::   ralpha_0  = 0.2_wp     ! (Eq.12) (ECHAM6 value)

   !!----------------------------------------------------------------------
CONTAINS


   FUNCTION CdN10_f_LU12( pfrice, pz0w,  pSc, phf, pDi  )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  CdN10_f_LU12  ***
      !!
      !!        GENERAL FORM OF EQUATION 22 of Lupkes et al. 2012
      !!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !!
      !! ** Purpose :    Computes the "form" contribution of the neutral air-ice
      !!                 drag referenced at 10m to make it dependent on edges at
      !!                 leads, melt ponds and flows (to be added to the "skin"
      !!                 contribution. After some
      !!                 approximations, this can be resumed to a dependency on
      !!                 ice concentration.
      !!
      !! ** References : Lupkes et al. JGR 2012 (theory)
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in)           :: pfrice ! ice concentration [fraction]  => at_i_b  ! NOT USED if pSc, phf and pDi all provided...
      REAL(wp), DIMENSION(:,:), INTENT(in)           :: pz0w   ! roughness length over water  [m]
      REAL(wp), DIMENSION(:,:), INTENT(in), OPTIONAL :: pSc    ! sheltering function [0-1] (Sc->1 for large distance between floes, ->0 for small distances)
      REAL(wp), DIMENSION(:,:), INTENT(in), OPTIONAL :: phf    ! mean freeboard of floes    [m]
      REAL(wp), DIMENSION(:,:), INTENT(in), OPTIONAL :: pDi    ! cross wind dimension of the floe (aka effective edge length for form drag)   [m]
      REAL(wp), DIMENSION(SIZE(pfrice,1),SIZE(pfrice,2))                       :: CdN10_f_LU12  ! neutral FORM drag coefficient contribution over sea-ice
      !!----------------------------------------------------------------------
      LOGICAL :: l_known_Sc=.FALSE., l_known_hf=.FALSE., l_known_Di=.FALSE.
      REAL(wp) :: ztmp, zrlog, zfri, zfrw, zSc, zhf, zDi
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------
      l_known_Sc = PRESENT(pSc)
      l_known_hf = PRESENT(phf)
      l_known_Di = PRESENT(pDi)

      DO jj = 1, SIZE(pfrice,2)
         DO ji = 1, SIZE(pfrice,1)

            zfri = pfrice(ji,jj)
            zfrw = (1._wp - zfri)

            IF(l_known_Sc) THEN
               zSc = pSc(ji,jj)
            ELSE
               !! Sc parameterized in terms of A (ice fraction):
               zSc = zfrw**(1._wp / ( 10._wp * rBeta_0 ))   ! Eq.(31)
               !PRINT *, 'LOLO: Sc PARAMETERIZED !!! =>', zSc
            END IF

            IF(l_known_hf) THEN
               zhf = phf(ji,jj)
            ELSE
               !! hf parameterized in terms of A (ice fraction):
               zhf = rhmax_0*zfri + rhmin_0*zfrw  ! Eq.(25)
               !PRINT *, 'LOLO: hf PARAMETERIZED !!! =>', zhf
            END IF

            IF(l_known_Di) THEN
               zDi = pDi(ji,jj)
            ELSE
               !! Di parameterized in terms of A (ice fraction):
               ztmp = 1._wp / ( 1._wp - (rDmin_0/rDmax_0)**(1._wp/rBeta_0) )   ! A* Eq.(27)
               zDi =  rDmin_0 * ( ztmp/(ztmp - zfri) )**rBeta_0                !    Eq.(26)
               !PRINT *, 'LOLO: Di PARAMETERIZED !!! =>', zDi
            END IF

            ztmp  = 1._wp/pz0w(ji,jj)
            zrlog = LOG(zhf*ztmp) / LOG(10._wp*ztmp)

            !#bug reported by V. Guemas: CdN10_f_LU12(:,:) = 0.5_wp* 0.3_wp * zrlog*zrlog * zSc*zSc * zhf/zDi * zfri  ! Eq.(22)
            CdN10_f_LU12(:,:) = 0.5_wp* 0.3_wp * zrlog*zrlog * zSc * zhf/zDi * zfri  ! Eq.(22)
            !!                   1/2      Ce

         END DO
      END DO
   END FUNCTION CdN10_f_LU12


   FUNCTION CdN_f_LU12_eq36( pzu, pfrice )
      REAL(wp),                 INTENT(in) :: pzu    ! reference height                       [m]
      REAL(wp), DIMENSION(:,:), INTENT(in) :: pfrice ! ice concentration [fraction]  => at_i_b  ! NOT USED if pSc, phf and pDi all provided...
      REAL(wp), DIMENSION(SIZE(pfrice,1),SIZE(pfrice,2))             :: CdN_f_LU12_eq36 ! neutral FORM drag coefficient contribution over sea-ice
      !!----------------------------------------------------------------------
      REAL(wp) :: ztmp, zrlog, zfri, zhf, zDi
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------
      !zhf   = 0.28   ! h_fc
      zhf   = 0.41   ! h_fc
      zDi   = rDmin_0

      ztmp  = 1._wp/rz0_w_0
      zrlog = LOG(zhf*ztmp) / LOG(pzu*ztmp)

      DO jj = 1, SIZE(pfrice,2)
         DO ji = 1, SIZE(pfrice,1)

            zfri = pfrice(ji,jj)

            CdN_f_LU12_eq36(:,:) = 0.5_wp* 0.3_wp * zrlog*zrlog * zhf/zDi  * (1._wp - zfri)**rBeta_0 ! Eq.(35) & (36)
            !!                        1/2      Ce

         END DO
      END DO
   END FUNCTION CdN_f_LU12_eq36




   FUNCTION CdN10_f_LU13( pfrice )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  CdN10_f_LU13  ***
      !!
      !! ** Purpose :    Computes the "form" contribution of the neutral air-ice
      !!                 drag referenced at 10m to make it dependent on edges at
      !!                 leads, melt ponds and flows (to be added to the "skin"
      !!                 contribution. After some
      !!                 approximations, this can be resumed to a dependency on
      !!                 ice concentration.
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
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pfrice           ! ice concentration [fraction]  => at_i_b
      REAL(wp), DIMENSION(SIZE(pfrice,1),SIZE(pfrice,2))              :: CdN10_f_LU13  ! neutral FORM drag coefficient contribution over sea-ice
      !!----------------------------------------------------------------------
      REAL(wp)            ::   zcoef
      !!----------------------------------------------------------------------
      zcoef = rNu_0 + 1._wp / ( 10._wp * rBeta_0 )

      !! We are not an AGCM, we are an OGCM!!! => we drop term "(1 - A)*Cd_w"
      !!  => so we keep only the last rhs terms of Eq.(1) of Lupkes et al, 2013 that we divide by "A":
      !! (we multiply Cd_i_s and Cd_i_f by A later, when applying ocean-ice partitioning...

      CdN10_f_LU13(:,:) = rCe_0 * pfrice(:,:)**(rMu_0 - 1._wp) * (1._wp - pfrice(:,:))**zcoef
      !! => seems okay for winter 100% sea-ice as second rhs term vanishes as pfrice == 1....

   END FUNCTION CdN10_f_LU13


   FUNCTION CdN_f_LG15( pzu, pfrice, pz0i,  pSc, phf, pDi  )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  CdN_f_LG15  ***
      !!
      !!        GENERAL FORM OF EQUATION 21 of Lupkes & Gryanik (2015)
      !!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !!
      !! ** Purpose :    Computes the "form" contribution of the neutral air-ice
      !!                 drag referenced at 10m to make it dependent on edges at
      !!                 leads, melt ponds and flows (to be added to the "skin"
      !!                 contribution. After some
      !!                 approximations, this can be resumed to a dependency on
      !!                 ice concentration.
      !!
      !! ** References : Lupkes & Gryanik (2015)
      !!
      !!----------------------------------------------------------------------
      REAL(wp),                     INTENT(in )          :: pzu    ! reference height                       [m]
      REAL(wp), DIMENSION(:,:), INTENT(in)           :: pfrice ! ice concentration [fraction]  => at_i_b  ! NOT USED if pSc, phf and pDi all provided...
      REAL(wp), DIMENSION(:,:), INTENT(in)           :: pz0i   ! roughness length over ICE  [m] (in LU12, it's over water ???)
      REAL(wp), DIMENSION(:,:), INTENT(in), OPTIONAL :: pSc    ! shletering function [0-1] (Sc->1 for large distance between floes, ->0 for small distances)
      REAL(wp), DIMENSION(:,:), INTENT(in), OPTIONAL :: phf    ! mean freeboard of floes    [m]
      REAL(wp), DIMENSION(:,:), INTENT(in), OPTIONAL :: pDi    ! cross wind dimension of the floe (aka effective edge length for form drag)   [m]
      REAL(wp), DIMENSION(SIZE(pfrice,1),SIZE(pfrice,2))                       :: CdN_f_LG15  ! neutral FORM drag coefficient contribution over sea-ice
      !!----------------------------------------------------------------------
      LOGICAL :: l_known_Sc=.FALSE., l_known_hf=.FALSE., l_known_Di=.FALSE.
      REAL(wp) :: ztmp, zrlog, zfri, zfrw, zSc, zhf, zDi
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------
      l_known_Sc    = PRESENT(pSc)
      l_known_hf    = PRESENT(phf)
      l_known_Di    = PRESENT(pDi)

      DO jj = 1, SIZE(pfrice,2)
         DO ji = 1, SIZE(pfrice,1)

            zfri = pfrice(ji,jj)
            zfrw = (1._wp - zfri)

            IF(l_known_Sc) THEN
               zSc = pSc(ji,jj)
            ELSE
               !! Sc parameterized in terms of A (ice fraction):
               zSc = zfrw**(1._wp / ( 10._wp * rBeta_0 ))   ! Eq.(31)
               !PRINT *, 'LOLO: Sc PARAMETERIZED !!! =>', zSc
            END IF

            IF(l_known_hf) THEN
               zhf = phf(ji,jj)
            ELSE
               !! hf parameterized in terms of A (ice fraction):
               zhf = rhmax_0*zfri + rhmin_0*zfrw  ! Eq.(25)
               !PRINT *, 'LOLO: hf PARAMETERIZED !!! =>', zhf
            END IF

            IF(l_known_Di) THEN
               zDi = pDi(ji,jj)
            ELSE
               !! Di parameterized in terms of A (ice fraction):
               ztmp = 1._wp / ( 1._wp - (rDmin_0/rDmax_0)**(1._wp/rBeta_0) )   ! A* Eq.(27)
               zDi =  rDmin_0 * ( ztmp/(ztmp - zfri) )**rBeta_0                !    Eq.(26)
               !PRINT *, 'LOLO: Di PARAMETERIZED !!! =>', zDi
            END IF

            ztmp  = 1._wp/pz0i(ji,jj)
            zrlog = LOG(zhf*ztmp/2.718_wp) / LOG(pzu*ztmp)  !LOLO: adding number "e" !!!

            !#bug reported by V. Guemas: CdN_f_LG15(:,:) = 0.5_wp* 0.4_wp * zrlog*zrlog * zSc*zSc * zhf/zDi * zfri  ! Eq.(21) Lukes & Gryanik (2015)
            CdN_f_LG15(:,:) = 0.5_wp* 0.4_wp * zrlog*zrlog * zSc * zhf/zDi * zfri  ! Eq.(21) Lukes & Gryanik (2015)
            !!                   1/2      Ce
         END DO
      END DO
   END FUNCTION CdN_f_LG15



   FUNCTION CdN_f_LG15_light( pzu, pfrice, pz0w )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  CdN_f_LG15_light  ***
      !!
      !! ** Purpose :    Computes the "form" contribution of the neutral air-ice
      !!                 drag referenced at 10m to make it dependent on edges at
      !!                 leads, melt ponds and flows (to be added to the "skin"
      !!                 contribution. After some
      !!                 approximations, this can be resumed to a dependency on
      !!                 ice concentration.
      !!
      !! ** References : Lupkes & Gryanik (2015)
      !!
      !!----------------------------------------------------------------------
      REAL(wp),                     INTENT(in) :: pzu    ! reference height [m]
      REAL(wp), DIMENSION(:,:), INTENT(in) :: pfrice ! ice concentration [fraction]  => at_i_b
      REAL(wp), DIMENSION(:,:), INTENT(in) :: pz0w   ! roughness length over water  [m]
      REAL(wp), DIMENSION(SIZE(pfrice,1),SIZE(pfrice,2)) :: CdN_f_LG15_light  ! neutral FORM drag coefficient contribution over sea-ice
      !!----------------------------------------------------------------------
      REAL(wp) :: ztmp, zrlog, zfri
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------
      DO jj = 1, SIZE(pfrice,2)
         DO ji = 1, SIZE(pfrice,1)

            zfri = pfrice(ji,jj)

            ztmp = 1._wp / pz0w(ji,jj)
            zrlog = LOG( 10._wp * ztmp ) / LOG( pzu * ztmp ) ! part of (Eq.46)

            CdN_f_LG15_light(:,:) = rce10_i_0 *zrlog*zrlog * zfri * (1._wp - zfri)**rbeta_0  ! (Eq.46)  [ index 1 is for ice, 2 for water ]

         END DO
      END DO
   END FUNCTION CdN_f_LG15_light



   !!======================================================================
END MODULE mod_cdn_form_ice
