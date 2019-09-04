! AeroBulk / 2016 / L. Brodeau
!
!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
!
!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
!   significant effects of some approximations in the bulk parameterizations of
!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
!
!
MODULE mod_wl_coare3p6
   !!====================================================================================
   !!       Warm-Layer correction of SST (if needed)
   !!
   !!       Routine "wl_coare3p6" maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, 2019
   !!
   !!====================================================================================
   USE mod_const   !: physical and othe constants
   USE mod_phymbl  !: thermodynamics

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: WL_COARE3P6

   REAL(wp), PARAMETER :: rich   = 0.65_wp   !: critical Richardson number
   REAL(wp), PARAMETER :: z_sst  = 1._wp     !: depth at which bulk SST is taken...
   REAL(wp), PARAMETER :: dz_max = 19._wp    !: maximum depth of warm layer (adjustable)
   REAL(wp), PARAMETER :: Qabs_thr = 50._wp  !: Threshold for heat flux absorbed in WL

   !LOGICAL, PUBLIC, SAVE :: l_wl_c36_never_called

CONTAINS

   SUBROUTINE WL_COARE3P6( pQsw, pQnsol, pTau, pSST, pdT, pTau_ac, pQ_ac, plon, isd, rdt, &
      &                    Hwl )
      !!---------------------------------------------------------------------
      !!
      !!  Cool-Skin Warm-Layer scheme according to COARE 3.6 (Fairall et al, 2019)
      !!     ------------------------------------------------------------------
      !!
      !!  **   INPUT:
      !!
      !!     *pQsw*       net solar radiative flux to the ocean
      !!     *pQnsol*     net non-solar heat flux to the ocean
      !!     *pTau*       surface wind stress
      !!     *pSST*       SST
      !!
      !!   **  INPUT/OUTPUT:
      !!     *pdT*  :
      !!
      !!------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pQsw     ! net solar radiation into the sea [W/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pQnsol   ! net non-solar heat flux into the sea [W/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pTau     ! wind stress [N/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pSST     ! bulk SST at depth z_sst [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pdT      ! dT due to warming at depth of pSST such that pSST_true = pSST + pdT
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pTau_ac  ! time integral / accumulated momentum Tauxdt => [N.s/m^2] (reset to zero every midnight)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pQ_ac    ! time integral / accumulated heat stored by the warm layer Qxdt => [J/m^2] (reset to zero every midnight)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: plon     ! longitude ! lolo
      INTEGER ,                     INTENT(in)    :: isd      ! current UTC time, counted in second since 00h of the current day
      REAL(wp),                     INTENT(in)    :: rdt      ! time step in seconds
      !!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out), OPTIONAL :: Hwl    ! depth of warm layer
      !
      !
      INTEGER :: ji,jj
      !
      REAL(wp) :: dT_wl, dz_wl, Qabs, zfs
      REAL(wp) :: zqac, ztac
      REAL(wp) :: Al, qjoule
      REAL(wp) :: ctd1, ctd2

      REAL(wp) :: rlag_gw_h  ! local solar time lag in hours   / Greenwich meridian (lon==0) => ex: ~ -10.47 hours for Hawai
      
      INTEGER  :: ilag_gw_s, &  ! local solar time LAG in seconds / Greenwich meridian (lon==0) => ex: ~ INT( -10.47*3600. ) seconds for Hawai
         &        isd_sol,   &  ! local solar time in number of seconds since local solar midnight
         &        jl
      INTEGER  :: jamset
      
      !! INITIALIZATION:

      !pTau=Bx(2)  !stress
      !!Qsen=Bx(3)   !sensible heat flux
      !!Qlat=Bx(7)   !latent heat flux - Webb corrected
      !!Qrain=Bx(34)  !rain heat flux
      ! => pQnsol == Qlw+Qsen+Qlat+Qrain
      !dT_skin=Bx(17)    !cool skin

      !pQ_ac   = 0._wp       ! accumulates heat from integral
      !pTau_ac = 0._wp     ! accumulates stress from integral

      ! Just in case not initialized to 0...
      !IF ( l_wl_c36_never_called ) THEN
      !   pTau_ac = 0.
      !   pQ_ac   = 0.
      !END IF

      pdT = 0._wp         ! dT initially set to 0._wp

      dT_wl = 0._wp       ! total warming (amplitude) in warm layer
      dz_wl = dz_max      ! initial depth set to max value
      Qabs  = 0._wp        ! total heat absorped in warm layer
      zfs   = 0.5_wp         ! initial value of solar flux absorption

      DO jj = 1, jpj
         DO ji = 1, jpi

            !PRINT *, ' *** jj =', jj

            !! Need to know the local solar time from longitude and isd
            !! *********************************************************
            PRINT *, ''
            rlag_gw_h = -1._wp * MODULO( ( 360._wp - MODULO(plon(ji,jj),360._wp) ) / 15._wp , 24._wp )
            rlag_gw_h = -1._wp * SIGN( MIN(ABS(rlag_gw_h) , ABS(MODULO(rlag_gw_h,24._wp))), rlag_gw_h + 12._wp )            
            !PRINT *, ' Lag in hours / Greenwich for local solar time =', rlag_gw_h            
            ilag_gw_s = INT( rlag_gw_h*3600._wp )
            isd_sol = MODULO( isd + ilag_gw_s , 24*3600 )
            !PRINT *, '     UTC     time in seconds:', isd
            !PRINT *, ' Local solar time in seconds:', isd_sol
            !PRINT *, '     UTC     time in hours:',   REAL(isd    ,wp)/3600._wp
            PRINT *, '  [WL_COARE3P6] Local solar time in hours:',   REAL(isd_sol,wp)/3600._wp            

            jamset = 0
            !jump = 1

            !**********************************************************
            !******************  setup read data loop  ****************



            !*****  variables for warm layer  ***
            Al   = 2.1e-5*(pSST(ji,jj)-rt0 + 3.2_wp)**0.79
            ctd1 = SQRT(2._wp*rich*rCp0_w/(Al*grav*rho0_w))        !mess-o-constants 1
            ctd2 = SQRT(2._wp*Al*grav/(rich*rho0_w))/(rCp0_w**1.5) !mess-o-constants 2

            !********************************************************
            !****  Compute apply warm layer  correction *************
            !********************************************************

            !intime  = yd-fix(yd)
            !loc     = (lonx+7.5)/15
            !chktime = loc+intime*24
            !IF (chktime>24) chktime = chktime-24
            !it_n = (chktime-24*fix(chktime/24))*3600

            !IF (icount>1) THEN                                  !not first time thru

            !IF ( (it_n <= 21600).OR.(jump == 0) ) THEN  ! (21600 == 6am)

            !   jump = 0

            IF (isd_sol < rdt ) THEN    !re-zero at midnight ! LOLO improve: risky if real midnight (00:00:00) is not a time in vtime...
               PRINT *, '  [WL_COARE3P6] MIDNIGHT RESET !!!!, isd_sol =>', isd_sol
               jamset         = 0
               zfs            = 0.5_wp
               dz_wl          = dz_max
               pTau_ac(ji,jj) = 0._wp
               ztac           = 0._wp
               pQ_ac(ji,jj)   = 0._wp
               zqac           = 0._wp
               dT_wl          = 0._wp
            END IF

            IF ( isd_sol >= 21600 ) THEN  ! (21600 == 6am)

               PRINT *, '  [WL_COARE3P6] WE DO WL !'
               PRINT *, '  [WL_COARE3P6] isd_sol, pTau, pSST, pdT =', isd_sol, REAL(pTau(ji,jj),4), REAL(pSST(ji,jj),4), REAL(pdT(ji,jj),4)

               !************************************
               !****   set warm layer constants  ***
               !************************************

               Qabs = zfs*pQsw(ji,jj) - pQnsol(ji,jj)       ! tot heat absorbed in warm layer

               PRINT *, '  [WL_COARE3P6] rdt,  pQsw, pQnsol, Qabs =', rdt, REAL(pQsw(ji,jj),4), REAL(pQnsol(ji,jj),4), REAL(Qabs,4)

               IF ( (Qabs >= Qabs_thr).OR.(jamset == 1) ) THEN         ! Check for threshold

                  PRINT *, '  [WL_COARE3P6] pTau_ac, pQ_ac =', REAL(pTau_ac(ji,jj),4), REAL(pQ_ac(ji,jj),4)

                  jamset = 1                                         ! indicates threshold crossed
                  !pTau_ac(ji,jj) = pTau_ac(ji,jj) + MAX(.002_wp , pTau(ji,jj))*rdt      ! momentum integral
                  ztac = pTau_ac(ji,jj) + MAX(.002_wp , pTau(ji,jj))*rdt      ! updated momentum integral

                  IF ( pQ_ac(ji,jj) + Qabs*rdt > 0._wp ) THEN         !check threshold for warm layer existence
                     !******************************************
                     ! Compute the absorption profile
                     !******************************************
                     DO jl = 1, 5                           !loop 5 times for zfs
                        zfs = 1. - ( 0.28*0.014*(1. - EXP(-dz_wl/0.014)) + 0.27*0.357*(1. - EXP(-dz_wl/0.357)) &
                           &        + 0.45*12.82*(1-EXP(-dz_wl/12.82)) ) / dz_wl
                        qjoule = (zfs*pQsw(ji,jj) - pQnsol(ji,jj))*rdt
                        zqac = pQ_ac(ji,jj) + qjoule ! updated heat absorbed
                        IF (zqac > 0._wp)  dz_wl = MIN( dz_max , ctd1*ztac/SQRT(zqac)) !Compute warm-layer depth
                     END DO

                  ELSE
                     !***********************
                     ! Warm layer wiped out
                     !***********************
                     zfs    = 0.75
                     dz_wl  = dz_max
                     qjoule = (zfs*pQsw(ji,jj) - pQnsol(ji,jj))*rdt
                     zqac   = pQ_ac(ji,jj) + qjoule ! updated heat absorbed
                  END IF


                  !*******  compute dt_warm  ******
                  !LOLOnew: dT_wl = ctd2*MAX(0._wp,zqac)**1.5/ztac
                  IF (zqac > 0._wp) THEN
                     dT_wl = ctd2*zqac**1.5/ztac
                  ELSE
                     dT_wl = 0._wp
                  END IF

               END IF ! IF ( (Qabs >= Qabs_thr).OR.(jamset == 1) )


               IF (dz_wl < z_sst) THEN           !Compute warm layer correction
                  pdT(ji,jj) = dT_wl
               ELSE
                  pdT(ji,jj) = dT_wl*z_sst/dz_wl
               END IF
                              
            END IF ! IF ( isd_sol >= 21600 ) THEN  ! (21600 == 6am)

            !END IF !  IF ( (isd_sol<=21600).OR.(jump==0) )  end 6am start first time thru


            IF ( (Qabs >= Qabs_thr).AND.(isd_sol >= 21600) ) THEN
               pQ_ac(ji,jj)   = zqac ! Updating pQ_ac, heat integral
               pTau_ac(ji,jj) = ztac !
            END IF

            
            !PRINT *, ' *** END => pdT =', pdT(ji,jj)
            
            IF ( PRESENT(Hwl) ) Hwl(ji,jj) = dz_wl

         END DO
         PRINT *, ''
      END DO
      
      !END IF  !  IF (icount>1)


      !it_b = isd_sol

      !************* output from routine  *****************************
      !Output:   A=[usr tau hsb hlb hbb hsbb hlwebb tsr qsr zo  zot zoq  Cd  Ch  Ce  L zet dT_skinx dqerx tkt Urf Trf Qrf RHrf UrfN Qlw  Le rhoa UN U10 U10N Cdn_10 Chn_10 Cen_10 RF Qs Evap T10 Q10 RH10 ug Whf Edis]
      !              1   2   3   4   5    6     7    8   9  10  11   12  13  14  15 16  17   18   19    20  21  22  23  24   25   26  27  28  29  30  31    32     33     34   35 36  37   38  39  40  41  42  43
      !pTs = pSST + pdT
      !Bx = coare36vn_zrf(u,zu(ibg),t,zt(ibg),rh,zq(ibg),P,ts,Rs,Rl,latx,zi,rain,Ss(ibg),cp(ibg),sigH(ibg),zrf_u(ibg),zrf_t(ibg),zrf_q(ibg))
      !pTau = Bx(2)            !hold stress
      !Qsen = Bx(3)             !hold shf
      !Qlat = Bx(4)             !hold lhf - use Webb corrected value
      !dT_skin = Bx(17)
      !Qrain = Bx(34)            !hold rain flux

      !B(ibg,1) = dT_wl   ! warming across entire warm layer deg.C
      !B(ibg,2) = dz_wl   ! warm layer thickness m
      !B(ibg,3) = pdT     ! heating at selected depth

      !icount = icount+1


      !l_wl_c36_never_called = .FALSE.

   END SUBROUTINE WL_COARE3P6

END MODULE mod_wl_coare3p6

!**************************************************
! Recompute fluxes with warm layer
!**************************************************
!clear Bx
!pSST=pSST+B(:,3)'
!Bx=coare36vn_zrf(Ur,zu,Tair,zt,RH,zq,Pair,pSST,Solar,IR,Lat,zi,Rainrate,Ss,cp,sigH,zrf_u,zrf_t,zrf_q)

!B=[Bx B]    !Add the warm layer variables

!************* output from routine  *****************************
!Output:   A=[usr tau hsb hlb hbb hsbb hlwebb tsr qsr zo  zot zoq  Cd  Ch  Ce  L zet dT_skinx dqerx tkt Urf Trf Qrf RHrf UrfN Qlw  Le rhoa UN U10 U10N Cdn_10 Chn_10 Cen_10 RF Qs Evap T10 Q10 RH10 ug Whf Edis]
!              1   2   3   4   5    6     7    8   9  10  11   12  13  14  15 16  17   18   19    20  21  22  23  24   25   26  27  28  29  30  31    32     33     34   35 36  37   38  39  40  41  42  43
!       dt_warm dz_wl pdT
!          44     45    46



!***********   input data **************
!   yday=       day-of-year
!       Ur=                     wind vector maginitude (m/s) relative to water at height zu
!       zu=                     height (m) of wind measurement
!       Tair=           air temp (degC)at height zt
!       zt=                     height (m) of air temperature measurement
!       RH=                     relative humidity (!) at height zq
!       zq=                     height (m) of air humidity measurement
!       Pair=           air pressure (mb)
!       pSST=       bulk surface sea temp (degC) at z_sst
!       Solar=          downward solar flux (w/m^2) defined positive down
!       IR=                     downward IR flux (w/m^2) defined positive down
!       Lat=            latitude (deg N=+)
!       Lon=            longitude (deg E=+)
!   zi=         inversion height (m)
!       Rainrate=       rain rate (mm/hr)
!       z_sst        depth (m) of water temperature measurement
!   Ss =        sea surface salinity (PSU)
!   cp =        phase speed of dominant waves (m/s)
!   sigH =      significant wave height (m)
!   zu, zt, zq heights of the observations (m)
!   zrf_u, zrf_t, zrf_q  reference height for profile.  Use this to compare observations at different heights
!
!
!********** output data  ***************
!Outputs
!From coare36vn_zrf
!Output:   A=[usr tau hsb hlb hbb hsbb hlwebb tsr qsr zo  zot zoq  Cd  Ch  Ce  L zet dT_skinx dqerx tkt Urf Trf Qrf RHrf UrfN Qlw  Le rhoa UN U10 U10N Cdn_10 Chn_10 Cen_10 RF Qs Evap T10 Q10 RH10 ug Whf Edis];
!              1   2   3   4   5    6     7    8   9  10  11   12  13  14  15 16  17   18   19    20  21  22  23  24   25   26  27  28  29  30  31    32     33     34   35 36  37   38  39  40  41  42  43
!where
!
!   usr = friction velocity that includes gustiness (m/s)
!   tau = wind stress (N/m^2)
!   hsb = sensible heat flux into ocean (W/m^2)
!   hlb = latent heat flux into ocean (W/m^2)
!   hbb = buoyany flux into ocean (W/m^2)
!   hsbb = "sonic" buoyancy flux measured directly by sonic anemometer
!   tsr = temperature scaling parameter (K)
!   qsr = specific humidity scaling parameter (g/Kg)
!   zot = thermal roughness length (m)
!   zoq = moisture roughness length (m)
!   Cd = wind stress transfer (drag) coefficient at height zu
!   Ch = sensible heat transfer coefficient (Stanton number) at height zu
!   Ce = latent heat transfer coefficient (Dalton number) at height zu
!    L = Obukhov length scale (m)
!  zet = Monin-Obukhov stability parameter zu/L
! dT_skin = cool-skin temperature depression (degC)
! dqer = cool-skin humidity depression (degC)
!  tkt = cool-skin thickness (m)
!  Urf = wind speed at reference height (user can select height below)
!  Tfr = temperature at reference height
!  Qfr = specific humidity at reference height
! RHfr = relative humidity at reference height
! UrfN = neutral value of wind speed at reference height
!  Qlw = Upwelling IR radiation computed by COARE
!   Le = latent heat of vaporization
! rhoa = density of air
!   UN = neutral value of wind speed at zu
!  U10 = wind speed adjusted to 10 m
! UN10 = neutral value of wind speed at 10m
!Cdn_10 = neutral value of drag coefficient at 10m
!Chn_10 = neutral value of Stanton number at 10m
!Cen_10 = neutral value of Dalton number at 10m
!Rf     = Rain heat flux (W/m^2)
!Qs     = surface specific humidity (g/g)
!Evap   = evaporation rate (mm/h)
!T10    = air temperature at 10m
!Q10    = air specific humidity at 10m
!RH10   = air relative humidity at 10m
!uq     = gustiness velocity (m/s)
!Whf    = whitecap fraction
!Edis   = energy dissipated by wave breaking (W/m^2)



!From WarmLayer
! 44: dT_wl - warming across entire warm layer degC
! 45: dz_wl - warm layer thickness m
! 46: pdT   - dT due to warming at depth of pSST such that pSST_true = pSST + pdT

!icount=1
!*********************  housekeep variables  ********
! Call coare35vn to get initial flux values
!Bx = coare36vn_zrf(Ur(1),zu(1),Tair(1),zt(1),RH(1),zq(1),Pair(1),pSST(1),Solar(1),IR(1),Lat(1),zi(1),Rainrate(1),Ss(1),cp(1),sigH(1),zrf_u(1),zrf_t(1),zrf_q(1))
