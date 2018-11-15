PROGRAM cx_vs_wind_test

   USE mod_const
   USE mod_thermo
   USE mod_blk_coare
   USE mod_blk_ncar
   USE mod_blk_ecmwf

   IMPLICIT NONE


   REAL(wp), PARAMETER :: &
      &                zt =  2. ,  &
      &                zu = 10. ,  &
      &             wind_max = 60.

   REAL(wp), DIMENSION(7), PARAMETER :: vrh = (/ 0.7 , 0.75 , 0.8 , 0.85 , 0.9 , 0.95 , 1. /)

   LOGICAL, PARAMETER :: &
      &   ldebug = .FALSE.

   REAL(wp), PARAMETER ::        &
      &   wdebug = 4.2

   INTEGER, PARAMETER :: &
      &   ndt    = 75,   &
      &   n_w    = 1201, &
      &   n_q    = 10
   !!
   CHARACTER(len=10) :: calgo = 'coare'
   CHARACTER(len=4)  :: csst, stab
   !!
   CHARACTER(len=100) :: &
      &   cextra = '', &
      &   cf_cd, cf_ce, cf_ch, cf_ac, cf_ublk, cf_us, cf_z0

   REAL(wp) :: Dtv, dt

   REAL(wp), DIMENSION(1,1) :: sst, sstk, qsat_sst, sst_v, slp, t10, w10, q10, U_bulk, &
      &  z0, z0t, z0q, us, ts, qs, &
      &  mCd, mCe, mCh, mzeta, mz0, mz0t, mz0q, mus, mts, mqs, &
      &  Cd, Ce, Ch, dw0, dw, zeta, &
      &  ac, mac, mublk, tmp


   INTEGER :: nrh, isst, jh, jdt, jw

   REAL(wp), DIMENSION(n_w) :: &
      &   t_w10

   REAL(wp), DIMENSION(2,1) :: v_tq
   !!
   REAL(wp), DIMENSION(:,:), ALLOCATABLE :: &
      &     v_ta_qa, &
      &     XT_a, XQ_a
   !!
   REAL(wp), DIMENSION(ndt) :: &
      &   t_dvt

   REAL(wp), DIMENSION(n_w) :: &
      &   t_zeta,  &
      &   t_cd,    &
      &   t_ch,    &
      &   t_ce,    &
      &   t_z0, t_z0t, t_z0q, &
      &   t_us, t_ts, t_qs, &
      &   t_ac, t_ublk

   IF ( iargc() /= 2 ) THEN
      PRINT *, 'USAGE: cx_vs_wind_test.x <algo (coare/coare35/ncar/ecmwf)> <SST (deg.C)>'
      STOP
   END IF
   
   CALL getarg(1,calgo)
   CALL getarg(2,csst) ; READ(csst,'(i2)') isst ; sst = REAL(isst)
   sstk = sst + rt0



   !! We want a large nb_itt !
   nb_itt = 20
   !nb_itt = 1



   nrh = size(vrh)

   PRINT *, 'nrh =', nrh


   ALLOCATE ( v_ta_qa(2,nrh) , XT_a(ndt,nrh), XQ_a(ndt,nrh) )



   jpi = 1 ; jpj = 1

   slp = Patm



   !!
   !! Table of wind speeds from 0 to wind_max  m/s :
   dw0 = wind_max/(n_w - 1)
   !FORALL (jw = 1:n_ew)  t_w10(jw) = (jw-1)*dw
   !PRINT *, 'dw0 =', dw0 ; STOP
   !!
   t_w10(1) = 0.0
   DO jw = 2, n_w
      dw = 0.25*dw0
      IF ( t_w10(jw-1) >= 5.  ) dw =    dw0
      IF ( t_w10(jw-1) >= 30. ) dw = 2.*dw0
      t_w10(jw) = t_w10(jw-1) + dw(1,1)
   END DO


   PRINT *, ''; PRINT *, 'Wind max =', t_w10(n_w); PRINT *, ''



   !! Table of air/sea virtual potential temperature differences to test Tpv_air - SST :
   t_dvt = (/ -15., -13., -12. , -11. , -10. , -9.5 , -9.  , -8.5 , -8. , -7.5 , -7. , -6.5 ,  &
      &      -6. , -5.5 , -5.  , -4.5 , -4.  , -3.5 , -3. , -2.5 , -2. , -1.75, -1.5 , -1.25, &
      &      -1. ,-0.75 , -0.5 , -0.4, -0.35, -0.3, -0.25 , -0.2, -0.15, -0.12, -0.1 , -0.07, -0.05, &
      &     -0.01, 0.0, 0.01, 0.05, 0.07, 0.1 , 0.12, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.75 ,  &
      &       1. , 1.25, 1.5, 1.75, 2.,  2.5 ,  3.  , 3.5  , 4.  , 4.5  , 5.  , 5.5  ,  &
      &       6. , 6.5  , 7.   , 7.5  ,  8.  , 8.5  , 9.  , 9.5  , 10. , 11.  , 12./)


   DO jdt = 1, ndt
      !!
      CALL FIND_COUPLES(sstk, t_dvt(jdt), nrh, vrh, v_ta_qa)
      !!
      XT_a(jdt,:) = v_ta_qa(1,:)    ! potential temperatures
      XQ_a(jdt,:) = v_ta_qa(2,:)    ! specific humidities
      !!
   END DO

   qsat_sst = 0.98*q_sat(sst + rt0,slp)

   sst_v = (sst + rt0)*(1. + rctv0*qsat_sst) ! virtual SST (K)
   PRINT *, ''; PRINT *, 'Virtual Sea Surface temperature =', REAL(sst_v - rt0,4); PRINT *, ''

   IF ( trim(stab) == 's' ) GOTO 201



   DO jdt = 1, ndt
      !!
      !!
      Dtv = t_dvt(jdt)
      !!
      IF ( ldebug ) THEN
         PRINT *, ''; PRINT *, ''; PRINT *, ''
         PRINT *, 'Air/sea difference in Virtual Potential Temp. =', REAL(Dtv,4)
         PRINT *, '====================================================================='
      END IF
      !!
      IF (Dtv >= 0.) THEN
         WRITE(cf_cd,'("dat/cd_dtv_+",i4.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
         PRINT *, 'cf_cd = ', trim(cf_cd)
         WRITE(cf_ch,'("dat/ch_dtv_+",i4.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
         WRITE(cf_ce,'("dat/ce_dtv_+",i4.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
         !WRITE(cf_zeta,'("dat/zeta_dtv_+",i4.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
         WRITE(cf_z0,'("dat/z0_dtv_+",i4.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
         !WRITE(cf_z0t,'("dat/z0t_dtv_+",i4.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
         !WRITE(cf_z0q,'("dat/z0q_dtv_+",i4.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
         WRITE(cf_ac,'("dat/ac_dtv_+",i4.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
         WRITE(cf_ublk,'("dat/ublk_dtv_+",i4.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
         WRITE(cf_us,'("dat/us_dtv_+",i4.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
         !WRITE(cf_ts,'("dat/ts_dtv_+",i4.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
         !WRITE(cf_qs,'("dat/qs_dtv_+",i4.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
      ELSE
         WRITE(cf_cd,'("dat/cd_dtv_",i5.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
         PRINT *, 'cf_cd = ', trim(cf_cd)
         WRITE(cf_ch,'("dat/ch_dtv_",i5.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
         WRITE(cf_ce,'("dat/ce_dtv_",i5.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
         !WRITE(cf_zeta,'("dat/zeta_dtv_",i5.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
         WRITE(cf_z0,'("dat/z0_dtv_",i5.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
         !WRITE(cf_z0t,'("dat/z0t_dtv_",i5.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
         !WRITE(cf_z0q,'("dat/z0q_dtv_",i5.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
         WRITE(cf_ac,'("dat/ac_dtv_",i5.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
         WRITE(cf_ublk,'("dat/ublk_dtv_",i5.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
         WRITE(cf_us,'("dat/us_dtv_",i5.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
         !WRITE(cf_ts,'("dat/ts_dtv_",i5.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
         !WRITE(cf_qs,'("dat/qs_dtv_",i5.4,"_sst_",i2.2,"_",a,a,".dat")') int(100*Dtv), int(sst), trim(calgo), trim(cextra)
      END IF
      !!
      !!
      !!
      !! Now working on a given Wind :
      DO jw = 1, n_w
         !!
         w10 = t_w10(jw)
         !!
         IF ( ldebug .and. (w10(1,1) == wdebug) ) THEN
            PRINT *, 'Pour dTv =', REAL(t_dvt(jdt),4), ' SST =', REAL(sst,4), &
               &   ' w10 =', REAL(w10,4)
            PRINT *, '     RH           T_a            q_a           Cd           Ce'
         END IF
         !!
         mCd = 0. ;  mCe = 0. ;  mCh = 0. ; mzeta = 0. ; mz0 = 0. ; mz0t = 0. ; mz0q = 0.
         mac = 0. ; mublk = 0.  ;  mus = 0. ;  mts = 0. ;  mqs = 0. ;
         !!
         DO jh = 1, nrh
            !!


            IF ( TRIM(calgo) == 'coare' ) &
               CALL TURB_COARE('3.0', zt, zu, sstk, XT_a(jdt,jh), qsat_sst, XQ_a(jdt,jh), w10, &
               &          Cd, Ch, Ce, t10, q10, U_bulk, xz0=z0, xu_star=us)
            
            IF ( TRIM(calgo) == 'coare35' ) &
               CALL TURB_COARE('3.5', zt, zu, sstk, XT_a(jdt,jh), qsat_sst, XQ_a(jdt,jh), w10, &
               &          Cd, Ch, Ce, t10, q10, U_bulk, xz0=z0, xu_star=us)

            IF ( TRIM(calgo) == 'ncar' ) &
               CALL TURB_NCAR( zt, zu, sstk, XT_a(jdt,jh), qsat_sst, XQ_a(jdt,jh), w10, &
               &          Cd, Ch, Ce, t10, q10, U_bulk, xz0=z0, xu_star=us)

            IF ( trim(calgo) == 'ecmwf' ) &
               CALL TURB_ECMWF(zt, zu, sstk, XT_a(jdt,jh), qsat_sst, XQ_a(jdt,jh), w10, &
               &          Cd, Ch, Ce, t10, q10, U_bulk, xz0=z0, xu_star=us)



            !!
            !! Recalculating eq. charnock parameter:
            tmp = us*us
            ac = z0*grav/tmp - 0.11*grav*visc_air(t10)/(tmp*us)

            !!
            mCd = mCd + Cd/nrh
            mCh = mCh + Ch/nrh
            mCe = mCe + Ce/nrh
            mzeta = mzeta + zeta/nrh
            mz0 = mz0 + z0/nrh
            mz0t = mz0t + z0t/nrh
            mz0q = mz0q + z0q/nrh
            mac = mac + ac/nrh
            mublk = mublk + U_bulk/nrh
            mus = mus + us/nrh
            mts = mts + ts/nrh
            mqs = mqs + qs/nrh
            !!
            IF ( ldebug .and. (w10(1,1) == wdebug) ) &
               PRINT *, REAL(vrh(jh),4), REAL(XT_a(jdt,jh)-rt0,4), &
               &   REAL(XQ_a(jdt,jh)*1000.,4), REAL(1000*Cd,4), REAL(1000*Ce,4)
            !!
         END DO
         !!
         !! Taking the mean over the nrh values (we checked that they only depend on
         !!    dTv for a given SST and wind
         t_cd(jw) = 1000.*mCd(1,1)
         t_ch(jw) = 1000.*mCh(1,1)
         t_ce(jw) = 1000.*mCe(1,1)
         t_zeta(jw) = mzeta(1,1)
         t_z0(jw) = mz0(1,1); t_z0t(jw) = mz0t(1,1); t_z0q(jw) = mz0q(1,1)
         t_ac(jw) = mac(1,1)
         t_ublk(jw) = mublk(1,1)
         t_us(jw) = mus(1,1) ; t_ts(jw) = mts(1,1) ; t_qs(jw) = mqs(1,1)
         !!
         !!
         !!
      END DO
      !!
      !!
      !!
      OPEN(unit = 11, file = cf_cd, status = 'unknown')
      WRITE(11,'("#     Wind               Cd")')
      DO jw = 1, n_w
         WRITE(11,'(1f16.8, " " ,1f16.8)') t_w10(jw), t_cd(jw)
      ENDDO
      CLOSE(11)
      !!
      OPEN(unit = 11, file = cf_ch, status = 'unknown')
      WRITE(11,'("#     Wind               Ch")')
      DO jw = 1, n_w
         WRITE(11,'(1f16.8, " " ,1f16.8)') t_w10(jw), t_ch(jw)
      ENDDO
      CLOSE(11)
      !!
      OPEN(unit = 11, file = cf_ce, status = 'unknown')
      WRITE(11,'("#     Wind               Ce")')
      DO jw = 1, n_w
         WRITE(11,'(1f16.8, " " ,1f16.8)') t_w10(jw), t_ce(jw)
      ENDDO
      CLOSE(11)
      !!
      OPEN(unit = 11, file = cf_us, status = 'unknown')
      WRITE(11,'("#     Wind               u*")')
      DO jw = 1, n_w
         WRITE(11,'(1f16.8, " " ,1f16.8)') t_w10(jw), t_us(jw)
      ENDDO
      CLOSE(11)
      !!
      OPEN(unit = 11, file = cf_z0, status = 'unknown')
      WRITE(11,'("#     Wind               z0")')
      DO jw = 1, n_w
         WRITE(11,'(1f16.8, " " ,1f16.8)') t_w10(jw), t_z0(jw)
      ENDDO
      CLOSE(11)
      !!
      OPEN(unit = 11, file = cf_ac, status = 'unknown')
      WRITE(11,'("#     Wind               u*")')
      DO jw = 1, n_w
         WRITE(11,'(1f16.8, " " ,1f16.8)') t_w10(jw), t_ac(jw)
      ENDDO
      CLOSE(11)
      !!
      OPEN(unit = 11, file = cf_ublk, status = 'unknown')
      WRITE(11,'("#     Wind               U_bulk")')
      DO jw = 1, n_w
         WRITE(11,'(1f16.8, " " ,1f16.8)') t_w10(jw), t_ublk(jw)
      ENDDO
      CLOSE(11)
      !!
   END DO
   !!
   !!
   !!
   !!
   PRINT *, ''; PRINT *, ''; PRINT *, ''
   !!
201 CONTINUE
   !!


   IF ( nb_itt < 20 ) THEN
      PRINT *, 'ACHTUNG!!!! nb_itt was:', nb_itt
   ELSE
      PRINT *, ' Number of itterations used:', nb_itt
   END IF



   STOP

   PRINT *, 'ALLO!!!!'


   !! Finding for what Dtv we get neutral stability, with current sst :
   !! =================================================================
   dt = 0.5
   w10 = 4.
   !!
   !! Starting from unstable state (zeta<0)
   Dtv = -5.   ! Air Dtv deg. cooler than sea
   !!
   DO WHILE ( ABS(dt) > 1.E-9 )
      !!
      Dtv = Dtv + dt
      !!
      !PRINT *, 'Dtv, dt =', Dtv, dt
      !! Couple t,q for rh=80% :
      CALL FIND_COUPLES(sstk, Dtv, 1, (/ 0.8_wp /), v_tq)
      !PRINT *, 't2, q2 =', v_tq(1,1)-rt0, 1000.*v_tq(2,1)
      !!

      IF ( trim(calgo) == 'coare' ) &
         CALL TURB_COARE('3.0', zt, zu, sstk, v_tq(1,1), qsat_sst, v_tq(2,1), w10, &
         &               Cd, Ch, Ce, t10, q10, U_bulk)

      IF ( TRIM(calgo) == 'coare35' ) &
         CALL TURB_COARE('3.5', zt, zu, sstk, v_tq(1,1), qsat_sst, v_tq(2,1), w10, &
         &               Cd, Ch, Ce, t10, q10, U_bulk)

      IF ( trim(calgo) == 'ncar' ) &
         CALL TURB_NCAR(zt, zu, sstk, v_tq(1,1), qsat_sst, v_tq(2,1), w10, &
         &             Cd, Ch, Ce, t10, q10, U_bulk)

      IF ( trim(calgo) == 'ecmwf' ) &
         CALL TURB_ECMWF(zt, zu, sstk, v_tq(1,1), qsat_sst, v_tq(2,1), w10, &
         &             Cd, Ch, Ce, t10, q10, U_bulk)


      !!
      !PRINT *, 'zeta =', zeta; PRINT *, ''
      IF( zeta(1,1) > 0. ) THEN ! On est passe du cote stable
         Dtv = Dtv - dt
         dt = dt/1.2  ! decreasing temperature increment
      END IF
      !!
   END DO
   !!
   PRINT *, 'Neutral state for Dtv =', REAL(Dtv,4)
   PRINT *, 'With sst  =', REAL(sst,4)
   PRINT *, 'With wind =', REAL(w10,4)
   !!
   CALL FIND_COUPLES(sst+rt0, Dtv, 1, (/ 0.8_wp /), v_tq)
   !!
   PRINT *, 'Corresponds (for RH=80%) to t2, q2   =', REAL(v_tq(1,1)-rt0,4), REAL(1000.*v_tq(2,1),4)
   PRINT *, 'Corresponds (for RH=80%) to t10, q10 =', REAL(t10-rt0,4), REAL(1000.*q10,4)





CONTAINS




   SUBROUTINE FIND_COUPLES(Ts, dvt, nh, trh, t_ta_qa)

      !!#######################################################################
      !!
      !!    Ts         surface temperature          (K)
      !!    dvt        difference virt. pot. (Tvp_a - Ts)       (K)
      !!    nh         number of different relative humidity values tested
      !!    trh        table containing the nh humidities     (ratio)
      !!    t_ta_qa    array contening the nh couples possible   (K,kg/kg)
      !!
      !!#######################################################################



      REAL(wp),    DIMENSION(jpi,jpj),  INTENT(in)  :: Ts
      REAL(wp),                     INTENT(in)  :: dvt
      INTEGER, INTENT(in)                   :: nh
      REAL(wp),    DIMENSION(nh)  , INTENT(in)  :: trh
      REAL(wp),    DIMENSION(2,nh), INTENT(out) :: t_ta_qa
      !!
      REAL(wp), PARAMETER :: reps = 1.E-7


      INTEGER :: jh, nitt

      REAL(wp) :: rdiff
      REAL(wp), DIMENSION(jpi,jpj)  :: RH_a, T_a, q_a, Tv_a, T_old, zqsat_sst, sstv, zslp




      zslp = Patm

      !! Virtual Sea Surface temperature:
      zqsat_sst = q_sat(Ts, zslp)
      sstv = Ts*(1. + rctv0*zqsat_sst) ! virtual SST (K)
      !!
      !!
      DO jh = 1, nh
         !!
         RH_a = trh(jh)
         !!
         !Tv_a = Ts + dvt
         Tv_a = sstv + dvt
         !!
         !! First guess of Ta :
         T_a = Tv_a
         !!
         !!
         nitt = 0 ;  rdiff = 10
         !!
         DO WHILE ( rdiff > reps )
            !!
            nitt = nitt + 1 ;  T_old = T_a
            !!
            q_a = q_air_rh(RH_a, T_a, zslp)
            !!
            !! Updating new air temperature:
            T_a = Tv_a/(1. + rctv0*q_a)
            !!
            rdiff = ABS(T_a(1,1) - T_old(1,1))
            !!
         END DO
         !!
         !!
         t_ta_qa(1,jh) = T_a(1,1)
         t_ta_qa(2,jh) = q_a(1,1)
         !!
      END DO
      !!
   END SUBROUTINE FIND_COUPLES
   !!
   !!
   !!
END PROGRAM cx_vs_wind_test
