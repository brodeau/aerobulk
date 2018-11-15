! AeroBulk / 2015 / L. Brodeau (brodeau@gmail.com)
! https://sourceforge.net/p/aerobulk

MODULE mod_blk_neutral_10m

   !!====================================================================================
   !!
   !!   Momentum, Latent and sensible heat exchange coefficients
   !!          In neutral conditions / at 10m
   !!
   !!            Author: Laurent Brodeau, 2016, brodeau@gmail.com
   !!
   !!====================================================================================

   USE mod_const   !: physical and othe constants

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: TURB_NEUTRAL_10M

   !! COARE own values for given constants:
   REAL(wp), PARAMETER :: &
      &   zu  = 10.,  &     ! we're at 10m !
      &  charn0_max = 0.028



CONTAINS

   SUBROUTINE turb_neutral_10m( calgo, U_N10, CdN10, ChN10, CeN10)
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_neutral_10m  ***
      !!
      !!            2015: L. Brodeau (brodeau@gmail.com)
      !!
      !! ** Purpose :   Computes turbulent neutra transfert coefficients at 10m
      !!                from the neutral wind speed at 10m
      !!                according to ...
      !!
      !!======================================================================================
      !!
      !! INPUT :
      !! -------
      !!    *  calgo : name of the algo to use
      !!    *  U_N10 : scalar wind speed at 10m                                [m/s]
      !!
      !! OUTPUT :
      !! --------
      !!    *  CdN10     : drag coefficient
      !!    *  ChN10     : sensible heat coefficient
      !!    *  CeN10     : evaporation coefficient
      !!
      !!============================================================================


      CHARACTER(len=*), INTENT(in)                ::   calgo
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   U_N10       ! relative wind module at zu            [m/s]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   CdN10       ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   ChN10       ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   CeN10       ! transfert coefficient for evaporation   (Q_lat)

      INTEGER :: j_itt

      REAL(wp), DIMENSION(:,:), ALLOCATABLE  ::  &
         &  u_star, &
         &  z0, z0t, z0q

      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   ztmp0, ztmp1, Ub

      !LOGICAL :: lreturn_z0 = .FALSE., lreturn_ustar = .FALSE.


      ALLOCATE ( u_star(jpi,jpj), &
         &     z0(jpi,jpj), z0t(jpi,jpj), z0q(jpi,jpj), Ub(jpi,jpj),         &
         &     ztmp0(jpi,jpj), ztmp1(jpi,jpj) )




      Ub = MAX(U_N10, 0.1_wp)


      IF ( (TRIM(calgo) == 'coare30').OR.(TRIM(calgo) == 'coare35').OR.(TRIM(calgo) == 'ecmwf') ) THEN

         !! First guess of CdN10
         CdN10 = 8.575E-5*Ub + 0.657E-3  ! from my curves


         !! ITERATION BLOCK
         DO j_itt = 1, nb_itt

            !! Need to know u*
            !! --------------
            u_star = Ub*SQRT(CdN10)


            !! Need to know z0 and z0t
            !! -----------------------

            IF ( TRIM(calgo) == 'coare35' ) THEN
               !! COARE 3.5: Charnock parameter is computed from the neutral wind speed at 10m: Eq. 13 (Edson al. 2013)
               ztmp0 = MIN( 0.0017_wp*Ub - 0.005_wp , charn0_max)  ! alpha Charnock parameter (Eq. 13 Edson al. 2013)
               !ztmp0 = MAX( ztmp0 , 0. )

            ELSEIF ( TRIM(calgo) == 'coare30' ) THEN
               ztmp0 = alfa_charn(Ub)

            ELSEIF ( TRIM(calgo) == 'ecmwf' ) THEN
               ztmp0 = 0.018

            ELSE

               PRINT *, TRIM(calgo)//' => Unknow algo !' ; STOP

            END IF

            !IF ( (jpi==1).AND.(jpj==1) ) PRINT *, '   *** alpha = ', ztmp0

            !! Roughness lengthes z0, momentum: z0t (z0q = z0t) :
            z0   = ztmp0*u_star*u_star/grav + 0.11*nu0_air/u_star ! Roughness length (eq.6)


            !! Compute drag coefficient at 10m
            !! -------------------------------

            ztmp0 = LOG(zu/z0)

            CdN10 = vkarmn*vkarmn / ( ztmp0*ztmp0 )

            !IF ( (jpi==1).AND.(jpj==1) ) PRINT *, '   *** CdN10 =', CdN10*1000.
            IF ( (jpi==1).AND.(jpj==1) ) PRINT *, ''

         END DO


         IF ( TRIM(calgo) == 'coare30' ) THEN

            ! Re_r: roughness Reynolds number:
            ztmp1 = z0*u_star/nu0_air

            !! Scalar roughness length z0t:
            z0t = MIN( 1.1E-4_wp , 5.5E-5_wp*ztmp1**(-0.6_wp) )     ! Scalar roughness for Theta and q (Fairall al 2003, eq.28)
            z0q = z0t

         END IF

         IF ( TRIM(calgo) == 'coare35' ) THEN

            ! Re_r: roughness Reynolds number:
            ztmp1 = z0*u_star/nu0_air

            !! Scalar roughness length z0t:
            !! Chris Fairall, Jim Edscon, private communication, March 2016 / COARE 3.5 :
            z0t   = MIN( 1.6e-4_wp , 5.8E-5_wp*ztmp1**(-0.72_wp) ) ! These thermal roughness lengths give Stanton and
            !!                                            ! Dalton numbers that closely approximate COARE3.0
            z0q = z0t

         END IF


         IF ( TRIM(calgo) == 'ecmwf' ) THEN
            ztmp1  = nu0_air/u_star
            z0t    = 0.40*ztmp1                              ! eq.3.26, Chap.3, p.34, IFS doc - Cy31r1
            z0q    = 0.62*ztmp1
         END IF



         
         PRINT *, ' *** Algo = ',trim(calgo)
         
         ztmp1 = LOG(zu/z0t)
         ChN10 = vkarmn*vkarmn / ( ztmp0*ztmp1 )     ! (ztmp0 = LOG(zu/z0))

         ztmp1 = LOG(zu/z0q)
         CeN10 = vkarmn*vkarmn / ( ztmp0*ztmp1 )     ! (ztmp0 = LOG(zu/z0))








      ELSEIF ( TRIM(calgo) == 'ncar' ) THEN

         Ub = MAX(U_N10, 0.5_wp)

         CdN10  = cd_neutral_10m( Ub )

         !Ch_n10  = 1.e-3*sqrt_Cd_n10*(18.*stab + 32.7*(1. - stab))  ! L&Y 2004 eq. (6c-6d)

         ztmp0 = SQRT(CdN10)

         ChN10  = 1.e-3*ztmp0*32.7  ! UNSTABLE CASE !!!      L&Y 2004 eq. (6c-6d)

         CeN10  = 1.e-3 * (34.6 * ztmp0)  ! L&Y 2004 eq. (6b)



      ELSE

         PRINT *, 'ERROR: algorithm '//TRIM(calgo)//' is not supported yet!'
         PRINT *, ''
         STOP

      END IF



      DEALLOCATE ( u_star, z0, z0t, Ub, ztmp0, ztmp1 )

   END SUBROUTINE turb_neutral_10m






   FUNCTION alfa_charn(dw)
      !!
      !! 2013: L. Brodeau
      !!
      !! Compute Charncok's constant depending on the wind speed
      !!
      !!  Wind below 10 m/s :  alfa = 0.011
      !!  Wind between 10 and 18 m/s : linear increase from 0.011 to 0.018
      !!  Wind greater than 18 m/s :  alfa = 0.018
      !!
      !!
      !!
      REAL(wp), DIMENSION(jpi,jpj) :: alfa_charn
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: dw   ! relative wind speed air/ocean
      !!
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: gt10, gt18
      !!
      ALLOCATE( gt10(jpi,jpj), gt18(jpi,jpj) )
      !!
      !! Charnock's constant, increases with the wind :
      gt10 = 0.5 + SIGN(0.5_wp,(dw - 10.)) ! If dw<10. --> 0, else --> 1
      gt18 = 0.5 + SIGN(0.5_wp,(dw - 18.)) ! If dw<18. --> 0, else --> 1
      !!
      alfa_charn =  (1. - gt10)*0.011    &    ! wind is lower than 10 m/s
         & + gt10*((1. - gt18)*(0.011 + (0.018 - 0.011) &
         & *(dw - 10.)/(18. - 10.)) + gt18*( 0.018 ) ) ! Hare et al. (1999)
      !!
      DEALLOCATE( gt10, gt18 )
      !!
   END FUNCTION alfa_charn






   FUNCTION cd_neutral_10m( pw10 )
      !!
      !!-----------------------------------------------------------------
      !! Estimate of the neutral drag coefficient at 10m as a function
      !! of neutral wind  speed at 10m
      !!
      !! Origin: Large & Yeager 2008 eq.(11a) and eq.(11b)
      !!
      !! Author: L. Brodeau, june 2016
      !!-----------------------------------------------------------------
      !!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pw10           ! scalar wind speed at 10m (m/s)
      REAL(wp), DIMENSION(jpi,jpj)             :: cd_neutral_10m
      !
      INTEGER  ::     ji, jj     ! dummy loop indices
      REAL(wp) :: zgt33, zw, zw6 ! local scalars
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zw  = pw10(ji,jj)
            zw6 = zw*zw*zw
            zw6 = zw6*zw6
            !
            ! When wind speed > 33 m/s => Cyclone conditions => special treatment
            zgt33 = 0.5 + SIGN( 0.5_wp, (zw - 33.) )   ! If pw10 < 33. => 0, else => 1
            !
            cd_neutral_10m(ji,jj) = 1.e-3 * ( &
               &       (1. - zgt33)*( 2.7/zw + 0.142 + zw/13.09 - 3.14807E-10*zw6) & ! wind <  33 m/s
               &      +    zgt33   *      2.34 )                                     ! wind >= 33 m/s
            !
            cd_neutral_10m(ji,jj) = MAX(cd_neutral_10m(ji,jj), 1e-6_wp)
            !
         END DO
      END DO
      !
   END FUNCTION cd_neutral_10m


END MODULE mod_blk_neutral_10m
