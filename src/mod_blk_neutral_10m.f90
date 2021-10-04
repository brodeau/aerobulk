! AeroBulk / 2015 / L. Brodeau

MODULE mod_blk_neutral_10m

   !!====================================================================================
   !!
   !!   Momentum, Latent and sensible heat exchange coefficients
   !!          In neutral conditions / at 10m
   !!
   !!            Author: Laurent Brodeau, 2016
   !!
   !!====================================================================================

   USE mod_const   !: physical and othe constants
   USE mod_phymbl,       ONLY: z0_from_Cd
   USE mod_blk_ncar,     ONLY: CD_N10_NCAR, CH_N10_NCAR, CE_N10_NCAR
   USE mod_blk_coare3p0, ONLY: charn_coare3p0
   USE mod_blk_coare3p6, ONLY: charn_coare3p6
   USE mod_blk_ecmwf,    ONLY: charn0_ecmwf
   USE mod_blk_andreas,  ONLY: u_star_andreas

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: TURB_NEUTRAL_10M

   REAL(wp), PARAMETER :: zu  = 10._wp     ! we're at 10m !


CONTAINS

   SUBROUTINE turb_neutral_10m( calgo, U_N10, CdN10, ChN10, CeN10, pz0 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_neutral_10m  ***
      !!
      !!            2015: L. Brodeau
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


      CHARACTER(len=*),         INTENT(in)    ::   calgo
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   U_N10       ! relative wind module at zu            [m/s]
      REAL(wp), DIMENSION(:,:), INTENT(  out) ::   CdN10       ! transfer coefficient for momentum         (tau)
      REAL(wp), DIMENSION(:,:), INTENT(  out) ::   ChN10       ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), DIMENSION(:,:), INTENT(  out) ::   CeN10       ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), DIMENSION(:,:), INTENT(  out) ::   pz0          ! roughness length [m]

      INTEGER :: Ni, Nj, jit

      REAL(wp), DIMENSION(:,:), ALLOCATABLE  :: u_star, z0t, z0q, Ub, ztmp0, ztmp1
      !!===================================================================================

      Ni = SIZE(U_N10,1)
      Nj = SIZE(U_N10,2)

      ALLOCATE ( u_star(Ni,Nj),   z0t(Ni,Nj), z0q(Ni,Nj), Ub(Ni,Nj),  &
         &        ztmp0(Ni,Nj), ztmp1(Ni,Nj) )

      Ub = MAX(U_N10, 0.1_wp)

      IF ( (TRIM(calgo) == 'coare3p0').OR.(TRIM(calgo) == 'coare3p6').OR.(TRIM(calgo) == 'ecmwf') ) THEN

         !! First guess of CdN10
         CdN10 = 8.575E-5*Ub + 0.657E-3  ! from my curves

         !! ITERATION BLOCK
         DO jit = 1, nb_iter

            !! Need to know u*
            !! --------------
            u_star = Ub*SQRT(CdN10)


            !! Charnock Parameter
            !! ------------------
            IF (     TRIM(calgo) == 'coare3p6' ) THEN
               ztmp0 = charn_coare3p6(Ub)

            ELSEIF ( TRIM(calgo) == 'coare3p0' ) THEN
               ztmp0 = charn_coare3p0(Ub)

            ELSEIF ( TRIM(calgo) == 'ecmwf' ) THEN
               ztmp0 = charn0_ecmwf

            ELSE
               PRINT *, TRIM(calgo)//' => Unknow algo !' ; STOP
            END IF

            !! Roughness lengthes z0, momentum: z0t (z0q = z0t) :
            pz0   = ztmp0*u_star*u_star/grav + 0.11*rnu0_air/u_star ! Roughness length (eq.6)


            !! Compute drag coefficient at 10m
            !! -------------------------------

            ztmp0 = LOG(zu/pz0)

            CdN10 = vkarmn2 / ( ztmp0*ztmp0 )

            !IF ( (Ni==1).AND.(Nj==1) ) PRINT *, '   *** CdN10 =', CdN10*1000.
            IF ( (Ni==1).AND.(Nj==1) ) PRINT *, ''

         END DO


         IF ( TRIM(calgo) == 'coare3p0' ) THEN

            ! Re_r: roughness Reynolds number:
            ztmp1 = pz0*u_star/rnu0_air

            !! Scalar roughness length z0t:
            z0t = MIN( 1.1E-4_wp , 5.5E-5_wp*ztmp1**(-0.6_wp) )     ! Scalar roughness for Theta and q (Fairall al 2003, eq.28)
            z0q = z0t

         END IF

         IF ( TRIM(calgo) == 'coare3p6' ) THEN

            ! Re_r: roughness Reynolds number:
            ztmp1 = pz0*u_star/rnu0_air

            !! Scalar roughness length z0t:
            !! Chris Fairall, Jim Edscon, private communication, March 2016 / COARE 3.5 :
            z0t   = MIN( 1.6e-4_wp , 5.8E-5_wp*ztmp1**(-0.72_wp) ) ! These thermal roughness lengths give Stanton and
            !!                                            ! Dalton numbers that closely approximate COARE3.0
            z0q = z0t

         END IF


         IF ( TRIM(calgo) == 'ecmwf' ) THEN
            ztmp1  = rnu0_air/u_star
            z0t    = 0.40*ztmp1                              ! eq.3.26, Chap.3, p.34, IFS doc - Cy31r1
            z0q    = 0.62*ztmp1
         END IF




         PRINT *, ' *** Algo = ',trim(calgo)

         ztmp1 = LOG(zu/z0t)
         ChN10 = vkarmn2 / ( ztmp0*ztmp1 )     ! (ztmp0 = LOG(zu/z0))

         ztmp1 = LOG(zu/z0q)
         CeN10 = vkarmn2 / ( ztmp0*ztmp1 )     ! (ztmp0 = LOG(zu/z0))



      ELSEIF ( TRIM(calgo) == 'ncar' ) THEN

         PRINT *, ' *** Algo = ',TRIM(calgo)

         Ub = MAX(U_N10, 0.5_wp)

         CdN10  = CD_N10_NCAR( Ub )

         ztmp0 = SQRT(CdN10)

         ChN10 = CH_N10_NCAR( ztmp0 , Ub*0. ) ! 0 => UNSTABLE CASE !!!      L&Y 2004 eq. (6c-6d)
         CeN10 = CE_N10_NCAR( ztmp0 )

         pz0    = MAX( z0_from_Cd( 10._wp , CdN10 ) , 0.0001_wp )
         pz0    = MIN( pz0, 0.1_wp )


      ELSEIF ( TRIM(calgo) == 'andreas' ) THEN

         PRINT *, ' *** Algo = ',TRIM(calgo)

         Ub = MAX(U_N10, 0.5_wp)
         ztmp0 = u_star_andreas( Ub ) ; ! => u*

         PRINT *, ' ===> YET TO BE CODED!!!!'
         STOP

      ELSE

         PRINT *, 'ERROR: algorithm '//TRIM(calgo)//' is not supported yet!'
         PRINT *, ''
         STOP

      END IF

      DEALLOCATE ( u_star, z0t, z0q, Ub, ztmp0, ztmp1 )

   END SUBROUTINE turb_neutral_10m

END MODULE mod_blk_neutral_10m
