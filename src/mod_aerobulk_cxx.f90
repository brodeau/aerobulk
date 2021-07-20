MODULE MOD_AEROBLUK_CXX

   USE mod_aerobulk
   USE iso_c_binding
   USE mod_const

   IMPLICIT NONE
   
   PRIVATE

   !PUBLIC :: l_vap_cxx
   PUBLIC :: aerobulk_cxx_skin, aerobulk_cxx_no_skin

   
CONTAINS

   !==== C++ interface for l_vap ====
   !SUBROUTINE l_vap_cxx( sst, m, l_vap_out ) BIND(c)
   !   USE mod_phymbl
   !   ! Arguments
   !   INTEGER(c_int),                    INTENT(in)  :: m
   !   REAL(c_double), DIMENSION(m,1),    INTENT(in)  :: sst
   !   REAL(c_double), DIMENSION(m,1),    INTENT(out) :: l_vap_out
   !   l_vap_out = L_vap(sst)      
   !END SUBROUTINE l_vap_cxx

   
   !==== C++ interface for using optional skin temperature shceme ====
   SUBROUTINE aerobulk_cxx_skin( jt, Nt, calgo, zt, zu, sst, t_zt, &
      &                          hum_zt, U_zu, V_zu, slp,          &
      &                          QL, QH, Tau_x, Tau_y, Evap,       &
      &                          Niter, l_skin, rad_sw, rad_lw, T_s,       &
      &                          l, m ) BIND(c)
      
      ! Arguments
      INTEGER(c_int),                    INTENT(in)  :: jt, Nt
      INTEGER(c_int),                    INTENT(in)  :: l, m
      CHARACTER(c_char), DIMENSION(l+1), INTENT(in)  :: calgo
      REAL(c_double),                    INTENT(in)  :: zt, zu
      REAL(c_double),    DIMENSION(m,1), INTENT(in)  :: sst, t_zt, hum_zt, U_zu, V_zu, slp
      REAL(c_double),    DIMENSION(m,1), INTENT(out) :: QL, QH, Tau_x, Tau_y, Evap, T_s
      INTEGER(c_int),                    INTENT(in)  :: Niter
      LOGICAL,                           INTENT(in)  :: l_skin
      REAL(c_double),    DIMENSION(m,1), INTENT(in)  :: rad_sw, rad_lw
      
      ! Locals
      CHARACTER(len=l) :: calgo_fort
      INTEGER :: i
      
      ! Loop over the input string to convert to Fortran string
      DO i=1,l
         calgo_fort(i:i) = calgo(i)
      ENDDO
      
      ! Call the actual routine
      ! We could/should transpose the arrays, but it's not neccesary
      CALL AEROBULK_MODEL( jt, Nt, calgo_fort, zt, zu, sst, t_zt,             &
         &                 hum_zt, U_zu, V_zu, slp,                           &
         &                 QL, QH, Tau_x, Tau_y, Evap,                        &
         &                 Niter=Niter, l_use_skin=l_skin, rad_sw=rad_sw, rad_lw=rad_lw, T_s=T_s )
      
   END SUBROUTINE aerobulk_cxx_skin

   
   !==== C++ interface without using optional skin temperature scheme ====
   SUBROUTINE aerobulk_cxx_no_skin( jt, Nt, calgo, zt, zu, sst, t_zt, &
      &                             hum_zt, U_zu, V_zu, slp,          &
      &                             QL, QH, Tau_x, Tau_y, Evap,       &
      &                             Niter, l, m ) BIND(c)
      
      ! Arguments
      INTEGER(c_int),                    INTENT(in)  :: jt, Nt
      INTEGER(c_int),                    INTENT(in)  :: l, m
      CHARACTER(c_char), DIMENSION(l+1), INTENT(in)  :: calgo
      REAL(c_double),                    INTENT(in)  :: zt, zu
      REAL(c_double),    DIMENSION(m,1), INTENT(in)  :: sst, t_zt, hum_zt, U_zu, V_zu, slp
      REAL(c_double),    DIMENSION(m,1), INTENT(out) :: QL, QH, Tau_x, Tau_y, Evap
      INTEGER(c_int),                    INTENT(in)  :: Niter
      
      ! Locals
      CHARACTER(len=l) :: calgo_fort
      INTEGER :: i
      
      ! Loop over the input string to convert to Fortran string
      DO i=1,l
         calgo_fort(i:i) = calgo(i)
      ENDDO
      
      ! Call the actual routine
      ! We could/should transpose the arrays, but it's not neccesary
      CALL AEROBULK_MODEL( jt, Nt, calgo_fort, zt, zu, sst, t_zt,  &
         &                 hum_zt, U_zu, V_zu, slp,                &
         &                 QL, QH, Tau_x, Tau_y, Evap, Niter=Niter )
      
   END SUBROUTINE aerobulk_cxx_no_skin
   
END MODULE MOD_AEROBLUK_CXX
