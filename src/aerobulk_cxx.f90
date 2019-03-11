module mod_aerobluk_cxx

  use mod_aerobulk
  use iso_c_binding
  use mod_const

  implicit none

  private

  public :: lvap_cxx, aerobulk_cxx_skin, aerobulk_cxx_no_skin

contains

!==== C++ interface for lvap ====
      subroutine lvap_cxx( sst, m, lvap_out ) bind(c)

      use mod_thermo

      ! Arguments
      INTEGER(c_int),                    INTENT(in)  :: m
      REAL(c_double), DIMENSION(m,1),    INTENT(in)  :: sst
      REAL(c_double), DIMENSION(m,1),    INTENT(out) :: lvap_out

      lvap_out = Lvap(sst)

      end subroutine lvap_cxx

!==== C++ interface for using optional skin temperature shceme ====
      subroutine aerobulk_cxx_skin( calgo, zt, zu, sst, t_zt, &
                                   q_zt, U_zu, V_zu, slp,    &
                                   QL, QH, Tau_x, Tau_y,     &
                                   rad_sw, rad_lw, T_s,      &
                                   Niter, l, m ) bind(c)

      ! Arguments
      INTEGER(c_int),                    INTENT(in)  :: l, m, Niter
      CHARACTER(c_char), DIMENSION(l+1), INTENT(in)  :: calgo
      REAL(c_double),                    INTENT(in)  :: zt, zu
      REAL(c_double), DIMENSION(m,1),    INTENT(in)  :: sst, t_zt, q_zt, U_zu, V_zu, slp
      REAL(c_double), DIMENSION(m,1),    INTENT(out) :: QL, QH, Tau_x, Tau_y, T_s
      REAL(c_double), DIMENSION(m,1),    INTENT(in)  :: rad_sw, rad_lw

      ! Locals
      character(len=l) :: calgo_fort
      integer :: i

      ! Loop over the input string to convert to Fortran string
      do i=1,l
        calgo_fort(i:i) = calgo(i)
      enddo

      ! Check that the jpi and jpj are right
      if ( jpi /= m .or. jpj /= 1 ) then
        l_first_call = .true.
      end if

      ! Call the actual routine
      ! We could/should transpose the arrays, but it's not neccesary
      CALL aerobulk_model( calgo_fort, zt, zu, sst, t_zt, &
                           q_zt, U_zu, V_zu, slp,         &
                           QL, QH, Tau_x, Tau_y,          &
                           Niter=Niter, rad_sw=rad_sw, rad_lw=rad_lw, T_s=T_s )

      end subroutine aerobulk_cxx_skin

!==== C++ interface without using optional skin temperature scheme ====
      subroutine aerobulk_cxx_no_skin( calgo, zt, zu, sst, t_zt, &
                                     q_zt, U_zu, V_zu, slp,    &
                                     QL, QH, Tau_x, Tau_y,     &
                                     Niter, l, m ) bind(c)

      ! Arguments
      INTEGER(c_int),                    INTENT(in)  :: l, m, Niter
      CHARACTER(c_char), DIMENSION(l+1), INTENT(in)  :: calgo
      REAL(c_double),                    INTENT(in)  :: zt, zu
      REAL(c_double), DIMENSION(m,1),    INTENT(in)  :: sst, t_zt, q_zt, U_zu, V_zu, slp
      REAL(c_double), DIMENSION(m,1),    INTENT(out) :: QL, QH, Tau_x, Tau_y

      ! Locals
      character(len=l) :: calgo_fort
      integer :: i

      ! Loop over the input string to convert to Fortran string
      do i=1,l
        calgo_fort(i:i) = calgo(i)
      enddo

      ! Check that the jpi and jpj are right
      if ( jpi /= m .or. jpj /= 1 ) then
        l_first_call = .true.
      end if

      ! Call the actual routine
      ! We could/should transpose the arrays, but it's not neccesary
      CALL aerobulk_model( calgo_fort, zt, zu, sst, t_zt, &
                           q_zt, U_zu, V_zu, slp,   &
                           QL, QH, Tau_x, Tau_y, Niter=Niter )

      end subroutine aerobulk_cxx_no_skin

end module mod_aerobluk_cxx

