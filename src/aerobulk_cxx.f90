      integer function aerobulk_cxx_all( calgo, zt, zu, sst, t_zt, &
     &                       q_zt, U_zu, V_zu, slp,    &
     &                       QL, QH, Tau_x, Tau_y,     &
     &                       Niter, rad_sw, rad_lw,    &
     &                       l, m ) bind(c)

      use iso_c_binding
      use mod_aerobulk

      implicit none

      ! Arguments
      INTEGER(c_int),                    INTENT(in)  :: l, m, Niter
      CHARACTER(c_char), DIMENSION(l+1), INTENT(in)  :: calgo
      REAL(c_double),                    INTENT(in)  :: zt, zu
      REAL(c_double), DIMENSION(m,1),    INTENT(in)  :: sst, t_zt, q_zt, U_zu, V_zu, slp
      REAL(c_double), DIMENSION(m,1),    INTENT(out) :: QL, QH, Tau_x, Tau_y
      REAL(c_double), DIMENSION(m,1),    INTENT(in)  :: rad_sw, rad_lw

      ! Locals
      character(len=l) :: calgo_fort
      integer :: i

      ! Loop over the input string to convert to Fortran string
      do i=1,l
        calgo_fort(i:i) = calgo(i)
      enddo

      ! Call the actual routine
      ! We could/should transpose the arrays, but it's not neccesary
      CALL aerobulk_model( calgo_fort, zt, zu, sst, t_zt, &
     &                 q_zt, U_zu, V_zu, slp,   &
     &                 QL, QH, Tau_x, Tau_y,    &
     &                 Niter=Niter, rad_sw=rad_sw, rad_lw=rad_lw )

      end function

      integer function aerobulk_cxx_rad( calgo, zt, zu, sst, t_zt, &
     &                       q_zt, U_zu, V_zu, slp,    &
     &                       QL, QH, Tau_x, Tau_y,     &
     &                       rad_sw, rad_lw,           &
     &                       l, m ) bind(c)

      use iso_c_binding
      use mod_aerobulk

      implicit none

      ! Arguments
      INTEGER(c_int),                    INTENT(in)  :: l, m
      CHARACTER(c_char), DIMENSION(l+1), INTENT(in)  :: calgo
      REAL(c_double),                    INTENT(in)  :: zt, zu
      REAL(c_double), DIMENSION(m,1),    INTENT(in)  :: sst, t_zt, q_zt, U_zu, V_zu, slp
      REAL(c_double), DIMENSION(m,1),    INTENT(out) :: QL, QH, Tau_x, Tau_y
      REAL(c_double), DIMENSION(m,1),    INTENT(in)  :: rad_sw, rad_lw

      ! Locals
      character(len=l) :: calgo_fort
      integer :: i

      ! Loop over the input string to convert to Fortran string
      do i=1,l
        calgo_fort(i:i) = calgo(i)
      enddo

      ! Call the actual routine
      ! We could/should transpose the arrays, but it's not neccesary
      CALL aerobulk_model( calgo_fort, zt, zu, sst, t_zt, &
     &                 q_zt, U_zu, V_zu, slp,   &
     &                 QL, QH, Tau_x, Tau_y,    &
     &                 rad_sw=rad_sw, rad_lw=rad_lw )

      end function

      integer function aerobulk_cxx_niter( calgo, zt, zu, sst, t_zt, &
     &                       q_zt, U_zu, V_zu, slp,    &
     &                       QL, QH, Tau_x, Tau_y,     &
     &                       Niter,                    &
     &                       l, m ) bind(c)

      use iso_c_binding
      use mod_aerobulk

      implicit none

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

      ! Call the actual routine
      ! We could/should transpose the arrays, but it's not neccesary
      CALL aerobulk_model( calgo_fort, zt, zu, sst, t_zt, &
     &                 q_zt, U_zu, V_zu, slp,   &
     &                 QL, QH, Tau_x, Tau_y,    &
     &                 Niter=Niter)

      end function

      integer function aerobulk_cxx_noopt( calgo, zt, zu, sst, t_zt, &
     &                       q_zt, U_zu, V_zu, slp,    &
     &                       QL, QH, Tau_x, Tau_y,     &
     &                       l, m ) bind(c)

      use iso_c_binding
      use mod_aerobulk

      implicit none

      ! Arguments
      INTEGER(c_int),                    INTENT(in)  :: l, m
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

      ! Call the actual routine
      ! We could/should transpose the arrays, but it's not neccesary
      CALL aerobulk_model( calgo_fort, zt, zu, sst, t_zt, &
     &                 q_zt, U_zu, V_zu, slp,   &
     &                 QL, QH, Tau_x, Tau_y )

      end function

