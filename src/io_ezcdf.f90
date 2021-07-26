MODULE io_ezcdf

   USE netcdf

   !! Netcdf input/output
   !!
   !! Author: Laurent Brodeau, 2010
   !!

   IMPLICIT NONE

   INTEGER, PARAMETER, PUBLIC :: nbatt_max=20

   TYPE, PUBLIC :: var_attr
      CHARACTER(LEN=128) :: cname
      INTEGER            :: itype
      INTEGER            :: ilength
      CHARACTER(LEN=128) :: val_char
      REAL, DIMENSION(9) :: val_num   ! assuming that numeric attributes are never a vector larger than 9...
   END TYPE var_attr


   TYPE, PUBLIC :: t_unit_t0
      CHARACTER(LEN=1)   :: unit
      INTEGER            :: year
      INTEGER            :: month
      INTEGER            :: day
      INTEGER            :: hour
      INTEGER            :: minute
      INTEGER            :: second
   END TYPE t_unit_t0

   TYPE, PUBLIC :: date
      INTEGER            :: year
      INTEGER            :: month
      INTEGER            :: day
      INTEGER            :: hour
      INTEGER            :: minute
      INTEGER            :: second
   END TYPE date



   PRIVATE


   INTERFACE GETVAR_1D
      MODULE PROCEDURE GETVAR_1D_R8, GETVAR_1D_R4, GETVAR_1D_R8_1x1
   END INTERFACE GETVAR_1D

   INTERFACE GETVAR_2D
      MODULE PROCEDURE GETVAR_2D_R8, GETVAR_2D_R4
   END INTERFACE GETVAR_2D

   INTERFACE GETVAR_3D
      MODULE PROCEDURE GETVAR_3D_R8, GETVAR_3D_R4
   END INTERFACE GETVAR_3D

   INTERFACE DUMP_FIELD
      MODULE PROCEDURE DUMP_2D_FIELD, DUMP_3D_FIELD
   END INTERFACE DUMP_FIELD



   !! List of public routines
   !! =======================
   PUBLIC :: dims,      &
      &    get_sf_ao,        &
      &    getvar_1d,        &
      &    getvar_2d,        &
      &    getvar_3d,        &
      &    getvar_attributes,&
      &    force_attr,       &
      &    getmask_2d,       & ! not time records (at least no more than 1!)
      &    getmask_3d,       & !        "              "
      &    pt_series,        &
      &    p2d_t,            &
      &    p3d_t,            &
      &    check_4_miss,     &
      &    get_var_info,     &
      &    dump_field,       &
      &    p2d_mapping_ab,   &
      &    rd_mapping_ab,    &
      &    phovmoller,       &
      &    get_time_unit_t0, &
      &    l_is_leap_year,   &
      &    test_xyz,         &
      &    time_to_date,     &
      &    to_epoch_time_scalar, to_epoch_time_vect, &
      &    coordinates_from_var_attr, &
      &    GETVAR_1D_R8_3x3_to_1x1
   !!===========================


   CHARACTER(len=80) :: cv_misc

   CHARACTER(len=2) :: cdt  ! '1d' or '2d'

   REAL(8), DIMENSION(3,2) :: vextrema

   INTEGER :: nd

   CHARACTER(len=400)    :: cu

   CHARACTER(len=8), PARAMETER :: cdum = 'dummy'

   CHARACTER(LEN=400), PARAMETER   ::     &
      &    cabout = 'Created with SOSIE interpolation environement => https://github.com/brodeau/sosie/'

   INTEGER :: ji, jj, jk

   !! About missing value attribute name:
   INTEGER, PARAMETER :: nmval = 4
   CHARACTER(len=14), DIMENSION(nmval), PARAMETER :: &
      &     c_nm_miss_val = (/  '_FillValue    ', 'missing_value ', 'FillValue     ', '_Fillvalue    ' /)
   CHARACTER(LEN=14), PARAMETER :: cmv0 = '_FillValue' ! Default name for Missing/Fill value ...

   INTEGER, DIMENSION(12), PARAMETER :: &
      &   tdmn = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /),        &
      &   tdml = (/ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /),        &
      &  tcdmn = (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /), &
      &  tcdml = (/ 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 /)


   REAL(8), SAVE :: rjs_t_old
   
CONTAINS



   SUBROUTINE DIMS(cf_in, cv_in, ni, nj, lz, Nt)

      !!-----------------------------------------------------------------------
      !!
      !! Updated August 2012, L. Brodeau
      !!
      !! This routine opens a netcdf file 'cf_in' to check the dimension
      !! of the variable 'cv_in'. It then gives the length of each of the dimension,
      !! if the length returns '-1' that means that the dimension does not exist
      !!
      !! example : if the variable has only 1 dimension, of length 132,
      !!           DIMS will return ni=132, nj=-1, lz=-1, Nt=-1
      !!
      !!
      !! INPUT :
      !! -------
      !!          * cf_in       : name of the input file          (character)
      !!          * cv_in       : name of the variable            (character)
      !!
      !! OUTPUT :
      !! --------
      !!          * ni      : first dimension                     (integer)
      !!          * nj      : second dimension  (-1 if none)      (integer)
      !!          * lz      : third dimension   (-1 if none)      (integer)
      !!          * Nt      : number of records (-1 if none)      (integer)
      !!
      !!------------------------------------------------------------------------
      INTEGER                         :: id_f, id_v
      CHARACTER(len=*),   INTENT(in)  :: cf_in, cv_in
      INTEGER,            INTENT(out) :: ni, nj, lz, Nt

      INTEGER, DIMENSION(:), ALLOCATABLE :: id_dim, nlen
      INTEGER :: jdim, id_unlim_dim
      CHARACTER(len=80), PARAMETER :: crtn = 'DIMS'

      ni = -1 ; nj = -1 ; lz = -1 ; Nt = -1

      CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, id_f),  crtn,cf_in,cv_in)

      ! Get ID of unlimited dimension
      CALL sherr(  NF90_INQUIRE(id_f, unlimitedDimId = id_unlim_dim), crtn,cf_in,cv_in)



      !! Getting variable ID:
      CALL sherr( NF90_INQ_VARID(id_f,cv_in,id_v),              crtn,cf_in,cv_in)

      !! nd => number of dimensions for the variable:
      CALL sherr( NF90_INQUIRE_VARIABLE(id_f, id_v, ndims=nd),  crtn,cf_in,cv_in)

      ALLOCATE ( id_dim(nd) , nlen(nd) )

      !! Vector containing the IDs of each dimension (id_dim):
      CALL sherr( NF90_INQUIRE_VARIABLE(id_f, id_v, dimids=id_dim),  crtn,cf_in,cv_in)

      DO jdim = 1, nd
         CALL sherr( NF90_INQUIRE_DIMENSION(id_f, id_dim(jdim), len=nlen(jdim)),   crtn,cf_in,cv_in)
      END DO

      IF ( (nd > 2).AND.(id_unlim_dim < 1) ) THEN
         WRITE(6,*) 'WARNING: DIMS of io_ezcdf.f90'
         WRITE(6,*) '   => variable '//TRIM(cv_in)//' in file:'
         WRITE(6,*) '      '//TRIM(cf_in)
         WRITE(6,*) '      does not have an UNLIMITED dimension !!!'
         WRITE(6,*) '   => if it is supposed to depend on a time-record, this time-record should be'
         WRITE(6,*) '      an unlimited DIMENSION into the netcdf file!!!'
         WRITE(6,*) '   => otherwize time dimension might be confused with a space dimension!!!'
         WRITE(6,*) '   => maybe SOSIE can overcome this, but at your own risks... ;)'
         WRITE(6,*) ''
      END IF

      SELECT CASE(nd)

      CASE(1)
         !! Purenj 1D
         ni = nlen(1)

      CASE(2)
         IF ( id_dim(2) == id_unlim_dim ) THEN
            !! 1D with time records
            ni = nlen(1) ; Nt = nlen(2)
         ELSE
            !! 2D with no time records
            ni = nlen(1) ; nj = nlen(2)
         END IF

      CASE(3)
         IF ( id_dim(3) == id_unlim_dim ) THEN
            !! 2D with time records
            ni = nlen(1) ; nj = nlen(2) ; Nt = nlen(3)
         ELSE
            !! 3D with no time records
            ni = nlen(1) ; nj = nlen(2) ; lz = nlen(3)
         END IF

      CASE(4)
         IF ( id_unlim_dim < 1 ) THEN
            WRITE(6,*) 'ERROR: file ',trim(cf_in),' doesnt have an unlimited dimension (time record)!'
         END IF

         ni = nlen(1) ; nj = nlen(2)
         IF ( id_dim(3) == id_unlim_dim ) THEN
            lz = nlen(4) ; Nt = nlen(3)   ! time record (unlimited dim) comes as 3rd dim and lz as 4th
         ELSE
            lz = nlen(3) ; Nt = nlen(4)   ! time record (unlimited dim) comes as last dim and lz as 3rd
         END IF

      CASE DEFAULT
         CALL print_err(crtn, 'the dimension is not realistic')

      END SELECT

      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_in,cv_in)

   END SUBROUTINE DIMS





   SUBROUTINE GETVAR_ATTRIBUTES(cf_in, cv_in,  Nb_att, v_att_list)
      !!-----------------------------------------------------------------------
      !! This routine gets all the attributes from a given variable (cv_in) from a given file (cf_in)
      !!
      !! INPUT :
      !! -------
      !!          * cf_in      : name of the input file             (character l=100)
      !!          * cv_in      : name of the variable               (character l=20)
      !!
      !! OUTPUT :
      !! --------
      !!          *
      !!
      !!------------------------------------------------------------------------
      CHARACTER(len=*),                   INTENT(in)  :: cf_in, cv_in
      INTEGER,                            INTENT(out) :: Nb_att
      TYPE(var_attr), DIMENSION(nbatt_max), INTENT(out) :: v_att_list
      !!
      !Local:
      INTEGER :: id_f, id_v
      CHARACTER(len=256) :: cname, cvalue
      REAL, DIMENSION(:), ALLOCATABLE :: rvalue ! will store attribute with numeric values
      INTEGER :: ierr, jatt, iwhat, ilg
      !!
      CHARACTER(len=80), PARAMETER :: crtn = 'GETVAR_ATTRIBUTES'

      CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, id_f),     crtn,cf_in,cv_in)
      CALL sherr( NF90_INQ_VARID(id_f, trim(cv_in), id_v),  crtn,cf_in,cv_in)

      v_att_list(:)%cname    = 'null'
      v_att_list(:)%val_char = 'null'

      DO jatt = 1, nbatt_max
         ierr = NF90_INQ_ATTNAME(id_f, id_v, jatt, cname)
         IF ( ierr == 0 ) THEN
            v_att_list(jatt)%cname = cname

            CALL sherr( NF90_INQUIRE_ATTRIBUTE(id_f, id_v, TRIM(cname), xtype=iwhat, len=ilg),  crtn,cf_in,cv_in)
            !PRINT *, 'LOLO: attribute = '//TRIM(cname)
            !PRINT *, '   iwhat = ', iwhat

            v_att_list(jatt)%itype   = iwhat
            !v_att_list(jatt)%ctype   = vtypes_def(iwhat)
            v_att_list(jatt)%ilength = ilg
            !! Getting value of attribute, depending on type!
            IF ( iwhat == 2 ) THEN
               CALL sherr( NF90_GET_ATT(id_f, id_v, cname, cvalue),  crtn,cf_in,cv_in)
               v_att_list(jatt)%val_char = TRIM(cvalue)
            ELSE
               ALLOCATE ( rvalue(ilg) )
               CALL sherr( NF90_GET_ATT(id_f, id_v, cname, rvalue),  crtn,cf_in,cv_in)
               v_att_list(jatt)%val_num(1:ilg) = rvalue(:)
               DEALLOCATE ( rvalue )
            END IF
         ELSE
            EXIT
         END IF
      END DO
      Nb_att = jatt-1
      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_in,cv_in)
   END SUBROUTINE GETVAR_ATTRIBUTES





   SUBROUTINE FORCE_ATTR(cattr, cval, v_att_list)
      !!------------------------------------------------------------------------
      !!------------------------------------------------------------------------
      CHARACTER(len=*),                     INTENT(in)    :: cattr, cval
      TYPE(var_attr), DIMENSION(nbatt_max), INTENT(inout) :: v_att_list
      !Local:
      INTEGER :: jatt
      !!
      CHARACTER(len=80), PARAMETER :: crtn = 'FORCE_ATTR'
      !!
      !! Find position of attribute to modify "cattr" and change its content if found!
      DO jatt = 0, nbatt_max-1
         IF ( TRIM(v_att_list(jatt+1)%cname) == 'null' ) EXIT
         IF ( TRIM(v_att_list(jatt+1)%cname) == TRIM(cattr) ) THEN
            v_att_list(jatt+1)%val_char = TRIM(cval) ! Setting value
            EXIT
         END IF
      END DO
   END SUBROUTINE FORCE_ATTR













   SUBROUTINE GETVAR_1D_R8(cf_in, cv_in, VX)
      !!-----------------------------------------------------------------------
      !! This routine extract a variable 1D from a netcdf file
      !!
      !! INPUT :
      !! -------
      !!          * cf_in      : name of the input file             (character l=100)
      !!          * cv_in      : name of the variable               (character l=20)
      !!
      !! OUTPUT :
      !! --------
      !!          * VX         : 1D array contening the variable   (double)
      !!
      !!------------------------------------------------------------------------
      INTEGER                             :: id_f, id_v
      CHARACTER(len=*),       INTENT(in)  :: cf_in, cv_in
      REAL(8), DIMENSION (:), INTENT(out) :: VX
      INTEGER :: ierr1, ierr2
      REAL(4) :: rsf, rao
      CHARACTER(len=80), PARAMETER :: crtn = 'GETVAR_1D_R8'
      CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, id_f),     crtn,cf_in,cv_in)
      CALL sherr( NF90_INQ_VARID(id_f, TRIM(cv_in), id_v),  crtn,cf_in,cv_in)
      ierr1 = NF90_GET_ATT(id_f, id_v, 'scale_factor', rsf)
      ierr2 = NF90_GET_ATT(id_f, id_v, 'add_offset',   rao)
      CALL sherr( NF90_GET_VAR(id_f, id_v, VX),              crtn,cf_in,cv_in)
      IF (ierr1 == NF90_NOERR) VX = rsf*VX
      IF (ierr2 == NF90_NOERR) VX = VX + rao
      CALL sherr( NF90_CLOSE(id_f),                         crtn,cf_in,cv_in)
   END SUBROUTINE GETVAR_1D_R8


   SUBROUTINE GETVAR_1D_R4(cf_in, cv_in, VX)
      !!-----------------------------------------------------------------------
      !! This routine extract a variable 1D from a netcdf file
      !!
      !! INPUT :
      !! -------
      !!          * cf_in      : name of the input file             (character l=100)
      !!          * cv_in      : name of the variable               (character l=20)
      !!
      !! OUTPUT :
      !! --------
      !!          * VX         : 1D array contening the variable   (double)
      !!
      !!------------------------------------------------------------------------
      INTEGER                             :: id_f, id_v
      CHARACTER(len=*),       INTENT(in)  :: cf_in, cv_in
      REAL(4), DIMENSION (:), INTENT(out) :: VX
      INTEGER :: ierr1, ierr2
      REAL(4) :: rsf, rao
      CHARACTER(len=80), PARAMETER :: crtn = 'GETVAR_1D_R4'
      CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, id_f),     crtn,cf_in,cv_in)
      CALL sherr( NF90_INQ_VARID(id_f, TRIM(cv_in), id_v),  crtn,cf_in,cv_in)
      ierr1 = NF90_GET_ATT(id_f, id_v, 'scale_factor', rsf)
      ierr2 = NF90_GET_ATT(id_f, id_v, 'add_offset',   rao)
      CALL sherr( NF90_GET_VAR(id_f, id_v, VX),              crtn,cf_in,cv_in)
      IF (ierr1 == NF90_NOERR) VX = rsf*VX
      IF (ierr2 == NF90_NOERR) VX = VX + rao
      CALL sherr( NF90_CLOSE(id_f),                         crtn,cf_in,cv_in)
   END SUBROUTINE GETVAR_1D_R4

   SUBROUTINE GETVAR_1D_R8_1x1(cf_in, cv_in, X)
      !!-----------------------------------------------------------------------
      !! This routine extract a variable 1D from a netcdf file
      !!
      !! INPUT :
      !! -------
      !!          * cf_in      : name of the input file             (character l=100)
      !!          * cv_in      : name of the variable               (character l=20)
      !!
      !! OUTPUT :
      !! --------
      !!          * X         : fake 3D array contening a 1D vector (field of shape 1 x 1 x Nx)  (double)
      !!
      !!------------------------------------------------------------------------
      INTEGER                             :: id_f, id_v
      CHARACTER(len=*),       INTENT(in)  :: cf_in, cv_in
      REAL(8), DIMENSION (:,:,:), INTENT(out) :: X
      INTEGER :: n1, n2, nx, ierr1, ierr2
      REAL(4) :: rsf, rao
      CHARACTER(len=80), PARAMETER :: crtn = 'GETVAR_1D_R8_1x1'
      n1 = SIZE(X,1)
      n2 = SIZE(X,2)
      IF ( (n1/=1).OR.(n2/=1) ) CALL print_err(crtn, ' PROBLEM #1 => 1st and 2nd dimmension /= 1 !!!', ivect=(/n1,n2/))
      nx = SIZE(X,3)
      CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, id_f),     crtn,cf_in,cv_in)
      CALL sherr( NF90_INQ_VARID(id_f, TRIM(cv_in), id_v),  crtn,cf_in,cv_in)
      ierr1 = NF90_GET_ATT(id_f, id_v, 'scale_factor', rsf)
      ierr2 = NF90_GET_ATT(id_f, id_v, 'add_offset',   rao)
      CALL sherr( NF90_GET_VAR(id_f, id_v, X, start=(/1,1,1/), count=(/1,1,nx/)), crtn,cf_in,cv_in)
      IF (ierr1 == NF90_NOERR) X = rsf*X
      IF (ierr2 == NF90_NOERR) X = X + rao
      CALL sherr( NF90_CLOSE(id_f),                         crtn,cf_in,cv_in)
   END SUBROUTINE GETVAR_1D_R8_1x1


   SUBROUTINE GETVAR_1D_R8_3x3_to_1x1(cf_in, cv_in, X)
      !!-----------------------------------------------------------------------
      !! This routine extract a variable 1D from a netcdf file
      !!
      !! INPUT :
      !! -------
      !!          * cf_in      : name of the input file             (character l=100)
      !!          * cv_in      : name of the variable               (character l=20)
      !!
      !! OUTPUT :
      !! --------
      !!          * X         : fake 3D array contening a 1D vector (field of shape 1 x 1 x Nx)  (double)
      !!
      !!------------------------------------------------------------------------
      !! The field in netcdf is 2D(3x3) + T, yet here X must be 2D(1x1)
      INTEGER                             :: id_f, id_v
      CHARACTER(len=*),       INTENT(in)  :: cf_in, cv_in
      REAL(8), DIMENSION (:,:,:), INTENT(out) :: X
      INTEGER :: n1, n2, nx, ierr1, ierr2
      REAL(4) :: rsf, rao
      CHARACTER(len=80), PARAMETER :: crtn = 'GETVAR_1D_R8_3x3_to_1x1'
      n1 = SIZE(X,1)
      n2 = SIZE(X,2)
      IF ( (n1/=1).OR.(n2/=1) ) CALL print_err(crtn, ' PROBLEM #1 => 1st and 2nd dimmension /= 3 !!!')
      nx = SIZE(X,3)
      CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, id_f),     crtn,cf_in,cv_in)
      CALL sherr( NF90_INQ_VARID(id_f, TRIM(cv_in), id_v),  crtn,cf_in,cv_in)
      ierr1 = NF90_GET_ATT(id_f, id_v, 'scale_factor', rsf)
      ierr2 = NF90_GET_ATT(id_f, id_v, 'add_offset',   rao)
      CALL sherr( NF90_GET_VAR(id_f, id_v, X, start=(/2,2,1/), count=(/1,1,nx/)), crtn,cf_in,cv_in) !lolo: 2,2,1 for star !!! the middle !!
      IF (ierr1 == NF90_NOERR) X = rsf*X
      IF (ierr2 == NF90_NOERR) X = X + rao
      CALL sherr( NF90_CLOSE(id_f),                         crtn,cf_in,cv_in)
   END SUBROUTINE GETVAR_1D_R8_3x3_to_1x1







   SUBROUTINE GETVAR_2D_R4(idx_f, idx_v, cf_in, cv_in, Nt, kz, jt, X, jt1, jt2, lz)
      !!-----------------------------------------------------------------------------
      !! This routine extract a 2D field from a netcdf file
      !! at a given time
      !!
      !! INPUT :
      !! -------
      !!          * idx_f    : ID of current file                  (integer)
      !!          * idx_v    : ID of current variable              (integer)
      !!          * cf_in      : name of the input file              (character)
      !!          * cv_in      : name of the variable                (character)
      !!          * Nt        : number of time records     (integer)
      !!          * kz        : level to extract                    (integer)
      !!                      0 => input file does not have levels ( 1 would work anyway...)
      !!          * jt        : time snapshot to extract            (integer)
      !!                      0 => input file does not have a time snapshot
      !! OUTPUT :
      !! --------
      !!          * X         : 2D array contening the variable     (real)
      !!
      !! OPTIONAL INPUT :
      !! ----------------
      !!          * jt1, jt2  : first and last time snapshot to extract
      !!          *       lz  : number of levels to know when they are all read
      !!                        so we can close the file
      !!------------------------------------------------------------------------
      INTEGER,                   INTENT(inout) :: idx_f, idx_v
      CHARACTER(len=*),          INTENT(in)    :: cf_in, cv_in
      INTEGER,                   INTENT(in)    :: Nt, kz, jt
      REAL(4),  DIMENSION (:,:), INTENT(out)   :: X
      INTEGER,       OPTIONAL,   INTENT(in)    :: jt1, jt2, lz

      INTEGER :: ni, nj, n1, n2, n3, n4, jlev, its, ite, kz_stop = 0, ierr1, ierr2
      REAL(4) :: rsf, rao
      LOGICAL :: l_okay
      CHARACTER(len=80), PARAMETER :: crtn = 'GETVAR_2D_R4'

      ni = size(X,1)
      nj = size(X,2)

      jlev = kz ! so we can modify jlev without affecting kz...

      CALL DIMS(cf_in, cv_in, n1, n2, n3, n4)
      IF ( (ni /= n1).OR.(nj /= n2) ) CALL print_err(crtn, ' PROBLEM #1 => '//TRIM(cv_in)//' in '//TRIM(cf_in), ivect=(/n1,n2,ni,nj/))

      its = 1 ; ite = Nt
      IF ( PRESENT(jt1) ) its = jt1
      IF ( PRESENT(jt2) ) ite = jt2

      IF ( present(lz) ) kz_stop = lz

      IF ( (jt == its).OR.(jt == 0) ) THEN   ! Opening file and defining variable :
         PRINT *, ''
         PRINT *, ' --- GETVAR_2D: opening file '//TRIM(cf_in)//' for '//TRIM(cv_in)//' !'
         CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, idx_f),  crtn,cf_in,cv_in)
         CALL sherr( NF90_INQ_VARID(idx_f, cv_in, idx_v),  crtn,cf_in,cv_in)
      END IF

      ierr1 = NF90_GET_ATT(idx_f, idx_v, 'scale_factor', rsf) ; !lolo, ugnj at each time...
      ierr2 = NF90_GET_ATT(idx_f, idx_v, 'add_offset',   rao)

      IF ( (idx_f==0).AND.(idx_v==0) ) CALL print_err(crtn, ' PROBLEM #2 file and variable handle never created => '//TRIM(cv_in)//' in '//TRIM(cf_in))

      l_okay = .FALSE.
      DO WHILE ( .NOT. l_okay )

         IF ( jlev == 0 ) THEN    ! No levels
            IF ( jt == 0 ) THEN
               CALL sherr( NF90_GET_VAR(idx_f, idx_v, X, start=(/1,1/), count=(/ni,nj/)), &
                  &      crtn,cf_in,cv_in)
            ELSEIF ( jt > 0 ) THEN
               CALL sherr( NF90_GET_VAR(idx_f, idx_v, X, start=(/1,1,jt/), count=(/ni,nj,1/)), &
                  &      crtn,cf_in,cv_in)
            END IF
            l_okay = .TRUE. ! we can exit the WHILE loop...

         ELSEIF ( jlev > 0 ) THEN

            !! User possibnj specified jlev = 1 and there is not an existing level dimension:
            IF ( n3 == -1 ) THEN
               PRINT *, ' *** warning: ',trim(crtn),' => there is actualnj no levels for ', trim(cv_in),' in ',trim(cf_in)
               PRINT *, '              => fixing it...'
               jlev = 0 ! => should be treated at next "while" loop...
            ELSE
               IF ( jlev >  n3 ) CALL print_err(crtn, ' you want extract a level greater than max value')
               IF ( jt == 0 ) THEN
                  CALL sherr( NF90_GET_VAR(idx_f, idx_v, X, start=(/1,1,jlev/), count=(/ni,nj,1/)), crtn,cf_in,cv_in)
               ELSEIF ( jt > 0 ) THEN
                  CALL sherr( NF90_GET_VAR(idx_f, idx_v, X, start=(/1,1,jlev,jt/), count=(/ni,nj,1,1/)), crtn,cf_in,cv_in)
               END IF
               l_okay = .TRUE.  ! we can exit the WHILE loop...
            END IF
         END IF
      END DO

      IF (ierr1 == NF90_NOERR) X = rsf*X
      IF (ierr2 == NF90_NOERR) X = X + rao

      IF ( ( (jt == ite ).OR.(jt == 0) ).AND.( (jlev == kz_stop).OR.(kz_stop == 0) ) )  THEN
         PRINT *, ' --- GETVAR_2D: closing file '//TRIM(cf_in)//' !'
         CALL sherr( NF90_CLOSE(idx_f),  crtn,cf_in,cv_in)
         idx_f = 0 ; idx_v = 0
         PRINT *, ''
      END IF

   END SUBROUTINE GETVAR_2D_R4

   SUBROUTINE GETVAR_2D_R8(idx_f, idx_v, cf_in, cv_in, Nt, kz, jt, X, jt1, jt2, lz)
      !!-----------------------------------------------------------------------------
      !! See GETVAR_2D_R4
      !!------------------------------------------------------------------------
      INTEGER,                   INTENT(inout) :: idx_f, idx_v
      CHARACTER(len=*),          INTENT(in)    :: cf_in, cv_in
      INTEGER,                   INTENT(in)    :: Nt, kz, jt
      REAL(8),  DIMENSION (:,:), INTENT(out)   :: X
      INTEGER,       OPTIONAL,   INTENT(in)    :: jt1, jt2, lz

      INTEGER :: ni, nj, n1, n2, n3, n4, jlev, its, ite, kz_stop = 0, ierr1, ierr2
      REAL(4) :: rsf, rao
      LOGICAL :: l_okay
      CHARACTER(len=80), PARAMETER :: crtn = 'GETVAR_2D_R8'

      ni = size(X,1)
      nj = size(X,2)

      jlev = kz ! so we can modify jlev without affecting kz...

      CALL DIMS(cf_in, cv_in, n1, n2, n3, n4)
      IF ( (ni /= n1).OR.(nj /= n2) ) CALL print_err(crtn, ' PROBLEM #1 => '//TRIM(cv_in)//' in '//TRIM(cf_in), ivect=(/n1,n2,ni,nj/))

      its = 1 ; ite = Nt
      IF ( PRESENT(jt1) ) its = jt1
      IF ( PRESENT(jt2) ) ite = jt2

      IF ( present(lz) ) kz_stop = lz

      IF ( (jt == its).OR.(jt == 0) ) THEN   ! Opening file and defining variable :
         PRINT *, ''
         PRINT *, ' --- GETVAR_2D: opening file '//TRIM(cf_in)//' for '//TRIM(cv_in)//' !'
         CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, idx_f),  crtn,cf_in,cv_in)
         CALL sherr( NF90_INQ_VARID(idx_f, cv_in, idx_v),  crtn,cf_in,cv_in)
      END IF

      ierr1 = NF90_GET_ATT(idx_f, idx_v, 'scale_factor', rsf) ; !lolo, ugnj at each time...
      ierr2 = NF90_GET_ATT(idx_f, idx_v, 'add_offset',   rao)

      IF ( (idx_f==0).AND.(idx_v==0) ) CALL print_err(crtn, ' PROBLEM #2 file and variable handle never created => '//TRIM(cv_in)//' in '//TRIM(cf_in))

      l_okay = .FALSE.
      DO WHILE ( .NOT. l_okay )

         IF ( jlev == 0 ) THEN    ! No levels
            IF ( jt == 0 ) THEN
               CALL sherr( NF90_GET_VAR(idx_f, idx_v, X, start=(/1,1/), count=(/ni,nj/)), &
                  &      crtn,cf_in,cv_in)
            ELSEIF ( jt > 0 ) THEN
               CALL sherr( NF90_GET_VAR(idx_f, idx_v, X, start=(/1,1,jt/), count=(/ni,nj,1/)), &
                  &      crtn,cf_in,cv_in)
            END IF
            l_okay = .TRUE. ! we can exit the WHILE loop...

         ELSEIF ( jlev > 0 ) THEN

            !! User possibnj specified jlev = 1 and there is not an existing level dimension:
            IF ( n3 == -1 ) THEN
               PRINT *, ' *** warning: ',trim(crtn),' => there is actualnj no levels for ', trim(cv_in),' in ',trim(cf_in)
               PRINT *, '              => fixing it...'
               jlev = 0 ! => should be treated at next "while" loop...
            ELSE
               IF ( jlev >  n3 ) CALL print_err(crtn, ' you want extract a level greater than max value')
               IF ( jt == 0 ) THEN
                  CALL sherr( NF90_GET_VAR(idx_f, idx_v, X, start=(/1,1,jlev/), count=(/ni,nj,1/)), crtn,cf_in,cv_in)
               ELSEIF ( jt > 0 ) THEN
                  CALL sherr( NF90_GET_VAR(idx_f, idx_v, X, start=(/1,1,jlev,jt/), count=(/ni,nj,1,1/)), crtn,cf_in,cv_in)
               END IF
               l_okay = .TRUE.  ! we can exit the WHILE loop...
            END IF
         END IF
      END DO

      IF (ierr1 == NF90_NOERR) X = rsf*X
      IF (ierr2 == NF90_NOERR) X = X + rao

      IF ( ( (jt == ite ).OR.(jt == 0) ).AND.( (jlev == kz_stop).OR.(kz_stop == 0) ) )  THEN
         PRINT *, ' --- GETVAR_2D: closing file '//TRIM(cf_in)//' !'
         CALL sherr( NF90_CLOSE(idx_f),  crtn,cf_in,cv_in)
         idx_f = 0 ; idx_v = 0
         PRINT *, ''
      END IF

   END SUBROUTINE GETVAR_2D_R8
















   SUBROUTINE GETVAR_3D_R4(idx_f, idx_v, cf_in, cv_in, Nt, jt, X, jt1, jt2, jz1, jz2)

      !!------------------------------------------------------------------
      !! This routine extract a 3D field from a netcdf file
      !! at a given time
      !!
      !! INPUT :
      !! -------
      !!          * idx_f    : ID of current file                  (integer)
      !!          * idx_v    : ID of current variable              (integer)
      !!          * cf_in      : name of the input file              (character)
      !!          * cv_in      : name of the variable                (character)
      !!          * Nt        : number of time records     (integer)
      !!
      !!          * jt        : time snapshot to extract            (integer)
      !!                      0 => input file does not have a time snapshot
      !!                           (= old GETVAR_2D_NOTIME)
      !!
      !! OUTPUT :
      !! --------
      !!          * X         : 3D array contening the variable     (real)
      !!
      !! OPTIONAL INPUT :
      !! ----------------
      !!          * jt1, jt2  : first and last time snapshot to extract
      !!          * jz1, jz2  : first and last levels to extract
      !!
      !!------------------------------------------------------------------------

      INTEGER,                    INTENT(inout) :: idx_f, idx_v
      CHARACTER(len=*),           INTENT(in)    :: cf_in, cv_in
      INTEGER,                    INTENT(in)    :: Nt, jt
      REAL(4), DIMENSION (:,:,:), INTENT(out)   :: X
      INTEGER,       OPTIONAL,    INTENT(in)    :: jt1, jt2
      INTEGER,       OPTIONAL,    INTENT(in)    :: jz1, jz2

      INTEGER ::  &
         & ni,  &        ! x dimension of the variable        (integer)
         & nj,  &        ! y dimension of the variable        (integer)
         & lz            ! z dimension of the variable        (integer)
      INTEGER :: n1, n2, n3, n4, its, ite, izs, ize
      INTEGER :: ierr1, ierr2
      REAL(4) :: rsf, rao
      CHARACTER(len=80), PARAMETER :: crtn = 'GETVAR_3D_R4'

      ni = size(X,1)
      nj = size(X,2)
      lz = size(X,3)

      CALL DIMS(cf_in, cv_in, n1, n2, n3, n4)

      IF ( (ni /= n1).OR.(nj /= n2) ) CALL print_err(crtn, ' PROBLEM #1 => '//TRIM(cv_in)//' in '//TRIM(cf_in))

      its = 1 ; ite = Nt
      IF ( PRESENT(jt1) ) its = jt1
      IF ( PRESENT(jt2) ) ite = jt2

      IF ( PRESENT(jz1).AND.PRESENT(jz2) ) THEN
         izs = jz1 ; ize = jz2
      ELSE
         izs = 1   ; ize = lz
      END IF

      IF ( (jt == its).OR.(jt == 0) ) THEN   ! Opening file and defining variable :
         PRINT *, ''
         PRINT *, ' --- GETVAR_3D_R4: opening file '//TRIM(cf_in)//' for '//TRIM(cv_in)//' !'
         CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE,  idx_f),  crtn,cf_in,cv_in)
         CALL sherr( NF90_INQ_VARID(idx_f, cv_in, idx_v),  crtn,cf_in,cv_in)
      END IF

      IF ( (idx_f==0).AND.(idx_v==0) ) CALL print_err(crtn, ' PROBLEM #2 file and variable handle never created => '//TRIM(cv_in)//' in '//TRIM(cf_in))

      ierr1 = NF90_GET_ATT(idx_f, idx_v, 'scale_factor', rsf) ; !lolo, ugnj at each time...
      ierr2 = NF90_GET_ATT(idx_f, idx_v, 'add_offset',   rao)

      IF ( jt == 0 ) THEN
         CALL sherr( NF90_GET_VAR(idx_f, idx_v, X, start=(/1,1,izs/), count=(/ni,nj,ize/)), &
            &      crtn,cf_in,cv_in)
      ELSEIF ( jt > 0 ) THEN
         CALL sherr( NF90_GET_VAR(idx_f, idx_v, X, start=(/1,1,izs,jt/), count=(/ni,nj,ize,1/)), &
            &      crtn,cf_in,cv_in)
      END IF

      IF (ierr1 == NF90_NOERR) X = rsf*X
      IF (ierr2 == NF90_NOERR) X = X + rao

      IF ( ( jt == ite ).OR.( jt == 0 ) )  THEN
         PRINT *, ' --- GETVAR_3D_R4: closing file '//TRIM(cf_in)//' !'
         CALL sherr( NF90_CLOSE(idx_f),  crtn,cf_in,cv_in)
         idx_f = 0 ; idx_v = 0
      END IF

   END SUBROUTINE GETVAR_3D_R4

   SUBROUTINE GETVAR_3D_R8(idx_f, idx_v, cf_in, cv_in, Nt, jt, X, jt1, jt2, jz1, jz2)
      INTEGER,                    INTENT(inout) :: idx_f, idx_v
      CHARACTER(len=*),           INTENT(in)    :: cf_in, cv_in
      INTEGER,                    INTENT(in)    :: Nt, jt
      REAL(8), DIMENSION (:,:,:), INTENT(out)   :: X
      INTEGER,       OPTIONAL,    INTENT(in)    :: jt1, jt2
      INTEGER,       OPTIONAL,    INTENT(in)    :: jz1, jz2
      INTEGER ::  &
         & ni,  &        ! x dimension of the variable        (integer)
         & nj,  &        ! y dimension of the variable        (integer)
         & lz            ! z dimension of the variable        (integer)
      INTEGER :: n1, n2, n3, n4, its, ite, izs, ize
      INTEGER :: ierr1, ierr2
      REAL(4) :: rsf, rao
      CHARACTER(len=80), PARAMETER :: crtn = 'GETVAR_3D_R8'
      ni = SIZE(X,1)
      nj = size(X,2)
      lz = size(X,3)
      CALL DIMS(cf_in, cv_in, n1, n2, n3, n4)
      IF ( (ni /= n1).OR.(nj /= n2) ) CALL print_err(crtn, ' PROBLEM #1 => '//TRIM(cv_in)//' in '//TRIM(cf_in))
      its = 1 ; ite = Nt
      IF ( PRESENT(jt1) ) its = jt1
      IF ( PRESENT(jt2) ) ite = jt2
      IF ( PRESENT(jz1).AND.PRESENT(jz2) ) THEN
         izs = jz1 ; ize = jz2
      ELSE
         izs = 1   ; ize = lz
      END IF
      IF ( (jt == its).OR.(jt == 0) ) THEN   ! Opening file and defining variable :
         PRINT *, ''
         PRINT *, ' --- GETVAR_3D_R8: opening file '//TRIM(cf_in)//' for '//TRIM(cv_in)//' !'
         CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE,  idx_f),  crtn,cf_in,cv_in)
         CALL sherr( NF90_INQ_VARID(idx_f, cv_in, idx_v),  crtn,cf_in,cv_in)
      END IF
      IF ( (idx_f==0).AND.(idx_v==0) ) CALL print_err(crtn, ' PROBLEM #2 file and variable handle never created => '//TRIM(cv_in)//' in '//TRIM(cf_in))
      ierr1 = NF90_GET_ATT(idx_f, idx_v, 'scale_factor', rsf) ; !lolo, ugnj at each time...
      ierr2 = NF90_GET_ATT(idx_f, idx_v, 'add_offset',   rao)
      IF ( jt == 0 ) THEN
         CALL sherr( NF90_GET_VAR(idx_f, idx_v, X, start=(/1,1,izs/), count=(/ni,nj,ize/)), &
            &      crtn,cf_in,cv_in)
      ELSEIF ( jt > 0 ) THEN
         CALL sherr( NF90_GET_VAR(idx_f, idx_v, X, start=(/1,1,izs,jt/), count=(/ni,nj,ize,1/)), &
            &      crtn,cf_in,cv_in)
      END IF
      IF (ierr1 == NF90_NOERR) X = rsf*X
      IF (ierr2 == NF90_NOERR) X = X + rao
      IF ( ( jt == ite ).OR.( jt == 0 ) )  THEN
         PRINT *, ' --- GETVAR_3D_R8: closing file '//TRIM(cf_in)//' !'
         CALL sherr( NF90_CLOSE(idx_f),  crtn,cf_in,cv_in)
         idx_f = 0 ; idx_v = 0
      END IF
   END SUBROUTINE GETVAR_3D_R8







   !lili
   SUBROUTINE GETMASK_2D(cf_in, cv_in, IX, jlev)

      !!-----------------------------------------------------------------------
      !!  Get mask (variable 'cv_in') from a netcdf file.
      !! - mask is stored in integer array IX
      !!
      !! INPUT :
      !! -------
      !!          * cf_in      : name of the input file          (character)
      !!          * cv_in      : name of mask variable           (character)
      !!
      !!          * jlev      : level to get (0 if no levels) |OPTIONAL|     (integer)
      !!
      !! OUTPUT :
      !! --------
      !!          * IX        :  2D array contening mask       (integer)
      !!
      !!------------------------------------------------------------------------
      INTEGER                              :: id_f, id_v
      CHARACTER(len=*),        INTENT(in)  :: cf_in, cv_in
      INTEGER, OPTIONAL,       INTENT(in)  :: jlev
      INTEGER(1), DIMENSION(:,:), INTENT(out) :: IX

      INTEGER :: &
         & ni, &    ! x dimension of the mask
         & nj       ! y dimension of the mask

      INTEGER :: nx, ny, nz, nt, icz

      CHARACTER(len=80), PARAMETER :: crtn = 'GETMASK_2D'

      ni = size(IX,1)
      nj = size(IX,2)

      icz = 1 ! getting mask at level 1 for default

      IF ( present(jlev) ) THEN
         IF ( jlev > 0 ) THEN
            icz = jlev ; WRITE(6,*) 'Getting mask at level', icz
         ELSE
            CALL print_err(crtn, 'you cannot specify a level jlev <= 0')
         END IF
      END IF

      CALL DIMS(cf_in, cv_in, nx, ny, nz, nt)

      IF ( (nx /= ni).OR.(ny /= nj) ) CALL print_err(crtn, 'data and mask file dont have same horizontal dimensions')



      !!    Opening MASK netcdf file
      CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE,    id_f),  crtn,cf_in,cv_in)
      CALL sherr( NF90_INQ_VARID(id_f, trim(cv_in), id_v),  crtn,cf_in,cv_in)


      IF ( nz > 0 ) THEN

         !! Mask is 3D
         !! ~~~~~~~~~~

         IF ( .NOT. present(jlev) ) THEN
            WRITE(6,*) TRIM(crtn)//': WARNING => mask is 3D! what level to read? => taking level #1 !'
         END IF

         !!  3D+T
         IF ( nt > 0 ) THEN
            CALL sherr( NF90_GET_VAR(id_f, id_v, IX, start=(/1,1,icz,1/), count=(/nx,ny,1,1/)),  &
               &      crtn,cf_in,cv_in)
         ELSE
            !! 3D
            CALL sherr( NF90_GET_VAR(id_f, id_v, IX, start=(/1,1,icz/), count=(/nx,ny,1/)), &
               &      crtn,cf_in,cv_in)
         END IF

      ELSE

         !! Mask is 2D
         !! ~~~~~~~~~~
         IF ( present(jlev) ) THEN
            IF (jlev > 1) CALL print_err(crtn, 'you want mask at a given level (jlev > 1) but mask is 2D!!!')
         END IF


         IF ( nt > 0 ) THEN
            !!  2D+T
            CALL sherr( NF90_GET_VAR(id_f, id_v, IX, start=(/1,1,1/), count=(/nx,ny,1/)), crtn,cf_in,cv_in)
         ELSE
            !! 2D
            CALL sherr( NF90_GET_VAR(id_f, id_v, IX),  crtn,cf_in,cv_in)
         END IF


      END IF

      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_in,cv_in)

   END SUBROUTINE GETMASK_2D





   SUBROUTINE GETMASK_3D(cf_in, cv_in, IX, jz1, jz2)

      !!-----------------------------------------------------------------------
      !!  Get mask (variable 'cv_in') from a netcdf file.
      !! - mask is stored in integer array IX
      !!
      !! INPUT :
      !! -------
      !!          * cf_in      : name of the input file          (character)
      !!          * cv_in      : name of mask variable           (character)
      !!
      !! OUTPUT :
      !! --------
      !!          * IX        :  3D array contening mask       (integer)
      !!
      !! OPTIONAL INPUT :
      !! ----------------
      !!          * jz1, jz2  : first and last levels to extract
      !!
      !!------------------------------------------------------------------------
      INTEGER                                :: id_f, id_v
      CHARACTER(len=*),          INTENT(in)  :: cf_in, cv_in
      INTEGER(1), DIMENSION(:,:,:), INTENT(out) :: IX
      INTEGER,       OPTIONAL,    INTENT(in)    :: jz1, jz2

      INTEGER :: &
         & ni, &    ! x dimension of the mask
         & nj, &    ! y dimension of the mask
         & lz       ! z dimension of the mask

      INTEGER :: nx, ny, nz, nt, izs, ize

      CHARACTER(len=80), PARAMETER :: crtn = 'GETMASK_3D'

      ni = size(IX,1)
      nj = size(IX,2)
      lz = size(IX,3)

      IF ( PRESENT(jz1).AND.PRESENT(jz2) ) THEN
         izs = jz1 ; ize = jz2
      ELSE
         izs = 1   ; ize = lz
      END IF

      CALL DIMS(cf_in, cv_in, nx, ny, nz, nt)

      IF ( nz < 1 ) THEN
         WRITE(6,*) 'mask 3D file => ', trim(cf_in)
         WRITE(6,*) 'mask 3D name => ', trim(cv_in)
         CALL print_err(crtn, 'mask is not 3D')
      END IF

      IF ( (nx /= ni).OR.(ny /= nj).OR.((ize-izs+1) /= lz) ) THEN
         !&   CALL print_err(crtn, 'data and mask file dont have same dimensions => '\\)
         PRINT *, ''
         WRITE(6,*) 'ERROR in ',TRIM(crtn),' (io_ezcdf.f90): '
         WRITE(6,*) 'data and mask file dont have same dimensions'
         WRITE(6,*) '  => nx, ny, nz =', nx, ny, nz
         WRITE(6,*) '  => ni, nj, lz =', ni, nj, lz
         PRINT *, ''
         STOP
      END IF

      !!    Opening MASK netcdf file
      CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, id_f),     crtn,cf_in,cv_in)
      CALL sherr( NF90_INQ_VARID(id_f, trim(cv_in), id_v),  crtn,cf_in,cv_in)

      IF ( nt > 0 ) THEN
         !! 3D+T
         CALL sherr( NF90_GET_VAR(id_f, id_v, IX, start=(/1,1,izs,1/), count=(/nx,ny,ize,1/)), &
            &      crtn,cf_in,cv_in)

      ELSE
         !! 3D
         CALL sherr( NF90_GET_VAR(id_f, id_v, IX, start=(/1,1,izs/), count=(/nx,ny,ize/)),  crtn,cf_in,cv_in)
      END IF

      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_in,cv_in)

   END SUBROUTINE GETMASK_3D



   
   SUBROUTINE PT_SERIES(vtime, vdt01, cf_in, cv_t, cv_dt01, cun01, cln01, vflag, &
      &                 ct_unit, ct_clnd,              &
      &                 vdt02, cv_dt02, cun02, cln02,  &
      &                 vdt03, cv_dt03, cun03, cln03,  &
      &                 vdt04, cv_dt04, cun04, cln04,  &
      &                 vdt05, cv_dt05, cun05, cln05,  &
      &                 vdt06, cv_dt06, cun06, cln06,  &
      &                 vdt07, cv_dt07, cun07, cln07,  &
      &                 vdt08, cv_dt08, cun08, cln08,  &
      &                 vdt09, cv_dt09, cun09, cln09,  &
      &                 vdt10, cv_dt10, cun10, cln10,  &
      &                 vdt11, cv_dt11, cun11, cln11,  &
      &                 vdt12, cv_dt12, cun12, cln12,  &
      &                 vdt13, cv_dt13, cun13, cln13,  &
      &                 vdt14, cv_dt14, cun14, cln14,  &
      &                 vdt15, cv_dt15, cun15, cln15,  &
      &                 vdt16, cv_dt16, cun16, cln16,  &
      &                 vdt17, cv_dt17, cun17, cln17,  &
      &                 vdt18, cv_dt18, cun18, cln18,  &
      &                 vdt19, cv_dt19, cun19, cln19,  &
      &                 vdt20, cv_dt20, cun20, cln20   )

      !! INPUT :
      !! -------
      !!        vtime  = time array                               [array 1D real8]
      !!        vdt01 = 1D array containing time-series         [array 1D real4]
      !!        cf_in  = name of the output file                  [character]
      !!        cv_t = name of time                               [character]
      !!        cv_dt01  = name of the variable                     [character]
      !!        cun01  = unit for treated variable                [character]
      !!        cln01 = long-name for treated variable              [character]
      !!        vflag = flag value or "0."                        [real]
      !!
      !!        ct_unit = time unit
      !!
      !!--------------------------------------------------------------------------

      REAL(8), DIMENSION(:),     INTENT(in)   :: vtime
      REAL(4), DIMENSION(:),      INTENT(in)  :: vdt01
      CHARACTER(len=*),           INTENT(in)  :: cf_in, cv_t, cv_dt01, cun01, cln01
      REAL(4),                    INTENT(in)  :: vflag
      CHARACTER(len=*), OPTIONAL, INTENT(in)  :: ct_unit, ct_clnd
      REAL(4), DIMENSION(:), OPTIONAL, INTENT(in)  :: vdt02, vdt03, vdt04, vdt05, vdt06, vdt07, vdt08, &
         &                                            vdt09, vdt10, vdt11, vdt12, vdt13, vdt14, vdt15, &
         &                                            vdt16, vdt17, vdt18, vdt19, vdt20
      CHARACTER(len=*),      OPTIONAL, INTENT(in)  :: cv_dt02, cv_dt03, cv_dt04, cv_dt05, cv_dt06, cv_dt07, cv_dt08, &
         &                                            cv_dt09, cv_dt10, cv_dt11, cv_dt12, cv_dt13, cv_dt14, cv_dt15, &
         &                                            cv_dt16, cv_dt17, cv_dt18, cv_dt19, cv_dt20
      CHARACTER(len=*),      OPTIONAL, INTENT(in)  :: cun02, cun03, cun04, cun05, cun06, cun07, cun08, &
         &                                            cun09, cun10, cun11, cun12, cun13, cun14, cun15, &
         &                                            cun16, cun17, cun18, cun19, cun20
      CHARACTER(len=*),      OPTIONAL, INTENT(in) ::  cln02, cln03, cln04, cln05, cln06, cln07, cln08,  &
         &                                            cln09, cln10, cln11, cln12, cln13, cln14, cln15,  &
         &                                            cln16, cln17, cln18, cln19, cln20
      !!
      INTEGER :: idf, idtd, idt, nbt, jt
      INTEGER :: idv01, idv02, idv03, idv04, idv05, idv06, idv07, idv08, idv09, idv10, &
         &       idv11, idv12, idv13, idv14, idv15, idv16, idv17, idv18, idv19, idv20      
      REAL(4) :: rmin, rmax
      LOGICAL :: ldv02=.FALSE.,ldv03=.FALSE.,ldv04=.FALSE.,ldv05=.FALSE.,ldv06=.FALSE.,ldv07=.FALSE.,ldv08=.FALSE.,ldv09=.FALSE., &
         &       ldv10=.FALSE.,ldv11=.FALSE.,ldv12=.FALSE.,ldv13=.FALSE.,ldv14=.FALSE.,ldv15=.FALSE., &
         &       ldv16=.FALSE.,ldv17=.FALSE.,ldv18=.FALSE.,ldv19=.FALSE.,ldv20=.FALSE.

      CHARACTER(len=80), PARAMETER :: crtn = 'PT_SERIES'

      IF (PRESENT(vdt02)) ldv02=.true.
      IF (PRESENT(vdt03)) ldv03=.true.
      IF (PRESENT(vdt04)) ldv04=.true.
      IF (PRESENT(vdt05)) ldv05=.true.
      IF (PRESENT(vdt06)) ldv06=.true.
      IF (PRESENT(vdt07)) ldv07=.true.
      IF (PRESENT(vdt08)) ldv08=.true.
      IF (PRESENT(vdt09)) ldv09=.true.
      IF (PRESENT(vdt10)) ldv10=.true.
      IF (PRESENT(vdt11)) ldv11=.true.
      IF (PRESENT(vdt12)) ldv12=.true.
      IF (PRESENT(vdt13)) ldv13=.true.
      IF (PRESENT(vdt14)) ldv14=.true.
      IF (PRESENT(vdt15)) ldv15=.true.
      IF (PRESENT(vdt16)) ldv16=.true.
      IF (PRESENT(vdt17)) ldv17=.true.
      IF (PRESENT(vdt18)) ldv18=.true.
      IF (PRESENT(vdt19)) ldv19=.true.
      IF (PRESENT(vdt20)) ldv20=.true.

      nbt = SIZE(vtime,1)
      
      IF (                     SIZE(vdt01,1)/=nbt ) CALL print_err(crtn, 'Time vec and series vec #01 dont agree in size! => '//TRIM(cv_dt01))
      IF ( ldv02 ) THEN ; IF ( SIZE(vdt02,1)/=nbt ) CALL print_err(crtn, 'Time vec and series vec #02 dont agree in size! => '//TRIM(cv_dt02)); ENDIF
      IF ( ldv03 ) THEN ; IF ( SIZE(vdt03,1)/=nbt ) CALL print_err(crtn, 'Time vec and series vec #03 dont agree in size! => '//TRIM(cv_dt03)); ENDIF
      IF ( ldv04 ) THEN ; IF ( SIZE(vdt04,1)/=nbt ) CALL print_err(crtn, 'Time vec and series vec #04 dont agree in size! => '//TRIM(cv_dt04)); ENDIF
      IF ( ldv05 ) THEN ; IF ( SIZE(vdt05,1)/=nbt ) CALL print_err(crtn, 'Time vec and series vec #05 dont agree in size! => '//TRIM(cv_dt05)); ENDIF
      IF ( ldv06 ) THEN ; IF ( SIZE(vdt06,1)/=nbt ) CALL print_err(crtn, 'Time vec and series vec #06 dont agree in size! => '//TRIM(cv_dt06)); ENDIF
      IF ( ldv07 ) THEN ; IF ( SIZE(vdt07,1)/=nbt ) CALL print_err(crtn, 'Time vec and series vec #07 dont agree in size! => '//TRIM(cv_dt07)); ENDIF
      IF ( ldv08 ) THEN ; IF ( SIZE(vdt08,1)/=nbt ) CALL print_err(crtn, 'Time vec and series vec #08 dont agree in size! => '//TRIM(cv_dt08)); ENDIF
      IF ( ldv09 ) THEN ; IF ( SIZE(vdt09,1)/=nbt ) CALL print_err(crtn, 'Time vec and series vec #09 dont agree in size! => '//TRIM(cv_dt09)); ENDIF
      IF ( ldv10 ) THEN ; IF ( SIZE(vdt10,1)/=nbt ) CALL print_err(crtn, 'Time vec and series vec #10 dont agree in size! => '//TRIM(cv_dt10)); ENDIF
      IF ( ldv11 ) THEN ; IF ( SIZE(vdt11,1)/=nbt ) CALL print_err(crtn, 'Time vec and series vec #11 dont agree in size! => '//TRIM(cv_dt11)); ENDIF
      IF ( ldv12 ) THEN ; IF ( SIZE(vdt12,1)/=nbt ) CALL print_err(crtn, 'Time vec and series vec #12 dont agree in size! => '//TRIM(cv_dt12)); ENDIF
      IF ( ldv13 ) THEN ; IF ( SIZE(vdt13,1)/=nbt ) CALL print_err(crtn, 'Time vec and series vec #13 dont agree in size! => '//TRIM(cv_dt13)); ENDIF
      IF ( ldv14 ) THEN ; IF ( SIZE(vdt14,1)/=nbt ) CALL print_err(crtn, 'Time vec and series vec #14 dont agree in size! => '//TRIM(cv_dt14)); ENDIF
      IF ( ldv15 ) THEN ; IF ( SIZE(vdt15,1)/=nbt ) CALL print_err(crtn, 'Time vec and series vec #15 dont agree in size! => '//TRIM(cv_dt15)); ENDIF
      IF ( ldv16 ) THEN ; IF ( SIZE(vdt16,1)/=nbt ) CALL print_err(crtn, 'Time vec and series vec #16 dont agree in size! => '//TRIM(cv_dt16)); ENDIF
      IF ( ldv17 ) THEN ; IF ( SIZE(vdt17,1)/=nbt ) CALL print_err(crtn, 'Time vec and series vec #17 dont agree in size! => '//TRIM(cv_dt17)); ENDIF
      IF ( ldv18 ) THEN ; IF ( SIZE(vdt18,1)/=nbt ) CALL print_err(crtn, 'Time vec and series vec #18 dont agree in size! => '//TRIM(cv_dt18)); ENDIF
      IF ( ldv19 ) THEN ; IF ( SIZE(vdt19,1)/=nbt ) CALL print_err(crtn, 'Time vec and series vec #19 dont agree in size! => '//TRIM(cv_dt19)); ENDIF
      IF ( ldv20 ) THEN ; IF ( SIZE(vdt20,1)/=nbt ) CALL print_err(crtn, 'Time vec and series vec #20 dont agree in size! => '//TRIM(cv_dt20)); ENDIF

      IF ( vflag /= 0.) THEN
         rmin =  1.E6 ; rmax = -1.E6
         DO jt = 1, nbt
            IF ((vdt01(jt) <= rmin).and.(vdt01(jt) /= vflag)) rmin = vdt01(jt)
            IF ((vdt01(jt) >= rmax).and.(vdt01(jt) /= vflag)) rmax = vdt01(jt)
         END DO
      ELSE
         rmin = minval(vdt01) ; rmax = maxval(vdt01)
      END IF
      vextrema(3,:) = (/MINVAL(vtime),MAXVAL(vtime)/)

      !!           CREATE NETCDF OUTPUT FILE :
      CALL sherr( NF90_CREATE(cf_in, NF90_NETCDF4, idf), crtn,cf_in,cv_dt01)

      !! Time
      CALL sherr( NF90_DEF_DIM(idf, TRIM(cv_t), NF90_UNLIMITED, idtd),                       crtn,cf_in,cv_t)
      CALL sherr( NF90_DEF_VAR(idf, TRIM(cv_t), NF90_DOUBLE, idtd, idt, deflate_level=9),    crtn,cf_in,cv_t)
      IF ( PRESENT(ct_unit) ) CALL sherr( NF90_PUT_ATT(idf, idt, 'units',    TRIM(ct_unit)), crtn,cf_in,cv_t)
      IF ( PRESENT(ct_clnd) ) CALL sherr( NF90_PUT_ATT(idf, idt, 'calendar', TRIM(ct_clnd)), crtn,cf_in,cv_t)
      CALL sherr( NF90_PUT_ATT(idf, idt, 'valid_min', vextrema(3,1)),                        crtn,cf_in,cv_t)
      CALL sherr( NF90_PUT_ATT(idf, idt, 'valid_max', vextrema(3,2)),                        crtn,cf_in,cv_t)

      !! Variable(s):
      CALL              sherr( NF90_DEF_VAR(idf, TRIM(cv_dt01), NF90_FLOAT, idtd, idv01, deflate_level=9), crtn,cf_in,cv_dt01 )
      IF ( ldv02 ) CALL sherr( NF90_DEF_VAR(idf, TRIM(cv_dt02), NF90_FLOAT, idtd, idv02, deflate_level=9), crtn,cf_in,cv_dt02 )
      IF ( ldv03 ) CALL sherr( NF90_DEF_VAR(idf, TRIM(cv_dt03), NF90_FLOAT, idtd, idv03, deflate_level=9), crtn,cf_in,cv_dt03 )
      IF ( ldv04 ) CALL sherr( NF90_DEF_VAR(idf, TRIM(cv_dt04), NF90_FLOAT, idtd, idv04, deflate_level=9), crtn,cf_in,cv_dt04 )
      IF ( ldv05 ) CALL sherr( NF90_DEF_VAR(idf, TRIM(cv_dt05), NF90_FLOAT, idtd, idv05, deflate_level=9), crtn,cf_in,cv_dt05 )
      IF ( ldv06 ) CALL sherr( NF90_DEF_VAR(idf, TRIM(cv_dt06), NF90_FLOAT, idtd, idv06, deflate_level=9), crtn,cf_in,cv_dt06 )
      IF ( ldv07 ) CALL sherr( NF90_DEF_VAR(idf, TRIM(cv_dt07), NF90_FLOAT, idtd, idv07, deflate_level=9), crtn,cf_in,cv_dt07 )
      IF ( ldv08 ) CALL sherr( NF90_DEF_VAR(idf, TRIM(cv_dt08), NF90_FLOAT, idtd, idv08, deflate_level=9), crtn,cf_in,cv_dt08 )
      IF ( ldv09 ) CALL sherr( NF90_DEF_VAR(idf, TRIM(cv_dt09), NF90_FLOAT, idtd, idv09, deflate_level=9), crtn,cf_in,cv_dt09 )
      IF ( ldv10 ) CALL sherr( NF90_DEF_VAR(idf, TRIM(cv_dt10), NF90_FLOAT, idtd, idv10, deflate_level=9), crtn,cf_in,cv_dt10 )
      IF ( ldv11 ) CALL sherr( NF90_DEF_VAR(idf, TRIM(cv_dt11), NF90_FLOAT, idtd, idv11, deflate_level=9), crtn,cf_in,cv_dt11 )
      IF ( ldv12 ) CALL sherr( NF90_DEF_VAR(idf, TRIM(cv_dt12), NF90_FLOAT, idtd, idv12, deflate_level=9), crtn,cf_in,cv_dt12 )
      IF ( ldv13 ) CALL sherr( NF90_DEF_VAR(idf, TRIM(cv_dt13), NF90_FLOAT, idtd, idv13, deflate_level=9), crtn,cf_in,cv_dt13 )
      IF ( ldv14 ) CALL sherr( NF90_DEF_VAR(idf, TRIM(cv_dt14), NF90_FLOAT, idtd, idv14, deflate_level=9), crtn,cf_in,cv_dt14 )
      IF ( ldv15 ) CALL sherr( NF90_DEF_VAR(idf, TRIM(cv_dt15), NF90_FLOAT, idtd, idv15, deflate_level=9), crtn,cf_in,cv_dt15 )
      IF ( ldv16 ) CALL sherr( NF90_DEF_VAR(idf, TRIM(cv_dt16), NF90_FLOAT, idtd, idv16, deflate_level=9), crtn,cf_in,cv_dt16 )
      IF ( ldv17 ) CALL sherr( NF90_DEF_VAR(idf, TRIM(cv_dt17), NF90_FLOAT, idtd, idv17, deflate_level=9), crtn,cf_in,cv_dt17 )
      IF ( ldv18 ) CALL sherr( NF90_DEF_VAR(idf, TRIM(cv_dt18), NF90_FLOAT, idtd, idv18, deflate_level=9), crtn,cf_in,cv_dt18 )
      IF ( ldv19 ) CALL sherr( NF90_DEF_VAR(idf, TRIM(cv_dt19), NF90_FLOAT, idtd, idv19, deflate_level=9), crtn,cf_in,cv_dt19 )
      IF ( ldv20 ) CALL sherr( NF90_DEF_VAR(idf, TRIM(cv_dt20), NF90_FLOAT, idtd, idv20, deflate_level=9), crtn,cf_in,cv_dt20 )

      !! V01:
      CALL sherr( NF90_PUT_ATT(idf, idv01, 'long_name', trim(cln01) ),  crtn,cf_in,cv_dt01)
      CALL sherr( NF90_PUT_ATT(idf, idv01, 'units',     trim(cun01) ),  crtn,cf_in,cv_dt01)
      IF ( vflag /= 0. ) CALL sherr( NF90_PUT_ATT(idf, idv01,trim(cmv0),vflag),  crtn,cf_in,cv_dt01)
      !CALL sherr( NF90_PUT_ATT(idf, idv01,'actual_range', (/rmin,rmax/)),  crtn,cf_in,cv_dt01)
      CALL sherr( NF90_PUT_ATT(idf, NF90_GLOBAL, 'About', trim(cabout)),  crtn,cf_in,cv_dt01)

      !! V02:
      IF ( ldv02 ) THEN
         CALL sherr( NF90_PUT_ATT(idf, idv02, 'long_name', TRIM(cln02) ),  crtn,cf_in,cv_dt02)
         CALL sherr( NF90_PUT_ATT(idf, idv02, 'units',     TRIM(cun02) ),  crtn,cf_in,cv_dt02)
         IF ( vflag /= 0. ) CALL sherr( NF90_PUT_ATT(idf, idv02,trim(cmv0),vflag),  crtn,cf_in,cv_dt02)
      END IF
      !! V03:
      IF ( ldv03 ) THEN
         CALL sherr( NF90_PUT_ATT(idf, idv03, 'long_name', TRIM(cln03) ),  crtn,cf_in,cv_dt03)
         CALL sherr( NF90_PUT_ATT(idf, idv03, 'units',     TRIM(cun03) ),  crtn,cf_in,cv_dt03)
         IF ( vflag /= 0. ) CALL sherr( NF90_PUT_ATT(idf, idv03,trim(cmv0),vflag),  crtn,cf_in,cv_dt03)
      END IF
      !! V04:
      IF ( ldv04 ) THEN
         CALL sherr( NF90_PUT_ATT(idf, idv04, 'long_name', TRIM(cln04) ),  crtn,cf_in,cv_dt04)
         CALL sherr( NF90_PUT_ATT(idf, idv04, 'units',     TRIM(cun04) ),  crtn,cf_in,cv_dt04)
         IF ( vflag /= 0. ) CALL sherr( NF90_PUT_ATT(idf, idv04,trim(cmv0),vflag),  crtn,cf_in,cv_dt04)
      END IF
      !! V05:
      IF ( ldv05 ) THEN
         CALL sherr( NF90_PUT_ATT(idf, idv05, 'long_name', TRIM(cln05) ),  crtn,cf_in,cv_dt05)
         CALL sherr( NF90_PUT_ATT(idf, idv05, 'units',     TRIM(cun05) ),  crtn,cf_in,cv_dt05)
         IF ( vflag /= 0. ) CALL sherr( NF90_PUT_ATT(idf, idv05,trim(cmv0),vflag),  crtn,cf_in,cv_dt05)
      END IF
      !! V06:
      IF ( ldv06 ) THEN
         CALL sherr( NF90_PUT_ATT(idf, idv06, 'long_name', TRIM(cln06) ),  crtn,cf_in,cv_dt06)
         CALL sherr( NF90_PUT_ATT(idf, idv06, 'units',     TRIM(cun06) ),  crtn,cf_in,cv_dt06)
         IF ( vflag /= 0. ) CALL sherr( NF90_PUT_ATT(idf, idv06,trim(cmv0),vflag),  crtn,cf_in,cv_dt06)
      END IF
      !! V07:
      IF ( ldv07 ) THEN
         CALL sherr( NF90_PUT_ATT(idf, idv07, 'long_name', TRIM(cln07) ),  crtn,cf_in,cv_dt07)
         CALL sherr( NF90_PUT_ATT(idf, idv07, 'units',     TRIM(cun07) ),  crtn,cf_in,cv_dt07)
         IF ( vflag /= 0. ) CALL sherr( NF90_PUT_ATT(idf, idv07,trim(cmv0),vflag),  crtn,cf_in,cv_dt07)
      END IF
      !! V08:
      IF ( ldv08 ) THEN
         CALL sherr( NF90_PUT_ATT(idf, idv08, 'long_name', TRIM(cln08) ),  crtn,cf_in,cv_dt08)
         CALL sherr( NF90_PUT_ATT(idf, idv08, 'units',     TRIM(cun08) ),  crtn,cf_in,cv_dt08)
         IF ( vflag /= 0. ) CALL sherr( NF90_PUT_ATT(idf, idv08,trim(cmv0),vflag),  crtn,cf_in,cv_dt08)
      END IF
      !! V09:
      IF ( ldv09 ) THEN
         CALL sherr( NF90_PUT_ATT(idf, idv09, 'long_name', TRIM(cln09) ),  crtn,cf_in,cv_dt09)
         CALL sherr( NF90_PUT_ATT(idf, idv09, 'units',     TRIM(cun09) ),  crtn,cf_in,cv_dt09)
         IF ( vflag /= 0. ) CALL sherr( NF90_PUT_ATT(idf, idv09,trim(cmv0),vflag),  crtn,cf_in,cv_dt09)
      END IF
      !! V10:
      IF ( ldv10 ) THEN
         CALL sherr( NF90_PUT_ATT(idf, idv10, 'long_name', TRIM(cln10) ),  crtn,cf_in,cv_dt10)
         CALL sherr( NF90_PUT_ATT(idf, idv10, 'units',     TRIM(cun10) ),  crtn,cf_in,cv_dt10)
         IF ( vflag /= 0. ) CALL sherr( NF90_PUT_ATT(idf, idv10,trim(cmv0),vflag),  crtn,cf_in,cv_dt10)
      END IF
      !! V11:
      IF ( ldv11 ) THEN
         CALL sherr( NF90_PUT_ATT(idf, idv11, 'long_name', TRIM(cln11) ),  crtn,cf_in,cv_dt11)
         CALL sherr( NF90_PUT_ATT(idf, idv11, 'units',     TRIM(cun11) ),  crtn,cf_in,cv_dt11)
         IF ( vflag /= 0. ) CALL sherr( NF90_PUT_ATT(idf, idv11,trim(cmv0),vflag),  crtn,cf_in,cv_dt11)
      END IF
      !! V12:
      IF ( ldv12 ) THEN
         CALL sherr( NF90_PUT_ATT(idf, idv12, 'long_name', TRIM(cln12) ),  crtn,cf_in,cv_dt12)
         CALL sherr( NF90_PUT_ATT(idf, idv12, 'units',     TRIM(cun12) ),  crtn,cf_in,cv_dt12)
         IF ( vflag /= 0. ) CALL sherr( NF90_PUT_ATT(idf, idv12,trim(cmv0),vflag),  crtn,cf_in,cv_dt12)
      END IF
      !! V13:
      IF ( ldv13 ) THEN
         CALL sherr( NF90_PUT_ATT(idf, idv13, 'long_name', TRIM(cln13) ),  crtn,cf_in,cv_dt13)
         CALL sherr( NF90_PUT_ATT(idf, idv13, 'units',     TRIM(cun13) ),  crtn,cf_in,cv_dt13)
         IF ( vflag /= 0. ) CALL sherr( NF90_PUT_ATT(idf, idv13,trim(cmv0),vflag),  crtn,cf_in,cv_dt13)
      END IF
      !! V14:
      IF ( ldv14 ) THEN
         CALL sherr( NF90_PUT_ATT(idf, idv14, 'long_name', TRIM(cln14) ),  crtn,cf_in,cv_dt14)
         CALL sherr( NF90_PUT_ATT(idf, idv14, 'units',     TRIM(cun14) ),  crtn,cf_in,cv_dt14)
         IF ( vflag /= 0. ) CALL sherr( NF90_PUT_ATT(idf, idv14,trim(cmv0),vflag),  crtn,cf_in,cv_dt14)
      END IF
      !! V15:
      IF ( ldv15 ) THEN
         CALL sherr( NF90_PUT_ATT(idf, idv15, 'long_name', TRIM(cln15) ),  crtn,cf_in,cv_dt15)
         CALL sherr( NF90_PUT_ATT(idf, idv15, 'units',     TRIM(cun15) ),  crtn,cf_in,cv_dt15)
         IF ( vflag /= 0. ) CALL sherr( NF90_PUT_ATT(idf, idv15,trim(cmv0),vflag),  crtn,cf_in,cv_dt15)
      END IF
      !! V16:
      IF ( ldv16 ) THEN
         CALL sherr( NF90_PUT_ATT(idf, idv16, 'long_name', TRIM(cln16) ),  crtn,cf_in,cv_dt16)
         CALL sherr( NF90_PUT_ATT(idf, idv16, 'units',     TRIM(cun16) ),  crtn,cf_in,cv_dt16)
         IF ( vflag /= 0. ) CALL sherr( NF90_PUT_ATT(idf, idv16,trim(cmv0),vflag),  crtn,cf_in,cv_dt16)
      END IF
      !! V17:
      IF ( ldv17 ) THEN
         CALL sherr( NF90_PUT_ATT(idf, idv17, 'long_name', TRIM(cln17) ),  crtn,cf_in,cv_dt17)
         CALL sherr( NF90_PUT_ATT(idf, idv17, 'units',     TRIM(cun17) ),  crtn,cf_in,cv_dt17)
         IF ( vflag /= 0. ) CALL sherr( NF90_PUT_ATT(idf, idv17,trim(cmv0),vflag),  crtn,cf_in,cv_dt17)
      END IF
      !! V18:
      IF ( ldv18 ) THEN
         CALL sherr( NF90_PUT_ATT(idf, idv18, 'long_name', TRIM(cln18) ),  crtn,cf_in,cv_dt18)
         CALL sherr( NF90_PUT_ATT(idf, idv18, 'units',     TRIM(cun18) ),  crtn,cf_in,cv_dt18)
         IF ( vflag /= 0. ) CALL sherr( NF90_PUT_ATT(idf, idv18,trim(cmv0),vflag),  crtn,cf_in,cv_dt18)
      END IF
      !! V19:
      IF ( ldv19 ) THEN
         CALL sherr( NF90_PUT_ATT(idf, idv19, 'long_name', TRIM(cln19) ),  crtn,cf_in,cv_dt19)
         CALL sherr( NF90_PUT_ATT(idf, idv19, 'units',     TRIM(cun19) ),  crtn,cf_in,cv_dt19)
         IF ( vflag /= 0. ) CALL sherr( NF90_PUT_ATT(idf, idv19,trim(cmv0),vflag),  crtn,cf_in,cv_dt19)
      END IF
      !! V20:
      IF ( ldv20 ) THEN
         CALL sherr( NF90_PUT_ATT(idf, idv20, 'long_name', TRIM(cln20) ),  crtn,cf_in,cv_dt20)
         CALL sherr( NF90_PUT_ATT(idf, idv20, 'units',     TRIM(cun20) ),  crtn,cf_in,cv_dt20)
         IF ( vflag /= 0. ) CALL sherr( NF90_PUT_ATT(idf, idv20,trim(cmv0),vflag),  crtn,cf_in,cv_dt20)
      END IF

      
      CALL sherr( NF90_ENDDEF(idf),  crtn,cf_in,cv_dt01)

      !!       Write time variable :
      CALL sherr( NF90_PUT_VAR(idf, idt, vtime),    crtn,cf_in,cv_dt01)

      !!      WRITE VARIABLE
      CALL sherr( NF90_PUT_VAR(idf, idv01, vdt01),  crtn,cf_in,cv_dt01)
      IF ( ldv02 ) CALL sherr( NF90_PUT_VAR(idf, idv02, vdt02),  crtn,cf_in,cv_dt02)
      IF ( ldv03 ) CALL sherr( NF90_PUT_VAR(idf, idv03, vdt03),  crtn,cf_in,cv_dt03)
      IF ( ldv04 ) CALL sherr( NF90_PUT_VAR(idf, idv04, vdt04),  crtn,cf_in,cv_dt04)
      IF ( ldv05 ) CALL sherr( NF90_PUT_VAR(idf, idv05, vdt05),  crtn,cf_in,cv_dt05)
      IF ( ldv06 ) CALL sherr( NF90_PUT_VAR(idf, idv06, vdt06),  crtn,cf_in,cv_dt06)
      IF ( ldv07 ) CALL sherr( NF90_PUT_VAR(idf, idv07, vdt07),  crtn,cf_in,cv_dt07)
      IF ( ldv08 ) CALL sherr( NF90_PUT_VAR(idf, idv08, vdt08),  crtn,cf_in,cv_dt08)
      IF ( ldv09 ) CALL sherr( NF90_PUT_VAR(idf, idv09, vdt09),  crtn,cf_in,cv_dt09)
      IF ( ldv10 ) CALL sherr( NF90_PUT_VAR(idf, idv10, vdt10),  crtn,cf_in,cv_dt10)
      IF ( ldv11 ) CALL sherr( NF90_PUT_VAR(idf, idv11, vdt11),  crtn,cf_in,cv_dt11)
      IF ( ldv12 ) CALL sherr( NF90_PUT_VAR(idf, idv12, vdt12),  crtn,cf_in,cv_dt12)
      IF ( ldv13 ) CALL sherr( NF90_PUT_VAR(idf, idv13, vdt13),  crtn,cf_in,cv_dt13)
      IF ( ldv14 ) CALL sherr( NF90_PUT_VAR(idf, idv14, vdt14),  crtn,cf_in,cv_dt14)
      IF ( ldv15 ) CALL sherr( NF90_PUT_VAR(idf, idv15, vdt15),  crtn,cf_in,cv_dt15)
      IF ( ldv16 ) CALL sherr( NF90_PUT_VAR(idf, idv16, vdt16),  crtn,cf_in,cv_dt16)
      IF ( ldv17 ) CALL sherr( NF90_PUT_VAR(idf, idv17, vdt17),  crtn,cf_in,cv_dt17)
      IF ( ldv18 ) CALL sherr( NF90_PUT_VAR(idf, idv18, vdt18),  crtn,cf_in,cv_dt18)
      IF ( ldv19 ) CALL sherr( NF90_PUT_VAR(idf, idv19, vdt19),  crtn,cf_in,cv_dt19)
      IF ( ldv20 ) CALL sherr( NF90_PUT_VAR(idf, idv20, vdt20),  crtn,cf_in,cv_dt20)

      CALL sherr( NF90_CLOSE(idf),  crtn,cf_in,cv_dt01)

   END SUBROUTINE PT_SERIES



   

   SUBROUTINE P2D_T(idx_f, idx_v, Nt, lct, xlon, xlat, vtime, x2d, cf_in, &
      &             cv_lo, cv_la, cv_t, cv_in, vflag,      &
      &             attr_lon, attr_lat, attr_time, attr_F, &
      &             cextrainfo, l_add_valid_min_max)
      !!
      !! INPUT :
      !! -------
      !!        idx_f = ID of the file (takes its value on the first call)
      !!        idx_v = ID of the variable //
      !!        Nt    = t dimension of array to plot              [integer]
      !!        lct   = current time step                         [integer]
      !!        xlon  = 2D array of longitude  (nx,ny) or (nx,1)  [double]
      !!        xlat  = 2D array of latitude   (nx,ny) or (ny,1)  [double]
      !!        vtime  = time array                               [array 1D]
      !!        x2d = 2D snap of 2D+T array at time jt to write   [real]
      !!        cf_in  = name of the output file                  [character]
      !!        cv_lo = name of longitude                         [character]
      !!        cv_la = name of latitude                          [character]
      !!        cv_t = name of time                               [character]
      !!        cv_in  = name of the variable                     [character]
      !!        vflag = flag value or "0."                        [real]
      !!
      !!        cextrainfo = extra information to go in "Info" of header of netcdf
      !!        l_add_valid_min_max = for each variable write valid_min and valid_max (default=.true.)
      !!
      !!--------------------------------------------------------------------------
      !!
      INTEGER,                    INTENT(inout) :: idx_f, idx_v
      INTEGER,                    INTENT(in)    :: Nt, lct
      REAL(8), DIMENSION(:,:),    INTENT(in)    :: xlat, xlon
      REAL(4), DIMENSION(:,:),    INTENT(in)    :: x2d
      REAL(8), DIMENSION(Nt),     INTENT(in)    :: vtime
      CHARACTER(len=*),           INTENT(in)    :: cf_in, cv_lo, cv_la, cv_t, cv_in
      REAL(4),                    INTENT(in)    :: vflag
      !! Optional:
      TYPE(var_attr), DIMENSION(nbatt_max), OPTIONAL, INTENT(in) :: attr_lon, attr_lat, attr_time, attr_F
      CHARACTER(len=*), OPTIONAL, INTENT(in)    :: cextrainfo
      LOGICAL, OPTIONAL, INTENT(in)             :: l_add_valid_min_max

      INTEGER  :: id_x, id_y, id_t
      INTEGER  :: id_lo, id_la
      INTEGER  :: id_tim
      INTEGER  :: ni, nj
      REAL(4)  :: rmin, rmax
      LOGICAL  :: lcopy_att_F = .FALSE., l_add_extrema=.true.
      INTEGER, DIMENSION(:), ALLOCATABLE :: vidim

      CHARACTER(len=80), PARAMETER :: crtn = 'P2D_T'

      IF ( PRESENT(attr_F) ) lcopy_att_F = .TRUE.

      IF ( PRESENT(l_add_valid_min_max) ) l_add_extrema = l_add_valid_min_max

      !! About dimensions of xlon, xlat and x2d:
      cdt = TEST_XYZ(xlon, xlat, x2d)
      ni = size(x2d,1) ; nj = size(x2d,2)

      IF ( lct == 1 ) THEN

         IF ( vflag /= 0.) THEN
            rmin =  1.E6 ; rmax = -1.E6
            DO jj=1, nj
               DO ji=1, ni
                  IF ((x2d(ji,jj) <= rmin).and.(x2d(ji,jj) /= vflag)) rmin = x2d(ji,jj)
                  IF ((x2d(ji,jj) >= rmax).and.(x2d(ji,jj) /= vflag)) rmax = x2d(ji,jj)
               END DO
            END DO
         ELSE
            rmin = minval(x2d) ; rmax = maxval(x2d)
         END IF

      END IF ! lct == 1

      IF ( lct == 1 ) THEN

         vextrema(1,:) = (/minval(xlon),maxval(xlon)/); vextrema(2,:) = (/minval(xlat),maxval(xlat)/)
         vextrema(3,:) = (/minval(vtime),maxval(vtime)/)

         !! Opening mesh file for grid quest :
         !! ----------------------------------
         !!
         !!           CREATE NETCDF OUTPUT FILE :

         CALL sherr( NF90_CREATE(cf_in, NF90_NETCDF4, idx_f), crtn,cf_in,cv_in)
         !PRINT *, 'cdt = ', TRIM(cdt),'---'
         !PRINT *, 'cv_lo = ', TRIM(cv_lo),'---'
         !PRINT *, 'cv_la = ', TRIM(cv_la),'---'
         !PRINT *, 'cv_t = ', TRIM(cv_t),'---'
         !PRINT *, 'crtn = ', TRIM(crtn),'---'
         !PRINT *, 'cf_in = ', TRIM(cf_in),'---'
         !PRINT *, 'cv_in = ', TRIM(cv_in),'---'
         !PRINT *, 'attr_lat = ', attr_lat,'---'
         !PRINT *, 'attr_time = ', attr_time,'---'

         CALL prepare_nc(idx_f, cdt, ni, nj, cv_lo, cv_la, cv_t, vextrema, &
            &            id_x, id_y, id_t, id_lo, id_la, id_tim, crtn,cf_in,cv_in, &
            &            attr_lon=attr_lon, attr_lat=attr_lat, attr_tim=attr_time, &
            &            l_add_valid_min_max=l_add_extrema)
         !!
         !! Variable
         IF ( TRIM(cv_t) /= '' ) THEN
            ALLOCATE (vidim(3))
            vidim = (/id_x,id_y,id_t/)
         ELSE
            ALLOCATE (vidim(2))
            vidim = (/id_x,id_y/)
         END IF
         CALL sherr( NF90_DEF_VAR(idx_f, TRIM(cv_in), NF90_FLOAT, vidim, idx_v, deflate_level=9), &
            &      crtn,cf_in,cv_in )
         DEALLOCATE ( vidim )

         !!  VARIABLE ATTRIBUTES
         IF ( lcopy_att_F ) CALL SET_ATTRIBUTES_TO_VAR(idx_f, idx_v, attr_F,  crtn,cf_in,cv_in)
         ! Forcing these attributes (given in namelist):
         IF (vflag/=0.) CALL sherr( NF90_PUT_ATT(idx_f, idx_v,trim(cmv0),           vflag),  crtn,cf_in,cv_in)
         CALL                sherr( NF90_PUT_ATT(idx_f, idx_v,'actual_range', (/rmin,rmax/)),  crtn,cf_in,cv_in)
         CALL                sherr( NF90_PUT_ATT(idx_f, idx_v,'coordinates', &
            &                                TRIM(cv_t)//" "//TRIM(cv_la)//" "//TRIM(cv_lo)),  crtn,cf_in,cv_in)

         !! Global attributes
         IF ( PRESENT(cextrainfo) ) &
            CALL sherr( NF90_PUT_ATT(idx_f, NF90_GLOBAL, 'Info', TRIM(cextrainfo)),  crtn,cf_in,cv_in)
         CALL sherr( NF90_PUT_ATT(idx_f, NF90_GLOBAL, 'About', trim(cabout)),  crtn,cf_in,cv_in)

         !!           END OF DEFINITION
         CALL sherr( NF90_ENDDEF(idx_f),  crtn,cf_in,cv_in)

         !!       Write longitude variable :
         CALL sherr( NF90_PUT_VAR(idx_f, id_lo, xlon),  crtn,cf_in,cv_in)
         !!
         !!       Write latitude variable :
         CALL sherr( NF90_PUT_VAR(idx_f, id_la, xlat),  crtn,cf_in,cv_in)
         !!
         !!       Write time variable :
         IF ( TRIM(cv_t) /= '' ) CALL sherr( NF90_PUT_VAR(idx_f, id_tim, vtime),  crtn,cf_in,cv_in)
         !!
      END IF

      !!               WRITE VARIABLE
      CALL sherr( NF90_PUT_VAR(idx_f, idx_v, x2d,  start=(/1,1,lct/), count=(/ni,nj,1/)),  crtn,cf_in,cv_in)

      IF ( lct == Nt ) CALL sherr( NF90_CLOSE(idx_f),  crtn,cf_in,cv_in)

   END SUBROUTINE P2D_T





   SUBROUTINE P3D_T(idx_f, idx_v, Nt, lct, xlon, xlat, vdpth, vtime, x3d, cf_in, &
      &             cv_lo, cv_la, cv_dpth, cv_t, cv_in, vflag, &
      &             attr_lon, attr_lat, attr_z, attr_time, attr_F, &
      &             cextrainfo, l_add_valid_min_max)

      !! INPUT :
      !! -------
      !!        idx_f = ID of the file (takes its value on the first call)
      !!        idx_v = ID of the variable //
      !!        Nt    = t dimension of array to plot              [integer]
      !!        lct   = current time step                         [integer]
      !!        xlon  = 2D array of longitude  (nx,ny) or (nx,1)  [double]
      !!        xlat  = 2D array of latitude   (nx,ny) or (ny,1)  [double]
      !!        vdpth = depth array                               [array 1D double]
      !!        vtime  = time array                               [array 1D double]
      !!        x3d = 3D snap of 3D+T array at time jt to write   [real]
      !!        cf_in  = name of the output file                  [character]
      !!        cv_lo = name of longitude                         [character]
      !!        cv_la = name of latitude                          [character]
      !!        cv_dpth = name of depth                           [character]
      !!        cv_t = name of time                               [character]
      !!        cv_in  = name of the variable                     [character]
      !!        vflag = flag value or "0."                        [real]
      !!
      !!        cextrainfo = extra information to go in "Info" of header of netcdf
      !!--------------------------------------------------------------------------

      INTEGER,                    INTENT(inout) :: idx_f, idx_v
      INTEGER                                   :: id_dpt
      INTEGER,                    INTENT(in)    :: Nt, lct
      REAL(4), DIMENSION(:,:,:),  INTENT(in)    :: x3d
      REAL(8), DIMENSION(:,:),    INTENT(in)    :: xlat, xlon
      REAL(8), DIMENSION(:),      INTENT(in)    :: vdpth
      REAL(8), DIMENSION(Nt),     INTENT(in)    :: vtime
      CHARACTER(len=*),           INTENT(in)    :: cf_in, cv_lo, cv_la, cv_dpth, cv_t, cv_in
      REAL(4),                    INTENT(in)    :: vflag
      !! Optional:
      TYPE(var_attr), DIMENSION(nbatt_max), OPTIONAL, INTENT(in) :: attr_lon, attr_lat, attr_z, &
         &                                                          attr_time, attr_F
      CHARACTER(len=*), OPTIONAL, INTENT(in)    :: cextrainfo
      LOGICAL, OPTIONAL, INTENT(in)             :: l_add_valid_min_max

      INTEGER          :: id_z
      INTEGER          :: id_x, id_y, id_t, id_lo, id_la, id_tim
      INTEGER          :: ni, nj, lz
      REAL(4)          :: dr, rmin, rmax
      LOGICAL :: lcopy_att_F = .FALSE., &
         &       lcopy_att_z = .FALSE., &
         &     l_add_extrema = .TRUE.
      INTEGER, DIMENSION(:), ALLOCATABLE :: vidim

      CHARACTER(len=80), PARAMETER :: crtn = 'P3D_T'

      IF ( PRESENT(attr_z) ) lcopy_att_z = .TRUE.
      IF ( PRESENT(attr_F) ) lcopy_att_F = .TRUE.

      IF ( PRESENT(l_add_valid_min_max) ) l_add_extrema = l_add_valid_min_max

      !! About dimensions of xlon, xlat, vdpth and x3d:
      cdt = TEST_XYZ(xlon, xlat, x3d(:,:,1))

      ni = size(x3d,1) ; nj = size(x3d,2) ; lz = size(vdpth)
      IF ( size(x3d,3) /= lz ) CALL print_err(crtn, 'depth array do not match data')

      IF ( lct == 1 ) THEN
         PRINT *, ''
         PRINT *, ' --- '//TRIM(crtn)//': creating file '//TRIM(cf_in)//' to write '//TRIM(cv_in)//'!'

         IF ( vflag /= 0.) THEN
            rmin =  1.E6 ; rmax = -1.E6
            DO jk=1, lz
               DO jj=1, nj

                  DO ji=1, ni
                     IF ((x3d(ji,jj,jk) <= rmin).and.(x3d(ji,jj,jk) /= vflag)) rmin = x3d(ji,jj,jk)
                     IF ((x3d(ji,jj,jk) >= rmax).and.(x3d(ji,jj,jk) /= vflag)) rmax = x3d(ji,jj,jk)
                  END DO
               END DO
            END DO
         ELSE
            rmin = minval(x3d) ; rmax = maxval(x3d)
         END IF
         dr = (rmax - rmin)/10.0 ; rmin = rmin - dr ; rmax = rmax + dr
         vextrema(1,:) = (/MINVAL(xlon),MAXVAL(xlon)/); vextrema(2,:) = (/MINVAL(xlat),MAXVAL(xlat)/)
         vextrema(3,:) = (/MINVAL(vtime),MAXVAL(vtime)/)

         !!           CREATE NETCDF OUTPUT FILE :
         CALL sherr( NF90_CREATE(cf_in, NF90_NETCDF4, idx_f),  crtn,cf_in,cv_in)
         !!
         CALL prepare_nc(idx_f, cdt, ni, nj, cv_lo, cv_la, cv_t, vextrema, &
            &            id_x, id_y, id_t, id_lo, id_la, id_tim, crtn,cf_in,cv_in, &
            &            attr_lon=attr_lon, attr_lat=attr_lat, attr_tim=attr_time, &
            &            l_add_valid_min_max=l_add_extrema)
         !! Depth vector:
         IF ( (TRIM(cv_dpth) == 'lev').OR.(TRIM(cv_dpth) == 'depth') ) THEN
            CALL sherr( NF90_DEF_DIM(idx_f, TRIM(cv_dpth), lz, id_z),  crtn,cf_in,cv_in)
         ELSE
            CALL sherr( NF90_DEF_DIM(idx_f, 'z', lz, id_z),  crtn,cf_in,cv_in)
         END IF
         CALL sherr( NF90_DEF_VAR(idx_f, TRIM(cv_dpth), NF90_DOUBLE, id_z,id_dpt, deflate_level=9),  crtn,cf_in,cv_in)
         IF ( lcopy_att_z )  CALL SET_ATTRIBUTES_TO_VAR(idx_f, id_dpt, attr_z, crtn,cf_in,cv_in)
         CALL sherr( NF90_PUT_ATT(idx_f, id_dpt, 'valid_min', MINVAL(vdpth)),  crtn,cf_in,cv_in)
         CALL sherr( NF90_PUT_ATT(idx_f, id_dpt, 'valid_max', MAXVAL(vdpth)),  crtn,cf_in,cv_in)

         !! Variable
         IF ( TRIM(cv_t) /= '' ) THEN
            ALLOCATE (vidim(4))
            vidim = (/id_x,id_y,id_z,id_t/)
         ELSE
            ALLOCATE (vidim(3))
            vidim = (/id_x,id_y,id_z/)
         END IF
         CALL sherr( NF90_DEF_VAR(idx_f, TRIM(cv_in), NF90_FLOAT, vidim, idx_v, deflate_level=9), &
            &       crtn,cf_in,cv_in)
         DEALLOCATE ( vidim )

         !!  VARIABLE ATTRIBUTES
         IF ( lcopy_att_F )  CALL SET_ATTRIBUTES_TO_VAR(idx_f, idx_v, attr_F,  crtn,cf_in,cv_in)
         ! Forcing these attributes:
         IF (vflag/=0.) CALL sherr( NF90_PUT_ATT(idx_f, idx_v,trim(cmv0),           vflag),  crtn,cf_in,cv_in)
         CALL                sherr( NF90_PUT_ATT(idx_f, idx_v,'actual_range', (/rmin,rmax/)),  crtn,cf_in,cv_in)
         CALL                sherr( NF90_PUT_ATT(idx_f, idx_v,'coordinates', &
            &                                TRIM(cv_t)//" "//TRIM(cv_dpth)//" "//TRIM(cv_la)//" "//TRIM(cv_lo)),  crtn,cf_in,cv_in)

         !! Global attributes
         IF ( PRESENT(cextrainfo) ) &
            CALL sherr( NF90_PUT_ATT(idx_f, NF90_GLOBAL, 'Info', TRIM(cextrainfo)),  crtn,cf_in,cv_in)
         CALL sherr( NF90_PUT_ATT(idx_f, NF90_GLOBAL, 'About', trim(cabout)),  crtn,cf_in,cv_in)
         !!
         !!           END OF DEFINITION
         CALL sherr( NF90_ENDDEF(idx_f),  crtn,cf_in,cv_in)
         !!
         !!       Write longitude variable :
         CALL sherr( NF90_PUT_VAR(idx_f, id_lo, xlon),  crtn,cf_in,cv_in)
         !!
         !!       Write latitude variable :
         CALL sherr( NF90_PUT_VAR(idx_f, id_la, xlat),  crtn,cf_in,cv_in)
         !!
         !!       Write depth variable :
         CALL sherr( NF90_PUT_VAR(idx_f, id_dpt, vdpth),  crtn,cf_in,cv_in)
         !!
         !!       Write time variable :
         IF ( TRIM(cv_t) /= '' ) CALL sherr( NF90_PUT_VAR(idx_f, id_tim, vtime),  crtn,cf_in,cv_in)
         !!
      END IF  !IF ( lct == 1 )

      !!                WRITE VARIABLE
      PRINT *, ' --- '//TRIM(crtn)//': writing '//TRIM(cv_in)//' in '//TRIM(cf_in)//', record #', lct
      CALL sherr( NF90_PUT_VAR(idx_f, idx_v,  x3d, start=(/1,1,1,lct/), count=(/ni,nj,lz,1/)),  crtn,cf_in,cv_in)

      !! Sync data from buffer to file
      IF ( lct /= Nt ) CALL sherr( NF90_SYNC (idx_f),  crtn,cf_in,cv_in)
      IF ( lct == Nt ) THEN
         PRINT *, ' --- '//TRIM(crtn)//': closing file '//TRIM(cf_in)//' .'
         CALL sherr( NF90_CLOSE(idx_f),  crtn,cf_in,cv_in)
      END IF
   END SUBROUTINE P3D_T




   SUBROUTINE CHECK_4_MISS(cf_in, cv_in, lmv, rmissv, cmissval)
      !!-
      !! o This routine looks for the presence of a missing value attribute for
      !!   variable "cv_in" in file "cf_in"
      !!          => returns "lmv" (true/false)
      !!          => returns the value in "rmissv" if "lmv==true"
      !!          => returns its attribute name: "cmissval"
      !!
      !! INPUT :
      !! -------
      !!         * cv_in    = variable                                [character]
      !!         * cf_in    = treated file                              [character]
      !!
      !! OUTPUT :
      !! --------
      !!         * imiss    = 0 -> no missing value, 1 -> missing value found [integer]
      !!         * rmissv = value of missing value                             [real]
      !!         * cmissval = name of the missing value attribute           [character]
      !!----------------------------------------------------------------------------
      !!
      INTEGER                       :: id_f, id_v
      CHARACTER(len=*), INTENT(in)  :: cf_in, cv_in
      LOGICAL,          INTENT(out) :: lmv
      REAL(4),          INTENT(out) :: rmissv
      CHARACTER(len=*), INTENT(out) :: cmissval
      !!
      INTEGER :: ierr, jm, ierr1, ierr2
      REAL(4) :: rsf, rao
      CHARACTER(len=80), PARAMETER :: crtn = 'CHECK_4_MISS'
      !!---------------------------------------------------------------------
      CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, id_f),  crtn,cf_in,cv_in) ! Opening file
      CALL sherr( NF90_INQ_VARID(id_f, cv_in, id_v),  crtn,cf_in,cv_in)    ! looking up variable
      ierr1 = NF90_GET_ATT(id_f, id_v, 'scale_factor', rsf)
      ierr2 = NF90_GET_ATT(id_f, id_v, 'add_offset',   rao)

      !! Scanning possible values until found:
      cmissval = '0'
      DO jm=1, nmval
         ierr = NF90_GET_ATT(id_f, id_v, TRIM(c_nm_miss_val(jm)), rmissv)
         IF ( ierr == NF90_NOERR ) THEN
            cmissval = TRIM(c_nm_miss_val(jm))
            EXIT
         END IF
      END DO
      !!
      IF (ierr1 == NF90_NOERR) rmissv = rsf*rmissv
      IF (ierr2 == NF90_NOERR) rmissv = rmissv + rao
      !!
      IF ( ierr == -43 ) THEN
         lmv = .FALSE.
      ELSE
         IF (ierr ==  NF90_NOERR) THEN
            lmv = .TRUE.
            PRINT *, ''
            PRINT *, '  *** CHECK_4_MISS: found missing value attribute '//TRIM(cmissval)//' for '//TRIM(cv_in)//' !'
            PRINT *, '      ( into '//TRIM(cf_in)//')'
            PRINT *, '      => value =', rmissv
            PRINT *, ''
         ELSE
            CALL print_err(crtn, 'problem getting missing_value attribute')
         END IF
      END IF
      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_in,cv_in)
   END SUBROUTINE CHECK_4_MISS



   SUBROUTINE GET_VAR_INFO(cf_in, cv_in, cunit, clnm,  clndr)
      !!
      !! o This routine returns the unit and longname of variable if they exist!
      !!
      !! INPUT :
      !! -------
      !!         * cv_in    = variable                                [character]
      !!         * cf_in    = treated file                            [character]
      !!
      !! OUTPUT :
      !! --------
      !!         * cunit = unit of cv_in                              [character]
      !!         * clnm  = name of the missing value arg.            [character]
      !! OPTIONAL OUTPUT :
      !! -----------------
      !!         * clndr = calendar of cv_in                              [character]
      !!
      !! Author : L. BRODEAU, 2008
      !!
      !!----------------------------------------------------------------------------
      !!
      INTEGER                       :: id_f, id_v
      CHARACTER(len=*), INTENT(in)  :: cf_in, cv_in
      CHARACTER(len=*) , INTENT(out) :: cunit
      CHARACTER(len=*), INTENT(out) :: clnm
      CHARACTER(len=*), INTENT(out), OPTIONAL :: clndr
      !!
      INTEGER :: ierr
      CHARACTER(len=400) :: c00
      CHARACTER(len=80), PARAMETER :: crtn = 'GET_VAR_INFO'
      LOGICAL :: lclndr = .FALSE.
      !!
      IF ( PRESENT(clndr) ) lclndr = .TRUE.
      !!
      !! Opening file :
      CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, id_f),  crtn,cf_in,cv_in)
      !!
      !! Chosing variable :
      CALL sherr( NF90_INQ_VARID(id_f, cv_in, id_v),  crtn,cf_in,cv_in)
      !!
      !!
      c00=''
      ierr = NF90_GET_ATT(id_f, id_v, 'units', c00)
      IF (ierr /= 0) c00 = 'UNKNOWN'
      cunit = TRIM(c00) ;
      !!
      c00=''
      ierr = NF90_GET_ATT(id_f, id_v, 'long_name', c00)
      IF (ierr /= 0) c00 = 'UNKNOWN'
      clnm = TRIM(c00)
      !!
      IF ( lclndr ) THEN
         c00=''
         ierr = NF90_GET_ATT(id_f, id_v, 'calendar', c00)
         IF (ierr /= 0) c00 = 'UNKNOWN'
         clndr = TRIM(c00) ;
      END IF
      !!
      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_in,cv_in)
      !!
   END SUBROUTINE GET_VAR_INFO



   SUBROUTINE DUMP_2D_FIELD(xfld, cf_in, cv_in,   xlon, xlat, cv_lo, cv_la,  rfill)
      !!-----------------------------------------------------------------------------
      !! This routine prints an 2D array into a netcdf file
      !!
      !! INPUT :
      !! -------
      !!        xfld  = 2D array (ni,nj) contening mask            [real4]
      !!        cf_in  = name of the output file                   [character]
      !!        cv_in  = name of the mask variable                 [character]
      !! OPTIONAL :
      !! ----------
      !!        xlon  = longitude array
      !!        xlat  = latitude array
      !!        cv_lo = longitude name
      !!        cv_la = latitude name
      !!        rfill = values of xfld to be masked in netcdf file
      !!
      !!------------------------------------------------------------------------------
      REAL(4),    DIMENSION(:,:), INTENT(in) :: xfld
      CHARACTER(len=*),           INTENT(in) :: cf_in, cv_in
      REAL(8), DIMENSION(:,:), INTENT(in), OPTIONAL :: xlon, xlat
      CHARACTER(len=*),        INTENT(in), OPTIONAL :: cv_lo, cv_la
      REAL(4),                 INTENT(in), OPTIONAL :: rfill
      !!
      INTEGER                                :: id_f, id_v
      INTEGER                                :: id_x, id_y, id_lo, id_la
      INTEGER     :: ni, nj, i01, i02
      LOGICAL :: lzcoord, l_mask
      CHARACTER(len=80), PARAMETER :: crtn = 'DUMP_2D_FIELD'
      ni = size(xfld,1) ; nj = size(xfld,2)
      lzcoord = .FALSE.
      l_mask  = .FALSE.
      IF ( PRESENT(xlon).AND.PRESENT(xlat) ) THEN
         IF ( PRESENT(cv_lo).AND.PRESENT(cv_la) ) THEN
            lzcoord = .TRUE.
            cdt = TEST_XYZ(xlon, xlat, xfld)
         ELSE
            CALL print_err(crtn, 'if you specify xlon and xlat, you must also specify cv_lo and cv_la')
         END IF
         vextrema(1,:) = (/MINVAL(xlon),MAXVAL(xlon)/); vextrema(2,:) = (/MINVAL(xlat),MAXVAL(xlat)/)
      END IF
      IF ( PRESENT(rfill) ) l_mask = .TRUE.
      CALL sherr( NF90_CREATE(cf_in, NF90_NETCDF4, id_f), crtn,cf_in,cv_in)
      IF ( lzcoord ) THEN
         CALL prepare_nc(id_f, cdt, ni, nj, cv_lo, cv_la, '',  vextrema, &
            &            id_x, id_y, i01, id_lo, id_la, i02, &
            &            crtn,cf_in,trim(cv_lo)//'+'//trim(cv_la))
      ELSE
         CALL prepare_nc(id_f, cdt, ni, nj, '', '', '',        vextrema, &
            &            id_x, id_y, i01, id_lo, id_la, i02, &
            &            crtn,cf_in,cv_in)
      END IF
      CALL sherr( NF90_DEF_VAR(id_f, TRIM(cv_in), NF90_FLOAT, (/id_x,id_y/), id_v, deflate_level=9), crtn,cf_in,cv_in)
      IF (l_mask)  CALL sherr( NF90_PUT_ATT(id_f, id_v,trim(cmv0),        rfill                 ), crtn,cf_in,cv_in)
      IF (lzcoord) CALL sherr( NF90_PUT_ATT(id_f, id_v,'coordinates',TRIM(cv_la)//" "//TRIM(cv_lo)), crtn,cf_in,cv_in)
      CALL sherr( NF90_ENDDEF(id_f),  crtn,cf_in,cv_in)
      IF ( lzcoord ) THEN
         CALL sherr( NF90_PUT_VAR(id_f, id_lo, xlon),  crtn,cf_in,cv_lo)
         CALL sherr( NF90_PUT_VAR(id_f, id_la, xlat),  crtn,cf_in,cv_la)
      END IF
      CALL sherr( NF90_PUT_VAR(id_f, id_v, xfld),  crtn,cf_in,cv_in)
      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_in,cv_in)
   END SUBROUTINE DUMP_2D_FIELD


   SUBROUTINE DUMP_3D_FIELD(xfld, cf_in, cv_in,   xlon, xlat, cv_lo, cv_la,  rfill)
      REAL(4),    DIMENSION(:,:,:), INTENT(in) :: xfld
      CHARACTER(len=*),           INTENT(in) :: cf_in, cv_in
      REAL(8), DIMENSION(:,:), INTENT(in), OPTIONAL :: xlon, xlat
      CHARACTER(len=*),        INTENT(in), OPTIONAL :: cv_lo, cv_la
      REAL(4),                 INTENT(in), OPTIONAL :: rfill
      !!
      INTEGER                                :: id_f, id_v
      INTEGER                                :: id_x, id_y, id_z, id_lo, id_la
      INTEGER     :: ni, nj, lz, i01, i02
      LOGICAL :: lzcoord, l_mask
      CHARACTER(len=80), PARAMETER :: crtn = 'DUMP_3D_FIELD'
      !! LOLO: should add levels!!!!
      !!
      ni = size(xfld,1) ; nj = size(xfld,2) ; lz = size(xfld,3)
      lzcoord = .FALSE.
      l_mask  = .FALSE.
      IF ( PRESENT(xlon).AND.PRESENT(xlat) ) THEN
         IF ( PRESENT(cv_lo).AND.PRESENT(cv_la) ) THEN
            lzcoord = .TRUE.
            cdt = TEST_XYZ(xlon, xlat, xfld(:,:,1))
         ELSE
            CALL print_err(crtn, 'if you specify xlon and xlat, you must also specify cv_lo and cv_la')
         END IF
         vextrema(1,:) = (/MINVAL(xlon),MAXVAL(xlon)/); vextrema(2,:) = (/MINVAL(xlat),MAXVAL(xlat)/)
      END IF
      IF ( PRESENT(rfill) ) l_mask = .TRUE.
      CALL sherr( NF90_CREATE(cf_in, NF90_NETCDF4, id_f), crtn,cf_in,cv_in)
      IF ( lzcoord ) THEN
         CALL prepare_nc(id_f, cdt, ni, nj, cv_lo, cv_la, '',  vextrema, &
            &            id_x, id_y, i01, id_lo, id_la, i02, &
            &            crtn,cf_in,trim(cv_lo)//'+'//trim(cv_la))
      ELSE
         CALL prepare_nc(id_f, cdt, ni, nj, '', '', '',        vextrema, &
            &            id_x, id_y, i01, id_lo, id_la, i02, &
            &            crtn,cf_in,cv_in)
      END IF

      CALL sherr( NF90_DEF_DIM(id_f, 'z', lz, id_z),  crtn,cf_in,cv_in)

      CALL sherr( NF90_DEF_VAR(id_f, TRIM(cv_in), NF90_FLOAT, (/id_x,id_y,id_z/), id_v, deflate_level=9), crtn,cf_in,cv_in)
      IF (l_mask)  CALL sherr( NF90_PUT_ATT(id_f, id_v,TRIM(cmv0),        rfill                 ), crtn,cf_in,cv_in)
      IF (lzcoord) CALL sherr( NF90_PUT_ATT(id_f, id_v,'coordinates',TRIM(cv_la)//" "//TRIM(cv_lo)), crtn,cf_in,cv_in)
      CALL sherr( NF90_ENDDEF(id_f),  crtn,cf_in,cv_in)
      IF ( lzcoord ) THEN
         CALL sherr( NF90_PUT_VAR(id_f, id_lo, xlon),  crtn,cf_in,cv_lo)
         CALL sherr( NF90_PUT_VAR(id_f, id_la, xlat),  crtn,cf_in,cv_la)
      END IF
      CALL sherr( NF90_PUT_VAR(id_f, id_v, xfld),  crtn,cf_in,cv_in)
      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_in,cv_in)
   END SUBROUTINE DUMP_3D_FIELD











   SUBROUTINE P2D_MAPPING_AB(cf_out, xlon, xlat, imtrcs, ralfbet, vflag, id_pb,  d2np)
      INTEGER                               :: id_f, id_x, id_y, id_lo, id_la
      CHARACTER(len=*),             INTENT(in) :: cf_out
      REAL(8),    DIMENSION(:,:),   INTENT(in) :: xlon, xlat
      INTEGER(4), DIMENSION(:,:,:), INTENT(in) :: imtrcs
      REAL(8),    DIMENSION(:,:,:), INTENT(in) :: ralfbet
      REAL(8),                      INTENT(in) :: vflag
      INTEGER(2), DIMENSION(:,:),   INTENT(in) :: id_pb
      REAL(4), DIMENSION(:,:), OPTIONAL, INTENT(in) :: d2np ! distance to nearest point (km)

      INTEGER          :: ni, nj, il0, id_n2, id_n3, id_v1, id_v2, id_v3, id_dnp
      LOGICAL          :: l_save_distance_to_np=.FALSE.
      CHARACTER(len=80), PARAMETER :: crtn = 'P2D_MAPPING_AB'

      IF ( PRESENT(d2np) ) l_save_distance_to_np=.TRUE.

      ni = size(ralfbet,1) ; nj = size(ralfbet,2)

      il0 = size(ralfbet,3)
      IF ( il0 /= 2 ) THEN
         PRINT *, 'ralfbet in P2D_MAPPING_AB of io_ezcdf.f90 has wrong shape:', ni, nj, il0
         STOP
      END IF

      il0 = size(imtrcs,3)
      IF ( il0 /= 3 ) THEN
         PRINT *, 'imtrcs in P2D_MAPPING_AB of io_ezcdf.f90 has wrong shape:', ni, nj, il0
         STOP
      END IF


      !!           CREATE NETCDF OUTPUT FILE :
      !!           ---------------------------
      CALL sherr( NF90_CREATE(cf_out, NF90_NETCDF4, id_f),  crtn,cf_out,cdum)

      CALL sherr( NF90_DEF_DIM(id_f, 'x',  ni, id_x), crtn,cf_out,cdum)
      CALL sherr( NF90_DEF_DIM(id_f, 'y',  nj, id_y), crtn,cf_out,cdum)
      CALL sherr( NF90_DEF_DIM(id_f, 'n2',  2, id_n2), crtn,cf_out,cdum)
      CALL sherr( NF90_DEF_DIM(id_f, 'n3',  3, id_n3), crtn,cf_out,cdum)

      CALL sherr( NF90_DEF_VAR(id_f, 'lon',       NF90_DOUBLE, (/id_x,id_y/),       id_lo, deflate_level=9), &
         &        crtn,cf_out,cdum)
      CALL sherr( NF90_DEF_VAR(id_f, 'lat',       NF90_DOUBLE, (/id_x,id_y/),       id_la, deflate_level=9), &
         &        crtn,cf_out,cdum)
      !CALL sherr( NF90_DEF_VAR(id_f, 'metrics',   NF90_INT64,    (/id_x,id_y,id_n3/), id_v1, deflate_level=9), &
      CALL sherr( NF90_DEF_VAR(id_f, 'metrics',   NF90_INT,    (/id_x,id_y,id_n3/), id_v1, deflate_level=9), &
         &        crtn,cf_out,cdum)
      CALL sherr( NF90_DEF_VAR(id_f, 'alphabeta', NF90_DOUBLE, (/id_x,id_y,id_n2/), id_v2, deflate_level=9), &
         &        crtn,cf_out,cdum)
      CALL sherr( NF90_DEF_VAR(id_f, 'iproblem',  NF90_INT,    (/id_x,id_y/),       id_v3, deflate_level=9), &
         &        crtn,cf_out,cdum)
      IF ( l_save_distance_to_np ) THEN
         CALL sherr( NF90_DEF_VAR(id_f, 'dist_np', NF90_REAL, (/id_x,id_y/), id_dnp, deflate_level=9),  crtn,cf_out,cdum)
         CALL sherr( NF90_PUT_ATT(id_f, id_dnp, 'long_name', 'Distance to nearest point'),  crtn,cf_out,cdum)
         CALL sherr( NF90_PUT_ATT(id_f, id_dnp, 'units'    , 'km'                       ),  crtn,cf_out,cdum)
      END IF

      IF ( vflag /= 0. ) THEN
         CALL sherr( NF90_PUT_ATT(id_f, id_v1,TRIM(cmv0),INT(vflag,8)), crtn,cf_out,'metrics (masking)')
         CALL sherr( NF90_PUT_ATT(id_f, id_v2,trim(cmv0),vflag),       crtn,cf_out,'alphabeta (masking)')
      END IF

      CALL sherr( NF90_PUT_ATT(id_f, NF90_GLOBAL, 'Info', 'File containing mapping/weight information for bilinear interpolation with SOSIE.'), &
         &      crtn,cf_out,cdum)
      CALL sherr( NF90_PUT_ATT(id_f, NF90_GLOBAL, 'About', trim(cabout)),  crtn,cf_out,cdum)

      CALL sherr( NF90_ENDDEF(id_f),  crtn,cf_out,cdum) ! END OF DEFINITION

      CALL sherr( NF90_PUT_VAR(id_f, id_lo, xlon),     crtn,cf_out,'lon')
      CALL sherr( NF90_PUT_VAR(id_f, id_la, xlat),     crtn,cf_out,'lat')

      CALL sherr( NF90_PUT_VAR(id_f, id_v1,  imtrcs),  crtn,cf_out,'metrics')
      CALL sherr( NF90_PUT_VAR(id_f, id_v2, ralfbet),  crtn,cf_out,'alphabeta')
      CALL sherr( NF90_PUT_VAR(id_f, id_v3,   id_pb),  crtn,cf_out,'iproblem')

      IF ( l_save_distance_to_np ) CALL sherr( NF90_PUT_VAR(id_f, id_dnp, d2np), crtn,cf_out,'lon' )

      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_out,cdum)

   END SUBROUTINE P2D_MAPPING_AB





   SUBROUTINE  RD_MAPPING_AB(cf_in, imtrcs, ralfbet, id_pb)
      CHARACTER(len=*),          INTENT(in)  :: cf_in
      !INTEGER(8), DIMENSION(:,:,:), INTENT(out) :: imtrcs
      INTEGER(4), DIMENSION(:,:,:), INTENT(out) :: imtrcs
      REAL(8),    DIMENSION(:,:,:), INTENT(out) :: ralfbet
      INTEGER(2), DIMENSION(:,:),   INTENT(out) :: id_pb

      INTEGER :: id_f
      INTEGER :: id_v1, id_v2, id_v3

      CHARACTER(len=80), PARAMETER :: crtn = 'RD_MAPPING_AB'

      CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, id_f),  crtn,cf_in,cdum)

      CALL sherr( NF90_INQ_VARID(id_f, 'metrics',   id_v1),  crtn,cf_in,cdum)
      CALL sherr( NF90_INQ_VARID(id_f, 'alphabeta', id_v2),  crtn,cf_in,cdum)
      CALL sherr( NF90_INQ_VARID(id_f, 'iproblem',  id_v3),  crtn,cf_in,cdum)

      CALL sherr( NF90_GET_VAR(id_f, id_v1, imtrcs),   crtn,cf_in,cdum)
      CALL sherr( NF90_GET_VAR(id_f, id_v2, ralfbet),  crtn,cf_in,cdum)
      CALL sherr( NF90_GET_VAR(id_f, id_v3, id_pb), crtn,cf_in,cdum)

      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_in,cdum)

   END SUBROUTINE RD_MAPPING_AB





   SUBROUTINE PHOVMOLLER(vx, vy, x2d, cf_in, cv_in, cv_x, cv_y, cunit, cln, cuX, cuY)
      !!
      !!-----------------------------------------------------------------------------
      !! This routine prints Hovmoller diagram from a 2D fiels and coordinates given as 2 vectors
      !!
      !! INPUT :
      !! -------
      !!        vx     = X array                                     [real8]
      !!        vy     = Y array                                     [real8]
      !!        x2d    = 2D array (ni,nj) contening mask             [real4]
      !!        cf_in  = name of the output file                   [character]
      !!        cv_in  = name of the mask variable                 [character]
      !!        cf_x   = name of X coordinates                     [character]
      !!        cf_y   = name of Y coordinates                     [character]
      !!        cunit  = unit for treated variable                [character]
      !!        cln = long-name for treated variable              [character]
      !!
      !! OPTIONAL :
      !! ----------
      !!        cuX  = unit X coordinate                          [character]
      !!        cuY  = unit Y coordinate                          [character]
      !!
      !!------------------------------------------------------------------------------
      INTEGER                                :: id_f, id_v
      INTEGER                                :: id_x, id_y, id_lo, id_la
      REAL(8),    DIMENSION(:)  , INTENT(in) :: vx, vy
      REAL(4),    DIMENSION(:,:), INTENT(in) :: x2d
      CHARACTER(len=*),           INTENT(in) :: cf_in, cv_in, cv_x, cv_y, cunit, cln

      CHARACTER(len=*), OPTIONAL, INTENT(in) :: cuX, cuY

      CHARACTER(len=64) :: cuXin, cuYin
      INTEGER :: ni, nj, i01, i02

      CHARACTER(len=80), PARAMETER :: crtn = 'PHOVMOLLER'

      ni = size(x2d,1) ; nj = size(x2d,2)

      cuXin = 'unknown'
      cuYin = 'unknown'
      IF ( present(cuX) ) cuXin = trim(cuX)
      IF ( present(cuY) ) cuYin = trim(cuY)



      vextrema(1,:) = (/minval(vx),maxval(vx)/); vextrema(2,:) = (/minval(vy),maxval(vy)/)

      !!           CREATE NETCDF OUTPUT FILE :
      CALL sherr( NF90_CREATE(CF_IN, NF90_NETCDF4, id_f),  crtn,cf_in,cv_in)

      CALL prepare_nc(id_f, '1d', ni, nj, cv_x, cv_y, '', vextrema, &
         &            id_x, id_y, i01, id_lo, id_la, i02, crtn,cf_in,cv_in)

      CALL sherr( NF90_DEF_VAR(id_f, TRIM(cv_in), NF90_FLOAT, (/id_x,id_y/), id_v, deflate_level=9), crtn,cf_in,cv_in)

      CALL sherr( NF90_PUT_ATT(id_f, id_v, 'long_name', trim(cln)),  crtn,cf_in,cv_in)
      CALL sherr( NF90_PUT_ATT(id_f, id_v, 'units',  trim(cunit) ),  crtn,cf_in,cv_in)
      !! Global attributes
      CALL sherr( NF90_PUT_ATT(id_f, NF90_GLOBAL, 'About', trim(cabout)),  crtn,cf_in,cv_in)

      CALL sherr( NF90_ENDDEF(id_f),  crtn,cf_in,cv_in) ! END OF DEFINITION


      !!          WRITE COORDINATES
      CALL sherr( NF90_PUT_VAR(id_f, id_lo, vx),  crtn,cf_in,cv_in)
      CALL sherr( NF90_PUT_VAR(id_f, id_la, vy),  crtn,cf_in,cv_in)

      !!          WRITE VARIABLE
      CALL sherr( NF90_PUT_VAR(id_f, id_v, x2d),  crtn,cf_in,cv_in)

      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_in,cv_in)

   END SUBROUTINE PHOVMOLLER











   SUBROUTINE GET_SF_AO(cf_in, cv_in, rsf, rao,  italk)
      !!
      !!-----------------------------------------------------------------------
      !! This routine extracts the 'scale_factor' and 'add_offset' of a given
      !! variable from a netcdf file
      !!
      !! INPUT :
      !! -------
      !!          * cf_in      : name of the input file              (character)
      !!          * cv_in      : name of the variable                (character)
      !!
      !! OUTPUT :
      !! --------
      !!          * rsf       : scale factor                        (real)
      !!          * rao       : add offset                          (real)
      !!
      !!  OPTIONAL (INPUT):
      !! ------------------
      !!          * italk     : verbose level: 0 nothing, 1 talks!    (integer)
      !!
      !!------------------------------------------------------------------------
      !!
      INTEGER                      :: id_f, id_v
      CHARACTER(len=*),INTENT(in) :: cf_in, cv_in
      REAL(4),         INTENT(out) :: rsf, rao
      INTEGER, OPTIONAL, INTENT(in) :: italk
      !!
      !! local :
      INTEGER :: ierr1, ierr2, itlk=0
      !!
      CHARACTER(len=80), PARAMETER :: crtn = 'GET_SF_AO'
      !!
      IF ( PRESENT(italk) ) itlk = italk
      !!
      CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, id_f),  crtn,cf_in,cv_in)
      !!
      CALL sherr( NF90_INQ_VARID(id_f, cv_in, id_v),  crtn,cf_in,cv_in)
      !!
      ierr1 = NF90_GET_ATT(id_f, id_v, 'scale_factor', rsf)
      ierr2 = NF90_GET_ATT(id_f, id_v, 'add_offset',   rao)
      !!
      IF ( itlk > 0 ) WRITE(6,*) ' --- GET_SF_AO: variable ', TRIM(cv_in), ' of file ',TRIM(cf_in), ' :'
      !!
      IF ( (ierr1 /= NF90_NOERR).OR.(ierr2 /= NF90_NOERR) ) THEN
         rsf = 1.
         rao = 0.
         IF ( itlk > 0 ) THEN
            WRITE(6,*) '       does not have a "scale_factor" and "add_offset" attributes'
         END IF
      END IF
      !!
      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_in,cv_in)
      !!
      IF ( itlk > 0 ) THEN
         WRITE(6,*) '       => scale_factor =', rsf
         WRITE(6,*) '       => add_offset =', rao
         PRINT *, ''
      END IF
      !!
   END SUBROUTINE GET_SF_AO


   SUBROUTINE sherr(ierr, croutine, ctf, ctv)
      !!
      !! To handle and display error messages
      !!
      INTEGER,            INTENT(in) :: ierr
      !!
      CHARACTER(len=*) , INTENT(in) :: &
         &            ctf,           &    !: treated file
         &            croutine,      &    !: routine name
         &            ctv                 !: treated varible
      !!
      !!
      IF ( ierr /= NF90_NOERR ) THEN
         PRINT *, ''
         WRITE(6,*) '************************************************'
         WRITE(6,*) 'Error occured in procedure ', trim(croutine),' !'
         PRINT *, ''
         WRITE(6,*) 'Treated file     = ', trim(ctf)
         WRITE(6,*) 'Treated variable = ', trim(ctv)
         PRINT *, ''
         WRITE(6,*) '--> aborting program'
         PRINT *, ''
         WRITE(6,*) 'Netcdf message was :'
         WRITE(6,*) trim(NF90_STRERROR(ierr))
         PRINT *, ''
         WRITE(6,*) '************************************************'
         PRINT *, ''
         STOP
      END IF
      !!
   END SUBROUTINE sherr


   FUNCTION TEST_XYZ(rx, ry, rz)
      !!
      !! Testing if 2D coordinates or 1D, and if match shape of data...
      !!
      CHARACTER(len=2) :: TEST_XYZ
      !!
      REAL(8), DIMENSION(:,:), INTENT(in) :: rx, ry
      REAL(4), DIMENSION(:,:), INTENT(in) :: rz
      !!
      INTEGER :: ix1, ix2, iy1, iy2, iz1, iz2
      !!
      ix1 = SIZE(rx,1) ; ix2 = SIZE(rx,2)
      iy1 = SIZE(ry,1) ; iy2 = SIZE(ry,2)
      iz1 = SIZE(rz,1) ; iz2 = SIZE(rz,2)
      !!
      IF ( (ix2 == 1).AND.(iy2 == 1) ) THEN
         !!
         IF ( (ix1 == iz1).AND.(iy1 == iz2) ) THEN
            TEST_XYZ = '1d'
         ELSEIF ( (ix1 == iy1).AND.(ix1 == iz1).AND.(iz2 == 1) ) THEN
            TEST_XYZ = 'y1'
            !! This is thechnicalnj 1D, yet in a 2D shape => [x=nx,y=1]
            !!  => may occur when 2D zonal sections
         ELSE
            PRINT *, 'ERROR, mod_manip.f90 = >TEST_XYZ 1: longitude and latitude array do not match data!'
            PRINT *, '  => ix1,ix2 / iy1,iy2 / iz1,iz2 ='
            PRINT *,      ix1,ix2 ,' /', iy1,iy2 ,' /', iz1,iz2
            PRINT *, ''; STOP
         END IF
         !!
      ELSE
         IF ( (ix1 == iz1).AND.(iy1 == iz1).AND.(ix2 == iz2).AND.(iy2 == iz2) ) THEN
            TEST_XYZ = '2d'
         ELSE
            PRINT *, 'ERROR, mod_manip.f90 = >TEST_XYZ 2: longitude and latitude array do not match data!'
            PRINT *, ''; STOP
         END IF
      END IF
      !!
   END FUNCTION TEST_XYZ


   !   SUBROUTINE TEST_XYZ(rx, ry, rd, cdm)
   !
   !      !! Testing if 2D coordinates or 1D, and if match shape of data...
   !
   !      REAL(8), DIMENSION(:,:), INTENT(in)  :: rx, ry
   !      REAL(4), DIMENSION(:,:), INTENT(in)  :: rd
   !      CHARACTER(len=2)       , INTENT(out) :: cdm
   !
   !      INTEGER :: ix1, ix2, iy1, iy2, id1, id2
   !
   !      ix1 = size(rx,1) ; ix2 = size(rx,2)
   !      iy1 = size(ry,1) ; iy2 = size(ry,2)
   !      id1 = size(rd,1) ; id2 = size(rd,2)
   !
   !      IF ( (ix2 == 1).AND.(iy2 == 1) ) THEN
   !
   !         IF ( (ix1 == id1).AND.(iy1 == id2) ) THEN
   !            cdm = '1d'
   !         ELSE
   !            CALL print_err('cdm', 'longitude and latitude array do not match data (1d)')
   !         END IF
   !
   !      ELSE!
   !
   !         IF ( (ix1 == id1).AND.(iy1 == id1).AND.(ix2 == id2).AND.(iy2 == id2) ) THEN
   !            cdm = '2d'
   !         ELSE
   !            CALL print_err('cdm', 'longitude and latitude array do not match data (2d)')
   !         END IF
   !
   !      END IF
   !
   !   END SUBROUTINE TEST_XYZ





   SUBROUTINE prepare_nc(id_file, cdt0, nx, ny, cv_lon, cv_lat, cv_time, vxtrm,     &
      &                  id_ji, id_jj, id_jt, id_lon, id_lat, id_time, cri,cfi,cvi, &
      &                  attr_lon, attr_lat, attr_tim, l_add_valid_min_max)
      !!----------------------------------------------------------------------------------
      !!----------------------------------------------------------------------------------
      INTEGER,                 INTENT(in)  :: id_file, nx, ny
      CHARACTER(len=2),        INTENT(in)  :: cdt0
      CHARACTER(len=*),        INTENT(in)  :: cv_lon, cv_lat, cv_time, cri,cfi,cvi
      REAL(8), DIMENSION(3,2), INTENT(in)  :: vxtrm
      INTEGER,                 INTENT(out) :: id_ji, id_jj, id_jt, id_lon, id_lat, id_time

      TYPE(var_attr), DIMENSION(nbatt_max), OPTIONAL, INTENT(in) :: attr_lon, attr_lat, attr_tim
      LOGICAL, OPTIONAL, INTENT(in)             :: l_add_valid_min_max

      LOGICAL ::  &
         &       lcopy_att_lon = .FALSE., &
         &       lcopy_att_lat = .FALSE., &
         &       lcopy_att_tim = .FALSE., &
         &       l_add_extrema = .TRUE.

      IF ( PRESENT(attr_lon).AND.(attr_lon(1)%itype>0) ) lcopy_att_lon = .TRUE.
      IF ( PRESENT(attr_lat).AND.(attr_lat(1)%itype>0) ) lcopy_att_lat = .TRUE.
      IF ( PRESENT(attr_tim).AND.(attr_tim(1)%itype>0) ) lcopy_att_tim = .TRUE.

      IF ( PRESENT(l_add_valid_min_max) ) l_add_extrema = l_add_valid_min_max


      !!    HORIZONTAL
      IF ( (TRIM(cv_lon) /= '').AND.(TRIM(cv_lat) /= '') ) THEN
         !!
         IF ( (cdt0 == '2d').OR.(cdt0 == 'y1') ) THEN
            CALL sherr( NF90_DEF_DIM(id_file, 'x', nx, id_ji), cri,cfi,cvi)
            CALL sherr( NF90_DEF_DIM(id_file, 'y', ny, id_jj), cri,cfi,cvi)
            CALL sherr( NF90_DEF_VAR(id_file, TRIM(cv_lon), NF90_DOUBLE, (/id_ji,id_jj/), id_lon, deflate_level=9), cri,cfi,cvi)
            CALL sherr( NF90_DEF_VAR(id_file, TRIM(cv_lat), NF90_DOUBLE, (/id_ji,id_jj/), id_lat, deflate_level=9), cri,cfi,cvi)
            !!
         ELSE IF ( cdt0 == '1d' ) THEN
            CALL sherr( NF90_DEF_DIM(id_file, TRIM(cv_lon), nx, id_ji), cri,cfi,cvi)
            CALL sherr( NF90_DEF_DIM(id_file, TRIM(cv_lat), ny, id_jj), cri,cfi,cvi)
            CALL sherr( NF90_DEF_VAR(id_file, TRIM(cv_lon), NF90_DOUBLE, id_ji, id_lon, deflate_level=9), cri,cfi,cvi)
            CALL sherr( NF90_DEF_VAR(id_file, TRIM(cv_lat), NF90_DOUBLE, id_jj, id_lat, deflate_level=9), cri,cfi,cvi)
         END IF
         !!
         IF ( lcopy_att_lon )  CALL SET_ATTRIBUTES_TO_VAR(id_file, id_lon, attr_lon,  cri,cfi,cvi)
         IF ( l_add_extrema ) THEN
            CALL sherr( NF90_PUT_ATT(id_file, id_lon, 'valid_min', vxtrm(1,1)), cri,cfi,cvi)
            CALL sherr( NF90_PUT_ATT(id_file, id_lon, 'valid_max', vxtrm(1,2)), cri,cfi,cvi)
         END IF
         !!
         IF ( lcopy_att_lat )  CALL SET_ATTRIBUTES_TO_VAR(id_file, id_lat, attr_lat,  cri,cfi,cvi)
         IF ( l_add_extrema ) THEN
            CALL sherr( NF90_PUT_ATT(id_file, id_lat, 'valid_min', vxtrm(2,1)), cri,cfi,cvi)
            CALL sherr( NF90_PUT_ATT(id_file, id_lat, 'valid_max', vxtrm(2,2)), cri,cfi,cvi)
         END IF
         !!
      ELSE
         CALL sherr( NF90_DEF_DIM(id_file, 'x', nx, id_ji), cri,cfi,cvi)
         CALL sherr( NF90_DEF_DIM(id_file, 'y', ny, id_jj), cri,cfi,cvi)
      END IF
      !!
      !!  TIME
      IF ( TRIM(cv_time) /= '' ) THEN
         CALL sherr( NF90_DEF_DIM(id_file, TRIM(cv_time), NF90_UNLIMITED, id_jt), cri,cfi,cvi)
         CALL sherr( NF90_DEF_VAR(id_file, TRIM(cv_time), NF90_DOUBLE, id_jt, id_time, deflate_level=9), cri,cfi,cvi)
         IF ( lcopy_att_tim )  CALL SET_ATTRIBUTES_TO_VAR(id_file, id_time, attr_tim,  cri,cfi,cvi)
         IF ( l_add_extrema ) THEN
            CALL sherr( NF90_PUT_ATT(id_file, id_time, 'valid_min',vxtrm(3,1)), cri,cfi,cvi)
            CALL sherr( NF90_PUT_ATT(id_file, id_time, 'valid_max',vxtrm(3,2)), cri,cfi,cvi)
         END IF
      END IF
      !!
   END SUBROUTINE prepare_nc


   SUBROUTINE SET_ATTRIBUTES_TO_VAR(idx_f, idx_v, vattr,  cri,cfi,cvi)
      !!
      INTEGER,                              INTENT(in) :: idx_f, idx_v
      TYPE(var_attr), DIMENSION(nbatt_max), INTENT(in) :: vattr
      CHARACTER(len=*),                     INTENT(in) :: cri,cfi,cvi
      !!
      INTEGER :: jat, il
      CHARACTER(len=128) :: cat
      !!
      DO jat = 1, nbatt_max
         cat = vattr(jat)%cname
         !PRINT *, 'LOLO: will set '//TRIM(cat)//'!'
         IF ( TRIM(cat) == 'null' ) EXIT
         IF ( (TRIM(cat)/='grid_type').AND.(TRIM(cat)/=trim(cmv0)).AND.(TRIM(cat)/='missing_value') &
            & .AND.(TRIM(cat)/='scale_factor').AND.(TRIM(cat)/='add_offset') ) THEN
            IF ( vattr(jat)%itype == 2 ) THEN
               CALL sherr( NF90_PUT_ATT(idx_f, idx_v, TRIM(cat), TRIM(vattr(jat)%val_char)), cri,cfi,cvi)
            ELSE
               il = vattr(jat)%ilength
               !PRINT *, 'LOLO: will set '//TRIM(cat)//' with:', vattr(jat)%val_num(:il)
               CALL sherr( NF90_PUT_ATT(idx_f, idx_v, TRIM(cat), vattr(jat)%val_num(:il)), cri,cfi,cvi)
            END IF
         END IF
      END DO
   END SUBROUTINE SET_ATTRIBUTES_TO_VAR



   FUNCTION GET_TIME_UNIT_T0( cstr )
      TYPE(t_unit_t0)              :: GET_TIME_UNIT_T0
      CHARACTER(len=*), INTENT(in) :: cstr  ! ex: "days since 1950-01-01 00:00:00" => 'd' 1950 1 1 0 0 0
      !!
      INTEGER :: i1, i2, is, ncday
      CHARACTER(len=32) :: cdum, chour, cday
      !!
      CHARACTER(len=80), PARAMETER :: crtn = 'GET_TIME_UNIT_T0'
      i1 = SCAN(cstr, ' ')
      i2 = SCAN(cstr, ' ', back=.TRUE.)
      !!
      chour = cstr(i2+1:)
      cdum =  cstr(1:i1-1)
      SELECT CASE(TRIM(cdum))
      CASE('days')
         GET_TIME_UNIT_T0%unit = 'd'
      CASE('hours')
         GET_TIME_UNIT_T0%unit = 'h'
      CASE('seconds')
         GET_TIME_UNIT_T0%unit = 's'
      CASE DEFAULT
         CALL print_err(crtn, 'the only time units we know are "seconds", "hours" and "days"')
      END SELECT
      !!
      cdum = cstr(1:i2-1)   ! => "days since 1950-01-01"
      i2 = SCAN(TRIM(cdum), ' ', back=.TRUE.)
      IF ( cstr(i1+1:i2-1) /= 'since' ) STOP 'Aborting GET_TIME_UNIT_T0!'
      !!
      !! Date:
      cdum = cstr(i2+1:)
      i1 = SCAN(cdum, '-')
      i2 = SCAN(cdum, '-', back=.TRUE.)

      is = SCAN(TRIM(cdum), ' ', back=.TRUE.)

      !! Day of calendar:
      cday = cdum(:is)
      ncday = len(TRIM(cday))
      IF ( i1 == 5 ) THEN
         READ(cdum(1:4),'(i4)')  GET_TIME_UNIT_T0%year
         IF ( i2 == 8 ) THEN
            READ(cdum(6:7),'(i2)')  GET_TIME_UNIT_T0%month
            IF ( ncday == 10 ) READ(cdum(9:10),'(i2)') GET_TIME_UNIT_T0%day
            IF ( ncday ==  9 ) READ(cdum(9:9) ,'(i2)') GET_TIME_UNIT_T0%day
         ELSEIF ( i2 == 7 ) THEN
            READ(cdum(6:6),'(i2)')  GET_TIME_UNIT_T0%month
            IF ( ncday == 8 ) READ(cdum(8:8),'(i2)') GET_TIME_UNIT_T0%day
            IF ( ncday == 9 ) READ(cdum(8:9),'(i2)') GET_TIME_UNIT_T0%day
         ELSE
            CALL print_err(crtn, 'origin date not recognized, must be "yyyy-mm-dd"')
         END IF
      ELSE
         CALL print_err(crtn, 'origin date not recognized, must be "yyyy-mm-dd"')
      END IF

      !! Time of day:
      i1 = SCAN(chour, ':')
      i2 = SCAN(chour, ':', back=.TRUE.)
      IF ( (i1 /= 3).OR.(i2 /= 6) ) CALL print_err(crtn, 'hour of origin date not recognized, must be "hh:mm:ss"')
      READ(chour(1:2),'(i2)')  GET_TIME_UNIT_T0%hour
      READ(chour(4:5),'(i2)')  GET_TIME_UNIT_T0%minute
      READ(chour(7:8),'(i2)') GET_TIME_UNIT_T0%second
      !!
   END FUNCTION GET_TIME_UNIT_T0

   FUNCTION L_IS_LEAP_YEAR( iyear )
      LOGICAL :: L_IS_LEAP_YEAR
      INTEGER, INTENT(in) :: iyear
      L_IS_LEAP_YEAR = .FALSE.
      !IF ( MOD(iyear,4)==0 )     L_IS_LEAP_YEAR = .TRUE.
      !IF ( MOD(iyear,100)==0 )   L_IS_LEAP_YEAR = .FALSE.
      !IF ( MOD(iyear,400)==0 )   L_IS_LEAP_YEAR = .TRUE.
      IF ( (MOD(iyear,4)==0).AND.(.NOT.((MOD(iyear,100)==0).AND.(MOD(iyear,400)/=0)) ) )  L_IS_LEAP_YEAR = .TRUE.
   END FUNCTION L_IS_LEAP_YEAR


   FUNCTION nbd_m( imonth, iyear )
      !! Number of days in month # imonth of year iyear
      INTEGER :: nbd_m
      INTEGER, INTENT(in) :: imonth, iyear
      IF (( imonth > 12 ).OR.( imonth < 1 ) ) THEN
         PRINT *, 'ERROR: nbd_m of io_ezcdf.f90 => 1 <= imonth cannot <= 12 !!!', imonth
         STOP
      END IF
      nbd_m = tdmn(imonth)
      IF ( L_IS_LEAP_YEAR(iyear) )  nbd_m = tdml(imonth)
   END FUNCTION nbd_m



   FUNCTION time_to_date( cal_unit_ref0, rt,  date_prev )
      !!
      !! Converts a time like "hours since 19XX" to "19XX-MM-DD-hh-mm-ss"
      !!
      TYPE(date)                  :: time_to_date
      TYPE(t_unit_t0), INTENT(in) :: cal_unit_ref0 ! date of the origin of the calendar ex: "'d',1950,1,1,0,0,0" for "days since 1950-01-01
      REAL(8)        , INTENT(in) :: rt ! time as specified as cal_unit_ref0
      TYPE(date),      INTENT(in), OPTIONAL :: date_prev 
      !!
      TYPE(date) :: date_start   ! date to start the itteration from...      
      REAL(8)    :: zt, zinc
      INTEGER    :: jy, jmn, jd, jm, jh, jd_old, js
      LOGICAL    :: l_control_start, lcontinue
      REAL(8)    :: rjs, rjs_t !, rjs_t_old
      CHARACTER(len=80), PARAMETER :: crtn = 'time_to_date'
      !!
      !! Takes values from cal_unit_ref0 !
      date_start%year   = cal_unit_ref0%year
      date_start%month  = cal_unit_ref0%month
      date_start%day    = cal_unit_ref0%day
      date_start%hour   = cal_unit_ref0%hour
      date_start%minute = cal_unit_ref0%minute
      date_start%second = cal_unit_ref0%second


      l_control_start = .FALSE.
      IF ( PRESENT(date_prev) ) THEN
         !! Only if reasonable value:
         IF ( date_prev%year >= cal_unit_ref0%year ) THEN
            l_control_start = .TRUE.
            date_start = date_prev
         END IF
      END IF
      
      zinc = 60. ! increment in seconds!
      !!
      jy  = date_start%year
      jmn = date_start%month
      jd  = date_start%day
      jh  = date_start%hour
      jm  = date_start%minute
      js  = date_start%second
      !!
      rjs = REAL(js, 8)

      IF ( l_control_start ) THEN
         rjs_t = rjs_t_old
      ELSE
         rjs_t = 0.
         rjs_t_old = 0.
      END IF
      !!
      
      IF ( rt > 0. ) THEN  !(if ==0, then date = to cal_unit_ref0 !)

         !! Converting rt to seconds
         SELECT CASE(cal_unit_ref0%unit)
         CASE('d')
            !PRINT *, ' Switching from days to seconds!'
            zt = rt*24.*3600.
         CASE('h')
            !PRINT *, ' Switching from hours to seconds!'
            zt = rt*3600.
         CASE('s')
            zt = rt
         CASE DEFAULT
            CALL print_err(crtn, 'the only time units we know are "s" (seconds), "h" (hours) and "d" (days)')
         END SELECT

         !WRITE(*,'(" *** start: ",i4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2," s cum =",i," d cum =",i)') jy, jmn, jd, jh, jm, js, NINT(rjs_t)

         lcontinue = .TRUE.
         DO WHILE ( lcontinue )
            jd_old = jd
            rjs_t  = rjs_t + zinc
            rjs    = rjs + zinc

            IF ( MOD(rjs_t,60.) == 0. ) THEN
               rjs = 0.
               jm  = jm+1
            END IF
            IF ( jm == 60 ) THEN
               jm = 0
               jh = jh+1
            END IF
            !IF ( MOD(rjs_t,3600) == 0 ) THEN
            !   jm = 0
            !   jh = jh+1
            !END IF
            !
            IF ( jh == 24 ) THEN
               jh = 0
               jd = jd + 1
            END IF
            IF ( jd == nbd_m(jmn,jy)+1 ) THEN
               jd = 1
               jmn = jmn + 1
            END IF
            IF ( jmn == 13 ) THEN
               jmn  = 1
               jy = jy+1
            END IF
            !
            !IF ( jd /= jd_old ) THEN
            !   WRITE(*,'(" ***  now : ",i4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2," s cum =",i," d cum =",i)') jy, jmn, jd, jh, jm, js,  rjs_t
            !END IF
            IF ( (zt <= rjs_t).AND.(zt > rjs_t_old) ) lcontinue = .FALSE.
            IF ( jy == 2021 ) THEN
               PRINT *, 'rjs_t =', rjs_t
               STOP 'ERROR: time_to_date => beyond 2020!'
            END IF
            rjs_t_old = rjs_t
         END DO
         !
      END IF
      
      time_to_date%year   = jy
      time_to_date%month  = jmn
      time_to_date%day    = jd
      time_to_date%hour   = jh
      time_to_date%minute = jm
      time_to_date%second = NINT(rjs)
      !
      !WRITE(*,'(" *** time_to_date => Date : ",i4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2)') jy, jmn, jd, jh, jm, NINT(rjs)
      !STOP'LOLO io_ezcdf.f90'
      !WRITE(6,*) 'LOLO io_ezcdf.f90: rjs, rjs_t =', rjs, rjs_t
      !STOP
      
   END FUNCTION time_to_date



   FUNCTION to_epoch_time_scalar( cal_unit_ref0, rt, dt )
      !!
      !! Tests: wcstools:
      !! => getdate fd2ts   2013-03-31T23:38:38  (to "seconds since 1950)
      !! => getdate fd2tsu  2013-03-31T23:38:38  (to "seconds since 1970" aka epoch unix...)
      !!
      INTEGER(8)                  :: to_epoch_time_scalar
      TYPE(t_unit_t0), INTENT(in) :: cal_unit_ref0 ! date of the origin of the calendar ex: "'d',1950,1,1,0,0,0" for "days since 1950-01-01
      REAL(8)        , INTENT(in) :: rt ! time as specified as cal_unit_ref0
      REAL(8)  , OPTIONAL      , INTENT(in) :: dt ! time as specified as cal_unit_ref0
      !!
      REAL(8)    :: zt, zinc
      INTEGER    :: jy, jmn, jd, jm, jh, jd_old, ipass, nb_pass, js
      LOGICAL    :: lcontinue
      REAL(8)    :: rjs_t, rjs0_epoch, rjs !, rjs_t_old
      CHARACTER(len=80), PARAMETER :: crtn = 'to_epoch_time_scalar'
      !!
      nb_pass = 1
      IF ( PRESENT(dt) ) THEN
         IF ( dt < 60. ) nb_pass = 2
      END IF
      !!
      SELECT CASE(cal_unit_ref0%unit)
      CASE('d')
         PRINT *, ' Switching from days to seconds!'
         zt = rt*24.*3600.
      CASE('h')
         PRINT *, ' Switching from hours to seconds!'
         zt = rt*3600.
      CASE('s')
         zt = rt
      CASE DEFAULT
         CALL print_err(crtn, 'the only time units we know are "s" (seconds), "h" (hours) and "d" (days)')
      END SELECT
      !!
      !!
      !! Starting with large time increment (in seconds):
      zinc = 60. ! increment in seconds!
      !!
      jy = cal_unit_ref0%year
      jmn = cal_unit_ref0%month
      jd = cal_unit_ref0%day
      js= cal_unit_ref0%second
      jm= cal_unit_ref0%minute
      jh= cal_unit_ref0%hour
      !!
      rjs = REAL(js, 8)
      rjs_t = 0. ; rjs_t_old = 0.
      !!
      DO ipass=1, nb_pass

         PRINT *, '' ; PRINT *, ' ipass = ', ipass

         IF ( ipass == 2 ) THEN
            !!
            PRINT *, ' after 1st pass:', jy, jmn, jd, jh, jm, rjs
            !! Tiny increment (sometime the time step is lower than 1 second on stelite ephem tracks!!!)
            !! Rewind 2 minutes backward:
            CALL REWIND_2MN(jy, jmn, jd, jh, jm)
            zinc = 0.1 ! seconds
            rjs_t     = rjs_t     - 120.
            rjs_t_old = rjs_t - zinc
            PRINT *, ' before 2nd pass:', jy, jmn, jd, jh, jm, rjs
         END IF


         !!
         !WRITE(*,'(" *** start: ",i4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2," s cum =",i," d cum =",i)') jy, jmn, jd, jh, jm, js,  js_t
         lcontinue = .TRUE.
         DO WHILE ( lcontinue )
            jd_old = jd
            rjs_t = rjs_t + zinc
            rjs = rjs + zinc
            !!

            IF ( ipass == 1 ) THEN
               IF ( MOD(rjs_t,60.) == 0. ) THEN
                  rjs = 0.
                  jm  = jm+1
               END IF
            ELSE
               IF ( rjs >= 60.) THEN
                  rjs = rjs - 60.
                  jm  = jm+1
               END IF
            END IF
            IF ( jm == 60 ) THEN
               jm = 0
               jh = jh+1
            END IF
            !IF ( MOD(rjs_t,3600) == 0 ) THEN
            !   jm = 0
            !   jh = jh+1
            !END IF
            !
            IF ( jh == 24 ) THEN
               jh = 0
               jd = jd + 1
            END IF
            IF ( jd == nbd_m(jmn,jy)+1 ) THEN
               jd = 1
               jmn = jmn + 1
            END IF
            IF ( jmn == 13 ) THEN
               jmn  = 1
               jy = jy+1
            END IF
            IF ( (jy==1970).AND.(jmn==1).AND.(jd==1).AND.(jh==0).AND.(jm==0).AND.(rjs==0.) ) rjs0_epoch = rjs_t
            !
            !IF ( jd /= jd_old ) THEN
            !   WRITE(*,'(" ***  now : ",i4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2," s cum =",i," d cum =",i)') jy, jmn, jd, jh, jm, js,  rjs_t
            !END IF
            IF ( (zt <= rjs_t).AND.(zt > rjs_t_old) ) lcontinue = .FALSE.
            IF ( jy == 2019 ) THEN
               PRINT *, 'rjs_t =', rjs_t
               STOP 'ERROR: to_epoch_time_scalar => beyond 2018!'
            END IF
            rjs_t_old = rjs_t
         END DO
         !
         to_epoch_time_scalar = rjs_t - rjs0_epoch
         !
         !PRINT *, 'Found !!!', zt, REAL(rjs_t)
         WRITE(*,'(" *** to_epoch_time_scalar => Date : ",i4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2)') jy, jmn, jd, jh, jm, NINT(rjs)
         !PRINT *, ' + rjs0_epoch =', rjs0_epoch
         !PRINT *, '  Date (epoch) =>', to_epoch_time_scalar

      END DO

   END FUNCTION to_epoch_time_scalar


   SUBROUTINE to_epoch_time_vect( cal_unit_ref0, vt,  l_dt_below_sec )
      !!
      TYPE(t_unit_t0), INTENT(in) :: cal_unit_ref0 ! date of the origin of the calendar ex: "'d',1950,1,1,0,0,0" for "days since 1950-01-01
      REAL(8)        , DIMENSION(:), INTENT(inout) :: vt           ! time as specified as cal_unit_ref0
      LOGICAL        , OPTIONAL, INTENT(in) :: l_dt_below_sec
      !!
      REAL(8)    :: zt !, dt_min
      REAL(8), DIMENSION(:), ALLOCATABLE :: vtmp

      INTEGER    :: ntr, jt, jd_old
      INTEGER    :: jy, jmn, jd, jh, jm, js, jx
      LOGICAL    :: lcontinue, l_be_accurate
      REAL(8)    :: rjs_t, rjs0_epoch, zinc, rjs !, rjs_t_old
      CHARACTER(len=80), PARAMETER :: crtn = 'to_epoch_time_vect'
      !!
      l_be_accurate = .FALSE.
      IF ( PRESENT(l_dt_below_sec) ) l_be_accurate = l_dt_below_sec
      !!
      SELECT CASE(cal_unit_ref0%unit)
      CASE('d')
         PRINT *, ' Switching from days to seconds!'
         vt = 24.*3600.*vt
      CASE('h')
         PRINT *, ' Switching from hours to seconds!'
         vt = 3600.*vt
      CASE('s')
         PRINT *, ' Already in seconds...'
      CASE DEFAULT
         CALL print_err(crtn, 'the only time units we know are "s" (seconds), "h" (hours) and "d" (days)')
      END SELECT
      !!
      IF ( l_be_accurate ) PRINT *, '    high accuracy turned on!'
      !!

      ntr = SIZE(vt,1)
      !PRINT *, ' ntr =', ntr
      !ALLOCATE ( vtmp(ntr-1) )
      !vtmp(:) = vt(2:ntr) - vt(1:ntr-1)
      !dt_min = MINVAL(vtmp)
      !PRINT *, ' * Minimum time-step (in seconds) => ', dt_min
      !dt_min = dt_min - dt_min/100.
      !PRINT *, vtmp(1:40)
      !DEALLOCATE ( vtmp )


      ALLOCATE ( vtmp(ntr) )
      !!
      !! Starting with large time increment (in seconds):

      !!
      jy = cal_unit_ref0%year
      jmn = cal_unit_ref0%month
      jd  = cal_unit_ref0%day
      jh  = cal_unit_ref0%hour
      jm  = cal_unit_ref0%minute
      js  = cal_unit_ref0%second
      jx  = 0
      !!
      !!
      rjs = REAL(js, 8)
      rjs_t = 0. ; rjs_t_old = 0.
      !!
      zinc = 60. ! seconds
      !!
      DO jt = 0, ntr ! 0 is the initial pass to find the start

         zt   = vt(MAX(jt,1))

         !PRINT *, '' ; PRINT *, ' jt, zt = ', jt, zt

         IF ( (l_be_accurate).AND.(jt > 0) ) zinc = 0.01 ! seconds

         !PRINT *, ' after previous jt:', jy, jmn, jd, jh, jm, js, jx
         IF ( jt == 1 ) THEN
            !! Rewind 2 minutes backward (to find again, more acuratenj the start date (accuracy of 1/100 of a second rather than 60 second!
            CALL REWIND_2MN(jy, jmn, jd, jh, jm)
            rjs_t     = rjs_t - 120.
            rjs_t_old = rjs_t - zinc
         END IF
         !PRINT *, ' before comming jt:', jy, jmn, jd, jh, jm, js, jx


         !WRITE(*,'(" *** start: ",i4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2," s cum =",i," d cum =",i)') jy, jmn, jd, jh, jm, js,  js_t
         lcontinue = .TRUE.
         DO WHILE ( lcontinue )
            jd_old = jd
            rjs_t = rjs_t + zinc

            !IF ( jt > 0 ) PRINT *, ' *** LOOP before zinc: jh, jm, js, jx =>', jh, jm, js, jx

            js = js + INT(zinc)
            IF ( l_be_accurate ) jx = jx + NINT(100.*MOD(zinc, 1.0))

            IF ( jt == 0 ) THEN
               IF ( MOD(rjs_t,60.) == 0. ) THEN
                  js = 0
                  jm = jm+1
               END IF
            ELSE
               IF ( jx >= 100) THEN
                  js = js + jx/100
                  jx = jx - 100
               END IF
               IF ( js >= 60 ) THEN
                  js  = js - 60
                  jm  = jm+1
               END IF
            END IF
            IF ( jm == 60 ) THEN
               jm = 0
               jh = jh+1
            END IF
            IF ( jh == 24 ) THEN
               jh = 0
               jd = jd + 1
            END IF
            IF ( jd == nbd_m(jmn,jy)+1 ) THEN
               jd = 1
               jmn = jmn + 1
            END IF
            IF ( jmn == 13 ) THEN
               jmn  = 1
               jy = jy+1
            END IF
            !IF ( jt > 0 ) PRINT *, ' *** LOOP after zinc: jh, jm, js, jx =>', jh, jm, js, jx
            !
            IF ( (jy==1970).AND.(jmn==1).AND.(jd==1).AND.(jh==0).AND.(jm==0).AND.(rjs==0.) ) rjs0_epoch = rjs_t
            !
            IF ( (zt <= rjs_t).AND.(zt > rjs_t_old) ) lcontinue = .FALSE.
            IF ( jy == 2019 ) THEN
               PRINT *, 'rjs_t =', rjs_t
               STOP 'ERROR: to_epoch_time_vect => beyond 2018!'
            END IF
            rjs_t_old = rjs_t
         END DO
         !
         vtmp(MAX(jt,1)) = rjs_t - rjs0_epoch
         !
         IF ( jt==1   ) WRITE(*,'(" *** to_epoch_time_vect => Start date : ",i4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2,":",i2.2," epoch: ",f15.4)') jy, jmn, jd, jh, jm, js,jx, vtmp(MAX(jt,1))
         IF ( jt==ntr ) WRITE(*,'(" *** to_epoch_time_vect =>   End date : ",i4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2,":",i2.2," epoch: ",f15.4)') jy, jmn, jd, jh, jm, js,jx, vtmp(MAX(jt,1))
         !
      END DO

      vt(:) = vtmp(:)

      DEALLOCATE ( vtmp )

   END SUBROUTINE to_epoch_time_vect



   SUBROUTINE REWIND_2MN(jy, jmn, jd, jh, jm)
      INTEGER, INTENT(inout) :: jy, jmn, jd, jh, jm
      !! jm:
      !! 3 => 1
      !! 2 => 0
      !! 1 => 59
      !! 0 => 58
      IF ( jm > 1 ) THEN
         jm = jm - 2
      ELSE
         jm = 58 + jm
         IF ( jh == 0 ) THEN
            jh = 23
            IF ( jd == 1 ) THEN
               jd = nbd_m(jmn-1,jy)
               IF ( jmn == 1 ) THEN
                  jmn = 12
                  jy  = jy - 1
               ELSE ! jmn
                  jmn = jmn - 1
               END IF
            ELSE ! jd
               jd = jd - 1
            END IF
         ELSE ! jh
            jh = jh - 1
         END IF
      END IF

   END SUBROUTINE REWIND_2MN









   SUBROUTINE print_err(crout, cmess,  ivect)
      CHARACTER(len=*),      INTENT(in) :: crout, cmess
      INTEGER, DIMENSION(:), INTENT(in), OPTIONAL :: ivect
      PRINT *, ''
      WRITE(6,*) 'ERROR in ',TRIM(crout),' (io_ezcdf.f90): '
      WRITE(6,*) TRIM(cmess)
      IF( PRESENT(ivect) ) WRITE(6,*) INT(ivect(:),2)
      PRINT *, ''
      STOP
   END SUBROUTINE print_err









   SUBROUTINE coordinates_from_var_attr( cf_in, cv_in, nb_coor, coord_list )
      !!
      !! Tests: wcstools:
      !! => getdate fd2ts   2013-03-31T23:38:38  (to "seconds since 1950)
      !! => getdate fd2tsu  2013-03-31T23:38:38  (to "seconds since 1970" aka epoch unix...)
      !!
      CHARACTER(len=*),                 INTENT(in)  :: cf_in, cv_in
      INTEGER,                          INTENT(out) :: nb_coor
      CHARACTER(len=128), DIMENSION(4), INTENT(out) :: coord_list
      !!
      INTEGER :: nl, ja, jp, jl, jlp, Natt
      TYPE(var_attr), DIMENSION(nbatt_max) :: v_att
      CHARACTER(len=128) :: csc

      nb_coor       = 0   ! 0 => no "coordinates" attribute !
      coord_list(:) = '0'

      CALL GETVAR_ATTRIBUTES(cf_in, cv_in,  Natt, v_att)
      PRINT *, '  => attributes of '//TRIM(cv_in)//' are:', v_att(:Natt)

      !! try to know coordinates the variable depends on from attribute coordinates:
      DO ja = 1, Natt

         IF ( TRIM(v_att(ja)%cname) == 'coordinates' ) THEN
            nb_coor = 1
            csc = TRIM(v_att(ja)%val_char)
            PRINT *, ' *** we found "coordinates" attribute for "'//TRIM(cv_in)//'":'
            PRINT *, '   => ', TRIM(csc); PRINT *, ''
            nl = LEN(TRIM(csc)) ; PRINT *, nl ; ! length of the string!!

            jp=0 ; jlp=1
            DO jl = 2, nl
               IF(csc(jl:jl)==' ') THEN
                  jp=jp+1
                  PRINT *, jl
                  coord_list(jp) = csc(jlp:jl)
                  jlp = jl + 1 ! jump the ' '
                  nb_coor = nb_coor + 1
               END IF
            END DO
            coord_list(jp+1) = csc(jlp:nl) ! the last one (had no ' ' at the end...)
            PRINT *, ''
            PRINT *, ' *** Based on argument "coordinates", "'//TRIM(cv_in)//'" depends on these:'
            DO jl=1, nb_coor
               PRINT *, ' - ', TRIM(coord_list(jl))
            END DO
            PRINT *, ''

         END IF
      END DO


   END SUBROUTINE coordinates_from_var_attr



END MODULE io_ezcdf
