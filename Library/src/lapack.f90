!============================================================================
! This is a collection of LAPACK routines that BATSRUS uses and independent
! of the real precision. 
! You are encouraged to use the local LAPACK library if available.
!
! subroutines: xerbla, lsame
!
! function:    ilaenv
! 
!============================================================================
      subroutine XERBLA( SRNAME, INFO )

!
!  -- LAPACK auxiliary routine (version 1.1) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      character*6        SRNAME
      integer            INFO
!     ..
!
!  Purpose
!  =======
!
!  XERBLA  is an error handler for the LAPACK routines.
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!  =========
!
!  SRNAME  (input) character*6
!          The name of the routine which called XERBLA.
!
!  INFO    (input) integer
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
!     .. Executable Statements ..
!
      WRITE( *, FMT = 9999 )SRNAME, INFO
!
      call CON_STOP_EXT('LAPACK::XERBLA')
!
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ', &
           'an illegal value' )
!
!     End of XERBLA
!
      end
!============================================================================
      logical          FUNCTION LSAME( CA, CB )
!
!  -- LAPACK auxiliary routine (version 1.1) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      character          CA, CB
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) character*1
!  CB      (input) character*1
!          CA and CB specify the single characters to be compared.
!
!     ..
!     .. Local Scalars ..
      integer            INTA, INTB, ZCODE
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
      LSAME = CA.eq.CB
      if( LSAME ) &
        RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR( 'Z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
!
      if( ZCODE.eq.90 .OR. ZCODE.eq.122 ) then
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
         if( INTA.GE.97 .and. INTA.LE.122 ) INTA = INTA - 32
         if( INTB.GE.97 .and. INTB.LE.122 ) INTB = INTB - 32
!
      else if( ZCODE.eq.233 .OR. ZCODE.eq.169 ) then
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
         if( INTA.GE.129 .and. INTA.LE.137 .OR. &
            INTA.GE.145 .and. INTA.LE.153 .OR. &
            INTA.GE.162 .and. INTA.LE.169 ) INTA = INTA + 64
         if( INTB.GE.129 .and. INTB.LE.137 .OR. &
            INTB.GE.145 .and. INTB.LE.153 .OR. &
            INTB.GE.162 .and. INTB.LE.169 ) INTB = INTB + 64
!
      else if( ZCODE.eq.218 .OR. ZCODE.eq.250 ) then
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
         if( INTA.GE.225 .and. INTA.LE.250 ) INTA = INTA - 32
         if( INTB.GE.225 .and. INTB.LE.250 ) INTB = INTB - 32
      end if
      LSAME = INTA.eq.INTB
!
!     RETURN
!
!     End of LSAME
!
      end
!=============================================================================
      integer FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      character*( * )    NAME, OPTS
      integer            ISPEC, N1, N2, N3, N4
!     ..
!
!  Purpose
!  =======
!
!  ILAENV is called from the LAPACK routines to choose problem-dependent
!  parameters for the local environment.  See ISPEC for a description of
!  the parameters.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPEC   (input) integer
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR and QZ methods
!               for nonsymmetric eigenvalue problems.
!
!  NAME    (input) character*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) character*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) integer
!  N2      (input) integer
!  N3      (input) integer
!  N4      (input) integer
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
! (ILAENV) (output) integer
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      if( NB.LE.1 ) NB = MAX( 1, N )
!
!  =====================================================================
!
!     .. Local Scalars ..
      logical            CNAME, SNAME
      character*1        C1
      character*2        C2, C4
      character*3        C3
      character*6        SUBNAM
      integer            I, IC, IZ, NB, NBMIN, NX
!     ..

!     .. Executable Statements ..
!
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
!
!     Invalid value for ISPEC
!
      ILAENV = -1
      RETURN
!
  100 CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      if( IZ.eq.90 .OR. IZ.eq.122 ) then
!
!        ASCII character set
!
         if( IC.GE.97 .and. IC.LE.122 ) then
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            do I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               if( IC.GE.97 .and. IC.LE.122 ) &
                 SUBNAM( I:I ) = CHAR( IC-32 )
            end do
         end if
!
      else if( IZ.eq.233 .OR. IZ.eq.169 ) then
!
!        EBCDIC character set
!
         if( ( IC.GE.129 .and. IC.LE.137 ) .OR. &
            ( IC.GE.145 .and. IC.LE.153 ) .OR. &
            ( IC.GE.162 .and. IC.LE.169 ) ) then
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            do I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               if( ( IC.GE.129 .and. IC.LE.137 ) .OR. &
                  ( IC.GE.145 .and. IC.LE.153 ) .OR. &
                  ( IC.GE.162 .and. IC.LE.169 ) ) &
                 SUBNAM( I:I ) = CHAR( IC+64 )
               end do 
         end if
!
      else if( IZ.eq.218 .OR. IZ.eq.250 ) then
!
!        Prime machines:  ASCII+128
!
         if( IC.GE.225 .and. IC.LE.250 ) then
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            do I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               if( IC.GE.225 .and. IC.LE.250 ) &
                 SUBNAM( I:I ) = CHAR( IC-32 )
            end do
         end if
      end if
!
      C1 = SUBNAM( 1:1 )
      SNAME = C1.eq.'S' .OR. C1.eq.'D'
      CNAME = C1.eq.'C' .OR. C1.eq.'Z'
      if( .not.( CNAME .OR. SNAME ) ) &
        RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
!
      GO TO ( 110, 200, 300 ) ISPEC
!
  110 CONTINUE
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
      NB = 1
!
      if( C2.eq.'GE' ) then
         if( C3.eq.'TRF' ) then
            if( SNAME ) then
               NB = 64
            else
               NB = 64
            end if
         else if( C3.eq.'QRF' .OR. C3.eq.'RQF' .OR. C3.eq.'LQF' .OR. &
                 C3.eq.'QLF' ) then
            if( SNAME ) then
               NB = 32
            else
               NB = 32
            end if
         else if( C3.eq.'HRD' ) then
            if( SNAME ) then
               NB = 32
            else
               NB = 32
            end if
         else if( C3.eq.'BRD' ) then
            if( SNAME ) then
               NB = 32
            else
               NB = 32
            end if
         else if( C3.eq.'TRI' ) then
            if( SNAME ) then
               NB = 64
            else
               NB = 64
            end if
         end if
      else if( C2.eq.'PO' ) then
         if( C3.eq.'TRF' ) then
            if( SNAME ) then
               NB = 64
            else
               NB = 64
            end if
         end if
      else if( C2.eq.'SY' ) then
         if( C3.eq.'TRF' ) then
            if( SNAME ) then
               NB = 64
            else
               NB = 64
            end if
         else if( SNAME .and. C3.eq.'TRD' ) then
            NB = 1
         else if( SNAME .and. C3.eq.'GST' ) then
            NB = 64
         end if
      else if( CNAME .and. C2.eq.'HE' ) then
         if( C3.eq.'TRF' ) then
            NB = 64
         else if( C3.eq.'TRD' ) then
            NB = 1
         else if( C3.eq.'GST' ) then
            NB = 64
         end if
      else if( SNAME .and. C2.eq.'OR' ) then
         if( C3( 1:1 ).eq.'G' ) then
            if( C4.eq.'QR' .OR. C4.eq.'RQ' .OR. C4.eq.'LQ' .OR. &
               C4.eq.'QL' .OR. C4.eq.'HR' .OR. C4.eq.'TR' .OR. &
               C4.eq.'BR' ) then
               NB = 32
            end if
         else if( C3( 1:1 ).eq.'M' ) then
            if( C4.eq.'QR' .OR. C4.eq.'RQ' .OR. C4.eq.'LQ' .OR. &
               C4.eq.'QL' .OR. C4.eq.'HR' .OR. C4.eq.'TR' .OR. &
               C4.eq.'BR' ) then
               NB = 32
            end if
         end if
      else if( CNAME .and. C2.eq.'UN' ) then
         if( C3( 1:1 ).eq.'G' ) then
            if( C4.eq.'QR' .OR. C4.eq.'RQ' .OR. C4.eq.'LQ' .OR. &
               C4.eq.'QL' .OR. C4.eq.'HR' .OR. C4.eq.'TR' .OR. &
               C4.eq.'BR' ) then
               NB = 32
            end if
         else if( C3( 1:1 ).eq.'M' ) then
            if( C4.eq.'QR' .OR. C4.eq.'RQ' .OR. C4.eq.'LQ' .OR. &
               C4.eq.'QL' .OR. C4.eq.'HR' .OR. C4.eq.'TR' .OR. &
               C4.eq.'BR' ) then
               NB = 32
            end if
         end if
      else if( C2.eq.'GB' ) then
         if( C3.eq.'TRF' ) then
            if( SNAME ) then
               if( N4.LE.64 ) then
                  NB = 1
               else
                  NB = 32
               end if
            else
               if( N4.LE.64 ) then
                  NB = 1
               else
                  NB = 32
               end if
            end if
         end if
      else if( C2.eq.'PB' ) then
         if( C3.eq.'TRF' ) then
            if( SNAME ) then
               if( N2.LE.64 ) then
                  NB = 1
               else
                  NB = 32
               end if
            else
               if( N2.LE.64 ) then
                  NB = 1
               else
                  NB = 32
               end if
            end if
         end if
      else if( C2.eq.'TR' ) then
         if( C3.eq.'TRI' ) then
            if( SNAME ) then
               NB = 64
            else
               NB = 64
            end if
         end if
      else if( C2.eq.'LA' ) then
         if( C3.eq.'UUM' ) then
            if( SNAME ) then
               NB = 64
            else
               NB = 64
            end if
         end if
      else if( SNAME .and. C2.eq.'ST' ) then
         if( C3.eq.'EBZ' ) then
            NB = 1
         end if
      end if
      ILAENV = NB
      RETURN
!
  200 CONTINUE
!
!     ISPEC = 2:  minimum block size
!
      NBMIN = 2
      if( C2.eq.'GE' ) then
         if( C3.eq.'QRF' .OR. C3.eq.'RQF' .OR. C3.eq.'LQF' .OR. &
            C3.eq.'QLF' ) then
            if( SNAME ) then
               NBMIN = 2
            else
               NBMIN = 2
            end if
         else if( C3.eq.'HRD' ) then
            if( SNAME ) then
               NBMIN = 2
            else
               NBMIN = 2
            end if
         else if( C3.eq.'BRD' ) then
            if( SNAME ) then
               NBMIN = 2
            else
               NBMIN = 2
            end if
         else if( C3.eq.'TRI' ) then
            if( SNAME ) then
               NBMIN = 2
            else
               NBMIN = 2
            end if
         end if
      else if( C2.eq.'SY' ) then
         if( C3.eq.'TRF' ) then
            if( SNAME ) then
               NBMIN = 8
            else
               NBMIN = 8
            end if
         else if( SNAME .and. C3.eq.'TRD' ) then
            NBMIN = 2
         end if
      else if( CNAME .and. C2.eq.'HE' ) then
         if( C3.eq.'TRD' ) then
            NBMIN = 2
         end if
      else if( SNAME .and. C2.eq.'OR' ) then
         if( C3( 1:1 ).eq.'G' ) then
            if( C4.eq.'QR' .OR. C4.eq.'RQ' .OR. C4.eq.'LQ' .OR. &
               C4.eq.'QL' .OR. C4.eq.'HR' .OR. C4.eq.'TR' .OR. &
               C4.eq.'BR' ) then
               NBMIN = 2
            end if
         else if( C3( 1:1 ).eq.'M' ) then
            if( C4.eq.'QR' .OR. C4.eq.'RQ' .OR. C4.eq.'LQ' .OR. &
               C4.eq.'QL' .OR. C4.eq.'HR' .OR. C4.eq.'TR' .OR. &
               C4.eq.'BR' ) then
               NBMIN = 2
            end if
         end if
      else if( CNAME .and. C2.eq.'UN' ) then
         if( C3( 1:1 ).eq.'G' ) then
            if( C4.eq.'QR' .OR. C4.eq.'RQ' .OR. C4.eq.'LQ' .OR. &
               C4.eq.'QL' .OR. C4.eq.'HR' .OR. C4.eq.'TR' .OR. &
               C4.eq.'BR' ) then
               NBMIN = 2
            end if
         else if( C3( 1:1 ).eq.'M' ) then
            if( C4.eq.'QR' .OR. C4.eq.'RQ' .OR. C4.eq.'LQ' .OR. &
               C4.eq.'QL' .OR. C4.eq.'HR' .OR. C4.eq.'TR' .OR. &
               C4.eq.'BR' ) then
               NBMIN = 2
            end if
         end if
      end if
      ILAENV = NBMIN
      RETURN
!
  300 CONTINUE
!
!     ISPEC = 3:  crossover point
!
      NX = 0
      if( C2.eq.'GE' ) then
         if( C3.eq.'QRF' .OR. C3.eq.'RQF' .OR. C3.eq.'LQF' .OR. &
            C3.eq.'QLF' ) then
            if( SNAME ) then
               NX = 128
            else
               NX = 128
            end if
         else if( C3.eq.'HRD' ) then
            if( SNAME ) then
               NX = 128
            else
               NX = 128
            end if
         else if( C3.eq.'BRD' ) then
            if( SNAME ) then
               NX = 128
            else
               NX = 128
            end if
         end if
      else if( C2.eq.'SY' ) then
         if( SNAME .and. C3.eq.'TRD' ) then
            NX = 1
         end if
      else if( CNAME .and. C2.eq.'HE' ) then
         if( C3.eq.'TRD' ) then
            NX = 1
         end if
      else if( SNAME .and. C2.eq.'OR' ) then
         if( C3( 1:1 ).eq.'G' ) then
            if( C4.eq.'QR' .OR. C4.eq.'RQ' .OR. C4.eq.'LQ' .OR. &
               C4.eq.'QL' .OR. C4.eq.'HR' .OR. C4.eq.'TR' .OR. &
               C4.eq.'BR' ) then
               NX = 128
            end if
         end if
      else if( CNAME .and. C2.eq.'UN' ) then
         if( C3( 1:1 ).eq.'G' ) then
            if( C4.eq.'QR' .OR. C4.eq.'RQ' .OR. C4.eq.'LQ' .OR. &
               C4.eq.'QL' .OR. C4.eq.'HR' .OR. C4.eq.'TR' .OR. &
               C4.eq.'BR' ) then
               NX = 128
            end if
         end if
      end if
      ILAENV = NX
      RETURN
!
  400 CONTINUE
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
      ILAENV = 6
      RETURN
!
  500 CONTINUE
!
!     ISPEC = 5:  minimum column dimension (not used)
!
      ILAENV = 2
      RETURN
!
  600 CONTINUE 
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
      ILAENV = INT( real( MIN( N1, N2 ) )*1.6E0 )
      RETURN
!
  700 CONTINUE
!
!     ISPEC = 7:  number of processors (not used)
!
      ILAENV = 1
      RETURN
!
  800 CONTINUE
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
      ILAENV = 50
      RETURN
!
!     End of ILAENV
!
      end
!============================================================================
!============================================================================
! This is a collection of single precision LAPACK routines that BATSRUS uses. 
! You are encouraged to use the local LAPACK library if available.
!
! subroutines: sgetrf, sgetrs, sgetf2, slaswp
! 
!============================================================================
!============================================================================
      subroutine SGETRF( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993 
!
!     .. Scalar Arguments ..
      integer            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      integer            IPIV( * )
      real*4             A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  SGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) integer
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) integer
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) real*4 array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) integer
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) integer array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) integer
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
!  =====================================================================
!
!     .. parameters ..
      real*4             ONE
      parameter          ( ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      integer            I, IINFO, J, JB, NB
!     ..
!     .. external subroutines ..
      external           SGEMM, SGETF2, SLASWP, STRSM, XERBLA
!     ..
!     .. external Functions ..
      integer            ILAENV
      external           ILAENV
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      if( M.lt.0 ) then
         INFO = -1
      else if( N.lt.0 ) then
         INFO = -2
      else if( LDA.lt.MAX( 1, M ) ) then
         INFO = -4
      end if
      if( INFO.ne.0 ) then
         CALL XERBLA( 'SGETRF', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if( M.eq.0 .OR. N.eq.0 ) &
        RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'SGETRF', ' ', M, N, -1, -1 )
      if( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) then
!
!        Use unblocked code.
!
         CALL SGETF2( M, N, A, LDA, IPIV, INFO )
      else
!
!        Use blocked code.
!
         do J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
            CALL SGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
            if( INFO.eq.0 .and. IINFO.GT.0 ) &
              INFO = IINFO + J - 1
            do I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
            end do
!
!           Apply interchanges to columns 1:J-1.
!
            CALL SLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!
            if( J+JB.LE.N ) then
!
!              Apply interchanges to columns J+JB:N.
!
               CALL SLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1, &
                           IPIV, 1 )
!
!              Compute block row of U.
!
               CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
                          N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), &
                          LDA )
               if( J+JB.LE.M ) then
!
!                 Update trailing submatrix.
!
                  CALL SGEMM( 'No transpose', 'No transpose', M-J-JB+1, &
                             N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA, &
                             A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ), &
                             LDA )
               end if
            end if
         end do
      end if
      RETURN
!
!     End of SGETRF
!
      end
!===========================================================================
      subroutine SGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993 
!
!     .. Scalar Arguments ..
      character          TRANS
      integer            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      integer            IPIV( * )
      real*4             A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  SGETRS solves a system of linear equations
!     A * X = B  or  A' * X = B
!  with a general N-by-N matrix A using the LU factorization computed
!  by SGETRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) character*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A'* X = B  (Transpose)
!          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!  N       (input) integer
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) integer
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) real*4 array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by SGETRF.
!
!  LDA     (input) integer
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) integer array, dimension (N)
!          The pivot indices from SGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) real*4 array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) integer
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) integer
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. parameters ..
      real*4             ONE
      parameter          ( ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      logical            NOTRAN
!     ..
!     .. external Functions ..
      logical            LSAME
      external           LSAME
!     ..
!     .. external subroutines ..
      external           SLASWP, STRSM, XERBLA
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      if( .not.NOTRAN .and. .not.LSAME( TRANS, 'T' ) .and. .not. &
         LSAME( TRANS, 'C' ) ) then
         INFO = -1
      else if( N.lt.0 ) then
         INFO = -2
      else if( NRHS.lt.0 ) then
         INFO = -3
      else if( LDA.lt.MAX( 1, N ) ) then
         INFO = -5
      else if( LDB.lt.MAX( 1, N ) ) then
         INFO = -8
      end if
      if( INFO.ne.0 ) then
         CALL XERBLA( 'SGETRS', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if( N.eq.0 .OR. NRHS.eq.0 ) &
        RETURN
!
      if( NOTRAN ) then
!
!        Solve A * X = B.
!
!        Apply row interchanges to the right hand sides.
!
         CALL SLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
!
!        Solve L*X = B, overwriting B with X.
!
         CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS, &
                    ONE, A, LDA, B, LDB )
!
!        Solve U*X = B, overwriting B with X.
!
         CALL STRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N, &
                    NRHS, ONE, A, LDA, B, LDB )
      else
!
!        Solve A' * X = B.
!
!        Solve U'*X = B, overwriting B with X.
!
         CALL STRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, &
                    ONE, A, LDA, B, LDB )
!
!        Solve L'*X = B, overwriting B with X.
!
         CALL STRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE, &
                    A, LDA, B, LDB )
!
!        Apply row interchanges to the solution vectors.
!
         CALL SLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      end if
!
      RETURN
!
!     End of SGETRS
!
      end
!===========================================================================
      subroutine SGETF2( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1992
!
!     .. Scalar Arguments ..
      integer            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      integer            IPIV( * )
      real*4             A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  SGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) integer
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) integer
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) real*4 array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) integer
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) integer array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) integer
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  =====================================================================
!
!     .. parameters ..
      real*4             ONE, ZERO
      parameter          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      integer            J, JP
!     ..
!     .. external Functions ..
      integer            ISAMAX
      external           ISAMAX
!     ..
!     .. external subroutines ..
      external           SGER, SSCAL, SSWAP, XERBLA
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      if( M.lt.0 ) then
         INFO = -1
      else if( N.lt.0 ) then
         INFO = -2
      else if( LDA.lt.MAX( 1, M ) ) then
         INFO = -4
      end if
      if( INFO.ne.0 ) then
         CALL XERBLA( 'SGETF2', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if( M.eq.0 .OR. N.eq.0 ) &
        RETURN
!
      do J = 1, MIN( M, N )
!
!        Find pivot and test for singularity.
!
         JP = J - 1 + ISAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         if( A( JP, J ).ne.ZERO ) then
!
!           Apply the interchange to columns 1:N.
!
            if( JP.ne.J ) &
              CALL SSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
!
!           Compute elements J+1:M of J-th column.
!
            if( J.lt.M ) &
              CALL SSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
!
         else if( INFO.eq.0 ) then
!
            INFO = J
         end if
!
         if( J.lt.MIN( M, N ) ) then
!
!           Update trailing submatrix.
!
            CALL SGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA, &
                      A( J+1, J+1 ), LDA )
         end if
      end do
      RETURN
!
!     End of SGETF2
!
      end
!===========================================================================

      subroutine SLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      integer            INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
      integer            IPIV( * )
      real*4             A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  SLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) integer
!          The number of columns of the matrix A.
!
!  A       (input/output) real*4 array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) integer
!          The leading dimension of the array A.
!
!  K1      (input) integer
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) integer
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) integer array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) integer
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
!  Further Details
!  ===============
!
!  Modified by
!   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!
! =====================================================================
!
!     .. Local Scalars ..
      integer            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      real*4             TEMP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
      if( INCX.GT.0 ) then
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      else if( INCX.lt.0 ) then
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      else
         RETURN
      end if
!
      N32 = ( N / 32 )*32
      if( N32.ne.0 ) then
         do J = 1, N32, 32
            IX = IX0
            do I = I1, I2, INC
               IP = IPIV( IX )
               if( IP.ne.I ) then
                  do K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
                     end do 
               end if
               IX = IX + INCX
            end do
         end do
      end if
      if( N32.ne.N ) then
         N32 = N32 + 1
         IX = IX0
         do I = I1, I2, INC
            IP = IPIV( IX )
            if( IP.ne.I ) then
               do K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
                  end do 
            end if
            IX = IX + INCX
            end do 
      end if
!
      RETURN
!
!     End of SLASWP
!
      end
!=============================================================================
!=============================================================================
! This is a collection of real*8 LAPACK routines that BATSRUS uses. 
! You are encouraged to use the local LAPACK library if available.
!
! subroutines: dgetrf, dgetrs, dgetf2, dlaswp
! 
!=============================================================================
!=============================================================================
      subroutine DGETRF( M, N, A, LDA, IPIV, INFO )

!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      integer            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      integer            IPIV( * )
      real*8         A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) integer
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) integer
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) real*8 array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) integer
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) integer array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) integer
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
!  =====================================================================
!
!     .. parameters ..
      real*8         ONE
      parameter          ( ONE = 1 )
!     ..
!     .. Local Scalars ..
      integer            I, IINFO, J, JB, NB
!     ..
!     .. external subroutines ..
      external           XERBLA, DGEMM, DGETF2, DLASWP, DTRSM
!     ..
!     .. external Functions ..
      integer            ILAENV
      external           ILAENV
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      if( M.lt.0 ) then
         INFO = -1
      else if( N.lt.0 ) then
         INFO = -2
      else if( LDA.lt.MAX( 1, M ) ) then
         INFO = -4
      end if
      if( INFO.ne.0 ) then
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if( M.eq.0 .OR. N.eq.0 ) &
        RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      if( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) then
!
!        Use unblocked code.
!
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      else
!
!        Use blocked code.
!
         do J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
            if( INFO.eq.0 .and. IINFO.GT.0 ) &
              INFO = IINFO + J - 1
            do I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
            end do
!
!           Apply interchanges to columns 1:J-1.
!
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!
            if( J+JB.LE.N ) then
!
!              Apply interchanges to columns J+JB:N.
!
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1, &
                           IPIV, 1 )
!
!              Compute block row of U.
!
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
                          N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), &
                          LDA )
               if( J+JB.LE.M ) then
!
!                 Update trailing submatrix.
!
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1, &
                             N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA, &
                             A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ), &
                             LDA )
               end if
            end if
         end do
      end if
      if (INFO.GT.0) then
        PRINT *,'LAPACK routine DGETRF:'
        PRINT *,'U(',INFO,INFO,') is exactly zero. The matrix'
        PRINT *,'is singular: the inverse cannot be computed.'
        call CON_STOP_EXT('LAPACK::DGETRF')
      endif
      RETURN
!
!     End of DGETRF
!
      end
!=============================================================================
      subroutine DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      character          TRANS
      integer            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      integer            IPIV( * )
      real*8         A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETRS solves a system of linear equations
!     A * X = B,  A**T * X = B,  or  A**H * X = B
!  with a general N-by-N matrix A using the LU factorization computed
!  by DGETRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) character*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B     (No transpose)
!          = 'T':  A**T * X = B  (Transpose)
!          = 'C':  A**H * X = B  (Conjugate transpose)
!
!  N       (input) integer
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) integer
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) real*8 array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by DGETRF.
!
!  LDA     (input) integer
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) integer array, dimension (N)
!          The pivot indices from DGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) real*8 array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) integer
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) integer
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. parameters ..
      real*8         ONE
      parameter          ( ONE = 1 )
!     ..
!     .. Local Scalars ..
      logical            NOTRAN
!     ..
!     .. external Functions ..
      logical            LSAME
      external           LSAME
!     ..
!     .. external subroutines ..
      external           XERBLA, DLASWP, DTRSM
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      if( .not.NOTRAN .and. .not.LSAME( TRANS, 'T' ) .and. .not. &
         LSAME( TRANS, 'C' ) ) then
         INFO = -1
      else if( N.lt.0 ) then
         INFO = -2
      else if( NRHS.lt.0 ) then
         INFO = -3
      else if( LDA.lt.MAX( 1, N ) ) then
         INFO = -5
      else if( LDB.lt.MAX( 1, N ) ) then
         INFO = -8
      end if
      if( INFO.ne.0 ) then
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if( N.eq.0 .OR. NRHS.eq.0 ) &
        RETURN
!
      if( NOTRAN ) then
!
!        Solve A * X = B.
!
!        Apply row interchanges to the right hand sides.
!
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
!
!        Solve L*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS, &
                    ONE, A, LDA, B, LDB )
!
!        Solve U*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N, &
                    NRHS, ONE, A, LDA, B, LDB )
      else
!
!        Solve A**T * X = B  or A**H * X = B.
!
!        Solve U'*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Upper', TRANS, 'Non-unit', N, NRHS, ONE, &
                    A, LDA, B, LDB )
!
!        Solve L'*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Lower', TRANS, 'Unit', N, NRHS, ONE, A, &
                    LDA, B, LDB )
!
!        Apply row interchanges to the solution vectors.
!
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      end if
!
      RETURN
!
!     End of DGETRS
!
      end
!=============================================================================
      subroutine DGETF2( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      integer            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      integer            IPIV( * )
      real*8         A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) integer
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) integer
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) real*8 array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) integer
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) integer array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) integer
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  =====================================================================
!
!     .. parameters ..
      real*8         ONE, ZERO
      parameter          ( ONE = 1, ZERO = 0 )
!     ..
!     .. Local Scalars ..
      integer            J, JP
!     ..
!     .. external Functions ..
      integer            IDAMAX
      external           IDAMAX
!     ..
!     .. external subroutines ..
      external           XERBLA, DGER, DSCAL, DSWAP
!     ..

!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      if( M.lt.0 ) then
         INFO = -1
      else if( N.lt.0 ) then
         INFO = -2
      else if( LDA.lt.MAX( 1, M ) ) then
         INFO = -4
      end if
      if( INFO.ne.0 ) then
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if( M.eq.0 .OR. N.eq.0 ) &
        RETURN
!
      do J = 1, MIN( M, N )
!
!        Find pivot and test for singularity.
!
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         if( A( JP, J ).ne.ZERO ) then
!
!           Apply the interchange to columns 1:N.
!
            if( JP.ne.J ) &
              CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
!
!           Compute elements J+1:M of J-th column.
!
            if( J.lt.M ) &
              CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
!
         else if( INFO.eq.0 ) then
!
            INFO = J
         end if
!
         if( J.lt.MIN( M, N ) ) then
!
!           Update trailing submatrix.
!
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), &
                       LDA, A( J+1, J+1 ), LDA )
         end if
      end do
      RETURN
!
!     End of DGETF2
!
      end
!=============================================================================
      subroutine DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      integer            INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
      integer            IPIV( * )
      real*8         A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) integer
!          The number of columns of the matrix A.
!
!  A       (input/output) real*8 array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) integer
!          The leading dimension of the array A.
!
!  K1      (input) integer
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) integer
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) integer array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) integer
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
! =====================================================================
!
!     .. Local Scalars ..
      integer            I, IP, IX
!     ..
!     .. external subroutines ..
      external           DSWAP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
      if( INCX.eq.0 ) &
        RETURN
      if( INCX.GT.0 ) then
         IX = K1
      else
         IX = 1 + ( 1-K2 )*INCX
      end if
      if( INCX.eq.1 ) then
         do I = K1, K2
            IP = IPIV( I )
            if( IP.ne.I ) &
              CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
         end do
      else if( INCX.GT.1 ) then
         do I = K1, K2
            IP = IPIV( IX )
            if( IP.ne.I ) &
              CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
            end do 
      else if( INCX.lt.0 ) then
         do I = K2, K1, -1
            IP = IPIV( IX )
            if( IP.ne.I ) &
              CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
         end do
      end if
!
      RETURN
!
!     End of DLASWP
!
      end
!============================================================================
