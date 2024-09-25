!============================================================================
! This is a collection of single precision BLAS routines that BATSRUS uses. 
! You are encouraged to use the local BLAS library if available.
!
! Subroutines: scopy, sgemm, sgemv, sger,  sscal, sswap, strsm
!
! Functions:   isamax
!=============================================================================
subroutine scopy(n,sx,incx,sy,incy)
  !
  !     copies a vector, x, to a vector, y.
  !     uses unrolled loops for increments equal to 1.
  !     jack dongarra, linpack, 3/11/78.
  !     modified 12/3/93, array(1) declarations changed to array(*)
  !
  real*4 sx(*),sy(*)
  integer i,incx,incy,ix,iy,m,mp1,n
  !
  if(n.le.0)return

  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do i =  1,n
     sy(iy) = sx(ix)
     ix = ix + incx
     iy = iy + incy
  end do
  return
end subroutine scopy
!=============================================================================
SUBROUTINE SGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
     BETA, C, LDC )
  !     .. Scalar Arguments ..
  CHARACTER*1        TRANSA, TRANSB
  INTEGER            M, N, K, LDA, LDB, LDC
  REAL*4             ALPHA, BETA
  !     .. Array Arguments ..
  REAL*4             A( LDA, * ), B( LDB, * ), C( LDC, * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  SGEMM  performs one of the matrix-matrix operations
  !
  !     C := alpha*op( A )*op( B ) + beta*C,
  !
  !  where  op( X ) is one of
  !
  !     op( X ) = X   or   op( X ) = X',
  !
  !  alpha and beta are scalars, and A, B and C are matrices, with op( A )
  !  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
  !
  !  Parameters
  !  ==========
  !
  !  TRANSA - CHARACTER*1.
  !           On entry, TRANSA specifies the form of op( A ) to be used in
  !           the matrix multiplication as follows:
  !
  !              TRANSA = 'N' or 'n',  op( A ) = A.
  !
  !              TRANSA = 'T' or 't',  op( A ) = A'.
  !
  !              TRANSA = 'C' or 'c',  op( A ) = A'.
  !
  !           Unchanged on exit.
  !
  !  TRANSB - CHARACTER*1.
  !           On entry, TRANSB specifies the form of op( B ) to be used in
  !           the matrix multiplication as follows:
  !
  !              TRANSB = 'N' or 'n',  op( B ) = B.
  !
  !              TRANSB = 'T' or 't',  op( B ) = B'.
  !
  !              TRANSB = 'C' or 'c',  op( B ) = B'.
  !
  !           Unchanged on exit.
  !
  !  M      - INTEGER.
  !           On entry,  M  specifies  the number  of rows  of the  matrix
  !           op( A )  and of the  matrix  C.  M  must  be at least  zero.
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry,  N  specifies the number  of columns of the matrix
  !           op( B ) and the number of columns of the matrix C. N must be
  !           at least zero.
  !           Unchanged on exit.
  !
  !  K      - INTEGER.
  !           On entry,  K  specifies  the number of columns of the matrix
  !           op( A ) and the number of rows of the matrix op( B ). K must
  !           be at least  zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - REAL*4          .
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  A      - REAL*4           array of DIMENSION ( LDA, ka ), where ka is
  !           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
  !           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
  !           part of the array  A  must contain the matrix  A,  otherwise
  !           the leading  k by m  part of the array  A  must contain  the
  !           matrix A.
  !           Unchanged on exit.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
  !           LDA must be at least  max( 1, m ), otherwise  LDA must be at
  !           least  max( 1, k ).
  !           Unchanged on exit.
  !
  !  B      - REAL*4           array of DIMENSION ( LDB, kb ), where kb is
  !           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
  !           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
  !           part of the array  B  must contain the matrix  B,  otherwise
  !           the leading  n by k  part of the array  B  must contain  the
  !           matrix B.
  !           Unchanged on exit.
  !
  !  LDB    - INTEGER.
  !           On entry, LDB specifies the first dimension of B as declared
  !           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
  !           LDB must be at least  max( 1, k ), otherwise  LDB must be at
  !           least  max( 1, n ).
  !           Unchanged on exit.
  !
  !  BETA   - REAL*4          .
  !           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
  !           supplied as zero then C need not be set on input.
  !           Unchanged on exit.
  !
  !  C      - REAL*4           array of DIMENSION ( LDC, n ).
  !           Before entry, the leading  m by n  part of the array  C must
  !           contain the matrix  C,  except when  beta  is zero, in which
  !           case C need not be set on entry.
  !           On exit, the array  C  is overwritten by the  m by n  matrix
  !           ( alpha*op( A )*op( B ) + beta*C ).
  !
  !  LDC    - INTEGER.
  !           On entry, LDC specifies the first dimension of C as declared
  !           in  the  calling  (sub)  program.   LDC  must  be  at  least
  !           max( 1, m ).
  !           Unchanged on exit.
  !
  !
  !  Level 3 Blas routine.
  !
  !  -- Written on 8-February-1989.
  !     Jack Dongarra, Argonne National Laboratory.
  !     Iain Duff, AERE Harwell.
  !     Jeremy Du Croz, Numerical Algorithms Group Ltd.
  !     Sven Hammarling, Numerical Algorithms Group Ltd.
  !
  !
  !     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
  !     .. External Subroutines ..
  EXTERNAL           XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX
  !     .. Local Scalars ..
  LOGICAL            NOTA, NOTB
  INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
  REAL*4             TEMP
  !     .. Parameters ..
  REAL*4             ONE         , ZERO
  PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
  !     ..
  !     .. Executable Statements ..
  !
  !     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
  !     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
  !     and  columns of  A  and the  number of  rows  of  B  respectively.
  !
  NOTA  = LSAME( TRANSA, 'N' )
  NOTB  = LSAME( TRANSB, 'N' )
  IF( NOTA )THEN
     NROWA = M
     NCOLA = K
  ELSE
     NROWA = K
     NCOLA = M
  END IF
  IF( NOTB )THEN
     NROWB = K
  ELSE
     NROWB = N
  END IF
  !
  !     Test the input parameters.
  !
  INFO = 0
  IF(      ( .NOT.NOTA                 ).AND. &
       ( .NOT.LSAME( TRANSA, 'C' ) ).AND. &
       ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
     INFO = 1
  ELSE IF( ( .NOT.NOTB                 ).AND. &
       ( .NOT.LSAME( TRANSB, 'C' ) ).AND. &
       ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
     INFO = 2
  ELSE IF( M  .LT.0               )THEN
     INFO = 3
  ELSE IF( N  .LT.0               )THEN
     INFO = 4
  ELSE IF( K  .LT.0               )THEN
     INFO = 5
  ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
     INFO = 8
  ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
     INFO = 10
  ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
     INFO = 13
  END IF
  IF( INFO.NE.0 )THEN
     CALL XERBLA( 'SGEMM ', INFO )
     RETURN
  END IF
  !
  !     Quick return if possible.
  !
  IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR. &
       ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) ) &
       RETURN
  !
  !     And if  alpha.eq.zero.
  !
  IF( ALPHA.EQ.ZERO )THEN
     IF( BETA.EQ.ZERO )THEN
        do J =  1, N
           do I =  1, M
              C( I, J ) = ZERO
           end do
        end do
     ELSE
        do J =  1, N
           do I =  1, M
              C( I, J ) = BETA*C( I, J )
           end do
        end do
     END IF
     RETURN
  END IF
  !
  !     Start the operations.
  !
  IF( NOTB )THEN
     IF( NOTA )THEN
        !
        !           Form  C := alpha*A*B + beta*C.
        !
        do J =  1, N
           IF( BETA.EQ.ZERO )THEN
              do I =  1, M
                 C( I, J ) = ZERO
              end do
           ELSE IF( BETA.NE.ONE )THEN
              do I =  1, M
                 C( I, J ) = BETA*C( I, J )
              end do
           END IF
           DO L = 1, K
              IF( B( L, J ).NE.ZERO )THEN
                 TEMP = ALPHA*B( L, J )
                 do I =  1, M
                    C( I, J ) = C( I, J ) + TEMP*A( I, L )
                 end do
              END IF
           end do
        end do
     ELSE
        !
        !           Form  C := alpha*A'*B + beta*C
        !
        do J =  1, N
           do I =  1, M
              TEMP = ZERO
              DO L = 1, K
                 TEMP = TEMP + A( L, I )*B( L, J )
              end do
              IF( BETA.EQ.ZERO )THEN
                 C( I, J ) = ALPHA*TEMP
              ELSE
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              END IF
           end do
        end do
     END IF
  ELSE
     IF( NOTA )THEN
        !
        !           Form  C := alpha*A*B' + beta*C
        !
        do J =  1, N
           IF( BETA.EQ.ZERO )THEN
              do I =  1, M
                 C( I, J ) = ZERO
              end do
           ELSE IF( BETA.NE.ONE )THEN
              do I =  1, M
                 C( I, J ) = BETA*C( I, J )
              end do
           END IF
           DO L = 1, K
              IF( B( J, L ).NE.ZERO )THEN
                 TEMP = ALPHA*B( J, L )
                 do I =  1, M
                    C( I, J ) = C( I, J ) + TEMP*A( I, L )
                 end do
              END IF
           end do
        end do
     ELSE
        !
        !           Form  C := alpha*A'*B' + beta*C
        !
        do J =  1, N
           do I =  1, M
              TEMP = ZERO
              DO L = 1, K
                 TEMP = TEMP + A( L, I )*B( J, L )
              end do
              IF( BETA.EQ.ZERO )THEN
                 C( I, J ) = ALPHA*TEMP
              ELSE
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              END IF
           end do
        end do
     END IF
  END IF
  !
  RETURN
  !
  !     End of SGEMM .
  !
END SUBROUTINE SGEMM
!=============================================================================
SUBROUTINE SGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, &
     BETA, Y, INCY )
  !     .. Scalar Arguments ..
  REAL*4             ALPHA, BETA
  INTEGER            INCX, INCY, LDA, M, N
  CHARACTER*1        TRANS
  !     .. Array Arguments ..
  REAL*4             A( LDA, * ), X( * ), Y( * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  SGEMV  performs one of the matrix-vector operations
  !
  !     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
  !
  !  where alpha and beta are scalars, x and y are vectors and A is an
  !  m by n matrix.
  !
  !  Parameters
  !  ==========
  !
  !  TRANS  - CHARACTER*1.
  !           On entry, TRANS specifies the operation to be performed as
  !           follows:
  !
  !              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
  !
  !              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
  !
  !              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
  !
  !           Unchanged on exit.
  !
  !  M      - INTEGER.
  !           On entry, M specifies the number of rows of the matrix A.
  !           M must be at least zero.
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the number of columns of the matrix A.
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - REAL*4          .
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  A      - REAL*4           array of DIMENSION ( LDA, n ).
  !           Before entry, the leading m by n part of the array A must
  !           contain the matrix of coefficients.
  !           Unchanged on exit.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. LDA must be at least
  !           max( 1, m ).
  !           Unchanged on exit.
  !
  !  X      - REAL*4           array of DIMENSION at least
  !           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
  !           and at least
  !           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
  !           Before entry, the incremented array X must contain the
  !           vector x.
  !           Unchanged on exit.
  !
  !  INCX   - INTEGER.
  !           On entry, INCX specifies the increment for the elements of
  !           X. INCX must not be zero.
  !           Unchanged on exit.
  !
  !  BETA   - REAL*4          .
  !           On entry, BETA specifies the scalar beta. When BETA is
  !           supplied as zero then Y need not be set on input.
  !           Unchanged on exit.
  !
  !  Y      - REAL*4           array of DIMENSION at least
  !           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
  !           and at least
  !           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
  !           Before entry with BETA non-zero, the incremented array Y
  !           must contain the vector y. On exit, Y is overwritten by the
  !           updated vector y.
  !
  !  INCY   - INTEGER.
  !           On entry, INCY specifies the increment for the elements of
  !           Y. INCY must not be zero.
  !           Unchanged on exit.
  !
  !
  !  Level 2 Blas routine.
  !
  !  -- Written on 22-October-1986.
  !     Jack Dongarra, Argonne National Lab.
  !     Jeremy Du Croz, Nag Central Office.
  !     Sven Hammarling, Nag Central Office.
  !     Richard Hanson, Sandia National Labs.
  !
  !
  !     .. Parameters ..
  REAL*4             ONE         , ZERO
  PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
  !     .. Local Scalars ..
  REAL*4             TEMP
  INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
  !     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
  !     .. External Subroutines ..
  EXTERNAL           XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  INFO = 0
  IF     ( .NOT.LSAME( TRANS, 'N' ).AND. &
       .NOT.LSAME( TRANS, 'T' ).AND. &
       .NOT.LSAME( TRANS, 'C' )      )THEN
     INFO = 1
  ELSE IF( M.LT.0 )THEN
     INFO = 2
  ELSE IF( N.LT.0 )THEN
     INFO = 3
  ELSE IF( LDA.LT.MAX( 1, M ) )THEN
     INFO = 6
  ELSE IF( INCX.EQ.0 )THEN
     INFO = 8
  ELSE IF( INCY.EQ.0 )THEN
     INFO = 11
  END IF
  IF( INFO.NE.0 )THEN
     CALL XERBLA( 'SGEMV ', INFO )
     RETURN
  END IF
  !
  !     Quick return if possible.
  !
  IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR. &
       ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) ) &
       RETURN
  !
  !     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
  !     up the start points in  X  and  Y.
  !
  IF( LSAME( TRANS, 'N' ) )THEN
     LENX = N
     LENY = M
  ELSE
     LENX = M
     LENY = N
  END IF
  IF( INCX.GT.0 )THEN
     KX = 1
  ELSE
     KX = 1 - ( LENX - 1 )*INCX
  END IF
  IF( INCY.GT.0 )THEN
     KY = 1
  ELSE
     KY = 1 - ( LENY - 1 )*INCY
  END IF
  !
  !     Start the operations. In this version the elements of A are
  !     accessed sequentially with one pass through A.
  !
  !     First form  y := beta*y.
  !
  IF( BETA.NE.ONE )THEN
     IF( INCY.EQ.1 )THEN
        IF( BETA.EQ.ZERO )THEN
           do I =  1, LENY
              Y( I ) = ZERO
           end do
        ELSE
           do I =  1, LENY
              Y( I ) = BETA*Y( I )
           end do
        END IF
     ELSE
        IY = KY
        IF( BETA.EQ.ZERO )THEN
           do I =  1, LENY
              Y( IY ) = ZERO
              IY      = IY   + INCY
           end do
        ELSE
           do I =  1, LENY
              Y( IY ) = BETA*Y( IY )
              IY      = IY           + INCY
           end do
        END IF
     END IF
  END IF
  IF( ALPHA.EQ.ZERO ) &
       RETURN
  IF( LSAME( TRANS, 'N' ) )THEN
     !
     !        Form  y := alpha*A*x + y.
     !
     JX = KX
     IF( INCY.EQ.1 )THEN
        do J =  1, N
           IF( X( JX ).NE.ZERO )THEN
              TEMP = ALPHA*X( JX )
              do I =  1, M
                 Y( I ) = Y( I ) + TEMP*A( I, J )
              end do
           END IF
           JX = JX + INCX
        end do
     ELSE
        do J =  1, N
           IF( X( JX ).NE.ZERO )THEN
              TEMP = ALPHA*X( JX )
              IY   = KY
              do I =  1, M
                 Y( IY ) = Y( IY ) + TEMP*A( I, J )
                 IY      = IY      + INCY
              end do
           END IF
           JX = JX + INCX
        end do
     END IF
  ELSE
     !
     !        Form  y := alpha*A'*x + y.
     !
     JY = KY
     IF( INCX.EQ.1 )THEN
        do J =  1, N
           TEMP = ZERO
           do I =  1, M
              TEMP = TEMP + A( I, J )*X( I )
           end do
           Y( JY ) = Y( JY ) + ALPHA*TEMP
           JY      = JY      + INCY
        end do
     ELSE
        do J =  1, N
           TEMP = ZERO
           IX   = KX
           do I =  1, M
              TEMP = TEMP + A( I, J )*X( IX )
              IX   = IX   + INCX
           end do
           Y( JY ) = Y( JY ) + ALPHA*TEMP
           JY      = JY      + INCY
        end do
     END IF
  END IF
  !
  RETURN
  !
  !     End of SGEMV .
  !
END SUBROUTINE SGEMV
!=============================================================================
SUBROUTINE SGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
  !     .. Scalar Arguments ..
  REAL*4             ALPHA
  INTEGER            INCX, INCY, LDA, M, N
  !     .. Array Arguments ..
  REAL*4             A( LDA, * ), X( * ), Y( * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  SGER   performs the rank 1 operation
  !
  !     A := alpha*x*y' + A,
  !
  !  where alpha is a scalar, x is an m element vector, y is an n element
  !  vector and A is an m by n matrix.
  !
  !  Parameters
  !  ==========
  !
  !  M      - INTEGER.
  !           On entry, M specifies the number of rows of the matrix A.
  !           M must be at least zero.
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the number of columns of the matrix A.
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - REAL*4          .
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  X      - REAL*4           array of dimension at least
  !           ( 1 + ( m - 1 )*abs( INCX ) ).
  !           Before entry, the incremented array X must contain the m
  !           element vector x.
  !           Unchanged on exit.
  !
  !  INCX   - INTEGER.
  !           On entry, INCX specifies the increment for the elements of
  !           X. INCX must not be zero.
  !           Unchanged on exit.
  !
  !  Y      - REAL*4           array of dimension at least
  !           ( 1 + ( n - 1 )*abs( INCY ) ).
  !           Before entry, the incremented array Y must contain the n
  !           element vector y.
  !           Unchanged on exit.
  !
  !  INCY   - INTEGER.
  !           On entry, INCY specifies the increment for the elements of
  !           Y. INCY must not be zero.
  !           Unchanged on exit.
  !
  !  A      - REAL*4           array of DIMENSION ( LDA, n ).
  !           Before entry, the leading m by n part of the array A must
  !           contain the matrix of coefficients. On exit, A is
  !           overwritten by the updated matrix.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. LDA must be at least
  !           max( 1, m ).
  !           Unchanged on exit.
  !
  !
  !  Level 2 Blas routine.
  !
  !  -- Written on 22-October-1986.
  !     Jack Dongarra, Argonne National Lab.
  !     Jeremy Du Croz, Nag Central Office.
  !     Sven Hammarling, Nag Central Office.
  !     Richard Hanson, Sandia National Labs.
  !
  !
  !     .. Parameters ..
  REAL*4             ZERO
  PARAMETER        ( ZERO = 0.0E+0 )
  !     .. Local Scalars ..
  REAL*4             TEMP
  INTEGER            I, INFO, IX, J, JY, KX
  !     .. External Subroutines ..
  EXTERNAL           XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  INFO = 0
  IF     ( M.LT.0 )THEN
     INFO = 1
  ELSE IF( N.LT.0 )THEN
     INFO = 2
  ELSE IF( INCX.EQ.0 )THEN
     INFO = 5
  ELSE IF( INCY.EQ.0 )THEN
     INFO = 7
  ELSE IF( LDA.LT.MAX( 1, M ) )THEN
     INFO = 9
  END IF
  IF( INFO.NE.0 )THEN
     CALL XERBLA( 'SGER  ', INFO )
     RETURN
  END IF
  !
  !     Quick return if possible.
  !
  IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) ) &
       RETURN
  !
  !     Start the operations. In this version the elements of A are
  !     accessed sequentially with one pass through A.
  !
  IF( INCY.GT.0 )THEN
     JY = 1
  ELSE
     JY = 1 - ( N - 1 )*INCY
  END IF
  IF( INCX.EQ.1 )THEN
     do J =  1, N
        IF( Y( JY ).NE.ZERO )THEN
           TEMP = ALPHA*Y( JY )
           do I =  1, M
              A( I, J ) = A( I, J ) + X( I )*TEMP
           end do
        END IF
        JY = JY + INCY
     end do
  ELSE
     IF( INCX.GT.0 )THEN
        KX = 1
     ELSE
        KX = 1 - ( M - 1 )*INCX
     END IF
     do J =  1, N
        IF( Y( JY ).NE.ZERO )THEN
           TEMP = ALPHA*Y( JY )
           IX   = KX
           do I =  1, M
              A( I, J ) = A( I, J ) + X( IX )*TEMP
              IX        = IX        + INCX
           end do
        END IF
        JY = JY + INCY
     end do
  END IF
  !
  RETURN
  !
  !     End of SGER  .
  !
END SUBROUTINE SGER
!=============================================================================
subroutine sscal(n,sa,sx,incx)
  !
  !     scales a vector by a constant.
  !     uses unrolled loops for increment equal to 1.
  !     jack dongarra, linpack, 3/11/78.
  !     modified 3/93 to return if incx .le. 0.
  !     modified 12/3/93, array(1) declarations changed to array(*)
  !
  real*4  sa,sx(*)
  integer i,incx,m,mp1,n,nincx
  !
  if( n.le.0 .or. incx.le.0 )return

  nincx = n*incx
  do i =  1,nincx,incx
     sx(i) = sa*sx(i)
  end do
end subroutine sscal
!=============================================================================
subroutine sswap (n,sx,incx,sy,incy)
  !
  !     interchanges two vectors.
  !     uses unrolled loops for increments equal to 1.
  !     jack dongarra, linpack, 3/11/78.
  !     modified 12/3/93, array(1) declarations changed to array(*)
  !
  real*4  sx(*),sy(*),stemp
  integer i,incx,incy,ix,iy,m,mp1,n
  !
  if(n.le.0)return

  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do i =  1,n
     stemp = sx(ix)
     sx(ix) = sy(iy)
     sy(iy) = stemp
     ix = ix + incx
     iy = iy + incy
  end do
end subroutine sswap
!=============================================================================
SUBROUTINE STRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, &
     B, LDB )
  !     .. Scalar Arguments ..
  CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
  INTEGER            M, N, LDA, LDB
  REAL*4             ALPHA
  !     .. Array Arguments ..
  REAL*4             A( LDA, * ), B( LDB, * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  STRSM  solves one of the matrix equations
  !
  !     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
  !
  !  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
  !  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
  !
  !     op( A ) = A   or   op( A ) = A'.
  !
  !  The matrix X is overwritten on B.
  !
  !  Parameters
  !  ==========
  !
  !  SIDE   - CHARACTER*1.
  !           On entry, SIDE specifies whether op( A ) appears on the left
  !           or right of X as follows:
  !
  !              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
  !
  !              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
  !
  !           Unchanged on exit.
  !
  !  UPLO   - CHARACTER*1.
  !           On entry, UPLO specifies whether the matrix A is an upper or
  !           lower triangular matrix as follows:
  !
  !              UPLO = 'U' or 'u'   A is an upper triangular matrix.
  !
  !              UPLO = 'L' or 'l'   A is a lower triangular matrix.
  !
  !           Unchanged on exit.
  !
  !  TRANSA - CHARACTER*1.
  !           On entry, TRANSA specifies the form of op( A ) to be used in
  !           the matrix multiplication as follows:
  !
  !              TRANSA = 'N' or 'n'   op( A ) = A.
  !
  !              TRANSA = 'T' or 't'   op( A ) = A'.
  !
  !              TRANSA = 'C' or 'c'   op( A ) = A'.
  !
  !           Unchanged on exit.
  !
  !  DIAG   - CHARACTER*1.
  !           On entry, DIAG specifies whether or not A is unit triangular
  !           as follows:
  !
  !              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
  !
  !              DIAG = 'N' or 'n'   A is not assumed to be unit
  !                                  triangular.
  !
  !           Unchanged on exit.
  !
  !  M      - INTEGER.
  !           On entry, M specifies the number of rows of B. M must be at
  !           least zero.
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the number of columns of B.  N must be
  !           at least zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - REAL*4          .
  !           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
  !           zero then  A is not referenced and  B need not be set before
  !           entry.
  !           Unchanged on exit.
  !
  !  A      - REAL*4           array of DIMENSION ( LDA, k ), where k is m
  !           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
  !           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
  !           upper triangular part of the array  A must contain the upper
  !           triangular matrix  and the strictly lower triangular part of
  !           A is not referenced.
  !           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
  !           lower triangular part of the array  A must contain the lower
  !           triangular matrix  and the strictly upper triangular part of
  !           A is not referenced.
  !           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
  !           A  are not referenced either,  but are assumed to be  unity.
  !           Unchanged on exit.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
  !           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
  !           then LDA must be at least max( 1, n ).
  !           Unchanged on exit.
  !
  !  B      - REAL*4           array of DIMENSION ( LDB, n ).
  !           Before entry,  the leading  m by n part of the array  B must
  !           contain  the  right-hand  side  matrix  B,  and  on exit  is
  !           overwritten by the solution matrix  X.
  !
  !  LDB    - INTEGER.
  !           On entry, LDB specifies the first dimension of B as declared
  !           in  the  calling  (sub)  program.   LDB  must  be  at  least
  !           max( 1, m ).
  !           Unchanged on exit.
  !
  !
  !  Level 3 Blas routine.
  !
  !
  !  -- Written on 8-February-1989.
  !     Jack Dongarra, Argonne National Laboratory.
  !     Iain Duff, AERE Harwell.
  !     Jeremy Du Croz, Numerical Algorithms Group Ltd.
  !     Sven Hammarling, Numerical Algorithms Group Ltd.
  !
  !
  !     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
  !     .. External Subroutines ..
  EXTERNAL           XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX
  !     .. Local Scalars ..
  LOGICAL            LSIDE, NOUNIT, UPPER
  INTEGER            I, INFO, J, K, NROWA
  REAL*4             TEMP
  !     .. Parameters ..
  REAL*4             ONE         , ZERO
  PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  LSIDE  = LSAME( SIDE  , 'L' )
  IF( LSIDE )THEN
     NROWA = M
  ELSE
     NROWA = N
  END IF
  NOUNIT = LSAME( DIAG  , 'N' )
  UPPER  = LSAME( UPLO  , 'U' )
  !
  INFO   = 0
  IF(      ( .NOT.LSIDE                ).AND. &
       ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
     INFO = 1
  ELSE IF( ( .NOT.UPPER                ).AND. &
       ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
     INFO = 2
  ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND. &
       ( .NOT.LSAME( TRANSA, 'T' ) ).AND. &
       ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
     INFO = 3
  ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND. &
       ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
     INFO = 4
  ELSE IF( M  .LT.0               )THEN
     INFO = 5
  ELSE IF( N  .LT.0               )THEN
     INFO = 6
  ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
     INFO = 9
  ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
     INFO = 11
  END IF
  IF( INFO.NE.0 )THEN
     CALL XERBLA( 'STRSM ', INFO )
     RETURN
  END IF
  !
  !     Quick return if possible.
  !
  IF( N.EQ.0 ) &
       RETURN
  !
  !     And when  alpha.eq.zero.
  !
  IF( ALPHA.EQ.ZERO )THEN
     do J =  1, N
        do I =  1, M
           B( I, J ) = ZERO
        end do
     end do
     RETURN
  END IF
  !
  !     Start the operations.
  !
  IF( LSIDE )THEN
     IF( LSAME( TRANSA, 'N' ) )THEN
        !
        !           Form  B := alpha*inv( A )*B.
        !
        IF( UPPER )THEN
           do J =  1, N
              IF( ALPHA.NE.ONE )THEN
                 do I =  1, M
                    B( I, J ) = ALPHA*B( I, J )
                 end do
              END IF
              do K =  M, 1, -1
                 IF( B( K, J ).NE.ZERO )THEN
                    IF( NOUNIT ) &
                         B( K, J ) = B( K, J )/A( K, K )
                    do I =  1, K - 1
                       B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
                    end do
                 END IF
              end do
           end do
        ELSE
           do J =  1, N
              IF( ALPHA.NE.ONE )THEN
                 do I =  1, M
                    B( I, J ) = ALPHA*B( I, J )
                 end do
              END IF
              do K =  1, M
                 IF( B( K, J ).NE.ZERO )THEN
                    IF( NOUNIT ) &
                         B( K, J ) = B( K, J )/A( K, K )
                    do I =  K + 1, M
                       B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
                    end do
                 END IF
              end do
           end do
        END IF
     ELSE
        !
        !           Form  B := alpha*inv( A' )*B.
        !
        IF( UPPER )THEN
           do J =  1, N
              do I =  1, M
                 TEMP = ALPHA*B( I, J )
                 do K =  1, I - 1
                    TEMP = TEMP - A( K, I )*B( K, J )
                 end do
                 IF( NOUNIT ) &
                      TEMP = TEMP/A( I, I )
                 B( I, J ) = TEMP
              end do
           end do
        ELSE
           do J =  1, N
              do I =  M, 1, -1
                 TEMP = ALPHA*B( I, J )
                 do K =  I + 1, M
                    TEMP = TEMP - A( K, I )*B( K, J )
                 end do
                 IF( NOUNIT ) &
                      TEMP = TEMP/A( I, I )
                 B( I, J ) = TEMP
              end do
           end do
        END IF
     END IF
  ELSE
     IF( LSAME( TRANSA, 'N' ) )THEN
        !
        !           Form  B := alpha*B*inv( A ).
        !
        IF( UPPER )THEN
           do J =  1, N
              IF( ALPHA.NE.ONE )THEN
                 do I =  1, M
                    B( I, J ) = ALPHA*B( I, J )
                 end do
              END IF
              do K =  1, J - 1
                 IF( A( K, J ).NE.ZERO )THEN
                    do I =  1, M
                       B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
                    end do
                 END IF
              end do
              IF( NOUNIT )THEN
                 TEMP = ONE/A( J, J )
                 do I =  1, M
                    B( I, J ) = TEMP*B( I, J )
                 end do
              END IF
           end do
        ELSE
           do J =  N, 1, -1
              IF( ALPHA.NE.ONE )THEN
                 do I =  1, M
                    B( I, J ) = ALPHA*B( I, J )
                 end do
              END IF
              do K =  J + 1, N
                 IF( A( K, J ).NE.ZERO )THEN
                    do I =  1, M
                       B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
                    end do
                 END IF
              end do
              IF( NOUNIT )THEN
                 TEMP = ONE/A( J, J )
                 do I =  1, M
                    B( I, J ) = TEMP*B( I, J )
                 end do
              END IF
           end do
        END IF
     ELSE
        !
        !           Form  B := alpha*B*inv( A' ).
        !
        IF( UPPER )THEN
           do K =  N, 1, -1
              IF( NOUNIT )THEN
                 TEMP = ONE/A( K, K )
                 do I =  1, M
                    B( I, K ) = TEMP*B( I, K )
                 end do
              END IF
              do J =  1, K - 1
                 IF( A( J, K ).NE.ZERO )THEN
                    TEMP = A( J, K )
                    do I =  1, M
                       B( I, J ) = B( I, J ) - TEMP*B( I, K )
                    end do
                 END IF
              end do
              IF( ALPHA.NE.ONE )THEN
                 do I =  1, M
                    B( I, K ) = ALPHA*B( I, K )
                 end do
              END IF
           end do
        ELSE
           do K =  1, N
              IF( NOUNIT )THEN
                 TEMP = ONE/A( K, K )
                 do I =  1, M
                    B( I, K ) = TEMP*B( I, K )
                 end do
              END IF
              do J =  K + 1, N
                 IF( A( J, K ).NE.ZERO )THEN
                    TEMP = A( J, K )
                    do I =  1, M
                       B( I, J ) = B( I, J ) - TEMP*B( I, K )
                    end do
                 END IF
              end do
              IF( ALPHA.NE.ONE )THEN
                 do I =  1, M
                    B( I, K ) = ALPHA*B( I, K )
                 end do
              END IF
           end do
        END IF
     END IF
  END IF
  !
  RETURN
  !
  !     End of STRSM .
  !
END SUBROUTINE STRSM
!=============================================================================
integer function isamax(n,sx,incx)
  !
  !     finds the index of element having max. absolute value.
  !     jack dongarra, linpack, 3/11/78.
  !     modified 3/93 to return if incx .le. 0.
  !     modified 12/3/93, array(1) declarations changed to array(*)
  !
  real*4  sx(*),smax
  integer i,incx,ix,n
  !
  isamax = 0
  if( n.lt.1 .or. incx.le.0 ) return
  isamax = 1
  if(n.eq.1)return

  ix = 1
  smax = abs(sx(1))
  ix = ix + incx
  do i =  2,n
     if(.not. (abs(sx(ix)).le.smax)) then
        isamax = i
        smax = abs(sx(ix))
     end if
     ix = ix + incx
  end do
  return
end function isamax
!=============================================================================
!============================================================================
! This is a collection of real*8 BLAS routines that BATSRUS uses. 
! You are encouraged to use the local BLAS library if available.
!
! Subroutines: dcopy, dgemm, dgemv, dger,  dscal, dswap, dtrsm
!
!=============================================================================
subroutine  dcopy(n,dx,incx,dy,incy)
  !
  !     copies a vector, x, to a vector, y.
  !     uses unrolled loops for increments equal to one.
  !     jack dongarra, linpack, 3/11/78.
  !
  real*8 dx(*),dy(*)
  integer i,incx,incy,ix,iy,m,mp1,n
  !
  if(n.le.0)return

  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do i =  1,n
     dy(iy) = dx(ix)
     ix = ix + incx
     iy = iy + incy
  end do
  return
end subroutine dcopy
!=============================================================================
SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
     BETA, C, LDC )
  !     .. Scalar Arguments ..
  CHARACTER*1        TRANSA, TRANSB
  INTEGER            M, N, K, LDA, LDB, LDC
  REAL*8   ALPHA, BETA
  !     .. Array Arguments ..
  REAL*8   A( LDA, * ), B( LDB, * ), C( LDC, * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  DGEMM  performs one of the matrix-matrix operations
  !
  !     C := alpha*op( A )*op( B ) + beta*C,
  !
  !  where  op( X ) is one of
  !
  !     op( X ) = X   or   op( X ) = X',
  !
  !  alpha and beta are scalars, and A, B and C are matrices, with op( A )
  !  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
  !
  !  Parameters
  !  ==========
  !
  !  TRANSA - CHARACTER*1.
  !           On entry, TRANSA specifies the form of op( A ) to be used in
  !           the matrix multiplication as follows:
  !
  !              TRANSA = 'N' or 'n',  op( A ) = A.
  !
  !              TRANSA = 'T' or 't',  op( A ) = A'.
  !
  !              TRANSA = 'C' or 'c',  op( A ) = A'.
  !
  !           Unchanged on exit.
  !
  !  TRANSB - CHARACTER*1.
  !           On entry, TRANSB specifies the form of op( B ) to be used in
  !           the matrix multiplication as follows:
  !
  !              TRANSB = 'N' or 'n',  op( B ) = B.
  !
  !              TRANSB = 'T' or 't',  op( B ) = B'.
  !
  !              TRANSB = 'C' or 'c',  op( B ) = B'.
  !
  !           Unchanged on exit.
  !
  !  M      - INTEGER.
  !           On entry,  M  specifies  the number  of rows  of the  matrix
  !           op( A )  and of the  matrix  C.  M  must  be at least  zero.
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry,  N  specifies the number  of columns of the matrix
  !           op( B ) and the number of columns of the matrix C. N must be
  !           at least zero.
  !           Unchanged on exit.
  !
  !  K      - INTEGER.
  !           On entry,  K  specifies  the number of columns of the matrix
  !           op( A ) and the number of rows of the matrix op( B ). K must
  !           be at least  zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - REAL*8.
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  A      - REAL*8 array of DIMENSION ( LDA, ka ), where ka is
  !           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
  !           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
  !           part of the array  A  must contain the matrix  A,  otherwise
  !           the leading  k by m  part of the array  A  must contain  the
  !           matrix A.
  !           Unchanged on exit.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
  !           LDA must be at least  max( 1, m ), otherwise  LDA must be at
  !           least  max( 1, k ).
  !           Unchanged on exit.
  !
  !  B      - REAL*8 array of DIMENSION ( LDB, kb ), where kb is
  !           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
  !           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
  !           part of the array  B  must contain the matrix  B,  otherwise
  !           the leading  n by k  part of the array  B  must contain  the
  !           matrix B.
  !           Unchanged on exit.
  !
  !  LDB    - INTEGER.
  !           On entry, LDB specifies the first dimension of B as declared
  !           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
  !           LDB must be at least  max( 1, k ), otherwise  LDB must be at
  !           least  max( 1, n ).
  !           Unchanged on exit.
  !
  !  BETA   - REAL*8.
  !           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
  !           supplied as zero then C need not be set on input.
  !           Unchanged on exit.
  !
  !  C      - REAL*8 array of DIMENSION ( LDC, n ).
  !           Before entry, the leading  m by n  part of the array  C must
  !           contain the matrix  C,  except when  beta  is zero, in which
  !           case C need not be set on entry.
  !           On exit, the array  C  is overwritten by the  m by n  matrix
  !           ( alpha*op( A )*op( B ) + beta*C ).
  !
  !  LDC    - INTEGER.
  !           On entry, LDC specifies the first dimension of C as declared
  !           in  the  calling  (sub)  program.   LDC  must  be  at  least
  !           max( 1, m ).
  !           Unchanged on exit.
  !
  !
  !  Level 3 Blas routine.
  !
  !  -- Written on 8-February-1989.
  !     Jack Dongarra, Argonne National Laboratory.
  !     Iain Duff, AERE Harwell.
  !     Jeremy Du Croz, Numerical Algorithms Group Ltd.
  !     Sven Hammarling, Numerical Algorithms Group Ltd.
  !
  !
  !     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
  !     .. External Subroutines ..
  EXTERNAL           XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX
  !     .. Local Scalars ..
  LOGICAL            NOTA, NOTB
  INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
  REAL*8   TEMP
  !     .. Parameters ..
  REAL*8   ONE    , ZERO
  PARAMETER        ( ONE = 1, ZERO = 0 )
  !     ..
  !     .. Executable Statements ..
  !
  !     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
  !     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
  !     and  columns of  A  and the  number of  rows  of  B  respectively.
  !
  NOTA  = LSAME( TRANSA, 'N' )
  NOTB  = LSAME( TRANSB, 'N' )
  IF( NOTA )THEN
     NROWA = M
     NCOLA = K
  ELSE
     NROWA = K
     NCOLA = M
  END IF
  IF( NOTB )THEN
     NROWB = K
  ELSE
     NROWB = N
  END IF
  !
  !     Test the input parameters.
  !
  INFO = 0
  IF(      ( .NOT.NOTA                 ).AND. &
       ( .NOT.LSAME( TRANSA, 'C' ) ).AND. &
       ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
     INFO = 1
  ELSE IF( ( .NOT.NOTB                 ).AND. &
       ( .NOT.LSAME( TRANSB, 'C' ) ).AND. &
       ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
     INFO = 2
  ELSE IF( M  .LT.0               )THEN
     INFO = 3
  ELSE IF( N  .LT.0               )THEN
     INFO = 4
  ELSE IF( K  .LT.0               )THEN
     INFO = 5
  ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
     INFO = 8
  ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
     INFO = 10
  ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
     INFO = 13
  END IF
  IF( INFO.NE.0 )THEN
     CALL XERBLA( 'DGEMM ', INFO )
     RETURN
  END IF
  !
  !     Quick return if possible.
  !
  IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR. &
       ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) ) &
       RETURN
  !
  !     And if  alpha.eq.zero.
  !
  IF( ALPHA.EQ.ZERO )THEN
     IF( BETA.EQ.ZERO )THEN
        do J =  1, N
           do I =  1, M
              C( I, J ) = ZERO
           end do
        end do
     ELSE
        do J =  1, N
           do I =  1, M
              C( I, J ) = BETA*C( I, J )
           end do
        end do
     END IF
     RETURN
  END IF
  !
  !     Start the operations.
  !
  IF( NOTB )THEN
     IF( NOTA )THEN
        !
        !           Form  C := alpha*A*B + beta*C.
        !
        do J =  1, N
           IF( BETA.EQ.ZERO )THEN
              do I =  1, M
                 C( I, J ) = ZERO
              end do
           ELSE IF( BETA.NE.ONE )THEN
              do I =  1, M
                 C( I, J ) = BETA*C( I, J )
              end do
           END IF
           DO  L = 1, K
              IF( B( L, J ).NE.ZERO )THEN
                 TEMP = ALPHA*B( L, J )
                 do I =  1, M
                    C( I, J ) = C( I, J ) + TEMP*A( I, L )
                 end do
              END IF
           end do
        end do
     ELSE
        !
        !           Form  C := alpha*A'*B + beta*C
        !
        do J =  1, N
           do I =  1, M
              TEMP = ZERO
              DO L = 1, K
                 TEMP = TEMP + A( L, I )*B( L, J )
              end do
              IF( BETA.EQ.ZERO )THEN
                 C( I, J ) = ALPHA*TEMP
              ELSE
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              END IF
           end do
        end do
     END IF
  ELSE
     IF( NOTA )THEN
        !
        !           Form  C := alpha*A*B' + beta*C
        !
        do J =  1, N
           IF( BETA.EQ.ZERO )THEN
              do I =  1, M
                 C( I, J ) = ZERO
              end do
           ELSE IF( BETA.NE.ONE )THEN
              do I =  1, M
                 C( I, J ) = BETA*C( I, J )
              end do
           END IF
           DO L = 1, K
              IF( B( J, L ).NE.ZERO )THEN
                 TEMP = ALPHA*B( J, L )
                 do I =  1, M
                    C( I, J ) = C( I, J ) + TEMP*A( I, L )
                 end do
              END IF
           end do
        end do
     ELSE
        !
        !           Form  C := alpha*A'*B' + beta*C
        !
        do J =  1, N
           do I =  1, M
              TEMP = ZERO
              DO L = 1, K
                 TEMP = TEMP + A( L, I )*B( J, L )
              end do
              IF( BETA.EQ.ZERO )THEN
                 C( I, J ) = ALPHA*TEMP
              ELSE
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              END IF
           end do
        end do
     END IF
  END IF
  !
  RETURN
  !
  !     End of DGEMM .
  !
END SUBROUTINE DGEMM
!============================================================================
SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, &
     BETA, Y, INCY )
  !     .. Scalar Arguments ..
  REAL*8   ALPHA, BETA
  INTEGER            INCX, INCY, LDA, M, N
  CHARACTER*1        TRANS
  !     .. Array Arguments ..
  REAL*8   A( LDA, * ), X( * ), Y( * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  DGEMV  performs one of the matrix-vector operations
  !
  !     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
  !
  !  where alpha and beta are scalars, x and y are vectors and A is an
  !  m by n matrix.
  !
  !  Parameters
  !  ==========
  !
  !  TRANS  - CHARACTER*1.
  !           On entry, TRANS specifies the operation to be performed as
  !           follows:
  !
  !              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
  !
  !              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
  !
  !              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
  !
  !           Unchanged on exit.
  !
  !  M      - INTEGER.
  !           On entry, M specifies the number of rows of the matrix A.
  !           M must be at least zero.
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the number of columns of the matrix A.
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - REAL*8.
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  A      - REAL*8 array of DIMENSION ( LDA, n ).
  !           Before entry, the leading m by n part of the array A must
  !           contain the matrix of coefficients.
  !           Unchanged on exit.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. LDA must be at least
  !           max( 1, m ).
  !           Unchanged on exit.
  !
  !  X      - REAL*8 array of DIMENSION at least
  !           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
  !           and at least
  !           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
  !           Before entry, the incremented array X must contain the
  !           vector x.
  !           Unchanged on exit.
  !
  !  INCX   - INTEGER.
  !           On entry, INCX specifies the increment for the elements of
  !           X. INCX must not be zero.
  !           Unchanged on exit.
  !
  !  BETA   - REAL*8.
  !           On entry, BETA specifies the scalar beta. When BETA is
  !           supplied as zero then Y need not be set on input.
  !           Unchanged on exit.
  !
  !  Y      - REAL*8 array of DIMENSION at least
  !           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
  !           and at least
  !           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
  !           Before entry with BETA non-zero, the incremented array Y
  !           must contain the vector y. On exit, Y is overwritten by the
  !           updated vector y.
  !
  !  INCY   - INTEGER.
  !           On entry, INCY specifies the increment for the elements of
  !           Y. INCY must not be zero.
  !           Unchanged on exit.
  !
  !
  !  Level 2 Blas routine.
  !
  !  -- Written on 22-October-1986.
  !     Jack Dongarra, Argonne National Lab.
  !     Jeremy Du Croz, Nag Central Office.
  !     Sven Hammarling, Nag Central Office.
  !     Richard Hanson, Sandia National Labs.
  !
  !
  !     .. Parameters ..
  REAL*8   ONE         , ZERO
  PARAMETER        ( ONE = 1, ZERO = 0 )
  !     .. Local Scalars ..
  REAL*8   TEMP
  INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
  !     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
  !     .. External Subroutines ..
  EXTERNAL           XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  INFO = 0
  IF     ( .NOT.LSAME( TRANS, 'N' ).AND. &
       .NOT.LSAME( TRANS, 'T' ).AND. &
       .NOT.LSAME( TRANS, 'C' )      )THEN
     INFO = 1
  ELSE IF( M.LT.0 )THEN
     INFO = 2
  ELSE IF( N.LT.0 )THEN
     INFO = 3
  ELSE IF( LDA.LT.MAX( 1, M ) )THEN
     INFO = 6
  ELSE IF( INCX.EQ.0 )THEN
     INFO = 8
  ELSE IF( INCY.EQ.0 )THEN
     INFO = 11
  END IF
  IF( INFO.NE.0 )THEN
     CALL XERBLA( 'DGEMV ', INFO )
     RETURN
  END IF
  !
  !     Quick return if possible.
  !
  IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR. &
       ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) ) &
       RETURN
  !
  !     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
  !     up the start points in  X  and  Y.
  !
  IF( LSAME( TRANS, 'N' ) )THEN
     LENX = N
     LENY = M
  ELSE
     LENX = M
     LENY = N
  END IF
  IF( INCX.GT.0 )THEN
     KX = 1
  ELSE
     KX = 1 - ( LENX - 1 )*INCX
  END IF
  IF( INCY.GT.0 )THEN
     KY = 1
  ELSE
     KY = 1 - ( LENY - 1 )*INCY
  END IF
  !
  !     Start the operations. In this version the elements of A are
  !     accessed sequentially with one pass through A.
  !
  !     First form  y := beta*y.
  !
  IF( BETA.NE.ONE )THEN
     IF( INCY.EQ.1 )THEN
        IF( BETA.EQ.ZERO )THEN
           do I =  1, LENY
              Y( I ) = ZERO
           end do
        ELSE
           do I =  1, LENY
              Y( I ) = BETA*Y( I )
           end do
        END IF
     ELSE
        IY = KY
        IF( BETA.EQ.ZERO )THEN
           do I =  1, LENY
              Y( IY ) = ZERO
              IY      = IY   + INCY
           end do
        ELSE
           do I =  1, LENY
              Y( IY ) = BETA*Y( IY )
              IY      = IY           + INCY
           end do
        END IF
     END IF
  END IF
  IF( ALPHA.EQ.ZERO ) &
       RETURN
  IF( LSAME( TRANS, 'N' ) )THEN
     !
     !        Form  y := alpha*A*x + y.
     !
     JX = KX
     IF( INCY.EQ.1 )THEN
        do J =  1, N
           IF( X( JX ).NE.ZERO )THEN
              TEMP = ALPHA*X( JX )
              do I =  1, M
                 Y( I ) = Y( I ) + TEMP*A( I, J )
              end do
           END IF
           JX = JX + INCX
        end do
     ELSE
        do J =  1, N
           IF( X( JX ).NE.ZERO )THEN
              TEMP = ALPHA*X( JX )
              IY   = KY
              do I =  1, M
                 Y( IY ) = Y( IY ) + TEMP*A( I, J )
                 IY      = IY      + INCY
              end do
           END IF
           JX = JX + INCX
        end do
     END IF
  ELSE
     !
     !        Form  y := alpha*A'*x + y.
     !
     JY = KY
     IF( INCX.EQ.1 )THEN
        do J =  1, N
           TEMP = ZERO
           do I =  1, M
              TEMP = TEMP + A( I, J )*X( I )
           end do
           Y( JY ) = Y( JY ) + ALPHA*TEMP
           JY      = JY      + INCY
        end do
     ELSE
        do J =  1, N
           TEMP = ZERO
           IX   = KX
           do I =  1, M
              TEMP = TEMP + A( I, J )*X( IX )
              IX   = IX   + INCX
           end do
           Y( JY ) = Y( JY ) + ALPHA*TEMP
           JY      = JY      + INCY
        end do
     END IF
  END IF
  !
  RETURN
  !
  !     End of DGEMV .
  !
END SUBROUTINE DGEMV
!============================================================================
SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
  !     .. Scalar Arguments ..
  REAL*8   ALPHA
  INTEGER            INCX, INCY, LDA, M, N
  !     .. Array Arguments ..
  REAL*8   A( LDA, * ), X( * ), Y( * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  DGER   performs the rank 1 operation
  !
  !     A := alpha*x*y' + A,
  !
  !  where alpha is a scalar, x is an m element vector, y is an n element
  !  vector and A is an m by n matrix.
  !
  !  Parameters
  !  ==========
  !
  !  M      - INTEGER.
  !           On entry, M specifies the number of rows of the matrix A.
  !           M must be at least zero.
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the number of columns of the matrix A.
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - REAL*8.
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  X      - REAL*8 array of dimension at least
  !           ( 1 + ( m - 1 )*abs( INCX ) ).
  !           Before entry, the incremented array X must contain the m
  !           element vector x.
  !           Unchanged on exit.
  !
  !  INCX   - INTEGER.
  !           On entry, INCX specifies the increment for the elements of
  !           X. INCX must not be zero.
  !           Unchanged on exit.
  !
  !  Y      - REAL*8 array of dimension at least
  !           ( 1 + ( n - 1 )*abs( INCY ) ).
  !           Before entry, the incremented array Y must contain the n
  !           element vector y.
  !           Unchanged on exit.
  !
  !  INCY   - INTEGER.
  !           On entry, INCY specifies the increment for the elements of
  !           Y. INCY must not be zero.
  !           Unchanged on exit.
  !
  !  A      - REAL*8 array of DIMENSION ( LDA, n ).
  !           Before entry, the leading m by n part of the array A must
  !           contain the matrix of coefficients. On exit, A is
  !           overwritten by the updated matrix.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. LDA must be at least
  !           max( 1, m ).
  !           Unchanged on exit.
  !
  !
  !  Level 2 Blas routine.
  !
  !  -- Written on 22-October-1986.
  !     Jack Dongarra, Argonne National Lab.
  !     Jeremy Du Croz, Nag Central Office.
  !     Sven Hammarling, Nag Central Office.
  !     Richard Hanson, Sandia National Labs.
  !
  !
  !     .. Parameters ..
  REAL*8   ZERO
  PARAMETER        ( ZERO = 0.0D+0 )
  !     .. Local Scalars ..
  REAL*8   TEMP
  INTEGER            I, INFO, IX, J, JY, KX
  !     .. External Subroutines ..
  EXTERNAL           XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  INFO = 0
  IF     ( M.LT.0 )THEN
     INFO = 1
  ELSE IF( N.LT.0 )THEN
     INFO = 2
  ELSE IF( INCX.EQ.0 )THEN
     INFO = 5
  ELSE IF( INCY.EQ.0 )THEN
     INFO = 7
  ELSE IF( LDA.LT.MAX( 1, M ) )THEN
     INFO = 9
  END IF
  IF( INFO.NE.0 )THEN
     CALL XERBLA( 'DGER  ', INFO )
     RETURN
  END IF
  !
  !     Quick return if possible.
  !
  IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) ) &
       RETURN
  !
  !     Start the operations. In this version the elements of A are
  !     accessed sequentially with one pass through A.
  !
  IF( INCY.GT.0 )THEN
     JY = 1
  ELSE
     JY = 1 - ( N - 1 )*INCY
  END IF
  IF( INCX.EQ.1 )THEN
     do J =  1, N
        IF( Y( JY ).NE.ZERO )THEN
           TEMP = ALPHA*Y( JY )
           do I =  1, M
              A( I, J ) = A( I, J ) + X( I )*TEMP
           end do
        END IF
        JY = JY + INCY
     end do
  ELSE
     IF( INCX.GT.0 )THEN
        KX = 1
     ELSE
        KX = 1 - ( M - 1 )*INCX
     END IF
     do J =  1, N
        IF( Y( JY ).NE.ZERO )THEN
           TEMP = ALPHA*Y( JY )
           IX   = KX
           do I =  1, M
              A( I, J ) = A( I, J ) + X( IX )*TEMP
              IX        = IX        + INCX
           end do
        END IF
        JY = JY + INCY
     end do
  END IF
  !
  RETURN
  !
  !     End of DGER  .
  !
END SUBROUTINE DGER
!============================================================================
subroutine  dscal(n,da,dx,incx)
  !
  !     scales a vector by a constant.
  !     uses unrolled loops for increment equal to one.
  !     jack dongarra, linpack, 3/11/78.
  !     modified 3/93 to return if incx .le. 0.
  !
  real*8 da,dx(*)
  integer i,incx,m,mp1,n,nincx
  !
  if( n.le.0 .or. incx.le.0 )return
  nincx = n*incx
  do i =  1,nincx,incx
     dx(i) = da*dx(i)
  end do
  return
end subroutine dscal
!=============================================================================
subroutine  dswap (n,dx,incx,dy,incy)
  !
  !     interchanges two vectors.
  !     uses unrolled loops for increments equal one.
  !     jack dongarra, linpack, 3/11/78.
  !
  real*8 dx(*),dy(*),dtemp
  integer i,incx,incy,ix,iy,m,mp1,n
  !
  if(n.le.0)return

  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do i =  1,n
     dtemp = dx(ix)
     dx(ix) = dy(iy)
     dy(iy) = dtemp
     ix = ix + incx
     iy = iy + incy
  end do
  return
end subroutine dswap
!============================================================================
SUBROUTINE DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, &
     B, LDB )
  !     .. Scalar Arguments ..
  CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
  INTEGER            M, N, LDA, LDB
  REAL*8   ALPHA
  !     .. Array Arguments ..
  REAL*8   A( LDA, * ), B( LDB, * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  DTRSM  solves one of the matrix equations
  !
  !     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
  !
  !  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
  !  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
  !
  !     op( A ) = A   or   op( A ) = A'.
  !
  !  The matrix X is overwritten on B.
  !
  !  Parameters
  !  ==========
  !
  !  SIDE   - CHARACTER*1.
  !           On entry, SIDE specifies whether op( A ) appears on the left
  !           or right of X as follows:
  !
  !              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
  !
  !              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
  !
  !           Unchanged on exit.
  !
  !  UPLO   - CHARACTER*1.
  !           On entry, UPLO specifies whether the matrix A is an upper or
  !           lower triangular matrix as follows:
  !
  !              UPLO = 'U' or 'u'   A is an upper triangular matrix.
  !
  !              UPLO = 'L' or 'l'   A is a lower triangular matrix.
  !
  !           Unchanged on exit.
  !
  !  TRANSA - CHARACTER*1.
  !           On entry, TRANSA specifies the form of op( A ) to be used in
  !           the matrix multiplication as follows:
  !
  !              TRANSA = 'N' or 'n'   op( A ) = A.
  !
  !              TRANSA = 'T' or 't'   op( A ) = A'.
  !
  !              TRANSA = 'C' or 'c'   op( A ) = A'.
  !
  !           Unchanged on exit.
  !
  !  DIAG   - CHARACTER*1.
  !           On entry, DIAG specifies whether or not A is unit triangular
  !           as follows:
  !
  !              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
  !
  !              DIAG = 'N' or 'n'   A is not assumed to be unit
  !                                  triangular.
  !
  !           Unchanged on exit.
  !
  !  M      - INTEGER.
  !           On entry, M specifies the number of rows of B. M must be at
  !           least zero.
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the number of columns of B.  N must be
  !           at least zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - REAL*8.
  !           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
  !           zero then  A is not referenced and  B need not be set before
  !           entry.
  !           Unchanged on exit.
  !
  !  A      - REAL*8 array of DIMENSION ( LDA, k ), where k is m
  !           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
  !           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
  !           upper triangular part of the array  A must contain the upper
  !           triangular matrix  and the strictly lower triangular part of
  !           A is not referenced.
  !           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
  !           lower triangular part of the array  A must contain the lower
  !           triangular matrix  and the strictly upper triangular part of
  !           A is not referenced.
  !           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
  !           A  are not referenced either,  but are assumed to be  unity.
  !           Unchanged on exit.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
  !           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
  !           then LDA must be at least max( 1, n ).
  !           Unchanged on exit.
  !
  !  B      - REAL*8 array of DIMENSION ( LDB, n ).
  !           Before entry,  the leading  m by n part of the array  B must
  !           contain  the  right-hand  side  matrix  B,  and  on exit  is
  !           overwritten by the solution matrix  X.
  !
  !  LDB    - INTEGER.
  !           On entry, LDB specifies the first dimension of B as declared
  !           in  the  calling  (sub)  program.   LDB  must  be  at  least
  !           max( 1, m ).
  !           Unchanged on exit.
  !
  !
  !  Level 3 Blas routine.
  !
  !
  !  -- Written on 8-February-1989.
  !     Jack Dongarra, Argonne National Laboratory.
  !     Iain Duff, AERE Harwell.
  !     Jeremy Du Croz, Numerical Algorithms Group Ltd.
  !     Sven Hammarling, Numerical Algorithms Group Ltd.
  !
  !
  !     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
  !     .. External Subroutines ..
  EXTERNAL           XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX
  !     .. Local Scalars ..
  LOGICAL            LSIDE, NOUNIT, UPPER
  INTEGER            I, INFO, J, K, NROWA
  REAL*8             TEMP
  !     .. Parameters ..
  REAL*8             ONE    , ZERO
  PARAMETER        ( ONE = 1, ZERO = 0 )
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  LSIDE  = LSAME( SIDE  , 'L' )
  IF( LSIDE )THEN
     NROWA = M
  ELSE
     NROWA = N
  END IF
  NOUNIT = LSAME( DIAG  , 'N' )
  UPPER  = LSAME( UPLO  , 'U' )
  !
  INFO   = 0
  IF(      ( .NOT.LSIDE                ).AND. &
       ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
     INFO = 1
  ELSE IF( ( .NOT.UPPER                ).AND. &
       ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
     INFO = 2
  ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND. &
       ( .NOT.LSAME( TRANSA, 'T' ) ).AND. &
       ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
     INFO = 3
  ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND. &
       ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
     INFO = 4
  ELSE IF( M  .LT.0               )THEN
     INFO = 5
  ELSE IF( N  .LT.0               )THEN
     INFO = 6
  ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
     INFO = 9
  ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
     INFO = 11
  END IF
  IF( INFO.NE.0 )THEN
     CALL XERBLA( 'DTRSM ', INFO )
     RETURN
  END IF
  !
  !     Quick return if possible.
  !
  IF( N.EQ.0 ) &
       RETURN
  !
  !     And when  alpha.eq.zero.
  !
  IF( ALPHA.EQ.ZERO )THEN
     do J =  1, N
        do I =  1, M
           B( I, J ) = ZERO
        end do
     end do
     RETURN
  END IF
  !
  !     Start the operations.
  !
  IF( LSIDE )THEN
     IF( LSAME( TRANSA, 'N' ) )THEN
        !
        !           Form  B := alpha*inv( A )*B.
        !
        IF( UPPER )THEN
           do J =  1, N
              IF( ALPHA.NE.ONE )THEN
                 do I =  1, M
                    B( I, J ) = ALPHA*B( I, J )
                 end do
              END IF
              do K =  M, 1, -1
                 IF( B( K, J ).NE.ZERO )THEN
                    IF( NOUNIT ) &
                         B( K, J ) = B( K, J )/A( K, K )
                    do I =  1, K - 1
                       B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
                    end do
                 END IF
              end do
           end do
        ELSE
           do J =  1, N
              IF( ALPHA.NE.ONE )THEN
                 do I =  1, M
                    B( I, J ) = ALPHA*B( I, J )
                 end do
              END IF
              do K =  1, M
                 IF( B( K, J ).NE.ZERO )THEN
                    IF( NOUNIT ) &
                         B( K, J ) = B( K, J )/A( K, K )
                    do I =  K + 1, M
                       B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
                    end do
                 END IF
              end do
           end do
        END IF
     ELSE
        !
        !           Form  B := alpha*inv( A' )*B.
        !
        IF( UPPER )THEN
           do J =  1, N
              do I =  1, M
                 TEMP = ALPHA*B( I, J )
                 do K =  1, I - 1
                    TEMP = TEMP - A( K, I )*B( K, J )
                 end do
                 IF( NOUNIT ) &
                      TEMP = TEMP/A( I, I )
                 B( I, J ) = TEMP
              end do
           end do
        ELSE
           do J =  1, N
              do I =  M, 1, -1
                 TEMP = ALPHA*B( I, J )
                 do K =  I + 1, M
                    TEMP = TEMP - A( K, I )*B( K, J )
                 end do
                 IF( NOUNIT ) &
                      TEMP = TEMP/A( I, I )
                 B( I, J ) = TEMP
              end do
           end do
        END IF
     END IF
  ELSE
     IF( LSAME( TRANSA, 'N' ) )THEN
        !
        !           Form  B := alpha*B*inv( A ).
        !
        IF( UPPER )THEN
           do J =  1, N
              IF( ALPHA.NE.ONE )THEN
                 do I =  1, M
                    B( I, J ) = ALPHA*B( I, J )
                 end do
              END IF
              do K =  1, J - 1
                 IF( A( K, J ).NE.ZERO )THEN
                    do I =  1, M
                       B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
                    end do
                 END IF
              end do
              IF( NOUNIT )THEN
                 TEMP = ONE/A( J, J )
                 do I =  1, M
                    B( I, J ) = TEMP*B( I, J )
                 end do
              END IF
           end do
        ELSE
           do J =  N, 1, -1
              IF( ALPHA.NE.ONE )THEN
                 do I =  1, M
                    B( I, J ) = ALPHA*B( I, J )
                 end do
              END IF
              do K =  J + 1, N
                 IF( A( K, J ).NE.ZERO )THEN
                    do I =  1, M
                       B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
                    end do
                 END IF
              end do
              IF( NOUNIT )THEN
                 TEMP = ONE/A( J, J )
                 do I =  1, M
                    B( I, J ) = TEMP*B( I, J )
                 end do
              END IF
           end do
        END IF
     ELSE
        !
        !           Form  B := alpha*B*inv( A' ).
        !
        IF( UPPER )THEN
           do K =  N, 1, -1
              IF( NOUNIT )THEN
                 TEMP = ONE/A( K, K )
                 do I =  1, M
                    B( I, K ) = TEMP*B( I, K )
                 end do
              END IF
              do J =  1, K - 1
                 IF( A( J, K ).NE.ZERO )THEN
                    TEMP = A( J, K )
                    do I =  1, M
                       B( I, J ) = B( I, J ) - TEMP*B( I, K )
                    end do
                 END IF
              end do
              IF( ALPHA.NE.ONE )THEN
                 do I =  1, M
                    B( I, K ) = ALPHA*B( I, K )
                 end do
              END IF
           end do
        ELSE
           do K =  1, N
              IF( NOUNIT )THEN
                 TEMP = ONE/A( K, K )
                 do I =  1, M
                    B( I, K ) = TEMP*B( I, K )
                 end do
              END IF
              do J =  K + 1, N
                 IF( A( J, K ).NE.ZERO )THEN
                    TEMP = A( J, K )
                    do I =  1, M
                       B( I, J ) = B( I, J ) - TEMP*B( I, K )
                    end do
                 END IF
              end do
              IF( ALPHA.NE.ONE )THEN
                 do I =  1, M
                    B( I, K ) = ALPHA*B( I, K )
                 end do
              END IF
           end do
        END IF
     END IF
  END IF
  !
  RETURN
  !
  !     End of DTRSM .
  !
END SUBROUTINE DTRSM
!============================================================================
integer function idamax(n,dx,incx)
  !
  !     finds the index of element having max. absolute value.
  !     jack dongarra, linpack, 3/11/78.
  !     modified 3/93 to return if incx .le. 0.
  !
  real*8  dx(*),dmax
  integer i,incx,ix,n
  !
  idamax = 0
  if( n.lt.1 .or. incx.le.0 ) return
  idamax = 1
  if(n.eq.1)return

  ix = 1
  dmax = dabs(dx(1))
  ix = ix + incx
  do i =  2,n
     if(.not. (dabs(dx(ix)).le.dmax) ) then 
        idamax = i
        dmax = dabs(dx(ix))
     end if
     ix = ix + incx
  end do
  return
end function idamax
