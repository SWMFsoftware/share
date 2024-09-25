!============================================================================
! This is a collection of single precision BLAS routines that BATSRUS uses. 
! You are encouraged to use the local BLAS library if available.
!
! subroutines: scopy, sgemm, sgemv, sger,  sscal, sswap, strsm
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
subroutine SGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
     BETA, C, LDC )
  !     .. Scalar Arguments ..
  character*1        TRANSA, TRANSB
  integer            M, N, K, LDA, LDB, LDC
  real*4             ALPHA, BETA
  !     .. Array Arguments ..
  real*4             A( LDA, * ), B( LDB, * ), C( LDC, * )
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
  !  parameters
  !  ==========
  !
  !  TRANSA - character*1.
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
  !  TRANSB - character*1.
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
  !  M      - integer.
  !           On entry,  M  specifies  the number  of rows  of the  matrix
  !           op( A )  and of the  matrix  C.  M  must  be at least  zero.
  !           Unchanged on exit.
  !
  !  N      - integer.
  !           On entry,  N  specifies the number  of columns of the matrix
  !           op( B ) and the number of columns of the matrix C. N must be
  !           at least zero.
  !           Unchanged on exit.
  !
  !  K      - integer.
  !           On entry,  K  specifies  the number of columns of the matrix
  !           op( A ) and the number of rows of the matrix op( B ). K must
  !           be at least  zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - real*4          .
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  A      - real*4           array of DIMENSION ( LDA, ka ), where ka is
  !           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
  !           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
  !           part of the array  A  must contain the matrix  A,  otherwise
  !           the leading  k by m  part of the array  A  must contain  the
  !           matrix A.
  !           Unchanged on exit.
  !
  !  LDA    - integer.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
  !           LDA must be at least  max( 1, m ), otherwise  LDA must be at
  !           least  max( 1, k ).
  !           Unchanged on exit.
  !
  !  B      - real*4           array of DIMENSION ( LDB, kb ), where kb is
  !           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
  !           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
  !           part of the array  B  must contain the matrix  B,  otherwise
  !           the leading  n by k  part of the array  B  must contain  the
  !           matrix B.
  !           Unchanged on exit.
  !
  !  LDB    - integer.
  !           On entry, LDB specifies the first dimension of B as declared
  !           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
  !           LDB must be at least  max( 1, k ), otherwise  LDB must be at
  !           least  max( 1, n ).
  !           Unchanged on exit.
  !
  !  BETA   - real*4          .
  !           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
  !           supplied as zero then C need not be set on input.
  !           Unchanged on exit.
  !
  !  C      - real*4           array of DIMENSION ( LDC, n ).
  !           Before entry, the leading  m by n  part of the array  C must
  !           contain the matrix  C,  except when  beta  is zero, in which
  !           case C need not be set on entry.
  !           On exit, the array  C  is overwritten by the  m by n  matrix
  !           ( alpha*op( A )*op( B ) + beta*C ).
  !
  !  LDC    - integer.
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
  !     .. external Functions ..
  logical            LSAME
  external           LSAME
  !     .. external subroutines ..
  external           XERBLA
  !     .. Local Scalars ..
  logical            NOTA, NOTB
  integer            I, INFO, J, L, NCOLA, NROWA, NROWB
  real*4             TEMP
  !     .. parameters ..
  real*4             ONE         , ZERO
  parameter        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
  !     ..
  !     .. Executable Statements ..
  !
  !     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
  !     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
  !     and  columns of  A  and the  number of  rows  of  B  respectively.
  !
  NOTA  = LSAME( TRANSA, 'N' )
  NOTB  = LSAME( TRANSB, 'N' )
  if( NOTA )then
     NROWA = M
     NCOLA = K
  else
     NROWA = K
     NCOLA = M
  end if
  if( NOTB )then
     NROWB = K
  else
     NROWB = N
  end if
  !
  !     Test the input parameters.
  !
  INFO = 0
  if(      ( .not.NOTA                 ).and. &
       ( .not.LSAME( TRANSA, 'C' ) ).and. &
       ( .not.LSAME( TRANSA, 'T' ) )      )then
     INFO = 1
  else if( ( .not.NOTB                 ).and. &
       ( .not.LSAME( TRANSB, 'C' ) ).and. &
       ( .not.LSAME( TRANSB, 'T' ) )      )then
     INFO = 2
  else if( M  .lt.0               )then
     INFO = 3
  else if( N  .lt.0               )then
     INFO = 4
  else if( K  .lt.0               )then
     INFO = 5
  else if( LDA.lt.MAX( 1, NROWA ) )then
     INFO = 8
  else if( LDB.lt.MAX( 1, NROWB ) )then
     INFO = 10
  else if( LDC.lt.MAX( 1, M     ) )then
     INFO = 13
  end if
  if( INFO.ne.0 )then
     CALL XERBLA( 'SGEMM ', INFO )
     RETURN
  end if
  !
  !     Quick return if possible.
  !
  if( ( M.eq.0 ).OR.( N.eq.0 ).OR. &
       ( ( ( ALPHA.eq.ZERO ).OR.( K.eq.0 ) ).and.( BETA.eq.ONE ) ) ) &
       RETURN
  !
  !     And if  alpha.eq.zero.
  !
  if( ALPHA.eq.ZERO )then
     if( BETA.eq.ZERO )then
        do J =  1, N
           do I =  1, M
              C( I, J ) = ZERO
           end do
        end do
     else
        do J =  1, N
           do I =  1, M
              C( I, J ) = BETA*C( I, J )
           end do
        end do
     end if
     RETURN
  end if
  !
  !     Start the operations.
  !
  if( NOTB )then
     if( NOTA )then
        !
        !           Form  C := alpha*A*B + beta*C.
        !
        do J =  1, N
           if( BETA.eq.ZERO )then
              do I =  1, M
                 C( I, J ) = ZERO
              end do
           else if( BETA.ne.ONE )then
              do I =  1, M
                 C( I, J ) = BETA*C( I, J )
              end do
           end if
           do L = 1, K
              if( B( L, J ).ne.ZERO )then
                 TEMP = ALPHA*B( L, J )
                 do I =  1, M
                    C( I, J ) = C( I, J ) + TEMP*A( I, L )
                 end do
              end if
           end do
        end do
     else
        !
        !           Form  C := alpha*A'*B + beta*C
        !
        do J =  1, N
           do I =  1, M
              TEMP = ZERO
              do L = 1, K
                 TEMP = TEMP + A( L, I )*B( L, J )
              end do
              if( BETA.eq.ZERO )then
                 C( I, J ) = ALPHA*TEMP
              else
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              end if
           end do
        end do
     end if
  else
     if( NOTA )then
        !
        !           Form  C := alpha*A*B' + beta*C
        !
        do J =  1, N
           if( BETA.eq.ZERO )then
              do I =  1, M
                 C( I, J ) = ZERO
              end do
           else if( BETA.ne.ONE )then
              do I =  1, M
                 C( I, J ) = BETA*C( I, J )
              end do
           end if
           do L = 1, K
              if( B( J, L ).ne.ZERO )then
                 TEMP = ALPHA*B( J, L )
                 do I =  1, M
                    C( I, J ) = C( I, J ) + TEMP*A( I, L )
                 end do
              end if
           end do
        end do
     else
        !
        !           Form  C := alpha*A'*B' + beta*C
        !
        do J =  1, N
           do I =  1, M
              TEMP = ZERO
              do L = 1, K
                 TEMP = TEMP + A( L, I )*B( J, L )
              end do
              if( BETA.eq.ZERO )then
                 C( I, J ) = ALPHA*TEMP
              else
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              end if
           end do
        end do
     end if
  end if
  !
  RETURN
  !
  !     End of SGEMM .
  !
end subroutine SGEMM
!=============================================================================
subroutine SGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, &
     BETA, Y, INCY )
  !     .. Scalar Arguments ..
  real*4             ALPHA, BETA
  integer            INCX, INCY, LDA, M, N
  character*1        TRANS
  !     .. Array Arguments ..
  real*4             A( LDA, * ), X( * ), Y( * )
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
  !  parameters
  !  ==========
  !
  !  TRANS  - character*1.
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
  !  M      - integer.
  !           On entry, M specifies the number of rows of the matrix A.
  !           M must be at least zero.
  !           Unchanged on exit.
  !
  !  N      - integer.
  !           On entry, N specifies the number of columns of the matrix A.
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - real*4          .
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  A      - real*4           array of DIMENSION ( LDA, n ).
  !           Before entry, the leading m by n part of the array A must
  !           contain the matrix of coefficients.
  !           Unchanged on exit.
  !
  !  LDA    - integer.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. LDA must be at least
  !           max( 1, m ).
  !           Unchanged on exit.
  !
  !  X      - real*4           array of DIMENSION at least
  !           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
  !           and at least
  !           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
  !           Before entry, the incremented array X must contain the
  !           vector x.
  !           Unchanged on exit.
  !
  !  INCX   - integer.
  !           On entry, INCX specifies the increment for the elements of
  !           X. INCX must not be zero.
  !           Unchanged on exit.
  !
  !  BETA   - real*4          .
  !           On entry, BETA specifies the scalar beta. When BETA is
  !           supplied as zero then Y need not be set on input.
  !           Unchanged on exit.
  !
  !  Y      - real*4           array of DIMENSION at least
  !           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
  !           and at least
  !           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
  !           Before entry with BETA non-zero, the incremented array Y
  !           must contain the vector y. On exit, Y is overwritten by the
  !           updated vector y.
  !
  !  INCY   - integer.
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
  !     .. parameters ..
  real*4             ONE         , ZERO
  parameter        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
  !     .. Local Scalars ..
  real*4             TEMP
  integer            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
  !     .. external Functions ..
  logical            LSAME
  external           LSAME
  !     .. external subroutines ..
  external           XERBLA

  INFO = 0
  if     ( .not.LSAME( TRANS, 'N' ).and. &
       .not.LSAME( TRANS, 'T' ).and. &
       .not.LSAME( TRANS, 'C' )      )then
     INFO = 1
  else if( M.lt.0 )then
     INFO = 2
  else if( N.lt.0 )then
     INFO = 3
  else if( LDA.lt.MAX( 1, M ) )then
     INFO = 6
  else if( INCX.eq.0 )then
     INFO = 8
  else if( INCY.eq.0 )then
     INFO = 11
  end if
  if( INFO.ne.0 )then
     CALL XERBLA( 'SGEMV ', INFO )
     RETURN
  end if
  !
  !     Quick return if possible.
  !
  if( ( M.eq.0 ).OR.( N.eq.0 ).OR. &
       ( ( ALPHA.eq.ZERO ).and.( BETA.eq.ONE ) ) ) &
       RETURN
  !
  !     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
  !     up the start points in  X  and  Y.
  !
  if( LSAME( TRANS, 'N' ) )then
     LENX = N
     LENY = M
  else
     LENX = M
     LENY = N
  end if
  if( INCX.GT.0 )then
     KX = 1
  else
     KX = 1 - ( LENX - 1 )*INCX
  end if
  if( INCY.GT.0 )then
     KY = 1
  else
     KY = 1 - ( LENY - 1 )*INCY
  end if
  !
  !     Start the operations. In this version the elements of A are
  !     accessed sequentially with one pass through A.
  !
  !     First form  y := beta*y.
  !
  if( BETA.ne.ONE )then
     if( INCY.eq.1 )then
        if( BETA.eq.ZERO )then
           do I =  1, LENY
              Y( I ) = ZERO
           end do
        else
           do I =  1, LENY
              Y( I ) = BETA*Y( I )
           end do
        end if
     else
        IY = KY
        if( BETA.eq.ZERO )then
           do I =  1, LENY
              Y( IY ) = ZERO
              IY      = IY   + INCY
           end do
        else
           do I =  1, LENY
              Y( IY ) = BETA*Y( IY )
              IY      = IY           + INCY
           end do
        end if
     end if
  end if
  if( ALPHA.eq.ZERO ) &
       RETURN
  if( LSAME( TRANS, 'N' ) )then
     !
     !        Form  y := alpha*A*x + y.
     !
     JX = KX
     if( INCY.eq.1 )then
        do J =  1, N
           if( X( JX ).ne.ZERO )then
              TEMP = ALPHA*X( JX )
              do I =  1, M
                 Y( I ) = Y( I ) + TEMP*A( I, J )
              end do
           end if
           JX = JX + INCX
        end do
     else
        do J =  1, N
           if( X( JX ).ne.ZERO )then
              TEMP = ALPHA*X( JX )
              IY   = KY
              do I =  1, M
                 Y( IY ) = Y( IY ) + TEMP*A( I, J )
                 IY      = IY      + INCY
              end do
           end if
           JX = JX + INCX
        end do
     end if
  else
     !
     !        Form  y := alpha*A'*x + y.
     !
     JY = KY
     if( INCX.eq.1 )then
        do J =  1, N
           TEMP = ZERO
           do I =  1, M
              TEMP = TEMP + A( I, J )*X( I )
           end do
           Y( JY ) = Y( JY ) + ALPHA*TEMP
           JY      = JY      + INCY
        end do
     else
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
     end if
  end if
  !
  RETURN
  !
  !     End of SGEMV .
  !
end subroutine SGEMV
!=============================================================================
subroutine SGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
  !     .. Scalar Arguments ..
  real*4             ALPHA
  integer            INCX, INCY, LDA, M, N
  !     .. Array Arguments ..
  real*4             A( LDA, * ), X( * ), Y( * )
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
  !  parameters
  !  ==========
  !
  !  M      - integer.
  !           On entry, M specifies the number of rows of the matrix A.
  !           M must be at least zero.
  !           Unchanged on exit.
  !
  !  N      - integer.
  !           On entry, N specifies the number of columns of the matrix A.
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - real*4          .
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  X      - real*4           array of dimension at least
  !           ( 1 + ( m - 1 )*abs( INCX ) ).
  !           Before entry, the incremented array X must contain the m
  !           element vector x.
  !           Unchanged on exit.
  !
  !  INCX   - integer.
  !           On entry, INCX specifies the increment for the elements of
  !           X. INCX must not be zero.
  !           Unchanged on exit.
  !
  !  Y      - real*4           array of dimension at least
  !           ( 1 + ( n - 1 )*abs( INCY ) ).
  !           Before entry, the incremented array Y must contain the n
  !           element vector y.
  !           Unchanged on exit.
  !
  !  INCY   - integer.
  !           On entry, INCY specifies the increment for the elements of
  !           Y. INCY must not be zero.
  !           Unchanged on exit.
  !
  !  A      - real*4           array of DIMENSION ( LDA, n ).
  !           Before entry, the leading m by n part of the array A must
  !           contain the matrix of coefficients. On exit, A is
  !           overwritten by the updated matrix.
  !
  !  LDA    - integer.
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
  !     .. parameters ..
  real*4             ZERO
  parameter        ( ZERO = 0.0E+0 )
  !     .. Local Scalars ..
  real*4             TEMP
  integer            I, INFO, IX, J, JY, KX
  !     .. external subroutines ..
  external           XERBLA

  INFO = 0
  if     ( M.lt.0 )then
     INFO = 1
  else if( N.lt.0 )then
     INFO = 2
  else if( INCX.eq.0 )then
     INFO = 5
  else if( INCY.eq.0 )then
     INFO = 7
  else if( LDA.lt.MAX( 1, M ) )then
     INFO = 9
  end if
  if( INFO.ne.0 )then
     CALL XERBLA( 'SGER  ', INFO )
     RETURN
  end if
  !
  !     Quick return if possible.
  !
  if( ( M.eq.0 ).OR.( N.eq.0 ).OR.( ALPHA.eq.ZERO ) ) &
       RETURN
  !
  !     Start the operations. In this version the elements of A are
  !     accessed sequentially with one pass through A.
  !
  if( INCY.GT.0 )then
     JY = 1
  else
     JY = 1 - ( N - 1 )*INCY
  end if
  if( INCX.eq.1 )then
     do J =  1, N
        if( Y( JY ).ne.ZERO )then
           TEMP = ALPHA*Y( JY )
           do I =  1, M
              A( I, J ) = A( I, J ) + X( I )*TEMP
           end do
        end if
        JY = JY + INCY
     end do
  else
     if( INCX.GT.0 )then
        KX = 1
     else
        KX = 1 - ( M - 1 )*INCX
     end if
     do J =  1, N
        if( Y( JY ).ne.ZERO )then
           TEMP = ALPHA*Y( JY )
           IX   = KX
           do I =  1, M
              A( I, J ) = A( I, J ) + X( IX )*TEMP
              IX        = IX        + INCX
           end do
        end if
        JY = JY + INCY
     end do
  end if
  !
  RETURN
  !
  !     End of SGER  .
  !
end subroutine SGER
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
subroutine STRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, &
     B, LDB )
  !     .. Scalar Arguments ..
  character*1        SIDE, UPLO, TRANSA, DIAG
  integer            M, N, LDA, LDB
  real*4             ALPHA
  !     .. Array Arguments ..
  real*4             A( LDA, * ), B( LDB, * )
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
  !  parameters
  !  ==========
  !
  !  SIDE   - character*1.
  !           On entry, SIDE specifies whether op( A ) appears on the left
  !           or right of X as follows:
  !
  !              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
  !
  !              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
  !
  !           Unchanged on exit.
  !
  !  UPLO   - character*1.
  !           On entry, UPLO specifies whether the matrix A is an upper or
  !           lower triangular matrix as follows:
  !
  !              UPLO = 'U' or 'u'   A is an upper triangular matrix.
  !
  !              UPLO = 'L' or 'l'   A is a lower triangular matrix.
  !
  !           Unchanged on exit.
  !
  !  TRANSA - character*1.
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
  !  DIAG   - character*1.
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
  !  M      - integer.
  !           On entry, M specifies the number of rows of B. M must be at
  !           least zero.
  !           Unchanged on exit.
  !
  !  N      - integer.
  !           On entry, N specifies the number of columns of B.  N must be
  !           at least zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - real*4          .
  !           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
  !           zero then  A is not referenced and  B need not be set before
  !           entry.
  !           Unchanged on exit.
  !
  !  A      - real*4           array of DIMENSION ( LDA, k ), where k is m
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
  !  LDA    - integer.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
  !           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
  !           then LDA must be at least max( 1, n ).
  !           Unchanged on exit.
  !
  !  B      - real*4           array of DIMENSION ( LDB, n ).
  !           Before entry,  the leading  m by n part of the array  B must
  !           contain  the  right-hand  side  matrix  B,  and  on exit  is
  !           overwritten by the solution matrix  X.
  !
  !  LDB    - integer.
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
  !     .. external Functions ..
  logical            LSAME
  external           LSAME
  !     .. external subroutines ..
  external           XERBLA

  logical            LSIDE, NOUNIT, UPPER
  integer            I, INFO, J, K, NROWA
  real*4             TEMP
  !     .. parameters ..
  real*4             ONE         , ZERO
  parameter        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  LSIDE  = LSAME( SIDE  , 'L' )
  if( LSIDE )then
     NROWA = M
  else
     NROWA = N
  end if
  NOUNIT = LSAME( DIAG  , 'N' )
  UPPER  = LSAME( UPLO  , 'U' )
  !
  INFO   = 0
  if(      ( .not.LSIDE                ).and. &
       ( .not.LSAME( SIDE  , 'R' ) )      )then
     INFO = 1
  else if( ( .not.UPPER                ).and. &
       ( .not.LSAME( UPLO  , 'L' ) )      )then
     INFO = 2
  else if( ( .not.LSAME( TRANSA, 'N' ) ).and. &
       ( .not.LSAME( TRANSA, 'T' ) ).and. &
       ( .not.LSAME( TRANSA, 'C' ) )      )then
     INFO = 3
  else if( ( .not.LSAME( DIAG  , 'U' ) ).and. &
       ( .not.LSAME( DIAG  , 'N' ) )      )then
     INFO = 4
  else if( M  .lt.0               )then
     INFO = 5
  else if( N  .lt.0               )then
     INFO = 6
  else if( LDA.lt.MAX( 1, NROWA ) )then
     INFO = 9
  else if( LDB.lt.MAX( 1, M     ) )then
     INFO = 11
  end if
  if( INFO.ne.0 )then
     CALL XERBLA( 'STRSM ', INFO )
     RETURN
  end if
  !
  !     Quick return if possible.
  !
  if( N.eq.0 ) &
       RETURN
  !
  !     And when  alpha.eq.zero.
  !
  if( ALPHA.eq.ZERO )then
     do J =  1, N
        do I =  1, M
           B( I, J ) = ZERO
        end do
     end do
     RETURN
  end if
  !
  !     Start the operations.
  !
  if( LSIDE )then
     if( LSAME( TRANSA, 'N' ) )then
        !
        !           Form  B := alpha*inv( A )*B.
        !
        if( UPPER )then
           do J =  1, N
              if( ALPHA.ne.ONE )then
                 do I =  1, M
                    B( I, J ) = ALPHA*B( I, J )
                 end do
              end if
              do K =  M, 1, -1
                 if( B( K, J ).ne.ZERO )then
                    if( NOUNIT ) &
                         B( K, J ) = B( K, J )/A( K, K )
                    do I =  1, K - 1
                       B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
                    end do
                 end if
              end do
           end do
        else
           do J =  1, N
              if( ALPHA.ne.ONE )then
                 do I =  1, M
                    B( I, J ) = ALPHA*B( I, J )
                 end do
              end if
              do K =  1, M
                 if( B( K, J ).ne.ZERO )then
                    if( NOUNIT ) &
                         B( K, J ) = B( K, J )/A( K, K )
                    do I =  K + 1, M
                       B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
                    end do
                 end if
              end do
           end do
        end if
     else
        !
        !           Form  B := alpha*inv( A' )*B.
        !
        if( UPPER )then
           do J =  1, N
              do I =  1, M
                 TEMP = ALPHA*B( I, J )
                 do K =  1, I - 1
                    TEMP = TEMP - A( K, I )*B( K, J )
                 end do
                 if( NOUNIT ) &
                      TEMP = TEMP/A( I, I )
                 B( I, J ) = TEMP
              end do
           end do
        else
           do J =  1, N
              do I =  M, 1, -1
                 TEMP = ALPHA*B( I, J )
                 do K =  I + 1, M
                    TEMP = TEMP - A( K, I )*B( K, J )
                 end do
                 if( NOUNIT ) &
                      TEMP = TEMP/A( I, I )
                 B( I, J ) = TEMP
              end do
           end do
        end if
     end if
  else
     if( LSAME( TRANSA, 'N' ) )then
        !
        !           Form  B := alpha*B*inv( A ).
        !
        if( UPPER )then
           do J =  1, N
              if( ALPHA.ne.ONE )then
                 do I =  1, M
                    B( I, J ) = ALPHA*B( I, J )
                 end do
              end if
              do K =  1, J - 1
                 if( A( K, J ).ne.ZERO )then
                    do I =  1, M
                       B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
                    end do
                 end if
              end do
              if( NOUNIT )then
                 TEMP = ONE/A( J, J )
                 do I =  1, M
                    B( I, J ) = TEMP*B( I, J )
                 end do
              end if
           end do
        else
           do J =  N, 1, -1
              if( ALPHA.ne.ONE )then
                 do I =  1, M
                    B( I, J ) = ALPHA*B( I, J )
                 end do
              end if
              do K =  J + 1, N
                 if( A( K, J ).ne.ZERO )then
                    do I =  1, M
                       B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
                    end do
                 end if
              end do
              if( NOUNIT )then
                 TEMP = ONE/A( J, J )
                 do I =  1, M
                    B( I, J ) = TEMP*B( I, J )
                 end do
              end if
           end do
        end if
     else
        !
        !           Form  B := alpha*B*inv( A' ).
        !
        if( UPPER )then
           do K =  N, 1, -1
              if( NOUNIT )then
                 TEMP = ONE/A( K, K )
                 do I =  1, M
                    B( I, K ) = TEMP*B( I, K )
                 end do
              end if
              do J =  1, K - 1
                 if( A( J, K ).ne.ZERO )then
                    TEMP = A( J, K )
                    do I =  1, M
                       B( I, J ) = B( I, J ) - TEMP*B( I, K )
                    end do
                 end if
              end do
              if( ALPHA.ne.ONE )then
                 do I =  1, M
                    B( I, K ) = ALPHA*B( I, K )
                 end do
              end if
           end do
        else
           do K =  1, N
              if( NOUNIT )then
                 TEMP = ONE/A( K, K )
                 do I =  1, M
                    B( I, K ) = TEMP*B( I, K )
                 end do
              end if
              do J =  K + 1, N
                 if( A( J, K ).ne.ZERO )then
                    TEMP = A( J, K )
                    do I =  1, M
                       B( I, J ) = B( I, J ) - TEMP*B( I, K )
                    end do
                 end if
              end do
              if( ALPHA.ne.ONE )then
                 do I =  1, M
                    B( I, K ) = ALPHA*B( I, K )
                 end do
              end if
           end do
        end if
     end if
  end if
  !
  RETURN
  !
  !     End of STRSM .
  !
end subroutine STRSM
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
! subroutines: dcopy, dgemm, dgemv, dger,  dscal, dswap, dtrsm
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
subroutine DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
     BETA, C, LDC )
  !     .. Scalar Arguments ..
  character*1        TRANSA, TRANSB
  integer            M, N, K, LDA, LDB, LDC
  real*8   ALPHA, BETA
  !     .. Array Arguments ..
  real*8   A( LDA, * ), B( LDB, * ), C( LDC, * )
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
  !  parameters
  !  ==========
  !
  !  TRANSA - character*1.
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
  !  TRANSB - character*1.
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
  !  M      - integer.
  !           On entry,  M  specifies  the number  of rows  of the  matrix
  !           op( A )  and of the  matrix  C.  M  must  be at least  zero.
  !           Unchanged on exit.
  !
  !  N      - integer.
  !           On entry,  N  specifies the number  of columns of the matrix
  !           op( B ) and the number of columns of the matrix C. N must be
  !           at least zero.
  !           Unchanged on exit.
  !
  !  K      - integer.
  !           On entry,  K  specifies  the number of columns of the matrix
  !           op( A ) and the number of rows of the matrix op( B ). K must
  !           be at least  zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - real*8.
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  A      - real*8 array of DIMENSION ( LDA, ka ), where ka is
  !           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
  !           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
  !           part of the array  A  must contain the matrix  A,  otherwise
  !           the leading  k by m  part of the array  A  must contain  the
  !           matrix A.
  !           Unchanged on exit.
  !
  !  LDA    - integer.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
  !           LDA must be at least  max( 1, m ), otherwise  LDA must be at
  !           least  max( 1, k ).
  !           Unchanged on exit.
  !
  !  B      - real*8 array of DIMENSION ( LDB, kb ), where kb is
  !           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
  !           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
  !           part of the array  B  must contain the matrix  B,  otherwise
  !           the leading  n by k  part of the array  B  must contain  the
  !           matrix B.
  !           Unchanged on exit.
  !
  !  LDB    - integer.
  !           On entry, LDB specifies the first dimension of B as declared
  !           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
  !           LDB must be at least  max( 1, k ), otherwise  LDB must be at
  !           least  max( 1, n ).
  !           Unchanged on exit.
  !
  !  BETA   - real*8.
  !           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
  !           supplied as zero then C need not be set on input.
  !           Unchanged on exit.
  !
  !  C      - real*8 array of DIMENSION ( LDC, n ).
  !           Before entry, the leading  m by n  part of the array  C must
  !           contain the matrix  C,  except when  beta  is zero, in which
  !           case C need not be set on entry.
  !           On exit, the array  C  is overwritten by the  m by n  matrix
  !           ( alpha*op( A )*op( B ) + beta*C ).
  !
  !  LDC    - integer.
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
  !     .. external Functions ..
  logical            LSAME
  external           LSAME
  !     .. external subroutines ..
  external           XERBLA
  !     .. Local Scalars ..
  logical            NOTA, NOTB
  integer            I, INFO, J, L, NCOLA, NROWA, NROWB
  real*8   TEMP
  !     .. parameters ..
  real*8   ONE    , ZERO
  parameter        ( ONE = 1, ZERO = 0 )
  !     ..
  !     .. Executable Statements ..
  !
  !     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
  !     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
  !     and  columns of  A  and the  number of  rows  of  B  respectively.
  !
  NOTA  = LSAME( TRANSA, 'N' )
  NOTB  = LSAME( TRANSB, 'N' )
  if( NOTA )then
     NROWA = M
     NCOLA = K
  else
     NROWA = K
     NCOLA = M
  end if
  if( NOTB )then
     NROWB = K
  else
     NROWB = N
  end if
  !
  !     Test the input parameters.
  !
  INFO = 0
  if(      ( .not.NOTA                 ).and. &
       ( .not.LSAME( TRANSA, 'C' ) ).and. &
       ( .not.LSAME( TRANSA, 'T' ) )      )then
     INFO = 1
  else if( ( .not.NOTB                 ).and. &
       ( .not.LSAME( TRANSB, 'C' ) ).and. &
       ( .not.LSAME( TRANSB, 'T' ) )      )then
     INFO = 2
  else if( M  .lt.0               )then
     INFO = 3
  else if( N  .lt.0               )then
     INFO = 4
  else if( K  .lt.0               )then
     INFO = 5
  else if( LDA.lt.MAX( 1, NROWA ) )then
     INFO = 8
  else if( LDB.lt.MAX( 1, NROWB ) )then
     INFO = 10
  else if( LDC.lt.MAX( 1, M     ) )then
     INFO = 13
  end if
  if( INFO.ne.0 )then
     CALL XERBLA( 'DGEMM ', INFO )
     RETURN
  end if
  !
  !     Quick return if possible.
  !
  if( ( M.eq.0 ).OR.( N.eq.0 ).OR. &
       ( ( ( ALPHA.eq.ZERO ).OR.( K.eq.0 ) ).and.( BETA.eq.ONE ) ) ) &
       RETURN
  !
  !     And if  alpha.eq.zero.
  !
  if( ALPHA.eq.ZERO )then
     if( BETA.eq.ZERO )then
        do J =  1, N
           do I =  1, M
              C( I, J ) = ZERO
           end do
        end do
     else
        do J =  1, N
           do I =  1, M
              C( I, J ) = BETA*C( I, J )
           end do
        end do
     end if
     RETURN
  end if
  !
  !     Start the operations.
  !
  if( NOTB )then
     if( NOTA )then
        !
        !           Form  C := alpha*A*B + beta*C.
        !
        do J =  1, N
           if( BETA.eq.ZERO )then
              do I =  1, M
                 C( I, J ) = ZERO
              end do
           else if( BETA.ne.ONE )then
              do I =  1, M
                 C( I, J ) = BETA*C( I, J )
              end do
           end if
           do  L = 1, K
              if( B( L, J ).ne.ZERO )then
                 TEMP = ALPHA*B( L, J )
                 do I =  1, M
                    C( I, J ) = C( I, J ) + TEMP*A( I, L )
                 end do
              end if
           end do
        end do
     else
        !
        !           Form  C := alpha*A'*B + beta*C
        !
        do J =  1, N
           do I =  1, M
              TEMP = ZERO
              do L = 1, K
                 TEMP = TEMP + A( L, I )*B( L, J )
              end do
              if( BETA.eq.ZERO )then
                 C( I, J ) = ALPHA*TEMP
              else
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              end if
           end do
        end do
     end if
  else
     if( NOTA )then
        !
        !           Form  C := alpha*A*B' + beta*C
        !
        do J =  1, N
           if( BETA.eq.ZERO )then
              do I =  1, M
                 C( I, J ) = ZERO
              end do
           else if( BETA.ne.ONE )then
              do I =  1, M
                 C( I, J ) = BETA*C( I, J )
              end do
           end if
           do L = 1, K
              if( B( J, L ).ne.ZERO )then
                 TEMP = ALPHA*B( J, L )
                 do I =  1, M
                    C( I, J ) = C( I, J ) + TEMP*A( I, L )
                 end do
              end if
           end do
        end do
     else
        !
        !           Form  C := alpha*A'*B' + beta*C
        !
        do J =  1, N
           do I =  1, M
              TEMP = ZERO
              do L = 1, K
                 TEMP = TEMP + A( L, I )*B( J, L )
              end do
              if( BETA.eq.ZERO )then
                 C( I, J ) = ALPHA*TEMP
              else
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              end if
           end do
        end do
     end if
  end if
  !
  RETURN
  !
  !     End of DGEMM .
  !
end subroutine DGEMM
!============================================================================
subroutine DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, &
     BETA, Y, INCY )
  !     .. Scalar Arguments ..
  real*8   ALPHA, BETA
  integer            INCX, INCY, LDA, M, N
  character*1        TRANS
  !     .. Array Arguments ..
  real*8   A( LDA, * ), X( * ), Y( * )
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
  !  parameters
  !  ==========
  !
  !  TRANS  - character*1.
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
  !  M      - integer.
  !           On entry, M specifies the number of rows of the matrix A.
  !           M must be at least zero.
  !           Unchanged on exit.
  !
  !  N      - integer.
  !           On entry, N specifies the number of columns of the matrix A.
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - real*8.
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  A      - real*8 array of DIMENSION ( LDA, n ).
  !           Before entry, the leading m by n part of the array A must
  !           contain the matrix of coefficients.
  !           Unchanged on exit.
  !
  !  LDA    - integer.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. LDA must be at least
  !           max( 1, m ).
  !           Unchanged on exit.
  !
  !  X      - real*8 array of DIMENSION at least
  !           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
  !           and at least
  !           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
  !           Before entry, the incremented array X must contain the
  !           vector x.
  !           Unchanged on exit.
  !
  !  INCX   - integer.
  !           On entry, INCX specifies the increment for the elements of
  !           X. INCX must not be zero.
  !           Unchanged on exit.
  !
  !  BETA   - real*8.
  !           On entry, BETA specifies the scalar beta. When BETA is
  !           supplied as zero then Y need not be set on input.
  !           Unchanged on exit.
  !
  !  Y      - real*8 array of DIMENSION at least
  !           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
  !           and at least
  !           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
  !           Before entry with BETA non-zero, the incremented array Y
  !           must contain the vector y. On exit, Y is overwritten by the
  !           updated vector y.
  !
  !  INCY   - integer.
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
  !     .. parameters ..
  real*8   ONE         , ZERO
  parameter        ( ONE = 1, ZERO = 0 )
  !     .. Local Scalars ..
  real*8   TEMP
  integer            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
  !     .. external Functions ..
  logical            LSAME
  external           LSAME
  !     .. external subroutines ..
  external           XERBLA

  INFO = 0
  if     ( .not.LSAME( TRANS, 'N' ).and. &
       .not.LSAME( TRANS, 'T' ).and. &
       .not.LSAME( TRANS, 'C' )      )then
     INFO = 1
  else if( M.lt.0 )then
     INFO = 2
  else if( N.lt.0 )then
     INFO = 3
  else if( LDA.lt.MAX( 1, M ) )then
     INFO = 6
  else if( INCX.eq.0 )then
     INFO = 8
  else if( INCY.eq.0 )then
     INFO = 11
  end if
  if( INFO.ne.0 )then
     CALL XERBLA( 'DGEMV ', INFO )
     RETURN
  end if
  !
  !     Quick return if possible.
  !
  if( ( M.eq.0 ).OR.( N.eq.0 ).OR. &
       ( ( ALPHA.eq.ZERO ).and.( BETA.eq.ONE ) ) ) &
       RETURN
  !
  !     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
  !     up the start points in  X  and  Y.
  !
  if( LSAME( TRANS, 'N' ) )then
     LENX = N
     LENY = M
  else
     LENX = M
     LENY = N
  end if
  if( INCX.GT.0 )then
     KX = 1
  else
     KX = 1 - ( LENX - 1 )*INCX
  end if
  if( INCY.GT.0 )then
     KY = 1
  else
     KY = 1 - ( LENY - 1 )*INCY
  end if
  !
  !     Start the operations. In this version the elements of A are
  !     accessed sequentially with one pass through A.
  !
  !     First form  y := beta*y.
  !
  if( BETA.ne.ONE )then
     if( INCY.eq.1 )then
        if( BETA.eq.ZERO )then
           do I =  1, LENY
              Y( I ) = ZERO
           end do
        else
           do I =  1, LENY
              Y( I ) = BETA*Y( I )
           end do
        end if
     else
        IY = KY
        if( BETA.eq.ZERO )then
           do I =  1, LENY
              Y( IY ) = ZERO
              IY      = IY   + INCY
           end do
        else
           do I =  1, LENY
              Y( IY ) = BETA*Y( IY )
              IY      = IY           + INCY
           end do
        end if
     end if
  end if
  if( ALPHA.eq.ZERO ) &
       RETURN
  if( LSAME( TRANS, 'N' ) )then
     !
     !        Form  y := alpha*A*x + y.
     !
     JX = KX
     if( INCY.eq.1 )then
        do J =  1, N
           if( X( JX ).ne.ZERO )then
              TEMP = ALPHA*X( JX )
              do I =  1, M
                 Y( I ) = Y( I ) + TEMP*A( I, J )
              end do
           end if
           JX = JX + INCX
        end do
     else
        do J =  1, N
           if( X( JX ).ne.ZERO )then
              TEMP = ALPHA*X( JX )
              IY   = KY
              do I =  1, M
                 Y( IY ) = Y( IY ) + TEMP*A( I, J )
                 IY      = IY      + INCY
              end do
           end if
           JX = JX + INCX
        end do
     end if
  else
     !
     !        Form  y := alpha*A'*x + y.
     !
     JY = KY
     if( INCX.eq.1 )then
        do J =  1, N
           TEMP = ZERO
           do I =  1, M
              TEMP = TEMP + A( I, J )*X( I )
           end do
           Y( JY ) = Y( JY ) + ALPHA*TEMP
           JY      = JY      + INCY
        end do
     else
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
     end if
  end if
  !
  RETURN
  !
  !     End of DGEMV .
  !
end subroutine DGEMV
!============================================================================
subroutine DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
  !     .. Scalar Arguments ..
  real*8   ALPHA
  integer            INCX, INCY, LDA, M, N
  !     .. Array Arguments ..
  real*8   A( LDA, * ), X( * ), Y( * )
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
  !  parameters
  !  ==========
  !
  !  M      - integer.
  !           On entry, M specifies the number of rows of the matrix A.
  !           M must be at least zero.
  !           Unchanged on exit.
  !
  !  N      - integer.
  !           On entry, N specifies the number of columns of the matrix A.
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - real*8.
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  X      - real*8 array of dimension at least
  !           ( 1 + ( m - 1 )*abs( INCX ) ).
  !           Before entry, the incremented array X must contain the m
  !           element vector x.
  !           Unchanged on exit.
  !
  !  INCX   - integer.
  !           On entry, INCX specifies the increment for the elements of
  !           X. INCX must not be zero.
  !           Unchanged on exit.
  !
  !  Y      - real*8 array of dimension at least
  !           ( 1 + ( n - 1 )*abs( INCY ) ).
  !           Before entry, the incremented array Y must contain the n
  !           element vector y.
  !           Unchanged on exit.
  !
  !  INCY   - integer.
  !           On entry, INCY specifies the increment for the elements of
  !           Y. INCY must not be zero.
  !           Unchanged on exit.
  !
  !  A      - real*8 array of DIMENSION ( LDA, n ).
  !           Before entry, the leading m by n part of the array A must
  !           contain the matrix of coefficients. On exit, A is
  !           overwritten by the updated matrix.
  !
  !  LDA    - integer.
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
  !     .. parameters ..
  real*8   ZERO
  parameter        ( ZERO = 0.0D+0 )
  !     .. Local Scalars ..
  real*8   TEMP
  integer            I, INFO, IX, J, JY, KX
  !     .. external subroutines ..
  external           XERBLA

  INFO = 0
  if     ( M.lt.0 )then
     INFO = 1
  else if( N.lt.0 )then
     INFO = 2
  else if( INCX.eq.0 )then
     INFO = 5
  else if( INCY.eq.0 )then
     INFO = 7
  else if( LDA.lt.MAX( 1, M ) )then
     INFO = 9
  end if
  if( INFO.ne.0 )then
     CALL XERBLA( 'DGER  ', INFO )
     RETURN
  end if
  !
  !     Quick return if possible.
  !
  if( ( M.eq.0 ).OR.( N.eq.0 ).OR.( ALPHA.eq.ZERO ) ) &
       RETURN
  !
  !     Start the operations. In this version the elements of A are
  !     accessed sequentially with one pass through A.
  !
  if( INCY.GT.0 )then
     JY = 1
  else
     JY = 1 - ( N - 1 )*INCY
  end if
  if( INCX.eq.1 )then
     do J =  1, N
        if( Y( JY ).ne.ZERO )then
           TEMP = ALPHA*Y( JY )
           do I =  1, M
              A( I, J ) = A( I, J ) + X( I )*TEMP
           end do
        end if
        JY = JY + INCY
     end do
  else
     if( INCX.GT.0 )then
        KX = 1
     else
        KX = 1 - ( M - 1 )*INCX
     end if
     do J =  1, N
        if( Y( JY ).ne.ZERO )then
           TEMP = ALPHA*Y( JY )
           IX   = KX
           do I =  1, M
              A( I, J ) = A( I, J ) + X( IX )*TEMP
              IX        = IX        + INCX
           end do
        end if
        JY = JY + INCY
     end do
  end if
  !
  RETURN
  !
  !     End of DGER  .
  !
end subroutine DGER
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
subroutine DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, &
     B, LDB )
  !     .. Scalar Arguments ..
  character*1        SIDE, UPLO, TRANSA, DIAG
  integer            M, N, LDA, LDB
  real*8   ALPHA
  !     .. Array Arguments ..
  real*8   A( LDA, * ), B( LDB, * )
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
  !  parameters
  !  ==========
  !
  !  SIDE   - character*1.
  !           On entry, SIDE specifies whether op( A ) appears on the left
  !           or right of X as follows:
  !
  !              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
  !
  !              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
  !
  !           Unchanged on exit.
  !
  !  UPLO   - character*1.
  !           On entry, UPLO specifies whether the matrix A is an upper or
  !           lower triangular matrix as follows:
  !
  !              UPLO = 'U' or 'u'   A is an upper triangular matrix.
  !
  !              UPLO = 'L' or 'l'   A is a lower triangular matrix.
  !
  !           Unchanged on exit.
  !
  !  TRANSA - character*1.
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
  !  DIAG   - character*1.
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
  !  M      - integer.
  !           On entry, M specifies the number of rows of B. M must be at
  !           least zero.
  !           Unchanged on exit.
  !
  !  N      - integer.
  !           On entry, N specifies the number of columns of B.  N must be
  !           at least zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - real*8.
  !           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
  !           zero then  A is not referenced and  B need not be set before
  !           entry.
  !           Unchanged on exit.
  !
  !  A      - real*8 array of DIMENSION ( LDA, k ), where k is m
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
  !  LDA    - integer.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
  !           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
  !           then LDA must be at least max( 1, n ).
  !           Unchanged on exit.
  !
  !  B      - real*8 array of DIMENSION ( LDB, n ).
  !           Before entry,  the leading  m by n part of the array  B must
  !           contain  the  right-hand  side  matrix  B,  and  on exit  is
  !           overwritten by the solution matrix  X.
  !
  !  LDB    - integer.
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
  !     .. external Functions ..
  logical            LSAME
  external           LSAME
  !     .. external subroutines ..
  external           XERBLA
  !     .. Local Scalars ..
  logical            LSIDE, NOUNIT, UPPER
  integer            I, INFO, J, K, NROWA
  real*8             TEMP
  !     .. parameters ..
  real*8             ONE    , ZERO
  parameter        ( ONE = 1, ZERO = 0 )
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  LSIDE  = LSAME( SIDE  , 'L' )
  if( LSIDE )then
     NROWA = M
  else
     NROWA = N
  end if
  NOUNIT = LSAME( DIAG  , 'N' )
  UPPER  = LSAME( UPLO  , 'U' )
  !
  INFO   = 0
  if(      ( .not.LSIDE                ).and. &
       ( .not.LSAME( SIDE  , 'R' ) )      )then
     INFO = 1
  else if( ( .not.UPPER                ).and. &
       ( .not.LSAME( UPLO  , 'L' ) )      )then
     INFO = 2
  else if( ( .not.LSAME( TRANSA, 'N' ) ).and. &
       ( .not.LSAME( TRANSA, 'T' ) ).and. &
       ( .not.LSAME( TRANSA, 'C' ) )      )then
     INFO = 3
  else if( ( .not.LSAME( DIAG  , 'U' ) ).and. &
       ( .not.LSAME( DIAG  , 'N' ) )      )then
     INFO = 4
  else if( M  .lt.0               )then
     INFO = 5
  else if( N  .lt.0               )then
     INFO = 6
  else if( LDA.lt.MAX( 1, NROWA ) )then
     INFO = 9
  else if( LDB.lt.MAX( 1, M     ) )then
     INFO = 11
  end if
  if( INFO.ne.0 )then
     CALL XERBLA( 'DTRSM ', INFO )
     RETURN
  end if
  !
  !     Quick return if possible.
  !
  if( N.eq.0 ) &
       RETURN
  !
  !     And when  alpha.eq.zero.
  !
  if( ALPHA.eq.ZERO )then
     do J =  1, N
        do I =  1, M
           B( I, J ) = ZERO
        end do
     end do
     RETURN
  end if
  !
  !     Start the operations.
  !
  if( LSIDE )then
     if( LSAME( TRANSA, 'N' ) )then
        !
        !           Form  B := alpha*inv( A )*B.
        !
        if( UPPER )then
           do J =  1, N
              if( ALPHA.ne.ONE )then
                 do I =  1, M
                    B( I, J ) = ALPHA*B( I, J )
                 end do
              end if
              do K =  M, 1, -1
                 if( B( K, J ).ne.ZERO )then
                    if( NOUNIT ) &
                         B( K, J ) = B( K, J )/A( K, K )
                    do I =  1, K - 1
                       B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
                    end do
                 end if
              end do
           end do
        else
           do J =  1, N
              if( ALPHA.ne.ONE )then
                 do I =  1, M
                    B( I, J ) = ALPHA*B( I, J )
                 end do
              end if
              do K =  1, M
                 if( B( K, J ).ne.ZERO )then
                    if( NOUNIT ) &
                         B( K, J ) = B( K, J )/A( K, K )
                    do I =  K + 1, M
                       B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
                    end do
                 end if
              end do
           end do
        end if
     else
        !
        !           Form  B := alpha*inv( A' )*B.
        !
        if( UPPER )then
           do J =  1, N
              do I =  1, M
                 TEMP = ALPHA*B( I, J )
                 do K =  1, I - 1
                    TEMP = TEMP - A( K, I )*B( K, J )
                 end do
                 if( NOUNIT ) &
                      TEMP = TEMP/A( I, I )
                 B( I, J ) = TEMP
              end do
           end do
        else
           do J =  1, N
              do I =  M, 1, -1
                 TEMP = ALPHA*B( I, J )
                 do K =  I + 1, M
                    TEMP = TEMP - A( K, I )*B( K, J )
                 end do
                 if( NOUNIT ) &
                      TEMP = TEMP/A( I, I )
                 B( I, J ) = TEMP
              end do
           end do
        end if
     end if
  else
     if( LSAME( TRANSA, 'N' ) )then
        !
        !           Form  B := alpha*B*inv( A ).
        !
        if( UPPER )then
           do J =  1, N
              if( ALPHA.ne.ONE )then
                 do I =  1, M
                    B( I, J ) = ALPHA*B( I, J )
                 end do
              end if
              do K =  1, J - 1
                 if( A( K, J ).ne.ZERO )then
                    do I =  1, M
                       B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
                    end do
                 end if
              end do
              if( NOUNIT )then
                 TEMP = ONE/A( J, J )
                 do I =  1, M
                    B( I, J ) = TEMP*B( I, J )
                 end do
              end if
           end do
        else
           do J =  N, 1, -1
              if( ALPHA.ne.ONE )then
                 do I =  1, M
                    B( I, J ) = ALPHA*B( I, J )
                 end do
              end if
              do K =  J + 1, N
                 if( A( K, J ).ne.ZERO )then
                    do I =  1, M
                       B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
                    end do
                 end if
              end do
              if( NOUNIT )then
                 TEMP = ONE/A( J, J )
                 do I =  1, M
                    B( I, J ) = TEMP*B( I, J )
                 end do
              end if
           end do
        end if
     else
        !
        !           Form  B := alpha*B*inv( A' ).
        !
        if( UPPER )then
           do K =  N, 1, -1
              if( NOUNIT )then
                 TEMP = ONE/A( K, K )
                 do I =  1, M
                    B( I, K ) = TEMP*B( I, K )
                 end do
              end if
              do J =  1, K - 1
                 if( A( J, K ).ne.ZERO )then
                    TEMP = A( J, K )
                    do I =  1, M
                       B( I, J ) = B( I, J ) - TEMP*B( I, K )
                    end do
                 end if
              end do
              if( ALPHA.ne.ONE )then
                 do I =  1, M
                    B( I, K ) = ALPHA*B( I, K )
                 end do
              end if
           end do
        else
           do K =  1, N
              if( NOUNIT )then
                 TEMP = ONE/A( K, K )
                 do I =  1, M
                    B( I, K ) = TEMP*B( I, K )
                 end do
              end if
              do J =  K + 1, N
                 if( A( J, K ).ne.ZERO )then
                    TEMP = A( J, K )
                    do I =  1, M
                       B( I, J ) = B( I, J ) - TEMP*B( I, K )
                    end do
                 end if
              end do
              if( ALPHA.ne.ONE )then
                 do I =  1, M
                    B( I, K ) = ALPHA*B( I, K )
                 end do
              end if
           end do
        end if
     end if
  end if
  !
  RETURN
  !
  !     End of DTRSM .
  !
end subroutine DTRSM
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
