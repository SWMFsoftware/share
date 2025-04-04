!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModLinearSolver

  ! Various methods to solve linear system of equations.
  ! There are both serial and parallel solvers, and direct and
  ! iterative solvers.

  ! GENERAL PRECONDITIONER FOR BLOCK HEPTADIAGONAL AND PENTADIAGONAL MATRICES.
  ! DIRECT SOLVER FOR BLOCK TRIDIAGONAL MATRIX.
  !
  ! G. Toth 2001
  ! (based on the original F77 block heptadiagonal preconditioner with
  !  Eisenstat trick implemented by Auke van der Ploeg, 1997)
  !
  ! subroutines:
  !          prehepta, Uhepta, Lhepta
  !
  ! Usage: call prehepta with the original matrix to obtain L and U.
  !        call Uhepta(.false.) to multiply a vector with U
  !        call Uhepta(.true.)  to multiply a vector with U^{-1}
  !        call Lhepta          to multiply a vector with L^{-1}
  !
  ! To solve a tridiagonal system A.x=rhs use these steps
  !
  !        call prehepta(A)        ! A --> LU
  !        x=rhs
  !        call Lhepta(x)          ! x --> L^{-1}.rhs
  !        call Uhepta(.true.,x)   ! x --> U^{-1}.L^{-1}.rhs = A^{-1}.rhs
  !
  ! To solve a penta- or hepta-diagonal problem with symmetric preconditioning
  !
  ! L^{-1}.A.U^{-1} U.x = L^{-1}.rhs
  !
  ! use the following steps:
  !
  !        call prehepta(A)                 ! A -> LU
  !        call Lhepta(rhs)                 ! rhs'=L^{-1}.rhs
  !        call Uhepta(.false.,x)           ! x'  = U.x (for initial guess)
  !        call bicgstab(matvec_prec,x,rhs) ! solve A'.x'=rhs'
  !        call Uhepta(.true.,x)            ! x = U^{-1}.x'
  !
  ! The preconditioned matrix vector multiplication is in
  ! subroutine matvec_prec(x,y,n), and it should calculate
  ! y = A'.x = L^{-1}.A.U^{-1}.x as
  !
  !        y=x
  !        call Uhepta(.true.,y)   ! multiply y with U^{-1}
  !        call matvec(y)          ! multiply y with A
  !        call Lhepta(y)          ! multiply y with L^{-1}

  use ModMpi
  use ModUtilities, ONLY: CON_stop, CON_stop_simple
  use ModBlasLapack, ONLY: BLAS_gemm, BLAS_copy, BLAS_gemv, &
       LAPACK_getrf, LAPACK_getrs
  use omp_lib

  ! BLAS_gemm, level 3 Matrix-Matrix Product.
  !
  ! BLAS_gemv, level 2 Matrix-Vector Product.
  !
  ! BLAS_copy, level 1 vector copy.
  !
  ! LAPACK_getrf, computes an LU factorization of a general
  !          M-by-N matrix A using partial pivoting with
  !          row interchanges. The factorization has the form
  !              A = P * L * U
  !          where P is a permutation matrix, L is lower triangular
  !          with unit diagonal elements (lower trapezoidal if m > n),
  !          and U is upper triangular (upper trapezoidal if m < n).
  !
  ! LAPACK_getrs, solves a system of linear equations
  !          A * X = B,  A**T * X = B,  or  A**H * X = B
  !          with a general N-by-N matrix A using the LU factorization
  !          computed by LAPACK_getrf

  implicit none
  SAVE

  private ! except

  public :: gmres           ! GMRES iterative solver
  public :: bicgstab        ! BiCGSTAB iterative solver
  public :: cg              ! CG iterative solver for symmetric positive matrix

  public :: prehepta        ! LU preconditioner for up to hepta block-diagonal
  public :: Uhepta          ! multiply with upper block triangular matrix
  public :: Lhepta          ! multiply with lower block triangular matrix
  public :: upper_hepta_scalar ! multiply with upper scalar triangular matrix
  public :: lower_hepta_scalar ! multiply with lower scalar triangular matrix
  public :: get_precond_matrix     ! get preconditioner matrix
  public :: multiply_left_precond  ! multiply with left preconditioner
  public :: multiply_right_precond ! multiply with right preconditioner
  public :: multiply_initial_guess ! multiply with initial guess with U

  public :: multiply_dilu         ! multiply with (LU)^-1
  !                                 assuming diagonal off-diag blocks
  public :: multiply_block_jacobi ! multiply with inverted diag blocks D^-1

  public :: implicit_solver       ! implicit solver in 1D with 3 point stencil

  public :: solve_linear_multiblock  ! solver for multiblock grid
  public :: precond_left_multiblock  ! left precond for multiblock grid
  public :: precond_right_multiblock ! right precond for multiblock grid

  public :: test_linear_solver    ! unit test

  ! To accurately calculate the quadratic form p.A.p in conjugate gradients
  real,    public:: pDotADotPPe = 0.0
  logical, public:: UsePDotADotP = .false.

  ! Use an effectively 16byte real accuracy for global sums
  logical, public:: UseAccurateSum = .false.

  ! Named indexes for various preconditioner options
  integer, public, parameter:: &
       Jacobi_=6, BlockJacobi_=5, GaussSeidel_=4, Dilu_=3, &
       Bilu1_=2, Bilu_=1, Mbilu_=0

  type LinearSolverParamType
     logical          :: DoPrecond        ! Do preconditioning
     character(len=10):: TypePrecondSide  ! Precondition left, right, symmetric
     character(len=11):: TypePrecond      ! Preconditioner type
     real             :: PrecondParam     ! Parameter (mostly for MBILU)
     character(len=10):: TypeKrylov       ! Krylov solver type
     character(len=3) :: TypeStop         ! Stopping criterion type (rel,abs)
     real             :: ErrorMax         ! Tolerance for solver
     integer          :: MaxMatvec        ! Maximum number of iterations
     integer          :: nKrylovVector    ! Number of vectors for GMRES
     logical          :: UseInitialGuess  ! non-zero initial guess
     real             :: Error            ! Actual accuracy achieved
     integer          :: nMatvec          ! Actual number of iterations
     integer          :: iError           ! Error code from the solver
  end type LinearSolverParamType

  public:: LinearSolverParamType

  ! revision history:
  ! 05Dec06 - Gabor Toth - initial prototype/prolog/code based on BATSRUS code
  ! 20Mar14 - Gabor Toth - lot of new code

  ! Used for tests
  integer, parameter :: rho_=1, rhou_=2, p_=3
  real, parameter :: Gamma=1.6, Dx=1.0, Inv2Dx=0.5/Dx
  integer, parameter  :: UnitTmp_=9
  integer :: iStep
  real    :: Time

  ! Number of unknowns
  integer:: n0 = 0
  ! Size of Krylov subspace
  integer:: nKrylov0 = 0
  ! Krylov subspace vectors
  real, allocatable :: Krylov_II(:,:)
  !$acc declare create(Krylov_II)

  ! Hessenberg matrix and some vectors
  real, allocatable :: hh(:,:)
  !$acc declare create(hh)
  real, allocatable :: c(:), s(:), rs(:)
  !$acc declare create(c, s, rs)
contains
  !============================================================================
  subroutine gmres(matvec, Rhs, Sol, IsInit, n, nKrylov, Tol, TypeStop, &
       Iter, info, DoTest, iCommIn)

    ! This code was initially written by Youcef Saad (May 23, 1985)
    ! then revised by Henk A. van der Vorst and Mike Botchev (Oct. 1996)
    ! Rewritten into F90 and parallelized for the BATSRUS code (May 2002)
    ! by Gabor Toth
    ! Moved into ModLinearSolver.f90 for SWMF by Gabor Toth (December 2006)

    ! subroutine for matrix vector multiplication
    interface
       subroutine matvec(a, b, n)
         implicit none
         ! Calculate b = M.a where M is the matrix
         integer, intent(in) :: n
         real, intent(in) ::  a(n)
         real, intent(out) :: b(n)
       end subroutine matvec
    end interface

    integer, intent(in) :: n         ! number of unknowns.
    integer, intent(in) :: nKrylov   ! size of krylov subspace
    real,    intent(in) :: Rhs(n)    ! right hand side vector
    real, intent(inout) :: Sol(n)    ! initial guess / solution vector
    logical, intent(in) :: IsInit    ! true if Sol contains initial guess
    real, intent(inout) :: Tol       ! required / achieved residual

    !        on input : required (relative) 2-norm or maximum norm of residual
    !        on output: achieved (relative) 2-norm or maximum norm of residual
    ! eps    == tolerance for stopping criterion. process is stopped
    !           as soon as ( ||.|| is the euclidean norm):
    !           || current residual||/||initial residual|| <= eps
    !           on OUTPUT: actual achieved norm residual (if iabs.ne.0)
    !           or achieved relative residual norm reduction

    character (len=3), intent(in) :: TypeStop
    !      Determine stopping criterion (||.|| denotes the 2-norm):
    !      typestop='rel'    -- relative stopping crit.:||res|| <= Tol*||res0||
    !      typestop='abs'    -- absolute stopping crit.: ||res|| <= Tol
    !      typestop='max'    -- maximum  stopping crit.: max(abs(res)) <= Tol

    integer, intent(inout) :: Iter    ! maximum/actual number of iterations

    integer, intent(out)   :: info    ! gives reason for returning:
    !     abs(info)=  0 - solution found satisfying given tolerance.
    !                 2 - no convergence within maximum number of iterations.
    !                 3 - initial guess satisfies the stopping criterion.
    !    sign(info)=  + - residual decreased
    !                 - - residual did not reduce

    logical, intent(in)    :: DoTest  !  write debug info if true
    integer, intent(in), optional :: iCommIn   ! MPI communicator

    ! Local variables
    integer :: iComm                           ! MPI communicator
    integer :: i, i1, its, j, k, k1
    real :: coeff, Tol1, epsmac, gam, ro, ro0, t, tmp
    !$acc declare create(coeff, ro, t, tmp)
    !--------------------------------------------------------------------------
    if(DoTest)write(*,*)'GMRES tol,iter:',Tol,Iter
    ! call timing_start('gmres')

    ! Assign the MPI communicator
    iComm = MPI_COMM_SELF
    if(present(iCommIn)) iComm = iCommIn

    ! If the size of arrays changes, dellocate the arrays
    if((n /= n0 .or. nKrylov /= nKrylov0) .and. allocated(Krylov_II)) then
       deallocate(Krylov_II, hh, c, s, rs)
    end if

    if(.not.allocated(Krylov_II)) then
       allocate(Krylov_II(n,nKrylov+2))
       allocate(hh(nKrylov+1,nKrylov))
       allocate(c(nKrylov), s(nKrylov), rs(nKrylov+1))
    end if

    if(range(1.0)>100)then
       epsmac=0.0000000000000001
    else
       epsmac=0.00000001
    endif

    its = 0

    ! **  outer loop starts here..
    !-------------- compute initial residual vector --------------

    RESTARTLOOP: do
       !
       !           Krylov_II(1):=A*Sol
       !
       if(IsInit .or. its>0)then
          call matvec(Sol, Krylov_II, n)
          !$acc parallel loop
          do i = 1, n
             Krylov_II(i,1) = Rhs(i) - Krylov_II(i,1)
          end do
       else
          ! Save a matvec when starting from zero initial condition
          !$acc parallel loop
          do i = 1, n
             Krylov_II(i,1) = Rhs(i)
          end do
       endif

       ro = sqrt( dot_product_mpi(n, Krylov_II(:,1), Krylov_II(:,1), iComm ))
       !$acc update device(ro)

       ! Since 'branching into or out of region is not allowed' by nvfortran,
       ! the following lines, which contain 'RETURN' can not be inside a
       ! acc region.
       if (ro == 0.0) then
          if(its == 0)then
             info=3
          else
             info = 0
          endif
          Tol = ro
          Iter = its
          RETURN
       end if

       ! set Tol1 for stopping criterion
       if (its == 0) then
          ro0 = ro
          if(DoTest) print *,'initial rnrm:',ro0
          if (TypeStop=='abs') then
             Tol1=Tol
             if (ro <= Tol1) then ! quit if accurate enough
                info = 3
                Tol  = ro
                Iter = its
                if(DoTest) print *,'GMRES: nothing to do. info = ',info
                RETURN
             end if
          else
             Tol1=Tol*ro
          end if
       end if

       !$acc parallel
       coeff = 1.0 / ro
       !$acc loop independent
       do i = 1, n
          Krylov_II(i,1)=coeff*Krylov_II(i,1)
       end do
       ! initialize 1-st term  of rhs of hessenberg system
       rs(1) = ro
       !$acc end parallel

       i = 0
       KRYLOVLOOP: do
          i = i + 1
          its = its + 1
          i1 = i + 1
          !
          !           Krylov_II(i1):=A*Krylov_II(i)
          !
          call matvec(Krylov_II(:,i),Krylov_II(:,i1),n)

          !-----------------------------------------
          !  modified gram - schmidt...
          !-----------------------------------------
          do j = 1, i

             t = dot_product_mpi(n, Krylov_II(:,j), Krylov_II(:,i1), iComm)
             !$acc update device(t)

             !$acc parallel
             hh(j,i) = t
             !$acc end parallel

             !$acc parallel loop gang vector independent
             do k = 1, n
                Krylov_II(k,i1) = Krylov_II(k,i1) - t*Krylov_II(k,j)
             end do
          end do

          t = sqrt(dot_product_mpi(n, Krylov_II(:,i1),Krylov_II(:,i1), iComm))
          !$acc update device(t)

          !$acc parallel
          hh(i1,i) = t
          if (t /= 0.0)then
             tmp = 1.0 / t
             !$acc loop gang vector independent
             do k = 1, n
                Krylov_II(k,i1) = tmp*Krylov_II(k,i1)
             end do
          endif
          !$acc end parallel

          !--------done with modified gram schmidt and arnoldi step

          !-------- now  update factorization of hh

          !-------- perform previous transformations  on i-th column of h
          !$acc serial
          do k = 2, i
             k1 = k-1
             t = hh(k1,i)
             hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
             hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
          end do

          gam = sqrt(hh(i,i)**2 + hh(i1,i)**2)
          if (gam == 0.0) gam = epsmac
          !-----------#  determine next plane rotation  #-------------------
          c(i) = hh(i,i)/gam
          s(i) = hh(i1,i)/gam
          rs(i1) = -s(i)*rs(i)
          rs(i) =  c(i)*rs(i)
          !---determine residual norm and test for convergence-
          hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
          ro = abs(rs(i1))
          !$acc end serial

          !$acc update host(ro)
          if(DoTest)then
             select case(TypeStop)
             case('rel')
                write(*,*) its,' matvecs, ',' ||rn||/||r0|| =',ro/ro0
             case('abs')
                write(*,*) its,' matvecs, ',' ||rn|| =',ro
             end select
          end if

          if (i >= nKrylov .or. (ro <= Tol1)) EXIT KRYLOVLOOP
       enddo KRYLOVLOOP

       !
       ! now compute solution. first solve upper triangular system.
       !
       ! rs := hh(1:i,1:i) ^-1 * rs

       !$acc serial
       do j = i, 1, -1
          if (rs(j)/=0.0) then
             rs(j)=rs(j)/hh(j,j)
             tmp = rs(j)
             do k=j-1,1,-1
                rs(k) = rs(k) - tmp*hh(k,j)
             enddo
          endif
       enddo
       !$acc end serial

       ! done with back substitution..
       ! now form linear combination to get solution
       !$acc parallel
       do j = 1, i
          !$acc loop gang vector independent
          do k = 1, n
             Sol(k) = Sol(k) + rs(j)*Krylov_II(k,j)
          end do
       end do
       !$acc end parallel

       ! exit from outer loop if converged or too many iterations
       if (ro <= Tol1 .or. its >= Iter) EXIT RESTARTLOOP
    end do RESTARTLOOP

    !    call cpu_time(finish)
    !    print '("TimeEndMainLoop = ",f6.3," seconds.")',finish-start

    Iter = its
    Tol = Tol/Tol1*ro ! (relative) tolerance achieved
    if(ro < Tol1)then
       info = 0
    elseif(ro < ro0)then
       info =  2
    else
       info = -2
    endif

    ! call timing_stop('gmres')

  end subroutine gmres
  !============================================================================
  subroutine bicgstab(matvec, rhs, qx, nonzero, n, tol, typestop, iter, info, &
       DoTest,iCommIn)

    ! Simple BiCGstab(\ell=1) iterative method
    ! Modified by G.Toth from the \ell<=2 version written
    ! by M.A.Botchev, Jan.'98.
    ! Parallelization for the BATS-R-US code (2000-2001) by G. Toth
    !
    ! This is the "vanilla" version of BiCGstab(\ell) as described
    ! in PhD thesis of D.R.Fokkema, Chapter 3.  It includes two enhancements
    ! to BiCGstab(\ell) proposed by G.Sleijpen and H.van der Vorst in
    ! 1) G.Sleijpen and H.van der Vorst "Maintaining convergence
    !    properties of BiCGstab methods in finite precision arithmetic",
    !    Numerical Algorithms, 10, 1995, pp.203-223
    ! 2) G.Sleijpen and H.van der Vorst "Reliable updated residuals in
    !    hybrid BiCG methods", Computing, 56, 1996, 141-163
    !
    ! {{ This code is based on:
    ! subroutine bistbl v1.0 1995
    !
    ! Copyright (c) 1995 by D.R. Fokkema.
    ! Permission to copy all or part of this work is granted,
    ! provided that the copies are not made or distributed
    ! for resale, and that the copyright notice and this
    ! notice are retained.  }}

    ! This subroutine determines the solution of A.QX=RHS, where
    ! the matrix-vector multiplication with A is performed by
    ! the subroutine 'matvec'. For symmetric matrix use the
    ! more efficient conjugate gradient (CG) algorithm!

    ! Arguments

    interface
       subroutine matvec(a,b,n)
         implicit none
         ! Calculate b = M.a where M is the matrix
         integer, intent(in) :: n
         real, intent(in) ::  a(n)
         real, intent(out) :: b(n)
       end subroutine matvec
    end interface
    !        subroutine for matrix vector multiplication

    integer, intent(in) :: n       ! number of unknowns
    real, intent(inout) :: rhs(n)  ! right-hand side / residual vector.
    real, intent(inout) :: qx(n)   ! initial guess /  solution vector.

    logical, intent(in):: nonzero  ! true  if qx contains initial guess
    real, intent(inout) :: tol     ! required / achieved norm of residual
    character (len=3), intent(in) :: typestop
    !  Determine stopping criterion (||.|| denotes the 2-norm):
    !  typestop='rel'    -- relative stopping crit.: ||res|| <= tol*||res0||
    !  typestop='abs'    -- absolute stopping crit.: ||res|| <= tol
    !  typestop='max'    -- maximum  stopping crit.: max(abs(res)) <= tol

    ! NOTE for typestop='rel' and 'abs':
    !   To save computational work, the value of
    !   residual norm used to check the convergence inside the main
    !   iterative loop is computed from
    !   projections, i.e. it can be smaller than the true residual norm
    !   (it may happen when e.g. the 'matrix-free' approach is used).
    !   Thus, it is possible that the true residual does NOT satisfy
    !   the stopping criterion ('rel' or 'abs').
    !   The true residual norm (or residual reduction) is reported on
    !   output in parameter TOL -- this can be changed to save 1 MATVEC
    !            (see comments at the end of the subroutine)

    integer, intent(inout) :: iter
    !       on input:  maximum number of iterations to be performed.
    !       on output: actual  number of iterations done.

    integer, intent(out)   :: info
    !       Gives reason for returning:
    !  abs(info)=  0 - solution found satisfying given tolerance.
    !              1 - iteration aborted due to division by very small value.
    !              2 - no convergence within maximum number of iterations.
    !              3 - initial guess satisfies the stopping criterion.
    !  sign(info)=  + - residual decreased
    !               - - residual did not reduce

    logical, intent(in)    :: DoTest        ! write debug info if true
    integer, intent(in), optional:: iCommIn ! MPI communicator

    ! Local parameters

    integer, parameter :: qz_=1,zz_=3,y0_=5,yl_=6,qy_=7

    ! Local variables (only 4 big vectors are needed):
    integer :: iComm             ! actual MPI communicator
    ! used to be automatic arrays
    real, allocatable :: bicg_r(:), bicg_u(:), bicg_r1(:), bicg_u1(:)

    ! allocatable array for initial guess
    real, allocatable :: qx0(:)

    real :: rwork(2,7)

    !--------------------------------------------------------------------------
    logical GoOn, rcmp, xpdt
    integer nmv
    real :: alpha, beta, omega, rho0, rho1, sigma
    real :: varrho, hatgamma
    real :: assumedzero, rnrm0, rnrm, rnrmMax0, rnrmMax
    real :: mxnrmx, mxnrmr, kappa0, kappal
    !--------------------------------------------------------------------------
    iComm = MPI_COMM_SELF
    if(present(iCommIn)) iComm = iCommIn

    ! Allocate arrays that used to be automatic
    allocate(bicg_r(n), bicg_u(n), bicg_r1(n), bicg_u1(n));

    if(DoTest)write(*,*)'BiCGSTAB tol,iter:',tol,iter

    info = 0

    if (tol<=0.0) call CON_stop('Error in BiCGSTAB: tolerance < 0')
    if (iter<=1)  call CON_stop('Error in BiCGSTAB: maxmatvec < 2')

    if(range(1.0)>100)then
       assumedzero=0.0000000000000001
    else
       assumedzero=0.00000001
    end if

    !
    !     --- Initialize first residual
    !
    ! Calculate initial residual
    if(nonzero)then
       ! Store initial guess into qx0
       allocate(qx0(n))
       qx0=qx
       call matvec(qx,bicg_r,n)
       bicg_r = rhs - bicg_r
    else
       bicg_r = rhs
    end if

    qx = 0.0
    bicg_u=0.0

    nmv = 0
    !
    !     --- Initialize iteration loop
    !

    rnrm0 = sqrt( dot_product_mpi(n, bicg_r,bicg_r,iComm))

    rnrm = rnrm0
    if(DoTest) print *,'initial rnrm:',rnrm

    mxnrmx = rnrm0
    mxnrmr = rnrm0
    rcmp = .false.
    xpdt = .false.

    alpha = 0.0
    omega = 1.0
    sigma = 1.0
    rho0 =  1.0
    !
    !     --- Iterate
    !
    select case(typestop)
    case('rel')
       GoOn = rnrm>tol*rnrm0 .and. nmv<iter
       assumedzero = assumedzero*rnrm0
       rnrmMax = 0
       rnrmMax0 = 0

    case('abs')
       GoOn = rnrm>tol       .and. nmv<iter
       assumedzero = assumedzero*rnrm0
       rnrmMax = 0
       rnrmMax0 = 0

    case('max')
       rnrmMax0 = maxval_abs_mpi(bicg_r,iComm)
       rnrmMax  = rnrmMax0
       if(DoTest) print *,'initial rnrmMax:',rnrmMax
       GoOn = rnrmMax>tol    .and. nmv<iter
       assumedzero = assumedzero*rnrmMax
    case default
       call CON_stop('Error in BiCGSTAB: unknown typestop value')
    end select

    if (.not.GoOn) then
       iter = nmv
       info = 3
       if(DoTest) print *,'BiCGSTAB: nothing to do. info = ',info
       deallocate(bicg_r, bicg_u, bicg_r1, bicg_u1)

       if(nonzero)then
          qx = qx0
          deallocate(qx0)
       end if

       RETURN
    end if

    do while (GoOn)
       !
       !     =====================
       !     --- The BiCG part ---
       !     =====================
       !
       rho0 = -omega*rho0

       rho1 = dot_product_mpi(n, rhs,bicg_r,iComm)

       if (abs(rho0)<assumedzero**2) then
          info = 1
          deallocate(bicg_r, bicg_u, bicg_r1, bicg_u1)
          if(nonzero)then
             qx = qx + qx0
             deallocate(qx0)
          end if
          RETURN
       endif
       beta = alpha*(rho1/rho0)
       rho0 = rho1
       bicg_u = bicg_r - beta*bicg_u

       ! DEBUG
       ! if(DoTest)&
       !     write(*,*)'rho0, rho1, alpha, beta:',rho0, rho1, alpha, beta

       call matvec(bicg_u,bicg_u1,n)
       nmv = nmv+1

       sigma=dot_product_mpi(n, rhs,bicg_u1,iComm)

       if (abs(sigma)<assumedzero**2) then
          info = 1
          deallocate(bicg_r, bicg_u, bicg_r1, bicg_u1)
          if(nonzero)then
             qx = qx + qx0
             deallocate(qx0)
          end if
          RETURN
       endif

       alpha = rho1/sigma

       qx     = qx     + alpha*bicg_u
       bicg_r = bicg_r - alpha*bicg_u1

       call matvec(bicg_r,bicg_r1,n)
       nmv = nmv+1

       rnrm = sqrt( dot_product_mpi(n, bicg_r,bicg_r,iComm) )

       mxnrmx = max (mxnrmx, rnrm)
       mxnrmr = max (mxnrmr, rnrm)

       ! DEBUG
       ! if(DoTest)&
       !     write(*,*)'rho0, rho1, beta, sigma, alpha, rnrm:',&
       !     rho0, rho1, beta, sigma, alpha, rnrm

       !
       !  ==================================
       !  --- The convex polynomial part ---
       !  ==================================
       !
       !    --- Z = R'R a 2 by 2 matrix
       ! i=1,j=0
       rwork(1,1) = dot_product_mpi(n, bicg_r,bicg_r,iComm)

       ! i=1,j=1
       rwork(2,1) = dot_product_mpi(n, bicg_r1,bicg_r,iComm)
       rwork(1,2) = rwork(2,1)

       ! i=2,j=1
       rwork(2,2) = dot_product_mpi(n, bicg_r1,bicg_r1,iComm)

       !
       !   --- tilde r0 and tilde rl (small vectors)
       !
       rwork(1:2,zz_:zz_+1)   = rwork(1:2,qz_:qz_+1)
       rwork(1,y0_) = -1.0
       rwork(2,y0_) = 0.0

       rwork(1,yl_) = 0.0
       rwork(2,yl_) = -1.0
       !
       !   --- Convex combination
       !
       rwork(1:2,qy_) = rwork(1,yl_)*rwork(1:2,qz_) + &
            rwork(2,yl_)*rwork(1:2,qz_+1)

       kappal = sqrt( sum( rwork(1:2,yl_)*rwork(1:2,qy_) ) )

       rwork(1:2,qy_) = rwork(1,y0_)*rwork(1:2,qz_) + &
            rwork(2,y0_)*rwork(1:2,qz_+1)

       kappa0 = sqrt( sum( rwork(1:2,y0_)*rwork(1:2,qy_) ) )

       varrho = sum( rwork(1:2,yl_)*rwork(1:2,qy_) )
       varrho = varrho / (kappa0*kappal)

       hatgamma = sign(1.0,varrho)*max(abs(varrho),0.7) * (kappa0/kappal)

       rwork(1:2,y0_) = -hatgamma*rwork(1:2,yl_) + rwork(1:2,y0_)

       !
       !    --- Update
       !
       omega = rwork(2,y0_)

       bicg_u = bicg_u - omega*bicg_u1

       qx     = qx     + omega*bicg_r

       bicg_r = bicg_r - omega*bicg_r1

       rwork(1:2,qy_) = rwork(1,y0_)*rwork(1:2,qz_) + &
            rwork(2,y0_)*rwork(1:2,qz_+1)

       rnrm = sqrt( max(0.0, sum( rwork(1:2,y0_)*rwork(1:2,qy_) )) )

       select case(typestop)
       case('rel')
          GoOn = rnrm>tol*rnrm0 .and. nmv<iter
          if(DoTest) print *, nmv,' matvecs, ', ' ||rn||/||r0|| =',rnrm/rnrm0
       case('abs')
          GoOn = rnrm>tol       .and. nmv<iter
          if(DoTest) print *, nmv,' matvecs, ||rn|| =',rnrm
       case('max')
          rnrmMax = maxval_abs_mpi(bicg_r,iComm)
          GoOn = rnrmMax>tol    .and. nmv<iter
          if(DoTest) print *, nmv,' matvecs, max(rn) =',rnrmMax
       end select

    end do
    !
    !     =========================
    !     --- End of iterations ---
    !     =========================

    ! Add initial guess if it is non-zero
    if(nonzero)then
       qx = qx + qx0
       deallocate(qx0)
    end if

    select case(typestop)
    case('rel')
       if (rnrm>tol*rnrm0) info = 2
       tol = rnrm/rnrm0

    case('abs')
       if (rnrm>tol) info = 2
       tol = rnrm

    case('max')
       if (rnrmMax>tol) info = 2
       tol = rnrmMax
    end select

    if((typestop/='max'.and.rnrm>rnrm0).or.(typestop=='max'.and.rnrmMax&
         >rnrmMax0)) info=-info

    iter = nmv
    deallocate(bicg_r, bicg_u, bicg_r1, bicg_u1)

  end subroutine bicgstab
  !============================================================================
  subroutine cg(matvec, Rhs_I, Sol_I, IsInit, n, Tol, TypeStop, &
       nIter, iError, DoTest, iCommIn, JacobiPrec_I, preconditioner)

    ! Based on an implementation by I. Sokolov

    ! Conjugated gradients, works if and only if the matrix A
    ! is symmetric and positive definite

    ! subroutine for matrix vector multiplication
    interface
       subroutine matvec(Vec_I, MatVec_I, n)
         ! Calculate MatVec_I = Mat_II.Vec_I where Mat_II is a symmetric
         ! positive definite matrix
         implicit none
         integer, intent(in) :: n
         real,    intent(in) :: Vec_I(n)
         real,    intent(out):: MatVec_I(n)
       end subroutine matvec
    end interface

    integer, intent(in) :: n         !  number of unknowns.
    real, intent(inout) :: Rhs_I(n)    ! right hand side vector
    real, intent(inout) :: Sol_I(n)  ! initial guess / solution vector
    logical, intent(in) :: IsInit    ! true  if Sol contains initial guess
    real,    intent(inout) :: Tol       ! required / achieved residual

    optional :: preconditioner
    interface
       subroutine preconditioner(Vec_I, PrecVec_I, n)
         ! Calculate PrecVec_I = Prec_II.Vec_I where Prec_II is the
         ! preconditione matrix that should be symmetric and positive definite
         implicit none
         integer, intent(in) :: n
         real, intent(in)    :: Vec_I(n)
         real, intent(out)   :: PrecVec_I(n)
       end subroutine preconditioner
    end interface

    character (len=3), intent(in) :: TypeStop
    !      Determine stopping criterion (||.|| denotes the 2-norm):
    !      typestop='rel'    -- relative stopping crit.:||res|| <= Tol*||res0||
    !      typestop='abs'    -- absolute stopping crit.: ||res|| <= Tol

    integer, intent(out):: iError    ! info about convergence
    logical, intent(in) :: DoTest
    integer, intent(inout) :: nIter  ! maximum/actual number of iterations

    integer, intent(in), optional :: iCommIn   ! MPI communicator
    real,    intent(in), optional :: JacobiPrec_I(n) ! Jacobi preconditioner

    ! Local variables
    integer :: iComm                           ! MPI communicator
    integer :: MaxIter
    real :: rDotR, pDotADotP , rDotR0, rDotRMax, Alpha, Beta

    ! These arrays used to be automatic
    real, allocatable :: Vec_I(:), aDotVec_I(:), PrecRhs_I(:)
    !--------------------------------------------------------------------------
    iComm = MPI_COMM_SELF ! Assign the MPI communicator
    if(present(iCommIn)) iComm = iCommIn

    ! Allocate the vectors needed for CG
    allocate(Vec_I(n), aDotVec_I(n))

    if(present(JacobiPrec_I) .or. present(preconditioner)) &
         allocate(PrecRhs_I(n))

    MaxIter = nIter

    ! compute initial residual vector

    if(IsInit)then
       call matvec(Sol_I, ADotVec_I, n)
       Rhs_I = Rhs_I - ADotVec_I
    else
       Sol_I = 0.0
    end if

    Vec_I = 0.0

    if(present(JacobiPrec_I))then
       PrecRhs_I = JacobiPrec_I*Rhs_I
       rDotR0 = dot_product_mpi(n, Rhs_I, PrecRhs_I, iComm)
    elseif(present(preconditioner))then
       call preconditioner(Rhs_I, PrecRhs_I, n)
       rDotR0 = dot_product_mpi(n, Rhs_I, PrecRhs_I, iComm)
    else
       rDotR0 = dot_product_mpi(n, Rhs_I, Rhs_I, iComm)
    end if

    if(TypeStop=='abs')then
       rDotRMax = Tol**2
    else
       rDotRMax = Tol**2 * max(rDotR0, 1e-99)
    end if

    if(DoTest)write(*,*)'CG rDotR0, rDotRMax=', rDotR0, rDotRMax

    rDotR = rDotR0

    ! Solve the problem iteratively
    do nIter = 1, MaxIter

       if( rDotR <= rDotRMax ) EXIT

       Alpha = 1.0/rDotR

       if(present(JacobiPrec_I) .or. present(preconditioner))then
          Vec_I = Vec_I + Alpha * PrecRhs_I
       else
          Vec_I = Vec_I + Alpha * Rhs_I
       end if

       UsePDotADotP = .false.

       call matvec(Vec_I, aDotVec_I, n)

       if(UsePDotADotP)then
          call MPI_ALLREDUCE(pDotADotPPe, pDotADotP, 1, MPI_REAL, MPI_SUM, &
               iComm, iError)
          UsePDotADotP = .false.
       else
          pDotADotP = dot_product_mpi(n, Vec_I, aDotVec_I, iComm)
       end if
       Beta = 1.0/pDotADotP

       Rhs_I = Rhs_I - Beta * aDotVec_I
       Sol_I = Sol_I + Beta * Vec_I

       if(present(JacobiPrec_I))then
          PrecRhs_I = JacobiPrec_I*Rhs_I
          rDotR = dot_product_mpi(n, Rhs_I, PrecRhs_I, iComm)
       elseif(present(preconditioner))then
          call preconditioner(Rhs_I, PrecRhs_I, n)
          rDotR = dot_product_mpi(n, Rhs_I, PrecRhs_I, iComm)
       else
          rDotR = dot_product_mpi(n, Rhs_I, Rhs_I, iComm)
       end if
       if(DoTest)write(*,*)'CG nIter, rDotR=',nIter, rDotR

    end do

    ! Deallocate the temporary vectors
    deallocate(Vec_I, aDotVec_I)
    if(allocated(PrecRhs_I))deallocate(PrecRhs_I)

    ! Calculate the achieved tolerance
    if(TypeStop=='abs')then
       Tol = sqrt(rDotR)
    else
       Tol = sqrt(rDotR/max(rDotR0, 1e-99))
    end if

    ! Set the error flag
    if(rDotR0 < rDotRMax)then
       nIter  = 0
       iError = 3
    elseif(nIter <= MaxIter)then
       iError = 0
    elseif(rDotR < rDotR0)then
       iError = 2
    else
       iError = -2
    end if

  end subroutine cg
  !============================================================================
  real function accurate_dot_product(a_I, b_I, LimitIn, iComm)

    ! Algorithm and implementation by G. Toth 2009

    ! arrays to be multiplied and added up
    real,    intent(in)           :: a_I(:), b_I(:)
    real,    intent(in), optional :: LimitIn  ! limit for splitting numbers
    integer, intent(in), optional :: iComm    ! MPI communicator

    ! This function returns the dot product of the a_I and b_I arrays
    ! for 1 or more processors.
    ! This is done more accurately than by simple sum(a_I*b_I) by splitting
    ! the values into two shorter parts in the binary representation.
    ! Eg. a=0.375353534 may be split into an upper part 0.375 (=3/8)
    ! and a lower part a - 3/8 = 0.000353534.
    ! The lower and uppper parts are added up independenty before they are
    ! added together, thus the round-off errors are somewhat mitigated.
    ! If LimitIn is provided, it should be an integer power of 2,
    ! and it provides the splitting limit. In the example the limit was 1/8.
    ! If LimitIn is missing the limit is based on maxval(abs(a_I*b_I)).
    ! If the optional MPI communicator argument iComm is provided,
    ! then the dot product (and max) are done for all processors belonging to
    ! iComm.

    ! Number of parts the binary numbers are cut into
    integer, parameter :: nPart = 2

    ! Limit is set to this fraction of the largest argument
    real, parameter :: Ratio = 1e-8

    integer :: i, n, iError, nProc
    real :: Limit, Maximum, a
    real :: Part_I(nPart), SumPart_I(nPart), SumAllPart_I(nPart)
    !--------------------------------------------------------------------------
    n = size(a_I)

    if(present(LimitIn))then
       Limit = LimitIn
    else
       Maximum = maxval(abs(a_I*b_I))
       if(present(iComm))then
          if(iComm /= MPI_COMM_SELF)then
             call MPI_COMM_SIZE(iComm, nProc, iError)
             if(nProc>1)call MPI_allreduce( &
                  MPI_IN_PLACE, Maximum, 1, MPI_REAL, MPI_MAX, iComm, iError)
          end if
       endif
       ! Formulate an integer power of 2 near Ratio*Maximum
       Limit = 2.0**exponent(Ratio*Maximum)
    end if

    SumPart_I = 0.0
    do i = 1, n
       a = a_I(i)*b_I(i)
       Part_I(1) = modulo(a, Limit)
       Part_I(2) = a - Part_I(1)
       SumPart_I = SumPart_I + Part_I
    end do

    if(present(iComm))then
       if(iComm /= MPI_COMM_SELF)then
          call MPI_allreduce(SumPart_I, SumAllPart_I, nPart, MPI_REAL, &
               MPI_SUM, iComm, iError)
          SumPart_I = SumAllPart_I
       end if
    endif

    accurate_dot_product = sum(SumPart_I)

  end function accurate_dot_product
  !============================================================================
  real function dot_product_mpi(n, a_I, b_I, iComm)

    integer, intent(in) :: n
    real, intent(in)    :: a_I(n), b_I(n)
    integer, intent(in) :: iComm

    integer :: i

    real :: DotProduct, DotProductMpi
    integer :: iError
    !--------------------------------------------------------------------------
#ifndef _OPENACC
    if(UseAccurateSum)then
       dot_product_mpi = accurate_dot_product(a_I, b_I, iComm=iComm)
       RETURN
    end if

    DotProduct = dot_product(a_I, b_I)
#else
    DotProduct = 0.0
    !$acc parallel loop gang vector reduction(+:DotProduct)
    do i = 1, n
       DotProduct = DotProduct + a_I(i)*b_I(i)
    end do
#endif

    if(iComm == MPI_COMM_SELF) then
       dot_product_mpi = DotProduct
       RETURN
    end if
    call MPI_allreduce(DotProduct, DotProductMpi, 1, MPI_REAL, MPI_SUM, &
         iComm, iError)
    dot_product_mpi = DotProductMpi
  end function dot_product_mpi
  !============================================================================
  real function maxval_abs_mpi(a_I, iComm)

    real, intent(in)    :: a_I(:)
    integer, intent(in) :: iComm

    real :: MaxvalAbs, MaxvalAbsMpi
    integer :: iError
    !--------------------------------------------------------------------------
    MaxvalAbs = maxval(abs(a_I))
    if(iComm == MPI_COMM_SELF)then
       maxval_abs_mpi = MaxvalAbs
       RETURN
    end if

    call MPI_allreduce(MaxvalAbs, MaxvalAbsMpi, 1, MPI_REAL, MPI_MAX, &
         iComm, iError)
    maxval_abs_mpi = MaxvalAbsMpi

  end function maxval_abs_mpi
  !============================================================================
  subroutine prehepta(nBlock, n, m1, m2, PrecondParam, d, e, f, e1, f1, e2, f2)
    !$acc routine vector

    ! This routine constructs an incomplete block LU-decomposition
    ! of a hepta- or penta-diagonal matrix in such a way that L+U has the
    ! same blockstructure as A. For block tri-diagonal matrix the subroutine
    ! provides a full LU decompostion.
    !
    ! For penta-diagonal matrix, set M2=nBlock and e2,f2 can be omitted.
    ! For tri-diagonal matrix set M1=M2=nBlock and e1,f1,e2,f2 can be omitted.
    !
    !     Gustafsson modification
    !
    ! It is possible to encorporate the so-called Gustafsson
    ! modification for the blocks on the main diagonal.
    ! In this appoach, a splitting A=LU+R, is constructed in such a way
    ! that the block row sums of R are zero. For systems of
    ! linear equations coming from problems with an elliptic character
    ! this can strongly improve the convergence behaviour of
    ! iterative methods. See page 22 of Phd thesis 'Preconditioning
    ! for sparse ..... ' for an illustration of this phenomenon.

    integer, intent(in):: n, M1, M2, nBlock
    real, intent(in)   :: PrecondParam
    real, intent(inout):: d(n,n,nBlock)
    real, intent(inout), dimension(n,n,nBlock), optional:: e, f, e1, f1, e2, f2

    !     Description of arguments:
    !
    ! nBlock:  Number of diagonal blocks.
    ! N:       The size of the blocks.
    ! M1:      Distance of outer blocks to the main diagonal blocks.
    ! M2:      Distance of outer-most blocks to main diagonal blocks.
    !           1 < M1 < M2.
    !           The matrix has a blockstructure corresponding to
    !           a seven-point stencil on a three-dimensional,
    !           rectangular grid. The blocks corresonding to the
    !           direction in which grid points are numbered
    !           first are the sub- and superdiagonal blocks.
    !           The blocks corresponding to the direction in which grid
    !           points are numbered secondly have distance M1 from the main
    !           diagonal blocks. Finally, the blocks corresponding to the
    !           direction in which grid points are numbered
    !           last have distance M2 from the main diagonal blocks.
    !
    ! PrecondParam:      The parameter for Gustafsson modification:
    !           +5 Jacobi       prec: multiply the diagonal elements
    !           +4 Block-Jacobi prec: invert the original diagonal blocks
    !           +3 Gauss-Seidel prec: invert the original diagonal blocks
    !                               and premultiply upper diagonal blocks
    !           +2 DILU  prec: LU for diagonal, keep off-diagonal blocks
    !           +1 BILU  prec: LU for diagonal, premultiply U with D^-1
    !           <0 MBILU prec: Gustaffson modification of diagonal blocks
    !                          using -1 <= PrecondParam < 0 parameter

    !
    ! d, e, f, e1, f1, e2, f2:
    !          on entrance: matrix A
    !          on exit: L + U - I (the diagonal of U is I)
    !           The matrix A and L+U are block heptadiagonal.
    !           The blocks are stored as follows:
    !           d(j): j=1..nBlock        main diagonal
    !           e(j): j=2..nBlock        sub diagonal blocks.
    !           f(j): j=1..nBlock-1      super diagonal blocks.
    !           e1(j): j=M1+1..nBlock    blocks in the lower-triangular
    !                                     part with distance M1 from
    !                                     the main diagonal.
    !           f1(j): j=1..nBlock-M1    blocks in the upper-triangular
    !                                     part with distance M1 from
    !                                     the main diagonal.
    !           e2(j): j=M2+1..nBlock    blocks in the lower-triangular
    !                                     part with distance M2 from
    !                                     the main diagonal.
    !           f2(j): j=1..nBlock-M2    blocks in the upper-triangular
    !                                     part with distance M2 from
    !                                     the main diagonal.
    !          It is assumed that the
    !          blocks are not very sparse, so the sparsity pattern
    !          within the separate blocks is not exploited.
    !          For example, the (i,k)-element of the j-th block on the
    !          main diagonal is stored in d(i,k,j).
    !
    !     Local variables:
    !
    ! dd:      a single block for manipulating diagonal block d(*,*,j)
    !
    ! pivot:   integer array which contains the sequence generated
    !          by partial pivoting in subroutine 'Lapack_getrf'.
    !

#ifndef _OPENACC
    ! these used to be automatic arrays
    real,    allocatable :: dd(:,:)
    integer, allocatable :: pivot(:)

    real, allocatable:: fOrig_VVI(:,:,:), f1Orig_VVI(:,:,:), f2Orig_VVI(:,:,:)

    ! info variable for lapack routines
    integer :: i, j, info, iPrecond
#endif
    ! call timing_start('precond')
    !--------------------------------------------------------------------------
    if(n == 1 .and. nint(PrecondParam) /= Dilu_) then
       if(nint(PrecondParam) == Bilu1_)then
          call prehepta_scalar_fast(nBlock, m1, d, e, f)
       else
          call prehepta_scalar(nBlock, m1, m2, PrecondParam, &
               d, e, f, e1, f1, e2, f2)
       end if
       RETURN
    else
#ifdef _OPENACC
       call CON_stop_simple('Error: unsupported on GPU!')
#endif
    endif

#ifndef _OPENACC
    ! Allocate arrays that used to be automatic
    allocate(dd(N,N), pivot(N))

    iPrecond = nint(PrecondParam)

    if(iPrecond == Dilu_)then
       allocate(fOrig_VVI(n,n,nBlock))
       fOrig_VVI = f
       if(present(f1))then
          allocate(f1Orig_VVI(n,n,nBlock))
          f1Orig_VVI = f1
          if(present(f2))then
             allocate(f2Orig_VVI(n,n,nBlock))
             f2Orig_VVI = f2
          end if
       end if
    end if

    do j = 1, nBlock

       dd = d(:,:,j)

       ! Modify D according to LU decomposition except for Jacobi/Gauss-Seidel
       if (iPrecond < GaussSeidel_)then
          ! D = D - E.F(j-1) - E1.F1(j-M1) - E2.F2(j-M2)
          if (j > 1 )call BLAS_gemm('n', 'n', n, n, n, -1.0, &
               e(:,:,j),  N, f(:,:,j- 1), n, 1.0, dd, n)
          if (j > m1) call BLAS_gemm('n', 'n', n, n, n, -1.0, &
               e1(:,:,j), N, f1(:,:,j-M1), n, 1.0, dd, n)
          if (j > m2) call BLAS_gemm('n','n', n, n, n, -1.0, &
               e2(:,:,j), N, f2(:,:,j-M2), n, 1.0, dd, n)
       end if
       if(iPrecond <= Mbilu_)then
          ! Relaxed Gustafsson modification for MBILU

          ! D = D + PrecondParam*
          !    (  E2.F(j-M2) + E2.F1(j-M2)
          !     + E1.F(j-M1) + E1.F2(j-M1)
          !     + E.F1(j-1)  + E .F2(j-1) )

          if (j > M2) then
             call BLAS_GEMM('n', 'n', N, N, N, PrecondParam, &
                  e2(:,:,j), N,  f(:,:,j-M2), N, 1.0, dd, N)
             call BLAS_GEMM('n', 'n', N, N, N, PrecondParam, &
                  e2(:,:,j), N, f1(:,:,j-M2), N, 1.0, dd, N)
          end if
          if (j > M1) call BLAS_GEMM('n','n',N,N,N,PrecondParam, &
               e1(:,:,j), N, f(:,:,j-M1), N, 1.0, dd, N)
          if (j > M1 .and. j-M1 <= nBlock-M2) &
               call BLAS_GEMM('n','n',N,N,N,PrecondParam, &
               e1(:,:,j), N, f2(:,:,j-M1), N, 1.0, dd, N)
          if (j > 1 .and. j-1 <= nBlock-M2) &
               call BLAS_GEMM('n', 'n', N, N, N, PrecondParam, &
               e(:,:,j), N, f2(:,:,j-1 ), N, 1.0, dd, N)
          if (j>1 .and. j-1 <= nBlock-M1) &
               call BLAS_GEMM('n', 'n', N, N, N, PrecondParam, &
               e(:,:,j), N, f1(:,:,j-1 ), N, 1.0, dd, N)
       end if

       ! Invert the diagonal block by first factorizing then solving D.D'=I
       call LAPACK_getrf(n, n, dd, n, pivot, info)
       ! Set the right hand side as identity matrix
       d(:,:,j) = 0.0
       do i = 1, N
          d(i,i,j) = 1.0
       end do
       ! Solve the problem, returns D^-1 into d(j)
       call LAPACK_getrs('n', n, n, dd, n, pivot, d(:,:,j), n, info)

       ! For Jacobi prec no need to do anything with upper diagonal blocks
       if (nint(PrecondParam) == BlockJacobi_) CYCLE

       ! Pre-multiply U with D^-1 to make the decomposition
       ! as well as multiplication with U^{-1} more efficient
       ! F2 = D^{-1}.F2, F1 = D^{-1}.F1, F = D^{-1}.F
       if (j   < nBlock)then
          dd = f(:,:,j)
          call BLAS_gemm('n', 'n', n, n, n, 1.0, d(:,:,j), n, dd, n, 0.0, &
               f(:,:,j), n)
       end if
       if (j+M1 <= nBlock)then
          dd = f1(:,:,j)
          call BLAS_gemm('n', 'n', n, n, n, 1.0, d(:,:,j), n, dd, n, 0.0, &
               f1(:,:,j), n)
       end if
       if (j+M2 <= nBlock)then
          dd = f2(:,:,j)
          call BLAS_gemm('n', 'n', n, n, n, 1.0, d(:,:,j), n, dd, n, 0.0, &
               f2(:,:,j), n)
       end if

    end do

    if(iPrecond == Dilu_)then
       ! DILU prec requires original diagonal blocks in U, so restore it
       f = fOrig_VVI
       deallocate(fOrig_VVI)
       if(present(f1))then
          f1 = f1Orig_VVI
          deallocate(f1Orig_VVI)
          if(present(f2))then
             f2 = f2Orig_VVI
             deallocate(f2Orig_VVI)
          end if
       end if
    end if

    ! Deallocate arrays
    deallocate(dd, pivot)

    ! call timing_stop('precond')
#endif
  end subroutine prehepta
  !============================================================================
  subroutine prehepta_scalar(nBlock, m1, m2, PrecondParam, &
       d, e, f, e1, f1, e2, f2)
    !$acc routine vector

    integer, intent(in):: M1, M2, nBlock
    real, intent(in)   :: PrecondParam
    real, intent(inout):: d(nBlock)
    real, intent(inout), dimension(nBlock), optional:: e, f, e1, f1, e2, f2

    real :: dd

    ! info variable for lapack routines
    integer :: j, iPrecond
    !--------------------------------------------------------------------------
    iPrecond = nint(PrecondParam)

    !$acc loop seq
    do j = 1, nBlock

       dd = d(j)
       if (iPrecond < GaussSeidel_)then
          if (j > 1 ) dd = dd - e(j)*f(j-1)
          if (j > m1) dd = dd - e1(j)*f1(j-M1)
          if (j > m2) dd = dd - e2(j)*f2(j-M2)
       end if
       if(iPrecond <= Mbilu_)then
          if (j > M2) then
             dd = dd + PrecondParam*(e2(j)*f(j-M2))
             dd = dd + PrecondParam*(e2(j)*f1(j-M2))
          end if

          if (j > M1) dd = dd + PrecondParam*(e1(j)*f(j-M1))

          if (j > M1 .and. j-M1 <= nBlock-M2) &
               dd = dd + PrecondParam*(e1(j)*f2(j-M1))

          if (j > 1 .and. j-1 <= nBlock-M2) &
               dd = dd + PrecondParam*(e(j)*f2(j-1))

          if (j>1 .and. j-1 <= nBlock-M1) &
               dd = dd + PrecondParam*(e(j)*f1(j-1))
       end if

       d(j) = 1./dd

       ! For Jacobi prec no need to do anything with upper diagonal blocks
       if (nint(PrecondParam) == BlockJacobi_) CYCLE

       if (j   < nBlock)then
          f(j) = d(j)*f(j)
       end if
       if (j+M1 <= nBlock)then
          f1(j) = d(j)*f1(j)
       end if
       if (j+M2 <= nBlock)then
          f2(j) = d(j)*f2(j)
       end if

    end do

  end subroutine prehepta_scalar
  !============================================================================
  subroutine prehepta_scalar_fast(nBlock, m1, d, e, f)
    !$acc routine vector

    ! Purely 1D preconditioning. Can be parallelized with nBlock/m1 threads,
    ! because every M1-s element of e and f are 0.

    integer, intent(in):: nBlock, m1
    real, intent(inout), dimension(nBlock):: d, e, f

    integer:: i, j
    !--------------------------------------------------------------------------
    !$acc loop vector
    do i = 1, nBlock, m1
       d(i) = 1/d(i)
       f(i) = d(i)*f(i)
       !$acc loop seq
       do j = i+1, i+m1-1
          d(j) = 1/(d(j) - e(j)*f(j-1))
          f(j) = d(j)*f(j)
       end do
    end do

  end subroutine prehepta_scalar_fast
  !============================================================================
  subroutine Uhepta(IsInverse, nBlock, n, M1, M2, x, f, f1, f2)
    !$acc routine vector

    ! G. Toth, 2001

    ! This routine multiplies x with the upper triagonal U or U^{-1}
    ! which must have been constructed in subroutine prehepta.
    !
    ! For penta-diagonal matrix, set M2=nBlock and f2 can be omitted.
    ! For tri-diagonal matrix set M1=M2=nBlock and f1,f2 can be omitted.

    logical, intent(in) :: IsInverse
    integer, intent(in) :: N, M1, M2, nBlock
    real, intent(inout) :: x(N,nBlock)
    real, intent(in)    :: f(N,N,nBlock)
    real, intent(in), optional :: f1(n,n,nBlock), f2(n,n,nBlock)

    !     Description of arguments:
    !
    ! IsInverse: logical switch
    !          Multiply by U^{-1} if true, otherwise multiply by U
    !
    ! nBlock:  Number of diagonal blocks.
    ! N:       the size of the blocks.
    ! M1:      distance of blocks to the main diagonal blocks.
    !          set M1=nBlock for block tri-diagonal matrices!
    ! M2:      distance of outer-most blocks to main diagonal blocks.
    !          set M2=nBlock for block tri- and penta-diagonal matrices.
    !
    ! x:       On input, the vector to be multiplied with U or U^{-1}.
    !          On output, the result of U.x or U^{-1}.x
    !
    ! f, f1, f2
    !           The matrix U is assumed to be in three block diagonals
    !
    !           The blocks are stored as follows:
    !           f(j): j=1..nBlock-1   super diagonal blocks.
    !
    !           f1(j): j=1..nBlock-M1 blocks in the upper-triangular part with
    !                                 distance M1 from the main diagonal.
    !                                 Omit for block tri-diagonal matrix!
    !
    !           f2(j): j=1..nBlock-M2 blocks in the upper-triangular part with
    !                                 distance M2 from the main diagonal.
    !                                 Omit for block tri/penta-diagonal matrix!
    !
    ! It is assumed that the blocks are not very sparse,
    ! so the sparsity pattern  within the separate blocks is not exploited.
    ! For example, the (i,k)-element
    ! of the j-th block on the super diagonal is stored in f(i,k,j).

    integer :: j
    !--------------------------------------------------------------------------
    ! call timing_start('Uhepta')

    if(n == 1)then
       if(.not.present(f1))then
          call upper_hepta_scalar_fast(IsInverse, nBlock, M1, x, f)
       else
          call upper_hepta_scalar(IsInverse, nBlock, M1, M2, x, f, f1, f2)
       end if
       RETURN
    end if

#ifndef _OPENACC
    if(n <= 20)then
       ! F90 VERSION
       if(IsInverse)then
          !  x' := U^{-1}.x = x - F.x'(j+1) - F1.x'(j+M1) - F2.x'(j+M2)
          do j=nBlock-1,1,-1
             !  x' := U^{-1}.x = x - F.x'(j+1) - F1.x'(j+M1) - F2.x'(j+M2)
             if (j+M2<=nBlock) then
                x(:,j) = x(:,j) - matmul( f(:,:,j),x(:,j+1 )) &
                     - matmul(f1(:,:,j),x(:,j+M1)) &
                     - matmul(f2(:,:,j),x(:,j+M2))
             else if(j+M1<=nBlock) then
                x(:,j) = x(:,j) - matmul( f(:,:,j),x(:,j+1 )) &
                     - matmul(f1(:,:,j),x(:,j+M1))
             else
                x(:,j) = x(:,j) - matmul(f(:,:,j),x(:,j+1))
             end if
          end do
       else
          !  x := U.x = x + F.x(j+1) + F1.x(j+M1) + F2.x(j+M2)
          do j = 1, nBlock-1
             if (j+M2<=nBlock) then
                x(:,j) = x(:,j) + matmul( f(:,:,j),x(:,j+1 )) &
                     + matmul(f1(:,:,j),x(:,j+M1)) &
                     + matmul(f2(:,:,j),x(:,j+M2))
             else if (j+M1<=nBlock) then
                x(:,j) = x(:,j) + matmul( f(:,:,j),x(:,j+1 )) &
                     + matmul(f1(:,:,j),x(:,j+M1))
             else
                x(:,j) = x(:,j) + matmul(f(:,:,j),x(:,j+1))
             end if
          end do
       end if
    else
       ! BLAS VERSION
       if(IsInverse)then
          !  x' := U^{-1}.x = x - F.x'(j+1) - F1.x'(j+M1) - F2.x'(j+M2)
          do j=nBlock-1,1,-1
             call BLAS_gemv('n', n, n, -1.0, &
                  f(:,:,j), n, x(:,j+1 ), 1, 1.0, x(:,j), 1)
             if(j+M1<=nBlock) call BLAS_gemv('n', n, n, -1.0, &
                  f1(:,:,j), n, x(:,j+M1), 1, 1.0, x(:,j), 1)
             if(j+M2<=nBlock) CALL BLAS_gemv('n', n, n, -1.0, &
                  f2(:,:,j), n, x(:,j+M2), 1, 1.0, x(:,j), 1)
          enddo
       else
          !  x := U.x = x + F.x(j+1) + F1.x(j+M1) + F2.x(j+M2)
          do j = 1, nBlock-1
             call BLAS_gemv( &
                  'n', n, n, 1.0, f(:,:,j), n, x(:,j+1 ), 1, 1.0, x(:,j), 1)
             if(j+M1<=nBlock) call BLAS_gemv( &
                  'n', n, n, 1.0, f1(:,:,j), n, x(:,j+M1), 1, 1.0, x(:,j), 1)
             if(j+M2<=nBlock) call BLAS_gemv( &
                  'n', n, n, 1.0, f2(:,:,j),n, x(:,j+M2), 1, 1.0, x(:,j), 1)
          end do
       end if
    end if
#endif
    ! call timing_stop('Uhepta')

  end subroutine Uhepta
  !============================================================================
  subroutine Lhepta(nBlock, N, M1, M2, x, d, e, e1, e2)
    !$acc routine vector

    ! This routine multiplies x with the lower triangular matrix L^{-1},
    ! which must have been constructed in subroutine prehepta:
    !
    ! x' = L^{-1}.x = D^{-1}.(x - E2.x'(j-M2) - E1.x'(j-M1) - E.x'(j-1))
    !
    ! For penta-diagonal matrix, set M2=nBlock and e2 can be omitted.
    ! For tri-diagonal matrix set M1=M2=nBlock and e1,e2 can be omitted.
    !
    ! For penta- or hepta-diagonal matrix, a missing e1 means that only
    ! the immediate diagonal is used in the preconditioner (optimal for GPU)

    integer, intent(in) :: n, M1, M2, nBlock
    real, intent(inout) :: x(n,nBlock)
    real, intent(in)    :: d(n,n,nBlock), e(n,n,nBlock)
    real, intent(in), optional :: e1(n,n,nBlock), e2(n,n,nBlock)

    !     Description of arguments:
    !
    ! nBlock:  Number of diagonal blocks.
    ! N:        the size of the blocks.
    ! M1:       distance of blocks to the main diagonal blocks.
    !           Set M1=nBlock for block tri-diagonal matrix!
    ! M2:       distance of outer-most blocks to main diagonal blocks.
    !           Set M2=nBlock for block tri- and penta-diagonal matrices!
    !
    ! x:        On input, the vector to be multiplied with L^{-1}.
    !           On output, the result of L^{-1}.x
    !
    ! d, e, e1, e2
    !           The matrix L is in four block diagonals.
    !
    !           The blocks are stored as follows:
    !           d(j): j=1..nBlock     Contains inverse of diagonal of L
    !                                 where L is from the incomplete LU
    !           e(j): j=2..nBlock     sub diagonal blocks.
    !           e1(j): j=M1+1..nBlock Blocks in the lower-triangular part with
    !                                 distance M1 from the main diagonal.
    !                                 Omit for block tri-diagonal matrix!
    !           e2(j): j=M2+1..nBlock Blocks in the lower-triangular part with
    !                                 distance M2 from the main diagonal.
    !                                 Omit for block tri/penta-diagonal matrix!
    ! External subroutine: DGEMV, BLAS level two Matrix-Vector Product.

    real, allocatable :: work(:)

    integer :: j
    !--------------------------------------------------------------------------
    ! call timing_start('Lhepta')
    if(n == 1)then
       if(.not.present(e1))then
          call lower_hepta_scalar_fast(nBlock, M1, x, d, e)
       else
          call lower_hepta_scalar(nBlock, M1, M2, x, d, e, e1, e2)
       end if
       RETURN
    end if

#ifndef _OPENACC
    ! Allocate arrays that used to be automatic
    allocate(work(N))
    if(n <= 20)then
       ! F90 version
       do j = 1, nBlock
          work = x(:,j)
          if (j>M2) then
             work = work                        &
                  - matmul( e(:,:,j),x(:,j-1 )) &
                  - matmul(e1(:,:,j),x(:,j-M1)) &
                  - matmul(e2(:,:,j),x(:,j-M2))
          else if (j>M1) then
             work = work                        &
                  - matmul( e(:,:,j),x(:,j-1 )) &
                  - matmul(e1(:,:,j),x(:,j-M1))
          else if (j>1) then
             work = work                        &
                  - matmul( e(:,:,j),x(:,j-1 ))
          end if
          x(:,j) = matmul( d(:,:,j),work)
       end do
    else
       ! BLAS VERSION
       do j = 1, nBlock

          call BLAS_gemv('n', N, N, 1.0, d(:,:,j) ,N, x(:,j), 1, 0.0, work, 1)
          if(j > 1 ) call BLAS_gemv( &
               'n', n, n, -1.0, e(:,:,j) ,n, x(:,j-1), 1, 1.0, work, 1)
          if(j > M1) call BLAS_gemv( &
               'n', n, n, -1.0, e1(:,:,j), n, x(:,j-M1), 1, 1.0, work, 1)
          if(j > M2) call BLAS_gemv( &
               'n', n, n, -1.0, e2(:,:,j), n, x(:,j-M2), 1, 1.0, work, 1)
          call BLAS_copy(n, work, 1, x(:,j), 1)

       enddo
    end if
    deallocate(work)
#endif
    ! call timing_stop('Lhepta')

  end subroutine Lhepta
  !============================================================================
  subroutine upper_hepta_scalar(IsInverse, nBlock, m1, m2, x, f, f1, f2)
    !$acc routine vector

    ! G. Toth, 2009

    ! This routine multiplies x with the upper triagonal U or U^{-1}
    ! which must have been constructed in subroutine prehepta.
    !
    ! For penta-diagonal matrix, set M2=nBlock and f2 can be omitted.
    ! For tri-diagonal matrix set M1=M2=nBlock and f1,f2 can be omitted.

    logical, intent(in) :: IsInverse
    integer, intent(in) :: m1, m2, nBlock
    real, intent(inout) :: x(nBlock)
    real, intent(in)    :: f(nBlock)
    real, intent(in), optional :: f1(nBlock), f2(nBlock)

    integer :: j
    !--------------------------------------------------------------------------
    if(IsInverse)then
       !  x' := U^{-1}.x = x - F.x'(j+1) - F1.x'(j+M1) - F2.x'(j+M2)
       !$acc loop seq
       do j = nBlock-1, 1, -1
          !  x' := U^{-1}.x = x - F.x'(j+1) - F1.x'(j+M1) - F2.x'(j+M2)
          if (j+M2<=nBlock) then
             x(j) = x(j) - f(j)*x(j+1) - f1(j)*x(j+M1) - f2(j)*x(j+M2)
          else if(j+M1 <= nBlock) then
             x(j) = x(j) - f(j)*x(j+1) - f1(j)*x(j+M1)
          else
             x(j) = x(j) - f(j)*x(j+1)
          end if
       end do
    else
       !  x := U.x = x + F.x(j+1) + F1.x(j+M1) + F2.x(j+M2)
       !$acc loop seq
       do j = 1, nBlock-1
          if (j + M2 <= nBlock) then
             x(j) = x(j) + f(j)*x(j+1) + f1(j)*x(j+M1) + f2(j)*x(j+M2)
          else if (j+M1 <= nBlock) then
             x(j) = x(j) + f(j)*x(j+1) + f1(j)*x(j+M1)
          else
             x(j) = x(j) + f(j)*x(j+1)
          end if
       end do
    end if

  end subroutine upper_hepta_scalar
  !============================================================================
  subroutine upper_hepta_scalar_fast(IsInverse, nBlock, m1, x, f)
    !$acc routine vector

    ! Multiply x with the upper triagonal U or U^{-1}
    ! restricted to immediate upper diagonal "f". This allows parallelizing
    ! nBlock/M1 threads for the case when M1 < nBlock as every M1-th element
    ! of f is 0 as there is no "lower neighbor" at the edge of the grid block.

    logical, intent(in) :: IsInverse
    integer, intent(in) :: m1, nBlock
    real, intent(inout) :: x(nBlock)
    real, intent(in)    :: f(nBlock)

    integer:: i, j
    !--------------------------------------------------------------------------
    if(IsInverse)then
       !  x' := U^{-1}.x = x - F.x'(j+1)
       !$acc loop vector
       do i = nBlock, 1, -M1
          !$acc loop seq
          do j = i-1, i-M1+1, -1
             x(j) = x(j) - f(j)*x(j+1)
          end do
       end do
    else
       !  x := U.x = x + F.x(j+1)
       !$acc loop vector
       do i = 1, nBlock, M1
          !$acc loop seq
          do j = i, i+M1-2
             x(j) = x(j) + f(j)*x(j+1)
          end do
       end do
    end if

  end subroutine upper_hepta_scalar_fast
  !============================================================================
  subroutine lower_hepta_scalar(nBlock, M1, M2, x, d, e, e1, e2)
    !$acc routine vector

    ! Multiply x with the lower triangular matrix L^{-1},
    ! which must have been constructed in subroutine prehepta.
    !
    ! For penta-diagonal matrix, set M2=nBlock and e2 can be omitted.
    ! For tri-diagonal matrix set M1=M2=nBlock and e1,e2 can be omitted.

    integer, intent(in):: m1, m2, nBlock
    real, intent(inout):: x(nBlock)
    real, intent(in)   :: d(nBlock), e(nBlock)
    real, intent(in), optional :: e1(nBlock), e2(nBlock)

    real:: Work1

    integer :: j
    !--------------------------------------------------------------------------
    ! x' = L^{-1}.x = D^{-1}.(x - E2.x'(j-M2) - E1.x'(j-M1) - E.x'(j-1))

    !$acc loop seq
    do j = 1, nBlock
       work1 = x(j)
       if (j > M2) then
          work1 = work1 - e(j)*x(j-1) - e1(j)*x(j-M1) - e2(j)*x(j-M2)
       else if (j > M1) then
          work1 = work1 - e(j)*x(j-1) - e1(j)*x(j-M1)
       else if (j > 1) then
          work1 = work1 - e(j)*x(j-1)
       end if
       x(j) = d(j)*work1
    end do

  end subroutine lower_hepta_scalar
  !============================================================================
  subroutine lower_hepta_scalar_fast(nBlock, M1, x, d, e)
    !$acc routine vector

    ! This routine multiplies x with the lower triangular matrix L^{-1},
    ! restricted to immediate lower diagonal "e". This allows parallelization
    ! nBlock/M1 threads for the case when M1 < nBlock as every M1-th element
    ! of e is 0 as there is no "lower neighbor" at the edge of the grid block.

    integer, intent(in):: nBlock, m1
    real, intent(inout):: x(nBlock)
    real, intent(in)   :: d(nBlock), e(nBlock)

    integer :: i, j
    !--------------------------------------------------------------------------
    ! x' = L^{-1}.x = D^{-1}.(x - E.x'(j-1)) = d.(x - e.x')
    ! Prehepta calculated d = D^{-1}

    !$acc loop vector
    do i = 1, nBlock, M1
       x(i) = d(i)*x(i)
       !$acc loop seq
       do j = i+1, i+M1-1
          x(j) = d(j)*(x(j) - e(j)*x(j-1))
       end do
    end do

  end subroutine lower_hepta_scalar_fast
  !============================================================================
  subroutine multiply_dilu(nBlock, n, m1, m2, x, d, e, f, e1, f1, e2, f2)

    ! G. Toth, 2009

    ! This routine multiplies x with U^{-1} L^{-1}
    ! which must have been constructed in subroutine prehepta.
    ! DILU makes the assumption that the off-diagonal blocks are diagonal!
    !
    ! For block penta-diagonal matrix, set M2=nBlock and omit e2, f2
    ! For block tri-diagonal matrix set M1=M2=nBlock and omit e1, f1, e2, f2

    integer, intent(in)                    :: nBlock, n, m1, m2
    real, intent(inout)                    :: x(n,nBlock)
    real, intent(in), dimension(n,n,nBlock):: d, e, f
    real, intent(in), dimension(n,n,nBlock), optional :: e1, f1, e2, f2

    real, allocatable :: x_V(:)
    integer :: i, j
    !--------------------------------------------------------------------------

    ! x' = L^{-1}.x = D^{-1}.(x - E2.x'(j-M2) - E1.x'(j-M1) - E.x'(j-1))
    allocate(x_V(n))

    do j = 1, nBlock
       x_V = x(:,j)
       if (j > M2) then
          do i = 1, n
             x_V(i) = x_V(i) - e(i,i,j)*x(i,j-1) - e1(i,i,j)*x(i,j-M1) &
                  - e2(i,i,j)*x(i,j-M2)
          end do
       else if (j > M1) then
          do i = 1, n
             x_V(i) = x_V(i) - e(i,i,j)*x(i,j-1) - e1(i,i,j)*x(i,j-M1)
          end do
       else if (j > 1) then
          do i = 1, n
             x_V(i) = x_V(i) - e(i,i,j)*x(i,j-1)
          end do
       end if
       x(:,j) = matmul(d(:,:,j), x_V)
    end do

    !  x' := U^{-1}.x = x - D^-1 [F.x'(j+1) - F1.x'(j+M1) - F2.x'(j+M2)]
    do j = nBlock-1, 1, -1
       if (j + M2 <= nBlock) then
          do i = 1, n
             x_V(i) = f(i,i,j)*x(i,j+1) + f1(i,i,j)*x(i,j+M1) &
                  + f2(i,i,j)*x(i,j+M2)
          end do
       else if(j + M1 <= nBlock) then
          do i = 1, n
             x_V(i) = f(i,i,j)*x(i,j+1) + f1(i,i,j)*x(i,j+M1)
          end do
       else
          do i = 1, n
             x_V(i) = f(i,i,j)*x(i,j+1)
          end do
       end if
       x(:,j) = x(:,j) - matmul(d(:,:,j), x_V)

    end do
    deallocate(x_V)

  end subroutine multiply_dilu
  !============================================================================
  subroutine multiply_block_jacobi(nBlock, nVar, x_VB, d_VVB)

    ! G. Toth, 2009

    ! This routine multiplies x with the already inverted diagonal blocks.

    integer, intent(in)   :: nBlock, nVar
    real,    intent(inout):: x_VB(nVar, nBlock)
    real,    intent(in)   :: d_VVB(nVar, nVar, nBlock)

    integer :: iBlock
    !--------------------------------------------------------------------------
    !$omp parallel do
    do iBlock = 1, nBlock
       x_VB(:,iBlock) = matmul(d_VVB(:,:,iBlock), x_VB(:,iBlock))
    end do
    !$omp end parallel do

  end subroutine multiply_block_jacobi
  !============================================================================
  subroutine get_precond_matrix(PrecondParam, nVar, nDim, nI, nJ, nK, nDiag, &
       a_II)
    !$acc routine vector

    ! Create approximate L-U decomposition of a_II matrix that corresponds to
    ! an nDim dimensional grid with nI*nJ*nK cells and nVar variables per cell.

    real, intent(in):: PrecondParam ! see description in subroutine prehepta
    integer, intent(in):: nVar   ! number of variables per cell
    integer, intent(in):: nDim   ! number of dimensions 1, 2 or 3
    integer, intent(in):: nI     ! number of cells in dim 1
    integer, intent(in):: nJ     ! number of cells in dim 2
    integer, intent(in):: nK     ! number of cells in dim 3
    integer, intent(in):: nDiag  ! number of diagonals
    real, intent(inout):: a_II(nVar*nVar*nI*nJ*nK,nDiag) ! Precond matrix

    character(len=*), parameter:: NameSub = 'get_precond_matrix'
    !--------------------------------------------------------------------------
    select case(nDim)
    case(1)
       ! Tridiagonal case
       call prehepta(nI, nVar, nI, nI, PrecondParam, &
            a_II(1,1), a_II(1,2), a_II(1,3))
    case(2)
       if(nDiag == 3)then
          ! BILU1 case
          call prehepta(nI*nJ, nVar, nI, nI*nJ, PrecondParam, &
               a_II(1,1), a_II(1,2), a_II(1,3))
       else
          ! Pentadiagonal case
          call prehepta(nI*nJ, nVar, nI, nI*nJ, PrecondParam, &
               a_II(1,1), a_II(1,2), a_II(1,3), a_II(1,4), a_II(1,5))
       end if
    case(3)
       if(nDiag == 3)then
          call prehepta(nI*nJ*nK, nVar, nI, nI*nJ, PrecondParam, &
               a_II(1,1), a_II(1,2), a_II(1,3))
       else
          ! Heptadiagonal case
          call prehepta(nI*nJ*nK, nVar, nI, nI*nJ, PrecondParam, &
               a_II(1,1), a_II(1,2), a_II(1,3), a_II(1,4), a_II(1,5), &
               a_II(1,6), a_II(1,7))
       end if
    case default
       write(*,*)'ERROR in ', NameSub, ' nDim=', nDim
#ifndef _OPENACC
       call CON_stop(NameSub//': invalid value for nDim')
#endif
    end select

  end subroutine get_precond_matrix
  !============================================================================
  subroutine multiply_left_precond(TypePrecond, TypePrecondSide, &
       nVar, nDim, nI, nJ, nK, nDiag, a_II, x_I)
    !$acc routine vector

    ! Multiply x_I with the left preconditioner matrix using the
    ! a_II matrix which was obtained with "get_precond_matrix"
    ! TypePrecond defines which type of preconditioner is used
    ! TypePrecondSide defines which side the preconditioner is applied

    character(len=*), intent(in):: TypePrecond     ! DILU, BILU, MBILU
    character(len=*), intent(in):: TypePrecondSide ! left, right, symm

    integer, intent(in):: nVar   ! number of variables per cell
    integer, intent(in):: nDim   ! number of dimensions 1, 2 or 3
    integer, intent(in):: nI     ! number of cells in dim 1
    integer, intent(in):: nJ     ! number of cells in dim 2
    integer, intent(in):: nK     ! number of cells in dim 3
    integer, intent(in):: nDiag  ! number of diagonals
    real,    intent(in):: a_II(nVar*nVar*nI*nJ*nK,nDiag) ! Precond matrix
    real, intent(inout):: x_I(nVar*nI*nJ*nK)             ! Vector of vars

    character(len=*), parameter:: NameSub = 'multiply_left_precond'
    !--------------------------------------------------------------------------
    if(TypePrecondSide == 'right') RETURN

#ifndef _OPENACC
    select case(TypePrecondSide)
    case('left', 'symmetric')
    case default
       call CON_stop(NameSub// &
            ': unknown value for TypePrecondSide='//TypePrecondSide)
    end select

    if(nDim < 1 .or. nDim > 3)then
       write(*,*)'ERROR in ', NameSub, ' nDim=', nDim
       call CON_stop(NameSub//': invalid value for nDim')
    end if
#endif

    select case(TypePrecond)
#ifndef _OPENACC
    case('BLOCKJACOBI')
       ! Multiply with the inverted diagonal blocks of the matrix
       call multiply_block_jacobi(nI*nJ*nK, nVar, x_I, &
            a_II(1,1))
    case('DILU')
       select case(nDim)
       case(1)
          ! Tridiagonal case
          call multiply_dilu(nI, nVar, nI, nI, x_I, &
               a_II(1,1), a_II(1,2), a_II(1,3))
       case(2)
          ! Pentadiagonal case
          call multiply_dilu(nI*nJ, nVar, nI, nI*nJ, x_I, &
               a_II(1,1), a_II(1,2), a_II(1,3), a_II(1,4), a_II(1,5))
       case(3)
          ! Heptadiagonal case
          call multiply_dilu(nI*nJ*nK, nVar, nI, nI*nJ, x_I, &
               a_II(1,1), a_II(1,2), a_II(1,3), a_II(1,4), a_II(1,5), &
               a_II(1,6), a_II(1,7))
       end select
#endif
    case('BILU1')
       ! Multiply with L^-1 restricted to main + one lower diagonal
       call lower_hepta_scalar_fast(nI*nJ*nK, nI, x_I, a_II(1,1), a_II(1,2))
       ! Multiply with U^-1 restricted to one upper diagonal
       if(TypePrecondSide == 'left')  &
            call upper_hepta_scalar_fast(.true., nI*nJ*nK, nI, x_I, a_II(1,3))
    case('BILU', 'MBILU')
       ! Multiply with L^-1 from the LU decomposition
       select case(nDim)
       case(1)
          ! Tridiagonal case
          call Lhepta(nI, nVar, nI, nI, x_I, &
               a_II(1,1), a_II(1,2))
       case(2)
          ! Pentadiagonal case
          call Lhepta(nI*nJ, nVar, nI, nI*nJ, x_I, &
               a_II(1,1), a_II(1,2), a_II(1,4))
       case(3)
          ! Heptadiagonal case
          call Lhepta(nI*nJ*nK, nVar, nI, nI*nJ, x_I, &
               a_II(1,1), a_II(1,2), a_II(1,4), a_II(1,6))
       end select

       if(TypePrecondSide == 'left')then
          ! Multiply with U^-1 from the LU decomposition
          select case(nDim)
          case(1)
             ! Tridiagonal case
             call Uhepta(.true., nI, nVar, nI, nI, x_I, &
                  a_II(1,3))
          case(2)
             ! Pentadiagonal case
             call Uhepta(.true., nI*nJ, nVar, nI, nI*nJ, x_I, &
                  a_II(1,3), a_II(1,5))
          case(3)
             ! Heptadiagonal case
             call Uhepta(.true., nI*nJ*nK, nVar, nI, nI*nJ, x_I, &
                  a_II(1,3), a_II(1,5), a_II(1,7))
          end select
       end if
    case default
#ifndef _OPENACC
       call CON_stop(NameSub//': unknown value for TypePrecond='//TypePrecond)
#endif
    end select

  end subroutine multiply_left_precond
  !============================================================================
  subroutine multiply_right_precond(TypePrecond, TypePrecondSide, &
       nVar, nDim, nI, nJ, nK, nDiag, a_II, x_I)
    !$acc routine vector

    ! Multiply x_I with the right preconditioner matrix using the
    ! a_II matrix which was obtained with "get_precond_matrix"
    ! TypePrecond defines which type of preconditioner is used
    ! TypePrecondSide defines which side the preconditioner is applied

    character(len=*), intent(in):: TypePrecond     ! DILU, BILU, MBILU
    character(len=*), intent(in):: TypePrecondSide ! left, right, symm

    integer, intent(in):: nVar   ! number of variables per cell
    integer, intent(in):: nDim   ! number of dimensions 1, 2 or 3
    integer, intent(in):: nI     ! number of cells in dim 1
    integer, intent(in):: nJ     ! number of cells in dim 2
    integer, intent(in):: nK     ! number of cells in dim 3
    integer, intent(in):: nDiag  ! number of diagonals
    real,    intent(in):: a_II(nVar*nVar*nI*nJ*nK,nDiag) ! Precond matrix
    real, intent(inout):: x_I(nVar*nI*nJ*nK)             ! Vector of vars

    character(len=*), parameter:: NameSub = 'multiply_right_precond'
    !--------------------------------------------------------------------------
    if(TypePrecondSide == 'left') RETURN

#ifndef _OPENACC
    select case(TypePrecondSide)
    case('right', 'symmetric')
    case default
       call CON_stop(NameSub// &
            ': unknown value for TypePrecondSide='//TypePrecondSide)
    end select

    if(nDim < 1 .or. nDim > 3)then
       write(*,*)'ERROR in ', NameSub, ' nDim=', nDim
       call CON_stop(NameSub//': invalid value for nDim')
    end if
#endif

    select case(TypePrecond)
    case('DILU')
       ! Currently DILU can only be applied from the left
       RETURN
    case('BILU', 'MBILU')
       if(TypePrecondSide == 'right')then
          ! Multiply with L^-1 from the LU decomposition
          select case(nDim)
          case(1)
             ! Tridiagonal case
             call Lhepta(nI, nVar, nI, nI, x_I, &
                  a_II(1,1), a_II(1,2))
          case(2)
             ! Pentadiagonal case
             call Lhepta(nI*nJ, nVar, nI, nI*nJ, x_I, &
                  a_II(1,1), a_II(1,2), a_II(1,4))
          case(3)
             ! Heptadiagonal case
             call Lhepta(nI*nJ*nK, nVar, nI, nI*nJ, x_I, &
                  a_II(1,1), a_II(1,2), a_II(1,4), a_II(1,6))
          end select
       end if
       ! Multiply with U^-1 from the LU decomposition
       select case(nDim)
       case(1)
          ! Tridiagonal case
          call Uhepta(.true., nI, nVar, nI, nI, x_I, &
               a_II(1,3))
       case(2)
          ! Pentadiagonal case
          call Uhepta(.true., nI*nJ, nVar, nI, nI*nJ, x_I, &
               a_II(1,3), a_II(1,5))
       case(3)
          ! Heptadiagonal case
          call Uhepta(.true., nI*nJ*nK, nVar, nI, nI*nJ, x_I, &
               a_II(1,3), a_II(1,5), a_II(1,7))
       end select
    case default
#ifndef _OPENACC
       call CON_stop(NameSub//': unknown value for TypePrecond='//TypePrecond)
#endif
    end select

  end subroutine multiply_right_precond
  !============================================================================
  subroutine precond_left_multiblock(Param, &
       nVar, nDim, nI, nJ, nK, nDiag, nBlock, Jac_VVCIB, x_I)

    ! Multiply x_I with the left preconditioner matrix

    type(LinearSolverParamType), intent(in):: Param

    integer, intent(in):: nVar   ! number of variables per cell
    integer, intent(in):: nDim   ! number of dimensions 1, 2 or 3
    integer, intent(in):: nI     ! number of cells in dim 1
    integer, intent(in):: nJ     ! number of cells in dim 2
    integer, intent(in):: nK     ! number of cells in dim 3
    integer, intent(in):: nDiag  ! number of diagonals
    integer, intent(in):: nBlock ! number of blocks

    ! Preconditioner matrix
    real, intent(in):: Jac_VVCIB(nVar,nVar,nI,nJ,nK,nDiag,nBlock)

    ! Vector of variables
    real, intent(inout):: x_I(nVar*nI*nJ*nK*nBlock)

    integer:: nVarIJK, iBlock

    character(len=*), parameter:: NameSub = 'precond_left_multiblock'
    !--------------------------------------------------------------------------
    if(.not. Param%DoPrecond) RETURN
    if(Param%TypePrecondSide == 'right') RETURN

    nVarIJK = nVar*nI*nJ*nK
    !$omp parallel do
    !$acc parallel loop gang independent
    do iBlock = 1, nBlock
       call multiply_left_precond( &
            Param%TypePrecond, Param%TypePrecondSide,&
            nVar, nDim, nI, nJ, nK, nDiag, Jac_VVCIB(1,1,1,1,1,1,iBlock), &
            x_I(nVarIJK*(iBlock-1) + 1))
    end do
    !$omp end parallel do

  end subroutine precond_left_multiblock
  !============================================================================
  subroutine precond_right_multiblock(Param, &
       nVar, nDim, nI, nJ, nK, nDiag, nBlock, Jac_VVCIB, x_I)

    ! Multiply x_I with the right preconditioner matrix using the

    type(LinearSolverParamType), intent(in):: Param

    integer, intent(in):: nVar   ! number of variables per cell
    integer, intent(in):: nDim   ! number of dimensions 1, 2 or 3
    integer, intent(in):: nI     ! number of cells in dim 1
    integer, intent(in):: nJ     ! number of cells in dim 2
    integer, intent(in):: nK     ! number of cells in dim 3
    integer, intent(in):: nDiag  ! number of diagonals
    integer, intent(in):: nBlock ! number of blocks

    ! Preconditioner matrix
    real, intent(in):: Jac_VVCIB(nVar,nVar,nI,nJ,nK,nDiag,nBlock)

    ! Vector of variables
    real, intent(inout):: x_I(nVar*nI*nJ*nK*nBlock)

    integer:: nVarIJK, iBlock

    character(len=*), parameter:: NameSub = 'precond_right_multiblock'
    !--------------------------------------------------------------------------
    if(.not. Param%DoPrecond) RETURN
    if(Param%TypePrecondSide == 'left') RETURN

    nVarIJK = nVar*nI*nJ*nK
    !$omp parallel do
    !$acc parallel loop gang independent
    do iBlock = 1, nBlock
       call multiply_right_precond( &
            Param%TypePrecond, Param%TypePrecondSide,&
            nVar, nDim, nI, nJ, nK, nDiag, Jac_VVCIB(1,1,1,1,1,1,iBlock), &
            x_I(nVarIJK*(iBlock-1) + 1))
    end do
    !$omp end parallel do

  end subroutine precond_right_multiblock
  !============================================================================
  subroutine multiply_initial_guess(nVar, nDim, nI, nJ, nK, nDiag, a_II, x_I)

    ! Multiply x_I with the upper triangular part of
    ! a_II matrix which was obtained with "get_precond_matrix"
    ! This is only needed if x_I is not zero initially and
    ! a symmetric preconditioning is used.
    ! Multiplying with L is not implemented, so non-zero initial guess
    ! cannot be combined with 'right' preconditioning.

    integer, intent(in):: nVar   ! number of variables per cell
    integer, intent(in):: nDim   ! number of dimensions 1, 2 or 3
    integer, intent(in):: nI     ! number of cells in dim 1
    integer, intent(in):: nJ     ! number of cells in dim 2
    integer, intent(in):: nK     ! number of cells in dim 3
    integer, intent(in):: nDiag  ! number of diagonals
    real,    intent(in):: a_II(nVar*nVar*nI*nJ*nK,nDiag) ! Precond matrix
    real, intent(inout):: x_I(nVar*nI*nJ*nK)             ! Vector of vars

    ! Multiply with U^-1 from the LU decomposition
    character(len=*), parameter:: NameSub = 'multiply_initial_guess'
    !--------------------------------------------------------------------------
    select case(nDim)
    case(1)
       ! Tridiagonal case
       call Uhepta(.false., nI, nVar, nI, nI, x_I, &
            a_II(1,3))
    case(2)
       if(nDiag == 3)then
          call Uhepta(.false., nI*nJ, nVar, nI, nI*nJ, x_I, &
               a_II(1,3))
       else
          ! Pentadiagonal case
          call Uhepta(.false., nI*nJ, nVar, nI, nI*nJ, x_I, &
               a_II(1,3), a_II(1,5))
       end if
    case(3)
       if(nDiag == 3)then
          call Uhepta(.false., nI*nJ*nK, nVar, nI, nI*nJ, x_I, &
               a_II(1,3))
       else
          ! Heptadiagonal case
          call Uhepta(.false., nI*nJ*nK, nVar, nI, nI*nJ, x_I, &
               a_II(1,3), a_II(1,5), a_II(1,7))
       end if
    case default
       write(*,*)'ERROR in ', NameSub, ' nDim=', nDim
       call CON_stop(NameSub//': invalid value for nDim')
    end select

  end subroutine multiply_initial_guess
  !============================================================================
  subroutine solve_linear_multiblock(Param, &
       nVar, nDim, nI, nJ, nK, nDiag, nBlock, iComm, impl_matvec, Rhs_I, x_I, &
       DoTest, Jac_VVCIB, JacobiPrec_I, cg_precond, hypre_precond)

    type(LinearSolverParamType), intent(inout):: Param
    integer, intent(in):: nVar       ! Number of impl. variables/cell
    integer, intent(in):: nDim       ! Number of spatial dimensions
    integer, intent(in):: nI, nJ, nK ! Number of cells in a block
    integer, intent(in):: nDiag      ! Number of diagonals in Jac_VVCIB
    integer, intent(in):: nBlock     ! Number of impl. grid blocks
    integer, intent(in):: iComm      ! MPI communicator for processors

    interface
       subroutine impl_matvec(Vec_I, MatVec_I, n)
         ! Calculate MatVec = Matrix.Vec
         implicit none
         integer, intent(in) :: n
         real,    intent(in) :: Vec_I(n)
         real,    intent(out):: MatVec_I(n)
       end subroutine impl_matvec
    end interface

    real, intent(inout):: Rhs_I(nVar*nI*nJ*nK*nBlock) ! RHS vector
    real, intent(inout):: x_I(nVar*nI*nJ*nK*nBlock)   ! Initial guess/solution

    logical, optional:: DoTest    ! show Krylov iterations and convergence

    real, intent(inout), optional:: &  ! Jacobian matrix --> preconditioner
         Jac_VVCIB(nVar,nVar,nI,nJ,nK,nDiag,nBlock)

    real, intent(inout), optional:: &  ! Point Jacobi preconditioner
         JacobiPrec_I(:)

    interface

       subroutine cg_precond(Vec_I, PrecVec_I, n)
         ! Preconditioner method for PCG
         implicit none
         integer, intent(in) :: n
         real,    intent(in) :: Vec_I(n)
         real,    intent(out):: PrecVec_I(n)
       end subroutine cg_precond

       subroutine hypre_precond(n, Vec_I)
         ! Preconditioner method for HYPRE AMG
         implicit none
         integer, intent(in) :: n
         real,    intent(inout) :: Vec_I(n)
       end subroutine hypre_precond

    end interface
    optional:: cg_precond, hypre_precond

    ! Local variables
    integer:: n, iBlock, i, j, k, iVar
    integer:: nVarIjk, nImpl
    !$acc declare create(nVarIjk, nImpl)

    character(len=*), parameter:: NameSub = 'solve_linear_multiblock'
    !--------------------------------------------------------------------------
    nVarIjk = nVar*nI*nJ*nK ! Number of variables per block
    !$acc update device(nVarIjk)

    ! Number of variables per processor
    nImpl   = nVarIjk*nBlock
    !$acc update device(nImpl)

    ! Make sure that left preconditioning is used when necessary
    select case(Param%TypePrecond)
    case('DILU', 'HYPRE', 'JACOBI', 'BLOCKJACOBI')
       Param%TypePrecondSide = 'left'
    end select

    if(Param%UseInitialGuess .and. Param%TypePrecondSide == 'right') &
         Param%TypePrecondSide = 'symmetric'

#ifdef _OPENACC
    ! So far, initial guess is not supported on GPU.
    !$acc parallel loop gang vector independent
    do i = 1, nImpl
       x_I(i) = 0
    end do
#else
    ! Initialize solution vector if needed
    if(.not.Param%UseInitialGuess) x_I = 0
#endif

    ! Get preconditioning matrix if required.
    ! Precondition RHS and initial guess (for symmetric prec only)
    if(Param%DoPrecond)then
       if(Param%TypePrecond == 'HYPRE')then
          call hypre_precond(nImpl, Rhs_I)
       elseif(Param%TypePrecond == 'JACOBI') then
          if(present(Jac_VVCIB))then
             n = 0
             !$omp parallel do private( n )
             do iBlock = 1, nBlock; do k = 1, nK; do j = 1, nJ; do i = 1, nI
                n = (iBlock-1)*nI*nJ*nk*nVar
                do iVar = 1, nVar
                   n = n + 1
                   JacobiPrec_I(n) = 1.0 / Jac_VVCIB(iVar,iVar,i,j,k,1,iBlock)
                end do
             end do; enddo; enddo; enddo
             !$omp end parallel do
          end if
          if(Param%TypeKrylov /= 'CG') &
               Rhs_I(1:nImpl) = JacobiPrec_I(1:nImpl)*Rhs_I(1:nImpl)
       else
          !$omp parallel do private( n )
          !$acc parallel loop gang private(n)
          do iBlock = 1, nBlock

             ! Preconditioning Jac_VVCIB matrix
             call get_precond_matrix( &
                  Param%PrecondParam, nVar, nDim, nI, nJ, nK, nDiag, &
                  Jac_VVCIB(1,1,1,1,1,1,iBlock))

             if(Param%TypeKrylov == 'CG') CYCLE

             ! Starting index in the linear arrays
             n = nVarIjk*(iBlock-1) + 1

             ! rhs --> P_L.rhs, where P_L=U^{-1}.L^{-1}, L^{-1}, or I
             ! for left, symmetric, and right preconditioning, respectively
             call multiply_left_precond(&
                  Param%TypePrecond, Param%TypePrecondSide, &
                  nVar, nDim, nI, nJ, nK, nDiag, &
                  Jac_VVCIB(1,1,1,1,1,1,iBlock), Rhs_I(n))

#ifndef _OPENACC
             ! Initial guess x --> P_R^{-1}.x where P_R^{-1} = I, U, LU for
             ! left, symmetric and right preconditioning, respectively
             ! Multiplication with LU is NOT implemented
             if(  Param%UseInitialGuess .and. &
                  Param%TypePrecondSide == 'symmetric') &
                  call multiply_initial_guess( &
                  nVar, nDim, nI, nJ, nK, nDiag, &
                  Jac_VVCIB(1,1,1,1,1,1,iBlock), x_I(n))
#endif
          end do
          !$omp end parallel do
       end if
    endif

    ! Initialize stopping conditions. Solver will return actual values.
    Param%nMatVec = Param%MaxMatvec
    Param%Error   = Param%ErrorMax

    if(DoTest)write(*,*)NameSub,': Before ', Param%TypeKrylov, &
         ' nMatVec, Error:', Param%nMatVec, Param%Error

    ! Solve linear problem
    ! call timing_start('krylov solver')
    select case(Param%TypeKrylov)
    case('BICGSTAB')
       call bicgstab(impl_matvec, Rhs_I, x_I, Param%UseInitialGuess, nImpl, &
            Param%Error, Param%TypeStop, Param%nMatvec, &
            Param%iError, DoTest, iComm)
    case('GMRES')
       call gmres(impl_matvec, Rhs_I, x_I, Param%UseInitialGuess, nImpl, &
            Param%nKrylovVector, &
            Param%Error, Param%TypeStop, Param%nMatvec, &
            Param%iError, DoTest, iComm)
    case('CG')
       if(.not. Param%DoPrecond)then
          call cg(impl_matvec, Rhs_I, x_I, Param%UseInitialGuess, nImpl,&
               Param%Error, Param%TypeStop, Param%nMatvec, &
               Param%iError, DoTest, iComm)
       elseif(Param%TypePrecond == 'JACOBI')then
          call cg(impl_matvec, Rhs_I, x_I, Param%UseInitialGuess, nImpl,&
               Param%Error, Param%TypeStop, Param%nMatvec, &
               Param%iError, DoTest, iComm, &
               JacobiPrec_I)
       else
          call cg(impl_matvec, Rhs_I, x_I, Param%UseInitialGuess, nImpl,&
               Param%Error, Param%TypeStop, Param%nMatvec, &
               Param%iError, DoTest, iComm, &
               preconditioner = cg_precond)
       end if
    case default
       call CON_stop(NameSub//': Unknown TypeKrylov='//Param%TypeKrylov)
    end select
    ! call timing_stop('krylov solver')

    ! Postprocessing: x = P_R.x where P_R = I, U^{-1}, U^{-1}L^{-1} for
    ! left, symmetric and right preconditioning, respectively
    if(Param%DoPrecond) call precond_right_multiblock(Param, &
         nVar, nDim, nI, nJ, nK, nDiag, nBlock, Jac_VVCIB, x_I)

    if(DoTest)write(*,*)NameSub,&
         ': After nMatVec, Error, iError=',&
         Param%nMatvec, Param%Error, Param%iError

    ! Converging without any iteration is not a real error, so set iError=0
    if(Param%iError==3) Param%iError=0

  end subroutine solve_linear_multiblock
  !============================================================================
  subroutine implicit_solver(ImplPar, DtImpl, DtExpl, nCell, nVar, State_GV, &
       calc_residual, update_boundary)

    ! Solve a linear system in 1D

    real,    intent(in) :: ImplPar, DtImpl, DtExpl
    integer, intent(in) :: nCell, nVar
    real, intent(inout) :: State_GV(-1:nCell+2, nVar)

    interface
       subroutine calc_residual(nOrder, Dt, nCell, nVar, State_GV, Resid_CV)
         implicit none
         integer, intent(in) :: nOrder, nCell, nVar
         real,    intent(in) :: Dt
         real,    intent(in) :: State_GV(-1:nCell+2, nVar)
         real,    intent(out):: Resid_CV(nCell, nVar)
       end subroutine calc_residual

       subroutine update_boundary(nCell, nVar, State_GV)
         implicit none
         integer, intent(in)    :: nCell, nVar
         real,    intent(inout) :: State_GV(-1:nCell+2, nVar)
       end subroutine update_boundary
    end interface

    real, parameter :: Eps = 1.e-6
    real, allocatable, dimension(:,:) :: &
         RightHand_CV, ResidOrig_CV, StateEps_GV, ResidEps_CV

    integer, parameter :: nDiag = 3
    real, allocatable :: Matrix_VVCI(:, :, :, :)

    real, allocatable :: Norm_V(:)  ! second norm of variables
    real, allocatable :: x_I(:)     ! linear vector of right hand side/unknowns
    real :: Coeff

    integer :: iVar, jVar, i, iStencil, iDiag, iX

    logical, parameter :: DoTest = .false.
    !--------------------------------------------------------------------------
    allocate(StateEps_GV(-1:nCell+2,nVar), &
         RightHand_CV(nCell,nVar), &
         ResidOrig_CV(nCell,nVar),  &
         ResidEps_CV(nCell,nVar), &
         Norm_V(nVar), &
         x_I(nCell*nVar), &
         Matrix_VVCI(nVar,nVar,nCell,nDiag))

    ! Make sure that ghost cells are up-to-date
    call update_boundary(nCell, nVar, State_GV)

    ! calculate right hand side
    call calc_residual(2, DtImpl, nCell, nVar, State_GV, RightHand_CV)

    ! Calculate the unperturbed residual with first order scheme
    call calc_residual(1, DtExpl, nCell, nVar, State_GV, ResidOrig_CV)

    ! Calculate the norm for the variables
    do iVar = 1, nVar
       Norm_V(iVar) = sqrt(sum(State_GV(1:nCell,iVar)**2)/nCell)
    end do

    if(DoTest)write(*,*)'Norm_V=',Norm_V,' Eps=',Eps

    ! Calculate the dR/dU matrix
    do jVar = 1, nVar
       do iStencil = 1, 3
          ! Get perturbed state
          StateEps_GV = State_GV
          do i = iStencil, nCell, 3
             StateEps_GV(i, jVar) = StateEps_GV(i, jVar) + Eps*Norm_V(jVar)
          end do
          call update_boundary(nCell, nVar, StateEps_GV)
          if(DoTest) call save_plot(iStep, Time, nCell, nVar, StateEps_GV)
          call calc_residual(1, DtExpl, nCell, nVar, StateEps_GV, ResidEps_CV)

          ! Jacobian is multiplied with -ImplPar*DtImpl
          Coeff= -Implpar*DtImpl/(Eps*Norm_V(jVar)*DtExpl)
          do i = 1, nCell
             iDiag = modulo(i-iStencil,3) + 1
             do iVar = 1, nVar
                Matrix_VVCI(iVar,jVar,i,iDiag) = &
                     Coeff*(ResidEps_CV(i,iVar) - ResidOrig_CV(i,iVar))
             end do
          end do

       end do
    end do

    if(DoTest)then
       write(*,'(a,/,3( 3f5.1,/ ))') 'Matrix(:,:,2,1)=', Matrix_VVCI(:,:,2,1)
       write(*,'(a,/,3( 3f5.1,/ ))') 'Matrix(:,:,2,2)=', Matrix_VVCI(:,:,2,2)
       write(*,'(a,/,3( 3f5.1,/ ))') 'Matrix(:,:,2,3)=', Matrix_VVCI(:,:,2,3)
    end if

    ! Add the diagonal part J = I - delta t*dR/dU
    do i = 1, nCell
       do iVar = 1, nVar
          Matrix_VVCI(iVar,iVar,i,1) = Matrix_VVCI(iVar,iVar,i,1) + 1.0
       end do
    end do

    ! L-U decomposition using block ILU (exact in 1D)
    call get_precond_matrix(real(bilu_), nVar, 1, nCell, 1, 1, 3, Matrix_VVCI)

    ! Put right hand side into a linear vector
    iX = 0
    do i = 1, nCell
       do iVar = 1, nVar
          iX=iX + 1
          x_I(iX) = RightHand_CV(i, iVar)
       end do
    end do

    ! x --> U^{-1}.L^{-1}.rhs = A^{-1}.rhs
    call multiply_left_precond('BILU', 'left', nVar, 1, nCell, 1, 1, &
         nDiag, Matrix_VVCI, x_I)

    ! Update the solution: U^n+1 = U^n + x
    iX = 0
    do i = 1, nCell
       do iVar = 1, nVar
          iX=iX + 1
          State_GV(i, iVar) = State_GV(i, iVar) + x_I(iX)
       end do
    end do

    call update_boundary(nCell, nVar, State_GV)

    deallocate(RightHand_CV, ResidOrig_CV, StateEps_GV, ResidEps_CV, &
         Norm_V, x_I, Matrix_VVCI)

  end subroutine implicit_solver
  !============================================================================
  subroutine save_plot(iStep, Time, nCell, nVar, State_GV)

    integer, intent(in) :: iStep
    real,    intent(in) :: Time
    integer, intent(in) :: nCell, nVar
    real,    intent(in) :: State_GV(-1:nCell+2, nVar)

    integer :: i
    !--------------------------------------------------------------------------
    write(UNITTMP_,'(a79)') 'linear solver test_hd11'
    write(UNITTMP_,'(i7,1pe13.5,3i3)') iStep,Time,1,1,nVar
    write(UNITTMP_,'(3i4)')            nCell
    write(UNITTMP_,'(100es13.5)')      Gamma
    write(UNITTMP_,'(a79)') 'x rho rhou p gamma'
    do i = 1, nCell
       write(UNITTMP_,'(100es18.10)') float(i), State_GV(i, rho_:p_)
    end do

  end subroutine save_plot
  !============================================================================
  subroutine calc_resid_test(nOrder, Dt, nCell, nVar, State_GV, Resid_CV)

    integer, intent(in) :: nOrder, nCell, nVar
    real,    intent(in) :: Dt
    real,    intent(in) :: State_GV(-1:nCell+2, nVar)
    real,    intent(out):: Resid_CV(nCell, nVar)

    integer :: i
    real :: uLeft, uRight
    !--------------------------------------------------------------------------
    do i = 1, nCell
       uRight = State_GV(i+1,rhou_)/State_GV(i+1,rho_)
       uLeft  = State_GV(i-1,rhou_)/State_GV(i-1,rho_)
       Resid_CV(i,rho_) = -Dt * Inv2Dx * &
            (State_GV(i+1,rhou_) - State_GV(i-1,rhou_))
       Resid_CV(i,rhou_) = -Dt * Inv2Dx * &
            ( uRight**2*State_GV(i+1,rho_) - uLeft**2*State_GV(i-1,rho_) &
            + State_GV(i+1,p_) - State_GV(i-1,p_) )
       ! dP/dt = -div(P*u) -(g-1)*P*Div(U)
       Resid_CV(i,p_) = -Dt * Inv2Dx * &
            ( uRight*State_GV(i+1,p_) - uLeft*State_GV(i-1,p_) &
            + (Gamma-1)*State_GV(i,p_)*(uRight - uLeft) )
    end do
  end subroutine calc_resid_test
  !============================================================================
  subroutine update_bound_test(nCell, nVar, State_GV)

    integer, intent(in)    :: nCell, nVar
    real,    intent(inout) :: State_GV(-1:nCell+2, nVar)

    ! periodic
    ! State_GV(-1,:)      = State_GV(nCell-1,:)
    ! State_GV( 0,:)      = State_GV(nCell,:)
    ! State_GV(nCell+1,:) = State_GV(1,:)
    ! State_GV(nCell+2,:) = State_GV(2,:)

    ! floating
    !--------------------------------------------------------------------------
    State_GV(-1,:)      = State_GV(1,:)
    State_GV( 0,:)      = State_GV(1,:)
    State_GV(nCell+1,:) = State_GV(nCell,:)
    State_GV(nCell+2,:) = State_GV(nCell,:)

  end subroutine update_bound_test
  !============================================================================
  subroutine test_linear_solver

    integer, parameter :: nStep=20
    integer, parameter :: nCell=51, nVar=3
    real, parameter :: DtExpl = 0.1, DtImpl = 1.0, ImplPar = 1.0

    real    :: State_GV(-1:nCell+2, nVar), StateOld_GV(-1:nCell+2, nVar)
    real    :: Resid_CV(nCell, nVar)

    ! initial condition
    !--------------------------------------------------------------------------
    State_GV(:,rho_) = 1.0; State_GV(5:10,rho_)=2
    State_GV(:,rhou_)= 2.0*State_GV(:,rho_);
    State_GV(:,p_)   = 3.0; State_GV(5:10,p_) = 6.0

    open(UNITTMP_, file='test_linear_solver.out', status='replace')
    Time = 0.0
    call save_plot(0, Time, nCell, nVar, State_GV)
    do iStep = 1, nStep

       StateOld_GV = State_GV

       Time  = Time + DtImpl

       if(.false.)then
          ! explicit
          call calc_resid_test(2, DtExpl, nCell, nVar, State_GV, Resid_CV)
          State_GV(1:nCell,:)=State_GV(1:nCell,:)+Resid_CV
          call update_bound_test(nCell, nVar, State_GV)
          call save_plot(iStep, Time, nCell, nVar, State_GV)
       else
          ! implicit
          call implicit_solver(ImplPar, DtImpl, DtExpl, nCell, nVar, &
               State_GV, calc_resid_test, update_bound_test)

          ! check
          call calc_resid_test(2, DtImpl, nCell, nVar, State_GV, Resid_CV)
          write(*,*)'iStep, Max error=',iStep, &
               maxval(abs(StateOld_GV(1:nCell,:)+Resid_CV-State_GV(1:nCell,:)))
          call save_plot(iStep, Time, nCell, nVar, State_GV)
       end if
    end do
    close(UNITTMP_)

  end subroutine test_linear_solver
  !============================================================================
end module ModLinearSolver
!==============================================================================
