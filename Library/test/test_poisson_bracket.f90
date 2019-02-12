!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!===========================TESTS============================================
module ModTestPoissonBracket
  use ModPoissonBracket, ONLY: explicit 
  use ModUtilities,      ONLY: CON_stop
  use ModNumConst,       ONLY: cTwoPi
  
  implicit none
contains
  subroutine test_poisson_bracket(tFinal)
    real, intent(in) :: tFinal
    !Misc:
    integer, parameter::  nPlot = 1,  nQ = 30,  nP = 1000
    integer           ::  iPlot, iQ, iP
    real, parameter :: qMax = 10.0, qMin = 0.01
    real:: MomentaRatio, MomentumMin, MomentumMax
    real:: VDF_G(-1:nQ+2, -1:nP+1), Momentum_I(nQ), Momentum2_I(0:nQ)
    !---------------------
    MomentaRatio = exp(log(qMax/qMin)/nQ)
    contains
      real function Hamiltonian(P2)
        real, intent(in) :: P2 ! momentum squared
        Hamiltonian = sqrt(1.0 + P2)
      end function Hamiltonian
  end subroutine test_poisson_bracket
end module ModTestPoissonBracket
!=============================================================================
!=============================================================================
program test_poisson_bracket
  use ModTestPoissonBracket, test => test_poisson_bracket
  implicit none
  call test(6.28318520)
end program test_poisson_bracket
