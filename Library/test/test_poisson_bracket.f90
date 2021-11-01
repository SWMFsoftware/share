!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!===========================TESTS============================================
module ModTestPoissonBracket
  use ModPoissonBracket, ONLY: explicit
  use ModUtilities,      ONLY: CON_stop
  use ModNumConst,       ONLY: cTwoPi
  use ModPlotFile,       ONLY: save_plot_file
  use ModConst
  implicit none

contains
  !============================================================================

  subroutine test_poisson_bracket(tFinal)
    real, intent(in) :: tFinal
    ! Misc:
    ! Gyration in a uniform magnetic field, nQ numper of points
    ! for the momentum grid, nP is number of point over azimuthal angle
    integer, parameter::  nQ = 30,  nP = 360
    ! Loop variables
    integer           ::  iQ, iStep
    ! Momentum max and min, in units of mc 
    real, parameter   :: qMax = 10.0, qMin = 0.01
    ! Mesh size in \phi
    real, parameter   :: DeltaPhi = cTwoPi/nP
    real :: MomentumRatio, MomentumMin, MomentumMax
    real :: VDF_G(-1:nQ+2, -1:nP+2), Volume_G(0:nQ+1, 0:nP+1)
    real :: Hamiltonian_N(-1:nQ+1, -1:nP+1)
    real :: LogMomentum_I(0:nQ+1), Momentum2_I(-1:nQ+1)
    real :: Time, Dt, Source_C(nQ,nP)
    !--------------------------------------------------------------------------
    MomentumRatio = exp(log(qMax/qMin)/nQ)
    MomentumMax =  qMin/MomentumRatio
    Momentum2_I(-1) = MomentumMax**2
    Hamiltonian_N(-1,:) = Hamiltonian(Momentum2_I(-1))
    do iQ = 0, nQ+1
       MomentumMin = MomentumMax
       MomentumMax = MomentumMin*MomentumRatio
       Momentum2_I( iQ) = MomentumMax**2
       Hamiltonian_N(iQ,:) = Hamiltonian(Momentum2_I(iQ))
       Volume_G(iQ,:) = 0.5*DeltaPhi*&
            (Momentum2_I(iQ) - Momentum2_I(iQ-1))
       LogMomentum_I(iQ)   = 0.50*log10(MomentumMin*MomentumMax)
    end do
    ! Initial distribution function
    VDF_G = 0.0; VDF_G(1:nQ, 172:189) = 1.0; Source_C = 0.0
    ! Compiutation
    Time = 0.0; iStep = 0
    do
       call explicit(nQ, nP, VDF_G, Volume_G, Source_C, &
            Hamiltonian_N,   &
            CFLIn=0.99, DtOut = Dt)
       iStep = iStep +1
       if(Time + Dt >= tFinal)then
          call explicit(nQ, nP, VDF_G, Volume_G, Source_C, &
               Hamiltonian_N,   &
               DtIn = tFinal - Time)
          VDF_G(1:nQ, 1:nP) = VDF_G(1:nQ, 1:nP) + Source_C
          EXIT
       else
          Time = Time + Dt
          VDF_G(1:nQ, 1:nP) = VDF_G(1:nQ, 1:nP) + Source_C
          ! Periodic boundary conditions
          VDF_G(1:nQ,-1:0 ) = VDF_G(1:nQ, nP-1:nP)
          VDF_G(1:nQ, nP+1:nP+2 ) = VDF_G(1:nQ, 1:2)
       end if
    end do
    call save_plot_file(NameFile='test_poisson.out', &
         TypeFileIn='ascii', TimeIn=tFinal, nStepIn = iStep, &
         NameVarIn='Phi LogMomentum VDF'  , &
         CoordMinIn_D=[0.50,      LogMomentum_I(1 )],&
         CoordMaxIn_D=[nP - 0.50, LogMomentum_I(nQ)],&
         StringFormatIn = '(4F10.3)'            ,&
         Coord2In_I = LogMomentum_I(1:nQ)           ,&
         VarIn_II = transpose(VDF_G(1:nQ,1:nP)))
  contains
    !==========================================================================
    real function Hamiltonian(P2)
      real, intent(in) :: P2 ! momentum squared
      !------------------------------------------------------------------------
      Hamiltonian = sqrt(1.0 + P2)
    end function Hamiltonian
    !==========================================================================
  end subroutine test_poisson_bracket
end module ModTestPoissonBracket
!==============================================================================
program test_program
  use ModTestPoissonBracket, ONLY: test_poisson_bracket
  use ModNumConst,           ONLY: cTwoPi
  ! use ModTestPoissonBracketAndScatter, ONLY: test_scatter
  implicit none

  !----------------------------------------------------------------------------
  call test_poisson_bracket(cTwoPi)
  ! call test_multipoisson_bracket(50.0)

  ! Test for comparing two diffusions is long. It should be repeated twice with
  ! Diffuornot = 0 or 1

  ! call test_scatter(50.0)
end program test_program
