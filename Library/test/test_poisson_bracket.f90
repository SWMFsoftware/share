!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!===========================TESTS==============================================
module ModDiffusion
  ! Solve the diffusion term in the Parker equation
  ! Revision history:
  ! Prototype:Sp/FLAMPA/src/SP_ModDiffusion.f90
  ! Adapted for the use in MFLAMPA (Dist_I is an input paramater,
  ! fixed contributions to M_I in the end points)-D.Borovikov, 2017
  ! Updated (identation, comments):  I.Sokolov, Dec.17, 2017

  implicit none

  PRIVATE
  public:: advance_diffusion1,tridiag

contains
  !============================================================================
  subroutine advance_diffusion1(Dt, n, Dist_I, F_I, DOuter_I, DInner_I)

    ! Solve the diffusion equation:
    !         f_t-D_outer(D_inner*f_x)_x=0,
    ! with zero Neumann boundary condition. The solution is advanced in time
    ! using fully implicit scheme.

    use ModNumConst, ONLY: cTiny

    real,   intent(in   ):: Dt     ! Time step
    integer,intent(in   ):: n      ! Number of meshes along the x-coordinate
    real,   intent(in   ):: Dist_I(n) ! Distance to the next mesh
    real,   intent(inout):: F_I(n) ! In:sol.to be advanced; Out:advanced sol

    ! Laplace multiplier and diffusion coefficient.
    real,   intent(in   ):: DOuter_I(n), DInner_I(n)

    ! Mesh spacing and face spacing.
    real                 :: DsMesh_I(2:n), DsFace_I(2:n-1)

    ! Main, upper, and lower diagonals.
    real, dimension(n)   :: Main_I,Upper_I,Lower_I, R_I
    integer:: i
    real:: Aux1,Aux2
    !--------------------------------------------------------------------------
    ! In M-FLAMPA D_I(i) is the distance between meshes i   and i+1
    ! while DsMesh_I(i) is the distance between centers of meshes
    ! i-1 and i. Therefore,
    DsMesh_I = 0.0
    DsFace_I = 0.0
    Lower_I = 0.0
    Main_I = 0.0
    Upper_I = 0.0
    R_I = 0.0
    do i=2,n
       DsMesh_I(i) = max(Dist_I(i-1),cTiny)
    end do
    ! Within the framework of finite volume method, the cell
    ! volume is used, which is proportional to  the distance between
    ! the faces bounding the volume with an index, i, which is half of
    ! sum of distance between meshes i-1 and i (i.e. D_I(i-1) and that
    ! between meshes i and i+1 (which is D_I(i)):
    do i=2,n-1
       DsFace_I(i) = max(0.5*(Dist_I(i) + Dist_I(i-1)),cTiny)
    end do

    ! In flux coordinates, the control volume associated with the
    ! given cell has a cross-section equal to (Magnetic Flux)/B,
    ! where the flux is a constant along the magnetic field line,
    ! set to one hereafter. Therefore, the diffusion equation has
    ! a following form:
    ! (DsFace_i/B_i)(f^(n+1) - f^n) = Flux_(i-1/2) - Flux_(i+1/2),
    ! where the particle density flux should be multiplied by the
    ! cross-section area too (magnetic flux factor is one!):
    !  Flux_(i-1/2) = (diffusion coefficient)_(i-1/2)/B_(i-1/2)*&
    !                 (f^(n+1)_(i-1) - f^(n+1)_i),
    !  The multiplier, DsFace_i/B_i, is denoted as DsFace_i/DOuter_i
    !  The face-centered combination,
    ! f^(n+1)_i-Dt*DOuter_I/DsFace_I*(&
    !     DInner_(i+1/2)*(f^(n+1)_(i+1)-f^(n+1)_i)/DsMesh_(i+1)-&
    !     DInner_(i-1/2)*(f^(n+1)_i -f^(n+1)_(i-1)/DsMesh_(i ))=f^n_i
    Main_I = 1.0
    ! For i=1:
    Aux1 = Dt*DOuter_I(1)*0.50*(DInner_I(1)+DInner_I(2))/&
         DsMesh_I(2)**2
    Main_I( 1) = Main_I(1)+Aux1
    Upper_I(1) = -Aux1
    ! For i=2,n-1:
    do i=2,n-1
       Aux1 = Dt*DOuter_I(i)*0.50*(DInner_I(i  ) + DInner_I(i+1))/&
            (DsMesh_I(i+1)*DsFace_I(i))
       Aux2 = Dt*DOuter_I(i)*0.50*(DInner_I(i-1) + DInner_I(i  ))/&
            (DsMesh_I(i  )*DsFace_I(i))
       Main_I(i)  = Main_I(i) + Aux1 + Aux2
       Upper_I(i) = -Aux1
       Lower_I(i) = -Aux2
    end do
    ! For i=n:
    Aux2 = Dt*DOuter_I(n)*0.50*(DInner_I(n-1) + DInner_I(n))/&
         DsMesh_I(n)**2
    Main_I( n) = Main_I(n) + Aux2
    Lower_I(n) = -Aux2
    ! Update the solution from f^(n) to f^(n+1):
    R_I = F_I
    call tridiag(n,Lower_I,Main_I,Upper_I,R_I,F_I)
  end subroutine advance_diffusion1
  !============================================================================
  subroutine tridiag(n, L_I, M_I, U_I, R_I, W_I)

    ! Solve tri-diagonal system of equations:
    !  ||m_1 u_1  0....        || ||w_1|| ||r_1||
    !  ||l_2 m_2 u_2...        || ||w_2|| ||r_2||
    !  || 0  l_3 m_3 u_3       ||.||w_3||=||r_3||
    !  ||...                   || ||...|| ||...||
    !  ||.............0 l_n m_n|| ||w_n|| ||r_n||
    ! From: Numerical Recipes, Chapter 2.6, p.40.

    ! input parameters
    integer,            intent(in):: n
    real, dimension(n), intent(in):: L_I, M_I ,U_I ,R_I
    ! Output parameters
    real, intent(out):: W_I(n)
    ! Misc
    integer:: j
    real:: Aux,Aux_I(2:n)
    !--------------------------------------------------------------------------
    if (M_I(1)==0.0) then
       write(*,*)' Error in tridiag: M_I(1)=0'
       stop
    end if
    Aux = M_I(1)
    W_I(1) = R_I(1)/Aux
    do j=2,n
       Aux_I(j) = U_I(j-1)/Aux
       Aux = M_I(j)-L_I(j)*Aux_I(j)
       if (Aux == 0.0) then
          write(*,*)'M_I(j), L_I(j), Aux_I(j) = ',&
               M_I(j),L_I(j),Aux_I(j)
          write(*,*)'  For j=',j
          write(*,*)'Tridiag failed'
          stop
       end if
       W_I(j) = (R_I(j)-L_I(j)*W_I(j-1))/Aux
    end do
    do j=n-1,1,-1
       W_I(j) = W_I(j)-Aux_I(j+1)*W_I(j+1)
    end do
  end subroutine tridiag
  !============================================================================
end module ModDiffusion
!==============================================================================
module ModTestPoissonBracket
  use ModPoissonBracket, ONLY: explicit
  use ModUtilities,      ONLY: CON_stop
  use ModNumConst,       ONLY: cTwoPi
  use ModPlotFile,       ONLY: save_plot_file
  use ModConst
  implicit none
  integer, parameter :: nX = 10000 !# of Lagrangian points
  integer, parameter :: nP = 40    ! 80  - for Fig5.Right Panel!# momentum bins
  real ::  VDFOld_G(-1:nX+2, -1:nP+2)
contains
  !============================================================================
  subroutine test_poisson_bracket(tFinal)

    real, intent(in) :: tFinal

    ! Misc:
    ! Gyration in a uniform magnetic field, nQ numper of points
    ! for the momentum grid, nP is number of point over azimuthal angle
    integer, parameter::  nQ = 300,  nP = 360

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
            CFLIn=0.99, DtOut = Dt, IsPeriodicIn_D=[.false.,.true.])
       iStep = iStep +1
       if(Time + Dt >= tFinal)then
          call explicit(nQ, nP, VDF_G, Volume_G, Source_C, &
               Hamiltonian_N,   &
               DtIn = tFinal - Time, IsPeriodicIn_D=[.false.,.true.])
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
         NameVarIn='-p_y  p_x VDF', &
         CoordMinIn_D=[10.0**LogMomentum_I(1 ), 0.50*DeltaPhi], &
         CoordMaxIn_D=[10.0**LogMomentum_I(nQ), cTwoPi - 0.50*DeltaPhi], &
         Coord1In_I = 10.0**LogMomentum_I(1:nQ), &
         VarIn_II = VDF_G(1:nQ,1:nP) )
  end subroutine test_poisson_bracket
  !==========================================================================
  real function Hamiltonian(P2)
    real, intent(in) :: P2 ! momentum squared
    !------------------------------------------------------------------------
    Hamiltonian = sqrt(1.0 + P2)
  end function Hamiltonian
  !============================================================================
  subroutine test_poisson_2d(tFinal)

    real, intent(in) :: tFinal

    ! Misc:
    ! Harmonic oscillators, nQ is the number of meshes over coordinate,
    ! nP is number of meshes overr generalized momentum
    integer, parameter::  nQ = 300,  nP = 300

    ! Loop variables
    integer           ::  iQ, iP, iStep

    ! Mesh size
    real, parameter   :: DeltaQ = 24.0/nQ, DeltaP = 24.0/nP
    real :: VDF_G(-1:nQ+2, -1:nP+2)
    real :: Volume_G(0:nQ+1, 0:nP+1)
    real :: Hamiltonian_N(-1:nQ+1, -1:nP+1)
    real :: Time, Dt, Source_C(nQ,nP)
    real :: NormL2Init, NormL2, EnergyInit, Energy, qNode, pNode, Q, P
    !--------------------------------------------------------------------------
    ! Control volume, for a uniform rectangular grid
    Volume_G = DeltaQ*DeltaP
    ! Hamiltonian at the nodes
    do iP = -1, nP+1
       pNode = DeltaP*(iP - nP/2)
       do iQ = -1, nQ+2
          qNode = DeltaQ*(iQ - nQ/2)
          Hamiltonian_N(iQ,iP) = hamiltonian(qNode**2 + pNode**2)
       end do
    end do
    ! Initial distribution function
    VDF_G = 0.0
    do iP = -1, nP+2
       P = DeltaP*(iP - nP/2 - 0.50)
       do iQ = -1, nQ+2
          Q = DeltaQ*(iQ - nQ/2 - 0.50)
          if(abs(P) < -Q*tan(cPi/20.0).and.P**2+Q**2<=100)VDF_G(iQ,iP) = 1.0
       end do
    end do
    Source_C = 0.0
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
          write(*,*)Time
       end if
    end do
    call save_plot_file(NameFile='test_poisson2d.out', &
         TypeFileIn='real8', TimeIn=tFinal, nStepIn = iStep, &
         NameVarIn='-Py    Px VDF', &
         CoordMinIn_D=[-12.0 + 0.50*DeltaQ, -12.0 +  0.50*DeltaP], &
         CoordMaxIn_D=[12.0 - 0.50*DeltaQ, 12.0 - 0.50*DeltaP], &
         VarIn_II = VDF_G(1:nQ,1:nP))

  end subroutine test_poisson_2d
  !===================================================================
  subroutine test_poisson_2d_smooth(tFinal)
    real, intent(in) :: tFinal

    ! Misc:
    ! Anharmonic oscillators, nQ is the number of meshes over coordinate,
    ! nP is number of meshes overr generalized momentum
    integer, parameter::  nQ = 300,  nP = 300

    ! Loop variables
    integer           ::  iQ, iP, iStep

    ! Mesh size
    real, parameter   :: DeltaQ = 24.0/nQ, DeltaP = 24.0/nP
    real :: VDF_G(-1:nQ+2, -1:nP+2), VDFFinal_G(-1:nQ+2,-1:nP+2)
    real :: Plot_VC(3,1:nQ,1:nP)
    real :: Volume_G(0:nQ+1, 0:nP+1)
    real :: Hamiltonian_N(-1:nQ+1, -1:nP+1)
    real :: Time, Dt, Source_C(nQ,nP), Error
    real :: NormL2Init, NormL2, EnergyInit, Energy, qNode, pNode, Q, P
    real, parameter :: WidthX = 6.0, WidthY = 6.0
    !--------------------------------------------------------------------------
    ! Control volume, for a uniform rectangular grid
    Volume_G = DeltaQ*DeltaP
    ! Hamiltonian at the nodes
    do iP = -1, nP+1
       pNode = DeltaP*(iP - nP/2)
       do iQ = -1, nQ+2
          qNode = DeltaQ*(iQ - nQ/2)
          Hamiltonian_N(iQ,iP) = hamiltonian(qNode**2 + pNode**2)
       end do
    end do
    ! Initial distribution function
    VDF_G = 0.0
    do iP = -1, nP+2
       P = DeltaP*(iP - nP/2 - 0.50)
       do iQ = -1, nQ+2
          Q = DeltaQ*(iQ - nQ/2 - 0.50)
          VDF_G(     iQ,iP) = initial_cap(Q,P)
          VDFFinal_G(iQ,iP) = final_cap(  Q,P)
       end do
    end do
    Plot_VC(1,:,:) = VDF_G(     1:nQ,1:nP)
    Plot_VC(2,:,:) = VDFFinal_G(1:nQ,1:nP)
    Source_C = 0.0
    ! Computation
    Time = 0.0; iStep = 0
    do
       call explicit(nQ, nP, VDF_G, Volume_G, Source_C, &
            Hamiltonian_N,   &
            CFLIn=0.50, DtOut = Dt)
       iStep = iStep +1
       if(Time + Dt >= tFinal)then
          call explicit(nQ, nP, VDF_G, Volume_G, Source_C, &
               Hamiltonian_N,   &
               DtIn = tFinal - Time)
          Time = tFinal
          write(*,*) Time
          VDF_G(1:nQ, 1:nP) = VDF_G(1:nQ, 1:nP) + Source_C
          EXIT
       else
          Time = Time + Dt
          VDF_G(1:nQ, 1:nP) = VDF_G(1:nQ, 1:nP) + Source_C
          write(*,*)Time
       end if
    end do
    Plot_VC(3,1:nQ,1:nP) = VDF_G(1:nQ, 1:nP)
    Error = sum(abs(VDF_G(1:nQ,1:nP) - VDFFinal_G(1:nQ,1:nP))*&
         Volume_G(1:nQ,1:nP))/(WidthX*WidthY)
    write(*,*)'Error=',Error
    call save_plot_file(NameFile='test_poisson_2d_smooth.out', &
         TypeFileIn='real8', TimeIn=tFinal, nStepIn = iStep, &
         NameVarIn='-Py    Px   VDFInit VDFAnal VDFFinal  Error', &
         CoordMinIn_D=[-12.0 + 0.50*DeltaQ, -12.0 +  0.50*DeltaP], &
         CoordMaxIn_D=[12.0 - 0.50*DeltaQ, 12.0 - 0.50*DeltaP], &
         ParamIn_I=[Error], &
         VarIn_VII = Plot_VC)
  contains
    !=============================
    real function initial_cap(x,y)
      real, intent(in) :: x,y
      real, parameter :: xCenter = -6.0, yCenter = 0.0
      !--------------------------------------------
      initial_cap = 0.0
      if(abs(x  - xCenter)<=WidthX/2.and.abs(y - yCenter)<=WidthY/2)then
         initial_cap = cos(cPi*(x - xCenter)/WidthX)**4*&
              cos(cPi*(y - yCenter)/WidthY)**4
      end if
    end function initial_cap
    !==========================================================================
    real function final_cap(x,y)
      real, intent(in) :: x,y
      real :: CosPhi,SinPhi,Phi, P, P2, OmegaT, xInit, yInit
      !----------------------
      P2  = x**2 + y**2
      OmegaT = tFinal/hamiltonian(P2)
      P  = sqrt(P2)
      CosPhi = x/P; SinPhi = y/P
      if(SinPhi>=0.0)then
         Phi =  acos(CosPhi)
      else
         Phi = -acos(CosPhi)
      end if
      xInit = P*cos(Phi + OmegaT)
      yInit = P*sin(Phi + OmegaT)
      final_cap = initial_cap(xInit,yInit)
    end function final_cap
    !==========================================================================
  end subroutine test_poisson_2d_smooth
  !============================================================================
  subroutine test_energy_conservation(tFinal)

    real, intent(in) :: tFinal

    ! Misc:
    ! Harmonic oscillators, nQ is the number of meshes over coordinate,
    ! nP is number of meshes overr generalized momentum
    integer, parameter::  nQ = 90,  nP = 90

    ! Loop variables
    integer           ::  iQ, iP, iStep

    ! Mesh size
    real, parameter   :: DeltaQ = 24.0/nQ, DeltaP = 24.0/nP
    real :: VDF_G(-1:nQ+2, -1:nP+2),VDFInitial_C(nQ,nP),PlotVar_VC(2,nQ,nP)
    real :: Volume_G(0:nQ+1, 0:nP+1)
    real :: Hamiltonian_N(-1:nQ+1, -1:nP+1)
    real :: Energy_C(nQ, nP)
    real :: Time, Dt, Source_C(nQ,nP)
    real :: NormL2Init, NormL2, EnergyInit, Energy, qNode, pNode, Q, P
    real, parameter :: CFL =0.5
    real, parameter:: pWidth = 10.0, qWidth = 2.0
    integer, parameter:: nPower = 4
    logical, parameter:: IsSmooth = .true.
    !--------------------------------------------------------------------------
    ! Control volume, for a uniform rectangular grid
    Volume_G = DeltaQ*DeltaP
    ! Hamiltonian at the nodes
    do iP = -1, nP+1
       pNode = DeltaP*(iP - nP/2)
       do iQ = -1, nQ+1
          qNode = DeltaQ*(iQ - nQ/2)
          Hamiltonian_N(iQ,iP) = 0.5*(qNode**2 + pNode**2)
       end do
    end do
    ! Energy at the cell centers
    do iP = 1, nP
       P = DeltaP*(iP - nP/2 - 0.50)
       do iQ = 1, nQ
          Q = DeltaQ*(iQ - nQ/2  - 0.50)
          Energy_C(iQ,iP) = 0.5*(Q**2 + P**2)
       end do
    end do
    ! Initial distribution function
    VDFInitial_C = 0.0
    do iP = 1, nP; do iQ = 1, nQ
       P = DeltaP*(iP - nP/2 - 0.5)
       Q = DeltaQ*(iQ - nQ/2 - 0.5)
       if(abs(Q) < qWidth .and. abs(P) < pWidth)then
          if(IsSmooth)then
             VDFInitial_C(iQ,iP) = cos(0.5*cPi*Q/qWidth)**nPower &
                  *                cos(0.5*cPi*P/pWidth)**nPower
          else
             VDFInitial_C(iQ,iP) = 1.0
          endif
       end if
    end do; end do

    VDF_G = 0.0; VDF_G(1:nQ, 1:nP) = VDFInitial_C
    Source_C = 0.0
    ! Initial NormL2 and Energy
    NormL2Init = sum(VDFInitial_C**2)
    EnergyInit = sum(VDFInitial_C*Energy_C)
    ! Compiutation
    Time = 0.0; iStep = 0
    do
       call explicit(nQ, nP, VDF_G, Volume_G, Source_C, &
            Hamiltonian_N,   &
            CFLIn=CFL, DtOut = Dt)
       iStep = iStep +1
       if(Time + Dt >= tFinal)then
          call explicit(nQ, nP, VDF_G, Volume_G, Source_C, &
               Hamiltonian_N,   &
               DtIn = tFinal - Time)
          VDF_G(1:nQ, 1:nP) = VDF_G(1:nQ, 1:nP) + Source_C
          NormL2 = sum( (VDF_G(1:nQ,1:nP) - VDFInitial_C)**2 )
          Energy = sum(VDF_G(1:nQ,1:nP)*Energy_C)
          write(*,*)tFinal, &
               NormL2/NormL2Init,&
               Energy/EnergyInit - 1.0
          EXIT
       else
          Time = Time + Dt
          VDF_G(1:nQ, 1:nP) = VDF_G(1:nQ, 1:nP) + Source_C
          NormL2 = sum( (VDF_G(1:nQ,1:nP) - VDFInitial_C)**2 )
          Energy = sum(VDF_G(1:nQ,1:nP)*Energy_C)
          write(*,*)Time, &
               NormL2/NormL2Init,&
               Energy/EnergyInit - 1.0
       end if
    end do
    PlotVar_VC(1,:,:) = VDF_G(1:nQ,1:nP)
    PlotVar_VC(2,:,:) = VDFInitial_C
    call save_plot_file(NameFile='test_energy.out', &
         TypeFileIn='real8', TimeIn=tFinal, nStepIn = iStep, &
         NameVarIn='Q    P VDF VDFIni ErrorL2 EnergyDefect', &
         ParamIn_I=[NormL2/NormL2Init, Energy/EnergyInit - 1.0],&
         CoordMinIn_D=[-12.0 + 0.5*DeltaQ, -12.0 + 0.5*DeltaP], &
         CoordMaxIn_D=[12.0 - 0.5*DeltaQ, 12.0 - 0.5*DeltaP], &
         VarIn_VII = PlotVar_VC)

  end subroutine test_energy_conservation
  !============================================================================
  subroutine test_in_action_angle(tFinal)

    real, intent(in) :: tFinal

    ! Misc:
    ! Harmonic oscillators, nQ is the number of meshes over coordinate,
    ! nP is number of meshes overr generalized momentum
    integer, parameter::  nJ = 300,  nPhi = 360

    ! Loop variables
    integer           ::  iJ, iPhi, iStep

    ! Mesh size
    ! Momentum max and min
    real, parameter   :: JMax = 12.0, JMin = 0.01

    ! Mesh size in \phi
    real, parameter   :: DeltaPhi = cTwoPi/nPhi
    real :: MomentumRatio, MomentumMin, MomentumMax
    real :: VDF_G(-1:nJ+2, -1:nPhi+2), Volume_G(0:nJ+1, 0:nPhi+1)
    real :: VDFInitial_C(nJ,nPhi), PlotVar_VC(2,nJ,nPhi)

    real :: Hamiltonian_N(-1:nJ+1, -1:nPhi+1)
    real :: LogMomentum_I(0:nJ+1), Momentum2_I(-1:nJ+1)

    real :: Energy_C(nJ, nPhi)
    real :: Time, Dt, Source_C(nJ,nPhi)
    real :: NormL2Init, NormL2, EnergyInit, Energy, J, Phi, Q, P
    !--------------------------------------------------------------------------
    MomentumRatio = exp(log(JMax/JMin)/nJ)
    MomentumMax =  JMin/MomentumRatio
    Momentum2_I(-1) = MomentumMax**2
    Hamiltonian_N(-1,:) = 0.50*Momentum2_I(-1)
    do iJ = 0, nJ+1
       MomentumMin = MomentumMax
       MomentumMax = MomentumMin*MomentumRatio
       Momentum2_I( iJ) = MomentumMax**2
       Hamiltonian_N(iJ,:) = 0.50*Momentum2_I(iJ)
       Volume_G(iJ,:) = 0.5*DeltaPhi*&
            (Momentum2_I(iJ) - Momentum2_I(iJ-1))
       LogMomentum_I(iJ)   = 0.50*log10(MomentumMin*MomentumMax)
    end do
    ! Energy at the cell centers,  and initial distribution
    VDFInitial_C = 0.0
    do iPhi = 1, nPhi
       Phi = DeltaPhi*(iPhi - 0.50)
       do iJ = 1, nJ
          J = (Momentum2_I(iJ)*Momentum2_I(iJ-1))**0.250
          Energy_C(iJ,iPhi) = 0.5*J**2
          Q = J*cos(Phi)
          P = J*sin(Phi)
          ! Initial distribution function
          if(Q < -1.0.or.Q > 1.0.or.P < -10.0.or.P > 10.0)CYCLE
          VDFInitial_C(iJ,iPhi) = 1.0
       end do
    end do
    VDF_G = 0.0; VDF_G(1:nJ, 1:nPhi) = VDFInitial_C
    ! Periodic boundary conditions
    VDF_G(1:nJ,-1:0 )   = VDF_G(1:nJ, nPhi-1:nPhi)
    VDF_G(1:nJ, nPhi+1:nPhi+2 ) = VDF_G(1:nJ, 1:2)

    Source_C = 0.0
    ! Initial NormL2 and Energy
    NormL2Init = sum(VDFInitial_C**2*Volume_G(1:nJ,1:nPhi))
    EnergyInit = sum(VDFInitial_C*Energy_C*Volume_G(1:nJ,1:nPhi))
    ! Compiutation
    Time = 0.0; iStep = 0
    ! write(*,*)'Time NormL2/NormL2Init Energy/EnergyInit-1'
    do
       call explicit(nJ, nPhi, VDF_G, Volume_G, Source_C, &
            Hamiltonian_N,   &
            CFLIn=0.99, DtOut = Dt)
       iStep = iStep +1
       if(Time + Dt >= tFinal)then
          call explicit(nJ, nPhi, VDF_G, Volume_G, Source_C, &
               Hamiltonian_N,   &
               DtIn = tFinal - Time)
          VDF_G(1:nJ, 1:nPhi) = VDF_G(1:nJ, 1:nPhi) + Source_C
          NormL2 = sum( (VDF_G(1:nJ,1:nPhi) - VDFInitial_C)**2 &
               *Volume_G(1:nJ,1:nPhi))
          Energy = sum(VDF_G(1:nJ,1:nPhi)*Energy_C &
               *Volume_G(1:nJ,1:nPhi))
          write(*,*)tFinal, &
               NormL2/NormL2Init,&
               Energy/EnergyInit - 1.0
          EXIT
       else
          Time = Time + Dt
          VDF_G(1:nJ, 1:nPhi) = VDF_G(1:nJ, 1:nPhi) + Source_C
          ! Periodic boundary conditions
          VDF_G(1:nJ,-1:0 )   = VDF_G(1:nJ, nPhi-1:nPhi)
          VDF_G(1:nJ, nPhi+1:nPhi+2 ) = VDF_G(1:nJ, 1:2)
          NormL2 = sum( (VDF_G(1:nJ,1:nPhi) - VDFInitial_C)**2 &
               *Volume_G(1:nJ,1:nPhi))
          Energy = sum(VDF_G(1:nJ,1:nPhi)*Energy_C &
               *Volume_G(1:nJ,1:nPhi))
          write(*,*)Time, &
               NormL2/NormL2Init,&
               Energy/EnergyInit - 1.0
       end if
    end do
    PlotVar_VC(1,:,:) = VDF_G(1:nJ,1:nPhi)
    PlotVar_VC(2,:,:) = VDFInitial_C
    call save_plot_file(NameFile='test_action_angle.out', &
         TypeFileIn='real8', TimeIn=tFinal, nStepIn = iStep, &
         NameVarIn='Q    P VDF VDFIni ErrorL2 EnergyDefect', &
         CoordMinIn_D=[10.0**LogMomentum_I(1 ), 0.50*DeltaPhi],&
         CoordMaxIn_D=[10.0**LogMomentum_I(nJ), cTwoPi - 0.50*DeltaPhi], &
         ParamIn_I=[NormL2/NormL2Init, Energy/EnergyInit - 1.0],&
         Coord1In_I = 10.0**LogMomentum_I(1:nJ), &
         VarIn_VII = PlotVar_VC)
  end subroutine test_in_action_angle
  !============================================================================
  subroutine test_dsa_poisson

    use ModDiffusion

    real,    parameter :: pMax = 100, pMin = 1
    real,    parameter :: tFinal = 6000.00
    real,    parameter :: DtTrial = 1.0
    real ::  DtInv = 1/DtTrial
    real ::  MomentumRatio, MomentumMin, MomentumMax
    real ::  VDF_G(-1:nX+2, -1:nP+2)
    real ::  Volume_G(0:nX+1, 0:nP+1), VolumeNew_G(0:nX+1, 0:nP+1)
    real ::  VolumeX_I(0:nX+1), VolumeNewX_I(0:nX+1)
    real ::  VolumeP_I(0:nP+1)
    real ::  dVolumeDt_G(0:nX+1, 0:nP+1), dVolumeXDt_G(0:nX+1)
    real ::  DOuter_I(nX) = 1.0, dInner_I(nX) = 100.0 ! 200.0 For Fig5
    real ::  dHamiltonian02_FY(0:nX+1, -1:nP+1)
    real ::  LogMomentum_I(0:nP+1)     ! Cell centered, for plots
    real ::  Momentum3_I(-1:nP+1)
    real ::  Time = 0.0, Dt, DtNext, DtAdjust, Source_C(nX,nP)
    real ::  Coord_I(-1:nX+2)  ! Time-dependent coordinate of a mesh
    real ::  Dist_I(-1:nX+1)   ! Distance from mesh i to mesh i+1

    ! Loop variables
    integer :: iX, iP, iStep
    !--------------------------------------------------------------------------
    Source_C(nX,nP) = 0.0

    ! Logarithmic grid in momentum. pMin is the value at the left face
    ! of the first physical cell, pMax is the value at the right face
    ! off the last cell. The momentum ratio at the faces of each cell is:
    MomentumRatio   = exp(log(pMax/pMin)/nP)

    ! The momentum at the left fface of 0th ghost cell is:
    MomentumMax     = pMin/MomentumRatio

    ! Calculate the generalized variable p^3/3 at each face, starting from
    ! the left fface of 0th cell
    Momentum3_I(-1) = MomentumMax**3/3

    ! For cells from 0 to nP+1 calculate face values and volume factor:
    do iP = 0, nP + 1
       MomentumMin = MomentumMax
       MomentumMax = MomentumMin*MomentumRatio
       Momentum3_I(iP) = MomentumMax**3/3
       VolumeP_I(iP)   = Momentum3_I(iP) - Momentum3_I(iP-1)
       ! Only for visualization log10 of the cell centered momentum
       LogMomentum_I(iP) = 0.50*log10(MomentumMin*MomentumMax)
    end do
    VDF_G = 1.0e-8; VDF_G(:,1) = 1/VolumeP_I(1)

    Time = 0.0; iStep = 0; DtNext = DtTrial
    call update_coords(Time)
    ! Figure 4, left panel
    ! call update_coords(6000.0)
    ! do iX = 5982,6005
    !   write(*,*)Coord_I(iX),1/VolumeNewX_I(iX)
    ! end do
    ! stop
    do
       Volume_G = VolumeNew_G
       VolumeX_I= VolumeNewX_I
       Dt = min(DtNext, tFinal - Time); DtInv = 1/Dt
       call update_coords(Time + Dt)
       dVolumeDt_G  = DtInv*(VolumeNew_G  - Volume_G)
       dVolumeXDt_G = DtInv*(VolumeNewX_I - VolumeX_I)
       do iP = -1, nP + 1
          dHamiltonian02_FY(:,iP) = - Momentum3_I(iP)*dVolumeXDt_G
       end do
       call explicit(nX, nP, VDF_G, Volume_G, Source_C, &
            dHamiltonian02_FY=dHamiltonian02_FY,        &
            dVolumeDt_G = dVolumeDt_G,                  &
            DtIn=Dt,           & ! Input time step
            CFLIn=0.98,        & ! Input CFL to calculate next time step
            DtOut=DtAdjust,    & ! Actual time step, which may be reduced
            DtRecommend=DtNext)  ! Calculation of next time step
       Dt = DtAdjust
       ! Correct final volumes, if Dt had been reduced
       VolumeNew_G  = Volume_G  + Dt*dVolumeDt_G
       VolumeNewX_I = VolumeX_I + Dt*dVolumeXDt_G
       iStep = iStep +1
       VDF_G(1:nX, 1:nP) = VDF_G(1:nX, 1:nP) + Source_C
       do iP =1, nP
          call advance_diffusion1(Dt,nX,Dist_I(1:nX),VDF_G(1:nX,&
               iP),DOuter_I(1:nX),DInner_I(1:nX))
       end do
       Time = Time + Dt
       if(Time> tFinal - 1.0e-8*DtNext)EXIT
       VDF_G(1:nX,-1:0 ) = 1.0e-8
       VDF_G(1:nX,nP+1:nP+2) = 1.0e-8
       VDF_G( 0,      :) = VDF_G(1,       :)
       VDF_G(-1,      :) = VDF_G(1,       :)
       VDF_G(nX+1,    :) = VDF_G(nX,      :)
       VDF_G(nX+2,    :) = VDF_G(nX,      :)
    end do
    ! do iP =1, nP
    !   write(*,*)LogMomentum_I(iP),alog10(VDF_G(5000,iP)*VolumeNew_G(5000,iP))
    ! end do
    call save_plot_file(NameFile='test_dsa_poisson.out', &
         TypeFileIn='ascii', TimeIn=tFinal, nStepIn = iStep, &
         NameVarIn='LogMomentum VDF'  ,                  &
         CoordMinIn_D=[ LogMomentum_I(1 ) ],             &
         CoordMaxIn_D=[ LogMomentum_I(nP) ],             &
         StringFormatIn = '(2F16.9)',                    &
         Coord1In_I = LogMomentum_I(1:nP),               &
         VarIn_I = alog10(VDF_G(5000,1:nP)))
    VDFOld_G = VDF_G
  contains
    !==========================================================================
    subroutine update_coords(Time)
      real, intent(in):: Time

      ! Effective width of the shock wave, in cell size
      integer, parameter:: nWidth = 1

      ! terms in the functions, describing the
      ! lagrangian point motion
      real, parameter:: cQuad = 0.375/nWidth, cConst = 0.375*nWidth
      ! Loop variables
      integer :: iX, iP
      real    :: CoordLagr
      ! Shock wave with the copmtression ratio of 4, with width ~1
      !------------------------------------------------------------------------
      do iX = -1, nX +2
         CoordLagr = iX - 0.50 -Time
         if(CoordLagr>=0.0)then
            Coord_I(iX) = CoordLagr + Time
         elseif(CoordLagr>=-real(nWidth))then
            Coord_I(iX) = CoordLagr + Time + cQuad*CoordLagr**2
         else
            Coord_I(iX) = 0.250*CoordLagr -cConst + Time
         end if
      end do
      Dist_I(-1:nX+1) = Coord_I(0:nX+2) - Coord_I(-1:nX+1)
      VolumeNewX_I(0:nX+1) = 0.50*(Dist_I(-1:nX) + Dist_I(0:nX+1))
      do iP = 0, nP+1
         VolumeNew_G(0:nX+1,iP) = VolumeP_I(iP)*VolumeNewX_I(0:nX+1)
      end do
    end subroutine update_coords
    !==========================================================================
  end subroutine test_dsa_poisson
  !============================================================================
  subroutine test_dsa_sa_mhd

    use ModDiffusion

    real,    parameter :: pMax = 100, pMin = 1
    real,    parameter :: tFinal = 6000.00
    real ::  MomentumRatio, MomentumMin, MomentumMax
    real ::  VDF_G(-1:nX+2, -1:nP+2)
    real ::  Hamiltonian_N(-1:nX+1, -1:nP+1)
    real ::  Volume_G(0:nX+1, 0:nP+1)
    real ::  VolumeX_I(-1:nX+2)
    real ::  VolumeP_I(0:nP+1)
    real ::  DOuter_I(nX) = 1.0, dInner_I(nX) = 200.0
    real ::  LogMomentum_I(0:nP+1)     ! Cell centered, for plots
    real ::  Momentum3_I(-1:nP+1)
    real ::  Time = 0.0, Dt, Source_C(nX,nP), CFL
    real ::  Coord_I(-1:nX+2)  ! Time-dependent coordinate of a mesh
    real ::  Dist_I(-1:nX+1)   ! Distance from mesh i to mesh i+1

    ! Loop variables
    integer :: iX, iP, iStep
    !--------------------------------------------------------------------------
    Source_C(nX,nP) = 0.0
    ! Logarithmic grid in momentum. pMin is the value at the left face
    ! of the first physical cell, pMax is the value at the right face
    ! off the last cell. The momentum ratio at the faces of each cell is:
    MomentumRatio   = exp(log(pMax/pMin)/nP)
    ! The momentum at the left fface of 0th ghost cell is:
    MomentumMax     = pMin/MomentumRatio
    ! Calculate the generalized variable p^3/3 at each face, starting from
    ! the left fface of 0th cell
    Momentum3_I(-1) = MomentumMax**3/3
    ! For cells from 0 to nP+1 calculate face values and volume factor:
    do iP = 0, nP + 1
       MomentumMin = MomentumMax
       MomentumMax = MomentumMin*MomentumRatio
       Momentum3_I(iP) = MomentumMax**3/3
       VolumeP_I(iP)   = Momentum3_I(iP) - Momentum3_I(iP-1)
       ! Only for visualization log10 of the cell centered momentum
       LogMomentum_I(iP) = 0.50*log10(MomentumMin*MomentumMax)
    end do
    VDF_G = VDFOld_G

    Time = 0.0; iStep = 0
    call update_coords(6000.0)
    ! Figure 5, left panel
    ! do iX = 5982,6005
    !   write(*,*)0.50*(Coord_I(iX)+Coord_I(iX+1)),&
    !        1/(Coord_I(iX+1)-Coord_I(iX))
    ! end do
    ! stop
    do
       VDFOld_G = VDF_G
       call explicit(nX, nP, VDF_G, Volume_G, Source_C, &
            Hamiltonian12_N=Hamiltonian_N,              &
            CFLIn=0.98, DtOut=Dt)
       iStep = iStep +1
       if(Time + Dt>= tFinal)then
          Dt = tFinal - Time
          call explicit(nX, nP, VDF_G, Volume_G, Source_C, &
               Hamiltonian12_N=Hamiltonian_N,              &
               DtIn=Dt)
          VDF_G(1:nX, 1:nP) = VDF_G(1:nX, 1:nP) + Source_C
          do iP =1, nP
             call advance_diffusion1(Dt,nX,Dist_I(1:nX),VDF_G(1:nX,&
                  iP),DOuter_I(1:nX),DInner_I(1:nX))
          end do
          Time = tFinal
          EXIT
       end if
       VDF_G(1:nX, 1:nP) = VDF_G(1:nX, 1:nP) + Source_C
       do iP =1, nP
          call advance_diffusion1(Dt,nX,Dist_I(1:nX),VDF_G(1:nX,&
               iP),DOuter_I(1:nX),DInner_I(1:nX))
       end do
       Time = Time + Dt
       VDF_G(1:nX,-1:0 ) = 1.0e-8
       VDF_G(1:nX,nP+1:nP+2) = 1.0e-8
       VDF_G( 0,      :) = VDF_G(1,       :)
       VDF_G(-1,      :) = VDF_G(1,       :)
       VDF_G(nX+1,    :) = max(VDF_G(nX,      :), 1.0e-8)
       VDF_G(nX+2,    :) = max(VDF_G(nX,      :), 1.0e-8)
       VDF_G(nX+1:nX+2,1) = 1/VolumeP_I(1)
    end do
    call save_plot_file(NameFile='test_dsa_sa_mhd.out', &
         TypeFileIn='ascii', TimeIn=tFinal, nStepIn = iStep, &
         NameVarIn='LogMomentum VDF'  ,                  &
         CoordMinIn_D=[ LogMomentum_I(1 ) ],             &
         CoordMaxIn_D=[ LogMomentum_I(nP) ],             &
         StringFormatIn = '(2F16.9)',                    &
         Coord1In_I = LogMomentum_I(1:nP),               &
         VarIn_I = alog10(VDF_G(5000,1:nP)))
  contains
    !==========================================================================
    subroutine update_coords(Time)

      real, intent(in):: Time

      ! Effective width of the shock wave, in cell size
      integer, parameter:: nWidth = 1

      ! terms in the functions, describing the
      ! lagrangian point motion
      real, parameter:: cQuad = 0.375/nWidth, cConst = 0.375*nWidth
      ! Loop variables
      integer :: iX, iP
      real    :: CoordLagr
      ! Shock wave with the copmtression ratio of 4, with width ~1
      !------------------------------------------------------------------------
      do iX = -1, nX +2
         CoordLagr = iX - 0.50 -Time
         if(CoordLagr>=0.0)then
            Coord_I(iX) = CoordLagr + Time
         elseif(CoordLagr>=-real(nWidth))then
            Coord_I(iX) = CoordLagr + Time + cQuad*CoordLagr**2
         else
            Coord_I(iX) = 0.250*CoordLagr -cConst + Time
         end if
      end do
      Dist_I(-1:nX+1) = Coord_I(0:nX+2) - Coord_I(-1:nX+1)
      VolumeX_I(0:nX+1) = 0.50*(Dist_I(-1:nX) + Dist_I(0:nX+1))
      VolumeX_I(-1)  = Dist_I(-1); VolumeX_I(nX+2) =  Dist_I(nX+1)
      do iP = 0, nP+1
         Volume_G(0:nX+1,iP) = VolumeP_I(iP)*VolumeX_I(0:nX+1)
      end do
      do iP = -1, nP+1
         Hamiltonian_N(-1:nX+1, iP)= -0.50*Momentum3_I(iP)*&
              (VolumeX_I(-1:nX+1) + VolumeX_I(0:nX+2))
      end do

    end subroutine update_coords
    !==========================================================================
  end subroutine test_dsa_sa_mhd
  !============================================================================
end module ModTestPoissonBracket
!==============================================================================
module ModHillVortex

  use ModPoissonBracket, ONLY: explicit
  use ModUtilities,      ONLY: CON_stop
  use ModPlotFile,       ONLY: save_plot_file
  use ModConst

  implicit none

  PRIVATE ! Except
  public :: test_hill_vortex

contains
  !============================================================================
  subroutine test_hill_vortex

    ! Grid in RSph-Theta variables
    ! Theta goes from 0 to 360 because it is mirrored over the pole for plots
    integer, parameter:: nR = 54, nTheta = 2*nR
    real,    parameter:: Dr = 4.50/nR, Cfl = 0.9

    ! Streamlines:
    integer, parameter:: r_ = 2, Ur_ = 2, z_ = 1, Uz_  = 1

    ! Position and width of the smooth bump or sharp ellipsoid
    real, parameter:: xCenter = -3.0, yCenter = 0.0, WidthX = 3.0, WidthY = 6.0
    logical:: IsSmooth = .true.

    real:: x, y, z, r, Uz, Ur, Vel_VC(Uz_:Ur_,-500:500,-500:500), RSph

    ! Loop variables
    integer:: i, j, iR, iTheta, iStep, iPlot

    ! Mesh size in Theta
    real, parameter   :: dTheta = cTwoPi/nTheta
    real :: VDF_G(-1:nR+2,-1:nTheta), Volume_G(0:nR+1, 0:nTheta/2+1)
    real :: R_I(-1:nR+1), Theta_I(-1:nTheta+1), CosTheta_I(-1:nTheta/2+1)
    real :: Hamiltonian_N(-1:nR+1, -1:nTheta/2+1)
    real :: Time, Dt, Source_C(nR,nTheta/2), tFinal
    real :: PlotVar_VC(5,nR,nTheta), DfDt_C(nR,nTheta), DfDt2_C(nR,nTheta)
    real :: Error, Error2, URSph, UTheta,ErrorTVD
    logical :: DoExit
    !--------------------------------------------------------------------------
    do j = -500,500
       r = 0.01*j
       do i = -500,500
          z = 0.01*i
          RSph = sqrt(z**2 + r**2)
          if(RSph < 1.0)then
             ! Inside the vortex
             Uz =  1.50*(-1.0 + r**2 + RSph**2)
             Ur = -1.50*r*z
          else
             Uz = 1.0 - 1.0/RSph**3 + 1.50*r**2/RSph**5
             Ur = -1.50*r*z/RSph**5
          end if
          Vel_VC(Uz_:Ur_,i,j) = [Uz, Ur]
       end do
    end do
    call save_plot_file( &
         NameFile = 'streamlines.out',                &
         StringHeaderIn='Hill vortex velocity field', &
         NameVarIn='z r Uz Ur'  ,                     &
         CoordMinIn_D=[-5.0, -5.0],                   &
         CoordMaxIn_D=[ 5.0,  5.0 ],                  &
         StringFormatIn = '(4F16.4)',                 &
         VarIn_VII = Vel_VC)

    ! Initialize coords
    ! These are cell face coordinates:
    ! r_I(iR) is the radius at iR+1/2, Theta_I(iTheta) is at iTheta+1/2.
    do iR = -1, nR + 1
       R_I(iR) = iR*Dr + 0.5
    end do
    do iTheta = -1, nTheta + 1
       Theta_I(iTheta) = iTheta*dTheta
    end do
    do iTheta = -1, nTheta/2+1
       CosTheta_I(iTheta) = cos(Theta_I(iTheta))
    end do
    ! Overwrite ghost cell values
    CosTheta_I(nTheta/2+1) = CosTheta_I(nTheta/2-1)
    CosTheta_I(-1) = CosTheta_I(1)

    ! Calculate Hamiltonian
    do iTheta = -1, nTheta/2+1
       ! Calculate theta-dependent factor
       Hamiltonian_N(:, iTheta) = cTwoPi*(1.0 - CosTheta_I(iTheta)**2)
       do iR = -1, nR/9 ! Unit radius
          Hamiltonian_N(iR,iTheta) = Hamiltonian_N(iR,iTheta)*0.750*&
               R_I(iR)**2*(R_I(iR)**2 - 1.0)
       end do
       do iR = nR/9+1, nR+1
          Hamiltonian_N(iR,iTheta) = Hamiltonian_N(iR,iTheta)*0.50*&
               R_I(iR)**2*(1.0 - 1.0/R_I(iR)**3)
       end do
    end do
    ! write(*,*)'Maxval(Hamiltonian)=',maxval(Hamiltonian_N)

    ! Calculate volume
    ! angle-dependent factor
    do iTheta = 1, nTheta/2
       Volume_G(0:nR+1,iTheta) = (CosTheta_I(iTheta-1) - CosTheta_I(iTheta))*&
            (R_I(0:nR+1)**3 - R_I(-1:nR)**3)*cTwoPi/3
    end do
    Volume_G(0:nR+1,0) = (CosTheta_I(0) - CosTheta_I(-1)) &
         *(R_I(0:nR+1)**3 - R_I(-1:nR)**3)*cTwoPi/3
    Volume_G(0:nR+1,nTheta/2+1) = &
         (CosTheta_I(nTheta/2+1) - CosTheta_I(nTheta/2)) &
         *(R_I(0:nR+1)**3 - R_I(-1:nR)**3)*cTwoPi/3
    ! write(*,*)'Minval(Volume)=',minval(Volume_G)
    ! Periodic boundary conditions
    Source_C = 0.0

    ! Initial conditions:
    VDF_G = 0.0; DfDt_C = 0.0; DfDt2_C = 0.0
    do iTheta = 1, nTheta; do iR = 1, nR
       x = cos(Theta_I(iTheta) - dTheta/2)*(R_I(iR) - dR/2) - xCenter
       y = sin(Theta_I(iTheta) - dTheta/2)*(R_I(iR) - dR/2) - yCenter
       if(IsSmooth)then
          ! inside a WidthX*WidthY rectangle centered on xCenter,yCenter
          ! cos^2(kx*(x-xCenter))*cos^2(ky*(y-yCenter))
          if(abs(x) < WidthX/2 .and. abs(y) < WidthY/2) then
             VDF_G(iR,iTheta) = cos(cPi*x/WidthX)**4 * cos(cPi*y/WidthY)**4
             z = x + xCenter; r = y + yCenter
             RSph = sqrt(z**2 + r**2)
             Uz = 1.0 - 1.0/RSph**3 + 1.50*r**2/RSph**5
             Ur = -1.50*r*z/RSph**5
             ! Anaylitical (df/dt=-u.grad f), to compare with numeric source
             DfDt_C(iR,iTheta) =cos(cPi*x/WidthX)**3 *cos(cPi*y/WidthY)**3&
                  *(Uz*cPi/WidthX*sin(cPi*x/WidthX)*cos(cPi*y/WidthY) + &
                  Ur*cPi/WidthY*cos(cPi*x/WidthX)*sin(cPi*y/WidthY))*4.0
          end if
       else
          if( (x/WidthX)**2 + (y/WidthY)**2 < 0.25) &
               VDF_G(iR,iTheta) = 1.0
       end if
    end do; end do
    if(IsSmooth)then
       do iTheta = 1, nTheta; do iR = 1, nR
          x = cos(Theta_I(iTheta) - dTheta/2)*(R_I(iR) - dR/2) - xCenter
          y = sin(Theta_I(iTheta) - dTheta/2)*(R_I(iR) - dR/2) - yCenter
          if(abs(x) < WidthX/2 .and. abs(y) < WidthY/2) then
             VDF_G(iR,iTheta) = cos(cPi*x/WidthX)**4 * cos(cPi*y/WidthY)**4
             z = x + xCenter; r = y + yCenter
             RSph = sqrt(z**2 + r**2)
             Uz = 1.0 - 1.0/RSph**3 + 1.50*r**2/RSph**5
             Ur = -1.50*r*z/RSph**5
             URsph = (Uz*z + Ur*r)/RSph
             UTheta = (-Uz*r + Ur*z)/RSph
             ! By differentiating  df/dt=-u.grad f, the second order derivative
             ! may be evaluated: d2f/dt2=-u.grad(df/dt), so that analytical
             ! df/dt may be differentiated numerically:
             DfDt2_C(iR,iTheta) =-0.50*(&
                  URsph*(DfDt_C(iR+1,iTheta) - DfDt_C(iR-1,iTheta))/Dr + &
                  UTheta*(DfDt_C(iR,iTheta+1) - DfDt_C(iR,iTheta-1))/&
                  (RSph*dTheta) )
          end if
       end do; end do
    end if

    Time = 0.0; iStep = 0

    ! Save initial conditions
    call save_plot_file('hill.outs', 'rewind', 'real8', &
         'Hill vortex', iStep, Time, &
         NameVarIn='r Theta Rho'  , &
         CoordMinIn_D=[0.5 + 0.5*Dr, 180.0/nTheta],&
         CoordMaxIn_D=[5.0 - 0.5*Dr, 360.0 - 180.0/nTheta],&
         VarIn_II = VDF_G(1:nR,1:nTheta) )

    call save_plot_file('hillcut.outs', 'rewind', 'real8', &
         'Hill vortex', iStep, Time, &
         NameVarIn='r Theta Rho'  , &
         CoordMinIn_D=[0.5 + 0.5*Dr, 180 + 180.0/nTheta],&
         CoordMaxIn_D=[5.0 - 0.5*Dr, 360 - 180.0/nTheta],&
         VarIn_II = VDF_G(1:nR,nTheta/2+1:nTheta) )
    if(IsSmooth)then
       ! Calculate source for small time step, to compare with
       ! analytical DfDt
       call explicit(nR, nTheta/2, VDF_G(-1:nR+2,-1:nTheta/2+2), Volume_G,&
            Source_C, Hamiltonian_N,   &
            CFLIn=1.0e-6, DtOut = Dt,ErrorTVD=ErrorTVD)
       PlotVar_VC(1,:,:) = VDF_G(1:nR,1:nTheta)
       PlotVar_VC(2,:,:) = DfDt_C(1:nR,1:nTheta)
       PlotVar_VC(3,:,1:nTheta/2) = Source_C/Dt
       Error = sum(abs(PlotVar_VC(2,:,1:nTheta/2) - &
            PlotVar_VC(3,:,1:nTheta/2))*&
            Volume_G(1:nR,1:nTheta/2))/(4/3.0*cPi*5.0**3)
       call explicit(nR, nTheta/2, VDF_G(-1:nR+2,-1:nTheta/2+2), Volume_G,&
            Source_C, Hamiltonian_N,   &
            CFLIn=CFL, DtOut = Dt)
       PlotVar_VC(4,:,:) = 0.50*Dt*DfDt2_C(1:nR,1:nTheta)
       PlotVar_VC(5,:,1:nTheta/2) = Source_C/Dt
       Error2 = sum(abs(&
            PlotVar_VC(2,:,1:nTheta/2) + PlotVar_VC(4,:,1:nTheta/2) - &
            PlotVar_VC(5,:,1:nTheta/2))*&
            Volume_G(1:nR,1:nTheta/2))/(4/3.0*cPi*5.0**3)
       do iTheta = nTheta/2 + 1, nTheta
          ! Symmetric prolongation for visualization
          PlotVar_VC(3,:, iTheta) = PlotVar_VC(3,:, 1+nTheta-iTheta)
          PlotVar_VC(5,:, iTheta) = PlotVar_VC(5,:, 1+nTheta-iTheta)
       end do
       call save_plot_file('InitialSource.out', 'rewind','real8', &
            'Hill vortex', iStep, Time, &
            NameVarIn=&
            'z r Rho DfDt Source DfDt2 Source2 Error ErrorTVD  Error2', &
            CoordMinIn_D=[0.5 + 0.5*Dr, 180.0/nTheta],&
            CoordMaxIn_D=[5.0 - 0.5*Dr, 360.0 - 180.0/nTheta],&
            ParamIn_I=[Error, ErrorTVD, Error2],&
            VarIn_VII = PlotVar_VC )
       write(*,*)'Error=', Error, ' ErrorTVD=',ErrorTVD,' Error2=', Error2
       !!! stop
    end if
    ! Computation
    do iPlot = 0, 99
       tFinal = (iPlot + 1)*0.1
       DoExit = .false.
       PLOT:do
          call explicit(nR, nTheta/2, VDF_G(-1:nR+2,-1:nTheta/2+2), Volume_G,&
               Source_C, Hamiltonian_N,   &
               CFLIn=Cfl, DtOut = Dt)
          iStep = iStep +1
          if(Time + Dt >= tFinal)then
             call explicit(nR, nTheta/2, VDF_G(-1:nR+2,-1:nTheta/2+2), &
                  Volume_G, Source_C, Hamiltonian_N, DtIn = tFinal - Time)
             VDF_G(1:nR, 1:nTheta/2) = VDF_G(1:nR, 1:nTheta/2) + Source_C
             Time = tFinal
             DoExit = .true.
          else
             Time = Time + Dt
             VDF_G(1:nR, 1:nTheta/2) = VDF_G(1:nR, 1:nTheta/2) + Source_C
          end if
          do iTheta = nTheta/2 + 1, nTheta
             ! Symmetric prolongation
             VDF_G(1:nR, iTheta) = VDF_G(1:nR, 1+nTheta-iTheta)
          end do
          ! Periodic prolongation
          VDF_G(1:nR,-1:0) = VDF_G(1:nR, nTheta-1:nTheta)
          write(*,*)'Time=',Time
          if(DoExit) EXIT PLOT
       end do PLOT
       call save_plot_file('hill.outs', 'append', 'real8', &
            'Hill vortex', iStep, tFinal, &
            NameVarIn='r Theta Rho'  , &
            CoordMinIn_D=[0.50 + 0.50*Dr, 180.0/nTheta],&
            CoordMaxIn_D=[5.0 - 0.50*Dr, 360.0 - 180.0/nTheta],&
            VarIn_II = VDF_G(1:nR,1:nTheta) )
       call save_plot_file('hillcut.outs', 'append', 'real8', &
            'Hill vortex', iStep, tFinal, &
            NameVarIn='r Theta Rho'  , &
            CoordMinIn_D=[0.5 + 0.5*Dr, 180 + 180.0/nTheta],&
            CoordMaxIn_D=[5.0 - 0.5*Dr, 360 - 180.0/nTheta],&
            VarIn_II = VDF_G(1:nR,nTheta/2+1:nTheta) )
    end do
  end subroutine test_hill_vortex
  !============================================================================
end module ModHillVortex
!==============================================================================
module ModStochastic
  use ModPoissonBracket, ONLY: explicit
  use ModUtilities,      ONLY: CON_stop
  use ModPlotFile,       ONLY: save_plot_file
  use ModConst
  implicit none
  PRIVATE ! Except
  public :: test_stochastic
  ! Streamlines:
  integer, parameter :: nPointPer2Pi = 180
  integer, parameter :: nJ = (4*nPointPer2Pi)/2, nTheta=nPointPer2Pi
  real, parameter    :: Delta = cTwoPi/nPointPer2Pi
  ! Loop variables
  integer           ::  iJ, iTheta, iStep, iPlot
  ! Mesh size in \Theta
  real :: VDF_G(  -1:nTheta+2, -1-nJ:nJ+2)
  real :: Volume_G(0:nTheta+1,   -nJ:nJ+1)
  real :: Theta_I(-1:nTheta+1), Action_I(-1-nJ:nJ+1), Action
  real :: HamiltonianFree_N(-1:nTheta+1, -1-nJ:nJ+1)
  real :: HamiltonianPush_N(-1:nTheta+1, -1-nJ:nJ+1)
  real :: Time, Dt, Source_C(nTheta,1-nJ:nJ), tFinal
  logical :: DoAgain
  character(LEN=11)::NameFile
  !----------------------------------------------------------------------------
contains
  !============================================================================
  subroutine test_stochastic(KChirikov)
    real, intent(in) :: KChirikov ! Left hand side in the Chhirikov criterion
    !--------------------------------------------------------------------------
    ! Uniform Delta*Delta grid:
    Volume_G = Delta**2
    do iJ = -1-nJ, nJ+1
       Action_I(iJ)  = iJ*Delta
    end do
    do iTheta = -1, nTheta+1
       Theta_I(iTheta) = iTheta*Delta
    end do
    ! Free rotator:
    do iJ = -1-nJ, nJ+1
       HamiltonianFree_N(-1:nTheta+1, iJ) = 0.50*Action_I(iJ)**2
    end do
    ! Sinusoidal perturbation
    do iTheta = -1, nTheta+1
       HamiltonianPush_N(iTheta, -1-nJ:nJ+1) = cos(Theta_I(iTheta))
    end do
    ! Account for:
    ! 1. Amplitude factor, KChirikov
    ! 2. Dirac delta-function, by dividing the amplitude by the short time
    !    interval and apply the perturbation only within this interval
    ! 3. Sum up both Hamiltonians:
    HamiltonianPush_N = (KChirikov/0.10)*HamiltonianPush_N + HamiltonianFree_N
    ! Initial distribution function: zero at high energies and in ghost cells:
    VDF_G(-1:nTheta+2,-1-nJ:nJ+2) =  exp(-12.50)
    do iJ  = 1, nJ
       Action = 0.50*(Action_I(iJ-1) + Action_I(iJ))
       if(Action < 1.0)then
          VDF_G(:,iJ) = exp(-12.50*Action**2)
       else
          EXIT
       end if
    end do
    do iJ  = 0, 1-nJ,-1
       Action = 0.50*(Action_I(iJ-1) + Action_I(iJ))
       if(abs(Action) < 1.0)then
          VDF_G(:,iJ) = exp(-12.50*Action**2)
       else
          EXIT
       end if
    end do
    Source_C = 0.0
    tFinal  = 0.0
    ! Compiutation
    Time = 0.0; iStep = 0
    call save_plot_file(NameFile='initial.out', &
         TypeFileIn='real8', TimeIn=tFinal, nStepIn = iStep, &
         NameVarIn='theta JT/2pi Log10VDF'  , &
         CoordMinIn_D=[         0.50*Delta, -2.0 + 0.50/nPointPer2Pi],&
         CoordMaxIn_D=[cTwoPi - 0.50*Delta,  2.0 - 0.50/nPointPer2Pi],&
         VarIn_II = log10(VDF_G(1:nTheta,1-nJ:nJ)) )
    do iPlot = 0, 999
       tFinal = tFinal + 0.1
       DoAgain = .true.
       do while(DoAgain)
          call explicit(nTheta,2*nJ, VDF_G, Volume_G,&
               Source_C, HamiltonianPush_N,   &
               CFLIn=0.5, DtOut = Dt,IsPeriodicIn_D=[.true.,.false.])
          iStep = iStep +1
          if(Time + Dt >= tFinal)then
             call explicit(nTheta,2*nJ,VDF_G,Volume_G,&
                  Source_C, HamiltonianPush_N,   &
                  DtIn = tFinal - Time,IsPeriodicIn_D=[.true.,.false.])
             VDF_G(1:nTheta,1-nJ:nJ) = VDF_G(1:nTheta,1-nJ:nJ) + Source_C
             Time = tFinal
             DoAgain = .false.
          else
             Time = Time + Dt
             VDF_G(1:nTheta,1-nJ:nJ) = VDF_G(1:nTheta,1-nJ:nJ) + Source_C
          end if
          ! Periodic prolongation
          VDF_G(-1:0,:) = VDF_G(nTheta-1:nTheta,:)
          VDF_G(nTheta+1:nTheta+2,:) = VDF_G(1:2,:)
          write(*,*)'Time=',Time
       end do
       tFinal = tFinal + 0.9
       DoAgain = .true.
       do while(DoAgain)
          call explicit(nTheta,2*nJ, VDF_G, Volume_G,&
               Source_C, HamiltonianFree_N,   &
               CFLIn=0.99, DtOut = Dt,IsPeriodicIn_D=[.true.,.false.])
          iStep = iStep +1
          if(Time + Dt >= tFinal)then
             call explicit(nTheta,2*nJ,VDF_G,Volume_G,&
                  Source_C, HamiltonianFree_N,   &
                  DtIn = tFinal - Time,IsPeriodicIn_D=[.true.,.false.])
             VDF_G(1:nTheta,1-nJ:nJ) = VDF_G(1:nTheta,1-nJ:nJ) + Source_C
             Time = tFinal
             DoAgain = .false.
          else
             Time = Time + Dt
             VDF_G(1:nTheta,1-nJ:nJ) = VDF_G(1:nTheta,1-nJ:nJ) + Source_C
          end if
          ! Periodic prolongation
          VDF_G(-1:0,:) = VDF_G(nTheta-1:nTheta,:)
          VDF_G(nTheta+1:nTheta+2,:) = VDF_G(1:2,:)
          write(*,*)'Time=',Time
       end do
       write(NameFile,'(a,i3.3,a)')'stoc',iPlot,'.out'
       call save_plot_file(NameFile=NameFile, &
            TypeFileIn='real8', TimeIn=tFinal, nStepIn = iStep, &
            NameVarIn='theta JT/2pi Log10VDF'  , &
            CoordMinIn_D=[         0.50*Delta, -2.0 + 0.50/nPointPer2Pi],&
            CoordMaxIn_D=[cTwoPi - 0.50*Delta,  2.0 - 0.50/nPointPer2Pi],&
            VarIn_II = log10(VDF_G(1:nTheta,1-nJ:nJ)) )
    end do
  end subroutine test_stochastic
  !============================================================================
end module ModStochastic
!==============================================================================
program test_program
  use ModTestPoissonBracket, ONLY: test_poisson_bracket, test_dsa_sa_mhd, &
       test_dsa_poisson, test_energy_conservation, test_in_action_angle,  &
       test_poisson_2d, test_poisson_2d_smooth
  use ModNumConst,           ONLY: cTwoPi
  use ModHillVortex, ONLY: test_hill_vortex
  use ModStochastic, ONLY: test_stochastic
  ! use ModTestPoissonBracketAndScatter, ONLY: test_scatter

  implicit none
  !----------------------------------------------------------------------------
  write(*,*)' Start poisson bracket tests'

  call test_poisson_bracket(cTwoPi)        ! nightly test1
  call test_dsa_poisson                    ! nightly test2
  ! call test_poisson_2d(cTwoPi)  ! Fig 1 (right panel)
  ! call test_poisson_2d_smooth(cTwoPi)
  ! call test_stochastic(1.2)
  ! call test_hill_vortex
  ! call test_energy_conservation(cTwoPi)
  ! call test_in_action_angle(cTwoPi)
  ! call test_dsa_sa_mhd ! for Fig5.Right Panel
  ! call test_multipoisson_bracket(50.0)

  ! Test for comparing two diffusions is long. It should be repeated twice with
  ! Diffuornot = 0 or 1

  ! call test_scatter(50.0)
  write(*,*)' Finish poisson bracket tests'
  
end program test_program
!==============================================================================
