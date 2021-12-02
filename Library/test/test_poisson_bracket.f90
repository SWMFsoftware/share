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
  ! This routine solves the diffusion equation:
  !         f_t-D_outer(D_inner*f_x)_x=0,
  ! with zero Neumann boundary condition. The solution is advanced in time
  ! using fully implicit scheme.
  subroutine advance_diffusion1(Dt,n,Dist_I,F_I,DOuter_I,DInner_I)
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
  ! This routine solves three-diagonal system of equations:
  !  ||m_1 u_1  0....        || ||w_1|| ||r_1||
  !  ||l_2 m_2 u_2...        || ||w_2|| ||r_2||
  !  || 0  l_3 m_3 u_3       ||.||w_3||=||r_3||
  !  ||...                   || ||...|| ||...||
  !  ||.............0 l_n m_n|| ||w_n|| ||r_n||
  ! From: Numerical Recipes, Chapter 2.6, p.40.
  subroutine tridiag(n, L_I, M_I, U_I, R_I, W_I)
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
    !   write(*,*)LogMomentum_I(iP),alog10(VDF_G(5000,iP))
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
program test_program
  use ModTestPoissonBracket, ONLY: test_poisson_bracket, test_dsa_sa_mhd, &
       test_dsa_poisson
  use ModNumConst,           ONLY: cTwoPi
  ! use ModTestPoissonBracketAndScatter, ONLY: test_scatter
  implicit none

  !----------------------------------------------------------------------------
  call test_poisson_bracket(cTwoPi)
  call test_dsa_poisson
  ! call test_dsa_sa_mhd ! for Fig5.Right Panel
  ! call test_multipoisson_bracket(50.0)

  ! Test for comparing two diffusions is long. It should be repeated twice with
  ! Diffuornot = 0 or 1

  ! call test_scatter(50.0)
end program test_program
!==============================================================================
