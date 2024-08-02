!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!
! For a contribution to the Liuville (kinetic) equation:
!
! df/dt + \sum_l{df/dq_l dH_l/dp_l - df/dp_l dH_l/dq_l}=0,
!
! in which f in the (unknown) distribution function,
! q_l and p_l are generalized coordinates,
! H_l are the (known!) Hamiltonian functions.
!
! See	arXiv:1910.12636 [physics.comp-ph]
!

module ModPoissonBracket

  use ModMpi
  use ModUtilities,    ONLY: CON_stop
  implicit none
  SAVE
  private ! except
  interface explicit
     module procedure explicit2
     module procedure explicit3
     module procedure explicit4
  end interface explicit
  !PUBLIC MEMBER FUNCTION:
  public :: explicit

  ! If DtIn results in CFL>CflMax, the time step is reduced
  real, parameter :: CflMax = 0.990
  character(LEN=*), parameter:: NameMod = 'ModPoissonBracket'
  logical, parameter :: UseLimiter = .true. ! false: switch off limiters
  ! The second order in time scheme is TVD only with  UseMinmodBeta=.true.,
  ! which scheme produces an extra dissipation  noticeable at 0.5<CFL<1,
  ! where \delta^-f changes sign. With UseMinmodBeta=.false. the TVD
  ! property is proved only at infinitesimal timestep.
  logical, parameter ::  UseMinmodBeta  = .true.
  logical, parameter ::  UseKoren = .false.
  real, parameter :: cTol = 1.0e-14

contains
  !============================================================================
  real function minmod(Arg1, Arg2)
    real, intent(in) :: Arg1, Arg2
    !--------------------------------------------------------------------------
    if(abs(Arg1)<=cTol)then
       minmod = 0.0
    elseif(.not.UseLimiter)then
       minmod = Arg1
    else
       minmod = (sign(0.5, Arg1) + sign(0.5, Arg2))*min(abs(Arg1), abs(Arg2))
    end if
  end function minmod
  !============================================================================
  real function half_beta(Cfl, DownwindDeltaMinusF, UpwindDeltaMinusF)
    !  To switch between two choices of beta-parameter:
    !  1. minmod, if UseMinmodBeta=.true.
    !  2. DownwindDeltaMinusF  otherwise
    real, intent(in) :: Cfl, DownwindDeltaMinusF, UpwindDeltaMinusF
    !--------------------------------------------------------------------------
    half_beta = 1.0 - 1.0e-14 - Cfl*(1.0 - minmod( &
         DownwindDeltaMinusF/UpwindDeltaMinusF, 1.0))
  end function half_beta
  !============================================================================
  real function betalimiter(DeltaF, UpwindDeltaF,&
       Cfl, DownwindDeltaMinusF, UpwindDeltaMinusF)
    !  \delta^-f in the neighboring cells
    real :: HalfBeta
    real, intent(in) :: DeltaF                ! f_j -f
    real, intent(in) :: UpwindDeltaF     ! f - f_j^\prime at the opposite face
    real, intent(in) :: Cfl, DownwindDeltaMinusF, UpwindDeltaMinusF
    real :: SignDeltaF, AbsDeltaF, AbsDeltaFLimited
    !--------------------------------------------------------------------------
    if(abs(DeltaF) <=cTol)then
       betalimiter = 0.0
    elseif(.not.UseLimiter)then
       betalimiter = 0.50*DeltaF
    else
       SignDeltaF = sign(1.0, DeltaF); AbsDeltaF = abs(DeltaF)
       if(UseKoren)then
          AbsDeltaFLimited = max(2*AbsDeltaF + &
               SignDeltaF*UpwindDeltaF,0.0)/6.0
       else
          AbsDeltaFLimited = 0.5*max(AbsDeltaF, SignDeltaF*UpwindDeltaF)
       end if
       if(UseMinmodBeta.and.CFL > 0.0.and.UpwindDeltaMinusF*DeltaF < 0.0)then
          HalfBeta = max(half_beta(Cfl, (1.0 - Cfl)*DownwindDeltaMinusF, &
               Cfl*UpwindDeltaMinusF), &
               1.0 - 1.0e-14 + CFL*UpwindDeltaMinusF/DeltaF)
       else
          HalfBeta = 1.0 - 1.0e-14
       end if
       betalimiter = SignDeltaF*min(AbsDeltaF*HalfBeta,AbsDeltaFLimited)
    end if
  end function betalimiter
  !============================================================================
  subroutine explicit2(nI, nJ, VDF_G, Volume_G, Source_C,    &
       Hamiltonian12_N, dHamiltonian01_FX, dHamiltonian02_FY,&
       DVolumeDt_G,                                          &
       IsSteadyState,                                        &
       DtIn, CFLIn, DtOut, DtOut_C, IsPeriodicIn_D)
    ! solve the contribution to the
    ! numerical flux from a single Poisson bracket,
    ! df/dq_l dH/dp_l - df/dp_l dH/dq_l
    integer, intent(in) :: nI    !# of cells along coordinate 1
    integer, intent(in) :: nJ    !# of cells along coord. 2 (momentum)
    ! Distribution function with gc. Two layers of face ghostcels
    ! and one level of corner ghost cells are used
    real, intent(in) :: VDF_G(-1:nI+2,-1:nJ+2)

    ! Hamiltonian function in nodes. One layer of ghost nodes is used
    ! Hamiltonian function for the Poisson bracket \{f,H_{12}}_{x,y}
    !    Node-centered at XY plane, cell-centered with respect to Z
    !    (In other words, Z-aligned-edge-centered)
    real, optional, intent(in) :: Hamiltonian12_N(-1:nI+1,-1:nJ+1)

    ! Increment in the Hamiltonian function for the Poisson bracket
    ! with respect to time, \{f,H_{01}\}_{t,x}. Is face-X centered.
    ! One layer of the ghost faces is needed.
    real, optional, intent(in) :: dHamiltonian01_FX(-1:nI+1, 0:nJ+1)
    ! Increment in the Hamiltonian function for the Poisson bracket
    ! with respect to time, \{f,H_{02}\}_{t,y}. Is face-Y centered.
    ! One layer of the ghost faces is needed.
    real, optional, intent(in) :: dHamiltonian02_FY( 0:nI+1,-1:nJ+1)

    ! Cell volume. One layer of face ghost cells is used
    real,           intent(in) :: Volume_G(0:nI+1,0:nJ+1)
    ! If non-canonical variables are used with time-dependent Jacobian,
    ! the cell volume changes in time. Need the volume derivative
    real, optional, intent(in) :: DVolumeDt_G(0:nI+1,0:nJ+1)
    ! Contribution to the conservative source (flux divergence) from
    ! the Poisson Bracket:
    !
    ! f(t+dt) - f(t) = Source_C
    real,           intent(out):: Source_C(1:nI,1:nJ)
    ! Logical to distinguish a time accurate mode from stedy-state one
    logical, optional, intent(in) :: IsSteadyState
    real, optional, intent(inout) :: DtIn   ! Options to set time step
    real, optional, intent(in)    :: CFLIn  ! Options to set time step
    real, optional, intent(out)   :: DtOut  ! Option to report time step
    real, optional, intent(out)   :: DtOut_C(nI,nJ)
    logical, optional, intent(in) :: IsPeriodicIn_D(:)
    character(len=*), parameter:: NameSub = 'explicit2'
    !--------------------------------------------------------------------------
    call explicit4(nI, nJ, 1, 1, VDF_G, Volume_G, Source_C,  &
         Hamiltonian12_N=Hamiltonian12_N,                    &
         dHamiltonian01_FX= dHamiltonian01_FX,               &
         dHamiltonian02_FY= dHamiltonian02_FY,               &
         DVolumeDt_G=DVolumeDt_G,                            &
         IsSteadyState=IsSteadyState,                        &
         DtIn=DtIn,                                          &
         CFLIn=CFLIn,                                        &
         DtOut=DtOut,                                        &
         DtOut_C=DtOut_C,                                    &
         IsPeriodicIn_D=IsPeriodicIn_D)
  end subroutine explicit2
  !============================================================================
  subroutine explicit3(nI, nJ, nK, VDF_G, Volume_G, Source_C,            &
       Hamiltonian12_N, Hamiltonian13_N, Hamiltonian23_N,                &
       dHamiltonian01_FX, dHamiltonian02_FY, dHamiltonian03_FZ,          &
       DVolumeDt_G,                                                      &
       IsSteadyState,                                                    &
       DtIn, CFLIn, DtOut, DtOut_C, IsPeriodicIn_D)

    ! solve the contribution to the numerical flux from multiple Poisson
    ! brackets, 1,2,3 enumerate phase coordinates,  0 relating to time.

    integer, intent(in) :: nI     !# of cells along coordinate 1
    integer, intent(in) :: nJ     !# of cells along coordinate 2
    integer, intent(in) :: nK     !# of cells along coordinate 3
    integer :: iKStart , iKLast

    ! Distribution function with gc. Two layers of face ghostcels
    ! and one level of corner ghost cells are used
    real, intent(in) :: VDF_G(-1:nI+2,-1:nJ+2,-1:nK+2)

    ! Hamiltonian functions in nodes. One layer of ghost nodes is used.
    ! 1. Hamiltonian function for the Poisson bracket \{f,H_{12}}_{x,y}
    !    Node-centered at XY plane, cell-centered with respect to Z
    !    (In other words, Z-aligned-edge-centered)
    real, optional, intent(in) :: Hamiltonian12_N(-1:nI+1,-1:nJ+1,0:nK+1)
    ! 2. Hamiltonian function for the Poisson bracket \{f,H_{13}}_{x,z}
    !    Node-centered at XZ plane, cell-centered with respect to Y
    !    (In other words, Y-aligned-edge-centered)
    real, optional, intent(in) :: Hamiltonian13_N(-1:nI+1,0:nJ+1,-1:nK+1)
    ! 3. Hamiltonian function for the Poisson bracket \{f,H_{23}}_{y,z}
    !    Node-centered at YZ plane, cell-centered with respect to X
    !    (In other words, X-aligned-edge-centered)
    real, optional, intent(in) :: Hamiltonian23_N(0:nI+1,-1:nJ+1,-1:nK+1)

    ! 1. Increment in the Hamiltonian function for the Poisson bracket
    ! with respect to time, \{f,H_{01}\}_{t,x}. Is face-X centered.
    ! One layer of the ghost faces is needed.
    real, optional, intent(in) :: dHamiltonian01_FX(-1:nI+1,0:nJ+1,0:nK+1)
    ! 2. Increment in the Hamiltonian function for the Poisson bracket
    ! with respect to time, \{f,H_{02}\}_{t,y}. Is face-Y centered.
    ! One layer of the ghost faces is needed.
    real, optional, intent(in) :: dHamiltonian02_FY(0:nI+1,-1:nJ+1,0:nK+1)
    ! 3. Increment in the Hamiltonian function for the Poisson bracket
    ! with respect to time, \{f,H_{03}\}_{t,z}. Is face-Z centered.
    ! One layer of the ghost faces is needed.
    real, optional, intent(in) :: dHamiltonian03_FZ(0:nI+1,0:nJ+1,-1:nK+1)

    ! Total Volume. One layer of face ghost cells is used
    real, intent(in) :: Volume_G(0:nI+1,0:nJ+1,0:nK+1)
    ! If non-canonical variables are used with time-dependent Jacobian,
    ! the cell volume changes in time. Need the volume derivative
    real, optional, intent(in) :: DVolumeDt_G(0:nI+1,0:nJ+1,0:nK+1)

    ! Contribution to the conservative source (flux divergence) for
    ! the Poisson Bracket:
    ! send the source_c back to the main code
    real, intent(out) :: Source_C(1:nI, 1:nJ, 1:nK)

    !OPTIONAL PARAMETERS:
    ! Logical to distinguish a time accurate mode from stedy-state one
    logical, optional, intent(in) :: IsSteadyState

    real, optional, intent(inout) :: DtIn    ! Options to set time step
    real, optional, intent(in)    :: CFLIn   ! Options to set time step
    real, optional, intent(out)   :: DtOut   ! Option to report time step
    real, optional, intent(out)   :: DtOut_C(nI,nJ,nK)
    logical, optional, intent(in) :: IsPeriodicIn_D(:)
    character(len=*), parameter:: NameSub = 'explicit3'
    !--------------------------------------------------------------------------
    call explicit4(nI, nJ, nK, 1, VDF_G, Volume_G, Source_C, &
         Hamiltonian12_N=Hamiltonian12_N,                    &
         Hamiltonian13_N=Hamiltonian13_N,                    &
         Hamiltonian23_N=Hamiltonian23_N,                    &
         dHamiltonian01_FX= dHamiltonian01_FX,               &
         dHamiltonian02_FY= dHamiltonian02_FY,               &
         dHamiltonian03_FZ= dHamiltonian03_FZ,               &
         DVolumeDt_G=DVolumeDt_G,                            &
         IsSteadyState=IsSteadyState,                        &
         DtIn=DtIn,                                          &
         CFLIn=CFLIn,                                        &
         DtOut=DtOut,                                        &
         DtOut_C=DtOut_C,                                    &
         IsPeriodicIn_D=IsPeriodicIn_D)
  end subroutine explicit3
  !============================================================================
  subroutine explicit4(nI, nJ, nK, nP, VDF_G, Volume_G, Source_C,        &
       Hamiltonian12_N, Hamiltonian13_N, Hamiltonian23_N,                &
       Hamiltonian14_N, Hamiltonian24_N, Hamiltonian34_N,                &
       dHamiltonian01_FX, dHamiltonian02_FY, dHamiltonian03_FZ,          &
       dHamiltonian04_FP, DVolumeDt_G,                                   &
       IsSteadyState,                                                    &
       DtIn, CFLIn, DtOut, DtOut_C, IsPeriodicIn_D)

    ! solve the contribution to the numerical flux from multiple Poisson
    ! brackets, 1,2,3 enumerate phase coordinates,  0 relating to time.

    integer, intent(in) :: nI     !# of cells along coordinate 1
    integer, intent(in) :: nJ     !# of cells along coordinate 2
    integer, intent(in) :: nK     !# of cells along coordinate 3
    integer, intent(in) :: nP     !# of cells along coordinate 4
    integer :: iKStart , iKLast, iPStart, iPLast

    ! Distribution function with gc. Two layers of face ghostcels
    ! and one level of corner ghost cells are used
    real, intent(in) :: VDF_G(-1:nI+2,-1:nJ+2,&
         -1 + 2*(1/nK):nK + 2*(1 - 1/nK),-1 + 2*(1/nP):nP +2*(1 - 1/nP))

    ! Hamiltonian functions in nodes. One layer of ghost nodes is used.
    ! 1. Hamiltonian function for the Poisson bracket \{f,H_{12}}_{x,y}
    !    Node-centered at XY plane, cell-centered with respect to Z
    !    (In other words, Z-aligned-edge-centered)
    real, optional, intent(in) :: Hamiltonian12_N(-1:nI+1,-1:nJ+1,&
         1/nK:nK+1-1/nK,1/nP:nP+1-1/nP)
    ! 2. Hamiltonian function for the Poisson bracket \{f,H_{13}}_{x,z}
    !    Node-centered at XZ plane, cell-centered with respect to Y
    !    (In other words, Y-aligned-edge-centered)
    real, optional, intent(in) :: Hamiltonian13_N(-1:nI+1, 0:nJ+1,&
         -1:nK+1,1/nP:nP+1-1/nP)
    ! 3. Hamiltonian function for the Poisson bracket \{f,H_{23}}_{y,z}
    !    Node-centered at YZ plane, cell-centered with respect to X
    !    (In other words, X-aligned-edge-centered)
    real, optional, intent(in) :: Hamiltonian23_N( 0:nI+1,-1:nJ+1,&
         -1:nK+1,1/nP:nP+1-1/nP)
    ! Hamiltonian functions in nodes. One layer of ghost nodes is used.
    ! 4. Hamiltonian function for the Poisson bracket \{f,H_{14}}_{x,P}
    !    Node-centered at XP plane, cell-centered with respect to y,z
    !    (In other words, YZ-aligned-edge-centered)
    real, optional, intent(in) :: Hamiltonian14_N(-1:nI+1, 0:nJ+1,&
         0:nK+1,-1:nP+1)
    ! 5. Hamiltonian function for the Poisson bracket \{f,H_{24}}_{y,P}
    !    Node-centered at YP plane, cell-centered with respect to XZ
    !    (In other words, XZ-aligned-edge-centered)
    real, optional, intent(in) :: Hamiltonian24_N(0:nI+1,-1:nJ+1,&
         0:nK+1,-1:nP+1)
    ! 6. Hamiltonian function for the Poisson bracket \{f,H_{34}}_{z,P}
    !    Node-centered at ZP plane, cell-centered with respect to XY
    !    (In other words, XY-aligned-edge-centered)
    real, optional, intent(in) :: Hamiltonian34_N( 0:nI+1, 0:nJ+1,&
         -1:nK+1,-1:nP+1)

    ! 1. Increment in the Hamiltonian function for the Poisson bracket
    ! with respect to time, \{f,H_{01}\}_{t,x}. Is face-X centered.
    ! One layer of the ghost faces is needed.
    real, optional, intent(in) :: dHamiltonian01_FX(-1:nI+1,0:nJ+1,&
         1/nK:nK+1-1/nK,1/nP:nP+1-1/nP)
    ! 2. Increment in the Hamiltonian function for the Poisson bracket
    ! with respect to time, \{f,H_{02}\}_{t,y}. Is face-Y centered.
    ! One layer of the ghost faces is needed.
    real, optional, intent(in) :: dHamiltonian02_FY(0:nI+1,-1:nJ+1,&
         1/nK:nK+1-1/nK,1/nP:nP+1-1/nP)
    ! 3. Increment in the Hamiltonian function for the Poisson bracket
    ! with respect to time, \{f,H_{03}\}_{t,z}. Is face-Z centered.
    ! One layer of the ghost faces is needed.
    real, optional, intent(in) :: dHamiltonian03_FZ( 0:nI+1,0:nJ+1,&
         -1:nK+1,1/nP:nP+1-1/nP)
    ! 4. Increment in the Hamiltonian function for the Poisson bracket
    ! with respect to time, \{f,H_{04}\}_{t,P}. Is face-P centered.
    ! One layer of the ghost faces is needed.
    real, optional, intent(in) :: dHamiltonian04_FP( 0:nI+1,0:nJ+1,&
         0:nK+1,-1:nP+1)

    ! Total Volume. One layer of face ghost cells is used
    real, intent(in) :: Volume_G(0:nI+1,0:nJ+1,1/nK:nK+1-1/nK,1/nP:nP+1-1/nP)
    ! If non-canonical variables are used with time-dependent Jacobian,
    ! the cell volume changes in time. Need the volume derivative
    real, optional, intent(in) :: DVolumeDt_G(0:nI+1,0:nJ+1,&
         1/nK:nK+1-1/nK,1/nP:nP+1-1/nP)

    ! Contribution to the conservative source (flux divergence) for
    ! the Poisson Bracket:
    ! send the source_c back to the main code
    real, intent(out) :: Source_C(1:nI, 1:nJ, 1:nK, 1:nP)

    !OPTIONAL PARAMETERS:
    ! Logical to distinguish a time accurate mode from stedy-state one
    logical, optional, intent(in) :: IsSteadyState

    real, optional, intent(inout) :: DtIn    ! Options to set time step
    real, optional, intent(in)    :: CFLIn   ! Options to set time step
    real, optional, intent(out)   :: DtOut   ! Option to report time step
    real, optional, intent(out)   :: DtOut_C(nI,nJ,nK,nP)
    logical, optional, intent(in) :: IsPeriodicIn_D(:)
    ! Local variables
    logical :: UseTimeDependentVolume = .false. !=present(DVolumeDt_G)

    ! Inverse volume. One layer of face ghost cells is used
    real :: vInv_G(0:nI+1,0:nJ+1,1/nK:nK+1-1/nK,1/nP:nP+1-1/nP)
    ! Boundary conddition:
    logical :: IsPeriodic_D(4 - 1/nK - 1/nP)
    ! Loop variables:
    integer :: i, j, k, iP, iDim, nDim
    ! Variations of VDF (one layer of ghost cell values):
    real,dimension(0:nI+1, 0:nJ+1, 1/nK:nK+1-1/nK, 1/nP:nP+1-1/nP) :: &
         DeltaMinusF_G, DeltaMinusH_G, SumDeltaHPlus_G, SumDeltaHMinus_G
    ! Face-centered vriations of Hamiltonian functions, with
    ! one layer of ghost faces
    real :: DeltaH_DG(4,-1:nI+1,-1:nJ+1,-1+2*(1/nK):nK+1-1/nK,&
         -1+2*(1/nP):nP+1-1/nP)
    ! Fluxes:
    ! Sum of \delta^+H fluxes, to be limited
    real :: SumFluxPlus
    ! Sum of all second order fluxes
    real :: SumFlux2_G(0:nI+1,0:nJ+1,1/nK:nK+1-1/nK,1/nP:nP+1-1/nP)
    integer, parameter :: DeltaH_ = 1, Flux_ = 2
    integer, parameter :: iShift_DS(4,8) = reshape(&
         [-1,  0,  0,  0,&
         1,  0,  0,  0,&
         0, -1,  0,  0,&
         0,  1,  0,  0,&
         0,  0, -1,  0,&
         0,  0,  1,  0,&
         0,  0,  0, -1,&
         0,  0,  0,  1], [4,8])
    real    :: Buff_VSG(DeltaH_:Flux_,8)
    integer :: nFlux, iFlux, iSide, iSide_SG(8), iD_D(4), iU_D(4)
    ! Local CFL number:
    real :: CFLCoef_G(0:nI+1,0:nJ+1,1/nK:nK+1-1/nK,1/nP:nP+1-1/nP), CFLLocal
    ! Time step (global, local)
    real :: Dt, TimeStep_G(0:nI+1,0:nJ+1,1/nK:nK+1-1/nK,1/nP:nP+1-1/nP)
    ! Misc:
    ! Sum of major contributions
    real    :: SumMajor
    ! Misc:
    real    :: VDF
    ! Gamma-limiter
    real    :: DeltaPlusFLimited, Gamma
    character(len=*), parameter:: NameSub = 'explicit4'
    !--------------------------------------------------------------------------
    if(present(DtIn))then
       Dt = DtIn
    else
       if(.not.present(CflIn))call CON_stop(&
            'Either CflIn or DtIn should be provided in '//NameSub)
    end if
    UseTimeDependentVolume = present(DVolumeDt_G)
    if(present(IsSteadyState))then
       if(.not.IsSteadyState)call CON_stop(&
            'If present IsSteadyState parameter, it must be be .true.')
       if(present(DtIn).or.present(DtOut))call CON_stop(&
            'DtIn and DtOut make no sense in the steady-state mode')
       if(.not.present(CflIn))call CON_stop(&
            'CflIn is required in the steady-state mode')
       if(UseTimeDependentVolume)call CON_stop(&
            'Volume cannot be time dependent in the steady-state mode')
    end if
    IsPeriodic_D = .false.
    if(present(IsPeriodicIn_D))&
         IsPeriodic_D = IsPeriodicIn_D(1:4-1/nK-1/nP)
    iKStart  = 1/nK ;  iKLast  = nK + 1 - 1/nK; nDim = 4 - 1/nK - 1/nP
    iPStart  = 1/nP ;  iPLast  = nP + 1 - 1/nP
    vInv_G = 1.0/Volume_G

    ! Nullify arrays:
    DeltaH_DG  = 0.0
    TimeStep_G = 0.0
    ! Bracket {F,H12}_{x,y}   Bracket {F,H13}_{x,z}    Bracket {F,H23}_{y,z}
    ! Hamiltonian 12 (xy)     Hamiltonian 13 (xz)      Hamiltonian 23 (yz)
    ! y                       z                        z
    ! 1---------<---------    1---------<----------    1----------<---------
    ! |                  |    |                   |    |                   |
    ! v       1,1,1      ^    v       1,1,1       |    v       1,1,1       ^
    ! |                  |    |                   |    |                   |
    ! |                  |    |                   |    |                   |
    ! 0--------->--------1x   0---------->--------1x   0---------->--------1y
    if (present(Hamiltonian12_N)) then
       DeltaH_DG(1,-1:nI+1, 0:nJ+1,iKStart:iKLast,iPStart:iPLast)          = &
            DeltaH_DG(1,-1:nI+1, 0:nJ+1,iKStart:iKLast,iPStart:iPLast)     + &
            Hamiltonian12_N(-1:nI+1, 0:nJ+1,iKStart:iKLast,iPStart:iPLast) - &
            Hamiltonian12_N(-1:nI+1,-1:nJ  ,iKStart:iKLast,iPStart:iPLast)
       DeltaH_DG(2, 0:nI+1,-1:nJ+1,iKStart:iKLast,iPStart:iPLast)          = &
            DeltaH_DG(2, 0:nI+1,-1:nJ+1,iKStart:iKLast,iPStart:iPLast)     + &
            Hamiltonian12_N(-1:nI  ,-1:nJ+1,iKStart:iKLast,iPStart:iPLast) - &
            Hamiltonian12_N(0:nI+1 ,-1:nJ+1,iKStart:iKLast,iPStart:iPLast)
    end if
    if (present(Hamiltonian13_N)) then
       DeltaH_DG(1,-1:nI+1, 0:nJ+1, 0:nK+1,iPStart:iPLast)                 = &
            DeltaH_DG(1,-1:nI+1, 0:nJ+1, 0:nK+1,iPStart:iPLast)            + &
            Hamiltonian13_N(-1:nI+1, 0:nJ+1, 0:nK+1,iPStart:iPLast)        - &
            Hamiltonian13_N(-1:nI+1, 0:nJ+1,-1:nK,iPStart:iPLast  )
       DeltaH_DG(3, 0:nI+1, 0:nJ+1,-1:nK+1,iPStart:iPLast)                 = &
            DeltaH_DG(3, 0:nI+1, 0:nJ+1,-1:nK+1,iPStart:iPLast)            + &
            Hamiltonian13_N(-1:nI  , 0:nJ+1,-1:nK+1,iPStart:iPLast)        - &
            Hamiltonian13_N( 0:nI+1, 0:nJ+1,-1:nK+1,iPStart:iPLast)
    end if
    if (present(Hamiltonian23_N)) then
       DeltaH_DG(2, 0:nI+1,-1:nJ+1, 0:nK+1,iPStart:iPLast)                 = &
            DeltaH_DG(2, 0:nI+1,-1:nJ+1, 0:nK+1,iPStart:iPLast)            + &
            Hamiltonian23_N( 0:nI+1,-1:nJ+1, 0:nK+1,iPStart:iPLast)        - &
            Hamiltonian23_N( 0:nI+1,-1:nJ+1,-1:nK,iPStart:iPLast  )
       DeltaH_DG(3, 0:nI+1, 0:nJ+1,-1:nK+1,iPStart:iPLast)                 = &
            DeltaH_DG(3, 0:nI+1, 0:nJ+1,-1:nK+1,iPStart:iPLast)            + &
            Hamiltonian23_N( 0:nI+1,-1:nJ  ,-1:nK+1,iPStart:iPLast)        - &
            Hamiltonian23_N( 0:nI+1, 0:nJ+1,-1:nK+1,iPStart:iPLast)
    end if
    if (present(Hamiltonian14_N)) then
       DeltaH_DG(1,-1:nI+1, 0:nJ+1,0:nK+1,0:nP+1)                          = &
            DeltaH_DG(1,-1:nI+1, 0:nJ+1,0:nK+1,0:nP+1)                     + &
            Hamiltonian14_N(-1:nI+1, 0:nJ+1,0:nK+1,0:nP+1)                 - &
            Hamiltonian14_N(-1:nI+1, 0:nJ+1,0:nK+1,-1:nP)
       DeltaH_DG(4, 0:nI+1,0:nJ+1,0:nK+1,-1:nP+1)                          = &
            DeltaH_DG(4, 0:nI+1,0:nJ+1,0:nK+1,-1:nP+1)                     + &
            Hamiltonian14_N(-1:nI  , 0:nJ+1, 0:nK+1,-1:nP+1)               - &
            Hamiltonian14_N(0:nI+1 , 0:nJ+1, 0:nK+1,-1:nP+1)
    end if
    if (present(Hamiltonian24_N)) then
       DeltaH_DG(2, 0:nI+1,-1:nJ+1, 0:nK+1,0:nP+1)                         = &
            DeltaH_DG(2,0:nI+1,-1:nJ+1,0:nK+1,0:nP+1)                      + &
            Hamiltonian24_N(0:nI+1,-1:nJ+1, 0:nK+1,0:nP+1)                 - &
            Hamiltonian24_N(0:nI+1,-1:nJ+1, 0:nK+1,-1:nP)
       DeltaH_DG(4, 0:nI+1, 0:nJ+1, 0:nK+1,-1:nP+1)                        = &
            DeltaH_DG(4, 0:nI+1, 0:nJ+1, 0:nK+1,-1:nP+1)                   + &
            Hamiltonian24_N( 0:nI+1,-1:nJ  ,0:nK+1,-1:nP+1)                - &
            Hamiltonian24_N( 0:nI+1, 0:nJ+1,0:nK+1,-1:nP+1)
    end if
    if (present(Hamiltonian34_N)) then
       DeltaH_DG(3,0:nI+1,0:nJ+1,-1:nK+1,0:nP+1)                           = &
            DeltaH_DG(3,0:nI+1,0:nJ+1,-1:nK+1,0:nP+1)                      + &
            Hamiltonian34_N( 0:nI+1,0:nJ+1,-1:nK+1,0:nP+1)                 - &
            Hamiltonian34_N( 0:nI+1,0:nJ+1,-1:nK+1,-1:nP)
       DeltaH_DG(4,0:nI+1,0:nJ+1,0:nK+1,-1:nP+1)                           = &
            DeltaH_DG(4,0:nI+1,0:nJ+1,0:nK+1,-1:nP+1)                      + &
            Hamiltonian34_N(0:nI+1,0:nJ+1,-1:nK,-1:nP+1)                   - &
            Hamiltonian34_N(0:nI+1,0:nJ+1,0:nK+1,-1:nP+1)
    end if
    ! Bracket {F,H01}_t,x    Bracket {F,H02}_t,y      Bracket {F,H03}_t,z
    ! Hamiltonian 01 (tx)    Hamiltonian 02 (ty)      Hamiltonian 03 (tz)
    ! x                      y                        z
    ! ----------<--------    ----------<---------     -----------<-------
    !
    !
    !
    !
    ! ---------->--------t   ---------->---------t    ----------->-------t
    if (present(dHamiltonian01_FX))&
         DeltaH_DG(1,-1:nI+1, 0:nJ+1,iKStart:iKLast,iPStart:iPLast)     = &
         DeltaH_DG(1,-1:nI+1, 0:nJ+1,iKStart:iKLast,iPStart:iPLast)     + &
         dHamiltonian01_FX
    if (present(dHamiltonian02_FY))&
         DeltaH_DG(2, 0:nI+1,-1:nJ+1,iKStart:iKLast,iPStart:iPLast)     = &
         DeltaH_DG(2, 0:nI+1,-1:nJ+1,iKStart:iKLast,iPStart:iPLast)     + &
         dHamiltonian02_FY
    if (present(dHamiltonian03_FZ))&
         DeltaH_DG(3, 0:nI+1, 0:nJ+1,-1:nK+1,iPStart:iPLast)            = &
         DeltaH_DG(3, 0:nI+1, 0:nJ+1,-1:nK+1,iPStart:iPLast)            + &
         dHamiltonian03_FZ
    if (present(dHamiltonian04_FP))&
         DeltaH_DG(4, 0:nI+1, 0:nJ+1, 0:nK+1,-1:nP+1)                   = &
         DeltaH_DG(4, 0:nI+1, 0:nJ+1, 0:nK+1,-1:nP+1)                   + &
         dHamiltonian04_FP
    ! Cleanup
    where(abs(DeltaH_DG)<=cTol) DeltaH_DG = 0.0

    ! Now, for each cell the value of DeltaH for face in positive
    ! directions of i and j may be found in the arrays, for
    ! negative directions the should be taken with opposite sign
    ! Nullify arrays:
    DeltaMinusF_G = 0.0
    DeltaMinusH_G = 0.0
    SumDeltaHPlus_G = 0.0
    SumDeltaHMinus_G = 0.0
    ! Calculate DeltaMinusF_G and SumDeltaHPlus_G
    do iDim = 1, nDim
       ! Calculate contributions from up faces to
       ! SumDeltaHPlus_G and DeltaMinusF_G
       ! Upper & lower bounds for +1 (may be with diffent meanings elsewhere)
       iU_D = [nI+1,nJ+1,iKLast,iPLast]; iU_D(iDim) = iU_D(iDim)+1
       iD_D = [0,0,iKStart,iPStart]; iD_D(iDim) = iD_D(iDim)+1
       DeltaMinusH_G = min(0.0, DeltaH_DG(iDim, &
            0:nI+1, 0:nJ+1, iKStart:iKLast, iPStart:iPLast))
       SumDeltaHMinus_G = SumDeltaHMinus_G + DeltaMinusH_G
       DeltaMinusF_G = DeltaMinusF_G + DeltaMinusH_G* &
            (VDF_G(iD_D(1):iU_D(1), iD_D(2):iU_D(2), &
            iD_D(3):iU_D(3), iD_D(4):iU_D(4)) - &
            VDF_G(0:nI+1, 0:nJ+1, iKStart:iKLast, iPStart:iPLast))
       SumDeltaHPlus_G = SumDeltaHPlus_G + max(0.0, DeltaH_DG(iDim, &
            0:nI+1, 0:nJ+1, iKStart:iKLast, iPStart:iPLast))

       ! Calculate contributions from down faces to
       ! SumDeltaHPlus and DeltaMinusF_G
       ! Upper & lower bounds for -1 (may be with diffent meanings elsewhere)
       iU_D = [nI+1,nJ+1,iKLast,iPLast]; iU_D(iDim) = iU_D(iDim)-1
       iD_D = [0,0,iKStart,iPStart]; iD_D(iDim) = iD_D(iDim)-1
       DeltaMinusH_G = min(0.0, -DeltaH_DG(iDim, &
            iD_D(1):iU_D(1), iD_D(2):iU_D(2), &
            iD_D(3):iU_D(3), iD_D(4):iU_D(4)))
       SumDeltaHMinus_G = SumDeltaHMinus_G + DeltaMinusH_G
       DeltaMinusF_G = DeltaMinusF_G + DeltaMinusH_G* &
            (VDF_G(iD_D(1):iU_D(1), iD_D(2):iU_D(2), &
            iD_D(3):iU_D(3), iD_D(4):iU_D(4)) - &
            VDF_G(0:nI+1, 0:nJ+1, iKStart:iKLast, iPStart:iPLast))
       SumDeltaHPlus_G = SumDeltaHPlus_G + max(0.0, &
            -DeltaH_DG(iDim, iD_D(1):iU_D(1), iD_D(2):iU_D(2), &
            iD_D(3):iU_D(3), iD_D(4):iU_D(4)))
    end do

    ! Get local CFLs
    where(SumDeltaHMinus_G == 0.0)
       DeltaMinusF_G = 0.0
       CFLCoef_G = 0.0
    elsewhere
       DeltaMinusF_G = -DeltaMinusF_G/SumDeltaHMinus_G
       ! UseTimeDependentVolume: Local CFLs are expressed via SumDeltaHMinus_G
       CFLCoef_G = -SumDeltaHMinus_G
    end where
    ! Not UseTimeDependentVolume: Local CFLs are expressed via *vInv_G
    if(.not.UseTimeDependentVolume) CFLCoef_G = CFLCoef_G*vInv_G

    ! Set CFL and time step
    if(UseTimeDependentVolume)then
       if(present(DtIn))then
          if(present(DtOut).and.present(CFLIn))&
               DtOut = CFLIn/maxval(vInv_G(1:nI,1:nJ,1:nK,1:nP)*&
               (CFLCoef_G(1:nI,1:nJ,1:nK,1:nP) - &
               CFLIn*DVolumeDt_G(1:nI,1:nJ,1:nK,1:nP)))
          ! Calculate the CFL factor with given Dt:
          vInv_G = 1.0/(Volume_G + Dt*DVolumeDt_G)
          CFLCoef_G = Dt*vInv_G*CFLCoef_G

          ! Check if the CFL satisfies the stability criterion
          if(maxval(CFLCoef_G(1:nI,1:nJ,1:nK,1:nP)) > CFLMax)then
             ! Restore CFLCoef_G and vInv
             CFLCoef_G = CFLCoef_G/(Dt*vInv_G)
             vInv_G = 1/Volume_G
             ! Reduce the time step using equation
             ! CFLMax = \Delta t*(-\sum\delta^-H)/(\Delta t*dV/dt + V)
             DtIn = CFLMax/maxval(vInv_G(1:nI,1:nJ,1:nK,1:nP)*&
                  ( CFLCoef_G(1:nI,1:nJ,1:nK,1:nP) &
                  - CFLMax*DVolumeDt_G(1:nI,1:nJ,1:nK,1:nP)))
             Dt = DtIn
             ! Calculate the CFL factor with given Dt:
             vInv_G = 1.0/(Volume_G + Dt*DVolumeDt_G)
             CFLCoef_G = Dt*vInv_G*CFLCoef_G
          end if
          if(present(DtOut).and..not.present(CFLIn))&
               DtOut = Dt
          ! Time accurate mode, equal time step everywhere
          TimeStep_G = Dt
       else
          ! Solve minimal time step from equation
          ! CFLIn = \Delta t*(-\sum\delta^-H)/(\Delta t*dV/dt + V)

          Dt = CFLIn*minval(  Volume_G(1:nI,1:nJ,1:nK,1:nP)/&
               max( CFLCoef_G(1:nI,1:nJ,1:nK,1:nP) -        &
               CFLIn*DVolumeDt_G(1:nI,1:nJ,1:nK,1:nP), 1.0D-30) )
          if(present(DtOut))DtOut = Dt
          ! Time accurate mode, equal time step everywhere
          TimeStep_G = Dt
          ! Calculate the volume at upper time level
          ! V(+\Delta t):
          vInv_G = 1.0/(Volume_G + TimeStep_G*DVolumeDt_G)
          CFLCoef_G = TimeStep_G*vInv_G*CFLCoef_G
       end if
    else
       if(present(DtIn))then
          if(present(DtOut ))&
               DtOut = CFLIn/maxval(CFLCoef_G(1:nI,1:nJ,1:nK,1:nP))
          CFLCoef_G = Dt*CFLCoef_G
          if(maxval(CFLCoef_G(1:nI,1:nJ,1:nK,1:nP)) > CFLMax)then
             ! Adjust time step and local CFL:
             Dt = Dt*CFLMax/maxval(CFLCoef_G(1:nI,1:nJ,1:nK,1:nP))
             CFLCoef_G = (Dt/DtIn)*CFLCoef_G
             DtIn = Dt
          end if
          ! Time accurate mode, equal time step everywhere
          TimeStep_G = Dt
       else
          where(CFLCoef_G/=0.0)TimeStep_G = CFLIn/CFLCoef_G
          if(.not.present(IsSteadyState))then
             Dt = minval(TimeStep_G(1:nI,1:nJ,1:nK,1:nP), &
                  MASK=TimeStep_G(1:nI,1:nJ,1:nK,1:nP)>0.0)
             if(present(DtOut))DtOut = Dt
             ! Time accurate mode, equal time step everywhere
             TimeStep_G = Dt
          end if
          CFLCoef_G = TimeStep_G*CFLCoef_G
       end if
    end if
    if(present(DtOut_C)) DtOut_C = TimeStep_G(1:nI,1:nJ,1:nK,1:nP)

    ! Calculate source = f(t+Dt) - f(t):
    ! First-order monotone scheme
    Source_C = -CFLCoef_G(1:nI,1:nJ,1:nK,1:nP)*&
         DeltaMinusF_G(1:nI,1:nJ,1:nK,1:nP)
    ! Second-order correction
    SumFlux2_G = 0.0
    do iP=1,nP; do k=1,nK; do j=1,nJ; do i=1,nI
       ! Limit and store fuxes across delta plus H faces from the given cell
       if(SumDeltaHPlus_G(i,j,k,iP)==0.0)CYCLE
       nFlux = 0
       do iDim = 1, nDim
          iSide = 2*iDim
          if(DeltaH_DG(iDim,i,j,k,iP) > 0.0)then
             nFlux = nFlux + 1; iSide_SG(nFlux) = iSide
             Buff_VSG(DeltaH_:Flux_,nFlux) = DeltaH_DG(iDim,i,j,k,iP)*&
                  [1.0, limiter(iSide,i,j,k,iP)]
          end if
          iSide = 2*iDim -1
          iD_D = [i,j,k,iP]; iD_D(iDim) = iD_D(iDim) - 1
          if(DeltaH_DG(iDim,iD_D(1),iD_D(2),iD_D(3),iD_D(4)) < 0.0 )then
             nFlux = nFlux + 1; iSide_SG(nFlux) = iSide
             Buff_VSG(DeltaH_:Flux_,nFlux) = &
                  (-DeltaH_DG(iDim,iD_D(1),iD_D(2),iD_D(3),iD_D(4)))*&
                  [1.0, limiter(iSide,i,j,k,iP)]
          end if
       end do
       SumFluxPlus = sum(Buff_VSG(Flux_,1:nFlux))
       ! The value of limited \deta^+H/2 to be achieved with gamma-limiter
       DeltaPlusFLimited = minmod(SumFluxPlus/SumDeltaHPlus_G(i,j,k,iP),&
            DeltaMinusF_G(i,j,k,iP))
       ! This is the major flux, having the same sign as SumFluxPlus
       SumMajor =  sum(Buff_VSG(Flux_,1:nFlux),&
            MASK= SumFluxPlus*Buff_VSG(Flux_,1:nFlux) > 0.0)
       if(abs(SumMajor) > 0.0)then
          Gamma = 1.0 + (DeltaPlusFLimited*SumDeltaHPlus_G(i,j,k,iP) - &
               SumFluxPlus)/SumMajor
          where(SumFluxPlus*Buff_VSG(Flux_,1:nFlux) > 0.0)&
               Buff_VSG(Flux_,1:nFlux) = Buff_VSG(Flux_,1:nFlux)*Gamma
       end if
       CFLLocal = CFLCoef_G(i,j,k,iP)
       do iFlux = 1, nFlux
          iD_D = [i,j,k,iP] + iShift_DS(:,iSide_SG(iFlux))
          SumFlux2_G(iD_D(1),iD_D(2),iD_D(3),iD_D(4)) =      &
               SumFlux2_G(iD_D(1),iD_D(2),iD_D(3),iD_D(4)) + &
               Buff_VSG(Flux_,iFlux) - &
               Buff_VSG(DeltaH_,iFlux)*CFLLocal*DeltaPlusFLimited
       end do
       SumFlux2_G(i,j,k,iP) = SumFlux2_G(i,j,k,iP) - &
            (1 - CFLLocal)*DeltaPlusFLimited*SumDeltaHPlus_G(i,j,k,iP)
    end do; end do; end do; end do
    ! Second order fluxes across the boundary
    if(IsPeriodic_D(1))then
       SumFlux2_G(1,1:nJ,1:nK,1:nP) = &
            SumFlux2_G(1,1:nJ,1:nK,1:nP) + SumFlux2_G(nI+1,1:nJ,1:nK,1:nP)
       SumFlux2_G(nI,1:nJ,1:nK,1:nP) = &
            SumFlux2_G(nI,1:nJ,1:nK,1:nP) + SumFlux2_G(0,1:nJ,1:nK,1:nP)
    else
       do iP=1,nP; do k=1,nK; do j = 1,nJ
          if(DeltaH_DG(1,0,j,k,iP) > 0.0)then
             SumFlux2_G(1,j,k,iP) = SumFlux2_G(1,j,k,iP) + &
                  DeltaH_DG(1,0,j,k,iP)*&
                  (1.0 - CFLCoef_G(0,j,k,iP))* minmod(&
                  limiter(2,0,j,k,iP), DeltaMinusF_G(0,j,k,iP))
          end if
          if(DeltaH_DG(1,nI,j,k,iP) < 0.0 )then
             SumFlux2_G(nI,j,k,iP)   = SumFlux2_G(nI,j,k,iP) - &
                  DeltaH_DG(1,nI,j,k,iP)*&
                  (1.0 - CFLCoef_G(nI+1,j,k,iP))*minmod(&
                  limiter(1,nI+1,j,k,iP), DeltaMinusF_G(nI+1,j,k,iP))
          end if
       end do; end do; end do
    end if
    if(IsPeriodic_D(2))then
       SumFlux2_G(1:nI,1,1:nK,1:nP) = &
            SumFlux2_G(1:nI,1,1:nK,1:nP) + SumFlux2_G(1:nI,nJ+1,1:nK,1:nP)
       SumFlux2_G(1:nI,nJ,1:nK,1:nP) = &
            SumFlux2_G(1:nI,nJ,1:nK,1:nP) + SumFlux2_G(1:nI,0,1:nK,1:nP)
    else
       do iP=1,nP; do k=1,nK; do i=1,nI
          if(DeltaH_DG(2,i,0,k,iP) > 0.0)then
             SumFlux2_G(i,1,k,iP) = SumFlux2_G(i,1,k,iP) + &
                  DeltaH_DG(2,i,0,k,iP)*&
                  (1.0 - CFLCoef_G(i,0,k,iP))*minmod( &
                  limiter(4,i,0,k,iP), DeltaMinusF_G(i,0,k,iP))
          end if
          if(DeltaH_DG(2,i,nJ,k,iP) < 0.0)then
             SumFlux2_G(i,nJ,k,iP)   = SumFlux2_G(i,nJ,k,iP) - &
                  DeltaH_DG(2,i,nJ,k,iP)*&
                  (1.0 - CFLCoef_G(i,nJ+1,k,iP))*minmod(  &
                  limiter(3,i,nJ+1,k,iP), DeltaMinusF_G(i,nJ+1,k,iP))
          end if
       end do; end do; end do
    end if
    if(nK > 1)then
       if(IsPeriodic_D(3))then
          SumFlux2_G(1:nI,1:nJ,1,1:nP) = &
               SumFlux2_G(1:nI,1:nJ,1,1:nP) + SumFlux2_G(1:nI,1:nJ,nK+1,1:nP)
          SumFlux2_G(1:nI,1:nJ,nK,1:nP) = &
               SumFlux2_G(1:nI,1:nJ,nK,1:nP) + SumFlux2_G(1:nI,1:nJ,0,1:nP)
       else
          do iP=1,nP; do j=1,nJ; do i=1,nI
             if(DeltaH_DG(3,i,j,0,iP) > 0.0)then
                SumFlux2_G(i,j,1,iP) = SumFlux2_G(i,j,1,iP) + &
                     DeltaH_DG(3,i,j,0,iP)*&
                     (1.0 - CFLCoef_G(i,j,0,iP))*minmod( &
                     limiter(6,i,j,0,iP), DeltaMinusF_G(i,j,0,iP))
             end if
             if(DeltaH_DG(3,i,j,nK,iP) < 0.0)then
                SumFlux2_G(i,j,nK,iP) = SumFlux2_G(i,j,nK,iP) - &
                     DeltaH_DG(3,i,j,nK,iP)*&
                     (1.0 - CFLCoef_G(i,j,nK+1,iP))*minmod( &
                     limiter(5,i,j,nK+1,iP), DeltaMinusF_G(i,j,nK+1,iP))
             end if
          end do; end do; end do
       end if
    end if
    Source_C = Source_C + TimeStep_G(1:nI,1:nJ,1:nK,1:nP)*         &
         vInv_G(1:nI,1:nJ,1:nK,1:nP)*SumFlux2_G(1:nI,1:nJ,1:nK,1:nP)
  contains
    !==========================================================================
    real function limiter(iSide,i,j,k,iP)
      integer, intent(in):: iSide,i,j,k,iP
      !------------------------------------------------------------------------
      iD_D = [i,j,k,iP] + iShift_DS(:,iSide)
      iU_D = [i,j,k,iP] - iShift_DS(:,iSide)
      limiter = betalimiter(                                                  &
           DeltaF=VDF_G(iD_D(1),iD_D(2),iD_D(3),iD_D(4))  - VDF_G(i,j,k,iP)  ,&
           UpwindDeltaF=VDF_G(i,j,k,iP) - &
           VDF_G(iU_D(1),iU_D(2),iU_D(3),iU_D(4))                            ,&
           Cfl = CFLCoef_G(i,j,k,iP)                                         ,&
           DownwindDeltaMinusF=DeltaMinusF_G(iD_D(1),iD_D(2),iD_D(3),iD_D(4)),&
           UpwindDeltaMinusF  =DeltaMinusF_G(i,j,k,iP))
    end function limiter
    !==========================================================================
  end subroutine explicit4
  !============================================================================
end module ModPoissonBracket
!==============================================================================
