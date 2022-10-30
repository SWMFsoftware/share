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
  end interface explicit
  !PUBLIC MEMBER FUNCTION:
  public :: explicit

  ! If DtIn results in CFL>CflMax, the time step is reduced
  real, parameter :: CflMax = 0.990
  character(LEN=*), parameter:: NameMod = 'ModPoissonBracket'
  logical, parameter :: UseLimiter = .true. ! false: switch off limiters
  real, parameter :: cTol = 1.0e-12

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
  real function limiter(ReductionCoef, DeltaF, UpwindDeltaF)
    !  \delta^-f in the neighboring cells
    real, intent(in) :: ReductionCoef
    real, intent(in) :: DeltaF                ! f_j -f
    real, intent(in) :: UpwindDeltaF     ! f - f_j^\prime at the opposite face
    real :: SignDeltaF, AbsDeltaF
    !--------------------------------------------------------------------------
    if(abs(DeltaF) <=cTol)then
       limiter = 0.0
    elseif(.not.UseLimiter)then
       limiter = 0.50*DeltaF
    else
       SignDeltaF = sign(1.0, DeltaF); AbsDeltaF = abs(DeltaF)
       limiter = SignDeltaF*min( 0.5*max(AbsDeltaF, SignDeltaF*UpwindDeltaF),&
            ReductionCoef*AbsDeltaF)
    end if
  end function limiter
  !============================================================================
  subroutine explicit2(nI, nJ, VDF_G, Volume_G, Source_C,    &
       Hamiltonian12_N, dHamiltonian01_FX, dHamiltonian02_FY,&
       DVolumeDt_G,                                          &
       DtIn, CFLIn, DtOut, CFLOut, DtRecommend, IsPeriodicIn_D)

    ! solve the contribution to the
    ! numerical flux from a single Poisson bracket,
    ! df/dq_l dH/dp_l - df/dp_l dH/dq_l

    integer, intent(in) :: nI    !# of cells along coordinate 1
    integer, intent(in) :: nJ    !# of cells along coord. 2 (momentum)

    ! Distribution function with gc. Two layers of face ghostcels
    ! and one level of corner ghost cells are used
    real, intent(in) :: VDF_G(-1:nI+2,-1:nJ+2)

    ! Hamiltonian function in nodes. One layer of ghost nodes is used
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

    real, optional, intent(in) :: DtIn, CFLIn   ! Options to set time step
    real, optional, intent(out):: DtOut, CFLOut ! Options to report time step
    real, optional, intent(out):: DtRecommend   ! Calculated for given CflIn
    logical, optional, intent(in) :: IsPeriodicIn_D(:)
    character(len=*), parameter:: NameSub = 'explicit2'
    !--------------------------------------------------------------------------
    call explicit3(nI, nJ, 1, VDF_G, Volume_G, Source_C,     &
       Hamiltonian12_N=Hamiltonian12_N,                      &
       dHamiltonian01_FX= dHamiltonian01_FX,                 &
       dHamiltonian02_FY= dHamiltonian02_FY,                 &
       DVolumeDt_G=DVolumeDt_G,                              &
       DtIn=DtIn,                                            &
       CFLIn=CFLIn,                                          &
       DtOut=DtOut,                                          &
       CFLOut=CFLOut,                                        &
       DtRecommend=DtRecommend,                              &
       IsPeriodicIn_D=IsPeriodicIn_D)
  end subroutine explicit2
  !============================================================================
  subroutine explicit3(nI, nJ, nK, VDF_G, Volume_G, Source_C,            &
       Hamiltonian12_N, Hamiltonian13_N, Hamiltonian23_N,                &
       dHamiltonian01_FX, dHamiltonian02_FY, dHamiltonian03_FZ,          &
       DVolumeDt_G,                                                      &
       DtIn, CFLIn, DtOut, CFLOut, DtRecommend, IsPeriodicIn_D)

    ! solve the contribution to the numerical flux from multiple Poisson
    ! brackets, 1,2,3 enumerate phase coordinates,  0 relating to time.

    integer, intent(in) :: nI     !# of cells along coordinate 1
    integer, intent(in) :: nJ     !# of cells along coordinate 2
    integer, intent(in) :: nK     !# of cells along coordinate 3
    integer :: iKStart , iKLast

    ! Distribution function with gc. Two layers of face ghostcels
    ! and one level of corner ghost cells are used
    real, intent(in) :: VDF_G(-1:nI+2,-1:nJ+2,&
         -1 + 2*(1/nK):nK + 2*(1 - 1/nK))

    ! Hamiltonian functions in nodes. One layer of ghost nodes is used.
    ! 1. Hamiltonian function for the Poisson bracket \{f,H_{12}}_{x,y}
    !    Node-centered at XY plane, cell-centered with respect to Z
    !    (In other words, Z-aligned-edge-centered)
    real, optional, intent(in) :: Hamiltonian12_N(-1:nI+1,-1:nJ+1,&
         1/nK:nK+1-1/nK)

    ! 2. Hamiltonian function for the Poisson bracket \{f,H_{13}}_{x,z}
    !    Node-centered at XZ plane, cell-centered with respect to Y
    !    (In other words, Y-aligned-edge-centered)
    real, optional, intent(in) :: Hamiltonian13_N(-1:nI+1,0:nJ+1,-1:nK+1)

    ! 3. Hamiltonian function for the Poisson bracket \{f,H_{23}}_{y,z}
    !    Node-centered at YZ plane, cell-centered with respect to X
    !    (In other words, X-aligned-edge-centered)
    real, optional, intent(in) :: Hamiltonian23_N(0:nI+1,-1:nJ+1,-1:nK+1)

    real, optional, intent(in) :: dHamiltonian01_FX(-1:nI+1,0:nJ+1,&
         1/nK:nK+1-1/nK)
    real, optional, intent(in) :: dHamiltonian02_FY(0:nI+1,-1:nJ+1,&
         1/nK:nK+1-1/nK)
    real, optional, intent(in) :: dHamiltonian03_FZ(0:nI+1,0:nJ+1,-1:nK+1)

    ! Total Volume. One layer of face ghost cells is used
    real, intent(in) :: Volume_G(0:nI+1,0:nJ+1,1/nK:nK+1-1/nK)

    ! If non-canonical variables are used with time-dependent Jacobian,
    ! the cell volume changes in time. Need the volume derivative
    real, optional, intent(in) :: DVolumeDt_G(0:nI+1,0:nJ+1,&
         1/nK:nK+1-1/nK)

    ! Contribution to the conservative source (flux divergence) for
    ! the Poisson Bracket:
    ! send the source_c back to the main code
    real, intent(out) :: Source_C(1:nI, 1:nJ, 1:nK)

    !OPTIONAL PARAMETERS:
    real, optional, intent(in) :: DtIn, CFLIn   ! Options to set time step
    real, optional, intent(out):: DtOut, CFLOut ! Options to report time step
    real, optional, intent(out):: DtRecommend   ! Calculated for given CflIn
    logical, optional, intent(in) :: IsPeriodicIn_D(:)
    ! Local variables
    logical :: UseTimeDependentVolume = .false. !=present(DVolumeDt_G)

    ! Inverse volume. One layer of face ghost cells is used
    real :: vInv_G(0:nI+1,0:nJ+1,1/nK:nK+1-1/nK)
    ! Boundary conddition:
    logical :: IsPeriodic_D(3) = .false.
    ! Loop variables:
    integer :: i, j, k
    ! Variations of VDF (one layer of ghost cell values):
    real,dimension(0:nI+1, 0:nJ+1, 0:nK+1) :: DeltaMinusF_G
    !
    ! face-centered vriations of Hamiltonian functions.
    ! one layer of ghost faces
    real :: DeltaH_FX(-1:nI+1,0:nJ+1,1/nK:nK+1-1/nK)
    real :: DeltaH_FY(0:nI+1,-1:nJ+1,1/nK:nK+1-1/nK)
    real :: DeltaH_FZ(0:nI+1,0:nJ+1,-1:nK+1)
    real :: SumDeltaHPlus_G(0:nI+1,0:nJ+1,1/nK:nK+1-1/nK), SumDeltaHMinus
    ! Fluxes:
    real :: Flux_FX(0:nI,1:nJ,1:nK)
    real :: Flux_FY(1:nI,0:nJ,1:nK)
    real :: Flux_FZ(1:nI,1:nJ,0:nK)
    ! Sum of \delta^+H fluxes, to be limited
    real :: SumFluxPlus_C(1:nI,1:nJ,1:nK)
    ! Sum of all second order fluxes
    real :: SumFlux2_G(0:nI+1,0:nJ+1,1/nK:nK+1-1/nK)
    ! Local CFL number:
    real :: CFLCoef_G(0:nI+1,0:nJ+1,1/nK:nK+1-1/nK), CFLLocal
    ! Time step
    real :: Dt, CFL
    ! Misc:
    ! Sum of major contributions
    real    :: SumMajor
    ! Downwind limiter:
    real    :: DownwindReduction_G(0:nI+1,0:nJ+1,1/nK:nK+1-1/nK)
    ! Misc: used to calculate downwind reduction coef:
    real    :: SignDeltaMinusF, VDF
    ! Limiter with downwind reduction coefficient:
    real    :: DeltaFLimited
    ! Upwind limitater
    real    :: DeltaMinusFLim_C(1:nI,1:nJ,1:nK)
    integer :: nMajorFlux, iFlux, iMajor_I(6), jMajor_I(6), kMajor_I(6)
    real    :: DeltaPlusFLimited, UpwindReduction, Flux, MajorFlux_I(6)
    character(len=*), parameter:: NameSub = 'explicit3'
    !--------------------------------------------------------------------------
    if(present(DtIn))then
       Dt = DtIn
    else
       if(.not.present(CflIn))call CON_stop(&
            'Either CflIn or DtIn should be provided in '//NameSub)
    end if
    IsPeriodic_D = .false.
    if(present(IsPeriodicIn_D))&
         IsPeriodic_D(1:size(IsPeriodicIn_D)) = IsPeriodicIn_D
    UseTimeDependentVolume = present(DVolumeDt_G)
    iKStart  = 1/nK ;  iKLast  = nK + 1 - 1/nK
    vInv_G = 1.0/Volume_G

    ! Nullify arrays:
    DeltaH_FX = 0.0; DeltaH_FY = 0.0; DeltaH_FZ = 0.0
    Flux_FX   = 0.0; Flux_FY   = 0.0; Flux_FZ   = 0.0
    SumDeltaHPlus_G = 0.0

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
       DeltaH_FX = DeltaH_FX(-1:nI+1, 0:nJ+1,iKStart:iKLast) + &
               Hamiltonian12_N(-1:nI+1, 0:nJ+1,iKStart:iKLast) - &
               Hamiltonian12_N(-1:nI+1,-1:nJ  ,iKStart:iKLast)
       DeltaH_FY = DeltaH_FY( 0:nI+1,-1:nJ+1,iKStart:iKLast) + &
               Hamiltonian12_N(-1:nI  ,-1:nJ+1,iKStart:iKLast) - &
               Hamiltonian12_N(0:nI+1 ,-1:nJ+1,iKStart:iKLast)
    end if
    if (present(Hamiltonian13_N)) then
       DeltaH_FX = DeltaH_FX(-1:nI+1, 0:nJ+1, 0:nK+1) + &
             Hamiltonian13_N(-1:nI+1, 0:nJ+1, 0:nK+1) - &
             Hamiltonian13_N(-1:nI+1, 0:nJ+1,-1:nK  )
       DeltaH_FZ = DeltaH_FZ( 0:nI+1, 0:nJ+1,-1:nK+1) + &
             Hamiltonian13_N(-1:nI  , 0:nJ+1,-1:nK+1) - &
             Hamiltonian13_N( 0:nI+1, 0:nJ+1,-1:nK+1)
    end if
    if (present(Hamiltonian23_N)) then
       DeltaH_FY = DeltaH_FY( 0:nI+1,-1:nJ+1, 0:nK+1) + &
             Hamiltonian23_N( 0:nI+1,-1:nJ+1, 0:nK+1) - &
             Hamiltonian23_N( 0:nI+1,-1:nJ+1,-1:nK  )
       DeltaH_FZ = DeltaH_FZ( 0:nI+1, 0:nJ+1,-1:nK+1) + &
             Hamiltonian23_N( 0:nI+1,-1:nJ  ,-1:nK+1) - &
             Hamiltonian23_N( 0:nI+1, 0:nJ+1,-1:nK+1)
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
         DeltaH_FX = DeltaH_FX + dHamiltonian01_FX
    if (present(dHamiltonian02_FY))&
         DeltaH_FY = DeltaH_FY + dHamiltonian02_FY
    if (present(dHamiltonian03_FZ))&
         DeltaH_FZ = DeltaH_FZ + dHamiltonian03_FZ
    ! Cleanup
    where(abs(DeltaH_FX)<cTol)DeltaH_FX = 0.0
    where(abs(DeltaH_FY)<cTol)DeltaH_FY = 0.0
    if(nK>1)then
       where(abs(DeltaH_FZ)<cTol)DeltaH_FZ = 0.0
    end if
    ! Now, for each cell the value of DeltaH for face in positive
    ! directions of i and j may be found in the arrays, for
    ! negative directions the should be taken with opposite sign
    ! Calculate DeltaMinusF and SumDeltaH
    do k=iKStart, iKLast; do j = 0, nJ+1; do i = 0, nI+1
       SumDeltaHPlus_G(i,j,k) = max(0.0, DeltaH_FX(i,  j,  k)) +&
                                max(0.0, DeltaH_FY(i,  j,  k)) +&
                                max(0.0,-DeltaH_FX(i-1,j,  k)) +&
                                max(0.0,-DeltaH_FY(i,j-1,  k))

       SumDeltaHMinus          =min(0.0, DeltaH_FX(i,  j,  k)) +&
                                min(0.0, DeltaH_FY(i,  j,  k)) +&
                                min(0.0,-DeltaH_FX(i-1,j,  k)) +&
                                min(0.0,-DeltaH_FY(i,j-1,  k))

       DeltaMinusF_G(i, j, k) = &
               min(0.0, DeltaH_FX(i,   j,   k))*VDF_G(i+1,j,k) +&
               min(0.0, DeltaH_FY(i,   j,   k))*VDF_G(i,j+1,k) +&
               min(0.0,-DeltaH_FX(i-1, j,   k))*VDF_G(i-1,j,k) +&
               min(0.0,-DeltaH_FY(i, j-1,   k))*VDF_G(i,j-1,k)
       if(nK>1)then
          ! Add three-dimensional effects.
          SumDeltaHPlus_G(i,j,k) = SumDeltaHPlus_G(i,j,k)      +&
                                    max(0.0, DeltaH_FZ(i,j,k)) +&
                                    max(0.0,-DeltaH_FZ(i,j,k-1))

          SumDeltaHMinus         = SumDeltaHMinus              +&
                                    min(0.0, DeltaH_FZ(i,j,k)) +&
                                    min(0.0,-DeltaH_FZ(i,j,k-1))

          DeltaMinusF_G(i, j, k) = DeltaMinusF_G(i, j, k)      +&
               min(0.0, DeltaH_FZ(i,   j,   k))*VDF_G(i,j,k+1) +&
               min(0.0,-DeltaH_FZ(i,   j, k-1))*VDF_G(i,j,k-1)
       end if
       if(abs(SumDeltaHPlus_G(i,j,k)) < cTol)SumDeltaHPlus_G(i,j,k) = 0.0
       VDF = VDF_G(i,j,k)
       if(SumDeltaHMinus==0.0)then
          DeltaMinusF_G(i,j,k) = 0.0
       else
          DeltaMinusF_G(i,j,k) = VDF - DeltaMinusF_G(i,j,k)/SumDeltaHMinus
       end if
       if(abs(DeltaMinusF_G(i,j,k))<cTol)then
          DeltaMinusF_G(i,j,k) = 0.0
          DownwindReduction_G(i,j,k) = 0.0
       else
          SignDeltaMinusF = sign(1.0,DeltaMinusF_G(i,j,k))
          SumMajor = &
               min(0.0, DeltaH_FX(i,   j,   k))*&
               max((VDF - VDF_G(i+1,j,k))*SignDeltaMinusF,0.0) +&
               min(0.0, DeltaH_FY(i,   j,   k))*&
               max((VDF -VDF_G(i,j+1,k))*SignDeltaMinusF,0.0) +&
               min(0.0,-DeltaH_FX(i-1, j,   k))*&
               max((VDF -VDF_G(i-1,j,k))*SignDeltaMinusF,0.0) +&
               min(0.0,-DeltaH_FY(i, j-1,   k))*&
               max((VDF -VDF_G(i,j-1,k))*SignDeltaMinusF,0.0)
          if(nK>1)then
             ! Add three-dimensional effects.
             SumMajor = SumMajor + &
                  min(0.0, DeltaH_FZ(i,   j,   k))*&
                  max((VDF - VDF_G(i,j,k+1))*SignDeltaMinusF,0.0) +&
                  min(0.0,-DeltaH_FZ(i,   j, k-1))*&
                  max((VDF - VDF_G(i,j,k-1))*SignDeltaMinusF,0.0)
          end if
          DownwindReduction_G(i,j,k) = SignDeltaMinusF*DeltaMinusF_G(i,j,k)*&
               SumDeltaHMinus/SumMajor
       end if
       if(UseTimeDependentVolume)then
          ! Local CFLs are expressed via SumDeltaHMinus
          CFLCoef_G(i,j,k) = -SumDeltaHMinus
       else
          CFLCoef_G(i,j,k) = vInv_G(i,j,k)*SumDeltaHPlus_G(i,j,k)
       end if
    end do; end do; end do
    ! Set CFL and time step
    if(UseTimeDependentVolume)then
       if(present(DtRecommend))then
          if(present(CFLIn))then
             CFL = CFLIn
          else
             CFL = CFLMax
          end if
          ! Solve time step from equation
          ! CFLIn = \Delta t*(-\sum\delta^-H)/(\Delta t*dV/dt + V)
          DtRecommend = CFL/maxval(vInv_G(1:nI,1:nJ,1:nK)*&
               (CFLCoef_G(1:nI,1:nJ,1:nK) - CFL*DVolumeDt_G(1:nI,1:nJ,1:nK)))
       end if
       if(present(DtIn))then

          ! Calculate the CFL factor with given Dt:
          vInv_G = 1.0/(Volume_G + Dt*DVolumeDt_G)
          CFLCoef_G = Dt*vInv_G*CFLCoef_G
          CFL = maxval(CFLCoef_G(1:nI,1:nJ,1:nK))

          ! Check if the CFL satisfies the stability criterion
          if(CFL>CFLMax)then

             ! Restore CFLCoef_G and vInv
             CFLCoef_G = CFLCoef_G/(Dt*vInv_G)
             vInv_G = 1/Volume_G

             ! Reduce the time step using equation
             ! CFLMax = \Delta t*(-\sum\delta^-H)/(\Delta t*dV/dt + V)
             CFL = CFLMax
             Dt = CFL/maxval(vInv_G(1:nI,1:nJ,1:nK)*&
                  ( CFLCoef_G(1:nI,1:nJ,1:nK) &
                  - CFL*DVolumeDt_G(1:nI,1:nJ,1:nK)))

             ! Calculate the CFL factor with given Dt:
             vInv_G = 1.0/(Volume_G + Dt*DVolumeDt_G)
             CFLCoef_G = Dt*vInv_G*CFLCoef_G
          end if
       else
          ! Solve time step from equation
          ! CFLIn = \Delta t*(-\sum\delta^-H)/(\Delta t*dV/dt + V)
          Dt = CFLIn/maxval(vInv_G(1:nI,1:nJ,1:nK)*&
               (CFLCoef_G(1:nI,1:nJ,1:nK) - CFLIn*DVolumeDt_G(1:nI,1:nJ,1:nK)))

          ! Calculate the volume at upper time level
          ! V(+\Delta t):
          vInv_G = 1.0/(Volume_G + Dt*DVolumeDt_G)
          CFLCoef_G = Dt*vInv_G*CFLCoef_G
       end if
    else
       if(.not.present(DtIn))&
            Dt = CFLIn/maxval(CFLCoef_G(1:nI,1:nJ,1:nK))
       CFLCoef_G = Dt*CFLCoef_G
    end if
    if(present(CFLOut))then

       ! The debugging version: normally, the CFL is known by now and there
       ! is no need to calculate it again
       CFLOut = maxval(CFLCoef_G(1:nI,1:nJ,1:nK))
       if(CFLOut > CFLMax*0.99 + 0.01)call CON_stop('CFL is too large')
    end if
    if(present(DtOut ))DtOut  = Dt

    ! Calculate source = f(t+Dt) - f(t):
    ! First order monotone scheme
    Source_C = -CFLCoef_G(1:nI,1:nJ,1:nK)*DeltaMinusF_G(1:nI,1:nJ,1:nK)

    ! Second order correction
    SumFlux2_G = 0.0; SumFluxPlus_C = 0.0
    DeltaMinusFLim_C = DeltaMinusF_G(1:nI,1:nJ,1:nK) 
    ! Calculate Face-X fluxes.
    do k=1, nK; do j = 1, nJ; do i = 0, nI
       if(DeltaH_FX(i,j,k) > 0.0)then
          DeltaFLimited = limiter(        &
               ReductionCoef=DownwindReduction_G(i+1,j,k),  &
               DeltaF=VDF_G(i+1,j,k)  - VDF_G(i  ,j,k),     &
               UpwindDeltaF=VDF_G(i,j,k) - VDF_G(i-1,j,k))
          if(i > 0)then
             Flux_FX(i,j,k) = DeltaH_FX(i,j,k)*DeltaFLimited
             SumFluxPlus_C(i,j,k) = SumFluxPlus_C(i,j,k) + Flux_FX(i,j,k)
             DeltaMinusFLim_C(i,j,k) = minmod(DeltaMinusFLim_C(i,j,k),&
                  DeltaMinusF_G(i+1,j,k))
          elseif(.not.IsPeriodic_D(1))then
             ! Limit both spatial and temporal corrections
             DeltaFLimited = minmod(DeltaFLimited, DeltaMinusF_G(i,j,k))
             SumFlux2_G(i+1,j,k) = SumFlux2_G(i+1,j,k) + &
                  DeltaH_FX(i,j,k)*(1.0 - CFLCoef_G(i,j,k))*&
                  minmod(DeltaFLimited,DeltaMinusF_G(i+1,j,k))
          end if
       elseif(DeltaH_FX(i,j,k) < 0.0 )then
          DeltaFLimited = limiter( &
               ReductionCoef=DownwindReduction_G(i,j,k),    &
               DeltaF=VDF_G(i,j,k)  - VDF_G(i+1,j,k),       &
               UpwindDeltaF=VDF_G(i+1,j,k) - VDF_G(i+2,j,k))
          if(i < nI)then
             Flux_FX(i,j,k) = DeltaH_FX(i,j,k)*DeltaFLimited
             SumFluxPlus_C(i+1,j,k) = SumFluxPlus_C(i+1,j,k) - Flux_FX(i,j,k)
             DeltaMinusFLim_C(i+1,j,k) = minmod(DeltaMinusFLim_C(i+1,j,k),&
                  DeltaMinusF_G(i,j,k))
          elseif(.not.IsPeriodic_D(1))then
             ! Limit both spatial and temporal corrections
             DeltaFLimited = minmod(DeltaFLimited, DeltaMinusF_G(i+1,j,k))
             ! Extra limitation for temporal correction
             SumFlux2_G(i,j,k) = SumFlux2_G(i,j,k) - DeltaH_FX(i,j,k)*&
                  (1.0 - CFLCoef_G(i+1,j,k))*minmod(DeltaFLimited,&
                  DeltaMinusF_G(i,j,k))
          end if
       end if
    end do; end do; end do
    ! Calculate Face-Y fluxes.
    do k=1, nK; do j = 0, nJ; do i = 1, nI
       if(DeltaH_FY(i,j,k) > 0.0)then
          DeltaFLimited = limiter( &
               ReductionCoef=DownwindReduction_G(i,j+1,k),  &
               DeltaF=VDF_G(i,j+1,k)  - VDF_G(i,j  ,k),     &
               UpwindDeltaF=VDF_G(i,j  ,k) - VDF_G(i,j-1,k))
          if(j > 0)then
             Flux_FY(i,j,k) = DeltaH_FY(i,j,k)*DeltaFLimited
             SumFluxPlus_C(i,j,k) = SumFluxPlus_C(i,j,k) + Flux_FY(i,j,k)
             DeltaMinusFLim_C(i,j,k) = minmod(DeltaMinusFLim_C(i,j,k),&
                  DeltaMinusF_G(i,j+1,k))
          elseif(.not.IsPeriodic_D(2))then
             ! Limit both spatial and temporal corrections
             DeltaFLimited = minmod(DeltaFLimited, DeltaMinusF_G(i,j,k))
             ! Extra limitation for temporal correction
             SumFlux2_G(i,j+1,k) = SumFlux2_G(i,j+1,k) + DeltaH_FY(i,j,k)*&
                  (1.0 - CFLCoef_G(i,j,k))*minmod( &
                  DeltaFLimited, DeltaMinusF_G(i,j+1,k))
          end if
       elseif(DeltaH_FY(i,j,k) < 0.0)then
          DeltaFLimited = limiter( &
               ReductionCoef=DownwindReduction_G(i,j,k),    &
               DeltaF=VDF_G(i,j,k)  - VDF_G(i,j+1 ,k),      &
               UpwindDeltaF=VDF_G(i,j+1,k) - VDF_G(i,j+2,k))
          if(j < nJ)then
             Flux_FY(i,j,k) = DeltaH_FY(i,j,k)*DeltaFLimited
             SumFluxPlus_C(i,j+1,k) = SumFluxPlus_C(i,j+1,k) - Flux_FY(i,j,k)
             DeltaMinusFLim_C(i,j+1,k) = minmod(DeltaMinusFLim_C(i,j+1,k),&
                  DeltaMinusF_G(i,j,k))
          elseif(.not.IsPeriodic_D(2))then
             ! Limit both spatial and temporal corrections
             DeltaFLimited = minmod(DeltaFLimited, DeltaMinusF_G(i,j+1,k))
             ! Extra limitation for temporal correction
             SumFlux2_G(i,j,k) = SumFlux2_G(i,j,k) - DeltaH_FY(i,j,k)*&
                  (1.0 - CFLCoef_G(i,j+1,k))*minmod(  &
                  DeltaFLimited, DeltaMinusF_G(i,j,k))
          end if
       end if
    end do; end do; end do
    if(nK>1)then
       ! Calculate Face-Z fluxes.
       do k=0, nK; do j = 1, nJ; do i = 1, nI
          if(DeltaH_FZ(i,j,k) > 0.0)then
             DeltaFLimited = limiter( &
                  ReductionCoef=DownwindReduction_G(i,j,k+1),  &
                  DeltaF=VDF_G(i,j,k+1)  - VDF_G(i,j  ,k),     &
                  UpwindDeltaF=VDF_G(i,j  ,k) - VDF_G(i,j,k-1))
             if(k > 0)then
                Flux_FZ(i,j,k) = DeltaH_FZ(i,j,k)*DeltaFLimited
                SumFluxPlus_C(i,j,k) = SumFluxPlus_C(i,j,k) + Flux_FZ(i,j,k)
                DeltaMinusFLim_C(i,j,k) = minmod(DeltaMinusFLim_C(i,j,k),&
                     DeltaMinusF_G(i,j,k+1))
             elseif(.not.IsPeriodic_D(3))then
                ! Limit both spatial and temporal corrections
                DeltaFLimited = minmod(DeltaFLimited, DeltaMinusF_G(i,j,k))
                ! Extra limitation for temporal correction
                SumFlux2_G(i,j,k+1) = SumFlux2_G(i,j,k+1) + DeltaH_FZ(i,j,k)*&
                     (1.0 - CFLCoef_G(i,j,k))*minmod( &
                     DeltaFLimited, DeltaMinusF_G(i,j,k+1))
             end if
          elseif(DeltaH_FZ(i,j,k) < 0.0)then
             DeltaFLimited = limiter( &
                  ReductionCoef=DownwindReduction_G(i,j,k),    &
                  DeltaF=VDF_G(i,j,k)  - VDF_G(i,j ,k+1),      &
                  UpwindDeltaF=VDF_G(i,j,k+1) - VDF_G(i,j,k+2))
             if(k < nK)then
                Flux_FZ(i,j,k) = DeltaH_FZ(i,j,k)*DeltaFLimited
                SumFluxPlus_C(i,j,k+1) = SumFluxPlus_C(i,j,k+1) - &
                     Flux_FZ(i,j,k)
                DeltaMinusFLim_C(i,j,k+1) = minmod(DeltaMinusFLim_C(i,j,k+1),&
                     DeltaMinusF_G(i,j,k))
             elseif(.not.IsPeriodic_D(3))then
                ! Limit both spatial and temporal corrections
                DeltaFLimited = minmod(DeltaFLimited, &
                     DeltaMinusF_G(i,j,k+1))
                ! Extra limitation for temporal correction
                SumFlux2_G(i,j,k) = SumFlux2_G(i,j,k) - DeltaH_FZ(i,j,k)*&
                     (1.0 - CFLCoef_G(i,j,k+1))*minmod( &
                     DeltaFLimited ,DeltaMinusF_G(i,j,k))
             end if
          end if
       end do; end do; end do
    end if
    do k=1, nK; do j = 1, nJ; do i = 1, nI
       if(SumDeltaHPlus_G(i,j,k)==0.0)CYCLE
       CFLLocal = CFLCoef_G(i,j,k)
       ! Calculate \delta^+\Psi(f_j-f)
       DeltaPlusFLimited  = SumFluxPlus_C(i,j,k)/SumDeltaHPlus_G(i,j,k)
       if( (DeltaPlusFLimited*DeltaMinusFLim_C(i,j,k)>=0.0.and.&
            abs(DeltaPlusFLimited)<=abs(DeltaMinusFLim_C(i,j,k)))&
            .or.(.not.UseLimiter))then
          ! There is no need to limit DeltaPlus with DeltaMinus_G(i,j,k)
          if(DeltaH_FX(i, j, k) > 0.0)&
               ! This is \delta H^+ face, contributing to second order flux
               SumFlux2_G(i+1,j,k) = SumFlux2_G(i+1,j,k) + Flux_FX(i,j,k)&
               - DeltaH_FX(i,j,k)*CFLLocal*DeltaPlusFLimited
          if(DeltaH_FY(i,j,k) > 0.0)&
               ! This is \delta H^+ face, contributing to second order flux
               SumFlux2_G(i,j+1,k) = SumFlux2_G(i,j+1,k) + Flux_FY(i,j,k)&
               - DeltaH_FY(i,j,k)*CFLLocal*DeltaPlusFLimited
          if(DeltaH_FX(i-1,j,k) < 0.0)&
               ! This is \delta H^+ face, contributing to second order flux
               SumFlux2_G(i-1,j,k) = SumFlux2_G(i-1,j,k) - Flux_FX(i-1,j,k)&
               + DeltaH_FX(i-1,j,k)* CFLLocal*DeltaPlusFLimited
          if(DeltaH_FY(i,j-1,k) < 0.0)&
               ! This is \delta H^+ face, contributing to second order flux
               SumFlux2_G(i,j-1,k) = SumFlux2_G(i,j-1,k) - Flux_FY(i,j-1,k)&
               + DeltaH_FY(i,j-1,k)*CFLLocal*DeltaPlusFLimited
          if(nK>1)then
             if(DeltaH_FZ(i,j,k) > 0.0)&
                  ! This is \delta H^+ face, contributing to second order flux
                  SumFlux2_G(i,j,k+1) = SumFlux2_G(i,j,k+1) + Flux_FZ(i,j,k)&
                  - DeltaH_FZ(i,j,k)*CFLLocal*DeltaPlusFLimited
             if(DeltaH_FZ(i,j,k-1) < 0.0)&
                  ! This is \delta H^+ face, contributing to second order flux
                  SumFlux2_G(i,j,k-1) = SumFlux2_G(i,j,k-1) - Flux_FX(i,j,k-1)&
                  + DeltaH_FZ(i,j,k-1)*CFLLocal*DeltaPlusFLimited
          end if
          SumFlux2_G(i,j,k) = SumFlux2_G(i,j,k) - &
               (1 - CFLLocal)*SumFluxPlus_C(i,j,k) 
       else
          ! Limit DeltaPlus
          DeltaPlusFLimited = minmod(DeltaPlusFLimited,DeltaMinusFLim_C(i,j,k))
          ! Att the total limited flux in the given cell
          SumFlux2_G(i,j,k) = SumFlux2_G(i,j,k) - &
               (1.0 - CFLLocal)*DeltaPlusFLimited*SumDeltaHPlus_G(i,j,k)
          nMajorFlux = 0; SumMajor = 0.0
          if(DeltaH_FX(i, j, k) > 0.0)then
             ! This is \delta H^+ face, contributing to second order flux
             SumFlux2_G(i+1,j,k) = SumFlux2_G(i+1,j,k)   &
                  - DeltaH_FX(i,j,k)*CFLLocal*DeltaPlusFLimited
             Flux = Flux_FX(i,j,k)
             if(Flux*SumFluxPlus_C(i,j,k) > 0.0)then
                ! This is the major flux, having the same sign as their sum
                SumMajor  = SumMajor  + Flux
                nMajorFlux = nMajorFlux + 1
                MajorFlux_I(nMajorFlux) = Flux
                iMajor_I(nMajorFlux) = i+1
                jMajor_I(nMajorFlux) = j
                kMajor_I(nMajorFlux) = k
             else
                SumFlux2_G(i+1,j,k) = SumFlux2_G(i+1,j,k) + Flux
             end if
          end if
          if(DeltaH_FY(i,   j,   k) > 0.0)then
             SumFlux2_G(i,j+1,k) = SumFlux2_G(i,j+1,k)  &
                  - DeltaH_FY(i,j,k)*CFLLocal*DeltaPlusFLimited
             Flux = Flux_FY(i,j,k)
             if(Flux*SumFluxPlus_C(i,j,k) > 0.0)then
                ! This is the major flux, having the same sign as their sum
                SumMajor  = SumMajor  + Flux
                nMajorFlux = nMajorFlux + 1
                MajorFlux_I(nMajorFlux) = Flux
                iMajor_I(nMajorFlux) = i
                jMajor_I(nMajorFlux) = j+1
                kMajor_I(nMajorFlux) = k
             else
                SumFlux2_G(i,j+1,k) = SumFlux2_G(i,j+1,k) + Flux
             end if
          end if
          if(DeltaH_FX(i-1, j, k) < 0.0)then
             ! This is \delta H^+ face, contributing to second order flux
             SumFlux2_G(i-1,j,k) = SumFlux2_G(i-1,j,k) &
                  + DeltaH_FX(i-1,j,k)*CFLLocal*DeltaPlusFLimited
             Flux = -Flux_FX(i-1,j,k)
             if(Flux*SumFluxPlus_C(i,j,k) > 0.0)then
                ! This is the major flux, having the same sign as their sum
                SumMajor  = SumMajor  + Flux
                nMajorFlux = nMajorFlux + 1
                MajorFlux_I(nMajorFlux) = Flux
                iMajor_I(nMajorFlux) = i-1
                jMajor_I(nMajorFlux) = j
                kMajor_I(nMajorFlux) = k
             else
                SumFlux2_G(i-1,j,k) = SumFlux2_G(i-1,j,k) + Flux
             end if
          end if
          if(DeltaH_FY(i,   j-1,   k) < 0.0)then
             ! This is \delta H^+ face, contributing to second order flux
             SumFlux2_G(i,j-1,k) = SumFlux2_G(i,j-1,k)  &
                  + DeltaH_FY(i,j-1,k)*CFLLocal*DeltaPlusFLimited
             Flux = -Flux_FY(i,j-1,k)
             if(Flux*SumFluxPlus_C(i,j,k) > 0.0)then
                ! This is the major flux, having the same sign as their sum
                SumMajor  = SumMajor  + Flux
                nMajorFlux = nMajorFlux + 1
                MajorFlux_I(nMajorFlux) = Flux
                iMajor_I(nMajorFlux) = i
                jMajor_I(nMajorFlux) = j-1
                kMajor_I(nMajorFlux) = k
             else
                SumFlux2_G(i,j-1,k) = SumFlux2_G(i,j-1,k) + Flux
             end if
          end if
          if(nK>1)then
             if(DeltaH_FZ(i,j,k) > 0.0)then
                ! This is \delta H^+ face, contributing to second order flux
                SumFlux2_G(i,j,k+1) = SumFlux2_G(i,j,k+1) &
                     - DeltaH_FZ(i,j,k)*CFLLocal*DeltaPlusFLimited
                Flux = Flux_FZ(i,j,k)
                if(Flux*SumFluxPlus_C(i,j,k) > 0.0)then
                   ! This is the major flux, having the same sign as their sum
                   SumMajor  = SumMajor  + Flux
                   nMajorFlux = nMajorFlux + 1
                   MajorFlux_I(nMajorFlux) = Flux
                   iMajor_I(nMajorFlux) = i
                   jMajor_I(nMajorFlux) = j
                   kMajor_I(nMajorFlux) = k+1
                else
                   SumFlux2_G(i,j,k+1) = SumFlux2_G(i,j,k+1) + Flux
                end if
             end if
             if(DeltaH_FZ(i,j,k-1) < 0.0)then
                ! This is \delta H^+ face, contributing to second order flux
                SumFlux2_G(i,j,k-1) = SumFlux2_G(i,j,k-1)  &
                     + DeltaH_FZ(i,j,k-1)*CFLLocal*DeltaPlusFLimited
                Flux = -Flux_FZ(i,j,k-1)
                if(Flux*SumFluxPlus_C(i,j,k) > 0.0)then
                   ! This is the major flux, having the same sign as their sum
                   SumMajor  = SumMajor  + Flux
                   nMajorFlux = nMajorFlux + 1
                   MajorFlux_I(nMajorFlux) = Flux
                   iMajor_I(nMajorFlux) = i
                   jMajor_I(nMajorFlux) = j
                   kMajor_I(nMajorFlux) = k-1
                else
                   SumFlux2_G(i,j,k-1) = SumFlux2_G(i,j,k-1) + Flux
                end if
             end if
          end if
          UpwindReduction = 1.0 + (DeltaPlusFLimited*SumDeltaHPlus_G(i,j,k) - &
               SumFluxPlus_C(i,j,k))/SumMajor
          do iFlux = 1, nMajorFlux
             SumFlux2_G(iMajor_I(iFlux),jMajor_I(iFlux),kMajor_I(iFlux)) =   &
                  SumFlux2_G(iMajor_I(iFlux),jMajor_I(iFlux),kMajor_I(iFlux))&
                  + MajorFlux_I(iFlux)*UpwindReduction
          end do
       end if
    end do; end do; end do
    if(IsPeriodic_D(1))then
       SumFlux2_G(1,1:nJ,1:nK) = &
            SumFlux2_G(1,1:nJ,1:nK) + SumFlux2_G(nI+1,1:nJ,1:nK)
       SumFlux2_G(nI,1:nJ,1:nK) = &
            SumFlux2_G(nI,1:nJ,1:nK) + SumFlux2_G(0,1:nJ,1:nK)
    end if
    if(IsPeriodic_D(2))then
       SumFlux2_G(1:nI,1,1:nK) = &
            SumFlux2_G(1:nI,1,1:nK) + SumFlux2_G(1:nI,nJ+1,1:nK)
       SumFlux2_G(1:nI,nJ,1:nK) = &
            SumFlux2_G(1:nI,nJ,1:nK) + SumFlux2_G(1:nI,0,1:nK)
    end if
    if(IsPeriodic_D(3))then
       SumFlux2_G(1:nI,1:nJ,1) = &
            SumFlux2_G(1:nI,1:nJ,1) + SumFlux2_G(1:nI,1:nJ,nK+1)
       SumFlux2_G(1:nI,1:nJ,nK) = &
            SumFlux2_G(1:nI,1:nJ,nK) + SumFlux2_G(1:nI,1:nJ,0)
    end if
    Source_C = Source_C + Dt*vInv_G(1:nI,1:nJ,1:nK)*SumFlux2_G(1:nI,1:nJ,1:nK)
  end subroutine explicit3
  !============================================================================
end module ModPoissonBracket
!==============================================================================
