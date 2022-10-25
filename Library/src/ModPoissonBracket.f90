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
  logical, parameter :: UseGroupSuperbee = .true.  .and. UseLimiter
  logical, parameter :: UseSimpleTvd = .false.     .or. .not. UseLimiter

  type RealPointer
     real, pointer :: Ptr
  end type RealPointer
  type(RealPointer) :: MajorFlux_I(6)

contains
  !============================================================================
  real function minmod(Arg1, Arg2)

    real, intent(in) :: Arg1, Arg2
    !--------------------------------------------------------------------------
    if(UseLimiter)then
       minmod = (sign(0.5, Arg1) + sign(0.5, Arg2))*min(abs(Arg1), abs(Arg2))
    else
       minmod = 0.5*(Arg1 + Arg2)
    end if
  end function minmod
  !============================================================================
  real function triple_superbee(DownwindDeltaMinusF,DeltaF, &
       UpwindDeltaF, UpwindDeltaMinusF)

    !  \delta^-f in the neighboring cells

    real, intent(in):: DownwindDeltaMinusF   ! Downwind \delta^-f
    real, intent(in):: DeltaF                ! f_j -f
    real, intent(in) :: UpwindDeltaF     ! f - f_j^\prime at the opposite face
    real, intent(in):: UpwindDeltaMinusF     ! Upwind \delta^-f
    real :: SignDeltaF, AbsDeltaF
    !--------------------------------------------------------------------------
    if(.not.UseLimiter)then
       triple_superbee = 0.5*DeltaF
       RETURN
    end if
    if(UseSimpleTvd)then
       triple_superbee = pair_superbee(DeltaF, UpwindDeltaF)
       RETURN
    end if
    if(DownwindDeltaMinusF*DeltaF<= 0.0 .or. DeltaF*UpwindDeltaMinusF<=0.0   &
         ! Nullify second order correction if at the level of machine zero
         .or.min(abs(DeltaF),abs(UpwindDeltaMinusF),abs(DownwindDeltaMinusF))&
         <1.0e-31)then
       triple_superbee = 0.0
       RETURN
    end if
    SignDeltaF = sign(1.0, DeltaF); AbsDeltaF = abs(DeltaF)
    triple_superbee = SignDeltaF*min(&
         0.50*max(AbsDeltaF, SignDeltaF*UpwindDeltaF),&
         abs(DownwindDeltaMinusF), abs(UpwindDeltaMinusF))
    if(abs(triple_superbee) < 1e-31)triple_superbee = 0.0

  end function triple_superbee
  !============================================================================
  real function group_superbee(DownwindDeltaMinusF, ReductionCoef, DeltaF, &
       UpwindDeltaF)

    !  \delta^-f in the neighboring cells

    real, intent(in):: DownwindDeltaMinusF   ! Downwind \delta^-f
    real, intent(in):: ReductionCoef
    real, intent(in):: DeltaF                ! f_j -f
    real, intent(in) :: UpwindDeltaF     ! f - f_j^\prime at the opposite face
    real :: SignDeltaF, AbsDeltaF
    !--------------------------------------------------------------------------
    ! if(DownwindDeltaMinusF*DeltaF<= 0.0   &
         ! Nullify second order correction if at the level of machine zero
    ! .or.min(abs(DeltaF),abs(DownwindDeltaMinusF)) < 1e-31)then
    if(abs(DeltaF) < 1e-31)then
       group_superbee = 0.0
       RETURN
    end if
    SignDeltaF = sign(1.0, DeltaF); AbsDeltaF = abs(DeltaF)
    group_superbee = SignDeltaF*min(&
         0.5*max(AbsDeltaF, SignDeltaF*UpwindDeltaF), ReductionCoef*AbsDeltaF)

  end function group_superbee
  !============================================================================
  real function pair_superbee(Arg1, Arg2)

    !  Used to limit only two distribution function variations
    !  \delta^-f in the neighboring cells

    real, intent(in):: Arg1, Arg2
    real :: AbsArg1, AbsArg2
    !--------------------------------------------------------------------------
    if(.not.UseLimiter)then
       pair_superbee = 0.5*(Arg1 + Arg2)
       RETURN
    end if
    if(Arg1*Arg2 <= 0.0)then
       pair_superbee = 0.0
       RETURN
    end if
    AbsArg1 = abs(Arg1); AbsArg2 = abs(Arg2)
    pair_superbee = &
         sign(min(AbsArg1, AbsArg2, 0.50*max(AbsArg1, AbsArg2)), Arg1)

  end function pair_superbee
  !============================================================================
  subroutine explicit2(nI, nJ, VDF_G, Volume_G, Source_C,    &
       Hamiltonian12_N, dHamiltonian01_FX, dHamiltonian02_FY,&
       DVolumeDt_G,                                          &
       DtIn, CFLIn, DtOut, CFLOut, DtRecommend)

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
       DtRecommend=DtRecommend)

  end subroutine explicit2
  !============================================================================
  subroutine explicit3(nI, nJ, nK, VDF_G, Volume_G, Source_C,            &
       Hamiltonian12_N, Hamiltonian13_N, Hamiltonian23_N,                &
       dHamiltonian01_FX, dHamiltonian02_FY, dHamiltonian03_FZ,          &
       DVolumeDt_G,                                                      &
       DtIn, CFLIn, DtOut, CFLOut, DtRecommend)

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

    ! Local variables
    logical :: UseTimeDependentVolume = .false. !=present(DVolumeDt_G)

    ! Inverse volume. One layer of face ghost cells is used
    real :: vInv_G(0:nI+1,0:nJ+1,1/nK:nK+1-1/nK)

    ! Loop variables:
    integer :: i, j, k
    ! Variations of VDF (one layer of ghost cell values):
    real :: DeltaMinusF_G(0:nI+1, 0:nJ+1, 0:nK+1)
    !
    ! face-centered vriations of Hamiltonian functions.
    ! one layer of ghost faces
    real :: DeltaH_FX(-1:nI+1,0:nJ+1,1/nK:nK+1-1/nK)
    real :: DeltaH_FY(0:nI+1,-1:nJ+1,1/nK:nK+1-1/nK)
    real :: DeltaH_FZ(0:nI+1,0:nJ+1,-1:nK+1)
    real, dimension(0:nI+1,0:nJ+1,1/nK:nK+1-1/nK) :: &
         SumDeltaHPlus_G, SumDeltaHMinus_G
    ! Fluxes:
    real,target :: Flux_FX(0:nI,1:nJ,1:nK)
    real,target :: Flux_FY(1:nI,0:nJ,1:nK)
    real,target :: Flux_FZ(1:nI,1:nJ,0:nK)
    real :: SumFlux_C(1:nI,1:nJ,1:nK)
    ! Local CFL number:
    real :: CFLCoef_G(0:nI+1,0:nJ+1,1/nK:nK+1-1/nK)
    ! Time step
    real :: Dt, CFL
    ! Group Limiter:
    real    :: SignDeltaMinusF, SumMajorDeltaF, VDF
    real    :: DownwindReduction_G(0:nI+1,0:nJ+1,1/nK:nK+1-1/nK)
    integer :: nMajorFlux, iFlux
    real    :: SumDeltaHPlusDeltaMinusF, SumMajorFlux, Limiter, Flux
    character(len=*), parameter:: NameSub = 'explicit3'
    !--------------------------------------------------------------------------
    if(present(DtIn))then
       Dt = DtIn
    else
       if(.not.present(CflIn))call CON_stop(&
            'Either CflIn or DtIn should be provided in '//NameSub)
    end if
    UseTimeDependentVolume = present(DVolumeDt_G)
    iKStart  = 1/nK ;  iKLast  = nK + 1 - 1/nK
    vInv_G = 1.0/Volume_G

    ! Nullify arrays:
    DeltaH_FX = 0.0; DeltaH_FY = 0.0; DeltaH_FZ = 0.0
    Flux_FX   = 0.0; Flux_FY   = 0.0; Flux_FZ   = 0.0
    SumDeltaHPlus_G = 0.0;     SumDeltaHMinus_G = 0.0
    if(UseGroupSuperbee)then
       do iFlux=1,6
          nullify(MajorFlux_I(iFlux)%Ptr)
       end do
    end if

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

    ! Now, for each cell the value of DeltaH for face in positive
    ! directions of i and j may be found in the arrays, for
    ! negative directions the should be taken with opposite sign
    ! Calculate DeltaMinusF and SumDeltaH
    do k=iKStart, iKLast; do j = 0, nJ+1; do i = 0, nI+1
       SumDeltaHPlus_G(i,j,k) = max(0.0, DeltaH_FX(i,  j,  k)) +&
                                max(0.0, DeltaH_FY(i,  j,  k)) +&
                                max(0.0,-DeltaH_FX(i-1,j,  k)) +&
                                max(0.0,-DeltaH_FY(i,j-1,  k))

       SumDeltaHMinus_G(i,j,k) =min(0.0, DeltaH_FX(i,  j,  k)) +&
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

          SumDeltaHMinus_G(i,j,k) = SumDeltaHMinus_G(i,j,k)    +&
                                    min(0.0, DeltaH_FZ(i,j,k)) +&
                                    min(0.0,-DeltaH_FZ(i,j,k-1))

          DeltaMinusF_G(i, j, k) = DeltaMinusF_G(i, j, k)      +&
               min(0.0, DeltaH_FZ(i,   j,   k))*VDF_G(i,j,k+1) +&
               min(0.0,-DeltaH_FZ(i,   j, k-1))*VDF_G(i,j,k-1)
       end if
       VDF = VDF_G(i,j,k)
       DeltaMinusF_G(i,j,k) = DeltaMinusF_G(i,j,k)           &
            /max(-SumDeltaHMinus_G(i,j,k), 1.0e-31) + VDF

       if(DeltaMinusF_G(i,j,k)==0.0)then
          DownwindReduction_G(i,j,k) = 0.0
       else
          SignDeltaMinusF = sign(1.0,DeltaMinusF_G(i,j,k))
          SumMajorDeltaF = &
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
             SumMajorDeltaF = SumMajorDeltaF + &
                  min(0.0, DeltaH_FZ(i,   j,   k))*&
                  max((VDF - VDF_G(i,j,k+1))*SignDeltaMinusF,0.0) +&
                  min(0.0,-DeltaH_FZ(i,   j, k-1))*&
                  max((VDF - VDF_G(i,j,k-1))*SignDeltaMinusF,0.0)
          end if
          DownwindReduction_G(i,j,k) = SignDeltaMinusF*DeltaMinusF_G(i,j,k)*&
               SumDeltaHMinus_G(i,j,k)/min(SumMajorDeltaF,-1.0e-31)
       end if
       if(UseTimeDependentVolume)then
          ! Local CFLs are expressed via SumDeltaHMinus
          CFLCoef_G(i,j,k) = -SumDeltaHMinus_G(i,j,k)
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

    ! By now, calculated are:
    ! 1. DeltaH_FX, DeltaH_FY, and, for nK > 1, Delta_FZ
    ! 2. Inverse volume. For time-dependent Jacobian, it is calculated at
    !    the end of time step (at t+Dt).
    ! 3. Time step, Dt.
    ! 4. CGL local, at each cell including one laayer of GC

    ! Calculate source = f(t+Dt) - f(t):
    ! First order monotone scheme
    Source_C = -CFLCoef_G(1:nI,1:nJ,1:nK)*DeltaMinusF_G(1:nI,1:nJ,1:nK)

    ! Second order correction
    if(UseGroupSuperbee)then
       SumFlux_C = 0.0
       ! Calculate Face-X fluxes.
       do k=1, nK; do j = 1, nJ; do i = 0, nI
          if(DeltaH_FX(i,j,k) > 0.0)then
             if(i > 0)then
                Flux_FX(i,j,k) = DeltaH_FX(i,j,k)*                &
                     (1.0 - CFLCoef_G(i,j,k))*group_superbee(     &
                     DownwindDeltaMinusF=DeltaMinusF_G(i+1,j,k),  &
                     ReductionCoef=DownwindReduction_G(i+1,j,k),  &
                     DeltaF=VDF_G(i+1,j,k)  - VDF_G(i  ,j,k),     &
                     UpwindDeltaF=VDF_G(i,j,k) - VDF_G(i-1,j,k))
                SumFlux_C(i,j,k) = SumFlux_C(i,j,k) +             &
                     Flux_FX(i,j,k)
             else
                Flux_FX(i,j,k) = DeltaH_FX(i,j,k)*                &
                     (1.0 - CFLCoef_G(i,j,k))*triple_superbee(    &
                     DownwindDeltaMinusF=DeltaMinusF_G(i+1,j,k),  &
                     DeltaF=VDF_G(i+1,j,k)  - VDF_G(i  ,j,k),     &
                     UpwindDeltaF=VDF_G(i,j,k) - VDF_G(i-1,j,k),  &
                     UpwindDeltaMinusF=DeltaMinusF_G(i,j,k))
             end if
          else
             if(i < nI)then
                Flux_FX(i,j,k) = DeltaH_FX(i,j,k)*                &
                     (1.0 - CFLCoef_G(i+1,j,k))*group_superbee(   &
                     DownwindDeltaMinusF=DeltaMinusF_G(i,j,k),    &
                     ReductionCoef=DownwindReduction_G(i,j,k),    &
                     DeltaF=VDF_G(i,j,k)  - VDF_G(i+1,j,k),       &
                     UpwindDeltaF=VDF_G(i+1,j,k) - VDF_G(i+2,j,k))
                SumFlux_C(i+1,j,k) = SumFlux_C(i+1,j,k) -         &
                     Flux_FX(i,j,k)
             else
                Flux_FX(i,j,k) = DeltaH_FX(i,j,k)*                &
                     (1.0 - CFLCoef_G(i+1,j,k))*triple_superbee(  &
                     DownwindDeltaMinusF=DeltaMinusF_G(i,j,k),    &
                     DeltaF=VDF_G(i,j,k)  - VDF_G(i+1,j,k),       &
                     UpwindDeltaF=VDF_G(i+1,j,k) - VDF_G(i+2,j,k),&
                     UpwindDeltaMinusF=DeltaMinusF_G(i+1  ,j,k))
             end if
          end if
       end do; end do; end do
       ! Calculate Face-Y fluxes.
       do k=1, nK; do j = 0, nJ; do i = 1, nI
          if(DeltaH_FY(i,j,k) > 0.0)then
             if(j > 0)then
                Flux_FY(i,j,k) = DeltaH_FY(i,j,k)*                &
                     (1.0 - CFLCoef_G(i,j,k))*group_superbee(     &
                     DownwindDeltaMinusF=DeltaMinusF_G(i,j+1,k),  &
                     ReductionCoef=DownwindReduction_G(i,j+1,k),  &
                     DeltaF=VDF_G(i,j+1,k)  - VDF_G(i,j  ,k),     &
                     UpwindDeltaF=VDF_G(i,j  ,k) - VDF_G(i,j-1,k))
                SumFlux_C(i,j,k) = SumFlux_C(i,j,k) +             &
                     Flux_FY(i,j,k)
             else
                Flux_FY(i,j,k) = DeltaH_FY(i,j,k)*                &
                     (1.0 - CFLCoef_G(i,j,k))*triple_superbee(    &
                     DownwindDeltaMinusF=DeltaMinusF_G(i,j+1,k),  &
                     DeltaF=VDF_G(i,j+1,k)  - VDF_G(i,j  ,k),     &
                     UpwindDeltaF=VDF_G(i,j  ,k) - VDF_G(i,j-1,k),&
                     UpwindDeltaMinusF=DeltaMinusF_G(i,j,k))
             end if
          else
             if(j < nJ)then
                Flux_FY(i,j,k) = DeltaH_FY(i,j,k)*                &
                     (1.0 - CFLCoef_G(i,j+1,k))*group_superbee(   &
                     DownwindDeltaMinusF=DeltaMinusF_G(i,j,k),    &
                     ReductionCoef=DownwindReduction_G(i,j,k),    &
                     DeltaF=VDF_G(i,j,k)  - VDF_G(i,j+1 ,k),      &
                     UpwindDeltaF=VDF_G(i,j+1,k) - VDF_G(i,j+2,k))
                SumFlux_C(i,j+1,k) = SumFlux_C(i,j+1,k) -         &
                     Flux_FY(i,j,k)
             else
                Flux_FY(i,j,k) = DeltaH_FY(i,j,k)*                &
                     (1.0 - CFLCoef_G(i,j+1,k))*triple_superbee(  &
                     DownwindDeltaMinusF=DeltaMinusF_G(i,j,k),    &
                     DeltaF=VDF_G(i,j,k)  - VDF_G(i,j+1 ,k),      &
                     UpwindDeltaF=VDF_G(i,j+1,k) - VDF_G(i,j+2,k),&
                     UpwindDeltaMinusF=DeltaMinusF_G(i,j+1,k))
             end if
          end if
       end do; end do; end do
       if(nK>1)then
          ! Calculate Face-Z fluxes.
          do k=0, nK; do j = 1, nJ; do i = 1, nI
             if(DeltaH_FZ(i,j,k) > 0.0)then
                if(k > 0)then
                   Flux_FZ(i,j,k) = DeltaH_FZ(i,j,k)*                &
                        (1.0 - CFLCoef_G(i,j,k))*group_superbee(     &
                        DownwindDeltaMinusF=DeltaMinusF_G(i,j,k+1),  &
                        ReductionCoef=DownwindReduction_G(i,j,k+1),  &
                        DeltaF=VDF_G(i,j,k+1)  - VDF_G(i,j  ,k),     &
                        UpwindDeltaF=VDF_G(i,j  ,k) - VDF_G(i,j,k-1))
                   SumFlux_C(i,j,k) = SumFlux_C(i,j,k) +             &
                        Flux_FZ(i,j,k)
                else
                   Flux_FZ(i,j,k) = DeltaH_FZ(i,j,k)*                &
                        (1.0 - CFLCoef_G(i,j,k))*triple_superbee(    &
                        DownwindDeltaMinusF=DeltaMinusF_G(i,j,k+1),  &
                        DeltaF=VDF_G(i,j,k+1)  - VDF_G(i,j  ,k),     &
                        UpwindDeltaF=VDF_G(i,j  ,k) - VDF_G(i,j,k-1),&
                        UpwindDeltaMinusF=DeltaMinusF_G(i,j  ,k))
                end if
             else
                if(k < nK)then
                   Flux_FZ(i,j,k) = DeltaH_FZ(i,j,k)*                &
                        (1.0 - CFLCoef_G(i,j,k+1))*group_superbee(   &
                        DownwindDeltaMinusF=DeltaMinusF_G(i,j,k),    &
                        ReductionCoef=DownwindReduction_G(i,j,k),    &
                        DeltaF=VDF_G(i,j,k)  - VDF_G(i,j ,k+1),      &
                        UpwindDeltaF=VDF_G(i,j,k+1) - VDF_G(i,j,k+2))
                   SumFlux_C(i,j,k+1) = SumFlux_C(i,j,k+1) -         &
                        Flux_FZ(i,j,k)
                else
                   Flux_FZ(i,j,k) = DeltaH_FZ(i,j,k)*                &
                        (1.0 - CFLCoef_G(i,j,k+1))*triple_superbee(  &
                        DownwindDeltaMinusF=DeltaMinusF_G(i,j,k),    &
                        DeltaF=VDF_G(i,j,k)  - VDF_G(i,j ,k+1),      &
                        UpwindDeltaF=VDF_G(i,j,k+1) - VDF_G(i,j,k+2),&
                        UpwindDeltaMinusF=DeltaMinusF_G(i,j,k+1))
                end if
             end if
          end do; end do; end do
       end if
       do k=1, nK; do j = 1, nJ; do i = 1, nI
          ! Compare sum of fluxes with \delta^-f
          SumDeltaHPlusDeltaMinusF = DeltaMinusF_G(i,j,k)*&
               SumDeltaHPlus_G(i,j,k)*(1.0 - CFLCoef_G(i,j,k))
          if(SumFlux_C(i,j,k)*DeltaMinusF_G(i,j,k)>=0.0.and.&
               abs(SumFlux_C(i,j,k))<=abs(SumDeltaHPlusDeltaMinusF))&
               CYCLE
          SumDeltaHPlusDeltaMinusF = minmod(SumDeltaHPlusDeltaMinusF, &
               SumFlux_C(i,j,k) )
          nMajorFlux = 0
          SumMajorFlux = 0.0
          if(DeltaH_FX(i, j, k) > 0.0)then
             ! This is \delta H^+ face, contributing to second order flux
             Flux = Flux_FX(i,j,k)
             if(Flux*SumFlux_C(i,j,k) > 0.0)then
                ! This is the major flux, having the same sign as their sum
                SumMajorFlux  = SumMajorFlux  + Flux
                nMajorFlux = nMajorFlux + 1
                MajorFlux_I(nMajorFlux)%Ptr=>Flux_FX(i,j,k)
             end if
          end if
          if(DeltaH_FY(i,   j,   k) > 0.0)then
             Flux = Flux_FY(i,j,k)
             if(Flux*SumFlux_C(i,j,k) > 0.0)then
                ! This is the major flux, having the same sign as their sum
                SumMajorFlux  = SumMajorFlux  + Flux
                nMajorFlux = nMajorFlux + 1
                MajorFlux_I(nMajorFlux)%Ptr=>Flux_FY(i,j,k)
             end if
          end if
          if(-DeltaH_FX(i-1, j, k) > 0.0)then
             ! This is \delta H^+ face, contributing to second order flux
             Flux = -Flux_FX(i-1,j,k)
             if(Flux*SumFlux_C(i,j,k) > 0.0)then
                ! This is the major flux, having the same sign as their sum
                SumMajorFlux  = SumMajorFlux  + Flux
                nMajorFlux = nMajorFlux + 1
                MajorFlux_I(nMajorFlux)%Ptr=>Flux_FX(i-1,j,k)
             end if
          end if
          if(-DeltaH_FY(i,   j-1,   k) > 0.0)then
             ! This is \delta H^+ face, contributing to second order flux
             Flux = -Flux_FY(i,j-1,k)
             if(Flux*SumFlux_C(i,j,k) > 0.0)then
                ! This is the major flux, having the same sign as their sum
                SumMajorFlux  = SumMajorFlux  + Flux
                nMajorFlux = nMajorFlux + 1
                MajorFlux_I(nMajorFlux)%Ptr=>Flux_FY(i,j-1,k)
             end if
          end if
          if(nK>1)then
             call CON_stop('group_superbee is not implemented for nK>1')
          end if
          Limiter = 1.0 + (SumDeltaHPlusDeltaMinusF - SumFlux_C(i,j,k))&
               /SumMajorFlux
          do iFlux = 1, nMajorFlux
             MajorFlux_I(iFlux)%Ptr = MajorFlux_I(iFlux)%Ptr *&
                  Limiter
             nullify(MajorFlux_I(iFlux)%Ptr)
          end do
       end do; end do; end do
    else
       ! Calculate Face-X fluxes.
       do k=1, nK; do j = 1, nJ; do i = 0, nI
          if(DeltaH_FX(i,j,k) > 0.0)then
             Flux_FX(i,j,k) = DeltaH_FX(i,j,k)*                &
                  (1.0 - CFLCoef_G(i,j,k))*triple_superbee(    &
                  DownwindDeltaMinusF=DeltaMinusF_G(i+1,j,k),  &
                  DeltaF=VDF_G(i+1,j,k)  - VDF_G(i  ,j,k),     &
                  UpwindDeltaF=VDF_G(i,j,k) - VDF_G(i-1,j,k),  &
                  UpwindDeltaMinusF=DeltaMinusF_G(i,j,k))
          else
             Flux_FX(i,j,k) = DeltaH_FX(i,j,k)*                &
                  (1.0 - CFLCoef_G(i+1,j,k))*triple_superbee(  &
                  DownwindDeltaMinusF=DeltaMinusF_G(i,j,k),    &
                  DeltaF=VDF_G(i,j,k)  - VDF_G(i+1,j,k),       &
                  UpwindDeltaF=VDF_G(i+1,j,k) - VDF_G(i+2,j,k),&
                  UpwindDeltaMinusF=DeltaMinusF_G(i+1  ,j,k))
          end if
       end do; end do; end do
       ! Calculate Face-Y fluxes.
       do k=1, nK; do j = 0, nJ; do i = 1, nI
          if(DeltaH_FY(i,j,k) > 0.0)then
             Flux_FY(i,j,k) = DeltaH_FY(i,j,k)*                &
                  (1.0 - CFLCoef_G(i,j,k))*triple_superbee(    &
                  DownwindDeltaMinusF=DeltaMinusF_G(i,j+1,k),  &
                  DeltaF=VDF_G(i,j+1,k)  - VDF_G(i,j  ,k),     &
                  UpwindDeltaF=VDF_G(i,j  ,k) - VDF_G(i,j-1,k),&
                  UpwindDeltaMinusF=DeltaMinusF_G(i,j,k))
          else
             Flux_FY(i,j,k) = DeltaH_FY(i,j,k)*                &
                  (1.0 - CFLCoef_G(i,j+1,k))*triple_superbee(  &
                  DownwindDeltaMinusF=DeltaMinusF_G(i,j,k),    &
                  DeltaF=VDF_G(i,j,k)  - VDF_G(i,j+1 ,k),      &
                  UpwindDeltaF=VDF_G(i,j+1,k) - VDF_G(i,j+2,k),&
                  UpwindDeltaMinusF=DeltaMinusF_G(i,j+1,k))
          end if
       end do; end do; end do
       if(nK>1)then
          ! Calculate Face-Z fluxes.
          do k=0, nK; do j = 1, nJ; do i = 1, nI
             if(DeltaH_FZ(i,j,k) > 0.0)then
                Flux_FZ(i,j,k) = DeltaH_FZ(i,j,k)*                &
                     (1.0 - CFLCoef_G(i,j,k))*triple_superbee(    &
                     DownwindDeltaMinusF=DeltaMinusF_G(i,j,k+1),  &
                     DeltaF=VDF_G(i,j,k+1)  - VDF_G(i,j  ,k),     &
                     UpwindDeltaF=VDF_G(i,j  ,k) - VDF_G(i,j,k-1),&
                     UpwindDeltaMinusF=DeltaMinusF_G(i,j  ,k))
             else
                Flux_FZ(i,j,k) = DeltaH_FZ(i,j,k)*                &
                     (1.0 - CFLCoef_G(i,j,k+1))*triple_superbee(  &
                     DownwindDeltaMinusF=DeltaMinusF_G(i,j,k),    &
                     DeltaF=VDF_G(i,j,k)  - VDF_G(i,j ,k+1),      &
                     UpwindDeltaF=VDF_G(i,j,k+1) - VDF_G(i,j,k+2),&
                     UpwindDeltaMinusF=DeltaMinusF_G(i,j,k+1))
             end if
          end do; end do; end do
       end if
    end if
    if(nK==1)then
       ! Two-dimensional formulation
       Source_C(1:nI,1:nJ,1) = Source_C(1:nI,1:nJ,1) + (&
            Flux_FX(0:nI-1, 1:nJ  , 1) - Flux_FX(1:nI, 1:nJ, 1)  + &
            Flux_FY(1:nI  , 0:nJ-1, 1) - Flux_FY(1:nI, 1:nJ, 1)) * &
            Dt*vInv_G(1:nI,1:nJ,1)
    else
       Source_C = Source_C + (&
            Flux_FX(0:nI-1, 1:nJ, 1:nK) - Flux_FX(1:nI, 1:nJ, 1:nK) + &
            Flux_FY(1:nI, 0:nJ-1, 1:nK) - Flux_FY(1:nI, 1:nJ, 1:nK) + &
            Flux_FZ(1:nI, 1:nJ, 0:nK-1) - Flux_FZ(1:nI, 1:nJ, 1:nK))* &
            Dt*vInv_G(1:nI,1:nJ,1:nK)
    end if
  end subroutine explicit3
  !============================================================================
end module ModPoissonBracket
!==============================================================================
