!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!BOP -------------------------------------------------------------------
!
!MODULE: ModPoissonBracket - a numerical flux from the Poisson bracket
!  
!
!DESCRIPTION:
!
! For a contribution to the Liuville (kinetic) equation:
! 
! df/dt + \sum{df/dq_l dH/dp_l - df/dp_l dH/dq_l}=0,
!
! in which f in the (unknown) distribution function, 
! q_l are generalized coordinates,
! p_l are canonically conjugate generalized momenta,
! H is the (known!) Hamiltonian function,
! 
! the methods in the module solve the contribution to the
! numerical flux from a single Poisson bracket,
! df/dq_l dH/dp_l - df/dp_l dH/dq_l
!
!INTERFACE:

module ModPoissonBracket

  !USES:
  use ModMpi
  use ModUtilities,    ONLY: CON_stop
  use ModLinearSolver, ONLY: bicgstab, prehepta, &
                             upper_hepta_scalar, lower_hepta_scalar, & 
                             get_precond_matrix
  !\     
  ! bicgstab             BiCGSTAB iterative solver
  ! prehepta             LU preconditioner for up to hepta block-diagonal 
  ! upper_hepta_scalar   multiply with upper scalar triangular matrix
  ! lower_hepta_scalar   multiply with lower scalar triangular matrix
  ! get_precond_matrix   get preconditioner matrix
  implicit none
  SAVE
  private ! except

  !PUBLIC MEMBER FUNCTION:
  public :: explicit
contains
  ! A choice of limiter functions
  real function pair_superbee(Arg1, Arg2)
    real, intent(in):: Arg1, Arg2
    real :: AbsArg1, AbsArg2
    !------------
    if(Arg1*Arg2 .le. 0.0)then
       pair_superbee = 0.0
       RETURN
    end if
    AbsArg1 = abs(Arg1)
    AbsArg2 = abs(Arg2)
    pair_superbee = sign(min(2*min(AbsArg1,AbsArg2),max(AbsArg1,AbsArg2)),&
         Arg1)
  end function pair_superbee
  !===============
  real function triple_superbee(Arg1, Arg2, Arg3)
    real, intent(in):: Arg1, Arg2, Arg3
    real :: AbsArg_I(3)
    !------------
    if(Arg1*Arg2 .le. 0.0 .or. Arg2*Arg3 .le. 0.0)then
       triple_superbee = 0.0
       RETURN
    end if
    AbsArg_I(1) = abs(Arg1)
    AbsArg_I(2) = abs(Arg2)
    AbsArg_I(3) = abs(Arg3)
    triple_superbee = sign(min(2*minval(AbsArg_I),maxval(AbsArg_I)),&
         Arg1)
  end function triple_superbee
  !==========================
  subroutine explicit(nQ, nP, VDF_G, Hamiltonian_N,   &
       Volume_G, Source_C, PreBracketFactor_G, DtIn, CFLIn, DtOut, CFLOut)
    !\
    ! solve the contribution to the
    ! numerical flux from a single Poisson bracket,
    ! df/dq_l dH/dp_l - df/dp_l dH/dq_l
    !/
    integer, intent(in) :: nQ           !number of cells along coordinate axis 
    integer, intent(in) :: nP           !number of cells along momentum   axis
    !\
    !Distribution function with gc. Two layers of face ghostcels and one level
    !of corner ghost cells are used
    !/
    real, intent(in) :: VDF_G(-1:nQ+2,-1:nP+2) 
    !\ 
    ! Hamiltonian function in nodes. One layer of ghost nodes is used
    !/
    real, intent(in) :: Hamiltonian_N(-1:nQ+1,-1:nP+1)  
    !\
    ! Inverse volume. One layer of face ghost cells is used
    !/  
    real, intent(in) :: Volume_G(0:nQ+1,0:nP+1)
    !\
    ! Inverse volume. One layer of face ghost cells is used
    !/
    real             ::vInv_G(0:nQ+1,0:nP+1)
    !\
    ! Contribution to the conservative source (flux divergence) for 
    ! the Poisson Bracket:
    !
    ! df/dt = sum(Source_C), sum is taken over all Poisson brackets
    !/                  
    real, intent(out):: Source_C(1:nQ, 1:nP)  
    !\
    !OPTIONAL:
    !/
    !\
    ! For the case of equation df/dt + Factor*{f, H} = 0
    !/
    real, optional, intent(in) :: PreBracketFactor_G(0:nQ+1,0:nP+1) 
    real, optional, intent(in) :: DtIn, CFLIn   !Options to set time step
    real, optional, intent(out):: DtOut, CFLOut !Options to report time step
    !\
    ! Local variables
    !/
    !Loop variables:
    integer :: iQ, iP
    
    ! Variations of VDF (one layer of ghost cell values):
    real :: DeltaMinusF_G(0:nQ+1, 0:nP+1)
    !\
    ! Variations of Hamiltonian functions (one layer of ghost faces, 
    ! for two directions:
    !/
    real :: DeltaH_FX(-1:nQ+1,0:nP+1), DeltaH_FY(0:nQ+1,-1:nP+1)
    real :: SumDeltaHPlus_G(0:nQ+1, 0:nP+1)
    !Fluxes
    real :: Flux_FX(0:nQ,1:nP), Flux_FY(1:nQ,0:nP)
    
    ! Local CFL number:
    real :: CFLCoef_G(0:nQ+1,0:nP+1)
    
    !Time step
    real :: Dt, CFL
    character(len=*), parameter:: NameSub = 'calc_poisson_bracket'
    !---------
    if(present(DtIn))then
       Dt = DtIn
    else
       if(.not.present(CflIn))call CON_stop(&
            'Either CflIn or DtIn should be provided in '//NameSub)
    end if
    VInv_G = 1/Volume_G
    !\
    ! Calculate DeltaH:
    !/
    DeltaH_FX(-1:nQ+1, 0:nP+1) = Hamiltonian_N(-1:nQ+1,   0:nP+1) - &
                            Hamiltonian_N(-1:nQ+1, -1:nP  )
    DeltaH_FY(0 :nQ+1,-1:nP+1) = Hamiltonian_N(-1:nQ,    -1:nP+1) - &
                            Hamiltonian_N(0:nQ+1 , -1:nP+1)
    ! Now, for each cell the value of DeltaH for face in positive 
    ! directions of iQ and iP may be found in the arrays, for
    ! negative directions the should be taken with opposite sign

    !\
    ! Calculate DeltaMinusF and SumDeltaH
    !/   
    do iP = 0, nP+1; do iQ = 0, nQ+1
       SumDeltaHPlus_G(iQ, iP) = max(0.0, DeltaH_FX(iQ,   iP)) +&
                                 max(0.0, DeltaH_FY(iQ,   iP)) +&
                                 max(0.0,-DeltaH_FX(iQ-1, iP)) +&
                                 max(0.0,-DeltaH_FY(iQ, iP-1))
       DeltaMinusF_G(iQ, iP) = (&
           min(0.0, DeltaH_FX(iQ,   iP))*VDF_G(iQ+1,iP) +&
           min(0.0, DeltaH_FY(iQ,   iP))*VDF_G(iQ,iP+1) +&
           min(0.0,-DeltaH_FX(iQ-1, iP))*VDF_G(iQ-1,iP) +&
           min(0.0,-DeltaH_FY(iQ, iP-1))*VDF_G(iQ,iP-1)  &
           )/SumDeltaHPlus_G(iQ, iP) + VDF_G(iQ,iP)
       CFLCoef_G(iQ, iP) = vInv_G(iQ, iP)*SumDeltaHPlus_G(iQ, iP)
    end do; end do
    !\
    ! Set Dt and construct array Dt*SumDeltaH/Volume_C
    !/
    if(present(DtIn))then
       CFLCoef_G = Dt*CFLCoef_G
       CFL = maxval(CFLCoef_G)
    else
       CFL = CFLIn
       Dt = CFL/maxval(CFLCoef_G)
       CFLCoef_G = Dt*CFLCoef_G
    end if 
    !\
    ! First order monotone scheme
    !/
    Source_C = -SumDeltaHPlus_G(1:nQ, 1:nP)*DeltaMinusF_G(1:nQ, 1:nP)
    !\
    ! Calculate Face-X fluxes.
    !/
    do iP = 1, nP; do iQ = 0, nQ
       if(DeltaH_FX(iQ, iP) > 0.0)then
          Flux_FX(iQ, iP) = DeltaH_FX(iQ, iP)*0.50*&
               (1.0 - CFLCoef_G(iQ,   iP))                           *&
               pair_superbee(DeltaMinusF_G(iQ, iP), DeltaMinusF_G(iQ+1, iP))
       else
          Flux_FX(iQ, iP) = DeltaH_FX(iQ, iP)*0.50*&
               (1.0 - CFLCoef_G(iQ+1, iP))                           *&
               pair_superbee(DeltaMinusF_G(iQ, iP), DeltaMinusF_G(iQ+1, iP))
       end if
    end do; end do
    !\
    ! Calculate Face-Y fluxes.
    !/
    do iP = 0, nP; do iQ = 1, nQ
       if(DeltaH_FY(iQ, iP) > 0.0)then
          Flux_FY(iQ, iP) = DeltaH_FY(iQ, iP)*0.50*&
               (1.0 - CFLCoef_G(iQ,   iP))                           *&
               pair_superbee(DeltaMinusF_G(iQ, iP), DeltaMinusF_G(iQ, iP+1))
       else
          Flux_FY(iQ, iP) = DeltaH_FY(iQ, iP)*0.50*&
               (1.0 - CFLCoef_G(iQ, iP+1))                           *&
               pair_superbee(DeltaMinusF_G(iQ, iP), DeltaMinusF_G(iQ, iP+1))
       end if
    end do; end do
    !\
    ! Finalize
    !/
    Source_C = (Source_C + Flux_FX(0:nQ-1, 1:nP) - Flux_FX(1:nQ, 1:nP) +  &
         Flux_FY(1:nQ, 0:nP-1) - Flux_FY(1:nQ, 1:nP) )*Dt*vInv_G(1:nQ,1:nP)
    if(present(CFLOut))CFLOut = CFL
    if(present(DtOut ))DtOut  = Dt
  end subroutine explicit
    
end module ModPoissonBracket
