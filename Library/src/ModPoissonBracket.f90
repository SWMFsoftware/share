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
  public :: explicit3
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
  !=========================================================================
  
  subroutine explicit(nQ, nP, VDF_G, Hamiltonian_N,   &
       Volume_G, Source_C, DtIn, CFLIn, DtOut, CFLOut)
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
    character(LEN=*), parameter:: NameSub = 'calc_poisson_bracket'
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

  
  
  !============================================================================================================================================================
  subroutine explicit3(nQ, nP, nR, VDF_G, Hamiltonian_N, Volume_G, Source_C, Hamiltonian01_N, Hamiltonian02_N, Hamiltonian03_N, Hamiltonian13_N, Hamiltonian23_N,   &
      DtIn, CFLIn, DtOut, CFLOut)
    !\
    ! solve the contribution to the
    ! numerical flux from three Poisson brackets with four subscripts, where "0" always refers to time!
    ! For our particular case: {f,H_1}_0,3+{f,H_2}_1,2+{f,H_3}_3,2=0
    ! In our particular cause: 0:t, 1: s_L, 2: \mu, 3:p^3/3,
    !/
    integer, intent(in) :: nQ           !number of cells along the first axis, which is s_L axis in our case
    integer, intent(in) :: nP           !number of cells along the second axis, which is \mu axis in our case
    integer, intent(in) :: nR           !number of cells along the third axis, which is p^3/3 axis in our case (note we actually take a "ln" for the p^3/3 axis)
    !\
    !Distribution function with gc. Two layers of face ghostcels and one level
    !of corner ghost cells are used
    !/
    real, intent(in) :: VDF_G(-1:nQ+2,-1:nP+2,-1:nR+2) 
    !\ 
    ! Hamiltonian functions in nodes. One layer of ghost nodes is used. At most three Hamiltonian functions are allowed
    ! The numbers in the name of Hamiltonian funcitons refer to the subscripts of Poisson brackets.
    ! Hamiltonian_N is necessary for this code. In any case, this Hamiltonian function always refer to the first non-time-dependent Hamiltonian in the equation.
    !/
    real, intent(in) :: Hamiltonian_N(-1:nQ+1,-1:nP+1,-1:nR+1)!The second Hamiltonian function in our equation. It's necessary to the code. In every case, the "1" and "2" in the array always refer to the first and second subscript.
    real, optional, intent(in) :: Hamiltonian01_N(-1:nQ+1,-1:nP+1,-1:nR+1,0:1), Hamiltonian02_N(-1:nQ+1,-1:nP+1,-1:nR+1,0:1), Hamiltonian03_N(-1:nQ+1,-1:nP+1,-1:nR+1,0:1)!The time dependent Hamiltonian functions. An optional part.
    real, optional, intent(in) :: Hamiltonian13_N(-1:nQ+1,-1:nR+1,-1:nR+1), Hamiltonian23_N(-1:nQ+1,-1:nR+1,-1:nR+1)!The Hamiltonian functions in the Poisson bracket with corresponding subscript. An optional part.
    
    !\
    ! Total Volume. One layer of face ghost cells is used
    !/  
    real, intent(in) :: Volume_G(0:nQ+1,0:nP+1,0:nR+1)
    !\
    ! Inverse volume. One layer of face ghost cells is used
    !/
    real :: vInv_G(0:nQ+1,0:nP+1,0:nR+1)
    !\
    ! Contribution to the conservative source (flux divergence) for 
    ! the Poisson Bracket:
    !/
    !send the source_c back to the main code
    real, intent(out) :: Source_C(1:nQ, 1:nP, 1:nR)  
    !\
    !OPTIONAL PARAMETERS:
    !/                                                                                                                                                                                                                                     
    real, optional, intent(in) :: DtIn, CFLIn   !Options to set time step
    real, optional, intent(out):: DtOut, CFLOut !Options to report time step
    !\
    ! Local variables
    !/
    !Loop variables:
    integer :: iQ, iP, iR
    
    ! Variations of VDF (one layer of ghost cell values):
    real :: DeltaMinusF_G(0:nQ+1, 0:nP+1, 0:nR+1)
    !\
    ! Variations of Hamiltonian functions (one layer of ghost faces, 
    ! for three directions:
    ! Note that the ranges of different subscripts in the following matrix are different because they depend on the specific direction.
    ! The "X", "Y" and "Z" in the following arrays refer to the corredponding direction of axis.
    !/
    real :: DeltaH_FX(-1:nQ+1,0:nP+1,0:nR+1)
    real :: DeltaH_FY(0:nQ+1,-1:nP+1,0:nR+1)
    real :: DeltaH_FZ(0:nQ+1,0:nP+1,-1:nR+1)
    real :: SumDeltaHMinus_G(0:nQ+1,0:nP+1,0:nR+1)
    !Fluxes:
    real :: Flux_FX(0:nQ,1:nP,1:nR), Flux_FY(1:nQ,0:nP,1:nR), Flux_FZ(1:nQ,1:nP,0:nR)
    
    ! Local CFL number:
    real :: CFLCoef_G(0:nQ+1,0:nP+1,0:nR+1)
    
    !Time step
    real :: Dt, CFL
    character(LEN=*), parameter:: NameSub = 'calc_poisson_bracket'
    !---------
    if(present(DtIn))then
       Dt = DtIn
    else
       if(.not.present(CflIn))call CON_stop(&
            'Either CflIn or DtIn should be provided in '//NameSub)
    end if
    vInv_G=1.0/Volume_G
    !\
    ! Calculate DeltaH:
    !/
    DeltaH_FX=0.0
    DeltaH_FY=0.0
    DeltaH_FZ=0.0
    !\
    !      Bracket {F,H12}_x,y         Bracket {F,H13}_x,z           Bracket {F,H23}_y,z
    !      Hamiltonian 12 (xy)         Hamiltonian 13 (xz)           Hamiltonian 23 (yz)
    !      y                           z                             z     
    !      ----------<--------         ----------<----------         -----------<----------- 
    !      |                  |        |                   |         |                     |
    !      v                  ^        v                   |         v                     ^ 
    !      |                  |        |                   |         |                     |
    !      |                  |        |                   |         |                     |
    !      ---------->---------  x     ---------->----------  x      ------------>----------  y
    !/
   
    DeltaH_FX(-1:nQ+1,0:nP+1,0:nR+1) = Hamiltonian_N(-1:nQ+1,0:nP+1,0:nR+1) - Hamiltonian_N(-1:nQ+1,-1:nP,0:nR+1)
    DeltaH_FY(0:nQ+1,-1:nP+1,0:nR+1) = Hamiltonian_N(-1:nQ,-1:nP+1,0:nR+1)  - Hamiltonian_N(0:nQ+1,-1:nP+1,0:nR+1)
    if (present(Hamiltonian13_N)) then
       DeltaH_FX(-1:nQ+1, 0:nP+1, 0:nR+1) = DeltaH_FX(-1:nQ+1, 0:nP+1, 0:nR+1) + Hamiltonian13_N(-1:nQ+1, 0:nP+1, 0:nR+1) - Hamiltonian13_N(-1:nQ+1, 0:nP+1, -1:nR)
       DeltaH_FZ(0:nQ+1, 0:nP+1, -1:nR+1) = DeltaH_FZ(0:nQ+1, 0:nP+1, -1:nR+1) + Hamiltonian13_N(-1:nQ, 0:nP+1, -1:nR+1) - Hamiltonian13_N(0:nQ+1, 0:nP+1, -1:nR+1)
    end if
    
    if (present(Hamiltonian23_N)) then
       DeltaH_FY(0:nQ+1, -1:nP+1, 0:nR+1) = DeltaH_FY(0:nQ+1, -1:nP+1, 0:nR+1) + Hamiltonian23_N(0:nQ+1, -1:nP+1, 0:nR+1) - Hamiltonian23_N(0:nQ+1, -1:nP+1, -1:nR)
       DeltaH_FZ(0:nQ+1, 0:nP+1, -1:nR+1) = DeltaH_FZ(0:nQ+1, 0:nP+1, -1:nR+1) + Hamiltonian23_N(0:nQ+1, -1:nP, -1:nR+1) - Hamiltonian23_N(0:nQ+1, 0:nP+1, -1:nR+1)
    end if

    !\
    !      Bracket {F,H01}_t,x         Bracket {F,H02}_t,y           Bracket {F,H03}_t,z
    !      Hamiltonian 01 (tx)         Hamiltonian 02 (ty)           Hamiltonian 03 (tz)                                                                                                                                                                                       
    !      x                           y                             z                                                                                                                                                                                                         
    !      ----------<----------       ----------<---------          -----------<-----------                                                                                                                                                                                   
    !                                                                                                                                                                                                                                                                          
    !                                                                                                                                                                                                                                                                          
    !                                                                                                                                                                                                                                                                          
    !                                                                                                                                                                                                                                                                          
    !      ---------->----------  t    ---------->----------  t      ----------->-----------  t                                                                                                                                                                                
    !/
    if (present(Hamiltonian01_N)) DeltaH_FX(-1:nQ+1, 0:nP+1, 0:nR+1) = DeltaH_FX(-1:nQ+1, 0:nP+1, 0:nR+1) + Hamiltonian01_N(-1:nQ+1, 0:nP+1, 0:nR+1, 0) - Hamiltonian01_N(-1:nQ+1, 0:nP+1, 0:nR+1, 1)
    if (present(Hamiltonian02_N)) DeltaH_FY(0:nQ+1, -1:nP+1, 0:nR+1) = DeltaH_FY(0:nQ+1, -1:nP+1, 0:nR+1) + Hamiltonian02_N(0:nQ+1, -1:nP+1, 0:nR+1, 0) - Hamiltonian02_N(0:nQ+1, -1:nP+1, 0:nR+1, 1)
    if (present(Hamiltonian03_N)) DeltaH_FZ(0:nQ+1, 0:nP+1, -1:nR+1) = DeltaH_FZ(0:nQ+1, 0:nP+1, -1:nR+1) + Hamiltonian03_N(0:nQ+1, 0:nP+1, -1:nR+1, 0) - Hamiltonian03_N(0:nQ+1, 0:nP+1, -1:nR+1, 1)
    
    
    ! Now, for each cell the value of DeltaH for face in positive 
    ! directions of iQ and iP may be found in the arrays, for
    ! negative directions the should be taken with opposite sign
    
    !\
    ! Calculate DeltaMinusF and SumDeltaH
    !/   
    do iR=0, nR+1; do iP = 0, nP+1; do iQ = 0, nQ+1
       
       SumDeltaHMinus_G(iQ, iP, iR) = min(0.0, DeltaH_FX(iQ,   iP,   iR)) +&
                                      min(0.0, DeltaH_FY(iQ,   iP,   iR)) +&
                                      min(0.0, DeltaH_FZ(iQ,   iP,   iR)) +&
                                      min(0.0,-DeltaH_FX(iQ-1, iP,   iR)) +&
                                      min(0.0,-DeltaH_FY(iQ, iP-1,   iR)) +&
                                      min(0.0,-DeltaH_FZ(iQ,   iP, iR-1))
                                      
       DeltaMinusF_G(iQ, iP, iR) = -(&
           min(0.0, DeltaH_FX(iQ,   iP,   iR))*VDF_G(iQ+1,iP,iR) +&
           min(0.0, DeltaH_FY(iQ,   iP,   iR))*VDF_G(iQ,iP+1,iR) +&
           min(0.0, DeltaH_FZ(iQ,   iP,   iR))*VDF_G(iQ,iP,iR+1) +&
           min(0.0,-DeltaH_FX(iQ-1, iP, iR))*VDF_G(iQ-1,iP,iR) +&
           min(0.0,-DeltaH_FY(iQ, iP-1, iR))*VDF_G(iQ,iP-1,iR) +&
           min(0.0,-DeltaH_FZ(iQ, iP, iR-1))*VDF_G(iQ,iP,iR-1)  &
           )/SumDeltaHMinus_G(iQ,iP,iR) + VDF_G(iQ,iP,iR)
       
       CFLCoef_G(iQ, iP, iR) = - vInv_G(iQ, iP, iR)*SumDeltaHMinus_G(iQ,iP,iR)
       
    end do; end do; end do
    
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
    Source_C = SumDeltaHMinus_G(1:nQ, 1:nP, 1:nR)*DeltaMinusF_G(1:nQ, 1:nP, 1:nR)
    
    !\                                                                                                                                                                                                                                  
    ! Calculate Face-X fluxes.                                                                                                                                                                                                                   
    !/                              
    do iR=1, nR; do iP = 1, nP; do iQ = 0, nQ
       if(DeltaH_FX(iQ, iP, iR) > 0.0)then
          Flux_FX(iQ, iP, iR) = DeltaH_FX(iQ, iP, iR)*0.50*&
               (1.0 - CFLCoef_G(iQ,   iP,  iR))                           *&
               pair_superbee(DeltaMinusF_G(iQ, iP, iR), DeltaMinusF_G(iQ+1, iP, iR))
       else
          Flux_FX(iQ, iP, iR) = DeltaH_FX(iQ, iP, iR)*0.50*&
               (1.0 - CFLCoef_G(iQ+1, iP,  iR))                           *&
               pair_superbee(DeltaMinusF_G(iQ, iP, iR), DeltaMinusF_G(iQ+1, iP, iR))
       end if
    end do; end do; end do
    
    !\
    ! Calculate Face-Y fluxes.
    !/
    do iR=1, nR; do iP = 0, nP; do iQ = 1, nQ
       if(DeltaH_FY(iQ, iP, iR) > 0.0)then
          Flux_FY(iQ, iP, iR) = DeltaH_FY(iQ, iP, iR)*0.50*&
               (1.0 - CFLCoef_G(iQ,   iP,  iR))                           *&
               pair_superbee(DeltaMinusF_G(iQ, iP, iR), DeltaMinusF_G(iQ, iP+1, iR))
       else
          Flux_FY(iQ, iP, iR) = DeltaH_FY(iQ, iP, iR)*0.50*&
               (1.0 - CFLCoef_G(iQ, iP+1,  iR))                           *&
               pair_superbee(DeltaMinusF_G(iQ, iP, iR), DeltaMinusF_G(iQ, iP+1, iR))
       end if
    end do; end do; end do
    
    !\
    ! Calculate Face-Z fluxes.
    !/
    do iR=0, nR; do iP = 1, nP; do iQ = 1, nQ
       if(DeltaH_FZ(iQ, iP, iR) > 0.0)then
          Flux_FZ(iQ, iP, iR) = DeltaH_FZ(iQ, iP, iR)*0.50*&
               (1.0 - CFLCoef_G(iQ,   iP,  iR))                           *&
               pair_superbee(DeltaMinusF_G(iQ, iP, iR), DeltaMinusF_G(iQ, iP, iR+1))
       else
          Flux_FZ(iQ, iP, iR) = DeltaH_FZ(iQ, iP, iR)*0.50*&
               (1.0 - CFLCoef_G(iQ, iP,  iR+1))                           *&
               pair_superbee(DeltaMinusF_G(iQ, iP, iR), DeltaMinusF_G(iQ, iP, iR+1))
       end if
    end do; end do; end do
    
    !\
    ! Finalize
    !/
    Source_C = ( Source_C + Flux_FX(0:nQ-1, 1:nP, 1:nR) - Flux_FX(1:nQ, 1:nP, 1:nR) +  &
         Flux_FY(1:nQ, 0:nP-1, 1:nR) - Flux_FY(1:nQ, 1:nP, 1:nR) +  Flux_FZ(1:nQ, 1:nP, 0:nR-1) - Flux_FZ(1:nQ, 1:nP, 1:nR) )*Dt*vInv_G(1:nQ,1:nP,1:nR)
    
    if(present(CFLOut))CFLOut = CFL
    if(present(DtOut ))DtOut  = Dt
  end subroutine explicit3
  
end module ModPoissonBracket
