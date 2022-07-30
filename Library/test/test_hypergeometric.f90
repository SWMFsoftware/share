!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!BOP
program test
  use ModHyperGeometric
  use ModPlotFile
  implicit none
  integer, parameter :: nStep = 400, SCR_ = 1, L0Ext_=2, Approx_ = 3
  integer, parameter :: AxialS_ = 1, PoloidalS_ = 2, &
       AxialU_ = 3, PoloidalU_ = 4, ToroidalU_ = 5
  real   , parameter :: DeltaKappaPrime = 0.0010
  integer            :: iLoop
  real               :: R0, a, Q0, Q1, CurrentE,         &
       KappaPrime, KappaPrime2, Kappa, Kappa2, Kappa3,   &
       KappaPrime0, KappaPrime02, Kappa0, Kappa02, Kappa03
  real :: Var_VI(SCR_:ToroidalU_,nStep), Coord_I(nSTep)
  !--------------------------
  !!!!!!!!!!!!!!!!!!!!!! Plot inductances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do iLoop = 1, nStep
     KappaPrime0  = iLoop*DeltaKappaPrime
     KappaPrime02 =  KappaPrime0**2 
     a  = 2*KappaPrime0/(1 - KappaPrime02)
     R0 = (1 + KappaPrime02)/(1 - KappaPrime02)
     Coord_I(iLoop) = a/R0
     Var_VI(SCR_,iLoop)   = scr_inductance(KappaPrime02)
     Var_VI(L0Ext_,iLoop) = l0_ext_inductance(KappaPrime02)
     Var_VI(Approx_,iLoop) = log(8.0/Coord_I(iLoop)) - 2.0
  end do
  call save_plot_file(NameFile='test_inductance.out', &
       TypeFileIn='ascii'                     ,&
       NameVarIn = &
       'a_over_R0 Inductance_SCR L_Ext_inductance Approximate_inductance', &
       StringFormatIn = '(4es18.10)'          ,&
       Coord1In_I = Coord_I                   ,&
       VarIn_VI = Var_VI(SCR_:Approx_,1:nStep) )
  !!!!!!!!!!!!!!!!!!!!! Inner field !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Parameters on the plasma filament surface at which \kappa^\prime_0 = 0.1
  KappaPrime0 = 0.10 
  KappaPrime02 = KappaPrime0**2
  Kappa02 = 1 - KappaPrime02 
  Kappa0 = sqrt(Kappa02); Kappa03 = Kappa0*Kappa02
  ! Inner field from the current over the toroid surface
  ! Eq. (31), constant field factor for surface current
  Q0 = 0.125*toroid_p(0, KappaPrime2In=KappaPrime02)/&
       toroid_q(0,KappaPrime2In=KappaPrime02)
  do iLoop = 1, nStep/2
     ! \kappa^\prime ranges from 0 to \kappa^\prime_0 = 0.1
     KappaPrime = iLoop*0.50*DeltaKappaPrime
     KappaPrime2 = KappaPrime**2
     Kappa2 = 1 - KappaPrime2; Kappa = sqrt(Kappa2); Kappa3 = Kappa*Kappa2
     Coord_I(iLoop) = KappaPrime
     ! Eq. 26 for surface curent
     Var_VI(AxialS_,iLoop) = Q0*Kappa3*toroid_q(0,&
          KappaPrime2In=KappaPrime2)
     Var_VI(PoloidalS_,iLoop) = 3*KappaPrime*Q0*Kappa3*toroid_q(1,&
          KappaPrime2In=KappaPrime2)
  end do
  ! Eq. 36, constant field factor for uniform current
  Q1 = 0.125*toroid_p(1,KappaPrime2In=KappaPrime02)/&
       (toroid_q(1,KappaPrime2In=KappaPrime02))
  ! Constant current I_E
  CurrentE = -1/(3*KappaPrime02*Kappa0*toroid_q(1,KappaPrime2In=KappaPrime02))
  ! Inner field from the uniform current
  do iLoop = 1, nStep/2
     KappaPrime = iLoop*0.50*DeltaKappaPrime
     KappaPrime2 = KappaPrime**2
     Kappa2 = 1 - KappaPrime2; Kappa = sqrt(Kappa2); Kappa3 = Kappa*Kappa2
     ! Eqs. 35
     Var_VI(AxialU_,iLoop) = Q1*Kappa3*toroid_q(0,KappaPrime2In = &
          KappaPrime2) + CurrentE
     Var_VI(PoloidalU_,iLoop) = 3*Q1*KappaPrime*Kappa3*toroid_q(1,    &
          KappaPrime2In=KappaPrime2)
     Var_VI(ToroidalU_,iLoop) = sqrt(&
          toroid_p(1,KappaPrime2In=KappaPrime02)  /            &
          (4*KappaPrime02*Kappa0*(toroid_q(1,KappaPrime2In=KappaPrime02))**2) &
          *(Kappa03*toroid_q(0,KappaPrime2In=KappaPrime02)     &
          - Kappa3*toroid_q(0,KappaPrime2In=KappaPrime2) )     )
  end do
  !  External magnetic field 
  do iLoop = 1 + nStep/2, nStep
     KappaPrime = iLoop*0.50*DeltaKappaPrime
     KappaPrime2 = KappaPrime**2
     Kappa2 = 1 - KappaPrime2; Kappa = sqrt(Kappa2); Kappa3 = Kappa*Kappa2
     Coord_I(iLoop) = KappaPrime
     ! Eqs. 27
     Var_VI(AxialS_,iLoop) = 0.125*Kappa3*toroid_p(0, KappaPrime2In = &
          KappaPrime2)
     Var_VI(AxialU_,iLoop) = 0.125*Kappa3*toroid_p(0, KappaPrime2In = &
          KappaPrime2)
     Var_VI(PoloidalS_,iLoop) = 0.375*Kappa3*KappaPrime* &
          toroid_p(1, KappaPrime2In = KappaPrime2)
     Var_VI(PoloidalU_,iLoop) = 0.375*Kappa3*KappaPrime* &
          toroid_p(1, KappaPrime2In = KappaPrime2)
     Var_VI(ToroidalU_,iLoop)  = 0.0
  end do
  call save_plot_file(NameFile='test_fields.out', &
       TypeFileIn='ascii'                     ,&
       NameVarIn=&
       'kappa_prime AxialS PoloidalS AxialU PoloidalU ToroidalU', &
       StringFormatIn = '(6es18.10)'          ,&
       Coord1In_I = Coord_I                   ,&
       VarIn_VI = Var_VI(AxialS_:ToroidalU_,1:nStep))
end program test

