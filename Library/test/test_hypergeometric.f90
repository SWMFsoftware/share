!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!BOP
program test
  use ModHyperGeometric
  use ModPlotFile
  implicit none
  integer, parameter :: nStep = 400, Exact_ = 1, Approx_ = 2
  integer, parameter :: AzimuthalS_ = 1, PoloidalS_ = 2, &
       AzimuthalU_ = 3, PoloidalU_ = 4, ToroidalU_ = 5
  real   , parameter :: DeltaKappaPrime = 0.0010
  integer            :: iLoop
  real               :: R0, a, ConstRatio, ConstLev,   &
       KappaPrime, KappaPrime2, Kappa, Kappa2, Kappa3, &
       KappaPrime0, KappaPrime02, Kappa0, Kappa02, Kappa03
  real :: Var_VI(Exact_:ToroidalU_,nStep), Coord_I(nSTep)
  !--------------------------
  do iLoop = 1, nStep
     KappaPrime0  = iLoop*DeltaKappaPrime
     KappaPrime02 =  KappaPrime0**2 
     a  = 2*KappaPrime0/(1 - KappaPrime02)
     R0 = (1 + KappaPrime02)/(1 - KappaPrime02)
     Coord_I(iLoop) = a/R0
     Var_VI(Exact_,iLoop)   = scr_inductance(KappaPrime02)
     Var_VI(Approx_,iLoop) = R0*(log(8.0/ Coord_I(iLoop)) - 2.0)
  end do
  call save_plot_file(NameFile='test_inductance.out', &
       TypeFileIn='ascii'                     ,&
       NameVarIn='a_over_R0 Exact_inductance Approximate_inductance', &
       StringFormatIn = '(3es18.10)'          ,&
       Coord1In_I = Coord_I                    ,&
       VarIn_VI = Var_VI(Exact_:Approx_,1:nStep) )
  ! test Green function (field from the current over the toroid surface
  KappaPrime0 = 0.10 
  KappaPrime02 = KappaPrime0**2
  Kappa02 = 1 - KappaPrime02 
  Kappa0 = sqrt(Kappa02); Kappa03 = Kappa0*Kappa02
  ConstRatio = 0.125*toroid_p(0, KappaPrime2In=KappaPrime02)/&
       toroid_q0(KappaPrime2In=KappaPrime02)
  do iLoop = 1, nStep/2
     KappaPrime = iLoop*0.50*DeltaKappaPrime
     KappaPrime2 = KappaPrime**2
     Kappa2 = 1 - KappaPrime2; Kappa = sqrt(Kappa2); Kappa3 = Kappa*Kappa2
     Coord_I(iLoop) = KappaPrime
     Var_VI(AzimuthalS_,iLoop) = ConstRatio*Kappa3*toroid_q0(&
          KappaPrime2In=KappaPrime2)
     Var_VI(PoloidalS_,iLoop) = KappaPrime*ConstRatio*Kappa3*&
          toroid_dq0du(KappaPrime2In=KappaPrime2)
  end do
  ConstRatio = 0.125*toroid_dpdu(0, KappaPrime2In=KappaPrime02)/&
       (toroid_dq0du(KappaPrime2In=KappaPrime02)*KappaPrime02)
  ConstLev = -1/( toroid_dq0du(KappaPrime2In=KappaPrime02)*KappaPrime02*&
       Kappa0)
  do iLoop = 1, nStep/2
     KappaPrime = iLoop*0.50*DeltaKappaPrime
     KappaPrime2 = KappaPrime**2
     Kappa2 = 1 - KappaPrime2; Kappa = sqrt(Kappa2); Kappa3 = Kappa*Kappa2
     Var_VI(AzimuthalU_,iLoop) = ConstRatio*Kappa3*toroid_q0(KappaPrime2In = &
          KappaPrime2) + ConstLev
     Var_VI(PoloidalU_,iLoop) = KappaPrime*ConstRatio*toroid_dq0du(&
          KappaPrime2In=KappaPrime2)*Kappa3
     Var_VI(ToroidalU_,iLoop) = sqrt(6*(-ConstRatio)*ConstLev*(&
          Kappa03*toroid_q0(KappaPrime2In=KappaPrime02)&
          - Kappa3*toroid_q0(KappaPrime2In=KappaPrime2)))
  end do
  do iLoop = 1 + nStep/2, nStep
     KappaPrime = iLoop*0.50*DeltaKappaPrime
     KappaPrime2 = KappaPrime**2
     Kappa2 = 1 - KappaPrime2; Kappa = sqrt(Kappa2); Kappa3 = Kappa*Kappa2
     Coord_I(iLoop) = KappaPrime
     Var_VI(AzimuthalS_,iLoop) = 0.125*Kappa3*toroid_p(0, KappaPrime2In = &
          KappaPrime2)
     Var_VI(AzimuthalU_,iLoop) = 0.125*Kappa3*toroid_p(0, KappaPrime2In = &
          KappaPrime2)
     Var_VI(PoloidalS_,iLoop) = 0.125*Kappa3*toroid_dpdu(0, KappaPrime2In = &
          KappaPrime2)/KappaPrime
     Var_VI(PoloidalU_,iLoop) = 0.125*Kappa3*toroid_dpdu(0, KappaPrime2In = &
          KappaPrime2)/KappaPrime
     Var_VI(ToroidalU_,iLoop)  = 0.0
  end do
  call save_plot_file(NameFile='test_fields.out', &
       TypeFileIn='ascii'                     ,&
       NameVarIn=&
       'kappa_prime AzimuthalS PoloidalS AzimuthalU_ PoloidalU ToroidalU', &
       StringFormatIn = '(6es18.10)'          ,&
       Coord1In_I = Coord_I                   ,&
       VarIn_VI = Var_VI(Exact_:ToroidalU_,1:nStep))
end program test

