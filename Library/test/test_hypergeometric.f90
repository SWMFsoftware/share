!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program test
  use ModHyperGeometric
  use ModPlotFile
  use ModCurrentFilament
  use ModFieldGS, ONLY: test_gs22=>test
  implicit none
  integer, parameter :: nStep = 400, SCR_ = 1, L0Ext_=2, Approx_ = 3
  integer, parameter :: AxialS_ = 1, PoloidalS_ = 2, &
       AxialU_ = 3, PoloidalU_ = 4, ToroidalU_ = 5
  real   , parameter :: DeltaKappaPrime = 0.0010
  integer            :: iLoop
  real               :: R0, a, Amplitude_I(Axial_:Toroidal_), &
       KappaPrime, KappaPrime2, Kappa2, Kappa3
  real :: Var_VI(SCR_:ToroidalU_,nStep), Coord_I(nSTep)
  !----------------------------------------------------------------------------
  !!!!!!!!!!!!!!!!!!!!!! Plot inductances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do iLoop = 1, nStep
     KappaPrime  = iLoop*DeltaKappaPrime
     KappaPrime2 =  KappaPrime**2
     a  = 2*KappaPrime/(1 - KappaPrime2)
     R0 = (1 + KappaPrime2)/(1 - KappaPrime2)
     Coord_I(iLoop) = a/R0
     Var_VI(SCR_,iLoop)   = scr_inductance(KappaPrime2)
     Var_VI(L0Ext_,iLoop) = l0_ext_inductance(KappaPrime2)
     Var_VI(Approx_,iLoop) = log(8.0/Coord_I(iLoop)) - 2.0
  end do
  call save_plot_file(NameFile='test_inductance.out', &
       TypeFileIn='ascii'                     ,&
       NameVarIn = &
       'a_over_R0 Inductance_SCR L_Ext_inductance Approximate_inductance', &
       StringFormatIn = '(4es18.10)'          ,&
       Coord1In_I = Coord_I                   ,&
       VarIn_VI = Var_VI(SCR_:Approx_,1:nStep) )
  !!!!!!!!!!!!!!!!!!!!!  Internal field  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Parameters on the plasma filament surface at which \kappa^\prime_0 = 0.1
  UseUniformCurrent = .true.;  UseSurfaceCurrent = .true.
  a  = 0.20/0.990;  R0 = 1.010/0.990
  call set_filament_geometry(a, R0)
  do iLoop = 1, nStep/2
     !
     ! \kappa^\prime ranges from 0 to \kappa^\prime_0 = 0.1
     KappaPrime     = iLoop*0.50*DeltaKappaPrime
     Coord_I(iLoop) = KappaPrime
     !
     ! For surface curent
     !
     call surface_current_field(KappaPrime2In = KappaPrime**2, &
          Amplitude_I=Amplitude_I)
     Var_VI(AxialS_,iLoop)    = Amplitude_I(Axial_)
     Var_VI(PoloidalS_,iLoop) = KappaPrime*Amplitude_I(Poloidal_)
     !
     ! Internal field from the uniform current
     !
     call uniform_current_field(KappaPrime2In = KappaPrime**2, &
          Amplitude_I=Amplitude_I)
     ! Eqs. 35
     Var_VI(AxialU_,iLoop)    = Amplitude_I(Axial_)
     Var_VI(PoloidalU_,iLoop) = KappaPrime*Amplitude_I(Poloidal_)
     Var_VI(ToroidalU_,iLoop) = Amplitude_I(Toroidal_)
  end do
  !!!!!!!!!!!!!!!!!!!!!  External magnetic field !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do iLoop = 1 + nStep/2, nStep
     KappaPrime = iLoop*0.50*DeltaKappaPrime
     Kappa2     = 1 - KappaPrime**2
     Kappa3 = sqrt(Kappa2)*Kappa2
     Coord_I(iLoop) = KappaPrime
     call  external_field(Kappa2In = Kappa2, Amplitude_I=Amplitude_I)
     ! Eqs. 27
     Var_VI(AxialS_,iLoop)    = 0.125*Kappa3*Amplitude_I(Axial_)
     Var_VI(AxialU_,iLoop)    = Var_VI(AxialS_,iLoop)
     Var_VI(PoloidalS_,iLoop) = 0.125*Kappa3*KappaPrime*Amplitude_I(Poloidal_)
     Var_VI(PoloidalU_,iLoop) = Var_VI(PoloidalS_,iLoop)
     Var_VI(ToroidalU_,iLoop) = 0.0
  end do
  call save_plot_file(NameFile='test_fields.out', &
       TypeFileIn='ascii'                     ,&
       NameVarIn=&
       'kappa_prime AxialS PoloidalS AxialU PoloidalU ToroidalU', &
       StringFormatIn = '(6es18.10)'          ,&
       Coord1In_I = Coord_I                   ,&
       VarIn_VI = Var_VI(AxialS_:ToroidalU_,1:nStep))
  call test_gs22
end program test
!==============================================================================

