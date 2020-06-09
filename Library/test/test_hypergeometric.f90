!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!BOP
program test
  use ModHyperGeometric
  use ModPlotFile
  implicit none
  integer, parameter :: nStep = 400
  real   , parameter :: DeltaKappaPrime = 0.0010
  integer            :: iLoop
  real               :: R0, a, KappaPrime
  real :: Exact_I(nStep), Approx_I(nStep), A2R0_I(nSTep)
  !--------------------------
  do iLoop = 1, nStep
     KappaPrime = iLoop*DeltaKappaPrime
     a  = 2*KappaPrime/(1 - KappaPrime**2)
     R0 = (1 + KappaPrime**2)/(1 - KappaPrime**2)
     A2R0_I(iLoop) = a/R0
     Exact_I(iLoop) = scr_inductance(KappaPrime**2)
     Approx_I(iLoop) = R0*(log(8.0/ A2R0_I(iLoop)) - 2.0)
     write(*,'(3es18.10)')A2R0_I(iLoop), Exact_I(iLoop), Approx_I(iLoop)
  end do
end program test

