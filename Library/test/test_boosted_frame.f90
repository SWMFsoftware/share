!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!===========================TESTS==============================================
module ModTestBoostedFrame

  ! Test the code in ModBoostedFrame and plot all figures in the paper
  ! Sokolov, I. V. and Gombosi, T. I. (2024),
  ! Physics-Based Forecasting of Tomorrow's Solar Wind at 1 AU
  ! Revision history:
  ! 12/18/2024 - Tests the schemes and plot figures for the paper
  !              by Sokolov&Gombosi (2004)
  use ModExactRs
  use ModPlotFile
  use ModBoostedFrame, ONLY: set_param, update_hydro_var, update_mhd8wave_var
  implicit none
  real, parameter ::  cThird = 1.0/3.0, &
       cTwoThird = 2.0/3.0
contains
  !============================================================================
  subroutine get_hydro_rs
    real, parameter :: Gamma = 5.0/3.0
    real :: Time = 0.04, TimeForecast = 0.022, TimeLimit = 0.1, X10 = 0.50
    real :: tOffsetPerR = 0.04
    real :: ExactRhoX_I(1000), ExactRhoForecastX_I(1000), X_I(1000)
    real :: Time_I(1000), ExactRhoT10_I(1000), ExactRhoForecastT10_I(1000)
    ! Misc:
    real :: Lambda, Uaux, Paux
    ! Loop variable:
    integer :: i
    !--------------------------------------------------------------------------

    call exact_rs_set_gamma(Gamma)
    ! set left and right state
    RhoL = 8
    pL   = 480
    UnL  = 0
    RhoR = 1
    pR   = 1
    UnR  = 0
    ! exact solver
    call exact_rs_pu_star
    ! Fill in coord arrays:
    do i = 1, 1000
       X_I(i)    = -0.5 + (i - 0.5)*1.0e-3
       Time_I(i) = TimeLimit*(i-0.5)*1.0e-3
    end do
    ! Distribution over x:
    open(11,file='SolutionsX.out',status='replace')
    do i = 1, 1000
       Lambda = X_I(i)/Time
       call exact_rs_sample(Lambda, ExactRhoX_I(i),Uaux, Paux)
       Lambda = X_I(i)/TimeForecast
       call exact_rs_sample(Lambda/(1 + max(Lambda, 0.0)*tOffsetPerR),&
            ExactRhoForecastX_I(i), Uaux, Paux)
       write(11,*)X_I(i) + 0.5, ExactRhoX_I(i), ExactRhoForecastX_I(i)
    end do
    close(11)
    ! Distribution over t:
    open(11,file='SolutionsT.out',status='replace')
    do i = 1, 1000
       Lambda = X10/Time_I(i)
       call exact_rs_sample(Lambda, ExactRhoT10_I(i), Uaux, Paux)
       call exact_rs_sample(Lambda/(1 + max(Lambda, 0.0)*tOffsetPerR),&
            ExactRhoForecastT10_I(i), Uaux, Paux)
       write(11,*)Time_I(i), ExactRhoT10_I(i), ExactRhoForecastT10_I(i)
    end do
    close(11)
  end subroutine get_hydro_rs
  !============================================================================
  subroutine test_hydro_x

    real, parameter :: Gamma = 5.0/3.0
    real :: TimeForecast = 0.022, tOffsetPerR = 0.04
    integer, parameter :: nCell = 100
    real, parameter :: Ds = 1.0/nCell
    integer, parameter :: Rho_ = 1, Ux_ = 2, Uz_= 4, RhoUx_ = 2, &
         RhoUz_ = 4, Energy_ = 5, P_ = 5
    real :: State_VC(Rho_:P_,nCell), StateOld_VC(Rho_:P_,nCell)
    integer, parameter :: iLogVar_V(2) = [Rho_, P_]
    real :: Flux_VF(Rho_:Energy_,0:nCell)
    real :: Cleft_F(0:nCell)
    real :: Cright_F(0:nCell)
    ! Interpolated left and right states at the face
    real    :: pLeft_VF(Rho_:P_,0:nCell)
    real    :: pRight_VF(Rho_:P_,0:nCell)
    real    :: dVarDown_V(Rho_:P_), dVarUp_V(Rho_:P_), Source_V(Rho_:Energy_)
    real    :: tSimulation, Dt
    integer :: nStep, i
    !--------------------------------------------------------------------------

    call exact_rs_set_gamma(Gamma)

    ! Init boosted frame
    call set_param(tOffsetPerRin=tOffsetPerR)

    ! Initial state for forecast numerical test
    State_VC(Ux_:Uz_,:) = 0
    State_VC(Rho_,:nCell/2) = 8
    State_VC(Rho_,nCell/2+1:) = 1
    State_VC(P_,:nCell/2) = 480
    State_VC(P_,nCell/2+1:) = 1
    tSimulation = 0; Dt = 0; nStep = 0
    do
       if(tSimulation >= TimeForecast)EXIT
       Dt = min(Dt,1.001*(TimeForecast - tSimulation))
       tSimulation = tSimulation + Dt
       nStep = nStep + 1
       write(*,*)'nStep, tSimulation, Dt=', nStep, tSimulation, Dt
       ! First Stage. Save old state
       StateOld_VC = State_VC
       call get_fluxes
       ! Advance to halfstep
       do i = 1, nCell/2
          call advance_no_offset(StateOld_VC(:,i)  , &
               0.5*Dt/Ds*(Flux_VF(:,i-1) - Flux_VF(:,i)), &
               State_VC(:,i))
       end do
       do i = nCell/2+1, nCell
          call update_hydro_var(Old_V=StateOld_VC(:,i),&
               Source_V=0.5*Dt/Ds*(Flux_VF(:,i-1) - Flux_VF(:,i)),&
               Radial_D=[1.0, 0.0, 0.0],&
               New_V=State_VC(:,i))
       end do
       ! Second stage
       call get_fluxes
       ! Advance to full time step
       do i = 1,nCell/2
          call advance_no_offset(StateOld_VC(:,i)  , &
               Dt/Ds*(Flux_VF(:,i-1) - Flux_VF(:,i)), &
               State_VC(:,i))
       end do
       do i = nCell/2+1, nCell
          call update_hydro_var(Old_V=StateOld_VC(:,i),&
               Source_V=Dt/Ds*(Flux_VF(:,i-1) - Flux_VF(:,i)),&
               Radial_D=[1.0, 0.0, 0.0],&
               New_V=State_VC(:,i))
       end do
       call get_time_step
       if(Dt<=0.0)then
          write(*,*)'Negative time step, reduce the time offset'
          stop
       end if
    end do
    call save_plot_file(NameFile='test_forecast.out', &
         TypeFileIn='ascii', TimeIn=tSimulation, nStepIn = nStep, &
         NameVarIn='X Rho ux uy uz p', &
         CoordMinIn_D=[0.5/nCell], &
         CoordMaxIn_D=[1 - 0.5/nCell], &
         VarIn_VI = State_VC )
  contains
    !==========================================================================
    subroutine get_fluxes
      integer :: iCell
      real, parameter :: Beta = 1.50
      ! Left boundary:
      !------------------------------------------------------------------------
      call get_godunov_flux(State_VC(:,1), State_VC(:,1), &
           Flux_VF(:,0), Cleft_F(0), Cright_F(0))
      State_VC(iLogVar_V,:) = log(State_VC(iLogVar_V,:))
      pLeft_VF(:,1) = State_VC(:,1)
      dVarUp_V = State_VC(:,2) - State_VC(:,1)
      do iCell = 2, nCell-1
         dVarDown_V = dVarUp_V
         ! Calculate Up slope
         dVarUp_V = State_VC(:,iCell+1) - State_VC(:,iCell)
         pRight_VF(:,iCell-1) = State_VC(:,iCell) - &
              (sign(0.25, dVarUp_V) + sign(0.25, dVarDown_V))* &
              min(Beta*abs(dVarUp_V), Beta*abs(dVarDown_V), &
              cThird*abs(2*dVarDown_V + dVarUp_V))
         pLeft_VF(iLogVar_V,iCell-1) = exp(pLeft_VF(iLogVar_V,iCell-1))
         pRight_VF(iLogVar_V,iCell-1) = exp(pRight_VF(iLogVar_V,iCell-1))
         call get_godunov_flux(pLeft_VF(:,iCell-1), pRight_VF(:,iCell-1), &
              Flux_VF(:,iCell-1), Cleft_F(iCell-1), Cright_F(iCell-1))
         pLeft_VF(:,iCell) = State_VC(:,iCell) + &
              (sign(0.25, dVarUp_V) + sign(0.25, dVarDown_V))* &
              min(Beta*abs(dVarUp_V), Beta*abs(dVarDown_V), &
              cThird*abs(dVarDown_V + 2*dVarUp_V))
         State_VC(iLogVar_V,iCell) = exp(State_VC(iLogVar_V,iCell))
      end do
      pLeft_VF(iLogVar_V,nCell-1) = exp(pLeft_VF(iLogVar_V,nCell-1))
      State_VC(iLogVar_V,nCell) = exp(State_VC(iLogVar_V,nCell))
      call get_godunov_flux(pLeft_VF(:,nCell-1), State_VC(:,nCell), &
           Flux_VF(:,nCell-1), Cleft_F(nCell-1), Cright_F(nCell-1))
      call get_godunov_flux(State_VC(:,nCell), State_VC(:,nCell), &
           Flux_VF(:,nCell), Cleft_F(nCell), Cright_F(nCell))
    end subroutine get_fluxes
    !==========================================================================
    subroutine get_time_step
      real :: CflLocal = 0.9, Dt_C(nCell)
      integer :: iCell
      !------------------------------------------------------------------------
      do iCell = 1, nCell/2
         ! No time offset in these cells
         Dt_C(iCell) = CflLocal*Ds/&
              max(Cright_F(iCell-1), -Cleft_F(iCell), 1e-30)
      end do
      do iCell = nCell/2+1, nCell
         ! Account for the time offset in these cells
         Dt_C(iCell) = CflLocal*Ds*&
              (1 - max(Cright_F(iCell-1), 0.0)*tOffsetPerR)/&
              max(Cright_F(iCell-1), -Cleft_F(iCell), 1e-30)
      end do
      Dt = minval(Dt_C)
    end subroutine get_time_step
    !==========================================================================
    subroutine advance_no_offset(StateIn_V, Source_V, StateOut_V)
      real, intent(in) :: StateIn_V(Rho_:P_), Source_V(Rho_:Energy_)
      real, intent(out):: StateOut_V(Rho_:P_)
      real :: RhoBar
      !------------------------------------------------------------------------
      StateOut_V = 0.0
      StateOut_V(Rho_) = StateIn_V(Rho_) + Source_V(Rho_)
      StateOut_V(Ux_)   = (StateIn_V(Rho_)*StateIn_V(Ux_) + Source_V(RhoUx_))/&
           StateOut_V(Rho_)
      StateOut_V(P_) = StateIn_V(P_) + (Gamma - 1)*&
           (Source_V(Energy_) + 0.5*(StateIn_V(Rho_)*StateIn_V(Ux_)**2 - &
           StateOut_V(Rho_)*StateOut_V(Ux_)**2) )
    end subroutine advance_no_offset
    !==========================================================================
    subroutine get_godunov_flux(pLeft_V, pRight_V, &
         Flux_V, Cleft, Cright)

      real,   intent(in) :: pLeft_V(Rho_:P_), pRight_V(Rho_:P_)
      real,  intent(out) :: Flux_V(Rho_:Energy_)
      real,  intent(out) :: Cleft, Cright

      real :: UnFace, RhoFace, pFace

      !------------------------------------------------------------------------
      RhoL = pLeft_V(Rho_ )
      UnL  = pLeft_V(Ux_)
      pL   = pLeft_V(P_)

      RhoR = pRight_V(Rho_ )
      UnR  = pRight_V(Ux_)
      pR   = pRight_V(P_)

      ! exact solver

      call exact_rs_pu_star
      call exact_rs_sample(0.0, RhoFace, UnFace, pFace)

      ! Store WL WR

      Cleft  = WL
      Cright = WR

      ! get the flux
      Flux_V(Rho_) = RhoFace * UnFace
      Flux_V(RhoUx_) = RhoFace*UnFace**2 + pFace
      Flux_V(RhoUx_+1:RhoUz_) = 0
      Flux_V(Energy_) = &
           (0.50*RhoFace*UnFace**2 + pFace*Gamma/(Gamma - 1))*UnFace
    end subroutine get_godunov_flux
    !==========================================================================
  end subroutine test_hydro_x
  !============================================================================
  subroutine test_hydro_t

    real, parameter :: Gamma = 5.0/3.0
    real :: TimeLimit = 0.1, tOffsetPerR = 0.04
    integer, parameter :: nCell = 200
    real, parameter :: Ds = 2.0/nCell
    integer, parameter :: Rho_ = 1, Ux_ = 2, Uz_= 4, RhoUx_ = 2, &
         RhoUz_ = 4, Energy_ = 5, P_ = 5
    real :: State_VC(Rho_:P_,nCell), StateOld_VC(Rho_:P_,nCell)
    integer, parameter :: iLogVar_V(2) = [Rho_, P_]
    real :: Flux_VF(Rho_:Energy_,0:nCell)
    real :: Cleft_F(0:nCell)
    real :: Cright_F(0:nCell)
    ! Interpolated left and right states at the face
    real    :: pLeft_VF(Rho_:P_,0:nCell)
    real    :: pRight_VF(Rho_:P_,0:nCell)
    real    :: dVarDown_V(Rho_:P_), dVarUp_V(Rho_:P_), Source_V(Rho_:Energy_)
    real    :: tSimulation, Dt
    integer :: nStep, i
    !--------------------------------------------------------------------------

    call exact_rs_set_gamma(Gamma)
    ! Init boosted frame
    call set_param(tOffsetPerRin=tOffsetPerR)

    open(11,file='TestForecastT.out',status='replace')
    ! Initial state for forecast numerical test
    State_VC(Ux_:Uz_,:) = 0
    State_VC(Rho_,:nCell/2) = 8
    State_VC(Rho_,nCell/2+1:) = 1
    State_VC(P_,:nCell/2) = 480
    State_VC(P_,nCell/2+1:) = 1
    tSimulation = 0; Dt = 0; nStep = 0
    do
       if(tSimulation >= 0.5*TimeLimit)EXIT
       Dt = min(Dt,1.001*(0.5*TimeLimit - tSimulation))
       tSimulation = tSimulation + Dt
       nStep = nStep + 1
       ! First Stage. Save old state
       StateOld_VC = State_VC
       call get_fluxes
       ! Advance to halfstep
       do i = 1, nCell/2
          call advance_no_offset(StateOld_VC(:,i)  , &
               0.5*Dt/Ds*(Flux_VF(:,i-1) - Flux_VF(:,i)), &
               State_VC(:,i))
       end do
       do i = nCell/2+1, nCell
          call update_hydro_var(Old_V=StateOld_VC(:,i),&
               Source_V=0.5*Dt/Ds*(Flux_VF(:,i-1) - Flux_VF(:,i)),&
               Radial_D=[1.0, 0.0, 0.0],&
               New_V=State_VC(:,i))
       end do
       ! Second stage
       call get_fluxes
       ! Advance to full time step
       do i = 1,nCell/2
          call advance_no_offset(StateOld_VC(:,i)  , &
               Dt/Ds*(Flux_VF(:,i-1) - Flux_VF(:,i)), &
               State_VC(:,i))
       end do
       do i = nCell/2+1, nCell
          call update_hydro_var(Old_V=StateOld_VC(:,i),&
               Source_V=Dt/Ds*(Flux_VF(:,i-1) - Flux_VF(:,i)),&
               Radial_D=[1.0, 0.0, 0.0],&
               New_V=State_VC(:,i))
       end do
       call get_time_step
       if(Dt<=0.0)then
          write(*,*)'Negative time step, reduce the time offset'
          stop
       end if
       write(11,*)tSimulation, 0.5*(State_VC(Rho_,3*nCell/4) + &
            State_VC(Rho_,3*nCell/4+1))
    end do
    close(11)
  contains
    !==========================================================================
    subroutine get_fluxes
      integer :: iCell
      real, parameter :: Beta = 1.50
      ! Left boundary:
      !------------------------------------------------------------------------
      call get_godunov_flux(State_VC(:,1), State_VC(:,1), &
           Flux_VF(:,0), Cleft_F(0), Cright_F(0))
      State_VC(iLogVar_V,:) = log(State_VC(iLogVar_V,:))
      pLeft_VF(:,1) = State_VC(:,1)
      dVarUp_V = State_VC(:,2) - State_VC(:,1)
      do iCell = 2, nCell-1
         dVarDown_V = dVarUp_V
         ! Calculate Up slope
         dVarUp_V = State_VC(:,iCell+1) - State_VC(:,iCell)
         pRight_VF(:,iCell-1) = State_VC(:,iCell) - &
              (sign(0.25, dVarUp_V) + sign(0.25, dVarDown_V))* &
              min(Beta*abs(dVarUp_V), Beta*abs(dVarDown_V), &
              cThird*abs(2*dVarDown_V + dVarUp_V))
         pLeft_VF(iLogVar_V,iCell-1) = exp(pLeft_VF(iLogVar_V,iCell-1))
         pRight_VF(iLogVar_V,iCell-1) = exp(pRight_VF(iLogVar_V,iCell-1))
         call get_godunov_flux(pLeft_VF(:,iCell-1), pRight_VF(:,iCell-1), &
              Flux_VF(:,iCell-1), Cleft_F(iCell-1), Cright_F(iCell-1))
         pLeft_VF(:,iCell) = State_VC(:,iCell) + &
              (sign(0.25, dVarUp_V) + sign(0.25, dVarDown_V))* &
              min(Beta*abs(dVarUp_V), Beta*abs(dVarDown_V), &
              cThird*abs(dVarDown_V + 2*dVarUp_V))
         State_VC(iLogVar_V,iCell) = exp(State_VC(iLogVar_V,iCell))
      end do
      pLeft_VF(iLogVar_V,nCell-1) = exp(pLeft_VF(iLogVar_V,nCell-1))
      State_VC(iLogVar_V,nCell) = exp(State_VC(iLogVar_V,nCell))
      call get_godunov_flux(pLeft_VF(:,nCell-1), State_VC(:,nCell), &
           Flux_VF(:,nCell-1), Cleft_F(nCell-1), Cright_F(nCell-1))
      call get_godunov_flux(State_VC(:,nCell), State_VC(:,nCell), &
           Flux_VF(:,nCell), Cleft_F(nCell), Cright_F(nCell))
    end subroutine get_fluxes
    !==========================================================================
    subroutine get_time_step
      real :: CflLocal = 0.9, Dt_C(nCell)
      integer :: iCell
      !------------------------------------------------------------------------
      do iCell = 1, nCell/2
         ! No time offset in these cells
         Dt_C(iCell) = CflLocal*Ds/&
              max(Cright_F(iCell-1), -Cleft_F(iCell), 1e-30)
      end do
      do iCell = nCell/2+1, nCell
         ! Account for the time offset in these cells
         Dt_C(iCell) = CflLocal*Ds*&
              (1 - max(Cright_F(iCell-1), 0.0)*tOffsetPerR)/&
              max(Cright_F(iCell-1), -Cleft_F(iCell), 1e-30)
      end do
      Dt = minval(Dt_C)
    end subroutine get_time_step
    !==========================================================================
    subroutine advance_no_offset(StateIn_V, Source_V, StateOut_V)
      real, intent(in) :: StateIn_V(Rho_:P_), Source_V(Rho_:Energy_)
      real, intent(out):: StateOut_V(Rho_:P_)
      real :: RhoBar
      !------------------------------------------------------------------------
      StateOut_V = 0.0
      StateOut_V(Rho_) = StateIn_V(Rho_) + Source_V(Rho_)
      StateOut_V(Ux_)   = (StateIn_V(Rho_)*StateIn_V(Ux_) + Source_V(RhoUx_))/&
           StateOut_V(Rho_)
      StateOut_V(P_) = StateIn_V(P_) + (Gamma - 1)*&
           (Source_V(Energy_) + 0.5*(StateIn_V(Rho_)*StateIn_V(Ux_)**2 - &
           StateOut_V(Rho_)*StateOut_V(Ux_)**2) )
    end subroutine advance_no_offset
    !==========================================================================
    subroutine get_godunov_flux(pLeft_V, pRight_V, &
         Flux_V, Cleft, Cright)

      real,   intent(in) :: pLeft_V(Rho_:P_), pRight_V(Rho_:P_)
      real,  intent(out) :: Flux_V(Rho_:Energy_)
      real,  intent(out) :: Cleft, Cright

      real :: UnFace, RhoFace, pFace

      !------------------------------------------------------------------------
      RhoL = pLeft_V(Rho_ )
      UnL  = pLeft_V(Ux_)
      pL   = pLeft_V(P_)

      RhoR = pRight_V(Rho_ )
      UnR  = pRight_V(Ux_)
      pR   = pRight_V(P_)

      ! exact solver

      call exact_rs_pu_star
      call exact_rs_sample(0.0, RhoFace, UnFace, pFace)

      ! Store WL WR

      Cleft  = WL
      Cright = WR

      ! get the flux
      Flux_V(Rho_) = RhoFace * UnFace
      Flux_V(RhoUx_) = RhoFace*UnFace**2 + pFace
      Flux_V(RhoUx_+1:RhoUz_) = 0
      Flux_V(Energy_) = &
           (0.50*RhoFace*UnFace**2 + pFace*Gamma/(Gamma - 1))*UnFace
    end subroutine get_godunov_flux
    !==========================================================================
  end subroutine test_hydro_t
  !============================================================================
  subroutine get_mhd_rs
    use ModInterpolate, ONLY: interpolate_vector, linear

    real, parameter :: Gamma = 1.40
    real :: Time = 80.0, TimeLimit = 200.0, X10 = 400.0, TimeForecast = 40.0
    real :: tOffsetPerR = 0.125
    ! Coords
    real :: X_I(4000), Lambda_I(4000)
    real :: ExactX_VI(8,4000), ExactForecastX_VI(8,4000)
    real :: Time_I(1000), ExactByT400_I(1000), ExactByForecastT400_I(1000)
    integer, parameter :: Rho_ = 1, Ux_ = 2, Uy_ = 3, Uz_= 4, RhoUx_ = 2, &
         RhoUz_ = 4, Bx_ = 5, By_ = 6, Bz_ = 7, Energy_ = 8, P_ = 8
    integer, parameter :: nCell = 4000
    real, parameter :: Ds = 800.0/nCell
    real :: State_VC(Rho_:P_,nCell), StateOld_VC(Rho_:P_,nCell)
    integer, parameter :: iLogVar_V(2) = [Rho_, P_]
    real :: Flux_VF(Rho_:Energy_,0:nCell)
    real :: Cleft_F(0:nCell)
    real :: Cright_F(0:nCell)
    ! Interpolated left and right states at the face
    real    :: pLeft_VF(Rho_:P_,0:nCell)
    real    :: pRight_VF(Rho_:P_,0:nCell)
    real    :: dVarDown_V(Rho_:P_), dVarUp_V(Rho_:P_), Source_V(Rho_:Energy_)
    real    :: tSimulation, Dt
    integer :: nStep, i
    ! Misc:
    real :: Lambda
    real, parameter :: Normal_D(3) = [1.0, 0.0, 0.0]
    ! Initial state for forecast numerical test
    !--------------------------------------------------------------------------
    State_VC(Ux_:Uz_,:) = 0; State_VC(Bz_,:) = 0.0; State_VC(Bx_,:) = 0.75
    State_VC(Rho_,:nCell/2) = 1
    State_VC(Rho_,nCell/2+1:) = 0.125
    State_VC(By_,1:nCell/2)   = 1.0
    State_VC(By_,nCell/2+1:nCell)   = -1.0
    State_VC(P_,:nCell/2) = 1
    State_VC(P_,nCell/2+1:) = 0.1
    tSimulation = 0; Dt = 0; nStep = 0
    do
       if(tSimulation >= Time)EXIT
       Dt = min(Dt,1.001*(Time - tSimulation))
       tSimulation = tSimulation + Dt
       nStep = nStep + 1
       ! First Stage. Save old state
       StateOld_VC = State_VC
       call get_fluxes
       if(nStep==1) call get_time_step
       ! Advance to halfstep
       do i = 1, nCell
          call advance_no_offset(StateOld_VC(:,i)  , &
               0.5*Dt/Ds*(Flux_VF(:,i-1) - Flux_VF(:,i)), &
               State_VC(:,i))
       end do
       ! Second stage
       call get_fluxes
       ! Advance to full time step
       do i = 1,nCell
          call advance_no_offset(StateOld_VC(:,i)  , &
               Dt/Ds*(Flux_VF(:,i-1) - Flux_VF(:,i)), &
               State_VC(:,i))
       end do
       call get_time_step
       if(nStep==1)State_VC = StateOld_VC
       if(Dt<=0.0)then
          write(*,*)'Negative time step, reduce the time offset'
          stop
       end if
    end do
    call save_plot_file(NameFile='mhd_rs_X_4000.out', &
         TypeFileIn='ascii', TimeIn=tSimulation, nStepIn = nStep, &
         NameVarIn='X Rho ux uy uz bx by bz p', &
         CoordMinIn_D=[400.0/nCell], &
         CoordMaxIn_D=[800 - 400.0/nCell], &
         VarIn_VI = State_VC )
    do i = 1, nCell
       X_I(i) = (i - 0.50)*Ds
       Lambda_I(i) = (X_I(i) - 400)/Time
       ExactX_VI(:,i) = State_VC(:,i)
    end do
    ! Interpolate solution in the boosted frame
    do i = 1, nCell
       Lambda = (X_I(i) - 400)/TimeForecast
       ExactForecastX_VI(:,i) = interpolate_vector(ExactX_VI, 8, 1, &
            [1], [nCell],&
            x_D = [Lambda/(1 + max(Lambda, 0.0)*tOffsetPerR)],      &
            x1_I = Lambda_I, DoExtrapolate = .false.)
    end do
    call save_plot_file(NameFile='mhd_forecast_rs_X_4000.out', &
         TypeFileIn='ascii', TimeIn=tSimulation, nStepIn = nStep, &
         NameVarIn='X Rho ux uy uz bx by bz p', &
         CoordMinIn_D=[400.0/nCell], &
         CoordMaxIn_D=[800 - 400.0/nCell], &
         VarIn_VI = ExactForecastX_VI )
    ! Timeline
    do i = 1, 1000
       Time_I(i) = TimeLimit*(i-0.5)*1.0e-3
    end do
    ! Distribution over t:
    open(11,file='mhd_rs_T_200.out',status='replace')
    do i = 1, 1000
       Lambda = X10/Time_I(i)
       ExactByT400_I(i) = linear(ExactX_VI(By_,:), 1, nCell,        &
            x = Lambda,      &
            x_I = Lambda_I, DoExtrapolate = .false.)
       ExactByForecastT400_I(i) = linear(ExactX_VI(By_,:), 1, nCell,        &
            x = Lambda/(1 + max(Lambda, 0.0)*tOffsetPerR),      &
            x_I = Lambda_I, DoExtrapolate = .false.)
       write(11,*)Time_I(i), ExactByT400_I(i), ExactByForecastT400_I(i)
    end do
    close(11)
  contains
    !==========================================================================
    subroutine get_fluxes
      integer :: iCell
      real, parameter :: Beta = 1.0
      ! Left boundary:
      !------------------------------------------------------------------------
      call get_aw_flux_mhd(State_VC(:,1), State_VC(:,1), &
           Flux_VF(:,0), Cleft_F(0), Cright_F(0))
      State_VC(iLogVar_V,:) = log(State_VC(iLogVar_V,:))
      pLeft_VF(:,1) = State_VC(:,1)
      dVarUp_V = State_VC(:,2) - State_VC(:,1)
      do iCell = 2, nCell-1
         dVarDown_V = dVarUp_V
         ! Calculate Up slope
         dVarUp_V = State_VC(:,iCell+1) - State_VC(:,iCell)
         pRight_VF(:,iCell-1) = State_VC(:,iCell) - &
              (sign(0.25, dVarUp_V) + sign(0.25, dVarDown_V))* &
              min(Beta*abs(dVarUp_V), Beta*abs(dVarDown_V), &
              cThird*abs(2*dVarDown_V + dVarUp_V))
         pLeft_VF(iLogVar_V,iCell-1) = exp(pLeft_VF(iLogVar_V,iCell-1))
         pRight_VF(iLogVar_V,iCell-1) = exp(pRight_VF(iLogVar_V,iCell-1))
         call get_aw_flux_mhd(pLeft_VF(:,iCell-1), pRight_VF(:,iCell-1), &
              Flux_VF(:,iCell-1), Cleft_F(iCell-1), Cright_F(iCell-1))
         pLeft_VF(:,iCell) = State_VC(:,iCell) + &
              (sign(0.25, dVarUp_V) + sign(0.25, dVarDown_V))* &
              min(Beta*abs(dVarUp_V), Beta*abs(dVarDown_V), &
              cThird*abs(dVarDown_V + 2*dVarUp_V))
         State_VC(iLogVar_V,iCell) = exp(State_VC(iLogVar_V,iCell))
      end do
      pLeft_VF(iLogVar_V,nCell-1) = exp(pLeft_VF(iLogVar_V,nCell-1))
      State_VC(iLogVar_V,nCell) = exp(State_VC(iLogVar_V,nCell))
      call get_aw_flux_mhd(pLeft_VF(:,nCell-1), State_VC(:,nCell), &
           Flux_VF(:,nCell-1), Cleft_F(nCell-1), Cright_F(nCell-1))
      call get_aw_flux_mhd(State_VC(:,nCell), State_VC(:,nCell), &
           Flux_VF(:,nCell), Cleft_F(nCell), Cright_F(nCell))
    end subroutine get_fluxes
    !==========================================================================
    subroutine get_aw_flux_mhd(pLeft_V, pRight_V, &
         Flux_V, Cleft, Cright)

      real,   intent(in) :: pLeft_V(Rho_:P_), pRight_V(Rho_:P_)
      real,  intent(out) :: Flux_V(Rho_:Energy_)
      real,  intent(out) :: Cleft, Cright
      ! Misc:
      real :: ConsL_V(Rho_:Energy_), ConsR_V(Rho_:Energy_)
      real :: FluxL_V(Rho_:Energy_), FluxR_V(Rho_:Energy_)
      real :: UnL, FastL, UnR, FastR
      !------------------------------------------------------------------------

      ConsL_V = conserved_vars(pLeft_V)
      call get_flux_function(pLeft_V, ConsL_V, FluxL_V, UnL)
      FastL = fast_magnetosonic_speed(pLeft_V)

      ConsR_V = conserved_vars(pRight_V)
      call get_flux_function(pRight_V, ConsR_V, FluxR_V, UnR)
      FastR = fast_magnetosonic_speed(pRight_V)

      Cright = max(UnL + FastL, UnR + FastR, 0.0)
      Cleft  = min(UnL - FastL, UnR - FastR, 0.0)
      Flux_V = (Cright*FluxL_V - Cleft*FluxR_V + &
           Cleft*Cright*(ConsR_V - ConsL_V) )/(Cright - Cleft)

    end subroutine get_aw_flux_mhd
    !==========================================================================
    function conserved_vars(Prim_V)

      real :: conserved_vars(Rho_:P_)
      real, intent(in) :: Prim_V(Rho_:P_)
      !------------------------------------------------------------------------
      conserved_vars(Rho_) = Prim_V(Rho_)
      conserved_vars(RhoUx_:RhoUz_) = Prim_V(Rho_)*Prim_V(Ux_:Uz_)
      conserved_vars(Bx_:Bz_) = Prim_V(Bx_:Bz_)
      conserved_vars(Energy_) = 0.5*(sum(conserved_vars(RhoUx_:RhoUz_)*&
           Prim_V(Ux_:Uz_)) + sum(Prim_V(Bx_:Bz_)**2)) + Prim_V(P_)/(Gamma - 1)
    end function conserved_vars
    !==========================================================================
    subroutine get_flux_function(Prime_V, Cons_V, Flux_V, Un)

      real, intent(in)  :: Prime_V(Rho_:P_), Cons_V(Rho_:P_)
      real, intent(out) :: Flux_V(Rho_:P_)
      real, intent(out) :: Un
      ! Mic:
      real :: Bn, pTot
      !------------------------------------------------------------------------

      Un = sum(Normal_D*Prime_V(Ux_:Uz_))
      Bn = sum(Normal_D*Prime_V(Bx_:Bz_))
      pTot = Prime_V(P_) + 0.5*sum(Prime_V(Bx_:Bz_)**2)
      Flux_V = Un*Cons_V
      Flux_V(RhoUx_:RhoUz_) = Flux_V(RhoUx_:RhoUz_) + &
           pTot*Normal_D - Bn*Prime_V(Bx_:Bz_)
      Flux_V(Bx_:Bz_) = Flux_V(Bx_:Bz_) - Bn*Prime_V(Ux_:Uz_)
      Flux_V(Energy_) = Flux_V(Energy_) + Un*pTot - &
           Bn*sum(Prime_V(Ux_:Uz_)*Prime_V(Bx_:Bz_))
    end subroutine get_flux_function
    !==========================================================================
    real function fast_magnetosonic_speed(Prim_V)
      real, intent(in) :: Prim_V(Rho_:P_)
      ! Square of speed of sound, of Alfven wave speed, of it projection on n
      real :: Cs2, Va2, VaN2
      ! Their total; inverse density:
      real :: V2Tot, InvRho
      !------------------------------------------------------------------------
      InvRho = 1/Prim_V(Rho_)
      ! Squares of speeds:
      Cs2 = Gamma*Prim_V(P_)*InvRho
      Va2 = sum(Prim_V(Bx_:Bz_)**2)*InvRho
      VaN2 = sum(Prim_V(Bx_:Bz_)*Normal_D)**2*InvRho
      V2Tot = Cs2 + Va2
      fast_magnetosonic_speed = sqrt(0.50*(V2Tot + &
           sqrt(V2Tot**2 - 4*Cs2*VaN2)))
    end function fast_magnetosonic_speed
    !==========================================================================
    subroutine get_time_step
      real :: CflLocal = 0.85, Dt_C(nCell)
      integer :: iCell
      !------------------------------------------------------------------------
      do iCell = 1, nCell
         ! No time offset in these cells
         Dt_C(iCell) = CflLocal*Ds/&
              max(Cright_F(iCell-1), -Cleft_F(iCell), 1e-30)
      end do
      Dt = minval(Dt_C)
    end subroutine get_time_step
    !==========================================================================
    subroutine advance_no_offset(StateIn_V, Source_V, StateOut_V)
      real, intent(in) :: StateIn_V(Rho_:P_), Source_V(Rho_:Energy_)
      real, intent(out):: StateOut_V(Rho_:P_)
      real :: RhoBar
      !------------------------------------------------------------------------
      StateOut_V(Rho_) = StateIn_V(Rho_) + Source_V(Rho_)
      StateOut_V(Ux_:Uz_)   = (StateIn_V(Rho_)*StateIn_V(Ux_:Uz_) + &
           Source_V(RhoUx_:RhoUz_))/StateOut_V(Rho_)
      StateOut_V(Bx_:Bz_)   = StateIn_V(Bx_:Bz_) + Source_V(Bx_:Bz_)
      StateOut_V(P_) = StateIn_V(P_) + (Gamma - 1)*&
           (Source_V(Energy_) + 0.5*(StateIn_V(Rho_)*sum(StateIn_V(Ux_:Uz_)**2&
           ) - StateOut_V(Rho_)*sum(StateOut_V(Ux_:Uz_)**2)) + 0.50*(&
           sum(StateIn_V(Bx_:Bz_)**2) - sum(StateOut_V(Bx_:Bz_)**2) ) )
    end subroutine advance_no_offset
    !==========================================================================
  end subroutine get_mhd_rs
  !============================================================================
  subroutine test_mhd_x

    real, parameter :: Gamma = 1.40
    real :: TimeForecast = 40.0, tOffsetPerR = 0.125
    integer, parameter :: Rho_ = 1, Ux_ = 2, Uy_ = 3, Uz_= 4, RhoUx_ = 2, &
         RhoUz_ = 4, Bx_ = 5, By_ = 6, Bz_ = 7, Energy_ = 8, P_ = 8
    integer, parameter :: nCell = 800
    real, parameter :: Ds = 800.0/nCell
    real :: State_VC(Rho_:P_,nCell), StateOld_VC(Rho_:P_,nCell)
    integer, parameter :: iLogVar_V(2) = [Rho_, P_]
    real :: Flux_VF(Rho_:Energy_,0:nCell)
    real :: Cleft_F(0:nCell)
    real :: Cright_F(0:nCell)
    ! Interpolated left and right states at the face
    real    :: pLeft_VF(Rho_:P_,0:nCell)
    real    :: pRight_VF(Rho_:P_,0:nCell)
    real    :: dVarDown_V(Rho_:P_), dVarUp_V(Rho_:P_), Source_V(Rho_:Energy_)
    real    :: tSimulation, Dt
    integer :: nStep, i

    real, parameter :: Normal_D(3) = [1.0, 0.0, 0.0]

    ! Init boosted frame
    !--------------------------------------------------------------------------
    call set_param(tOffsetPerRin=tOffsetPerR, Gamma = Gamma)

    ! Initial state for forecast numerical test
    State_VC(Ux_:Uz_,:) = 0; State_VC(Bz_,:) = 0.0; State_VC(Bx_,:) = 0.75
    State_VC(Rho_,:nCell/2) = 1
    State_VC(Rho_,nCell/2+1:) = 0.125
    State_VC(By_,1:nCell/2)   = 1.0
    State_VC(By_,nCell/2+1:nCell)   = -1.0
    State_VC(P_,:nCell/2) = 1
    State_VC(P_,nCell/2+1:) = 0.1
    open(11,file='test_mhd_T.out',status='replace')
    tSimulation = 0; Dt = 0; nStep = 0
    do
       if(tSimulation >= TimeForecast)EXIT
       Dt = min(Dt,1.001*(TimeForecast - tSimulation))
       tSimulation = tSimulation + Dt
       nStep = nStep + 1
       write(*,*)'nStep, tSimulation, Dt=', nStep, tSimulation, Dt
       ! First Stage. Save old state
       StateOld_VC = State_VC
       call get_fluxes
       ! Advance to halfstep
       do i = 1, nCell/2
          call advance_no_offset(StateOld_VC(:,i)  , &
               0.5*Dt/Ds*(Flux_VF(:,i-1) - Flux_VF(:,i)), &
               State_VC(:,i))
       end do
       do i = nCell/2+1, nCell
          call update_mhd8wave_var(Old_V=StateOld_VC(:,i),&
               Source_V=0.5*Dt/Ds*(Flux_VF(:,i-1) - Flux_VF(:,i)),&
               Radial_D=[1.0, 0.0, 0.0],&
               New_V=State_VC(:,i))
       end do
       ! Second stage
       call get_fluxes
       ! Advance to full time step
       do i = 1,nCell/2
          call advance_no_offset(StateOld_VC(:,i)  , &
               Dt/Ds*(Flux_VF(:,i-1) - Flux_VF(:,i)), &
               State_VC(:,i))
       end do
       do i = nCell/2+1, nCell
          call update_mhd8wave_var(Old_V=StateOld_VC(:,i),&
               Source_V=Dt/Ds*(Flux_VF(:,i-1) - Flux_VF(:,i)),&
               Radial_D=[1.0, 0.0, 0.0],&
               New_V=State_VC(:,i))
       end do
       call get_time_step
       if(Dt<=0.0)then
          write(*,*)'Negative time step, reduce the time offset'
          stop
       end if
    end do
    call save_plot_file(NameFile='test_mhd_x.out', &
         TypeFileIn='ascii', TimeIn=tSimulation, nStepIn = nStep, &
         NameVarIn='X Rho ux uy uz bx by bz p', &
         CoordMinIn_D=[400.0/nCell], &
         CoordMaxIn_D=[800.0 - 400.0/nCell], &
         VarIn_VI = State_VC )
  contains
    !==========================================================================
    subroutine get_fluxes
      integer :: iCell
      real, parameter :: Beta = 1.0
      ! Left boundary:
      !------------------------------------------------------------------------
      call get_aw_flux_mhd(State_VC(:,1), State_VC(:,1), &
           Flux_VF(:,0), Cleft_F(0), Cright_F(0))
      State_VC(iLogVar_V,:) = log(State_VC(iLogVar_V,:))
      pLeft_VF(:,1) = State_VC(:,1)
      dVarUp_V = State_VC(:,2) - State_VC(:,1)
      do iCell = 2, nCell-1
         dVarDown_V = dVarUp_V
         ! Calculate Up slope
         dVarUp_V = State_VC(:,iCell+1) - State_VC(:,iCell)
         pRight_VF(:,iCell-1) = State_VC(:,iCell) - &
              (sign(0.25, dVarUp_V) + sign(0.25, dVarDown_V))* &
              min(Beta*abs(dVarUp_V), Beta*abs(dVarDown_V), &
              cThird*abs(2*dVarDown_V + dVarUp_V))
         pLeft_VF(iLogVar_V,iCell-1) = exp(pLeft_VF(iLogVar_V,iCell-1))
         pRight_VF(iLogVar_V,iCell-1) = exp(pRight_VF(iLogVar_V,iCell-1))
         call get_aw_flux_mhd(pLeft_VF(:,iCell-1), pRight_VF(:,iCell-1), &
              Flux_VF(:,iCell-1), Cleft_F(iCell-1), Cright_F(iCell-1))
         pLeft_VF(:,iCell) = State_VC(:,iCell) + &
              (sign(0.25, dVarUp_V) + sign(0.25, dVarDown_V))* &
              min(Beta*abs(dVarUp_V), Beta*abs(dVarDown_V), &
              cThird*abs(dVarDown_V + 2*dVarUp_V))
         State_VC(iLogVar_V,iCell) = exp(State_VC(iLogVar_V,iCell))
      end do
      pLeft_VF(iLogVar_V,nCell-1) = exp(pLeft_VF(iLogVar_V,nCell-1))
      State_VC(iLogVar_V,nCell) = exp(State_VC(iLogVar_V,nCell))
      call get_aw_flux_mhd(pLeft_VF(:,nCell-1), State_VC(:,nCell), &
           Flux_VF(:,nCell-1), Cleft_F(nCell-1), Cright_F(nCell-1))
      call get_aw_flux_mhd(State_VC(:,nCell), State_VC(:,nCell), &
           Flux_VF(:,nCell), Cleft_F(nCell), Cright_F(nCell))
    end subroutine get_fluxes
    !==========================================================================
    subroutine get_aw_flux_mhd(pLeft_V, pRight_V, &
         Flux_V, Cleft, Cright)

      real,   intent(in) :: pLeft_V(Rho_:P_), pRight_V(Rho_:P_)
      real,  intent(out) :: Flux_V(Rho_:Energy_)
      real,  intent(out) :: Cleft, Cright
      ! Misc:
      real :: ConsL_V(Rho_:Energy_), ConsR_V(Rho_:Energy_)
      real :: FluxL_V(Rho_:Energy_), FluxR_V(Rho_:Energy_)
      real :: UnL, FastL, UnR, FastR
      !------------------------------------------------------------------------

      ConsL_V = conserved_vars(pLeft_V)
      call get_flux_function(pLeft_V, ConsL_V, FluxL_V, UnL)
      FastL = fast_magnetosonic_speed(pLeft_V)

      ConsR_V = conserved_vars(pRight_V)
      call get_flux_function(pRight_V, ConsR_V, FluxR_V, UnR)
      FastR = fast_magnetosonic_speed(pRight_V)

      Cright = max(UnL + FastL, UnR + FastR, 0.0)
      Cleft  = min(UnL - FastL, UnR - FastR, 0.0)
      Flux_V = (Cright*FluxL_V - Cleft*FluxR_V + &
           Cleft*Cright*(ConsR_V - ConsL_V) )/(Cright - Cleft)

    end subroutine get_aw_flux_mhd
    !==========================================================================
    function conserved_vars(Prim_V)

      real :: conserved_vars(Rho_:P_)
      real, intent(in) :: Prim_V(Rho_:P_)
      !------------------------------------------------------------------------
      conserved_vars(Rho_) = Prim_V(Rho_)
      conserved_vars(RhoUx_:RhoUz_) = Prim_V(Rho_)*Prim_V(Ux_:Uz_)
      conserved_vars(Bx_:Bz_) = Prim_V(Bx_:Bz_)
      conserved_vars(Energy_) = 0.5*(sum(conserved_vars(RhoUx_:RhoUz_)*&
           Prim_V(Ux_:Uz_)) + sum(Prim_V(Bx_:Bz_)**2)) + Prim_V(P_)/(Gamma - 1)
    end function conserved_vars
    !==========================================================================
    subroutine get_flux_function(Prime_V, Cons_V, Flux_V, Un)

      real, intent(in)  :: Prime_V(Rho_:P_), Cons_V(Rho_:P_)
      real, intent(out) :: Flux_V(Rho_:P_)
      real, intent(out) :: Un
      ! Mic:
      real :: Bn, pTot
      !------------------------------------------------------------------------

      Un = sum(Normal_D*Prime_V(Ux_:Uz_))
      Bn = sum(Normal_D*Prime_V(Bx_:Bz_))
      pTot = Prime_V(P_) + 0.5*sum(Prime_V(Bx_:Bz_)**2)
      Flux_V = Un*Cons_V
      Flux_V(RhoUx_:RhoUz_) = Flux_V(RhoUx_:RhoUz_) + &
           pTot*Normal_D - Bn*Prime_V(Bx_:Bz_)
      Flux_V(Bx_:Bz_) = Flux_V(Bx_:Bz_) - Bn*Prime_V(Ux_:Uz_)
      Flux_V(Energy_) = Flux_V(Energy_) + Un*pTot - &
           Bn*sum(Prime_V(Ux_:Uz_)*Prime_V(Bx_:Bz_))
    end subroutine get_flux_function
    !==========================================================================
    real function fast_magnetosonic_speed(Prim_V)
      real, intent(in) :: Prim_V(Rho_:P_)
      ! Square of speed of sound, of Alfven wave speed, of it projection on n
      real :: Cs2, Va2, VaN2
      ! Their total; inverse density:
      real :: V2Tot, InvRho
      !------------------------------------------------------------------------
      InvRho = 1/Prim_V(Rho_)
      ! Squares of speeds:
      Cs2 = Gamma*Prim_V(P_)*InvRho
      Va2 = sum(Prim_V(Bx_:Bz_)**2)*InvRho
      VaN2 = sum(Prim_V(Bx_:Bz_)*Normal_D)**2*InvRho
      V2Tot = Cs2 + Va2
      fast_magnetosonic_speed = sqrt(0.50*(V2Tot + &
           sqrt(V2Tot**2 - 4*Cs2*VaN2)))
    end function fast_magnetosonic_speed
    !==========================================================================
    subroutine get_time_step
      real :: CflLocal = 0.85, Dt_C(nCell)
      integer :: iCell
      !------------------------------------------------------------------------
      do iCell = 1, nCell
         ! No time offset in these cells
         Dt_C(iCell) = CflLocal*Ds/&
              max(Cright_F(iCell-1), -Cleft_F(iCell), 1e-30)
      end do
      do iCell = nCell/2 + 1, nCell
         Dt_C(iCell) = Dt_C(iCell)*(1 - tOffsetPerR*&
              (sum(State_VC(Ux_:Uz_,iCell)*Normal_D) + &
              fast_magnetosonic_speed(State_VC(:,iCell))))
      end do
      Dt = minval(Dt_C)
    end subroutine get_time_step
    !==========================================================================
    subroutine advance_no_offset(StateIn_V, Source_V, StateOut_V)
      real, intent(in) :: StateIn_V(Rho_:P_), Source_V(Rho_:Energy_)
      real, intent(out):: StateOut_V(Rho_:P_)
      real :: RhoBar
      !------------------------------------------------------------------------
      StateOut_V(Rho_) = StateIn_V(Rho_) + Source_V(Rho_)
      StateOut_V(Ux_:Uz_)   = (StateIn_V(Rho_)*StateIn_V(Ux_:Uz_) + &
           Source_V(RhoUx_:RhoUz_))/StateOut_V(Rho_)
      StateOut_V(Bx_:Bz_)   = StateIn_V(Bx_:Bz_) + Source_V(Bx_:Bz_)
      StateOut_V(P_) = StateIn_V(P_) + (Gamma - 1)*&
           (Source_V(Energy_) + 0.5*(StateIn_V(Rho_)*sum(StateIn_V(Ux_:Uz_)**2&
           ) - StateOut_V(Rho_)*sum(StateOut_V(Ux_:Uz_)**2)) + 0.50*(&
           sum(StateIn_V(Bx_:Bz_)**2) - sum(StateOut_V(Bx_:Bz_)**2) ) )
    end subroutine advance_no_offset
    !==========================================================================
  end subroutine test_mhd_x
  !============================================================================
  subroutine test_mhd_t

    real, parameter :: Gamma = 1.40
    real :: TimeLimit = 200.0, tOffsetPerR = 0.125
    integer, parameter :: Rho_ = 1, Ux_ = 2, Uy_ = 3, Uz_= 4, RhoUx_ = 2, &
         RhoUz_ = 4, Bx_ = 5, By_ = 6, Bz_ = 7, Energy_ = 8, P_ = 8
    integer, parameter :: nCell = 1600
    real, parameter :: Ds = 1600.0/nCell
    real :: State_VC(Rho_:P_,nCell), StateOld_VC(Rho_:P_,nCell)
    integer, parameter :: iLogVar_V(2) = [Rho_, P_]
    real :: Flux_VF(Rho_:Energy_,0:nCell)
    real :: Cleft_F(0:nCell)
    real :: Cright_F(0:nCell)
    ! Interpolated left and right states at the face
    real    :: pLeft_VF(Rho_:P_,0:nCell)
    real    :: pRight_VF(Rho_:P_,0:nCell)
    real    :: dVarDown_V(Rho_:P_), dVarUp_V(Rho_:P_), Source_V(Rho_:Energy_)
    real    :: tSimulation, Dt
    integer :: nStep, i

    real, parameter :: Normal_D(3) = [1.0, 0.0, 0.0]

    ! Init boosted frame
    !--------------------------------------------------------------------------
    call set_param(tOffsetPerRin=tOffsetPerR, Gamma = Gamma)

    ! Initial state for forecast numerical test
    State_VC(Ux_:Uz_,:) = 0; State_VC(Bz_,:) = 0.0; State_VC(Bx_,:) = 0.75
    State_VC(Rho_,:nCell/2) = 1
    State_VC(Rho_,nCell/2+1:) = 0.125
    State_VC(By_,1:nCell/2)   = 1.0
    State_VC(By_,nCell/2+1:nCell)   = -1.0
    State_VC(P_,:nCell/2) = 1
    State_VC(P_,nCell/2+1:) = 0.1
    open(11,file='test_mhd_t.out',status='replace')
    tSimulation = 0; Dt = 0; nStep = 0
    do
       if(tSimulation >= 0.50*TimeLimit)EXIT
       Dt = min(Dt,1.001*(0.50*TimeLimit - tSimulation))
       tSimulation = tSimulation + Dt
       nStep = nStep + 1
       ! First Stage. Save old state
       StateOld_VC = State_VC
       call get_fluxes
       ! Advance to halfstep
       do i = 1, nCell/2
          call advance_no_offset(StateOld_VC(:,i)  , &
               0.5*Dt/Ds*(Flux_VF(:,i-1) - Flux_VF(:,i)), &
               State_VC(:,i))
       end do
       do i = nCell/2+1, nCell
          call update_mhd8wave_var(Old_V=StateOld_VC(:,i),&
               Source_V=0.5*Dt/Ds*(Flux_VF(:,i-1) - Flux_VF(:,i)),&
               Radial_D=[1.0, 0.0, 0.0],&
               New_V=State_VC(:,i))
       end do
       ! Second stage
       call get_fluxes
       ! Advance to full time step
       do i = 1,nCell/2
          call advance_no_offset(StateOld_VC(:,i)  , &
               Dt/Ds*(Flux_VF(:,i-1) - Flux_VF(:,i)), &
               State_VC(:,i))
       end do
       do i = nCell/2+1, nCell
          call update_mhd8wave_var(Old_V=StateOld_VC(:,i),&
               Source_V=Dt/Ds*(Flux_VF(:,i-1) - Flux_VF(:,i)),&
               Radial_D=[1.0, 0.0, 0.0],&
               New_V=State_VC(:,i))
       end do
       call get_time_step
       if(Dt<=0.0)then
          write(*,*)'Negative time step, reduce the time offset'
          stop
       end if
       write(11,*)tSimulation, 0.5*(State_VC(By_,3*nCell/4) + &
            State_VC(By_,3*nCell/4+1))
    end do
    close(11)
  contains
    !==========================================================================
    subroutine get_fluxes
      integer :: iCell
      real, parameter :: Beta = 1.0
      ! Left boundary:
      !------------------------------------------------------------------------
      call get_aw_flux_mhd(State_VC(:,1), State_VC(:,1), &
           Flux_VF(:,0), Cleft_F(0), Cright_F(0))
      State_VC(iLogVar_V,:) = log(State_VC(iLogVar_V,:))
      pLeft_VF(:,1) = State_VC(:,1)
      dVarUp_V = State_VC(:,2) - State_VC(:,1)
      do iCell = 2, nCell-1
         dVarDown_V = dVarUp_V
         ! Calculate Up slope
         dVarUp_V = State_VC(:,iCell+1) - State_VC(:,iCell)
         pRight_VF(:,iCell-1) = State_VC(:,iCell) - &
              (sign(0.25, dVarUp_V) + sign(0.25, dVarDown_V))* &
              min(Beta*abs(dVarUp_V), Beta*abs(dVarDown_V), &
              cThird*abs(2*dVarDown_V + dVarUp_V))
         pLeft_VF(iLogVar_V,iCell-1) = exp(pLeft_VF(iLogVar_V,iCell-1))
         pRight_VF(iLogVar_V,iCell-1) = exp(pRight_VF(iLogVar_V,iCell-1))
         call get_aw_flux_mhd(pLeft_VF(:,iCell-1), pRight_VF(:,iCell-1), &
              Flux_VF(:,iCell-1), Cleft_F(iCell-1), Cright_F(iCell-1))
         pLeft_VF(:,iCell) = State_VC(:,iCell) + &
              (sign(0.25, dVarUp_V) + sign(0.25, dVarDown_V))* &
              min(Beta*abs(dVarUp_V), Beta*abs(dVarDown_V), &
              cThird*abs(dVarDown_V + 2*dVarUp_V))
         State_VC(iLogVar_V,iCell) = exp(State_VC(iLogVar_V,iCell))
      end do
      pLeft_VF(iLogVar_V,nCell-1) = exp(pLeft_VF(iLogVar_V,nCell-1))
      State_VC(iLogVar_V,nCell) = exp(State_VC(iLogVar_V,nCell))
      call get_aw_flux_mhd(pLeft_VF(:,nCell-1), State_VC(:,nCell), &
           Flux_VF(:,nCell-1), Cleft_F(nCell-1), Cright_F(nCell-1))
      call get_aw_flux_mhd(State_VC(:,nCell), State_VC(:,nCell), &
           Flux_VF(:,nCell), Cleft_F(nCell), Cright_F(nCell))
    end subroutine get_fluxes
    !==========================================================================
    subroutine get_aw_flux_mhd(pLeft_V, pRight_V, &
         Flux_V, Cleft, Cright)

      real,   intent(in) :: pLeft_V(Rho_:P_), pRight_V(Rho_:P_)
      real,  intent(out) :: Flux_V(Rho_:Energy_)
      real,  intent(out) :: Cleft, Cright
      ! Misc:
      real :: ConsL_V(Rho_:Energy_), ConsR_V(Rho_:Energy_)
      real :: FluxL_V(Rho_:Energy_), FluxR_V(Rho_:Energy_)
      real :: UnL, FastL, UnR, FastR
      !------------------------------------------------------------------------

      ConsL_V = conserved_vars(pLeft_V)
      call get_flux_function(pLeft_V, ConsL_V, FluxL_V, UnL)
      FastL = fast_magnetosonic_speed(pLeft_V)

      ConsR_V = conserved_vars(pRight_V)
      call get_flux_function(pRight_V, ConsR_V, FluxR_V, UnR)
      FastR = fast_magnetosonic_speed(pRight_V)

      Cright = max(UnL + FastL, UnR + FastR, 0.0)
      Cleft  = min(UnL - FastL, UnR - FastR, 0.0)
      Flux_V = (Cright*FluxL_V - Cleft*FluxR_V + &
           Cleft*Cright*(ConsR_V - ConsL_V) )/(Cright - Cleft)

    end subroutine get_aw_flux_mhd
    !==========================================================================
    function conserved_vars(Prim_V)

      real :: conserved_vars(Rho_:P_)
      real, intent(in) :: Prim_V(Rho_:P_)
      !------------------------------------------------------------------------
      conserved_vars(Rho_) = Prim_V(Rho_)
      conserved_vars(RhoUx_:RhoUz_) = Prim_V(Rho_)*Prim_V(Ux_:Uz_)
      conserved_vars(Bx_:Bz_) = Prim_V(Bx_:Bz_)
      conserved_vars(Energy_) = 0.5*(sum(conserved_vars(RhoUx_:RhoUz_)*&
           Prim_V(Ux_:Uz_)) + sum(Prim_V(Bx_:Bz_)**2)) + Prim_V(P_)/(Gamma - 1)
    end function conserved_vars
    !==========================================================================
    subroutine get_flux_function(Prime_V, Cons_V, Flux_V, Un)

      real, intent(in)  :: Prime_V(Rho_:P_), Cons_V(Rho_:P_)
      real, intent(out) :: Flux_V(Rho_:P_)
      real, intent(out) :: Un
      ! Mic:
      real :: Bn, pTot
      !------------------------------------------------------------------------

      Un = sum(Normal_D*Prime_V(Ux_:Uz_))
      Bn = sum(Normal_D*Prime_V(Bx_:Bz_))
      pTot = Prime_V(P_) + 0.5*sum(Prime_V(Bx_:Bz_)**2)
      Flux_V = Un*Cons_V
      Flux_V(RhoUx_:RhoUz_) = Flux_V(RhoUx_:RhoUz_) + &
           pTot*Normal_D - Bn*Prime_V(Bx_:Bz_)
      Flux_V(Bx_:Bz_) = Flux_V(Bx_:Bz_) - Bn*Prime_V(Ux_:Uz_)
      Flux_V(Energy_) = Flux_V(Energy_) + Un*pTot - &
           Bn*sum(Prime_V(Ux_:Uz_)*Prime_V(Bx_:Bz_))
    end subroutine get_flux_function
    !==========================================================================
    real function fast_magnetosonic_speed(Prim_V)
      real, intent(in) :: Prim_V(Rho_:P_)
      ! Square of speed of sound, of Alfven wave speed, of it projection on n
      real :: Cs2, Va2, VaN2
      ! Their total; inverse density:
      real :: V2Tot, InvRho
      !------------------------------------------------------------------------
      InvRho = 1/Prim_V(Rho_)
      ! Squares of speeds:
      Cs2 = Gamma*Prim_V(P_)*InvRho
      Va2 = sum(Prim_V(Bx_:Bz_)**2)*InvRho
      VaN2 = sum(Prim_V(Bx_:Bz_)*Normal_D)**2*InvRho
      V2Tot = Cs2 + Va2
      fast_magnetosonic_speed = sqrt(0.50*(V2Tot + &
           sqrt(V2Tot**2 - 4*Cs2*VaN2)))
    end function fast_magnetosonic_speed
    !==========================================================================
    subroutine get_time_step
      real :: CflLocal = 0.85, Dt_C(nCell)
      integer :: iCell
      !------------------------------------------------------------------------
      do iCell = 1, nCell
         ! No time offset in these cells
         Dt_C(iCell) = CflLocal*Ds/&
              max(Cright_F(iCell-1), -Cleft_F(iCell), 1e-30)
      end do
      do iCell = nCell/2 + 1, nCell
         Dt_C(iCell) = Dt_C(iCell)*(1 - tOffsetPerR*&
              (sum(State_VC(Ux_:Uz_,iCell)*Normal_D) + &
              fast_magnetosonic_speed(State_VC(:,iCell))))
      end do
      Dt = minval(Dt_C)
    end subroutine get_time_step
    !==========================================================================
    subroutine advance_no_offset(StateIn_V, Source_V, StateOut_V)
      real, intent(in) :: StateIn_V(Rho_:P_), Source_V(Rho_:Energy_)
      real, intent(out):: StateOut_V(Rho_:P_)
      real :: RhoBar
      !------------------------------------------------------------------------
      StateOut_V(Rho_) = StateIn_V(Rho_) + Source_V(Rho_)
      StateOut_V(Ux_:Uz_)   = (StateIn_V(Rho_)*StateIn_V(Ux_:Uz_) + &
           Source_V(RhoUx_:RhoUz_))/StateOut_V(Rho_)
      StateOut_V(Bx_:Bz_)   = StateIn_V(Bx_:Bz_) + Source_V(Bx_:Bz_)
      StateOut_V(P_) = StateIn_V(P_) + (Gamma - 1)*&
           (Source_V(Energy_) + 0.5*(StateIn_V(Rho_)*sum(StateIn_V(Ux_:Uz_)**2&
           ) - StateOut_V(Rho_)*sum(StateOut_V(Ux_:Uz_)**2)) + 0.50*(&
           sum(StateIn_V(Bx_:Bz_)**2) - sum(StateOut_V(Bx_:Bz_)**2) ) )
    end subroutine advance_no_offset
    !==========================================================================
  end subroutine test_mhd_t
  !============================================================================
end module ModTestBoostedFrame
!==============================================================================
program test_program
  use ModTestBoostedFrame, ONLY: get_hydro_rs, test_hydro_x, test_hydro_t, &
       get_mhd_rs, test_mhd_x, test_mhd_t
  implicit none
  !----------------------------------------------------------------------------
  write(*,*)' Start exact rs and forecast'
  call get_hydro_rs

  write(*,*)' Start test_hydro_x'
  call test_hydro_x

  write(*,*)' Start test_hydro_t'
  call test_hydro_t

  write(*,*)' Start get_mhd_rs'
  call get_mhd_rs

  write(*,*)' Start test_mhd_x'
  call test_mhd_x

  write(*,*)' Start test_mhd_t'
  call test_mhd_t

end program test_program
!==============================================================================
