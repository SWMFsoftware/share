!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!===========================TESTS============================================
!Solve energetic particle transport via the new numerical scheme
!Dimensionless parameters: Energy: Mev, Speed: c Mass: Mev/c^2 Length: solar radius
!In the code, label Q: s_L  P: \mu  R:P^3/3  NLoop: t
module ModTestPoissonBracket
  use ModPoissonBracket, ONLY: explicit!, explicit3 
  use ModUtilities,      ONLY: CON_stop
  use ModNumConst,       ONLY: cTwoPi
  use ModPlotFile,       ONLY: save_plot_file
  use ModConst
  implicit none
contains
  
    subroutine test_poisson_bracket(tFinal)
    real, intent(in) :: tFinal
    !Misc:
    integer, parameter::  nQ = 30,  nP = 360
    integer           ::  iQ, iP, iStep
    real, parameter   :: qMax = 10.0, qMin = 0.01 
    real, parameter   :: DeltaPhi = cTwoPi/nP 
    real :: MomentumRatio, MomentumMin, MomentumMax
    real :: VDF_G(-1:nQ+2, -1:nP+2), Volume_G(0:nQ+1, 0:nP+1)
    real :: Hamiltonian_N(-1:nQ+1, -1:nP+1)
    real :: LogMomentum_I(0:nQ+1), Momentum2_I(-1:nQ+1)
    real :: Time, Dt, Source_C(nQ,nP)
    !---------------------
    MomentumRatio = exp(log(qMax/qMin)/nQ)
    MomentumMax =  qMin/MomentumRatio
    Momentum2_I(-1) = MomentumMax**2
    Hamiltonian_N(-1,:) = Hamiltonian(Momentum2_I(-1))
    do iQ = 0, nQ+1
       MomentumMin = MomentumMax
       MomentumMax = MomentumMin*MomentumRatio
       Momentum2_I( iQ) = MomentumMax**2
       Hamiltonian_N(iQ,:) = Hamiltonian(Momentum2_I(iQ))
       Volume_G(iQ,:) = 0.5*DeltaPhi*&
            (Momentum2_I(iQ) - Momentum2_I(iQ-1))
       LogMomentum_I(iQ)   = 0.50*log10(MomentumMin*MomentumMax)       
    end do
    VDF_G = 0.0; VDF_G(1:nQ, 172:189) = 1.0; Source_C = 0.0
    !\
    ! Compiutation
    !/
    Time = 0.0; iStep = 0
    do 
       call explicit(nQ, nP, VDF_G, Volume_G, Source_C, &
            Hamiltonian_N,   &
            CFLIn=0.99, DtOut = Dt)
       iStep = iStep +1
       if(Time + Dt >= tFinal)then
          call explicit(nQ, nP, VDF_G, Volume_G, Source_C, &
               Hamiltonian_N,   &
               DtIn = tFinal - Time)
          VDF_G(1:nQ, 1:nP) = VDF_G(1:nQ, 1:nP) + Source_C
          EXIT
       else
          Time = Time + Dt
          VDF_G(1:nQ, 1:nP) = VDF_G(1:nQ, 1:nP) + Source_C
          VDF_G(1:nQ,-1:0 ) = VDF_G(1:nQ, nP-1:nP)
          VDF_G(1:nQ, nP+1:nP+2 ) = VDF_G(1:nQ, 1:2)
       end if
    end do
    call save_plot_file(NameFile='test_poisson.out', &
         TypeFileIn='ascii', TimeIn=tFinal, nStepIn = iStep, &
         NameVarIn='Phi LogMomentum VDF'  , &
         CoordMinIn_D=(/0.50,      LogMomentum_I(1 )/),&
         CoordMaxIn_D=(/nP - 0.50, LogMomentum_I(nQ)/),&
         StringFormatIn = '(4F10.3)'            ,&
         Coord2In_I = LogMomentum_I(1:nQ)           ,&
         VarIn_II = transpose(VDF_G(1:nQ,1:nP)))
    contains
      real function Hamiltonian(P2)
        real, intent(in) :: P2 ! momentum squared
        Hamiltonian = sqrt(1.0 + P2)
      end function Hamiltonian
  end subroutine test_poisson_bracket  
!Solve energetic particle transport via the new numerical scheme with multi Poisson bracket
!{f,H_1}_t,p^3/3 + {f,H_2}_s_L,\mu + {f,H_3}_p^3/3,\mu = 0
!Dimensionless parameters: Energy: Mev, Speed: c Mass: Mev/c^2 Length: c*s
!In the code, label Q: s_L  P: \mu  R: ln(P^3/3)  NLoop: t
  subroutine test_multipoisson_bracket(TimeOut)
    real,intent(in) :: TimeOut
    !Misc:
    real,parameter :: tEachFile=100.0                                   !time range of each input file
    real           :: Time=0.0, Dt=0.0                                  !simulation time, time step
    real,parameter :: ProtonmassInGeV=cProtonmass*clightspeed**2/&
         (1.0e9*cElectroncharge)!proton mass in unit of GeV/c^2
    real,parameter :: ParticleRange1=2.0, ParticleRange2=4.0            !parameter of the place of the initial particles
    integer,parameter :: nQ = 130,  nP = 20,  nR = 20                   !number of grid of s_L, \mu, ln(p^3/3) axis      
    integer,parameter :: nVar = 14                                      !number of vatiales in each line of the input files               
    integer,parameter :: x_ = 2,y_ = 3,z_ = 4,ux_ = 7,uy_ = 8,uz_ = 9,&
                        Bx_ = 10, By_= 11, Bz_ = 12                     !Position of the variables in the input line
    integer,parameter :: nStride = 25                                   !stride of the data reading
    integer,parameter :: iLagrIDSkipped = 20                            !Lagr index that is gaped for field line at each time
    real :: DeltaMu = 2.0/nP                                            !\delta\mu for each grid along the \mu axis
    real :: ParticleEnergyMin, ParticleEnergyMax                        !Min and Max value for particle energy
    real :: LnP3min, LnP3max                                            !Ln(P^3/3)_min, Ln(P^3/3)_max
    real :: DeltaP3_I(nR), DeltaLnP3                                    ! \delta(P^3/3)  \deltaLn(P^3/3)
    real :: Volume_G(0:nQ+1, 0:nP+1, 0:nR+1)                            !control volume 
    real :: VDF_G(-1:nQ+2, -1:nP+2, -1:nR+2)                            !vdf
    real :: VDFOutput_II(1:nQ, 1:nR)                                    !vdf output
    real :: Hamiltonian2_N(-1:nQ+1,-1:nP+1,0:nR+1)                      !the Poisson bracket with regard to the first and second vars 
    real :: Hamiltonian3_N(0:nQ+1,-1:nP+1,-1:nR+1)                      !consider the case when there are more than one Poisson bracket
    real :: DeltaHamiltonian1_N(0:nQ+1,0:nP+1,-1:nR+1)                  !and when there is a Poisson bracket with respect to the time
    real :: Source_C(nQ,nP,nR)                                          !Source at each time step
    real :: RawData1_II(nQ,nVar), RawData2_II(nQ,nVar)                  !the raw data that is imported from MHD
    real :: GapData                                                     !Read the redundant data
    real :: DeltaSOverB_C(1:nQ), DDeltaSOverBDt_C(1:nQ), InvB_C(1:nQ), &
         bDuDt_C(1:nQ), DLnBDeltaSSquaredDt_C(1:nQ)                     !The data calculated directly from RawData
    real :: tOutput                                                     !Time of each output file
    real :: HelioDist_I(nQ)                                             !Heliocentric distance, for visualization
    real :: Energy_I(nR)                                                !Energy in unit of KeV, for output
    integer :: m, iQ, iP, iR, iStep=0, iStride                          !Loop integers
    integer :: nOutputFile = 2                                          !The number of output files
    integer :: IndexOfFile                                              !Number Index of File
    character(LEN=3) :: NameFileSuffix                                  !The number index of output files in character
    
    
       
    !\
    ! Import data from files
    !/
    open(unit=13,file='FieldLineOld.in')
    open(unit=14,file='FieldLineNew.in')
    
    !\ 
    ! Get the lines that we need 
    !/
    GapData = 0.0
    do while(nint(GapData) /= iLagrIDSkipped); read(13,*)GapData; end do
    do iQ = 1, nQ
       read(13,*)(RawData1_II(iQ,m), m = 1,nVar)
       do iStride = 1, nStride -1
          read(13,*)
       end do
    end do
    close(13)
    
    GapData = 0.0
    do while(nint(GapData) /= iLagrIDSkipped); read(14,*)GapData; end do
    do iQ = 1, nQ
       read(14,*)(RawData2_II(iQ,m), m = 1,nVar)
       do iStride = 1, nStride -1
          read(14,*)
       end do
    end do
    close(14)
    
    !\
    ! Normalize process: This is because the unit we use in this code is different from the unit of input data 
    !/
    RawData1_II(:,x_ :z_ ) = RawData1_II(:,x_ :z_ )*rSun/clightspeed
    RawData1_II(:,ux_:uz_) = RawData1_II(:,ux_:uz_)/clightspeed
    RawData2_II(:,x_ :z_ ) = RawData2_II(:,x_ :z_ )*rSun/clightspeed
    RawData2_II(:,ux_:uz_) = RawData2_II(:,ux_:uz_)/clightspeed
    
    
    !\
    ! Initialize the Ln(p^3/3) axis
    !/
    ParticleEnergyMin = 1.0e-3       !Min energy in unit of Gev
    ParticleEnergyMax = 1.0          !Max energy in unit of Gev
    !\
    !give the min value of Ln(P^3/3). We set the minimum of energy 
    !to be 1MeV, then Ln(p^3/3) can be calculated
    !/
    LnP3min = log(sqrt((ParticleEnergyMin + ProtonmassInGeV)**2 - ProtonmassInGeV**2)**3/3)
    !\ 
    !give the max value of Ln(P^3/3). We set the maximum of energy to be 1GeV, 
    !then Ln(p^3/3) can be calculated 
    !/
    LnP3max = log(sqrt((ParticleEnergyMax + ProtonmassInGeV)**2 - ProtonmassInGeV**2)**3/3) 
    !\
    !get  \DeltaLn(P^3/3)
    !/
    DeltaLnP3 = (LnP3max - LnP3min)/nR 
    do iR=1,nR
       !get \deltaP^3/3 for each grid
       DeltaP3_I(iR) = exp(LnP3min + iR*DeltaLnP3) - exp(LnP3min + (iR - 1)*DeltaLnP3)
    end do
    
    !\                                       
    ! Initialize the distribution function VDF_G
    !/
    VDF_G = 1.0e-7
    !Assign the VDF at positive \mu
    do iR = 1,nR
       do iQ = 1,nQ
          if (iQ <= ParticleRange1) then
             !Assume the momentum of particls to be power distribution: f ~ p^-5
             VDF_G(iQ,nP/2:nP,iR) = max(0.1*exp(-5.0/3.0*(LnP3min+iR*DeltaLnP3)),1.0e-7)
          elseif (iQ <= ParticleRange2)then
             !using a linear slope for the middle part
             VDF_G(iQ,nP/2:nP,iR) = max(0.1*exp(-5.0/3.0*(LnP3min+iR*DeltaLnP3))*&
                  (ParticleRange2 - real(iQ))/ &
                  (ParticleRange2 - ParticleRange1),1.0e-7)
          end if
       end do
    end do
    !\
    !boundary conditions
    !/
    !use reflective boundary condition for \mu axis
    VDF_G(:,    0, :) = VDF_G(:,    1, :)
    VDF_G(:,   -1, :) = VDF_G(:,    2, :)
    VDF_G(:, nP+1, :) = VDF_G(:,   nP, :)
    VDF_G(:, nP+2, :) = VDF_G(:, nP-1, :)
    !use outflow boundary condition for s_L axis
    VDF_G(-1,   :, :) = VDF_G(1,    :, :)
    VDF_G(0,    :, :) = VDF_G(1,    :, :)
    VDF_G(nQ+1, :, :) = VDF_G(nQ,   :, :)
    VDF_G(nQ+2, :, :) = VDF_G(nQ,   :, :)
    !use outflow boundary condition for Ln(p^3/3) axis
    VDF_G(:, :,   -1) = VDF_G(:, :,    1)
    VDF_G(:, :,    0) = VDF_G(:, :,    1)
    VDF_G(:, :, nR+1) = VDF_G(:, :,   nR)
    VDF_G(:, :, nR+2) = VDF_G(:, :,   nR)
    
    !\
    ! Initialize output data
    !/
    Source_C = 0.0
    
    
    !\ 
    ! Calculate the Lagrangian particle positions
    !/   
    do iQ = 1, nQ
       HelioDist_I(iQ) = sqrt(RawData1_II(iQ,x_)**2 + RawData1_II(iQ,y_)**2 + RawData1_II(iQ,z_)**2 )*clightspeed/rSun
    end do
    
    !\
    ! Calculate the position for ln(p^3/3) axis: note that the energy is in the unit of KeV now
    !/
    do iR=1,nR
       Energy_I(iR) = 6.0 + log10(sqrt((3.0*exp(LnP3min + (iR - 0.50)*DeltaLnP3))**(2.0/3.0) + ProtonmassInGeV**2) &
            - ProtonmassInGeV)
    end do
    
    !\
    ! Calculate the VDF for output: here, we integrate over the \mu axis
    !/                   
    VDFOutput_II = 0.0
    do iP = 1, nP
       VDFOutput_II = VDFOutput_II + VDF_G(1:nQ, iP, 1:nR)
    end do
    VDFOutput_II = VDFOutput_II*DeltaMu
    
    !\
    ! Give the initial outpur of VDF
    !/
    call save_plot_file(NameFile='test_multipoisson_000.out', &
         TypeFileIn='ascii', TimeIn=Time, nStepIn = iStep,    &
         NameVarIn='sL logE VDF'  ,                           &
         CoordMinIn_D=(/HelioDist_I(1) ,      Energy_I(1)/),  &
         CoordMaxIn_D=(/HelioDist_I(nQ),      Energy_I(nR)/), &
         StringFormatIn = '(3F15.7)',                         &  
         Coord1In_I = HelioDist_I,                            &
         Coord2In_I = Energy_I,                               &
         VarIn_II = log10(VDFOutput_II(1:nQ,1:nR)))
    
    
    !\
    ! ***Start the main simulation loop***
    ! In the loop, we output several snapshot of the distribution function VDF 
    !/
    do IndexofFile=1,nOutputFile
       tOutput = TimeOut/real(nOutputFile)*real(IndexofFile)
       !\
       ! In this loop, we first update the RawData and calculate 3 Hamiltonian functions and total volume 
       ! Then we use subroutine explicit3 to calculate Source_C and update VDF
       !/
       do
          
          call calc_initial_data(&
               RawData1_II,     &
               RawData2_II,     &
               Time,            &
               DeltaSOverB_C = DeltaSOverB_C,       &
               DDeltaSOverBDt_C = DDeltaSOverBDt_C, &
               InvB_C = InvB_C,                     &
               bDuDt_C = bDuDt_C,                   &
               DLnBDeltaSSquaredDt_C = DLnBDeltaSSquaredDt_C)!import MHD data at all the calculation time
          DeltaHamiltonian1_N = 0.0
          Hamiltonian2_N = 0.0
          Hamiltonian3_N = 0.0
          !calculate Hamiltonian function which is in the time-dependent poisson bracket
          call calc_hamiltonian_1(DDeltaSOverBDt_C,   &
                                  DeltaHamiltonian1_N,&
                                  DeltaP3_I,          &
                                  DeltaLnP3,          &     
                                  LnP3min,LnP3max)
          !calculate the second Hamiltonian function 
          call calc_hamiltonian_2(InvB_C,             &
                                  Hamiltonian2_N     ,&
                                  DeltaP3_I,DeltaLnP3,&
                                  LnP3min,LnP3max)
          !\
          !calculate the third Hamiltonian function
          !/ 
          call calc_hamiltonian_3(DeltaSOverB_C,      &
                                  bDuDt_C,            &
                                  DLnBDeltaSSquaredDt_C,&
                                  Hamiltonian3_N     ,&
                                  DeltaP3_I          ,&
                                  DeltaLnP3,LnP3min,LnP3max)
          
          !\
          !calculate the total volume
          !/
          do iQ=1,nQ
             do iR=1,nR
                Volume_G(iQ,0:nP+1,iR)=DeltaSOverB_C(iQ)*DeltaMu*DeltaP3_I(iR)
             end do
          end do
          !boundary conditions for total volume
          Volume_G(0,   :,   :)=Volume_G(1 ,:, :)
          Volume_G(nQ+1,:,   :)=Volume_G(nQ,:, :)
          Volume_G(:,   :,   0)=Volume_G(: ,:, 1)
          Volume_G(:,   :,nR+1)=Volume_G(: ,:,nR)
          
          
          !\
          ! Begin the calculation of Source_C
          !/
          call explicit(nQ, nP, nR, VDF_G,  Volume_G, Source_C,  &
               Hamiltonian2_N,                                   &
               Hamiltonian23_N = -Hamiltonian3_N,                &
               dHamiltonian03_FZ = DeltaHamiltonian1_N,          &
               DtOut = Dt, CFLIn=0.99)
          
          iStep = iStep + 1
          if(Time + Dt >= tOutput)then
             call explicit(nQ, nP, nR, VDF_G,  Volume_G, Source_C, &
                  Hamiltonian2_N,                                  &
                  Hamiltonian23_N = -Hamiltonian3_N,               &
                  dHamiltonian03_FZ = DeltaHamiltonian1_N,         &
                  DtIn = tOutput - Time)
             VDF_G(1:nQ, 1:nP, 1:nR) = VDF_G(1:nQ, 1:nP, 1:nR) + Source_C/Volume_G(1:nQ,1:nP,1:nR)
             Time = tOutput
             Source_C = 0.0
             EXIT
          else
             Time = Time + Dt
             !\
             ! Update DeltaSOverB_C for the calculation of volume
             !/
             call calc_initial_data(RawData1_II,RawData2_II,Time,DeltaSOverB_C = DeltaSOverB_C)
             !\                                                                                                                              
             ! Calculate the total volume                                        
             !/
             do iQ=1,nQ
                do iR=1,nR
                   Volume_G(iQ,0:nP+1,iR)=DeltaSOverB_C(iQ)*DeltaMu*DeltaP3_I(iR)
                end do
             end do
             !boundary conditions for total volume
             Volume_G(0,   :,   :)=Volume_G(1 ,:, :)
             Volume_G(nQ+1,:,   :)=Volume_G(nQ,:, :)
             Volume_G(:,   :,   0)=Volume_G(: ,:, 1)
             Volume_G(:,   :,nR+1)=Volume_G(: ,:,nR)
             
             VDF_G(1:nQ, 1:nP, 1:nR) = VDF_G(1:nQ, 1:nP, 1:nR) + Source_C/Volume_G(1:nQ,1:nP,1:nR)
             Source_C=0.0
             !use reflective boundary condition for \mu axis
             VDF_G(:,    0, :) = VDF_G(:,    1, :)
             VDF_G(:,   -1, :) = VDF_G(:,    2, :)
             VDF_G(:, nP+1, :) = VDF_G(:,   nP, :)
             VDF_G(:, nP+2, :) = VDF_G(:, nP-1, :)
             !use outflow boundary condition for s_L axis
             VDF_G(-1,   :, :) = VDF_G(1,    :, :)
             VDF_G(0,    :, :) = VDF_G(1,    :, :)
             VDF_G(nQ+1, :, :) = VDF_G(nQ,   :, :)
             VDF_G(nQ+2, :, :) = VDF_G(nQ,   :, :)
             !use outflow boundary condition for ln(p^3/3) axis
             VDF_G(:, :,   -1) = VDF_G(:, :,    1)
             VDF_G(:, :,    0) = VDF_G(:, :,    1)
             VDF_G(:, :, nR+1) = VDF_G(:, :,   nR)
             VDF_G(:, :, nR+2) = VDF_G(:, :,   nR)
          end if
       end do
       
       !\
       ! To avoid some mistakes when taking log10 to VDF_G, we here set a min value for VDF_G
       !/
       VDF_G = max(VDF_G, 1.0e-7)
       
       !\
       ! Set the name of files, here the largest number of files is 999
       !/
       write(NameFileSuffix,'(I3.3)')IndexOfFile       
       
       !\
       ! Calculate the VDF for output: here, we integrate over the \mu axis
       !/                   
       VDFOutput_II = 0.0
       do iP = 1, nP
          VDFOutput_II = VDFOutput_II + VDF_G(1:nQ, iP, 1:nR)
       end do
       VDFOutput_II = VDFOutput_II*DeltaMu
       
       !\
       ! Save the files
       !/
       call save_plot_file(NameFile='test_multipoisson_'//NameFileSuffix//'.out', &
         TypeFileIn='ascii', TimeIn=Time, nStepIn = iStep,       &
         NameVarIn='sL logE VDF'  ,                              &
         CoordMinIn_D=(/HelioDist_I(1) ,      Energy_I(1)/),     &
         CoordMaxIn_D=(/HelioDist_I(nQ),      Energy_I(nR)/),    &
         StringFormatIn = '(3F15.7)',                            &
         Coord1In_I = HelioDist_I,                               &
         Coord2In_I = Energy_I,                                  &
         VarIn_II = log10(VDFOutput_II(1:nQ,1:nR)))
   
    end do
    
  contains

    
    subroutine calc_initial_data(&
         RawData1_II, RawData2_II,Time,DeltaSOverB_C,DDeltaSOverBDt_C,InvB_C,bDuDt_C,DLnBDeltaSSquaredDt_C)!calculate data from the input files
      real, optional, intent(inout) :: DeltaSOverB_C(1:nQ), DDeltaSOverBDt_C(1:nQ), InvB_C(1:nQ), bDuDt_C(1:nQ), DLnBdeltaSSquaredDt_C(1:nQ)
      real, intent(in)    :: RawData1_II(nQ,nVar), RawData2_II(nQ,nVar)
      real, intent(in)    :: Time
      real :: DeltaSOverBOld_C(1:nQ), InvBOld_C(1:nQ), LnBDeltaSSquaredOld_C(1:nQ)!\deltas/b, 1/(2B), ln(B\deltas^2) at Old time
      real :: DeltaSOverBNew_C(1:nQ), InvBNew_C(1:nQ), LnBDeltaSSquaredNew_C(1:nQ)!\deltas/b, 1/(2B), ln(B\deltas^2) at New time 
      real :: B_C(1:nQ,3)!magnetic field B at cell center
      real :: MidPoint_ID(nQ-1,x_:z_), DeltaS_I(nQ)!midpoint for to consecutive points, \deltas
      integer :: iQ!loop integer
      
      
      !\
      ! Calculate values at Old time first
      !/
      !calculate midpoints
      do iQ=1,nQ-1
         MidPoint_ID(iQ,:) = (RawData1_II(iQ+1,x_:z_) + RawData1_II(iQ,x_:z_))*0.5
      end do
      !calculate deltas
      do iQ=2,nQ-1
         DeltaS_I(iQ) = sqrt((MidPoint_ID(iQ,x_) - MidPoint_ID(iQ-1,x_))**2 + &
              (MidPoint_ID(iQ,y_) - MidPoint_ID(iQ-1,y_))**2 + (MidPoint_ID(iQ,z_) &
              - MidPoint_ID(iQ-1,z_))**2)
      end do
      !Linear interpolate the delats such that there will be nQ deltas
      DeltaS_I(1)  = 2*sqrt((MidPoint_ID(1,   x_) - RawData1_II(1, x_))**2 +       &
           (MidPoint_ID(1,   y_) - RawData1_II(1, y_))**2 + (MidPoint_ID(1,   z_)&
           - RawData1_II(1, z_))**2)
      DeltaS_I(nQ) = 2*sqrt((MidPoint_ID(nQ-1,x_) - RawData1_II(nQ,x_))**2 + &
           (MidPoint_ID(nQ-1,y_) - RawData1_II(nQ,y_))**2 + (MidPoint_ID(nQ-1,z_) &
           - RawData1_II(nQ,z_))**2)
      
      !\
      ! calculate part of hamiltonian function (1/2B) and \DeltaS/B at grid center
      !/
      do iQ=1,nQ
         InvBOld_C(iQ) = 1.0/(2.0*sqrt(RawData1_II(iQ, Bx_)**2 + &
              RawData1_II(iQ, By_)**2 + RawData1_II(iQ, Bz_)**2))
         DeltaSOverBOld_C(iQ) = DeltaS_I(iQ)*2.0*InvBOld_C(iQ)
         LnBDeltaSSquaredOld_C(iQ) = log(sqrt(RawData1_II(iQ, Bx_)**2 + &
              RawData1_II(iQ, By_)**2 + RawData1_II(iQ, Bz_)**2)*DeltaS_I(iQ)**2)
      end do
      
      
      !\
      ! Calculate values at New time
      !/
      do iQ=1,nQ-1
         MidPoint_ID(iQ,:) = (RawData2_II(iQ+1,x_:z_) + RawData2_II(iQ,x_:z_))*0.5
      end do
      !calculate deltas
      do iQ=2,nQ-1
         DeltaS_I(iQ) = sqrt((MidPoint_ID(iQ,x_) - MidPoint_ID(iQ-1,x_))**2 +   &
         (MidPoint_ID(iQ,y_) - MidPoint_ID(iQ-1,y_))**2 + (MidPoint_ID(iQ,z_) - &
              MidPoint_ID(iQ-1,z_))**2)
      end do
      !Linear interpolate the delats such that there will be nQ deltas
      DeltaS_I(1)  = 2*sqrt((MidPoint_ID(1,   x_) - RawData2_II(1, x_))**2 + &
           (MidPoint_ID(1,   y_) - RawData2_II(1, y_))**2 + (MidPoint_ID(1,   z_) - &
           RawData2_II(1, z_))**2)
      DeltaS_I(nQ) = 2*sqrt((MidPoint_ID(nQ-1,x_) - RawData2_II(nQ,x_))**2 + &
           (MidPoint_ID(nQ-1,y_) - RawData2_II(nQ,y_))**2 + (MidPoint_ID(nQ-1,z_) - &
           RawData2_II(nQ,z_))**2)
      
      !\
      ! calculate part of hamiltonian function (1/2B) and \DeltaS/B at grid center
      !/
      do iQ=1,nQ
         InvBNew_C(iQ) = 1.0/(2.0*sqrt(RawData2_II(iQ, Bx_)**2 + RawData2_II(iQ, By_)**2 + &
              RawData2_II(iQ, Bz_)**2))
         DeltaSOverBNew_C(iQ) = DeltaS_I(iQ)*2.0*InvBNew_C(iQ)
         LnBDeltaSSquaredNew_C(iQ) = log(sqrt(RawData2_II(iQ, Bx_)**2 + &
              RawData2_II(iQ, By_)**2 + RawData2_II(iQ, Bz_)**2)*DeltaS_I(iQ)**2)
      end do
      
      
      !\
      ! Calculate values for Current time
      !/
      !Calculate B at cell center
      if (present(bDuDt_C)) B_C = RawData1_II(:,Bx_:Bz_) + (RawData2_II(:,Bx_:Bz_) - &
           RawData1_II(:,Bx_:Bz_))/tEachFile*Time
      !Calculate  (1/2B) 
      if (present(InvB_C)) InvB_C = InvBOld_C + (InvBNew_C - InvBOld_C)/tEachFile*Time
      !Calculate \deltas/B at current time
      if (present(DeltaSOverB_C)) DeltaSOverB_C = DeltaSOverBOld_C + (DeltaSOverBNew_C &
           - DeltaSOverBOld_C)/tEachFile*Time
      !Calculate \deltas/B at the next time
      if (present(DdeltaSOverBDt_C)) DDeltaSOverBDt_C = (DeltaSOverBNew_C &
           - DeltaSOverBOld_C)/tEachFile
      !Calculate Dln(B\deltas^2)/Dt
      if (present(DLnBdeltaSSquaredDt_C)) DLnBdeltaSSquaredDt_C = (LnBDeltaSSquaredNew_C&
           - LnBDeltaSSquaredOld_C)/tEachfile
      !Calculate b*Du/Dt
      if (present(bDuDt_C)) then
         do iQ=1,nQ
            bDuDt_C(iQ) = (B_C(iQ,1)*(RawData2_II(iQ,Bx_) - RawData1_II(iQ,Bx_)) + &
                 B_C(iQ,2) * (RawData2_II(iQ,By_) - RawData1_II(iQ,By_)) + &
                 B_C(iQ,3) * (RawData2_II(iQ,Bz_) - RawData1_II(iQ,Bz_)))  &
                 *2.0*InvB_C(iQ)/tEachFile
         end do
      end if
      !\
      !Here we use linear interpolation to get the data of each time step from every to consecutive files
      !/
      
    end subroutine calc_initial_data
    
    
    subroutine calc_hamiltonian_1(DdeltaSOverBDt_C,DeltaHamiltonian1_N,DeltaP3_I,DeltaLnP3,LnP3min,LnP3max)
      !\
      !calculate the first Hamiltonian function with time: p^3/3*\deltas/B\
      !/
      real, intent(out) :: DeltaHamiltonian1_N(0:nQ+1,0:nP+1,-1:nR+1) 
      real, intent(in) :: DeltaP3_I(nR), DeltaLnP3, LnP3min, LnP3max
      real, intent(in) :: DDeltaSOverBDt_C(1:nQ)
      integer :: iP,iR!loop integers
      
      !\
      ! calculate the first hamiltonian function
      !/
      do iP=0,nP+1
         do iR=0,nR
            DeltaHamiltonian1_N(1:nQ,iP,iR) = -exp(LnP3min + iR*DeltaLnP3)*DDeltaSOverBDt_C(1:nQ)
         end do
      end do
      
      !\
      ! calculate the Hamiltonian function that is to be used in actual calculation, which is \deltaH\tuta!!
      !/
      
      DeltaHamiltonian1_N(1:nQ,0:nP+1,0:nR) = DeltaHamiltonian1_N(1:nQ,0:nP+1,0:nR)*DeltaMu
      
      !\
      ! Boundary condition of Hamiltonian function
      !/
      DeltaHamiltonian1_N(0   ,:,   :) = DeltaHamiltonian1_N(1 ,:, :)
      DeltaHamiltonian1_N(nQ+1,:,   :) = DeltaHamiltonian1_N(nQ,:, :)
      DeltaHamiltonian1_N(:   ,:,  -1) = DeltaHamiltonian1_N(: ,:, 0)
      DeltaHamiltonian1_N(:   ,:,nR+1) = DeltaHamiltonian1_N(: ,:,nR)
      
    end subroutine calc_hamiltonian_1
    
    !calculate the second Hamiltonian function at each fixed time: (1-mu^2)v/2B 
    subroutine calc_hamiltonian_2(InvB_C,Hamiltonian2_N,DeltaP3_I,DeltaLnP3,LnP3min,LnP3max)
      real, intent(out) :: Hamiltonian2_N(-1:nQ+1, -1:nP+1, 0:nR+1)                                                                                                                                        
      real, intent(in) :: DeltaP3_I(nR),DeltaLnP3,LnP3min,LnP3max
      real, intent(in) :: InvB_C(nQ)
      real :: InvB_F(0:nQ)!intermediate matrix or arrays  
      real :: Momentum_I(1:nR)!momentum at each center of the cell
      integer :: iQ,iP,iR!loop integers
      
      !\
      ! calculate the value of part of the hamiltonian function (1/2B) on the boundary of grid
      !/                     
      do iQ=1,nQ-1
         InvB_F(iQ) = (InvB_C(iQ+1) + InvB_C(iQ))*0.5
      end do
      InvB_F(0)  = InvB_C(1 ) - (InvB_C(2 ) - InvB_C(1   ))*0.5
      InvB_F(nQ) = InvB_C(nQ) + (InvB_C(nQ) - InvB_C(nQ-1))*0.5
      
      !\
      ! calculate the real hamiltonian function, mutiply (1-mu^2) and v (nocite that v is a 
      ! function of P^3/3)!!!
      !/                   
      do iR=1,nR
         Momentum_I(iR) = (3.0*exp(LnP3min + (iR - 0.5)*DeltaLnP3))**(1.0/3.0)
      end do
      do iQ=0,nQ
         do iP=0,nP
            do iR=1,nR
               !Consider the law of relativity, v=1/sqrt(1+m^2*c^2/p^2), we can calculate v as 
               !a function of p. Note that light speed is the unit of speed in our code, 
               !so we do not need to multiply a c^2 in the following expression
               Hamiltonian2_N(iQ,iP,iR) = InvB_F(iQ)*(1.0 - (-1.0 + real(iP)*DeltaMu)**2)&
                    /sqrt(1 + ProtonmassInGeV**2/Momentum_I(iR)**2)
            end do
         end do
      end do
      
      !\
      ! calculate the Hamiltonian function that is to be used in actual calculation, which is \deltaH\tuta!!
      !/
      do iR=1,nR
         Hamiltonian2_N(0:nQ,0:nP,iR) = Hamiltonian2_N(0:nQ,0:nP,iR)*DeltaP3_I(iR)
      end do

      !\
      ! Boundary condition of Hamiltonian function
      !/
      Hamiltonian2_N(-1  ,:   ,:   ) = Hamiltonian2_N(0 ,:   ,: )
      Hamiltonian2_N(nQ+1,:   ,:   ) = Hamiltonian2_N(nQ,:   ,: )
      Hamiltonian2_N(:   ,-1  ,:   ) = Hamiltonian2_N(: ,1   ,: )
      Hamiltonian2_N(:   ,nP+1,:   ) = Hamiltonian2_N(: ,nP-1,: )
      Hamiltonian2_N(:   ,:   ,0   ) = Hamiltonian2_N(: ,:   ,1 )
      Hamiltonian2_N(:   ,:   ,nR+1) = Hamiltonian2_N(: ,:   ,nR)

      !\
      ! because of the requirement of equation, the Hamiltonian with a minus poisson bracket must be its opposite
      !/
      Hamiltonian2_N = -Hamiltonian2_N       
    end subroutine calc_hamiltonian_2
    
    
    subroutine calc_hamiltonian_3(DeltaSOverB_C,bDuDt_C,DLnBDeltaSSquaredDt_C,Hamiltonian3_N,DeltaP3_I,DeltaLnP3,LnP3min,LnP3max)
      real, intent(out) :: Hamiltonian3_N(0:nQ+1,-1:nP+1,-1:nR+1)
      real, intent(in) :: DeltaP3_I(nR), DeltaLnP3, LnP3min, LnP3max
      real, intent(in) :: DeltaSOverB_C(1:nQ), bDuDt_C(1:nQ), DLnBDeltaSSquaredDt_C(1:nQ)
      real :: Momentum_I(0:nR)!momentum at each cell face
      integer :: iQ,iP,iR!loop integers
      
      
      !\
      !calculate the momentum for convenience
      !/
      do iR=0,nR
         Momentum_I(iR)=(3.0*exp(LnP3min + iR*DeltaLnP3))**(1.0/3.0)
      end do

      !\
      !calculate the third hamiltonian function
      !/
      do iQ = 1, nQ
         do iP = 0, nP
            do iR = 0, nR
               Hamiltonian3_N(iQ,iP,iR) = (1.0 - (-1.0 + real(iP)*DeltaMu)**2)/2.0*&
                    ((-1.0 + real(iP)*DeltaMu)*exp(LnP3min + iR*DeltaLnP3)*DLnBDeltaSSquaredDt_C(iQ) &
                    + ProtonmassInGeV*Momentum_I(iR)**2*bDuDt_C(iQ))
            end do
         end do
      end do
      
      !\
      !calculate the Hamiltonian function that is to be used in actual calculation, which is \deltaH\tuta!!
      !/
      do iQ = 1, nQ
         Hamiltonian3_N(iQ,0:nP,0:nR) = Hamiltonian3_N(iQ,0:nP,0:nR)*DeltaSOverB_C(iQ)
      end do
      
      !\
      ! Boundary condition of Hamiltonian function
      !/
      Hamiltonian3_N(0   ,:   ,:   ) = Hamiltonian3_N(1 ,:   ,: )
      Hamiltonian3_N(nQ+1,:   ,:   ) = Hamiltonian3_N(nQ,:   ,: )
      Hamiltonian3_N(:   ,-1  ,:   ) = Hamiltonian3_N(: ,1   ,: )
      Hamiltonian3_N(:   ,nP+1,:   ) = Hamiltonian3_N(: ,nP-1,: )
      Hamiltonian3_N(:   ,:   ,-1  ) = Hamiltonian3_N(: ,:   ,0 )
      Hamiltonian3_N(:   ,:   ,nR+1) = Hamiltonian3_N(: ,:   ,nR)
    end subroutine calc_hamiltonian_3
    
    
  end subroutine test_multipoisson_bracket
end module ModTestPoissonBracket
!=============================================================================
!Solve energetic particle transport via the new numerical scheme
!Dimensionless parameters: Energy: Mev, Speed: c Mass: Mev/c^2 Length: solar radius
!In the code, label Q: s_L  P: \mu  R:P^3/3  NLoop: t
module ModDiffusion
  !Solve the diffusion term in the Parker equation
  !Revision history:
  !Prototype:Sp/FLAMPA/src/SP_ModDiffusion.f90 
  !Adapted for the use in MFLAMPA (Dist_I is an input paramater, 
  !fixed contributions to M_I in the end points)-D.Borovikov, 2017
  !Updated (identation, comments):  I.Sokolov, Dec.17, 2017
  implicit none
  PRIVATE
  public:: advance_diffusion1,tridiag
contains
  !=========================================================================
  ! This routine solves the diffusion equation:
  !         f_t-D_outer(D_inner*f_x)_x=0,
  ! with zero Neumann boundary condition. The solution is advanced in time
  ! using fully implicit scheme.
  subroutine advance_diffusion1(Dt,n,Dist_I,F_I,DOuter_I,DInner_I)
    use ModNumConst, ONLY: cTiny

    real,   intent(in   ):: Dt     !Time step                       
    integer,intent(in   ):: n      !Number of meshes along the x-coordinate
    real,   intent(in   ):: Dist_I(n) !Distance to the next mesh
    real,   intent(inout):: F_I(n) !In:sol.to be advanced; Out:advanced sol
    !Laplace multiplier and diffusion coefficient.
    real,   intent(in   ):: DOuter_I(n), DInner_I(n)
   
    !Mesh spacing and face spacing.
    real                 :: DsMesh_I(2:n), DsFace_I(2:n-1)
    !Main, upper, and lower diagonals.
    real, dimension(n)   :: Main_I,Upper_I,Lower_I, R_I 
    integer:: i
    real:: Aux1,Aux2
    !-----------------------------------------------------------------
    !\
    !In M-FLAMPA D_I(i) is the distance between meshes i   and i+1
    !while DsMesh_I(i) is the distance between centers of meshes  
    !i-1 and i. Therefore,
    DsMesh_I = 0.0
    DsFace_I = 0.0
    Lower_I = 0.0
    Main_I = 0.0
    Upper_I = 0.0
    R_I = 0.0
    do i=2,n
       DsMesh_I(i) = max(Dist_I(i-1),cTiny)
    end do
    !/
    !\
    ! Within the framework of finite volume method, the cell 
    ! volume is used, which is proportional to  the distance between 
    ! the faces bounding the volume with an index, i, which is half of
    ! sum of distance between meshes i-1 and i (i.e. D_I(i-1) and that 
    ! between meshes i and i+1 (which is D_I(i)):  
    do i=2,n-1
       DsFace_I(i) = max(0.5*(Dist_I(i) + Dist_I(i-1)),cTiny)
    end do
    
    !/
    ! In flux coordinates, the control volume associated with the 
    ! given cell has a cross-section equal to (Magnetic Flux)/B, 
    ! where the flux is a constant along the magnetic field line, 
    ! set to one hereafter. Therefore, the diffusion equation has 
    ! a following form:
    ! (DsFace_i/B_i)(f^(n+1) - f^n) = Flux_(i-1/2) - Flux_(i+1/2),
    ! where the particle density flux should be multiplied by the
    ! cross-section area too (magnetic flux factor is one!):  
    !  Flux_(i-1/2) = (diffusion coefficient)_(i-1/2)/B_(i-1/2)*&
    !                 (f^(n+1)_(i-1) - f^(n+1)_i),
    !  The multiplier, DsFace_i/B_i, is denoted as DsFace_i/DOuter_i
    !  The face-centered combination, 
    !\
    ! f^(n+1)_i-Dt*DOuter_I/DsFace_I*(&
    !     DInner_(i+1/2)*(f^(n+1)_(i+1)-f^(n+1)_i)/DsMesh_(i+1)-&
    !     DInner_(i-1/2)*(f^(n+1)_i -f^(n+1)_(i-1)/DsMesh_(i ))=f^n_i
    !/
    Main_I = 1.0
    !\
    ! For i=1:
    !/
    Aux1 = Dt*DOuter_I(1)*0.50*(DInner_I(1)+DInner_I(2))/&
         DsMesh_I(2)**2
    Main_I( 1) = Main_I(1)+Aux1
    Upper_I(1) = -Aux1
    !\
    ! For i=2,n-1:
    !/
    do i=2,n-1
       Aux1 = Dt*DOuter_I(i)*0.50*(DInner_I(i  ) + DInner_I(i+1))/&
            (DsMesh_I(i+1)*DsFace_I(i))
       Aux2 = Dt*DOuter_I(i)*0.50*(DInner_I(i-1) + DInner_I(i  ))/&
            (DsMesh_I(i  )*DsFace_I(i))
       Main_I(i)  = Main_I(i) + Aux1 + Aux2
       Upper_I(i) = -Aux1
       Lower_I(i) = -Aux2
    end do
    !\
    ! For i=n:
    !/
    Aux2 = Dt*DOuter_I(n)*0.50*(DInner_I(n-1) + DInner_I(n))/&
         DsMesh_I(n)**2
    Main_I( n) = Main_I(n) + Aux2
    Lower_I(n) = -Aux2
    !\
    ! Update the solution from f^(n) to f^(n+1):
    !/
    R_I = F_I
    call tridiag(n,Lower_I,Main_I,Upper_I,R_I,F_I)
  end subroutine advance_diffusion1
  !================================================================
  ! This routine solves three-diagonal system of equations: 
  !  ||m_1 u_1  0....        || ||w_1|| ||r_1||
  !  ||l_2 m_2 u_2...        || ||w_2|| ||r_2||
  !  || 0  l_3 m_3 u_3       ||.||w_3||=||r_3||
  !  ||...                   || ||...|| ||...||
  !  ||.............0 l_n m_n|| ||w_n|| ||r_n||
  ! From: Numerical Recipes, Chapter 2.6, p.40.                               
  subroutine tridiag(n, L_I, M_I, U_I, R_I, W_I)
    !\
    ! input parameters
    !/
    integer,            intent(in):: n
    real, dimension(n), intent(in):: L_I, M_I ,U_I ,R_I
    !\
    ! Output parameters
    !/
    real, intent(out):: W_I(n)
    !\
    ! Misc
    !/
    integer:: j
    real:: Aux,Aux_I(2:n)
    !-------------------!
    if (M_I(1)==0.0) then
       write(*,*)' Error in tridiag: M_I(1)=0'
       stop
    end if
    Aux = M_I(1)
    W_I(1) = R_I(1)/Aux
    do j=2,n
       Aux_I(j) = U_I(j-1)/Aux
       Aux = M_I(j)-L_I(j)*Aux_I(j)
       if (Aux.eq.0.0) then
          write(*,*)'M_I(j), L_I(j), Aux_I(j) = ',&
               M_I(j),L_I(j),Aux_I(j)
          write(*,*)'  For j=',j
          write(*,*)'Tridiag failed'
          stop
       end if
       W_I(j) = (R_I(j)-L_I(j)*W_I(j-1))/Aux
    end do
    do j=n-1,1,-1
       W_I(j) = W_I(j)-Aux_I(j+1)*W_I(j+1)
    end do
  end subroutine tridiag
end module ModDiffusion
!=======================
module ModTestPoissonBracketAndScatter
  use ModPoissonBracket!, ONLY: explicit3, explicit => explicit3 
  use ModUtilities,      ONLY: CON_stop
  use ModNumConst,       ONLY: cTwoPi
  use ModPlotFile,       ONLY: save_plot_file
  use ModConst
  implicit none
contains  
!Solve energetic particle transport via the new numerical scheme with multi Poisson bracket
!{f,H_1}_t,p^3/3 + {f,H_2}_s_L,\mu + {f,H_3}_p^3/3,\mu = 0
!Dimensionless parameters: Energy: Mev, Speed: c Mass: Mev/c^2 Length: c*s
!In the code, label Q: s_L  P: \mu  R: ln(P^3/3)  NLoop: t
!this code consider only the most simple equation
  subroutine test_multipoisson_bracket(TimeOut)
    real,intent(in) :: TimeOut
    !Misc:
    integer,parameter :: Diffuornot=0                                                   !whether using dissusive approximation or not
    real           :: DtFixed = 0.03                                                    !optional fixed time step
    real           :: Time=0.0, Dt=0.0                                                  !simulation time, time step
    real,parameter :: ProtonmassInGeV=cProtonmass*clightspeed**2/(10**9*cElectroncharge)!proton mass in unit of GeV/c^2
    real,parameter :: ParticleRange1=2.0, ParticleRange2=4.0                            !parameter of the place of the initial particles
    integer,parameter :: nQ = 130,  nP = 20,  nR = 20                                   !number of grid of s_L, \mu, ln(p^3/3) axis      
    integer,parameter :: nVar = 14                                                      !number of vatiales in each line of the input files               
    integer,parameter :: x_ = 2, y_ = 3, z_ = 4, rho_=5 ,&
                        Bx_ = 10, By_= 11, Bz_ = 12, Wave1_ = 13, Wave2_ = 14           !Position of the variables in the input line
    integer,parameter :: nStride = 25                                                   !stride of line of the data reading
    integer,parameter :: iLagrIDSkipped = 20                                            !Lagr index that is gaped for field line at each time
    real :: DeltaMu = 2.0/(nP)                                                          !\delta\mu for each grid along the \mu axis
    real :: ParticleEnergyMin, ParticleEnergyMax                                        !Min and Max value for particle energy
    real :: LnP3min, LnP3max                                                            !Ln(P^3/3)_min, Ln(P^3/3)_max
    real :: DeltaP3_I(nR), DeltaLnP3                                                    ! \delta(P^3/3)  \deltaLn(P^3/3)
    real :: Volume_G(0:nQ+1, 0:nP+1, 0:nR+1)                                            !control volume 
    real :: VDF_G(-1:nQ+2, -1:nP+2, -1:nR+2)                                            !vdf
    real :: VDF_Final(1:nQ, 1:nR)                                                       !vdf for output
    real :: Hamiltonian2_N(-1:nQ+1,-1:nP+1,0:nR+1)                                      !the Poisson bracket with regard to the first and second vars 
    real :: Source_C(nQ,nP,nR)                                                          !Source of each time step
    real :: RawData1_II(nQ,nVar)                                                        !the raw data that is imported from MHD
    real :: GapData                                                                     !Store the gaped data 
    real :: InvB_C(1:nQ)                                                                !1/(2B)
    real :: DeltaSOverB_C(1:nQ)                                                         !Cell centered Delta S over B
    real :: DeltaS0_I(1:nQ-1)                                                           !Distance between the cell centers
    real :: LambdaMuMu_II(1:nQ,1:nR)                                                    !The data that is calculated directly from RawData
    real :: tOutput                                                                     !Time of each output file
    real :: HelioDist_I(nQ)                                                             !Heliocentric distance in the unit of rSun, for visualization
    real :: Energy_I(nR)                                                                !Energy in unit of KeV, for output
    real :: VDFValue                                                                    ! VDFSpectrum
    integer :: m, iQ, iP, iR, iStep=0, iStride, iHeaderInform, iGap, iLineBegin         !Loop integers
    integer :: nOutputFile=2                                                            !The number of output files
    integer :: IndexofFile                                                              !Number Index of File
    character(LEN=3) :: NameFileSuffix                                                  !The number index of output files in character
    character(LEN=1) :: RemainData                                                      !Read the remaining data afger the read of RawData
    logical :: DoExit
    real :: CoeffMuToxx = 27.0/14.0
    !momentum at the cell center                                                                                                                                                                                                                                                                                                                                        
    real :: Momentum                !Cell centered particle momentum                                                                                                                                                                                                                                                                                                    
    real :: V_I(1:nR)               !Cell centered particle velocity
    !\
    ! Open file that contain data
    !/
    open(unit=13,file='FieldLineOld.in')
    rewind(13)
    !\
    ! The data that is gaped
    !/
    GapData = 0.0
    do while(nint(GapData) /= iLagrIDSkipped); read(13,*)GapData; end do
    do iQ = 1, nQ
       read(13,*)(RawData1_II(iQ,m), m = 1,nVar)
       do iStride = 1, nStride -1
          read(13,*)
       end do
    end do
    close(13)
    !\
    ! Normalize process: This is because the unit we use in this code is different from the unit of input data 
    !/
    RawData1_II(:,x_ :z_ ) = RawData1_II(:,x_ :z_ )*rSun/clightspeed
    !\
    ! Initialize the Ln(p^3/3) axis
    !/
    ParticleEnergyMin = 1.0e-3       !Min energy in unit of Gev
    ParticleEnergyMax = 1.0          !Max energy in unit of Gev
    LnP3min = log(sqrt((ParticleEnergyMin+ProtonmassInGeV)**2-ProtonmassInGeV**2)**3.0/3.0) !give the min value of Ln(P^3/3). We set the minimum of energy to be 1MeV, then Ln(p^3/3) can be calculated
    LnP3max = log(sqrt((ParticleEnergyMax+ProtonmassInGeV)**2-ProtonmassInGeV**2)**3.0/3.0) !give the max value of Ln(P^3/3). We set the maximum of energy to be 1GeV, then Ln(p^3/3) can be calculated 
    DeltaLnP3 = (LnP3max-LnP3min)/nR !give the \DeltaLn(P^3/3)
    do iR=1,nR
       DeltaP3_I(iR) = exp(LnP3min+iR*DeltaLnP3) - exp(LnP3min+(iR-1)*DeltaLnP3)!give \deltaP^3/3 for each grid
    end do
    do iR = 1, nR
       !\
       ! Calculate momentum at cell centers, in unit of GeV/c
       !/
       Momentum = (3.0*exp(LnP3min + (iR - 0.5)*DeltaLnP3))**(1.0/3.0)
       ! v/ (pc/1Gev)**(1/3)
       V_I(iR) = 1.0/sqrt(1.0 + ProtonmassInGeV**2/Momentum**2)
    end do

    !\                                       
    ! Initialize the distribution function VDF_G
    !/
    VDF_G = 1.0e-7
    do iR = 1, nR
       VDFValue = max(exp(-5.0/3.0*(LnP3min+iR*DeltaLnP3)), 1.0e-7)!Assume the momentum of particls to be power distribution: f ~ p^-5 
       do iQ = 1, nQ
          if (iQ<=ParticleRange1)then
             VDF_G(iQ,1:nP,iR) = VDFValue
          else
             VDF_G(iQ,1:nP,iR) = max(VDFValue * &!using a linear slope for the middle part
                  (ParticleRange2-real(iQ)) / &
                  (ParticleRange2-ParticleRange1), 1.0e-7)
          end if
       end do
    end do
    
    !\
    !boundary conditions
    !/
    !use reflective boundary condition for \mu axis
    VDF_G(:,    0, :) = VDF_G(:,    1, :)
    VDF_G(:,   -1, :) = VDF_G(:,    2, :)
    VDF_G(:, nP+1, :) = VDF_G(:,   nP, :)
    VDF_G(:, nP+2, :) = VDF_G(:, nP-1, :)
    !use outflow boundary condition for s_L axis
    VDF_G(-1,   :, :) = VDF_G(1,    :, :)
    VDF_G(0,    :, :) = VDF_G(1,    :, :)
    VDF_G(nQ+1, :, :) = VDF_G(nQ,   :, :)
    VDF_G(nQ+2, :, :) = VDF_G(nQ,   :, :)
    !use outflow boundary condition for Ln(p^3/3) axis
    VDF_G(:, :,   -1) = VDF_G(:, :,    1)
    VDF_G(:, :,    0) = VDF_G(:, :,    1)
    VDF_G(:, :, nR+1) = VDF_G(:, :,   nR)
    VDF_G(:, :, nR+2) = VDF_G(:, :,   nR)
    
    !\
    ! Initialize output data
    !/
    Source_C = 0.0
    
    !\
    ! Calculate the position for the s_L axis
    !/
    do iQ = 1, nQ
       HelioDist_I(iQ) = sqrt(RawData1_II(iQ,x_)**2 + RawData1_II(iQ,y_)**2 + RawData1_II(iQ,z_)**2) * clightspeed / rSun
    end do
    
    !\
    ! Calculate the position for ln(p^3/3) axis: note that the energy is in the unit of KeV now
    !/
    do iR=1,nR
       Energy_I(iR) = 6.0 + log10(sqrt((3.0*exp(LnP3min + (iR - 0.5)*DeltaLnP3))**(2.0/3.0) + ProtonmassInGeV**2) - ProtonmassInGeV)
    end do
    
    !\
    ! Calculate the VDF for output
    !/
    VDF_Final = 0.0
    do iP = 1, nP
       VDF_Final = VDF_Final + VDF_G(1:nQ, iP, 1:nR)*DeltaMu
    end do
    
    !\
    ! Give the initial outpur of VDF
    !/
    call save_plot_file(NameFile='test_multipoissonscatter_000.out', &
         TypeFileIn='ascii', TimeIn=Time, nStepIn = iStep,    &
         NameVarIn='sL[rSun] log10_E[KeV] VDF'  ,             &
         CoordMinIn_D=(/HelioDist_I(1) ,      Energy_I(1)/),  &
         CoordMaxIn_D=(/HelioDist_I(nQ),      Energy_I(nR)/), &
         StringFormatIn = '(3F15.7)',                         &  
         Coord1In_I = HelioDist_I,                            &
         Coord2In_I = Energy_I,                               &
         VarIn_II = log10(VDF_Final(1:nQ,1:nR)))
    
    !\
    ! ***Start the main simulation loop***
    ! In the loop, we output several snapshot of the distribution function VDF 
    !/
    Hamiltonian2_N = 0.0
    call calc_initial_data
    do iQ=1,nQ
       do iR=1,nR
          Volume_G(iQ,0:nP+1,iR)=DeltaSOverB_C(iQ)*DeltaMu*DeltaP3_I(iR)
       end do
    end do
    !boundary conditions for total volume
    Volume_G(0,   :,   :) = Volume_G(1 ,:, :)
    Volume_G(nQ+1,:,   :) = Volume_G(nQ,:, :)
    Volume_G(:,   :,   0) = Volume_G(: ,:, 1)
    Volume_G(:,   :,nR+1) = Volume_G(: ,:,nR)

    call explicit(nQ, nP, nR, VDF_G, Volume_G, Source_C, &
                  Hamiltonian2_N, DtOut = DtFixed, CFLIn=0.99)
    
    do IndexofFile=1,nOutputFile
       tOutput = TimeOut/real(nOutputFile)*real(IndexofFile)
       !\
       ! In this loop, we first update the RawData and calculate 3 Hamiltonian functions and total volume 
       ! Then we use subroutine explicit3 to calculate Source_C and update VDF
       !/
       do
          
          iStep = iStep + 1
          if(Time + DtFixed >= tOutput)then
             Dt = tOutput - Time
             DoExit = .true.
          else
             DoExit = .false.
             Dt = DtFixed
          end if
          Time = Time + Dt
          if (Diffuornot == 0) then
             call explicit(nQ, nP, nR, VDF_G, Volume_G, Source_C,&
                  Hamiltonian2_N, DtIn = Dt)!Note that we use a fixed time step here
             VDF_G(1:nQ, 1:nP, 1:nR) = VDF_G(1:nQ, 1:nP, 1:nR) + Source_C/Volume_G(1:nQ,1:nP,1:nR)
             Source_C=0.0
             call calc_scatter(Dt)
          else
             call calc_scatterDiffu(Dt)
          end if   
          !use reflective boundary condition for \mu axis
          VDF_G(:,    0, :) = VDF_G(:,    1, :)
          VDF_G(:,   -1, :) = VDF_G(:,    2, :)
          VDF_G(:, nP+1, :) = VDF_G(:,   nP, :)
          VDF_G(:, nP+2, :) = VDF_G(:, nP-1, :)
          !use outflow boundary condition for s_L axis
          VDF_G(-1,   :, :) = VDF_G(1,    :, :)
          VDF_G(0,    :, :) = VDF_G(1,    :, :)
          VDF_G(nQ+1, :, :) = VDF_G(nQ,   :, :)
          VDF_G(nQ+2, :, :) = VDF_G(nQ,   :, :)
          !use outflow boundary condition for ln(p^3/3) axis
          VDF_G(:, :,   -1) = VDF_G(:, :,    1)
          VDF_G(:, :,    0) = VDF_G(:, :,    1)
          VDF_G(:, :, nR+1) = VDF_G(:, :,   nR)
          VDF_G(:, :, nR+2) = VDF_G(:, :,   nR)
          if(DoExit) EXIT
       end do
       
       !\
       ! To avoid some mistakes when taking log10 to VDF_G, we here set a min value for VDF_G
       !/
       VDF_G = max(VDF_G, 1.0e-7)
       
       !\
       ! Calculate the final position
       !/
       do iQ = 1, nQ
          HelioDist_I(iQ) = sqrt(RawData1_II(iQ,x_)**2 + RawData1_II(iQ,y_)**2 + RawData1_II(iQ,z_)**2) * clightspeed / rSun
       end do
       
       !\                                                                                
       ! Calculate the VDF for output                                                                                                            
       !/
       VDF_Final = 0.0
       do iP = 1, nP
          VDF_Final = VDF_Final + VDF_G(1:nQ, iP, 1:nR)*DeltaMu
       end do
       
       !\
       ! Set the name of files, here the largest number of files is 999
       !/
       write(NameFileSuffix,'(I3.3)')IndexOfFile
       
       !\
       ! Save the files
       !/
       call save_plot_file(NameFile='test_multipoissonscatter_'//NameFileSuffix//'.out', &
         TypeFileIn='ascii', TimeIn=Time, nStepIn = iStep,       &
         NameVarIn='sL[rSun] log10_E[KeV] VDF'  ,                &
         CoordMinIn_D=(/HelioDist_I(1) ,      Energy_I(1)/),     &
         CoordMaxIn_D=(/HelioDist_I(nQ),      Energy_I(nR)/),    &
         StringFormatIn = '(3F15.7)',                            &
         Coord1In_I = HelioDist_I,                               &
         Coord2In_I = Energy_I,                                  &
         VarIn_II = log10(VDF_Final(1:nQ,1:nR)))
       
       
    end do
    
  contains
    
    
    subroutine calc_initial_data
      real :: MagneticOverWave_C(1:nQ), LarmorRadius13_C(1:nQ), Lmax_C(1:nQ)
      real :: MidPoint_ID(nQ-1,x_:z_), DeltaS_I(nQ)!midpoint for to consecutive points, \deltas
      real :: Momentum, B
      real :: CoeffMu
      integer :: iQ, iR!loop integer
      real :: InvB_F(0:nQ)!intermediate matrix or arrays  
      
      
      !\
      ! Calculate values at Old time first
      !/
      !calculate midpoints
      do iQ=1,nQ-1
         MidPoint_ID(iQ,:) = (RawData1_II(iQ+1,x_:z_) + RawData1_II(iQ,x_:z_))*0.5
      end do
      !calculate deltas
      do iQ=2,nQ-1
         DeltaS_I(iQ) = sqrt(sum((MidPoint_ID(iQ,x_:z_) - MidPoint_ID(iQ-1,x_:z_))**2))
      end do
      !Linear interpolate the delats such that there will be nQ deltas
      DeltaS_I(1)  = 2*sqrt(sum((MidPoint_ID(1,   x_:z_) - RawData1_II(1, x_:z_))**2))
      DeltaS_I(nQ) = 2*sqrt(sum((MidPoint_ID(nQ-1,x_:z_) - RawData1_II(nQ,x_:z_))**2))
      !\
      ! calculate part of hamiltonian function (1/2B) and \DeltaS/B at grid center
      !/
      do iQ=1,nQ
         B = sqrt(sum(RawData1_II(iQ, Bx_:Bz_)**2))
         InvB_C(iQ) = 0.50/B
         DeltaSOverB_C(iQ) = DeltaS_I(iQ)*2.0*InvB_C(iQ)
         MagneticOverWave_C(iQ) = (cMu*sum(RawData1_II(iQ, wave1_:wave2_)) + B**2)/&
                                  (cMu*sum(RawData1_II(iQ, wave1_:wave2_)))
         LarmorRadius13_C(iQ) = (1.0e9/(cLightSpeed**2*B))**(1.0/3.0)
         Lmax_C(iQ) = 0.03*sqrt(sum(RawData1_II(iQ,x_:z_)**2))
      end do
      do iQ=1,nQ-1
         DeltaS0_I(iQ) = sqrt(sum((RawData1_II(iQ+1,x_:z_)-RawData1_II(iQ,x_:z_))**2))
      end do
      
      !\
      ! Calculate \lambda_mumu and \lambda_xx
      !/
      CoeffMu = 6.0/(2.0**(2.0/3.0)*cPi**(5.0/3.0))
      do iR = 1, nR
         Momentum = (3.0*exp(LnP3min + (iR - 0.5)*DeltaLnP3))**(1.0/3.0)
         do iQ = 1, nQ
            LambdaMuMu_II(iQ, iR) = CoeffMu*MagneticOverWave_C(iQ)*LarmorRadius13_C(iQ)*Lmax_C(iQ)**(2.0/3.0)*Momentum**(1.0/3.0)
         end do
      end do
      
      !\
      ! calculate the value of part of the hamiltonian function (1/2B) on the boundary of grid
      !/                     
      do iQ=1,nQ-1
         InvB_F(iQ) = (InvB_C(iQ+1) + InvB_C(iQ))*0.5
      end do
      InvB_F(0)  = InvB_C(1 ) - (InvB_C(2 ) - InvB_C(1   ))*0.5
      InvB_F(nQ) = InvB_C(nQ) + (InvB_C(nQ) - InvB_C(nQ-1))*0.5
      
      do iQ=0,nQ
         do iP=0,nP
            do iR=1,nR
               !Consider the law of relativity, v=1/sqrt(1+m^2*c^2/p^2), we can calculate v as a function of p
               !Note that light speed is the unit of speed in our code, so we do not need to multiply a c^2 in the following expression
               Hamiltonian2_N(iQ,iP,iR) = InvB_F(iQ)*(1.0-(-1.0+real(iP)*DeltaMu)**2)*V_I(iR)
            end do
         end do
      end do
      
      !\
      ! calculate the Hamiltonian function that is to be used in actual calculation, which is \deltaH\tuta!!
      !/
      do iR=1,nR
         Hamiltonian2_N(0:nQ,0:nP,iR) = Hamiltonian2_N(0:nQ,0:nP,iR)*DeltaP3_I(iR)
      end do
      
      !\
      ! Boundary condition of Hamiltonian function
      !/
      Hamiltonian2_N(-1  ,:   ,:   ) = Hamiltonian2_N(0 ,:   ,: )
      Hamiltonian2_N(nQ+1,:   ,:   ) = Hamiltonian2_N(nQ,:   ,: )
      Hamiltonian2_N(:   ,-1  ,:   ) = Hamiltonian2_N(: ,1   ,: )
      Hamiltonian2_N(:   ,nP+1,:   ) = Hamiltonian2_N(: ,nP-1,: )
      Hamiltonian2_N(:   ,:   ,0   ) = Hamiltonian2_N(: ,:   ,1 )
      Hamiltonian2_N(:   ,:   ,nR+1) = Hamiltonian2_N(: ,:   ,nR)
      
      !\
      ! because of the requirement of equation, the Hamiltonian with a minus poisson bracket must be its opposite
      !/
      Hamiltonian2_N = - Hamiltonian2_N       
    end subroutine calc_initial_data

    
    !=================================
    !calculate scatter: \deltaf/\deltat = (Dmumu*f_mu)_mu
    subroutine calc_scatter(Dt)
      use ModDiffusion, ONLY: advance_diffusion1,tridiag
      real, intent(in) :: Dt
      
      !Dmumu for each fixed iQ and iR. Dmumu is the coefficient of diffusion along the pitch angle, \mu
      real :: DMuMu_I(0:nP)
      
      !Factorize DMuMu, each factor being only a function of s_L, \mu, ln(p^3/3) respectively    
      real :: FactorMU_F(0:nP)
      
      !lower, main, upper diagonal, output value of the calculation of scatter 
      real :: L_I(1:nP), M_I(1:nP), U_I(1:nP), W_I(1:nP)
      
      !Alfven wave speed
      real :: AlfvenSpeed

      !Face centered value of the pitch angle
      real :: Mu_I(0:nP)
      
      real :: DtOverDMu2              !\Deltat/\Delta\mu^2
      real :: LowerLimit              !Lower limit of Factor_mu
      integer :: iQ, iP, iR           !loop integers
      integer :: switch1, Pnumber     !control parameter
      !----------------------------------      
      
      !\
      ! Calculate factorized diffusion coefficient
      !/
      LowerLimit = (1.0 - DeltaMu**2)*abs(DeltaMu)**(2.0/3.0)
      
      do iP = 0, nP
         Mu_I(iP) = -1.0 + (real(iP))*DeltaMu
         FactorMU_F(iP) = (1.0 - Mu_I(iP)**2)*abs(Mu_I(iP))**(2.0/3.0)
      end do
      DtOverDMu2 = Dt/DeltaMu**2
      !\
      ! In this loop, for each fixed iQ and iR, we calculate the effect of scatter along \mu axis
      !/
      do iQ = 1, nQ
         AlfvenSpeed = 0.5/(InvB_C(iQ)*clightspeed*sqrt(RawData1_II(iQ,rho_)*cAtomicMass*cMu))
         do iR = 1, nR
            DMuMu_I = 0.0
            switch1 = 1
            
            do iP = 0, nP
               !use the selection to control whether we will floor the value of D_\mu\mu or not
               if (V_I(iR)*abs(Mu_I(iP))/AlfvenSpeed>=10.0) then
                  DMuMu_I(iP) = V_I(iR)/LambdaMuMu_II(iQ,iR)*FactorMU_F(iP)*DtOverDMu2
               else
                  if (switch1==1) Pnumber = iP
                  switch1 = 0
                  !set the lower limit of D_\mu\mu when |\mu| is close to zero
                  DMuMu_I(iP) = V_I(iR)/LambdaMuMu_II(iQ,iR)*max(FactorMU_F(Pnumber), LowerLimit)*DtOverDMu2
               end if
            end do
            L_I = - DMuMu_I(0:nP - 1)
            U_I = - DMuMu_I(1:nP)
            M_I = 1 - L_I - U_I
            W_I = VDF_G(iQ, 1:nP, iR)
            call tridiag(nP, L_I, M_I, U_I, W_I, VDF_G(iQ, 1:nP, iR))
         end do
      end do
    end subroutine calc_scatter    
    !=================================
    !calculate scatter: \deltaf/\deltat = (Dmumu*f_mu)_mu
    subroutine calc_scatterDiffu(Dt)
      use ModDiffusion, ONLY: advance_diffusion1
      real, intent(in) :: Dt
      
      !Inner and Outer coffecient of diffusion equation
      real :: DInner_I(nQ), DOuter_I(nQ), Dist_I(nQ)
      real :: dxx(nQ)
      integer :: iQ, iP, iR           !loop integers
      !----------------------------------      
      
      !\
      ! Calculate factorized diffusion coefficient
      !/
      
      !\
      ! Calculate Dist
      !/
      Dist_I = 1.0
      Dist_I(1:nQ - 1) = DeltaS0_I
      !\
      ! In this loop, for each fixed iP and iR, we calculate the effect of scatter along s_L axis
      !/
      do iP = 1, nP   
         do iR = 1, nR
            DInner_I = 2.0/3.0*V_I(iR)*CoeffMuToxx*LambdaMuMu_II(:,iR)*InvB_C
            DOuter_I = 0.50/(InvB_C)!there should be no Dt here!!!
            call advance_diffusion1(Dt, nQ, Dist_I, VDF_G(1:nQ, iP, iR), DOuter_I, DInner_I)
         end do
      end do
      dxx = dinner_I/2.0/InvB_C
    end subroutine calc_scatterDiffu
  end subroutine test_multipoisson_bracket
end module ModTestPoissonBracketAndScatter
!=============================================================================
!=============================================================================
program test_poisson_bracket
  use ModTestPoissonBracket, test => test_poisson_bracket, test3 => test_multipoisson_bracket
  use ModNumConst,       ONLY: cTwoPi
  use ModTestPoissonBracketAndScatter, test_scatter => test_multipoisson_bracket
  implicit none
  call test(cTwoPi)
  call test3(50.0)
  !\
  ! Test for comparing two diffusions is long. It should be repeated twice with
  ! Diffuornot = 0 or 1
  !/
  !call test_scatter(50.0)
end program test_poisson_bracket
