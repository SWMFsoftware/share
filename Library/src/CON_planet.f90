!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!
!
module CON_planet

  ! Physical information about the planet. The planet is described
  ! with its name. Default values can be set with {\bf planet\_init}.
  ! Simplifying assumptions, such as no rotation, aligned magnetic
  ! and rotational axes etc. can be made.
  !
  ! This is a public class. The variables should be modified by CON only.
  ! Components can only access the data through the inquiry methods
  ! via the {\bf CON\_physics} class.

  use ModNumConst, ONLY: cTwoPi, cDegToRad, cPi
  use ModPlanetConst
  use ModTimeConvert, ONLY: TimeType, time_int_to_real
  use ModUtilities, ONLY: CON_stop

  ! revision history:
  ! 01Aug03 - Aaron Ridly <ridley@umich.edu> and
  !           Gabor Toth <gtoth@umich.edu>   - initial prototype/prolog/code
  ! 23Mar03 - added get_planet subroutine for OO type access
  ! 06May04 - K.C. Hansen and G. Toth added Saturn
  !           G.Toth fixed bugs in degree to radian conversions

  implicit none

  save

  character(len=*), parameter, private :: NameMod='CON_planet'

  character (len=lNamePlanet) :: NamePlanet = ''

  ! Define variables
  real           :: RadiusPlanet
  real           :: MassPlanet
  real           :: TiltRotation
  real           :: IonosphereHeight
  real           :: OmegaPlanet     ! Rotation + Orbit
  real           :: OmegaRotation    ! Rotation
  real           :: RotPeriodPlanet ! Rotation
  real           :: OmegaOrbit      ! Orbit
  real           :: AngleEquinox
  type(TimeType) :: TimeEquinox = TimeType(2000, 1, 1, 0, 0, 0, &
       0.0, 0.0_REAL8_, '20000101000000')
  !$acc declare create(OmegaPlanet, OmegaRotation, AngleEquinox, TimeEquinox)

  ! Magnetic field type and strength in teslas
  character (len=lTypeBField) :: TypeBField = 'DIPOLE'
  real                        :: DipoleStrength
  !$acc declare create(DipoleStrength)
  real    :: MagAxisThetaGeo  ! Permanent theta  in GEO
  real    :: MagAxisPhiGeo    ! Permanent phi    in GEO

  ! Orientation of the axes
  real    :: RotAxisTheta      ! Permanent theta angle in GSE
  real    :: RotAxisPhi        ! Permanent phi   angle in GSE
  real    :: MagAxisTheta      ! Current   theta angla in GSE
  real    :: MagAxisPhi        ! Current   phi   angla in GSE
  !$acc declare create(MagAxisTheta, MagAxisPhi)

  ! Offset of the magnetic field center
  real    :: MagCenter_D(3) = [ 0.0, 0.0, 0.0]
  !$acc declare create(MagCenter_D)

  ! Optional changes relative to the "real" planet
  logical :: UseRotation     = .true.
  logical :: UseAlignedAxes  = .false.
  logical :: UseRealRotAxis  = .true.
  logical :: UseSetRotAxis   = .false.
  logical :: UseRealMagAxis  = .true.
  logical :: UseSetMagAxis   = .false.
  !$acc declare create(UseRotation)

  ! Frequency of updating the magnetic field information
  logical :: DoUpdateB0      = .true.
  real    :: DtUpdateB0      = 0.0001
  !$acc declare create(DoUpdateB0, DtUpdateB0)

  ! A primary axis is set to the true value
  ! a secondary axis is aligned with the primary axis
  logical :: IsRotAxisPrimary = .true., IsMagAxisPrimary = .true.

  ! A logical to indicate it the default parameters have been
  ! modified for this planet
  logical :: IsPlanetModified = .false.

  ! Variables added for the multipole option for the magnetic field
  logical :: UseMultipoleB0 = .false.
  integer :: MaxHarmonicDegree = -1
  character(len=200) :: NamePlanetHarmonicsFile
  real, allocatable :: gPlanet_II(:, :), hPlanet_II(:, :)

  ! Orbit parameters
  ! Logical specifying if we use these elements
  ! rather thhan the harwired in geopack and more
  ! accurate description applicable to the Earth only
  logical :: UseOrbitElements = .false.
  real :: rOrbitPlanet     ! [m]
  real :: Excentricity    ! Dimless
  ! Euler angles of orbit (in degs):
  real :: RightAscension, Inclination, ArgPeriapsis
  ! Derived orbital elements
  ! 1. Matrix to transform from orbital coordinates
  !    xOrb, yOrb to Xyz in HGI:
  !    XyzHgi = matmul(HgiOrb_DD, [xOrb, yOrb, 0.0])
  real :: HgiOrb_DD(3,3)
  ! 2. Major and minor semiaxes
  real :: SemiMinorAxis, SemiMajorAxis
  !$acc declare create(rOrbitPlanet, Excentricity)
  !$acc declare create(RightAscension, Inclination, ArgPeriapsis)
  !$acc declare create(HgiOrb_DD,SemiMinorAxis, SemiMajorAxis)
contains
  !============================================================================

  subroutine set_planet_defaults

    ! Initialize parameters for Earth as the default planet.  This is in case
    ! there is no \#PLANET command.

    character(len=*), parameter:: NameSub = 'set_planet_defaults'
    !--------------------------------------------------------------------------
    Planet_          = Earth_
    NamePlanet       = NamePlanet_I(Earth_)
    RadiusPlanet     = rPlanet_I(Earth_)
    MassPlanet       = mPlanet_I(Earth_)
    RotPeriodPlanet  = RotationPeriodPlanet_I(Earth_)
    TiltRotation     = TiltPlanet_I(Earth_)
    IonosphereHeight = IonoHeightPlanet_I(Earth_)
    OmegaOrbit       = cTwoPi/OrbitalPeriodPlanet_I(Earth_)
    OmegaPlanet      = OmegaOrbit + cTwoPi/RotationPeriodPlanet_I(Earth_)
    AngleEquinox     = cTwoPi * &
         ( iHourEquinoxPlanet_I(Earth_) * 3600 &
         + iMinuteEquinoxPlanet_I(Earth_) * 60 &
         + iSecondEquinoxPlanet_I(Earth_) &
         + FracSecondEquinoxPlanet_I(Earth_) &
         ) / (24 * 3600)
    TimeEquinox      = TimeType(&
         iYearEquinoxPlanet_I(Earth_), &
         iMonthEquinoxPlanet_I(Earth_), &
         iDayEquinoxPlanet_I(Earth_), &
         iHourEquinoxPlanet_I(Earth_),   &
         iMinuteEquinoxPlanet_I(Earth_), &
         iSecondEquinoxPlanet_I(Earth_), &
         FracSecondEquinoxPlanet_I(Earth_), &
         0.0_Real8_, '20000320073500')
    call time_int_to_real(TimeEquinox)
    TypeBField       = TypeBFieldPlanet_I(Earth_)
    DipoleStrength   = DipoleStrengthPlanet_I(Earth_)
    MagAxisThetaGeo  = bAxisThetaPlanet_I(Earth_)
    MagAxisPhiGeo    = bAxisPhiPlanet_I(Earth_)
    rOrbitPlanet     = rOrbitPlanet_I(Earth_)  ! [m]
    Excentricity     = Excentricity_I(Earth_) ! dimless
    ! Euler angles of orbit (read in degs, convert to rads):
    RightAscension   = RightAscension_I(Earth_)
    RightAscension   = RightAscension*cDegToRad
    Inclination      = Inclination_I(Earth_)
    Inclination      = Inclination*cDegToRad
    ArgPeriapsis     = ArgPeriapsis_I(Earth_)
    ArgPeriapsis     = ArgPeriapsis*cDegToRad
    call get_orbit_elements
    !$acc update device(OmegaPlanet, AngleEquinox, OmegaRotation,             &
    !$acc TimeEquinox, MagAxisPhi, MagAxisTheta,  rOrbitPlanet, Excentricity, &
    !$acc RightAscension, Inclination, ArgPeriapsis)

  end subroutine set_planet_defaults
  !============================================================================
  function is_planet_init(NamePlanetIn) result(IsKnown)

    character(len=*), intent(in) :: NamePlanetIn

    ! return value
    logical :: IsKnown

    ! Initialize parameters for the planet identified by its name and
    ! return true if the planet is known. If it is not known return false.
    ! Store the name in either case. The planet data can be initialized at most
    ! once.
    integer :: i

    logical :: IsInitialized = .false.

    character(len=*), parameter:: NameSub = 'is_planet_init'
    !--------------------------------------------------------------------------
    IsKnown = .true.
    if(IsInitialized)then
       if(NamePlanet == NamePlanetIn) RETURN
       call CON_stop(NameSub//&
            ' ERROR: attempt to change planet name from '// &
            trim(NamePlanet)//' to '//NamePlanetIn)
    end if

    NamePlanet    = NamePlanetIn
    IsInitialized = .true.

    IsKnown       = .false.
    do i = NoPlanet_, MaxPlanet
      if (NamePlanet == NamePlanet_I(i)) then
         IsKnown = .true.
         Planet_ = i
         EXIT
      end if
    end do

    if (.not. IsKnown)  then
       Planet_    = NewPlanet_
       NamePlanet = NamePlanetIN
    end if

    ! Set all values for the selected planet
    RadiusPlanet     = rPlanet_I(Planet_)
    MassPlanet       = mPlanet_I(Planet_)
    RotPeriodPlanet  = RotationPeriodPlanet_I(Planet_)
    TiltRotation     = TiltPlanet_I(Planet_)
    IonosphereHeight = IonoHeightPlanet_I(Planet_)
    if (RotationPeriodPlanet_I(Planet_) == 0.0) then
       OmegaRotation = 0.0
    else
       OmegaRotation = cTwoPi/RotationPeriodPlanet_I(Planet_)
    end if
    if (OrbitalPeriodPlanet_I(Planet_) == 0.0) then
       OmegaOrbit = 0.0
    else
       OmegaOrbit = cTwoPi/OrbitalPeriodPlanet_I(Planet_)
    end if
    OmegaPlanet     = OmegaRotation + OmegaOrbit
    AngleEquinox =  &
         cTwoPi*(iHourEquinoxPlanet_I(Planet_)*3600 &
         +       iMinuteEquinoxPlanet_I(Planet_)*60 &
         +       iSecondEquinoxPlanet_I(Planet_)  &
         +       FracSecondEquinoxPlanet_I(Planet_))/(24*3600)
    TimeEquinox  = TimeType(&
         iYearEquinoxPlanet_I(Planet_), &
         iMonthEquinoxPlanet_I(Planet_), &
         iDayEquinoxPlanet_I(Planet_), &
         iHourEquinoxPlanet_I(Planet_), &
         iMinuteEquinoxPlanet_I(Planet_), &
         iSecondEquinoxPlanet_I(Planet_), &
         FracSecondEquinoxPlanet_I(Planet_), &
         0.0_Real8_, '')

    ! Magnetic field type and strength in teslas
    TypeBField        = TypeBFieldPlanet_I(Planet_)
    DipoleStrength    = DipoleStrengthPlanet_I(Planet_)
    MagAxisThetaGeo   = bAxisThetaPlanet_I(Planet_)  ! Permanent theta  in GEO
    MagAxisPhiGeo     = bAxisPhiPlanet_I(Planet_)    ! Permanent phi    in GEO

    rOrbitPlanet     = rOrbitPlanet_I(Planet_)  ! [m]
    Excentricity     = Excentricity_I(Planet_) ! dimless
    ! Euler angles of orbit (read in degs, convert to rads):
    RightAscension   = RightAscension_I(Planet_)*cDegToRad
    Inclination      = Inclination_I(Planet_)*cDegToRad
    ArgPeriapsis     = ArgPeriapsis_I(Planet_)*cDegToRad
    call get_orbit_elements
    ! For Enceladus the dipole is at Saturn's center
    if(Planet_==Enceladus_) MagCenter_D(2)    = 944.23

    !$acc update device(OmegaPlanet, AngleEquinox, OmegaRotation, TimeEquinox,&
    !$acc MagAxisPhi, MagAxisTheta,  rOrbitPlanet, Excentricity,              &
    !$acc RightAscension, Inclination, ArgPeriapsis)

  end function is_planet_init
  !============================================================================
  subroutine read_planet_var(NameCommand)

    use ModUtilities, ONLY: upper_case
    use ModReadParam, ONLY: read_var, lStringLine
    use ModIoUnit, ONLY: UnitTmp_

    character (len=*), intent(in) :: NameCommand

    ! Planet related temporary variables
    character (len=lNamePlanet) :: NamePlanetIn
    character (len=lStringLine) :: NamePlanetCommands=''
    logical :: UseNonDipole

    integer :: m, n, m1, n1
    character(len=100) :: StringHeader

    character(len=*), parameter:: NameSub = 'read_planet_var'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#PLANET", "#MOON", "#COMET")

       call read_var('NamePlanet',NamePlanetIn)
       call upper_case(NamePlanetIn)

       ! Check if planet has been already initialized
       if(NamePlanet == NamePlanetIn) RETURN

       if (NamePlanetCommands /= '') &
            call CON_stop(NameSub// &
            ' ERROR: #PLANET should precede '// &
            NamePlanetCommands)

       ! If planet has not been initialized ...
       if ( .not. is_planet_init(NamePlanetIn) ) then
          call read_var('RadiusPlanet', RadiusPlanet)
          call read_var('MassPlanet',   MassPlanet)
          call read_var('OmegaPlanet',  OmegaPlanet)
          if (OmegaPlanet /= 0.0) then
             RotPeriodPlanet = cTwoPi/OmegaPlanet
          else
             RotPeriodPlanet = 0.0
          end if
          call read_var('TiltRotation', TiltRotation)
          TiltRotation = TiltRotation * cDegToRad
          call read_var('TypeBField',   TypeBField)
          call upper_case(TypeBField)

          select case(TypeBField)
          case('NONE')
            MagAxisTheta   = 0.0
            MagAxisPhi     = 0.0
            UseSetMagAxis  = .true.
            UseRealMagAxis = .false.
            DipoleStrength = 0.0

          case('DIPOLE','QUADRUPOLE','OCTUPOLE')
            call read_var('MagAxisThetaGeo', MagAxisThetaGeo)
            MagAxisThetaGeo = MagAxisThetaGeo * cDegToRad
            call read_var('MagAxisPhiGeo',   MagAxisPhiGeo)
            MagAxisPhiGeo = MagAxisPhiGeo * cDegToRad
            call read_var('DipoleStrength',DipoleStrength)

            if (TypeBField == 'QUADRUPOLE') then
              call CON_stop(NameSub// &
              ' ERROR: quadrupole field unimplemented')
            endif

            if (TypeBField == 'OCTUPOLE') then
              call CON_stop(NameSub// &
                  ' ERROR: octupole field unimplemented')
            endif

            case default
              call CON_stop(NameSub// &
              ' ERROR: TypeBfield not specified for planet.'//TypeBField)

            end select
       end if
       !$acc update device(OmegaPlanet)
    case('#IDEALAXES')
       ! This is a short version of setting one axis parallel with Z
       ! and the other one aligned with it

       NamePlanetCommands = '#IDEALAXES ' // NamePlanetCommands
       IsPlanetModified = .true.

       UseRealRotAxis   = .false.
       IsRotAxisPrimary = .true.
       UseSetRotAxis    = .true.
       RotAxisTheta     = 0.0
       RotAxisPhi       = 0.0
       UseRealMagAxis   = .false.
       IsMagAxisPrimary = .false.
       UseSetMagAxis    = .false.

    case('#ROTATIONAXIS')

       NamePlanetCommands = '#ROTATIONAXIS ' // NamePlanetCommands
       IsPlanetModified = .true.
       UseRealRotAxis = .false.

       call read_var('IsRotAxisPrimary', IsRotAxisPrimary)
       if (IsRotAxisPrimary) then

          UseSetRotAxis = .true.

          call read_var('RotAxisTheta', RotAxisTheta)
          if(RotAxisTheta < 0.0)call CON_stop(NameSub// &
               ' ERROR: negative tilt should be entered as Phi=180.0')
          RotAxisTheta = cDegToRad * RotAxisTheta

          call read_var('RotAxisPhi', RotAxisPhi)
          RotAxisPhi = cDegToRad * RotAxisPhi
       else
          if(.not.IsMagAxisPrimary)call CON_stop(NameSub// &
               ' ERROR: either rotation or magnetic axis must be primary')
       end if

    case('#MAGNETICAXIS')

       NamePlanetCommands = '#MAGNETICAXIS ' // NamePlanetCommands
       IsPlanetModified = .true.
       UseRealMagAxis = .false.

       call read_var('IsMagAxisPrimary', IsMagAxisPrimary)
       if (IsMagAxisPrimary) then

          UseSetMagAxis = .true.

          call read_var('MagAxisTheta', MagAxisTheta)
          if(MagAxisTheta < 0.0)call CON_stop(NameSub// &
               ' ERROR: negative tilt should be entered as Phi=180.0')
          MagAxisTheta = cDegToRad * MagAxisTheta

          call read_var('MagAxisPhi', MagAxisPhi)
          MagAxisPhi = cDegToRad * MagAxisPhi
       else
          if(.not.IsRotAxisPrimary)call CON_stop(NameSub//&
               ' ERROR: either rotation or magnetic axis must be primary')
       end if

    case('#MAGNETICCENTER')

       NamePlanetCommands = '#MAGNETICCENTER ' // NamePlanetCommands
       IsPlanetModified = .true.
       call read_var('MagneticCenterX',MagCenter_D(1))
       call read_var('MagneticCenterY',MagCenter_D(2))
       call read_var('MagneticCenterZ',MagCenter_D(3))

    case('#ROTATION')

       NamePlanetCommands = '#ROTATION ' // NamePlanetCommands
       IsPlanetModified = .true.

       call read_var('UseRotation', UseRotation)
       if (.not.UseRotation) then
          OmegaPlanet     = 0.0
          RotPeriodPlanet = 0.0
       else
          call read_var('Rotation period [hours]',  RotPeriodPlanet)
          RotPeriodPlanet = RotPeriodPlanet * 3600
          OmegaPlanet = cTwoPi / RotPeriodPlanet
       endif
       !$acc update device(UseRotation, OmegaPlanet)
    case('#NONDIPOLE')

       NamePlanetCommands = '#NONDIPOLE ' // NamePlanetCommands
       IsPlanetModified = .true.

       call read_var('UseNonDipole',UseNonDipole)
       if (.not.UseNonDipole) then
          TypeBField = 'DIPOLE'
       else
          call CON_stop(NameSub// &
               ' ERROR: nondipole magnetic field unimplemented')
       endif

    case('#DIPOLE')

       NamePlanetCommands = '#DIPOLE ' // NamePlanetCommands
       IsPlanetModified = .true.

       call read_var('DipoleStrength',DipoleStrength)
       if(DipoleStrength ==0.0)then
          TypeBField="NONE"
       else
          TypeBField="DIPOLE"
       end if

    case('#UPDATEB0')

       call read_var('DtUpdateB0',DtUpdateB0)

    case('#MULTIPOLEB0')

       call read_var('UseMultipoleB0', UseMultipoleB0)
       if(UseMultipoleB0) then

          ! If multipole is used, the magnetic axis is aligned
          ! with the rotation axis. Rotation axis may be real or
          ! user-specified using #ROTATIONAXIS command.
          IsPlanetModified = .true.
          UseRealRotAxis   = .true.
          IsRotAxisPrimary = .true.
          UseRealMagAxis   = .false.
          IsMagAxisPrimary = .false.
          UseSetMagAxis    = .false.

          ! Set multipole properties
          TypeBField="MULTIPOLE"
          call read_var('MaxHarmonicDegree', MaxHarmonicDegree)
          call read_var('NamePlanetHarmonicsFile', NamePlanetHarmonicsFile)

          if(.not.allocated(gPlanet_II) .and. .not.allocated(hPlanet_II)) then
             allocate(gPlanet_II(0:MaxHarmonicDegree, 0:MaxHarmonicDegree))
             allocate(hPlanet_II(0:MaxHarmonicDegree, 0:MaxHarmonicDegree))
          end if

          ! Read in the coefficients from the file
          open(UnitTmp_, file=NamePlanetHarmonicsFile, status='old')
          read(UnitTmp_, *) StringHeader
          do n = 0, MaxHarmonicDegree
             do m = 0, n
                ! Read data from the file, default fortran format
                read(UnitTmp_, *) &
                     n1, m1, gPlanet_II(n, m), hPlanet_II(n, m)
             end do
          end do
          close(UnitTmp_)

          call normalize_schmidt_coefficients
       end if
    case('#ORBIT')
       call read_var('OrbitalPeriodPlanet', OrbitalPeriodPlanet_I(Planet_))
       OmegaOrbit = cTwoPi/OrbitalPeriodPlanet_I(Planet_)
       ! Correct omega planet?
       call read_var('rOrbitPlanet',rOrbitPlanet) ! [m]
       call read_var('Excentricity',Excentricity) ! dimless
       ! read Euler angles of orbit (read in degs, convert to rads):
       call read_var('RightAscension', RightAscension)
       RightAscension = RightAscension*cDegToRad
       call read_var('Inclination', Inclination)
       Inclination = Inclination*cDegToRad
       call read_var('ArgPeriapsis',ArgPeriapsis)
       ArgPeriapsis = ArgPeriapsis*cDegToRad
       !$acc update device(OmegaOrbit,  rOrbitPlanet, Excentricity,    &
       !$acc RightAscension, Inclination, ArgPeriapsis)
       call get_orbit_elements
    case('#TIMEEQUINOX')
       call read_var('iYear',  TimeEquinox%iYear)
       call read_var('iMonth', TimeEquinox%iMonth)
       call read_var('iDay',   TimeEquinox%iDay)
       call read_var('iHour',  TimeEquinox%iHour)
       call read_var('iMinute',TimeEquinox%iMinute)
       call time_int_to_real(TimeEquinox)
    end select

  end subroutine read_planet_var
  !============================================================================
  subroutine check_planet_var(IsProc0, DoTimeAccurate)

    logical, intent(in) :: IsProc0, DoTimeAccurate

    ! The rotation and magnetic axes are aligned if any of them is not a
    ! primary axis or if multipoleB0 is used (implicitly).
    character(len=*), parameter:: NameSub = 'check_planet_var'
    !--------------------------------------------------------------------------
    UseAlignedAxes = (.not. IsRotAxisPrimary) .or. (.not. IsMagAxisPrimary)

    ! Warn if setting is unphysical
    if(UseSetMagAxis .and. UseRealRotAxis .and. IsProc0) &
         write(*,*)NameSub,' WARNING: magnetic axis is explicitly set ',&
         'while rotation axis is calculated from real time ?!'

    ! Check if there is a need to update the magnetic field
    DoUpdateB0 = DtUpdateB0 > 0.0 .and. DoTimeAccurate .and. UseRotation &
         .and. .not.(UseAlignedAxes .and. .not.UseMultipoleB0)

    ! Multipole B0 uses IdealAxes to set the proxy mag-rot axes,
    ! so it needs DoUpdateB0 to be True even when IdealAxes is used.
    ! If dipole is used with IdealAxes, no update is needed.

    !$acc update device(DoUpdateB0, DtUpdateB0)

  end subroutine check_planet_var
  !============================================================================
  subroutine get_planet( &
       NamePlanetOut, RadiusPlanetOut, MassPlanetOut, OmegaPlanetOut, &
       RotationPeriodOut, IonosphereHeightOut, &
       UseRotationOut, DipoleStrengthOut, DoUpdateB0Out, DtUpdateB0Out)

    character(len=*), optional, intent(out) :: NamePlanetOut
    real,             optional, intent(out) :: RadiusPlanetOut
    real,             optional, intent(out) :: MassPlanetOut
    real,             optional, intent(out) :: OmegaPlanetOut
    real,             optional, intent(out) :: RotationPeriodOut
    real,             optional, intent(out) :: IonosphereHeightOut
    logical,          optional, intent(out) :: UseRotationOut
    real,             optional, intent(out) :: DipoleStrengthOut
    logical,          optional, intent(out) :: DoUpdateB0Out
    real,             optional, intent(out) :: DtUpdateB0Out
    !--------------------------------------------------------------------------
    if(present(NamePlanetOut))      NamePlanetOut       = NamePlanet
    if(present(RadiusPlanetOut))    RadiusPlanetOut     = RadiusPlanet
    if(present(MassPlanetOut))      MassPlanetOut       = MassPlanet
    if(present(OmegaPlanetOut))     OmegaPlanetOut      = OmegaPlanet
    if(present(RotationPeriodOut))  RotationPeriodOut   = RotPeriodPlanet
    if(present(IonosphereHeightOut))IonosphereHeightOut = IonosphereHeight
    if(present(UseRotationOut))     UseRotationOut      = UseRotation
    if(present(DipoleStrengthOut))  DipoleStrengthOut   = DipoleStrength
    if(present(DoUpdateB0Out))      DoUpdateB0Out       = DoUpdateB0
    if(present(DtUpdateB0Out))      DtUpdateB0Out       = DtUpdateB0

  end subroutine get_planet
  !============================================================================
  subroutine normalize_schmidt_coefficients
    ! Instead of normalizing the Legendre polynomials, from Gauss to Schmidt
    ! -semi-normalized, we normalize the Schmidt coefficients instead for
    ! a more efficient calculation. We only need to do this once at the
    ! beginning of the run.
    integer :: m, n
    real, allocatable :: s_II(:,:)
    !--------------------------------------------------------------------------
    if(.not.allocated(s_II)) &
         allocate(s_II(0:MaxHarmonicDegree,0:MaxHarmonicDegree))

    s_II = 1.0

    do n = 1, MaxHarmonicDegree
      s_II(n,0) = s_II(n-1,0)*(2*n - 1.0)/n
      s_II(n,1) = s_II(n,0)*sqrt(2*n/(n + 1.0))
      gPlanet_II(n,0) = gPlanet_II(n,0)*s_II(n,0)
      hPlanet_II(n,0) = hPlanet_II(n,0)*s_II(n,0)
      gPlanet_II(n,1) = gPlanet_II(n,1)*s_II(n,1)
      hPlanet_II(n,1) = hPlanet_II(n,1)*s_II(n,1)
      do m = 2, n
        s_II(n,m) = s_II(n,m-1) * sqrt((n - m + 1.0)/(n + m))
        gPlanet_II(n,m) = gPlanet_II(n,m)*s_II(n,m)
        hPlanet_II(n,m) = hPlanet_II(n,m)*s_II(n,m)
      end do
    end do

    deallocate(s_II)

  end subroutine normalize_schmidt_coefficients
  !============================================================================
  subroutine get_orbit_elements

    use ModCoordTransform, ONLY: rot_matrix_x, rot_matrix_z

    ! The true formula for transformation from the 'orbital' coordinates
    ! (for the Earth, in the ecliptic plane) to an arbitrary reference frame)
    ! (any solar equatorial plane), the true formula is
    ! XyzStarEquator_D = matmul(HgiOrb_DD,XyzOrb_D)
    ! where
    ! HgiOrb_DD = rot_matrix_z(RightAscension).rot_matrix_x(Inclination).&
    !             rot_matrix_z(ArgPeriapsis),
    ! where RightAscension is the angle between x-axis of the solar
    ! equatorial coordinate system, and ascending node,
    ! Inclination is for the orbital plane with
    ! respect to the solar equatorial one and ArgPeriapsis is the angle
    ! between the direction to periapsis with that to ascending node.
    ! HOWEVER:
    ! 1. For HGI coordinate system thus defined RightAscension = cPi
    ! 2. and ArgPeriapsis and RightAscension are usually defined with regard
    !    to the equinox, so what we need here is their difference
    ! In this way the formula works for the Sun and the Earth or the star
    ! and exoplanet, (in which cases the definition of HGI is based on the
    ! particular choice of planet, otherwise the general formula shoud be
    ! used with properly re-defined RightAscension and ArgPeriapsis.

    !--------------------------------------------------------------------------
    HgiOrb_DD = matmul(rot_matrix_x(Inclination),&
         rot_matrix_z(ArgPeriapsis - RightAscension) )
    HgiOrb_DD = matmul(rot_matrix_z(cPi),HgiOrb_DD)
    SemiMajorAxis = rOrbitPlanet   ! Check acculacy!
    SemiMinorAxis= SemiMajorAxis*sqrt(1 - Excentricity**2)
    UseOrbitElements = NamePlanet /= 'EARTH'
    !$acc update device(HgiOrb_DD, SemiMajorAxis, SemiMinorAxis)

  end subroutine get_orbit_elements
  !============================================================================
  subroutine orbit_in_hgi(Time, XyzHgi_D, vHgi_D)

    use ModKind
    real(real8_), intent(in)  :: Time
    real, intent(out)           :: XyzHgi_D(3)  ! Coordinates in HGI
    real, optional, intent(out) :: vHgi_D(3)    ! Velocity in HGI
    ! True (in fact, not true) anomaly:
    real :: TrueAnomaly
    ! Coordinates and velocity in the orbital plane
    real              :: XyzOrbit_D(3),  vOrbit_D(3)
    !--------------------------------------------------------------------------
    ! Check accuracy!!!
    TrueAnomaly = OmegaOrbit*(Time - TimeEquinox%Time) - ArgPeriapsis
    TrueAnomaly= modulo(TrueAnomaly,cTwoPi)
    XyzOrbit_D = [SemiMajorAxis*(cos(TrueAnomaly) - Excentricity), &
         SemiMinorAxis*sin(TrueAnomaly), 0.0]
    XyzHgi_D = matmul(HgiOrb_DD, XyzOrbit_D)
    if(present(vHgi_D))then
       vOrbit_D = [-SemiMajorAxis*sin(TrueAnomaly), &
            SemiMinorAxis*cos(TrueAnomaly), 0.0]*OmegaOrbit
       vHgi_D =  matmul(HgiOrb_DD, vOrbit_D)
    end if
  end subroutine orbit_in_hgi
  !============================================================================
end module CON_planet
!==============================================================================
