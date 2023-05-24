!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!
!
module CON_star

  ! Physical information about the planet. The planet is described
  ! with its name. Default values can be set with {\bf planet\_init}.
  ! Simplifying assumptions, such as no rotation, aligned magnetic
  ! and rotational axes etc. can be made.
  !
  ! This is a public class. The variables should be modified by CON only.
  ! Components can only access the data through the inquiry methods
  ! via the {\bf CON\_physics} class.

  use ModConst, ONLY: cTwoPi, rSun, mSun, RotationPeriodSun, &
       cSecondPerDay
  use ModTimeConvert, ONLY: TimeType, time_int_to_real
  use ModUtilities, ONLY: CON_stop
  implicit none

  save

  character(len=*), parameter, private :: NameMod='CON_star'

  ! Define variables
  real           :: RadiusStar = rSun
  real           :: MassStar   = mSun
  real           :: OmegaStar  = cTwoPi/RotationPeriodSun
  !$acc declare create(OmegaStar)
  real           :: RotPeriodStar  = RotationPeriodSun
  ! Logical, claiming if the star is not the Sun 
  logical        :: UseStar
contains
  !============================================================================

  subroutine read_star_var(NameCommand)

    use ModReadParam, ONLY: read_var, lStringLine
    use ModIoUnit,    ONLY: UnitTmp_

    character (len=*), intent(in) :: NameCommand

    ! Planet related temporary variables
    character (len=3) :: NameStarIn

    character(len=*), parameter:: NameSub = 'read_star_var'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#STAR")
       UseStar=.true.
       call read_var('RadiusStar',         RadiusStar)  ! In rSun
       RadiusStar = RadiusStar*rSun
       call read_var('MassStar',           MassStar)
       MassStar = MassStar*mSun
       call read_var('RotPeriodStar',      RotPeriodStar)
       if(RotPeriodStar == 0.0)then
          OmegaStar = 0.0
       else
          RotPeriodStar = RotPeriodStar*cSecondPerDay
          OmegaStar = cTwoPi/RotPeriodStar
       end if
       !$acc update device(OmegaStar)
    case default
       call CON_stop('Unknwn NameCommand='//NameCommand//' in '//NameSub)
    end select
  end subroutine read_star_var
  !============================================================================

end module CON_star
!==============================================================================
