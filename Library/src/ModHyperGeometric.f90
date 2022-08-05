!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModHyperGeometric
  use ModMpi, ONLY: &
       iRealPrec ! 1, if the code is compiled with double precision
  use ModNumConst, ONLY: cPi
  use ModUtilities, ONLY: CON_stop
  implicit none
  PRIVATE ! Except
  public :: toroid_p ! tilde{P}_n function, related to k**3 * (k^\prime)**n
  public :: toroid_q ! tilde{Q}_n function, related to k**3 * (k^\prime)**n
  public :: scr_inductance    ! inductance of a superconducting ring
  public :: l0_ext_inductance ! the external inductance of a ring current
  public :: calc_elliptic_int_1kind  ! Elliptic integral, of the first kind
  public :: calc_elliptic_int_2kind  ! Elliptic integral, of the second kind
  real, parameter:: cTolerance_I(0:1) = [1.0e-7, 1.0e-15]
  real, parameter:: cEiler    = 0.57721566490
  real, parameter:: cSqrtPi   = 1.7724538509055159
contains
  !============================================================================
  real function psi_semi(n)
    integer, intent(in) :: n

    ! Definition of psi function: psi(z) = d\log(\Gamma(z))/dz
    ! For semiinteger z:
    ! psi(0.5 + n) = -C - log 4 +\sum_{k=1}^n{1/(k - 0.5)}
    integer :: k

    !--------------------------------------------------------------------------
    psi_semi = -cEiler -log(4.0)
    do k=1, abs(n)
       psi_semi = psi_semi +1.0/(real(k) - 0.50)
    end do
  end function psi_semi
  !============================================================================
  real function psi_int(n)
    integer, intent(in) :: n

    ! Definition of psi function: psi(z) = d\log(\Gamma(z))/dz
    ! For integer z:
    ! psi(1 + n) = -C +\sum_{k=1}^{n-1}{1.0/k}
    integer :: k

    !--------------------------------------------------------------------------
    psi_int = -cEiler
    do k = 1, n -1
       psi_int = psi_int +1.0/real(k)
    end do
  end function psi_int
  !============================================================================
  real function factorial(n)
    integer, intent(in) :: n
    ! n! = \Gamma(n+1)
    integer :: k

    !--------------------------------------------------------------------------
    factorial = 1.0
    do k=2, n
       factorial = factorial*real(k)
    end do
  end function factorial
  !============================================================================
  real function gamma_semi(n)
    integer, intent(in) :: n
    ! \Gamma(0.5 + n)
    integer :: k
    ! n >= -1!
    !--------------------------------------------------------------------------
    gamma_semi = -2.0*cSqrtPi
    do k = 0, n
       gamma_semi = gamma_semi * (real(k) - 0.50)
    end do
  end function gamma_semi
  !============================================================================

  ! Hypergeometric series
  real function hypergeom(A, B, C, Z)

    ! Input parameters and argument
    real, intent(in) :: A, B, C, Z

    ! Loop variable
    integer:: i

    ! Misc
    real :: aPlusI, bPlusI, cPlusI, rMember
    character(len=*), parameter:: NameSub = 'hypergeom'
    !--------------------------------------------------------------------------
    hypergeom = 1.0
    rMember   = 1.0
    aPlusI = a -1.0
    bPlusI = b -1.0
    cPlusI = c -1.0
    i = 0
    do while(abs(rMember) >= cTolerance_I(iRealPrec))
       i = i +1
       aPlusI = aPlusI + 1.0
       bPlusI = bPlusI + 1.0
       cPlusI = cPlusI + 1.0
       rMember = rMember*aPlusI*bPlusI*z/(cPlusI*i)
       hypergeom = hypergeom + rMember
    end do
  end function hypergeom
  !============================================================================
  recursive function hyper_semi_semi_int(nA, nB, nC, ZIn, OneMinusZIn) &
       RESULT(FuncVal)
    integer,        intent(in) :: nA, nB, nC
    real, optional, intent(in) :: ZIn, OneMinusZIn
    real    :: Z, OneMinusZ
    real    :: FuncVal
    ! Calculate hypergeometric series F(a, b, c, z), if
    ! semiinteger a = 0.5 + nA, nA = 0, 1, 2
    ! semiinteger b = 0.5 + nB, nB = 0, 1, 2
    ! integer     c = nC
    real :: A, B, C
    integer :: n ! Discriminator= c - a -b

    ! Loop variable
    integer :: i

    ! Misc
    real :: aPlusI, bPlusI, cPlusI, rMember, LogFactor

    character(len=*), parameter:: NameSub = 'hyper_semi_semi_int'
    !--------------------------------------------------------------------------
    if(present(ZIn))then
       Z = ZIn; OneMinusZ = 1.0 - Z
    elseif(present(OneMinusZIn))then
       OneMinusZ = OneMinusZIn; Z = 1.0 - OneMinusZ
    else
       call CON_stop(&
            NameSub//': ZIn or OneMinusZIn should be present')
    end if
    ! Real arguments of the hypergeometric function:
    A =  0.50 + real(nA); B = 0.50 + real(nB); C = real(nC)
    if (abs(z) < 0.50) then

       ! Direct summation of the hypergeometric series, well withing the
       ! convergence radius
       FuncVal = hypergeom(&
            A=A,        &
            B=B,        &
            C=C,        &
            Z=Z)
       RETURN
    end if

    ! Use the analytic extension to the singular point z=1
    ! OneMinusZ = 1.0 - z
    ! The difference C - (A+B) is integer. Calculate this.
    n = nC - (1 + nA + nB)

    ! The formulae for the "logarithmic case" (integer n)
    ! https://dlmf.nist.gov
    ! strongly depend on the sign of n. Consider case-by-case
    if(n==0)then
       ! Apply Eq. 15.8.10 in https://dlmf.nist.gov:
       LogFactor    = -log(OneMinusZ) + 2.0*psi_int(1) - &
            (psi_semi(nA) + psi_semi(nB))
       rMember      = factorial(nA + nB)/(gamma_semi(nA)*gamma_semi(nB))
       FuncVal = rMember*LogFactor
       aPlusI       = A - 1.0
       bPlusI       = B - 1.0
       i = 0
       do while(abs(rMember) >= cTolerance_I(iRealPrec))
          i = i + 1
          aPlusI = aPlusI + 1.0
          bPlusI = bPlusI + 1.0
          rMember = rMember*aPlusI*bPlusI/i**2*OneMinusZ
          LogFactor = LogFactor + 2.0/real(i) - 1.0/aPlusI - 1.0/bPlusI
          FuncVal = FuncVal + rMember*LogFactor
       end do
    elseif(n<0)then
       ! Apply Eq. 15.8.12 in https://dlmf.nist.gov:
       n = -n
       FuncVal = hyper_semi_semi_int(nA - n, nB - n, nC, ZIn, OneMinusZIn)/&
            OneMinusZ**n
    else
       ! Apply Eq. 15.8.10 in https://dlmf.nist.gov:
       rMember      = factorial(nC-1)*factorial(n-1)/(&
            gamma_semi(nA + n)*gamma_semi(nB + n))
       FuncVal = rMember
       aPlusI       = A - 1.0
       bPlusI       = B - 1.0
       do i = 1, n-1
          aPlusI = aPlusI + 1.0
          bPlusI = bPlusI + 1.0
          rMember = rMember*aPlusI*bPlusI/(real(i)*real(n - i))*(-OneMinusZ)
          FuncVal = FuncVal + rMember
       end do
       aPlusI = aPlusI + 1.0
       bPlusI = bPlusI + 1.0
       cPlusI = real(n)
       rMember = rMember*aPlusI*bPlusI/cPlusI*(-OneMinusZ)
       LogFactor    = -log(OneMinusZ) + psi_int(1) + psi_int(1 + n) - &
            (psi_semi(nA + n) + psi_semi(nB + n))
       FuncVal = FuncVal + rMember*LogFactor
       i = 0
       do while(abs(rMember) >= cTolerance_I(iRealPrec))
          i = i + 1
          aPlusI = aPlusI + 1.0
          bPlusI = bPlusI + 1.0
          cPlusI = cPlusI + 1.0
          rMember = rMember*aPlusI*bPlusI/(real(i)*cPlusI)*OneMinusZ
          LogFactor = LogFactor + 1.0/real(i) + 1.0/cPlusI - &
               1.0/aPlusI - 1.0/bPlusI
          FuncVal = FuncVal + rMember*LogFactor
       end do

    end if
  end function hyper_semi_semi_int
  !============================================================================
  !=======================|Toroidal functions|=================================
  !                       v                  v
  !
  ! Calculate functions
  !
  ! sqrt(2\sinh u)P^{-1}_{n - 1/2}(\cosh u)/(k^3(k^\prime)^n)
  ! and
  ! sqrt(2\sinh u)Q^{-1}_{n - 1/2}(\cosh u)/(k^3(k^\prime)^n)
  !
  ! The multiplier, k^3, is present in a definition
  ! of the toroid functions and while calculated detivatives
  ! this multiplier is taken into account, however, it is not
  ! calculated here, to avoid dividing zero by zero at k=0
  !
  ! Herewith, cosh u = (2 - k^2)/(2k^\prime);
  ! sinh u = k^2/(2k^\prime)
  ! k^\prime=sqrt(1 - k^2)=exp(-u)
  !
  real function toroid_p(n, Kappa2In, KappaPrime2In)
    ! tilde{P} function, integer n > 0 or n=0, related to k**3 * (k^\prime)**n
    integer, intent(in):: n
    real, optional, intent(in) :: Kappa2In, KappaPrime2In
    real :: Kappa2, KappaPrime2
    character(len=*), parameter:: NameSub = 'toroid_p'
    !--------------------------------------------------------------------------
    if(n < 0)call CON_stop(&
         NameSub//': argument n should be non-negative')
    if(present(Kappa2In))then
       toroid_p = 0.250*hyper_semi_semi_int(nA=1, nB=1+n, nC=3, &
            ZIn = Kappa2In)
    elseif(present(KappaPrime2In))then
       toroid_p = 0.250*hyper_semi_semi_int(nA=1, nB=1+n, nC=3, &
            OneMinusZIn = KappaPrime2In)
    else
       call CON_stop(&
            NameSub//': Kappa2In or KappaPrime2InIn should be present')
    end if
  end function toroid_p
  !============================================================================
  real function toroid_q(n, KappaPrime2In, Kappa2In)
    ! tilde{Q}_n function, integer n > 0, related to k**3 * (k^\prime)**n
    integer, intent(in):: n
    real, optional, intent(in) :: KappaPrime2In, Kappa2In
    character(len=*), parameter:: NameSub = 'toroid_q'
    !--------------------------------------------------------------------------
    if(n < 0)call CON_stop(&
         NameSub//': argument n should be non-negative')
    if(present(KappaPrime2In))then
       toroid_q = hyper_semi_semi_int(nA=1, nB=1+n, nC=n+1, &
            ZIn = KappaPrime2In)
    elseif(present(Kappa2In))then
       toroid_q = hyper_semi_semi_int(nA=1, nB=1+n, nC=n+1, &
            OneMinusZIn = Kappa2In)
    else
       call CON_stop(&
            NameSub//': Kappa2In or KappaPrime2InIn should be present')
    end if
    toroid_q = -toroid_q*cSqrtPi*gamma_semi(n-1)/factorial(n)
  end function toroid_q
  !============================================================================
  !                           ^                  ^
  !===========================|Toroidal functions|=============================
  !==============================|Inductances|=================================
  !                              v           v
  real function scr_inductance(KappaPrime2)
    real, intent(in):: KappaPrime2

    ! for a superconducting ring calculate the ratio of inductance
    ! to \mu_0R_\infty as a function of (k^prime)^2. For a given R0 and a,
    ! R_\infty = sqrt(R_0^2-a^2) and k^\prime=a/(R_0+R_\infty)

    ! Inverse of inductances: from a given harmonic harmonic and total
    real:: InvInductance, InvInductanceTotal
    real:: Tolerance

    ! Loop variable
    integer::n
    !--------------------------------------------------------------------------
    InvInductance = toroid_q(0,KappaPrime2In=KappaPrime2)/&
         (0.25*toroid_p(0, KappaPrime2In=KappaPrime2))
    Tolerance = InvInductance*cTolerance_I(iRealPrec)
    InvInductanceTotal = InvInductance
    n = 0
    do while(InvInductance > Tolerance)
       n = n +1
       InvInductance = toroid_q(n, KappaPrime2In=KappaPrime2)/&
            ((0.25 - n*n)*toroid_p(n, KappaPrime2In=KappaPrime2))
       InvInductanceTotal = InvInductanceTotal + 2.0*InvInductance
    end do

    ! Invert and accunf for the common multiplier
    scr_inductance = 2.0*cPi**2/InvInductanceTotal
  end function scr_inductance
  !============================================================================
  real function l0_ext_inductance(KappaPrime2)
    real, intent(in):: KappaPrime2

    ! calculate the ratio of external inductance
    ! to \mu_0R_\infty as a function of (k^prime_0)^2.
    !--------------------------------------------------------------------------
    l0_ext_inductance = 0.50*toroid_p(0, KappaPrime2In=KappaPrime2)*cPi**2/&
         toroid_q(0,KappaPrime2In=KappaPrime2)
  end function l0_ext_inductance
  !============================================================================
  subroutine calc_elliptic_int_1kind(Z, KElliptic)
    real, intent(in):: Z
    real, intent(out):: KElliptic
    character(len=*), parameter:: NameSub = 'calc_elliptic_int_1kind'
    !--------------------------------------------------------------------------
    ! Calculate 2F1(0.5 +0, 0.5 + 0; 1; Z)
    KElliptic = 0.50*cPi*hyper_semi_semi_int(nA=0, nB=0, nC=1, ZIn=Z**2)
  end subroutine calc_elliptic_int_1kind
  !============================================================================
  subroutine calc_elliptic_int_2kind(ArgK,EElliptic)

    real, intent(in):: ArgK
    real, intent(out):: EElliptic

    ! Loop variable
    integer:: i

    ! Misc
    real :: aPlusI, bPlusI, rMember, LogFactor, ArgKPrime2
    real :: OneOver2n2nMinus1
    character(len=*), parameter:: NameSub = 'calc_elliptic_int_2kind'
    !--------------------------------------------------------------------------
    EElliptic = 0.50*cPi*hyper_semi_semi_int(nA=-1, nB=0, nC=1, ZIn=ArgK**2)
  end subroutine calc_elliptic_int_2kind
  !============================================================================
end module ModHyperGeometric
!==============================================================================
module  ModUniformCurrentFilament
  use ModHyperGeometric, ONLY: toroid_p, toroid_q
  implicit none
  PRIVATE  ! Except
  ! Named indexes for the field amplitudes
  integer, parameter  ::  Axial_  = 1, Poloidal_ =  2, Toroidal_  = 3
  ! \kappa^\prime at the  boundary
  real :: KappaPrime0
  !
  ! Constant  factors  to calculate the internal field
  !
  real :: Q1            ! Eq. 36, constant field factor for uniform current
  real :: CurrentE      ! Constant current I_E
  !
  ! Constants determining toroidal field:
  !
  real  :: Q2           ! Eq. 49 constant  torroidl field factor
  real  :: ToroidQ0AtU0 ! Q^{-1}_{-1/2}(u_0)
  public:: set_kappaprime0  ! Set constant coefficients for given KappaPrime0
  public:: get_amplitude_int! Internal field amplitudes
contains
  !============================================================================
  subroutine set_kappaprime0(KappaPrime0In, Kappa02)
    real, intent(in)  :: KappaPrime0In
    real, intent(out) :: Kappa02
    real :: KappaPrime02, Kappa0,  Kappa03
    !--------------------------------------------------------------------------
    KappaPrime0 = KappaPrime0In
    KappaPrime02 = KappaPrime0**2
    Kappa02 = 1 - KappaPrime02
    Kappa0 = sqrt(Kappa02); Kappa03 = Kappa0*Kappa02
    ! Eq. 36, constant field factor for uniform current
    Q1 = 0.125*toroid_p(1,KappaPrime2In=KappaPrime02)/&
         (toroid_q(1,KappaPrime2In=KappaPrime02))
    ! Constant current I_E
    CurrentE = -1/(3*KappaPrime02*Kappa0*toroid_q(1,KappaPrime2In=KappaPrime02))
    ! Constants determining toroidal field:
    ! Eq. 49 constant  torroidl field factor
    Q2 = toroid_p(1,KappaPrime2In=KappaPrime02)  / &
         (4*KappaPrime02*Kappa0*(toroid_q(1,KappaPrime2In=KappaPrime02))**2)
    ! Q^{-1}_{-1/2}(u_0):
    ToroidQ0AtU0 = Kappa03*toroid_q(0,KappaPrime2In=KappaPrime02)
  end subroutine set_kappaprime0
  !============================================================================
  subroutine get_amplitude_int(KappaPrime2In, Amplitude_I)
    real, intent(in)    :: KappaPrime2In
    real, intent(out)   :: Amplitude_I(Axial_:Toroidal_)
    !
    ! Misc
    !
    real :: Kappa, Kappa2, Kappa3
    !--------------------------------------------------------------------------
    Kappa2 = 1 - KappaPrime2In; Kappa = sqrt(Kappa2); Kappa3 = Kappa2*Kappa

    ! Eqs. 35
    Amplitude_I(Axial_)    = Q1*Kappa3*toroid_q(0,KappaPrime2In = &
         KappaPrime2In) + CurrentE
    Amplitude_I(Poloidal_) = 3*Q1*Kappa3*toroid_q(1,    &
         KappaPrime2In=KappaPrime2In)
    Amplitude_I(Toroidal_) = sqrt(Q2*&
         (ToroidQ0AtU0 - Kappa3*toroid_q(0,KappaPrime2In=KappaPrime2In) )   )
  end subroutine get_amplitude_int
  !============================================================================
end module ModUniformCurrentFilament
!==============================================================================
module  ModSurfaceCurrentFilament
  use ModHyperGeometric, ONLY: toroid_p, toroid_q
  implicit none
  PRIVATE  ! Except
  ! Named indexes for the field amplitudes
  integer, parameter  ::  Axial_  = 1, Poloidal_ =  2, Toroidal_  = 3
  ! \kappa^\prime at the  boundary
  real :: KappaPrime0
  !
  ! Constant  factor  to calculate the internal field
  !
  real, public :: Q0        ! Eq. 31, constant field factor for surface current
  public:: set_kappaprime0  ! Set constant coefficients for given KappaPrime0
  public:: get_amplitude_int! Internal field amplitudes
contains
  !============================================================================
  subroutine set_kappaprime0(KappaPrime0In, Kappa02)
    real, intent(in)  :: KappaPrime0In
    real, intent(out) :: Kappa02
    real :: KappaPrime02
    !--------------------------------------------------------------------------
    KappaPrime0  = KappaPrime0In
    KappaPrime02 = KappaPrime0**2
    Kappa02      = 1 - KappaPrime02
    ! Eq. 31, constant field factor for surface current
    Q0 = 0.125*toroid_p(0, KappaPrime2In=KappaPrime02)/&
       toroid_q(0,KappaPrime2In=KappaPrime02)
  end subroutine set_kappaprime0
  !============================================================================
  subroutine get_amplitude_int(KappaPrime2In, Amplitude_I)
    real, intent(in)    :: KappaPrime2In
    real, intent(out)   :: Amplitude_I(Axial_:Toroidal_)
    !
    ! Misc
    !
    real :: Kappa2, Kappa3
    !--------------------------------------------------------------------------
    Kappa2 = 1 - KappaPrime2In; Kappa3 = sqrt(Kappa2)*Kappa2

    ! Eqs. 35
    Amplitude_I(Axial_)    =   Q0*Kappa3*toroid_q(0, KappaPrime2In = &
         KappaPrime2In)
    Amplitude_I(Poloidal_) = 3*Q0*Kappa3*toroid_q(1, KappaPrime2In=  &
         KappaPrime2In)
    Amplitude_I(Toroidal_) = 0.0
  end subroutine get_amplitude_int
  !============================================================================
end module ModSurfaceCurrentFilament
!==============================================================================
module ModCurrentFilament
  use ModUniformCurrentFilament, ONLY: uniform_current_field=>get_amplitude_int
  use ModSurfaceCurrentFilament, ONLY: surface_current_field=>get_amplitude_int
  implicit none
  !
  ! Radius of torooidal  maagnetic axis
  !
  real  :: rInfty, rInfty2
  ! Named indexes for the field amplitudes
  integer, parameter  ::  Axial_  = 1, Poloidal_ =  2, Toroidal_  = 3
  ! \kappa^2 and \kappa^Prime at the  boundary
  real :: Kappa02, KappaPrime0
  !
  ! Logicals determining if we use manufactured  current distrributions
  ! 1. Uniform current form-factor:
  logical :: UseUniformCurrent = .false.
  ! 2. Surface current form-factor
  logical :: UseSurfaceCurrent = .false.
  ! Inductunce  coefficient; in the SI units this is the ratio of the total
  ! inductance of the filament normaalized by \mu_0 R_\infty
  real :: Inductance
contains
  !============================================================================
  subroutine set_filament_geometry(rMinor, rMajor)
    use ModUniformCurrentFilament, ONLY: set_uniform_current=>set_kappaprime0
    use ModSurfaceCurrentFilament, ONLY: set_surface_current=>set_kappaprime0
    use ModHypergeometric, ONLY: l0_ext_inductance

    real, intent(in) :: rMinor, rMajor
    !--------------------------------------------------------------------------
    rInfty2 = rMajor**2 -  rMinor**2; rInfty = sqrt(rInfty2)
    KappaPrime0 = rMinor/(rMajor + rInfty)
    Kappa02  = -1.0
    if(UseUniformCurrent)then
       call set_uniform_current(KappaPrime0, Kappa02)
       ! Calculate inductance depending on the choice of the current form-factor
       ! Inductance includes external and internal field iductances as well as
       ! the torooidal field inductance:
       !                      external               internal  toroidal
       Inductance = l0_ext_inductance(KappaPrime0**2) + 0.250 + 0.50
    end if
    if(UseSurfaceCurrent)then
       call set_surface_current(KappaPrime0, Kappa02)
       ! Only external inductance matters:
       Inductance = l0_ext_inductance(KappaPrime0**2)
    end if
    if(Kappa02 <= 0.0) then
       ! With no foorm-factor, the field is always external, since for any
       ! \kappa one has \kappa^2 < Kappa_0^2=1
       Kappa02 = 1.0; rInfty = rMajor; rInfty2  = rInfty**2
       ! Approximate formula expressed in terms of a/R0 ratio
       Inductance =  log(8.0*rMajor/rMinor) - 1.250
       ! Set \kappa^\prime = a/(2R_0)
       KappaPrime0 = 0.5* rMinor / rMajor
    end if
  end subroutine set_filament_geometry
  !============================================================================
  subroutine external_field(Kappa2In, Amplitude_I)
    use ModHyperGeometric, ONLY: toroid_p, toroid_q
    real, intent(in)    :: Kappa2In
    real, intent(out)   :: Amplitude_I(Axial_:Toroidal_)
    !--------------------------------------------------------------------------
    ! Eqs. 35
    Amplitude_I(Axial_)    =   toroid_p(0, Kappa2In = Kappa2In)
    Amplitude_I(Poloidal_) = 3*toroid_p(1, Kappa2In = Kappa2In)
    Amplitude_I(Toroidal_) = 0.0
  end subroutine external_field
  !============================================================================
end module ModCurrentFilament
!==============================================================================
module ModFieldGS
  use ModCurrentFilament
  implicit none
  ! Magnitude of the magnetic field at the center of filament:
  real :: Bc         = 1.0   ! In SI, Bc = \mu_0 I / (2 R_\infty)
  ! Magnetic field vector at the center of filament: Bc times
  ! the unit vector of the axis direction
  real :: Bc_D(3)    = [1.0, 0.0, 0.0]
  ! The sign of ratio between always positive toroidal current density  and
  ! the toridal magnetic field, which may be either positive or negative
  ! everywhere
  real :: Helicity   = 1.0
  !  Strapping field. Its action on the filament current balances
  !  the hoop force.
  real :: BStrap_D(3)
contains
  !============================================================================
  subroutine set_filament_field(BcIn, BDir_D)
    use ModNumConst,       ONLY:  cTwoPi
    ! Inputs:
    ! Combined paraameter: its magnititude is Bc, the sign being helicity.
    real, intent(in)  :: BcIn
    ! Unit direction vector of the axis of symmetry, equaal to (-1, 0, 0) in
    ! the coordinate frame, z, x, y  (or z, r, \phi) with z-axis colinear
    ! with the direction of axis of symmetry.
    real, intent(in) :: BDir_D(3)
    !--------------------------------------------------------------------------
    Bc = abs(BcIn); Helicity = sign(1.0, BcIn)
    Bc_D = BDir_D*Bc
    ! With the earlier calculated inductance, determine the strapping field:
    BStrap_D = - (Inductance / cTwoPi)*Bc_D ! Eq. 57 in GS22
  end subroutine set_filament_field
  !============================================================================
  subroutine get_filament_field(R_D, B_D, BetaIn, BPhiOut)
    ! Caalculates:
    ! 1. With preset
    !      UseUniformCurrent = .true.
    !      call set_filament_geometry(a, R0)
    !      call set_filament_field(Bc, BDir_D)
    !    this subroutine calculates the field vector B_D at the point R_D for
    !    the unoform current form-factor (GS22)
    ! 2. With preset UseUniformCurrent = .false. (default) this subroutine
    !    calculates the external field in the TD solution. It should not be
    !    applied inside the filament
    use ModCoordTransform, ONLY: cross_product
    real, intent(in ) :: R_D(3)
    real, intent(out) :: B_D(3)
    real, optional, intent(in )  ::  BetaIn
    real, optional, intent(out)  ::  BPhiOut
    !
    ! Misc: geometry
    !
    real :: R2, Z2, rPerp2, rPerp, RPlus2, Kappa2, KappaPrime2, RMinus2
    !
    ! Msic: field charaacteristics
    !
    real :: BcDotR, CommonFactor, BPhi
    real :: Poloidal_D(3)  !  Not a unit vector, its length is Bc*KappaPrime
    real :: Amplitude_I(Axial_:Toroidal_)
    !--------------------------------------------------------------------------
    !
    ! Geometry
    !
    R2 = sum(R_D**2)
    BcDotR  = sum(Bc_D*R_d)
    Z2= (BcDotR/Bc)**2
    rPerp2  = R2 - Z2;   rPerp =  sqrt(rPerp2)
    RPlus2 = Z2 +  (rInfty + rPerp)**2
    Poloidal_D = ( 2*BcDotR*R_D + (rInfty2 - R2)*Bc_D )/RPlus2
    Kappa2 = 4*rPerp*rInfty/RPlus2
    if(Kappa2 <= Kappa02) then
       ! external field
       call external_field(Kappa2, Amplitude_I)
       CommonFactor = (sqrt(rInfty2/RPlus2))**3
       B_D = CommonFactor*(Bc_D*Amplitude_I(Axial_) + &
            Poloidal_D*Amplitude_I(Poloidal_))
       ! No toroidal external field:
       if(present(BPhiOut)) BPhiOut = 0.0
    else
       ! internal field
       !
       ! Calculate (\kappa^\prime)^2:
       ! 1. Calculate R_-.
       RMinus2 = Z2 +  (rInfty - rPerp)**2
       ! Exprress (\kappa^\prime)^2 via R_- (not via Kappa2 = 1 - KappaPrime2):
       KappaPrime2 = RMinus2 / RPlus2
       ! Field harmonics
       call uniform_current_field(KappaPrime2, Amplitude_I)
       ! Calculate separately the toroidal field:
       BPhi = Helicity*Amplitude_I(Toroidal_)
       if(present(BetaIn)) BPhi = BPhi/sqrt(1 + BetaIn)
       CommonFactor = (sqrt(rInfty/rPerp))**3
       B_D = CommonFactor*(Bc_D*Amplitude_I(Axial_)    + &
            Poloidal_D*Amplitude_I(Poloidal_) + &
            (BPhi / rPerp)*cross_product(Bc_D, R_D) )
       if(present(BPhiOut)) BPhiOut = CommonFactor*BPhi
    end if
  end subroutine get_filament_field
  !============================================================================
  subroutine test
    use ModPlotFile
    ! Loop variables
    integer:: i, j
    ! Number of points = (2N+1)*(2N+1)
    integer, parameter :: N = 300
    ! Field to show:
    real :: Field_DII(Axial_:Toroidal_, -N:N, -N:N)
    ! Domain size (-ZMax to ZMax)
    real, parameter:: ZMax = 2.0
    ! Grid size
    real, parameter :: DGrid = ZMax/N
    !--------------------------------------------------------------------------
    UseUniformCurrent = .true.;  UseSurfaceCurrent = .false.
    ! Set rInfty = 1 and KappaPrime  = 0.1
    call set_filament_geometry(0.20 / 0.990, 1.01 / 0.990)
    call set_filament_field(1.0, [1.0, 0.0, 0.0])
    do j = -N, N
       do i = -N, N
          call get_filament_field(R_D=[DGrid*i, DGrid*j, 0.0],&
               B_D = Field_DII(:, i, j))
       end do
    end do
    call save_plot_file(NameFile='field_gs22.out', &
         TypeFileIn='ascii', &
         NameVarIn='z r Bz Br BPhi'  , &
         CoordMinIn_D=[-ZMax, -ZMax],&
         CoordMaxIn_D=[ ZMax,  ZMax],&
         StringFormatIn = '(6F12.5)',&
         VarIn_VII = Field_DII )
    ! Test the strapped field:
    do j = -N, N
       do i = -N, N
           Field_DII(:, i, j) = Field_DII(:, i, j) + BStrap_D
       end do
    end do
    call save_plot_file(NameFile='strapped_field.out', &
         TypeFileIn='ascii', &
         NameVarIn='z r Bz Br BPhi'  , &
         CoordMinIn_D=[-ZMax, -ZMax],&
         CoordMaxIn_D=[ ZMax,  ZMax],&
         StringFormatIn = '(6F12.5)',&
         VarIn_VII = Field_DII )

  end subroutine test
  !============================================================================
end module ModFieldGS
!==============================================================================
