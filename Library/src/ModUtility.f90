!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModUtilities

  ! Simple methods which are used by CON and can be used
  ! by the components of the SWMF too.
  !
  ! F77 and C++ codes need an F90 interface to access these utilities.

  ! These are all the dependencies:
  use ModKind, ONLY: Real4_, Real8_
  use ModMpi,  ONLY: MPI_wtime, MPI_initialized, MPI_abort, &
       MPI_comm_rank, MPI_comm_world
  use ModIoUnit, ONLY: UnitTmp_, io_unit_clean
  use iso_c_binding, ONLY: c_null_char

  implicit none

  private ! except

  public:: make_dir
  public:: check_dir
  public:: fix_dir_name
  public:: open_file
  public:: close_file
  public:: remove_file
  public:: touch_file
  public:: flush_unit
  public:: split_string
  public:: join_string
  public:: upper_case
  public:: lower_case
  public:: string_to_char_array
  public:: char_array_to_string
  public:: write_string_tabs_name
  public:: sleep
  public:: check_allocate
  public:: greatest_common_divisor
  public:: CON_stop
  public:: CON_set_do_test
  public:: test_mod_utility
  public:: norm2
  public:: find_cell
  interface find_cell
     ! Single and double precision coordinates
     module procedure find_cell4, find_cell8
  end interface find_cell

  logical, public :: DoFlush = .true. ! parameter for flush_unit
  logical, public :: DoWriteCallSequence = .false. ! parameter for CON_stop

  ! Parameters for CON_set_do_test
  character(len=200), public :: StringTest = ''
  integer, public            :: iProcTest = 0
  integer, public            :: lVerbose = 1

  character, public:: cTab = char(9)

  interface split_string
     module procedure split_string, split_string_simple
  end interface

  interface join_string
     module procedure join_string, join_string_simple
  end interface

contains
  !============================================================================
  subroutine write_string_tabs_name(String, Name)

    character(len=*), intent(in):: String, Name

    ! Write String and Name with 2 or 3 tabs in between (PARAM.in format)

    !--------------------------------------------------------------------------
    if(len_trim(String) < 8)then
       write(UnitTmp_,'(a)') trim(String)//cTab//cTab//cTab//trim(Name)
    else
       write(UnitTmp_,'(a)') trim(String)//cTab//cTab//trim(Name)
    end if

  end subroutine write_string_tabs_name
  !============================================================================
  subroutine CON_set_do_test(String, DoTest, DoTestMe)

    ! Set DoTest to true if " String " can be found in " StringTest "
    ! If the optional DoTestMe variable is present, it is set to
    ! DoTestMe = DoTest .and. iProc == iProcTest,
    ! where iProcTest is a public variable set to the test processor index.
    ! If only DoTest is present, it behaves like DoTestMe.
    !
    ! Depending on the value of the public variable lVerbose,
    ! different amount of output is printed:
    ! If String matches StringTest or lVerbose==10, the test processor prints
    !    "String CALLED"
    ! If lVerbose==100, all processors print "String CALLED by iProc=",iProc

    character(len=*),  intent(in) :: String
    logical,           intent(out):: DoTest
    logical, optional, intent(out):: DoTestMe

    logical:: IsMpiInitialized
    integer:: iProc=0, iError
    !--------------------------------------------------------------------------
    call MPI_initialized(IsMpiInitialized, iError)

    if(IsMpiInitialized) &
         call MPI_comm_rank(MPI_COMM_WORLD, iProc, iError)

    if(present(DoTestMe))then
       DoTest   = index(' '//StringTest,' '//String//' ') > 0
       DoTestMe = DoTest .and. iProc==iProcTest
    else
       if(iProc==iProcTest)then
          DoTest = index(' '//StringTest,' '//String//' ') > 0
       else
          DoTest = .false.
       end if
    end if
    if((DoTest .or. lVerbose>=10) .and. iProc==iProcTest)then
       write(*,*)String,' CALLED'
    else if(lVerbose>=100)then
       write(*,*)String,' CALLED by iProc=', iProc
    end if

  end subroutine CON_set_do_test
  !============================================================================
  subroutine CON_stop(String, Value1, Value2, Value3, Value4)

    character(len=*), intent(in):: String
    class(*), optional, intent(in):: Value1, Value2, Value3, Value4

    ! Stop execution after the following actions:
    !
    ! Write out error message with processor rank and String.
    ! Write out optional arguments Value1 ... Value4 of arbitrary type.
    ! Close open files.
    ! If public variable DoCreateCallSequence is true, then make a floating
    ! exception to produce call sequence (with NAG compiler in debugging mode).
    ! Abort execution with MPI_abort and stop.

    logical:: IsMpiInitialized
    integer:: iProc=0, nError, iError
    !--------------------------------------------------------------------------
    call MPI_initialized(IsMpiInitialized, iError)

    if(IsMpiInitialized) call MPI_comm_rank(MPI_COMM_WORLD, iProc, iError)

    write(*,*) 'ERROR: aborting execution on processor', iProc, ' with message:'
    write(*,'(a)') 'ERROR: '//String
    if(present(Value1)) call write_value(Value1)
    if(present(Value2)) call write_value(Value2)
    if(present(Value3)) call write_value(Value3)
    if(present(Value4)) call write_value(Value4)

    ! Try closing all open IO units
    call io_unit_clean

    ! Create call sequence if requested
    if(DoWriteCallSequence)then
       write(*,*)'Making floating point exception to write call sequence!'
       write(*,*) sqrt(-1.0-iProc)
    end if

    ! Stop execution
    if(IsMpiInitialized) call MPI_abort(MPI_COMM_WORLD, nError, iError)
    stop

  end subroutine CON_stop
  !============================================================================
  subroutine write_value(Value, StringBefore, Advance)

    ! Write out value of type real, integer, logical, or character.
    ! If StringBefore is present, write it before the value with no advance.
    ! If Advance is present and set to ADVANCE="NO", then the Value is
    ! written out formatted with ADVANCE="NO".

    class(*),         intent(in):: Value
    character(len=*), optional, intent(in):: StringBefore
    character(len=*), optional, intent(in):: Advance

    logical:: DoAdvance
    !--------------------------------------------------------------------------
    DoAdvance = .true.
    if(present(Advance)) DoAdvance = Advance /= "NO"

    if(present(StringBefore)) write(*,'(a)',ADVANCE="NO") StringBefore
    select type (Value)
    type is (real)
       if(DoAdvance)then
          write(*,*) Value
       else
          write(*,'(es18.5)',ADVANCE="NO") Value
       end if
    type is (integer)
       if(DoAdvance)then
          write(*,*) Value
       else
          write(*,'(i10)',ADVANCE="NO") Value
       end if
    type is (logical)
       if(DoAdvance)then
          write(*,*) Value
       else
          write(*,'(l1)',ADVANCE="NO") Value
       end if
    type is (character(*))
       if(DoAdvance)then
          write(*,'(a)') Value
       else
           write(*,'(a)',ADVANCE="NO") Value
        end if
    end select

  end subroutine write_value
  !============================================================================

  subroutine make_dir(NameDir, iErrorOut)

    character(len=*), intent(in):: NameDir ! Directory name
    integer, intent(out), optional:: iErrorOut    ! 0 on success, -1 otherwise

    ! Create the directory using "mkdir -p NameDir"
    ! If the directory already exists, this function does nothing.
    ! The permission can be set by including it into
    ! NameDir, for example NameDir = '-m777 Public'.
    !
    ! If the optional iErrorOut is present, it is set to the error code,
    ! which is 0 if no error occurred. Otherwise the code stops with an
    ! error message if an error occurs.

    integer:: iError = -1             ! Return value as retrieved from shell

    character(len=100):: StringCommand

    ! Check if directory exists
    character(len=*), parameter:: NameSub = 'make_dir'
    !--------------------------------------------------------------------------
    call check_dir(NameDir, iError)
    if(iError == 0)then
       if(present(iErrorOut)) iErrorOut = 0
       RETURN
    end if

    ! Create directory
    StringCommand = 'mkdir -p '//trim(NameDir)
    call execute_command_line(trim(StringCommand), EXITSTAT=iError)

    if(present(iErrorOut)) then
       iErrorOut = iError
       RETURN
    elseif(iError /= 0)then
       write(*,*) NameSub,' iError =', iError
       call CON_stop(NameSub//' failed to make directory '//trim(NameDir))
    end if

  end subroutine make_dir
  !============================================================================

  subroutine check_dir(NameDir, iErrorOut)

    character(len=*), intent(in) :: NameDir
    integer, optional, intent(out):: iErrorOut

    ! Check if a directory exists by trying to open a file in it.
    ! Die with an error message if the directory does not exist.
    ! Directory names are cached, so multiple calls with the same
    ! name do not repeat the check.
    !
    ! {\bf This subroutine should be called by the root PE of the component
    ! only!}
    ! Calling the subroutine from multiple PE-s may result in a fatal error,
    ! namely one PE may delete the file written by the other PE, so the
    ! other PE thinks that the directory does not exist.

    integer:: iError

    ! Try to open a file in this directory
    character(len=*), parameter:: NameSub = 'check_dir'
    !--------------------------------------------------------------------------
    open(UnitTmp_, file=trim(NameDir)//'/.test', status='unknown', &
         iostat = iError)

    if (iError /= 0) then
       if(present(iErrorOut))then
          iErrorOut = iError
          RETURN
       else
          call CON_stop(NameSub//' ERROR: Cannot find/write into directory '&
               //trim(NameDir))
       end if
    else
       close(UnitTmp_, status = 'DELETE')
    endif

  end subroutine check_dir
  !============================================================================

  subroutine fix_dir_name(NameDir)

    character(len=*), intent(inout) :: NameDir

    ! Append a '/' at the end of the directory name if it is not there
    ! and the directory name is not zero length (empty string).
    !
    ! {\bf This subroutine should be called by all PE-s of the component!}

    integer :: i
    character(len=*), parameter:: NameSub = 'fix_dir_name'
    !--------------------------------------------------------------------------
    i = len_trim(NameDir)
    if(i == 0) RETURN
    if(NameDir(i:i) == '/') RETURN

    if(i >= len(NameDir)) call CON_stop(NameSub// &
         "ERROR cannot append / to directory name "//NameDir)

    NameDir(i+1:i+1) = '/'

  end subroutine fix_dir_name
  !============================================================================

  subroutine open_file(iUnitIn, File, Form, Status, Position, Access, Recl, &
       iComm, NameCaller)

    ! Interface for the Fortran open statement with error checking.
    ! If an error occurs, the code stops and writes out the unit number,
    ! the error code and the name of the file.
    ! If NameCaller is present, it is also shown.
    ! If no unit number is passed, open UnitTmp_.
    ! Default format is 'formatted' as in the open statement.
    ! Default status is 'replace' (not unknown) as it is well defined.
    ! Default position is 'rewind' (not asis) as it is well defined.
    ! Default access is 'sequential' as in the open statement.
    ! There is no default record length Recl.
    ! If the MPI communicator iComm is present together with Recl,
    ! the file will be opened with status='replace' on processor 0,
    ! and with status='old' on other processors with an MPI_barrier
    ! in between.

    integer, optional, intent(in):: iUnitIn
    character(len=*), optional, intent(in):: File
    character(len=*), optional, intent(in):: Form
    character(len=*), optional, intent(in):: Status
    character(len=*), optional, intent(in):: Position
    character(len=*), optional, intent(in):: Access
    integer,          optional, intent(in):: Recl
    integer,          optional, intent(in):: iComm
    character(len=*), optional, intent(in):: NameCaller

    character(len=20):: TypeForm, TypeStatus, TypePosition, TypeAccess

    integer:: iUnit
    integer:: iError, iProc, nProc

    character(len=*), parameter:: NameSub = 'open_file'
    !--------------------------------------------------------------------------
    iUnit = UnitTmp_
    if(present(iUnitIn)) iUnit = iUnitIn

    TypeForm = 'formatted'
    if(present(Form)) TypeForm = Form

    TypeStatus = 'replace'
    if(present(Status)) TypeStatus = Status

    TypePosition = 'rewind'
    if(present(Position)) TypePosition = Position

    TypeAccess = 'sequential'
    if(present(Access)) TypeAccess = Access

    if(present(Recl))then
       if(present(iComm))then
          ! Get iProc and nProc
          call MPI_comm_size(iComm, nProc, iError)
          call MPI_comm_rank(iComm, iProc, iError)
          ! Open file with status "replace" on processor 0
          if(iProc == 0) &
               open(iUnit, FILE=File, FORM=TypeForm, STATUS='replace', &
               ACCESS=TypeAccess, RECL=Recl, IOSTAT=iError)
          if(nProc > 1)then
             ! Check if open worked on processor 0
             if(iProc == 0 .and. iError /= 0)then
                write(*,*) NameSub,' iUnit, iError=', iUnit, iError
                if(present(NameCaller)) write(*,*) 'NameCaller=', NameCaller
                call CON_stop(NameSub// &
                     ' processor 0 could not open file='//trim(File))
             end if
             ! Make sure all processors wait until proc 0 has opened file
             call MPI_barrier(iComm, iError)
             ! Other processors open with status "old"
             if(iProc > 0) &
                  open(iUnit, FILE=File, FORM=TypeForm, STATUS='old', &
                  ACCESS=TypeAccess, RECL=Recl, IOSTAT=iError)
          end if
       else
          open(iUnit, FILE=File, FORM=TypeForm, STATUS=TypeStatus, &
               ACCESS=TypeAccess, RECL=Recl, IOSTAT=iError)
       end if
    else
       open(iUnit, FILE=File, FORM=TypeForm, STATUS=TypeStatus, &
            POSITION=TypePosition, ACCESS=TypeAccess, IOSTAT=iError)
    end if

    if(iError /= 0)then
       write(*,*) NameSub,' iUnit, iError=', iUnit, iError
       if(present(NameCaller)) write(*,*) 'NameCaller=', NameCaller
       call CON_stop(NameSub//' could not open file='//trim(File))
    end if

  end subroutine open_file
  !============================================================================

  subroutine close_file(iUnitIn, Status, NameCaller)

    ! Interface for the Fortran close statement with error checking
    ! If an error occurs, the code stops and writes out the unit number,
    ! the error code and the name of the file.
    ! If NameCaller is present, it is also shown.
    ! If no unit number is passed, close UnitTmp_

    integer, optional, intent(in):: iUnitIn
    character(len=*), optional, intent(in):: Status
    character(len=*), optional, intent(in):: NameCaller

    integer:: iUnit
    integer:: iError

    character(len=*), parameter:: NameSub = 'close_file'
    !--------------------------------------------------------------------------
    iUnit = UnitTmp_
    if(present(iUnitIn)) iUnit = iUnitIn

    if(present(Status))then
       close(iUnit, STATUS=Status, IOSTAT=iError)
    else
       close(iUnit, IOSTAT=iError)
    end if

    if(iError /= 0)then
       write(*,*) NameSub,' iUnit, iError=', iUnit, iError
       if(present(NameCaller)) write(*,*) 'NameCaller=', NameCaller
       call CON_stop(NameSub//' could not close unit')
    end if

  end subroutine close_file
  !============================================================================

  subroutine remove_file(NameFile, NameCaller)

    ! Remove file NameFile if it exists
    ! Pass NameCaller to open_file and close_file in case of errors

    character(len=*), intent(in):: NameFile
    character(len=*), optional, intent(in):: NameCaller

    logical:: IsFound
    !--------------------------------------------------------------------------

    inquire(FILE=NameFile, EXIST=IsFound)
    if(.not.IsFound) RETURN

    call open_file(FILE=NameFile, NameCaller=NameCaller)
    call close_file(STATUS='DELETE', NameCaller=NameCaller)

  end subroutine remove_file
  !============================================================================

  subroutine touch_file(NameFile, NameCaller)

    ! Create file NameFile if it does not exist
    ! Pass NameCaller to open_file and close_file in case of errors

    character(len=*), intent(in):: NameFile
    character(len=*), optional, intent(in):: NameCaller

    !--------------------------------------------------------------------------
    call open_file(FILE=NameFile, NameCaller=NameCaller)
    call close_file(NameCaller=NameCaller)

  end subroutine touch_file
  !============================================================================

  subroutine flush_unit(iUnit)

    integer, intent(in) :: iUnit

    ! Flush the I/O unit iUnit if DoFlush is true
    !--------------------------------------------------------------------------
    if(DoFlush) flush(iUnit)

  end subroutine flush_unit
  !============================================================================

  subroutine split_string_simple(String, String_I, nStringOut, &
       StringSepIn, UseArraySyntaxIn, DoAddSeparator)

    character(len=*),    intent(in):: String    ! string to be split
    character (len=*), intent(out):: String_I(:) ! array of substrings
    integer,          optional, intent(out):: nStringOut ! number of substrings
    character(len=*), optional, intent(in):: StringSepIn ! separator string
    logical, optional, intent(in):: UseArraySyntaxIn     ! expand Var(10:20:2)
    logical, optional, intent(in):: DoAddSeparator ! add sep. to substrings
    !--------------------------------------------------------------------------
    call split_string(String, size(String_I), String_I, nStringOut, &
         StringSepIn, UseArraySyntaxIn, DoAddSeparator)

  end subroutine split_string_simple
  !============================================================================

  subroutine split_string(String, MaxString, String_I, nStringOut, &
       StringSepIn, UseArraySyntaxIn, DoAddSeparator)

    character(len=*),    intent(in):: String    ! string to be split
    integer,             intent(in):: MaxString ! maximum array size
    character (len=*), intent(out):: String_I(MaxString) ! array of substrings
    integer,          optional, intent(out):: nStringOut ! number of substrings
    character(len=*), optional, intent(in):: StringSepIn ! separator string
    logical, optional, intent(in):: UseArraySyntaxIn     ! expand Var(10:20:2)
    logical, optional, intent(in):: DoAddSeparator ! add sep. to substrings

    ! Cut the input string into an array of substrings. The separator
    ! string is either StringSepIn or a single space (default).
    ! Multiple consecutive separator strings are treated as one.
    ! Leading and trailing spaces are ignored. For example
    !\begin{verbatim}
    ! ' IE  GM ' --> nString=2, String\_I=(/'IE','GM'/)
    !\end{verbatim}
    ! When UseArraySyntax is present, then expand strings containing
    ! parens into an array of substrings ending with numbers, e.g.
    !\begin{verbatim}
    ! 'Var(4)'      --> nString=4,  String\_I=(/'Var1','Var2','Var3','Var4'/)
    ! 'Var(11)'     --> nString=11, String\_I=(/'Var01','Var02',...,'Var11'/)
    ! 'Var(3:5)'    --> nString=3,  String\_I=(/'Var3','Var4','Var5'/)
    ! 'Var(7:11:2)' --> nString=3,  String\_I=(/'Var07','Var09','Var11'/)
    !\end{verbatim}

    integer          :: nString
    character(len=10):: StringSep
    logical          :: UseArraySyntax

    character(len=len(String)+10) :: StringTmp
    character(:), allocatable :: StringTmp2

    integer :: i, l, lSep, lKeep

    character(len=*), parameter:: NameSub = 'split_string'
    !--------------------------------------------------------------------------
    if(present(StringSepIn))then
       StringSep = StringSepIn
       lSep      = len(StringSepIn)
    else
       StringSep = ' '
       lSep      = 1
    end if

    UseArraySyntax = .false.
    if(present(UseArraySyntaxIn)) UseArraySyntax = UseArraySyntaxIn

    lKeep = 0
    if(present(DoAddSeparator))then
       if(DoAddSeparator) lKeep = lSep
    end if

    nString   = 0
    StringTmp = String
    l         = len_trim(StringTmp)
    StringTmp = trim(StringTmp) // StringSep(1:lSep) ! Add separator to the end
    do
       StringTmp = adjustl(StringTmp)          ! Remove leading spaces
       i = index(StringTmp, StringSep(1:lSep)) ! Find end of first part
       if(i <= 1) EXIT                         ! Nothing before the separator
       nString = nString + 1                   ! Count parts

       if(lKeep>0)then                         ! Do not keep added separator
          if(i+lKeep >= len_trim(StringTmp)) lKeep = 0
       end if

       String_I(nString) = StringTmp(1:i-1+lKeep) ! Put part into string array
       StringTmp2= StringTmp(i+lSep:l+lSep) ! Delete part+separator from string
       StringTmp = StringTmp2               ! Delete part+separator from string

       if(UseArraySyntax) call expand_array(String_I(nString))

       if(nString == MaxString) EXIT        ! Check for maximum number of parts
    end do

    if(present(nStringOut)) nStringOut = nString

  contains
    !==========================================================================
    subroutine expand_array(String1)

      ! Expand String1 if it contains array syntax, e.g.
      ! "Var(04)"     to   "Var01", "Var02", "Var03", "Var04"
      ! "Var(2:4)"    to   "Var2", "Var3", "Var4"
      ! "Var(8:12:2)" to   "Var08", "Var10", "Var12"

      character(len=*), intent(inout):: String1

      character(len=len(String1)) :: String2
      character(len=6):: StringFormat

      integer:: j, k, l, m, lNum, iFirst, iLast, Di, iNum, iError
      !------------------------------------------------------------------------
      ! Find the opening paren if any
      j = index(String1,'(')
      if(j < 1) RETURN
      k = index(String1,')')

      if(k < j) &
           call CON_stop(NameSub//' missing closing paren in String='//String)

      ! Check for colon
      l = index(String1,':')
      if(l > j)then
         ! read initial index value before the first colon
         read(String1(j+1:l-1),*,IOSTAT=iError) iFirst
         if(iError /= 0 .or. iFirst < 1) call CON_stop(NameSub// &
              ' invalid initial index value in String='//String)
      else
         iFirst = 1
         l = j
      end if

      ! Check for a second colon
      m = index(String1,':',back=.true.)
      if(m > l)then
         ! read index stride value after the seecond colon
         read(String1(m+1:k-1),*,IOSTAT=iError) Di
         if(iError /= 0 .or. Di < 1) call CON_stop(NameSub// &
              ' invalid index stride value in String='//String)
      else
         Di = 1
         m  = k
      end if

      ! read the last index value between the l and m positions
      read(String1(l+1:m-1),*,IOSTAT=iError) iLast
      if(iError /= 0 .or. iLast < iFirst) call CON_stop(NameSub// &
           ' invalid maximum index value in String='//String)

      ! Set length of numerical string and the format string
      lNum = m - l - 1
      write(StringFormat,'(a,i1,a,i1,a)') "(i",lNum,".",lNum,")"

      ! Set the beginning part of the string to the variable name
      String2 = ' '
      String2(1:j-1) = String1(1:j-1)

      ! Expand variable names by repating name and adding numerical value
      nString = nString - 1
      do iNum = iFirst, iLast, Di

         write(String2(j:j+lNum),StringFormat) iNum
         nString = nString + 1
         String_I(nString) = String2

         if(nString == MaxString) RETURN

      end do
    end subroutine expand_array
    !==========================================================================

  end subroutine split_string
  !============================================================================

  subroutine join_string_simple(String_I, String, StringSepIn)

    character(len=*), intent(in) :: String_I(:)
    character(len=*), intent(out):: String
    character(len=*), optional, intent(in):: StringSepIn
    !--------------------------------------------------------------------------
    call join_string(size(String_I), String_I, String, StringSepIn)

  end subroutine join_string_simple
  !============================================================================

  subroutine join_string(nString, String_I, String, StringSepIn)

    integer,          intent(in) :: nString
    character(len=*), intent(in) :: String_I(nString)
    character(len=*), intent(out):: String
    character(len=*), optional, intent(in):: StringSepIn

    ! Join the input string array into one string. The parts are joined
    ! with spaces or the optional StringSepIn string (up to 10 characters).

    character(len=10):: StringSep

    integer :: i, l

    character(len=*), parameter:: NameSub = 'join_string'
    !--------------------------------------------------------------------------
    if(present(StringSepIn))then
       StringSep = StringSepIn
       l = len(StringSepIn)
    else
       StringSep = ' '
       l = 1
    endif
    String = String_I(1)

    do i = 2, nString
      String = trim(String) // StringSep(1:l) // String_I(i)
    end do

  end subroutine join_string
  !============================================================================

  subroutine upper_case(String)

    character (len=*), intent(inout) :: String

    ! Change characters to upper case in String

    integer, parameter :: iA=ichar('a'), iZ=ichar('z'), Di=ichar('A')-iA
    integer :: i, iC
    !--------------------------------------------------------------------------
    do i = 1, len_trim(String)
       iC = ichar(String(i:i))
       if(iC >= iA .and. iC <= iZ) String(i:i) = char(iC+Di)
    end do

  end subroutine upper_case
  !============================================================================

  subroutine lower_case(String)

    character (len=*), intent(inout) :: String

    ! Change characters to lower case in String

    integer, parameter :: iA=ichar('A'), iZ=ichar('Z'), Di=ichar('a')-iA
    integer :: i, iC
    !--------------------------------------------------------------------------
    do i = 1, len_trim(String)
       iC = ichar(String(i:i))
       if(iC >= iA .and. iC <= iZ) String(i:i) = char(iC+Di)
    end do

  end subroutine lower_case
  !============================================================================

  subroutine string_to_char_array(String, String_I)

    ! Convert Fortran string into a C-style character array
    ! Ignore trailing spaces.
    ! Add null character to the end. Return length if needed.

    character(len=*),  intent(in) :: String
    character,         intent(out):: String_I(*)

    integer:: i, n
    !--------------------------------------------------------------------------
    n = len_trim(String)
    do i = 1, n
       String_I(i) = String(i:i)
    end do
    String_I(n+1) = c_null_char

  end subroutine string_to_char_array
  !============================================================================

  subroutine char_array_to_string(String_I, String)

    ! Convert C-style character array into a Fortran string.

    character,         intent(in) :: String_I(*)
    character(len=*),  intent(out):: String

    integer:: i
    !--------------------------------------------------------------------------
    String = ' '
    do i = 1, len(String)
       if(String_I(i) == c_null_char) EXIT
       String(i:i) = String_I(i)
    end do

  end subroutine char_array_to_string
  !============================================================================

  subroutine sleep(DtCpu)

    real, intent(in) :: DtCpu  ! CPU time to sleep (in seconds)

    real(Real8_) :: tCpu0

    ! This subroutine returns after the number of seconds
    ! given in its argument.
    !--------------------------------------------------------------------------
    tCpu0 = MPI_WTIME()
    do
       if(MPI_WTIME() > tCpu0 + DtCpu) RETURN
    end do

  end subroutine sleep
  !============================================================================

  subroutine check_allocate(iError, NameArray)

    integer,intent(in)::iError
    character(LEN=*),intent(in)::NameArray

    !--------------------------------------------------------------------------
    if (iError > 0) call CON_stop('check_allocate F90_ERROR '// &
         'Could not allocate array '//NameArray)

  end subroutine check_allocate
  !============================================================================

  recursive function greatest_common_divisor(i, j) result(kGCD)
    ! Calculate the greatest common divisor of i and j
    ! with Euclid's algorithm

    integer, intent(in):: i, j
    integer:: kGCD
    !--------------------------------------------------------------------------
    if(j == 0)then
       kGCD = i
    else
       kGCD = greatest_common_divisor(j, mod(i, j))
    end if

  end function greatest_common_divisor
  !============================================================================
  subroutine test_mod_utility

    ! Unit tests

    character(len=500):: String
    character:: StringC_I(501)
    integer, parameter :: MaxString = 20
    integer :: nString
    character(len=30) :: String_I(MaxString)
    integer :: iString, l
    integer:: iError

    character(len=*), parameter:: NameSub = 'test_mod_utility'
    !--------------------------------------------------------------------------
    write(*,'(a)') 'testing check_dir'
    write(*,'(a)') 'check directory "."'
    call check_dir('.')
    write(*,'(a)') 'check_dir returned successfully'
    write(*,'(a)') 'check directory "xxx/"'
    call check_dir('xxx/', iErrorOut=iError)
    if(iError==0) write(*,*) 'ERROR: iError should not be zero!'

    write(*,'(a)') 'testing make_dir'
    write(*,'(a)') 'make directory "xxx/"'
    call make_dir('xxx')
    call check_dir('xxx/')
    write(*,'(a)') 'make xxx/ directory again (should not produce an error)'
    call make_dir('xxx', iErrorOut=iError)
    write(*,*) iError

    write(*,'(a)') 'testing open_file and close_file'

    ! Use defaults
    call open_file(FILE='xxx/testfile.dat')
    write(UnitTmp_,'(a)') 'Some text'
    call close_file
    ! Create an error message by passing incorrect filename
    ! Since the error code varies by compiler, this is commented out
    ! call open_file(FILE='xxx/testfile.bad', STATUS='old', &
    !     NameCaller=NameSub)

    ! Use all arguments
    call open_file(UnitTmp_, FILE='xxx/testfile.dat', FORM='formatted', &
         STATUS='old', POSITION='append', NameCaller=NameSub)
    write(UnitTmp_,'(a)') 'More text'
    call close_file(UnitTmp_, STATUS='delete', NameCaller=NameSub)

    write(*,'(a)') 'testing touch_file and remove_file'
    call touch_file('xxx/touched', NameCaller=NameSub)
    call touch_file('xxx/touched', NameCaller=NameSub)
    call remove_file('xxx/touched', NameCaller=NameSub)
    call remove_file('xxx/touched', NameCaller=NameSub)

    write(*,'(/,a)') 'testing fix_dir_name'
    String = ' '
    call fix_dir_name(String)
    write(*,'(a)') 'fixed empty string=' // trim(String)
    String = 'GM/BATSRUS/data'
    write(*,'(a)') 'original    string=' // trim(String)
    call fix_dir_name(String)
    write(*,'(a)') 'fixed first string=' // trim(String)
    call fix_dir_name(String)
    write(*,'(a)') 'fixed again string=' // trim(String)

    write(*,'(/,a)') 'testing split_string'
    String = '  a(3)  bb(04:06) c,ddd ee,ff gg(8:12:2) '
    write(*,'(a)') 'String=' // trim(String)

    call split_string(String, MaxString, String_I, nString)
    write(*,'(a,i3,a)') 'with space separator split to', nString, ' parts:'
    do iString = 1, nString
       write(*,'(a)') trim(String_I(iString))
    end do

    call split_string(String, String_I, nString, ",")
    write(*,'(a,i3,a)') 'with "," separator split to', nString, ' parts:'
    do iString = 1, nString
       write(*,'(a)') trim(String_I(iString))
    end do

    call split_string(String, String_I, nString, ") ", DoAddSeparator=.true.)
    write(*,'(a,i3,a)') 'with ") " separator split to', nString, ' parts:'
    do iString = 1, nString
       write(*,'(a)') trim(String_I(iString))
    end do

    call split_string(String, String_I, nString, UseArraySyntaxIn=.true.)
    write(*,'(a,i3,a)') 'with UseArraySyntax split to', nString,' parts:'
    do iString = 1, nString
       write(*,'(a)') trim(String_I(iString))
    end do

    write(*,'(/,a)') 'testing join_string'
    call join_string(nString, String_I, String)
    write(*,'(a)') 'joined string='//trim(String)

    call join_string(String_I(1:4), String, '; ')
    write(*,'(a)') 'joined string='//trim(String)

    write(*,'(/,a)') 'testing upper_case and lower_case'
    String = 'abCD 123:'
    write(*,'(a)') 'mixed case string='//trim(String)
    call upper_case(String)
    write(*,'(a)') 'upper case string='//trim(String)
    call lower_case(String)
    write(*,'(a)') 'lower case string='//trim(String)

    write(*,'(/,a)') 'testing string_to_char_array'
    String = "it's a string"
    call string_to_char_array(String, StringC_I)
    l = len_trim(String)
    write(*,'(a,100a1)')'C character array: ', StringC_I(1:l)
    if(StringC_I(l+1) /= c_null_char) &
         write(*,*)'Error: null terminator is missing'

    write(*,'(/,a)') 'testing char_array_to_string'
    call char_array_to_string(StringC_I, String)
    write(*,'(a,a)')    'Fortran string   : ', trim(String)
    if(String /= "it's a string") &
         write(*,*)'Error: incorrect conversion to Fortran String'

    write(*,'(/,a)') 'testing greatest_common_divisor'
    l = greatest_common_divisor(36, 26)
    if(l /= 2) &
         write(*,*)'Error: greatest_common_divisor(36,26)=', l,' should be 2'

    l = greatest_common_divisor(26, 36)
    if(l /= 2) &
         write(*,*)'Error: greatest_common_divisor(26,36)=', l,' should be 2'

    l = greatest_common_divisor(36, 12)
    if(l /= 12) &
         write(*,*)'Error: greatest_common_divisor(36,12)=', l,' should be 12'

    l = greatest_common_divisor(12, 36)
    if(l /= 12) &
         write(*,*)'Error: greatest_common_divisor(12,36)=', l,' should be 12'

  end subroutine test_mod_utility
  !============================================================================
  function norm2(x_I) result(norm)
    !$acc routine seq

    ! Does the same as the intrinsic norm2 function in Fortran 2008
    ! It is needed because nvfortran+ACC does not have norm2 implemented

    real, intent(in) :: x_I(:)
    real:: norm
    !--------------------------------------------------------------------------
    norm = sqrt(sum(x_I**2))

  end function norm2
  !============================================================================
  subroutine find_cell4(MinCoord, MaxCoord, Coord, iCoord, dCoord, &
       Coord_I, DoExtrapolate, StringError, IsInside)

    ! Single precision coordinates.
    !
    ! Find cell index and distance from cell for either
    ! - a uniform grid with normalized coordinate (Coord_I is NOT present)
    ! - a nonuniform grid with monotone coordinates (Coord_I is present)
    !
    ! For sake of easy usage the returned coordinate index iCoord always
    ! satisfies MinCoord <= iCoord < MaxCoord.
    !
    ! If the coordinate is out of bounds, and DoExtrapolate is not present,
    ! the code stops with an error message. If DoExtrapolate is present and
    ! false, dCoord is modified to 0 or 1 so that the last grid cell is used.
    ! If DoExtrapolate is true, iCoord and dCoord are set
    ! corresponding to linear extrapolation from the last two grid values.
    !
    ! For interpolation the normalized distance dCoord measured from
    ! coordinate iCoord satisfies 0.0 <= dCoord <= 1.0
    ! but for extrapolation dCoord < 0.0 or > 1.0 is also possible.
    !---
    ! FOR THE UNIFORM CASE the normalized coordinate Coord should be equal to
    ! the index at the cell centers, therefore:
    !
    ! iCoord = max(MinCoord, min(MaxCoord-1, floor(Coord)))
    ! dCoord = Coord - iCoord
    !
    ! The optional IsInside = MinCoord <= Coord <= MaxCoord

    ! IN THE NON-UNIFORM CASE the cell iCoord that is left to coordinate Coord
    ! is found with a binary search in the Coord_I coordinates.
    !
    ! The normalized distance is set to
    !    dCoord = (Coord-Coord_I(iCoord))/(Coord_I(iCoord+1)-Coord_I(iCoord))
    !
    ! The optional IsInside = Coord_I(1) <= Coord <= Coord_I(nCoord).

    ! Example for linear interpolation on a 1D uniform grid of nX cells,
    ! DeltaX cell size and the first cell center is at DeltaX/2:
    !
    !   call find_cell(1, nX, x/DeltaX+0.5, iX, d)
    !   State_V = (1.0 - d)*State_VC(:,iX) + d*State_VC(:,iX+1)
    !
    ! Example for linear interpolation on a 1D non-uniform grid
    ! with 2 ghost cells:
    !
    !   call find_cell(-1, nI+2, x, iX, d, x_G)
    !   State_V = (1.0 - d)*State_VG(:,iX) + d*State_VG(:,iX+1)

    integer,          intent(in)           :: MinCoord, MaxCoord
    real(Real4_),     intent(in)           :: Coord
    integer,          intent(out)          :: iCoord
    real,             intent(out), optional:: dCoord
    real(Real4_),     intent(in),  optional:: Coord_I(MinCoord:)
    logical,          intent(in),  optional:: DoExtrapolate
    character(len=*), intent(in),  optional:: StringError
    logical,          intent(out), optional:: IsInside

    integer:: i, Di

    logical:: IsUniform

    character(len=*), parameter:: NameSub = 'find_cell4'
    !--------------------------------------------------------------------------

    if(present(Coord_I)) then
       IsUniform = size(Coord_I) < 2
    else
       IsUniform = .true.
    endif

    if(IsUniform)then
       ! Uniform grid case with normalized coordinate

       iCoord = min(MaxCoord-1, max(MinCoord, floor(Coord)))
       dCoord = Coord - iCoord

       if(Coord < MinCoord)then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) NameSub, ': ', StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord=', Coord
             call CON_stop(NameSub//': normalized coordinate is too small!')
          elseif(.not.DoExtrapolate)then
             ! Use lefttmost cell (first order accurate)
             dCoord = 0.0
          end if
          if(present(IsInside)) IsInside = .false.
       elseif(Coord > MaxCoord)then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord=', Coord
             call CON_stop(NameSub//': normalized coordinate is too large!')
          elseif(.not.DoExtrapolate)then
             ! Use rightmost cell (first order accurate)
             dCoord = 1.0
          endif
          if(present(IsInside)) IsInside = .false.
       else
          if(present(IsInside)) IsInside = .true.
       end if

    elseif(Coord_I(MinCoord) < Coord_I(MaxCoord))then

       ! Monotone increasing coordinates

       if(Coord < Coord_I(MinCoord))then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord, CoordMin, CoordMax=', &
                  Coord, Coord_I(MinCoord), Coord_I(MaxCoord)
             call CON_stop(NameSub//': coordinate is too small!')
          elseif(DoExtrapolate)then
             iCoord = MinCoord
             dCoord = (Coord - Coord_I(iCoord)) &
                  /   (Coord_I(iCoord+1) - Coord_I(iCoord))
          else
             iCoord = MinCoord
             dCoord = 0.0
          end if
          if(present(IsInside)) IsInside = .false.
          RETURN
       end if

       if(Coord > Coord_I(MaxCoord))then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord, CoordMin, CoordMax=', &
                  Coord, Coord_I(MinCoord), Coord_I(MaxCoord)
             call CON_stop(NameSub//': coordinate is too large!')
          elseif(DoExtrapolate)then
             iCoord = MaxCoord - 1
             dCoord = (Coord - Coord_I(iCoord))  &
                  /   (Coord_I(iCoord+1) - Coord_I(iCoord))
          else
             iCoord = MaxCoord - 1
             dCoord = 1.0
          end if
          if(present(IsInside)) IsInside = .false.
          RETURN
       end if

       if(present(IsInside)) IsInside = .true.

       ! binary search
       i  = (MinCoord + MaxCoord)/2
       Di = (MaxCoord - MinCoord)/2
       do
          Di = (Di + 1)/2
          if(Coord < Coord_I(i)) then
             i = max(MinCoord, i - Di)
          elseif(Coord > Coord_I(i+1))then
             i = min(MaxCoord-1, i + Di)
          else
             EXIT
          end if
       end do
       iCoord = i
       if(Coord_I(iCoord+1) == Coord_I(iCoord))then
          dCoord = 0.0
       else
          dCoord = (Coord             - Coord_I(iCoord)) &
               /   (Coord_I(iCoord+1) - Coord_I(iCoord))
       end if
    else

       ! Monotone decreasing coordinates

       if(Coord < Coord_I(MaxCoord))then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord, CoordMin, CoordMax=', &
                  Coord, Coord_I(MaxCoord), Coord_I(MinCoord)
             call CON_stop(NameSub//': coordinate is too small!')
          elseif(DoExtrapolate)then
             iCoord = MaxCoord - 1
             dCoord = (Coord_I(iCoord) - Coord) &
                  /   (Coord_I(iCoord) - Coord_I(iCoord+1))
          else
             iCoord = MaxCoord - 1
             dCoord = 1.0
          end if
          if(present(IsInside)) IsInside = .false.
          RETURN
       end if

       if(Coord > Coord_I(MinCoord))then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord, CoordMin, CoordMax=', &
                  Coord, Coord_I(MaxCoord), Coord_I(MinCoord)
             call CON_stop(NameSub//': coordinate is too large!')
          elseif(DoExtrapolate)then
             iCoord = MinCoord
             dCoord = (Coord_I(iCoord) - Coord)  &
                  /   (Coord_I(iCoord) - Coord_I(iCoord+1))
          else
             iCoord = MinCoord
             dCoord = 0.0
          end if
          if(present(IsInside)) IsInside = .false.
          RETURN
       end if

       if(present(IsInside)) IsInside = .true.

       ! binary search
       i  = (MinCoord + MaxCoord)/2
       Di = (MaxCoord - MinCoord)/2
       do
          Di = (Di + 1)/2
          if(Coord > Coord_I(i)) then
             i = max(MinCoord, i - Di)
          elseif(Coord < Coord_I(i+1))then
             i = min(MaxCoord-1, i + Di)
          else
             EXIT
          end if
       end do
       iCoord = i
       if(Coord_I(iCoord+1) == Coord_I(iCoord))then
          dCoord = 0.0
       else
          dCoord = (Coord_I(iCoord) - Coord  ) &
               /   (Coord_I(iCoord) - Coord_I(iCoord+1))
       end if
    end if

  end subroutine find_cell4
  !============================================================================
  subroutine find_cell8(MinCoord, MaxCoord, Coord, iCoord, dCoord, &
       Coord_I, DoExtrapolate, StringError, IsInside)

    ! double precision version of the subroutine above

    integer,          intent(in)           :: MinCoord, MaxCoord
    real(Real8_),     intent(in)           :: Coord
    integer,          intent(out)          :: iCoord
    real,             intent(out), optional:: dCoord
    real(Real8_),     intent(in),  optional:: Coord_I(MinCoord:)
    logical,          intent(in),  optional:: DoExtrapolate
    character(len=*), intent(in),  optional:: StringError
    logical,          intent(out), optional:: IsInside

    integer:: i, Di

    logical:: IsUniform

    character(len=*), parameter:: NameSub = 'find_cell8'
    !--------------------------------------------------------------------------

    if(present(Coord_I)) then
       IsUniform = size(Coord_I) < 2
    else
       IsUniform = .true.
    endif

    if(IsUniform)then
       ! Uniform grid case with normalized coordinate

       iCoord = min(MaxCoord-1, max(MinCoord, floor(Coord)))
       dCoord = Coord - iCoord

       if(Coord < MinCoord)then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) NameSub, ': ', StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord=', Coord
             call CON_stop(NameSub//': normalized coordinate is too small!')
          elseif(.not.DoExtrapolate)then
             ! Use lefttmost cell (first order accurate)
             dCoord = 0.0
          end if
          if(present(IsInside)) IsInside = .false.
       elseif(Coord > MaxCoord)then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord=', Coord
             call CON_stop(NameSub//': normalized coordinate is too large!')
          elseif(.not.DoExtrapolate)then
             ! Use rightmost cell (first order accurate)
             dCoord = 1.0
          endif
          if(present(IsInside)) IsInside = .false.
       else
          if(present(IsInside)) IsInside = .true.
       end if

    elseif(Coord_I(MinCoord) < Coord_I(MaxCoord))then

       ! Monotone increasing coordinates

       if(Coord < Coord_I(MinCoord))then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord, CoordMin, CoordMax=', &
                  Coord, Coord_I(MinCoord), Coord_I(MaxCoord)
             call CON_stop(NameSub//': coordinate is too small!')
          elseif(DoExtrapolate)then
             iCoord = MinCoord
             dCoord = (Coord - Coord_I(iCoord)) &
                  /   (Coord_I(iCoord+1) - Coord_I(iCoord))
          else
             iCoord = MinCoord
             dCoord = 0.0
          end if
          if(present(IsInside)) IsInside = .false.
          RETURN
       end if

       if(Coord > Coord_I(MaxCoord))then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord, CoordMin, CoordMax=', &
                  Coord, Coord_I(MinCoord), Coord_I(MaxCoord)
             call CON_stop(NameSub//': coordinate is too large!')
          elseif(DoExtrapolate)then
             iCoord = MaxCoord - 1
             dCoord = (Coord - Coord_I(iCoord))  &
                  /   (Coord_I(iCoord+1) - Coord_I(iCoord))
          else
             iCoord = MaxCoord - 1
             dCoord = 1.0
          end if
          if(present(IsInside)) IsInside = .false.
          RETURN
       end if

       if(present(IsInside)) IsInside = .true.

       ! binary search
       i  = (MinCoord + MaxCoord)/2
       Di = (MaxCoord - MinCoord)/2
       do
          Di = (Di + 1)/2
          if(Coord < Coord_I(i)) then
             i = max(MinCoord, i - Di)
          elseif(Coord > Coord_I(i+1))then
             i = min(MaxCoord-1, i + Di)
          else
             EXIT
          end if
       end do
       iCoord = i
       if(Coord_I(iCoord+1) == Coord_I(iCoord))then
          dCoord = 0.0
       else
          dCoord = (Coord             - Coord_I(iCoord)) &
               /   (Coord_I(iCoord+1) - Coord_I(iCoord))
       end if
    else

       ! Monotone decreasing coordinates

       if(Coord < Coord_I(MaxCoord))then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord, CoordMin, CoordMax=', &
                  Coord, Coord_I(MaxCoord), Coord_I(MinCoord)
             call CON_stop(NameSub//': coordinate is too small!')
          elseif(DoExtrapolate)then
             iCoord = MaxCoord - 1
             dCoord = (Coord_I(iCoord) - Coord) &
                  /   (Coord_I(iCoord) - Coord_I(iCoord+1))
          else
             iCoord = MaxCoord - 1
             dCoord = 1.0
          end if
          if(present(IsInside)) IsInside = .false.
          RETURN
       end if

       if(Coord > Coord_I(MinCoord))then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord, CoordMin, CoordMax=', &
                  Coord, Coord_I(MaxCoord), Coord_I(MinCoord)
             call CON_stop(NameSub//': coordinate is too large!')
          elseif(DoExtrapolate)then
             iCoord = MinCoord
             dCoord = (Coord_I(iCoord) - Coord)  &
                  /   (Coord_I(iCoord) - Coord_I(iCoord+1))
          else
             iCoord = MinCoord
             dCoord = 0.0
          end if
          if(present(IsInside)) IsInside = .false.
          RETURN
       end if

       if(present(IsInside)) IsInside = .true.

       ! binary search
       i  = (MinCoord + MaxCoord)/2
       Di = (MaxCoord - MinCoord)/2
       do
          Di = (Di + 1)/2
          if(Coord > Coord_I(i)) then
             i = max(MinCoord, i - Di)
          elseif(Coord < Coord_I(i+1))then
             i = min(MaxCoord-1, i + Di)
          else
             EXIT
          end if
       end do
       iCoord = i
       if(Coord_I(iCoord+1) == Coord_I(iCoord))then
          dCoord = 0.0
       else
          dCoord = (Coord_I(iCoord) - Coord  ) &
               /   (Coord_I(iCoord) - Coord_I(iCoord+1))
       end if
    end if

  end subroutine find_cell8
  !============================================================================

end module ModUtilities
!==============================================================================
