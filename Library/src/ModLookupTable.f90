!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModLookupTable

  ! Lookup tables define a function on a discrete grid so that the
  ! function can be interpolated within the range of the table.
  ! The function can have one or more arguments.
  ! The function value can be scalar or an array of reals.
  ! The lookup table is stored in the standard IDL format file, so it
  ! can be easily visualized.
  !
  ! Use lookup tables to calculate lookup properties, such as
  ! equation of state, opacities, ionization level etc.
  ! For example interpolate pressure as a function of the logarithm of
  ! density and logarithm of internal energy. All variables are in SI units.

  use ModReadParam,   ONLY: read_var
  use ModPlotFile,    ONLY: read_plot_file, save_plot_file
  use ModInterpolate, ONLY: interpolate_vector, interpolate_vector4, find_cell
  use ModUtilities,   ONLY: split_string, lower_case, CON_stop
  use ModIoUnit,      ONLY: UnitTmp_
  use ModKind,        ONLY: Real4_
  use ModMpi

  implicit none
  SAVE

  private ! except

  public:: init_lookup_table        ! set parameters of  the lookup table(s)
  public:: read_lookup_table_param  ! read parameters of the lookup table(s)
  public:: i_lookup_table           ! function returning the index of table
  public:: make_lookup_table        ! create/save 2D table from calculations
  public:: make_lookup_table_1d     ! create/save 1d table from calculations
  public:: make_lookup_table_3d     ! create/save 3D table from calculations
  public:: interpolate_lookup_table ! interpolate from lookup table
  public:: get_lookup_table         ! get information from a lookup table
  public:: test_lookup_table        ! unit test
  public:: copy_lookup_table_to_gpu ! copy Table_I to GPU

  integer, public, parameter:: MaxTable = 40 ! maximum number of tables
  integer, public :: nTable = 0     ! actual number of tables
  !$acc declare create(nTable)

  public:: TableType
  type TableType
     character(len=100):: NameTable        ! unique name for identification
     character(len=4)  :: NameCommand      ! command: load, make, save
     character(len=100):: NameFile         ! file name containing the table
     character(len=10) :: TypeFile         ! file type (ascii, real4, real8)
     character(len=500):: StringDescription! description of table
     character(len=500):: NameVar          ! name of indexes and values
     integer:: nIndex                      ! number of function arguments
     integer:: nValue                      ! number of values in each element
     integer, allocatable:: nIndex_I(:)    ! number of elements per dim
     real   , allocatable:: IndexMin_I(:)  ! minimum values for indexes
     real   , allocatable:: IndexMax_I(:)  ! maximum values for indexes
     real   , allocatable:: dIndex_I(:)    ! increment of indexes
     logical, allocatable:: IsLogIndex_I(:)! true if arguments are logarithmic
     real,         allocatable :: Value_VC(:,:,:,:,:,:) ! array of values
     real(Real4_), allocatable:: Value4_VC(:,:,:,:,:,:) ! real4 values
     integer:: nParam                      ! number of extra parameters
     real, allocatable :: Param_I(:)       ! parameter values
     logical, allocatable :: IsUniform_I(:)! false if index is non-uniform
     real, allocatable :: Index1_I(:)      ! non-uniform index 1
     real, allocatable :: Index2_I(:)      ! non-uniform index 2
     real, allocatable :: Index3_I(:)      ! non-uniform index 3
     real, allocatable :: Index4_I(:)      ! non-uniform index 4
     real, allocatable :: Index5_I(:)      ! non-uniform index 5
     real :: Time                          ! time, if applicable
  end type TableType

  ! The array of tables
  type(TableType), public, target :: Table_I(MaxTable)
  !$acc declare create(Table_I)

  ! private variables

  interface interpolate_lookup_table
     module procedure &
          interpolate_arg1, interpolate_arg1_scalar, &
          interpolate_arg2, interpolate_arg2_scalar, &
          interpolate_arg3, interpolate_arg3_scalar, &
          interpolate_arg4, interpolate_arg4_scalar, &
          interpolate_arg5, interpolate_arg5_scalar, &
          interpolate_arg_array
     module procedure interpolate_with_known_val  ! Table value is given
  end interface interpolate_lookup_table

  ! Array for variable names
  integer, parameter:: MaxVar = 200
  character(len=20):: NameVar_I(MaxVar)

contains
  !============================================================================

  subroutine deallocate_table_index(Ptr)

    type(TableType), pointer:: Ptr
    !--------------------------------------------------------------------------
    if(allocated(Ptr%nIndex_I))     deallocate(Ptr%nIndex_I)
    if(allocated(Ptr%IndexMin_I))   deallocate(Ptr%IndexMin_I)
    if(allocated(Ptr%IndexMax_I))   deallocate(Ptr%IndexMax_I)
    if(allocated(Ptr%IsLogIndex_I)) deallocate(Ptr%IsLogIndex_I)
    if(allocated(Ptr%IsUniform_I))  deallocate(Ptr%IsUniform_I)
    if(allocated(Ptr%dIndex_I))     deallocate(Ptr%dIndex_I)
    if(allocated(Ptr%Index1_I))     deallocate(Ptr%Index1_I)
    if(allocated(Ptr%Index2_I))     deallocate(Ptr%Index2_I)
    if(allocated(Ptr%Index3_I))     deallocate(Ptr%Index3_I)
    if(allocated(Ptr%Index4_I))     deallocate(Ptr%Index4_I)
    if(allocated(Ptr%Index5_I))     deallocate(Ptr%Index5_I)

  end subroutine deallocate_table_index
  !============================================================================
  subroutine allocate_table_index(Ptr)

    type(TableType), pointer:: Ptr
    !--------------------------------------------------------------------------

    call deallocate_table_index(Ptr)
    allocate( &
         Ptr%nIndex_I(Ptr%nIndex),     &
         Ptr%IndexMin_I(Ptr%nIndex),   &
         Ptr%IndexMax_I(Ptr%nIndex),   &
         Ptr%IsLogIndex_I(Ptr%nIndex), &
         Ptr%IsUniform_I(Ptr%nIndex),  &
         Ptr%dIndex_I(Ptr%nIndex)      )

  end subroutine allocate_table_index
  !============================================================================
  subroutine init_lookup_table(NameTable, NameCommand, NameVar, &
       nIndex_I, IndexMin_I, IndexMax_I, &
       NameFile, TypeFile, StringDescription, nParam, Param_I, &
       Index1_I, Index2_I, Index3_I, Index4_I, Index5_I, &
       Value1d_VC, Value2d_VC, Value3d_VC,Time)

    ! Initialize a lookup table
    !
    ! NameTable is a unique string identifier for the table.
    !
    ! NameCommand can be "load", "save", "make", "use".
    ! For "load" the table is loaded from file NameFile.
    ! For "make" the table is produced in memory, so its size nIndex_I
    ! and range of indexes from IndexMin_I to IndexMax_I must be given.
    ! The "save" command makes the table and then saves it into NameFile.
    ! The "use" command will load the table if the file already exists,
    ! otherwise it makes and saves the table.
    !
    ! TypeFile defines the file type to "ascii" (text file)
    ! "real4" (single precision binary) or "real8" (double precision binary).
    !
    ! StringDescription should describe the content of the table and
    ! it is saved into the first line of the file.
    !
    ! The file can contain nParam reals with values Param_I

    character(len=*), intent(in):: NameTable, NameCommand
    character(len=*), optional, intent(in):: NameVar
    integer,          optional, intent(in):: nIndex_I(:)
    real,             optional, intent(in):: IndexMin_I(:), IndexMax_I(:)
    character(len=*), optional, intent(in):: &
         NameFile, TypeFile, StringDescription
    integer,          optional, intent(in):: nParam
    real,             optional, intent(in):: Param_I(:)
    real,             optional, intent(in):: Index1_I(:), Index2_I(:)
    real,             optional, intent(in):: Index3_I(:), Index4_I(:)
    real,             optional, intent(in):: Index5_I(:)
    real,             optional, intent(in):: Value1d_VC(:,:)
    real,             optional, intent(in):: Value2d_VC(:,:,:)
    real,             optional, intent(in):: Value3d_VC(:,:,:,:)
    real,             optional, intent(in):: Time

    integer :: iTable, iError
    type(TableType), pointer:: Ptr

    character(len=*), parameter:: NameSub = 'init_lookup_table'
    !--------------------------------------------------------------------------

    ! Check if the table has been set already (say in a previous session)
    iTable = i_lookup_table(NameTable)
    if(iTable < 0)then
       ! new table
       nTable = nTable + 1

       if(nTable > MaxTable)then
          write(*,*)NameSub,' MaxTable =',MaxTable
          call CON_stop(NameSub//': number of tables exceeded MaxTable')
       end if

       iTable = nTable
    end if

    ! For sake of more concise source code, use a pointer to the table
    Ptr => Table_I(iTable)
    Ptr%NameTable = NameTable

    Ptr%TypeFile = 'real4'

    Ptr%NameCommand = NameCommand
    call lower_case(Ptr%NameCommand)

    select case(Ptr%NameCommand)
    case("load", "save", "use")
       Ptr%NameFile = NameFile
       Ptr%TypeFile = TypeFile
    case("make")
       Ptr%TypeFile = TypeFile
    case default
       call CON_stop(NameSub//': unknown command='//Ptr%NameCommand)
    end select

    if(ptr%NameCommand == "use")then
       ! Check if the file is already there or not
       open(UnitTmp_, FILE=NameFile, STATUS="old", IOSTAT=iError)
       if(iError == 0)then
          close(UnitTmp_)
          ptr%NameCommand = "load"
       else
          ptr%NameCommand = "save"
       end if
    end if

    if(Ptr%NameCommand == "load" )then
       call load_lookup_table(iTable)
       RETURN
    end if

    ! Save size of table and index ranges
    Ptr%nIndex = size(nIndex_I)
    call allocate_table_index(Ptr)

    Ptr%nIndex_I   = nIndex_I
    Ptr%IndexMin_I = IndexMin_I
    Ptr%IndexMax_I = IndexMax_I

    if(present(Param_I))then
       Ptr%nParam = size(Param_I)
       allocate(Ptr%Param_I(Ptr%nParam))
       Ptr%Param_I = Param_I
    elseif(present(nParam))then
       Ptr%nParam = nParam
       if(nParam > 0) allocate(Ptr%Param_I(nParam))
    else
       Ptr%nParam = 0
    end if

    if(present(StringDescription))then
       Ptr%StringDescription = StringDescription
    else
       Ptr%StringDescription = NameTable
    endif

    Ptr%NameVar = NameVar

    call split_string(Ptr%NameVar, MaxVar, NameVar_I, Ptr%nValue, &
         UseArraySyntaxIn=.true.)

    ! Do not count the names of the indexes and parameters
    Ptr%nValue = Ptr%nValue - Ptr%nIndex - Ptr%nParam

    ! Figure out which index is logarithmic
    Ptr%IsLogIndex_I = index(NameVar_I(1:Ptr%nIndex), "log") == 1

    ! Take logarithm of the ranges if logarithmic
    where(Ptr%IsLogIndex_I)
       Ptr%IndexMin_I = log10(Ptr%IndexMin_I)
       Ptr%IndexMax_I = log10(Ptr%IndexMax_I)
    end where

    ! Assign values if provided by optional arguments.
    ! Replaces the need for simplistic make_lookup_table calls.
    !  -> used for imf files in 1d
    if(present(Value1d_VC) .or. present(Value2d_VC) .or. present(Value3d_VC))then
       if(allocated(Ptr%Value4_VC)) deallocate(Ptr%Value4_VC)
       if(allocated(Ptr%Value_VC))  deallocate(Ptr%Value_VC)
    end if
    if(ptr%TypeFile == 'real4')then
       if(present(Value1d_VC))then
          allocate( &
               Ptr%Value4_VC(Ptr%nValue, Ptr%nIndex_I(1),1,1,1,1))
          Ptr%Value4_VC(:,:,1,1,1,1) = Value1d_VC
       elseif(present(Value2d_VC))then
          allocate( &
               Ptr%Value4_VC(Ptr%nValue, Ptr%nIndex_I(1), &
               Ptr%nIndex_I(2), 1, 1, 1))
          Ptr%Value4_VC(:,:,:,1,1,1) = Value2d_VC
       elseif(present(Value3d_VC))then
          allocate( &
               Ptr%Value4_VC(Ptr%nValue, Ptr%nIndex_I(1), &
               Ptr%nIndex_I(2), Ptr%nIndex_I(3), 1, 1))
          Ptr%Value4_VC(:,:,:,:,1,1) = Value3d_VC
       endif
    else
       if(present(Value1d_VC))then
          allocate( &
               Ptr%Value_VC(Ptr%nValue, Ptr%nIndex_I(1),1,1,1,1))
          Ptr%Value_VC(:,:,1,1,1,1) = Value1d_VC
       elseif(present(Value2d_VC))then
          allocate( &
               Ptr%Value_VC(Ptr%nValue, Ptr%nIndex_I(1), &
               Ptr%nIndex_I(2), 1, 1, 1))
          Ptr%Value_VC(:,:,:,1,1,1) = Value2d_VC
       elseif(present(Value3d_VC))then
          allocate( &
               Ptr%Value_VC(Ptr%nValue, Ptr%nIndex_I(1), &
               Ptr%nIndex_I(2), Ptr%nIndex_I(3), 1, 1))
          Ptr%Value_VC(:,:,:,:,1,1) = Value3d_VC
       endif
    end if

    ! Initialize uniform indices.
    Ptr%IsUniform_I = .true.

    ! Calculate increments
    Ptr%dIndex_I = (Ptr%IndexMax_I - Ptr%IndexMin_I)/(Ptr%nIndex_I - 1)

    ! Set indices to non-uniform if provided.
    if(present(Index1_I))then
       allocate(Ptr%Index1_I(Ptr%nIndex_I(1)))
       Ptr%Index1_I = Index1_I
       call check_index(iTable, Ptr%Index1_I, 1)
    endif
    if(present(Index2_I))then
       allocate(Ptr%Index2_I(Ptr%nIndex_I(2)))
       Ptr%Index2_I = Index2_I
       call check_index(iTable, Ptr%Index2_I, 2)
    endif
    if(present(Index3_I))then
       allocate(Ptr%Index3_I(Ptr%nIndex_I(3)))
       Ptr%Index3_I = Index3_I
       call check_index(iTable, Ptr%Index3_I, 3)
    endif
    if(present(Index4_I))then
       allocate(Ptr%Index4_I(Ptr%nIndex_I(4)))
       Ptr%Index4_I = Index4_I
       call check_index(iTable, Ptr%Index4_I, 4)
    endif
    if(present(Index5_I))then
       allocate(Ptr%Index5_I(Ptr%nIndex_I(5)))
       Ptr%Index5_I = Index5_I
       call check_index(iTable, Ptr%Index5_I, 5)
    endif
    if(present(Time))then
       Ptr%Time = Time
    else
       Ptr%Time = 0.0
    end if
  end subroutine init_lookup_table
  !============================================================================
  subroutine read_lookup_table_param

    ! Read parameters for one or more tables from an input parameter file.
    ! The table is identified by a name string which may contain a list
    ! inside curly brackets, e.g. {Xe Be Pl}_eos is expanded to three
    ! table names: "Xe_eos", "Be_eos", "Pl_eos".
    !
    ! The table maybe loaded

    integer, parameter:: MaxTableName = 20

    character(len=100):: NameCommand, NameTable, NameFile, TypeFile, &
         NameTable_I(MaxTableName), NameFile_I(MaxTableName)
    integer :: iTable, iIndex, nTableName, iTableName, nFileName
    type(TableType), pointer:: Ptr

    logical           :: DoReadTableParam
    integer           :: nTableParam, iTableParam, iError
    real, allocatable :: TableParam_I(:)
    character(len=100):: NameParam

    character(len=*), parameter:: NameSub = 'read_lookup_table_param'
    !--------------------------------------------------------------------------
    call read_var('NameTable', NameTable)

    call read_var('NameCommand', NameCommand, IsLowerCase = .true.)

    DoReadTableParam = index(NameCommand, "para") > 0
    NameCommand = NameCommand(1:4)
    ! If table is loaded, the table parameters are loaded too
    if(NameCommand == "load") DoReadTableParam = .false.

    ! Expand name
    call check_braces(NameTable, NameTable_I, nTableName)

    if(nTableName /= 1 .and. &
         NameCommand(1:4) /= "load" .and. NameCommand /= 'para') &
         call CON_stop( NameSub// &
         ': multiple table names can be used only for load and param')

    do iTableName = 1, nTableName

       NameTable = NameTable_I(iTableName)

       ! Check if the table has been set already (say in a previous session)

       iTable = i_lookup_table(NameTable)
       if(iTable < 0)then
          ! new table
          nTable = nTable + 1

          if(nTable > MaxTable)then
             write(*,*)NameSub,' MaxTable =',MaxTable
             call CON_stop(NameSub//': number of tables exceeded MaxTable')
          end if

          iTable = nTable
       end if

       ! For sake of more concise source code, use a pointer to the table
       Ptr => Table_I(iTable)
       Ptr%NameTable = NameTable
       Ptr%TypeFile = 'real4'
       if(iTableName == 1)then

          ! Read the parameters for this table
          select case(NameCommand)
          case("load", "save", "use")
             call read_var('NameFile', NameFile)
             call check_braces(NameFile, NameFile_I, nFileName)
             if(nTableName /= nFileName)call CON_stop(NameSub // &
                  ': the number of table names is not equal'// &
                  ' to the number of files names')
             call read_var('TypeFile', TypeFile)
          case("make", "para")
             call read_var('TypeFile', TypeFile)
          case default
             call CON_stop(NameSub//': unknown command='//Ptr%NameCommand)
          end select

          if(NameCommand == "use")then
             ! Check if table is already there or not
             open(UnitTmp_, FILE=NameFile_I(iTableName), STATUS="old", &
                  IOSTAT=iError)
             if(iError == 0)then
                close(UnitTmp_)
                NameCommand = "load"
                DoReadTableParam = .false.
             else
                NameCommand = "save"
             end if
          end if

          ! The table parameters have to be the same for all tables
          if(DoReadTableParam)then
             call read_var('NameParam', NameParam)
             call split_string(NameParam, MaxVar, NameVar_I, nTableParam, &
                  UseArraySyntaxIn=.true.)
             allocate(TableParam_I(nTableParam))
             do iTableParam = 1, nTableParam
                call read_var('TableParam', TableParam_I(iTableParam))
             end do
          end if

       end if

       if(NameCommand /= "para")then
          Ptr%NameCommand = NameCommand
          Ptr%TypeFile = TypeFile
       end if

       if(NameCommand == "load" .or. NameCommand == "save") &
            Ptr%NameFile = NameFile_I(iTableName)

       ! Load table on all processors
       if(NameCommand == "load") call load_lookup_table(iTable)

       ! Set table parameters (also allow overwriting loaded parameters)
       if(DoReadTableParam)then
          Ptr%nParam = nTableParam
          if(allocated(Ptr%Param_I)) deallocate(Ptr%Param_I)
          allocate(Ptr%Param_I(nTableParam))
          Ptr%Param_I = TableParam_I
          if(NameCommand == "load" .or. NameCommand == "para") &
               Ptr%NameVar = trim(Ptr%NameVar) // ' ' // trim(NameParam)
       end if

    end do

    if(DoReadTableParam) deallocate(TableParam_I)

    if(NameCommand == "load" .or. NameCommand == "para") RETURN

    call read_var('StringDescription', Ptr%StringDescription)
    call read_var('NameVar',           Ptr%NameVar)
    call read_var('nIndex',            Ptr%nIndex)
    call allocate_table_index(Ptr)

    call split_string(Ptr%NameVar, MaxVar, NameVar_I, Ptr%nValue, &
         UseArraySyntaxIn=.true.)

    ! Do not count the names of the indexes
    Ptr%nValue = Ptr%nValue - Ptr%nIndex

    ! Append parameter names if needed
    if(DoReadTableParam) &
         Ptr%NameVar = trim(Ptr%NameVar) // ' ' // trim(NameParam)

    ! Figure out which index is logarithmic
    Ptr%IsLogIndex_I = index(NameVar_I(1:Ptr%nIndex), "log") == 1

    do iIndex = 1, Ptr%nIndex
       call read_var('nIndex_I',   Ptr%nIndex_I(iIndex))
       call read_var('IndexMin',   Ptr%IndexMin_I(iIndex))
       call read_var('IndexMax',   Ptr%IndexMax_I(iIndex))

       ! Take logarithm of the ranges if logarithmic
       if(Ptr%IsLogIndex_I(iIndex)) then
          Ptr%IndexMin_I(iIndex) = log10(Ptr%IndexMin_I(iIndex))
          Ptr%IndexMax_I(iIndex) = log10(Ptr%IndexMax_I(iIndex))
       end if
    end do

    ! Calculate increments
    Ptr%dIndex_I = (Ptr%IndexMax_I - Ptr%IndexMin_I)/(Ptr % nIndex_I - 1)

    ! Tables defined by "make" are always uniform
    Ptr%IsUniform_I = .true.

  contains
    !==========================================================================

    subroutine check_braces(Name, Name_I, nString)

      character(LEN=*), intent(in) :: Name
      character(LEN=*), intent(out):: Name_I(MaxTableName)
      integer, intent(out) :: nString

      integer:: iBracePosition1, iBracePosition2, iString
      !------------------------------------------------------------------------

      iBracePosition1 = index(Name,'{')
      iBracePosition2 = index(Name,'}')

      if(iBracePosition1 < 1 .or. iBracePosition2 < 1 .or. &
           iBracePosition2 < iBracePosition1 )then
         nString = 1
         Name_I(1) = Name
         RETURN
      end if

      call split_string(Name(iBracePosition1+1:iBracePosition2-1),&
           MaxTableName, Name_I, nString)

      if(iBracePosition1 > 1)then
         do iString = 1, nString
            Name_I(iString) = Name(1:iBracePosition1-1)//trim(Name_I(iString))
         end do
      end if

      if(iBracePosition2 < len_trim(Name))then
         do iString = 1, nString
            Name_I(iString) =trim(Name_I(iString))// &
                 Name(iBracePosition2 + 1: len_trim(Name))
         end do
      end if

    end subroutine check_braces
    !==========================================================================

  end subroutine read_lookup_table_param
  !============================================================================

  integer function i_lookup_table(Name)

    ! return the index of the lookup table based on its name
    ! return -1 if the table was not found

    character(len=*), intent(in):: Name
    integer :: iTable

    !--------------------------------------------------------------------------
    do iTable = 1, nTable
       if(Table_I(iTable)%NameTable == Name) then
          i_lookup_table = iTable
          RETURN
       end if
    end do
    i_lookup_table = -1

  end function i_lookup_table
  !============================================================================

  subroutine load_lookup_table(iTable)

    integer, intent(in) :: iTable

    type(TableType), pointer:: Ptr
    integer :: nVar, iIndex

    ! since number of parameters is not known in advance, this array is needed
    real, allocatable:: TableParam_I(:)
    integer, allocatable:: nIndex_I(:)

    character(len=*), parameter:: NameSub = 'load_lookup_table'
    !--------------------------------------------------------------------------

    if(iTable > nTable) call CON_stop(NameSub//' iTable larger than nTable')

    Ptr => Table_I(iTable)

    allocate(TableParam_I(1000), nIndex_I(5))

    call read_plot_file( Ptr%NameFile,            &
         TypeFileIn      = Ptr%TypeFile,          &
         StringHeaderOut = Ptr%StringDescription, &
         nDimOut         = Ptr%nIndex,            &
         nOut_D          = nIndex_I,              &
         nVarOut         = Ptr%nValue,            &
         NameVarOut      = Ptr%NameVar,           &
         nParamOut       = Ptr%nParam,            &
         TimeOut         = Ptr%Time,              &
         ParamOut_I      = TableParam_I)

    call allocate_table_index(Ptr)

    Ptr%nIndex_I = nIndex_I(1:Ptr%nIndex)

    if(Ptr%nParam > 0)then
       if(allocated(Ptr%Param_I)) deallocate(Ptr%Param_I)
       allocate(Ptr%Param_I(Ptr%nParam))
       Ptr%Param_I = TableParam_I(1:Ptr%nParam)
    end if

    ! Figure out which index is logarithmic
    call split_string(Ptr%NameVar, MaxVar, NameVar_I, nVar)
    Ptr%IsLogIndex_I = index(NameVar_I(1:Ptr%nIndex), "log") == 1

    allocate( &
         Ptr%Index1_I(nIndex_I(1)), &
         Ptr%Index2_I(nIndex_I(2)), &
         Ptr%Index3_I(nIndex_I(3)), &
         Ptr%Index4_I(nIndex_I(4)), &
         Ptr%Index5_I(nIndex_I(5)) )

    if(allocated(Ptr%Value_VC))  deallocate(Ptr%Value_VC)
    if(allocated(Ptr%Value4_VC)) deallocate(Ptr%Value4_VC)
    if(Ptr%TypeFile == 'real4')then
       allocate(Ptr%Value4_VC(Ptr%nValue, &
            nIndex_I(1),nIndex_I(2),nIndex_I(3),nIndex_I(4),nIndex_I(5)))
       call read_plot_file( Ptr%NameFile,   &
            TypeFileIn    = Ptr%TypeFile,   &
            CoordMinOut_D = Ptr%IndexMin_I, &
            CoordMaxOut_D = Ptr%IndexMax_I, &
            Var4Out_VI5   = Ptr%Value4_VC,  &
            Coord1Out_I   = Ptr%Index1_I,   &
            Coord2Out_I   = Ptr%Index2_I,   &
            Coord3Out_I   = Ptr%Index3_I,   &
            Coord4Out_I   = Ptr%Index4_I,   &
            Coord5Out_I   = Ptr%Index5_I    )
    else
       allocate(Ptr%Value_VC(Ptr%nValue, &
            nIndex_I(1),nIndex_I(2),nIndex_I(3),nIndex_I(4),nIndex_I(5)))
       call read_plot_file( Ptr%NameFile,   &
            TypeFileIn    = Ptr%TypeFile,   &
            CoordMinOut_D = Ptr%IndexMin_I, &
            CoordMaxOut_D = Ptr%IndexMax_I, &
            VarOut_VI5    = Ptr%Value_VC,   &
            Coord1Out_I   = Ptr%Index1_I,   &
            Coord2Out_I   = Ptr%Index2_I,   &
            Coord3Out_I   = Ptr%Index3_I,   &
            Coord4Out_I   = Ptr%Index4_I,   &
            Coord5Out_I   = Ptr%Index5_I    )
    end if

    ! Calculate increments assuming uniform grid
    Ptr%dIndex_I = (Ptr%IndexMax_I - Ptr%IndexMin_I)/(Ptr%nIndex_I - 1)

    ! Check monotonicity and uniformity of 1..Ptr%nIndex index arrays
    do iIndex = 1, Ptr%nIndex
       select case(iIndex)
       case(1)
          call check_index(iTable, Ptr%Index1_I, iIndex)
       case(2)
          call check_index(iTable, Ptr%Index2_I, iIndex)
       case(3)
          call check_index(iTable, Ptr%Index3_I, iIndex)
       case(4)
          call check_index(iTable, Ptr%Index4_I, iIndex)
       case(5)
          call check_index(iTable, Ptr%Index5_I, iIndex)
       end select
    end do

    deallocate(TableParam_I, nIndex_I)

  end subroutine load_lookup_table
  !============================================================================

  subroutine check_index(iTable, Index_I, iIndex)

    ! Check monotonicity and uniformity of index arrays
    ! Uniform indexes are reallocated to a one-element array

    integer, intent(in)              :: iTable, iIndex
    real, allocatable, intent(inout) :: Index_I(:)

    integer :: i, n
    type(TableType), pointer:: Ptr

    character(len=*), parameter:: NameSub = 'check_index'
    !--------------------------------------------------------------------------
    Ptr => Table_I(iTable)
    n = size(Index_I)
    if(n < 3) RETURN

    if(any(Index_I(2:n) < Index_I(1:n-1)))then
       do i = 2, n
          if(Index_I(i-1) < Index_I(i)) CYCLE
          write(*,*) '!!! ERROR at i, Index_I(i-1:i)=', i, Index_I(i-1:i)
       end do
       call CON_stop(NameSub// &
            ' ERROR: decreasing index in '//trim(Ptr%NameFile)//':', iIndex)
    end if
    ! Figure out which indices are non-uniform.
    if (maxval(Index_I(2:n) - Index_I(1:n-1)) > 1.3*Ptr%dIndex_I(iIndex))then
       Ptr%IsUniform_I(iIndex) = .false.
    else
       Ptr%IsUniform_I(iIndex) = .true.
       deallocate(Index_I)
       allocate(Index_I(1))
    end if

  end subroutine check_index
  !============================================================================

  subroutine make_lookup_table_1d(iTable, calc_table_var, iComm, UseRealIndex)

    ! Fill in table iTable using the subroutine calc_table_var
    ! The optional communicator allows for parallel execution

    integer, intent(in):: iTable  ! table index
    interface
       subroutine calc_table_var(iTable, Arg1, Value_V)
         integer, intent(in) :: iTable
         real,    intent(in) :: Arg1
         real,    intent(out):: Value_V(:)
       end subroutine calc_table_var
    end interface

    integer, optional, intent(in):: iComm
    logical, optional, intent(in):: UseRealIndex

    integer:: iProc, nProc, iError
    integer:: i1, n1, nValue
    logical:: IsLog1, IsUniform1
    real   :: Index1Min, dIndex1, Index1
    real, allocatable:: Value_VC(:,:)
    type(TableType), pointer:: Ptr

    character(len=*), parameter:: NameSub = 'make_lookup_table_1d'
    !--------------------------------------------------------------------------
    Ptr => Table_I(iTable)

    if(Ptr%NameCommand /= "make" .and. Ptr%NameCommand /= "save") &
         RETURN

    if(Ptr%nIndex /= 1)call CON_stop(NameSub//': table '&
         //trim(Ptr%NameTable)//' should be 1D')

    ! Use simple scalars for sake of legibility
    n1        = Ptr%nIndex_I(1)
    IsLog1    = Ptr%IsLogIndex_I(1)
    IsUniform1= Ptr%IsUniform_I(1)
    Index1Min = Ptr%IndexMin_I(1)
    dIndex1   = Ptr%dIndex_I(1)

    nValue    = Ptr%nValue

    ! Get processor index and total number of processors
    if(present(iComm))then
       call MPI_comm_rank(iComm,iProc,iError)
       call MPI_comm_size(iComm,nProc,iError)
    else
       iProc = 0
       nProc = 1
    end if

    if(iProc == 0)write(*,'(3a)')NameSub,' is creating table ', &
         trim(Ptr%NameTable)

    ! Allocate Value_VC array
    if(allocated(Ptr%Value4_VC)) deallocate(Ptr%Value4_VC)
    if(allocated(Ptr%Value_VC))  deallocate(Ptr%Value_VC)
    allocate(Value_VC(nValue,n1))
    Value_VC = 0.0

    ! Fill up lookup table in parallel
    do i1 = iProc+1, n1, nProc
       if(IsUniform1) then
          Index1 = Index1Min + (i1 - 1)*dIndex1
          if(IsLog1) Index1 = 10**Index1
       else
          if(UseRealIndex) then
             Index1 = Ptr%Index1_I(i1)
          else
             Index1 = i1
          endif
       endif
       call calc_table_var(iTable, Index1, Value_VC(:,i1))
    end do

    ! Collect and copy result into table
    if(nProc > 1)call MPI_allreduce(MPI_IN_PLACE, Value_VC, nValue*n1, &
         MPI_REAL, MPI_SUM, iComm, iError)

    if(Ptr%TypeFile == 'real4')then
       allocate(Ptr%Value4_VC(nValue,n1,1,1,1,1))
       Ptr%Value4_VC(:,:,1,1,1,1) = Value_VC
    else
       allocate(Ptr%Value_VC(nValue,n1,1,1,1,1))
       Ptr%Value_VC(:,:,1,1,1,1) = Value_VC
    end if

    deallocate(Value_VC)

    if(Ptr%NameCommand == "save" .and. iProc == 0)then
       ! Cannot pass unallocated arrays so check Ptr%Param_I
       if( allocated(Ptr%Param_I) ) then
          if(Ptr%TypeFile == 'real4')then
             call save_plot_file( &
                  Ptr%NameFile,                                 &
                  TypeFileIn     = Ptr%TypeFile,                &
                  StringHeaderIn = Ptr%StringDescription,       &
                  NameVarIn      = Ptr%NameVar,                 &
                  nDimIn         = Ptr%nIndex,                  &
                  CoordMinIn_D   = Ptr%IndexMin_I,              &
                  CoordMaxIn_D   = Ptr%IndexMax_I,              &
                  Var4In_VI      = Ptr%Value4_VC(:,:,1,1,1,1),  &
                  ParamIn_I      = Ptr%Param_I,                 &
                  TimeIn         = Ptr%Time)
          else
             call save_plot_file( &
                  Ptr%NameFile,                                 &
                  TypeFileIn     = Ptr%TypeFile,                &
                  StringHeaderIn = Ptr%StringDescription,       &
                  NameVarIn      = Ptr%NameVar,                 &
                  nDimIn         = Ptr%nIndex,                  &
                  CoordMinIn_D   = Ptr%IndexMin_I,              &
                  CoordMaxIn_D   = Ptr%IndexMax_I,              &
                  VarIn_VI      = Ptr%Value_VC(:,:,1,1,1,1),    &
                  ParamIn_I      = Ptr%Param_I,                 &
                  TimeIn         = Ptr%Time)
          end if
       else
          if(Ptr%TypeFile == 'real4')then
             call save_plot_file( &
                  Ptr%NameFile,                                 &
                  TypeFileIn     = Ptr%TypeFile,                &
                  StringHeaderIn = Ptr%StringDescription,       &
                  NameVarIn      = Ptr%NameVar,                 &
                  nDimIn         = Ptr%nIndex,                  &
                  CoordMinIn_D   = Ptr%IndexMin_I,              &
                  CoordMaxIn_D   = Ptr%IndexMax_I,              &
                  Var4In_VI      = Ptr%Value4_VC(:,:,1,1,1,1),  &
                  TimeIn         = Ptr%Time)
          else
             call save_plot_file( &
                  Ptr%NameFile,                                 &
                  TypeFileIn     = Ptr%TypeFile,                &
                  StringHeaderIn = Ptr%StringDescription,       &
                  NameVarIn      = Ptr%NameVar,                 &
                  nDimIn         = Ptr%nIndex,                  &
                  CoordMinIn_D   = Ptr%IndexMin_I,              &
                  CoordMaxIn_D   = Ptr%IndexMax_I,              &
                  VarIn_VI       = Ptr%Value_VC(:,:,1,1,1,1),   &
                  TimeIn         = Ptr%Time)
          end if
       end if
    end if

    ! Make sure that all processors are done
    if(nProc>1)call MPI_barrier(iComm, iError)

    ! Make sure it is not saved again
    Ptr%NameCommand = "done"

  end subroutine make_lookup_table_1d
  !============================================================================

  subroutine make_lookup_table(iTable, calc_table_var, iComm)

    ! Fill in table iTable using the subroutine calc_table_var
    ! The optional communicator allows for parallel execution

    integer, intent(in):: iTable  ! table index
    interface
       subroutine calc_table_var(iTable, Arg1, Arg2, Value_V)
         integer, intent(in) :: iTable
         real,    intent(in) :: Arg1, Arg2
         real,    intent(out):: Value_V(:)
       end subroutine calc_table_var
    end interface

    integer, optional, intent(in):: iComm

    integer:: iProc, nProc, iError
    integer:: i1, i2, n1, n2, nValue
    logical:: IsLog1, IsLog2
    real   :: Index1Min, Index2Min, dIndex1, dIndex2, Index1, Index2
    real, allocatable:: Value_VC(:,:,:)
    type(TableType), pointer:: Ptr

    character(len=*), parameter:: NameSub = 'make_lookup_table'
    !--------------------------------------------------------------------------

    Ptr => Table_I(iTable)

    if(Ptr%NameCommand /= "make" .and. Ptr%NameCommand /= "save") &
         RETURN

    if(Ptr%nIndex /= 2)call CON_stop(NameSub//': table '&
         //trim(Ptr%NameTable)//' should be 2D')

    ! Use simple scalars for sake of legibility
    n1        = Ptr%nIndex_I(1)
    IsLog1    = Ptr%IsLogIndex_I(1)
    Index1Min = Ptr%IndexMin_I(1)
    dIndex1   = Ptr%dIndex_I(1)

    n2        = Ptr%nIndex_I(2)
    IsLog2    = Ptr%IsLogIndex_I(2)
    Index2Min = Ptr%IndexMin_I(2)
    dIndex2   = Ptr%dIndex_I(2)

    nValue    = Ptr%nValue

    ! Get processor index and total number of processors
    if(present(iComm))then
       call MPI_comm_rank(iComm,iProc,iError)
       call MPI_comm_size(iComm,nProc,iError)
    else
       iProc = 0
       nProc = 1
    end if

    if(iProc == 0)write(*,'(3a)')NameSub,' is creating table ', &
         trim(Ptr%NameTable)

    ! Allocate Value_VC arrays
    if(allocated(Ptr%Value4_VC)) deallocate(Ptr%Value4_VC)
    if(allocated(Ptr%Value_VC))  deallocate(Ptr%Value_VC)
    allocate(Value_VC(nValue,n1,n2))
    Value_VC = 0.0

    ! Fill up lookup table in parallel
    do i2 = iProc+1, n2, nProc
       Index2 = Index2Min + (i2 - 1)*dIndex2
       if(IsLog2) Index2 = 10**Index2
       do i1 = 1, n1
          Index1 = Index1Min + (i1 - 1)*dIndex1
          if(IsLog1) Index1 = 10**Index1
          call calc_table_var(iTable, Index1, Index2, Value_VC(:,i1,i2))
       end do
    end do

    ! Collect and copy result into table
    if(nProc > 1) call MPI_allreduce(MPI_IN_PLACE, Value_VC, nValue*n1*n2, &
         MPI_REAL, MPI_SUM, iComm, iError)

    if(Ptr%TypeFile == 'real4')then
       allocate(Ptr%Value4_VC(nValue,n1,n2,1,1,1))
       Ptr%Value4_VC(:,:,:,1,1,1) = Value_VC
    else
       allocate(Ptr%Value_VC(nValue,n1,n2,1,1,1))
       Ptr%Value_VC(:,:,:,1,1,1) = Value_VC
    end if

    deallocate(Value_VC)

    if(Ptr%NameCommand == "save" .and. iProc == 0)then
       ! Cannot pass unallocated arrays so check Ptr%Param_I
       if( allocated(Ptr%Param_I) ) then
          if(Ptr%TypeFile == 'real4')then
             call save_plot_file( &
                  Ptr%NameFile,                                 &
                  TypeFileIn     = Ptr%TypeFile,                &
                  StringHeaderIn = Ptr%StringDescription,       &
                  NameVarIn      = Ptr%NameVar,                 &
                  nDimIn         = Ptr%nIndex,                  &
                  CoordMinIn_D   = Ptr%IndexMin_I,              &
                  CoordMaxIn_D   = Ptr%IndexMax_I,              &
                  Var4In_VII     = Ptr%Value4_VC(:,:,:,1,1,1),  &
                  ParamIn_I      = Ptr%Param_I,                 &
                  TimeIn         = Ptr%Time)
          else
             call save_plot_file( &
                  Ptr%NameFile,                                 &
                  TypeFileIn     = Ptr%TypeFile,                &
                  StringHeaderIn = Ptr%StringDescription,       &
                  NameVarIn      = Ptr%NameVar,                 &
                  nDimIn         = Ptr%nIndex,                  &
                  CoordMinIn_D   = Ptr%IndexMin_I,              &
                  CoordMaxIn_D   = Ptr%IndexMax_I,              &
                  VarIn_VII      = Ptr%Value_VC(:,:,:,1,1,1),   &
                  ParamIn_I      = Ptr%Param_I,                 &
                  TimeIn         = Ptr%Time)
          end if
       else
          if(Ptr%TypeFile == 'real4')then
             call save_plot_file( &
                  Ptr%NameFile,                                 &
                  TypeFileIn     = Ptr%TypeFile,                &
                  StringHeaderIn = Ptr%StringDescription,       &
                  NameVarIn      = Ptr%NameVar,                 &
                  nDimIn         = Ptr%nIndex,                  &
                  CoordMinIn_D   = Ptr%IndexMin_I,              &
                  CoordMaxIn_D   = Ptr%IndexMax_I,              &
                  Var4In_VII     = Ptr%Value4_VC(:,:,:,1,1,1),  &
                  TimeIn         = Ptr%Time)
          else
             call save_plot_file( &
                  Ptr%NameFile,                                 &
                  TypeFileIn     = Ptr%TypeFile,                &
                  StringHeaderIn = Ptr%StringDescription,       &
                  NameVarIn      = Ptr%NameVar,                 &
                  nDimIn         = Ptr%nIndex,                  &
                  CoordMinIn_D   = Ptr%IndexMin_I,              &
                  CoordMaxIn_D   = Ptr%IndexMax_I,              &
                  VarIn_VII      = Ptr%Value_VC(:,:,:,1,1,1),   &
                  TimeIn         = Ptr%Time)
          end if
       end if
    end if

    ! Make sure that all processors are done
    if(nProc > 1)call MPI_barrier(iComm, iError)

    ! Make sure it is not saved again
    Ptr%NameCommand = "done"

  end subroutine make_lookup_table
  !============================================================================

  subroutine make_lookup_table_3d(iTable, calc_table_var, iComm)

    ! Fill in table iTable using the subroutine calc_table_var
    ! The optional communicator allows for parallel execution

    integer, intent(in):: iTable  ! table index
    interface
       subroutine calc_table_var(iTable, Arg1, Arg2, Arg3, Value_V)
         integer, intent(in) :: iTable
         real,    intent(in) :: Arg1, Arg2, Arg3
         real,    intent(out):: Value_V(:)
       end subroutine calc_table_var
    end interface

    integer, optional, intent(in):: iComm

    integer:: iProc, nProc, iError
    integer:: i1, i2, i3, n1, n2, n3, nValue
    logical:: IsLog1, IsLog2, IsLog3
    real   :: Index1Min, Index2Min, Index3Min, dIndex1, dIndex2, dIndex3
    real   :: Index1, Index2, Index3
    real, allocatable:: Value_VC(:,:,:,:)
    type(TableType), pointer:: Ptr

    character(len=*), parameter:: NameSub = 'make_lookup_table_3d'
    !--------------------------------------------------------------------------

    Ptr => Table_I(iTable)

    if(Ptr%NameCommand /= "make" .and. Ptr%NameCommand /= "save") &
         RETURN

    if(Ptr%nIndex /= 3)call CON_stop(NameSub//': table '&
         //trim(Ptr%NameTable)//' should be 3D')

    ! Use simple scalars for sake of legibility
    n1        = Ptr%nIndex_I(1)
    IsLog1    = Ptr%IsLogIndex_I(1)
    Index1Min = Ptr%IndexMin_I(1)
    dIndex1   = Ptr%dIndex_I(1)

    n2        = Ptr%nIndex_I(2)
    IsLog2    = Ptr%IsLogIndex_I(2)
    Index2Min = Ptr%IndexMin_I(2)
    dIndex2   = Ptr%dIndex_I(2)

    n3        = Ptr%nIndex_I(3)
    IsLog3    = Ptr%IsLogIndex_I(3)
    Index3Min = Ptr%IndexMin_I(3)
    dIndex3   = Ptr%dIndex_I(3)

    nValue    = Ptr%nValue

    ! Get processor index and total number of processors
    if(present(iComm))then
       call MPI_comm_rank(iComm,iProc,iError)
       call MPI_comm_size(iComm,nProc,iError)
    else
       iProc = 0
       nProc = 1
    end if

    if(iProc == 0)write(*,'(3a)')NameSub,' is creating table ', &
         trim(Ptr%NameTable)

    ! Allocate Value_VC array
    if(allocated(Ptr%Value4_VC)) deallocate(Ptr%Value4_VC)
    if(allocated(Ptr%Value_VC))  deallocate(Ptr%Value_VC)
    allocate(Value_VC(nValue,n1,n2,n3))
    Value_VC = 0.0

    ! Fill up lookup table in parallel
    do i3 = 1, n3
       Index3 = Index3Min + (i3 - 1)*dIndex3
       if(IsLog3) Index3 = 10**Index3
       do i2 = 1, n2
          ! Parallel processing in the last 2 indexes
          if(modulo(i2-1 + (i3-1)*n2, nProc) /= iProc) CYCLE

          Index2 = Index2Min + (i2 - 1)*dIndex2
          if(IsLog2) Index2 = 10**Index2
          do i1 = 1, n1
             Index1 = Index1Min + (i1 - 1)*dIndex1
             if(IsLog1) Index1 = 10**Index1
             call calc_table_var(iTable, Index1, Index2, Index3, &
                  Value_VC(:,i1,i2,i3))
          end do
       end do
    end do

    ! Collect and copy result into table
    if(nProc > 1) call MPI_allreduce(MPI_IN_PLACE, Value_VC, nValue*n1*n2*n3, &
         MPI_REAL, MPI_SUM, iComm, iError)

    if(Ptr%TypeFile == 'real4')then
       allocate(Ptr%Value4_VC(nValue,n1,n2,n3,1,1))
       Ptr%Value4_VC(:,:,:,:,1,1) = Value_VC
    else
       allocate(Ptr%Value_VC(nValue,n1,n2,n3,1,1))
       Ptr%Value_VC(:,:,:,:,1,1) = Value_VC
    end if

    deallocate(Value_VC)

    if(Ptr%NameCommand == "save" .and. iProc == 0)then
       ! Cannot pass unallocated arrays so check Ptr%Param_I
       if( allocated(Ptr%Param_I) ) then
          if(Ptr%TypeFile == 'real4')then
             call save_plot_file( &
                  Ptr%NameFile,                                 &
                  TypeFileIn     = Ptr%TypeFile,                &
                  StringHeaderIn = Ptr%StringDescription,       &
                  NameVarIn      = Ptr%NameVar,                 &
                  nDimIn         = Ptr%nIndex,                  &
                  CoordMinIn_D   = Ptr%IndexMin_I,              &
                  CoordMaxIn_D   = Ptr%IndexMax_I,              &
                  Var4In_VIII    = Ptr%Value4_VC(:,:,:,:,1,1),  &
                  ParamIn_I      = Ptr%Param_I,                 &
                  TimeIn         = Ptr%Time)
          else
             call save_plot_file( &
                  Ptr%NameFile,                                 &
                  TypeFileIn     = Ptr%TypeFile,                &
                  StringHeaderIn = Ptr%StringDescription,       &
                  NameVarIn      = Ptr%NameVar,                 &
                  nDimIn         = Ptr%nIndex,                  &
                  CoordMinIn_D   = Ptr%IndexMin_I,              &
                  CoordMaxIn_D   = Ptr%IndexMax_I,              &
                  VarIn_VIII     = Ptr%Value_VC(:,:,:,:,1,1),   &
                  ParamIn_I      = Ptr%Param_I,                 &
                  TimeIn         = Ptr%Time)
          end if
       else
          if(Ptr%TypeFile == 'real4')then
             call save_plot_file( &
                  Ptr%NameFile,                                 &
                  TypeFileIn     = Ptr%TypeFile,                &
                  StringHeaderIn = Ptr%StringDescription,       &
                  NameVarIn      = Ptr%NameVar,                 &
                  nDimIn         = Ptr%nIndex,                  &
                  CoordMinIn_D   = Ptr%IndexMin_I,              &
                  CoordMaxIn_D   = Ptr%IndexMax_I,              &
                  Var4In_VIII    = Ptr%Value4_VC(:,:,:,:,1,1),  &
                  TimeIn         = Ptr%Time)
          else
             call save_plot_file( &
                  Ptr%NameFile,                                 &
                  TypeFileIn     = Ptr%TypeFile,                &
                  StringHeaderIn = Ptr%StringDescription,       &
                  NameVarIn      = Ptr%NameVar,                 &
                  nDimIn         = Ptr%nIndex,                  &
                  CoordMinIn_D   = Ptr%IndexMin_I,              &
                  CoordMaxIn_D   = Ptr%IndexMax_I,              &
                  VarIn_VIII     = Ptr%Value_VC(:,:,:,:,1,1),   &
                  TimeIn         = Ptr%Time)
          end if
       end if
    end if

    ! Make sure that all processors are done
    if(nProc > 1)call MPI_barrier(iComm, iError)

    ! Make sure it is not saved again
    Ptr%NameCommand = "done"

  end subroutine make_lookup_table_3d
  !============================================================================

  subroutine interpolate_arg1(iTable, Arg1, Value_V, DoExtrapolate)
    !$acc routine seq

    ! Return the array of values Value_V corresponding to argument Arg1.
    ! If DoExtrapolate is not present, stop with an error if the arguments
    ! are out of range. If it is present and false, return the value of
    ! the closest element in the table. If it is present and true, do a
    ! linear extrapolation.

    integer, intent(in) :: iTable            ! table
    real,    intent(in) :: Arg1              ! input arguments
    real,    intent(out):: Value_V(:)        ! output values
    logical, optional, intent(in):: DoExtrapolate ! optional extrapolation

    !--------------------------------------------------------------------------
    call interpolate_arg_array(iTable, [Arg1], Value_V, DoExtrapolate)

  end subroutine interpolate_arg1
  !============================================================================

  subroutine interpolate_arg2(iTable, Arg1, Arg2, Value_V, DoExtrapolate)
    !$acc routine seq

    ! Return the array of values Value_V corresponding to arguments
    ! Arg1 and Arg2.
    ! If DoExtrapolate is not present, stop with an error if the arguments
    ! are out of range. If it is present and false, return the value of
    ! the closest element in the table. If it is present and true, do a
    ! linear extrapolation.

    integer, intent(in) :: iTable            ! table
    real,    intent(in) :: Arg1, Arg2        ! input arguments
    real,    intent(out):: Value_V(:)        ! output values
    logical, optional, intent(in):: DoExtrapolate ! optional extrapolation

    !--------------------------------------------------------------------------
    call interpolate_arg_array(iTable, [Arg1, Arg2], Value_V, &
         DoExtrapolate)

  end subroutine interpolate_arg2
  !============================================================================

  subroutine interpolate_arg3(iTable, Arg1, Arg2, Arg3, Value_V, DoExtrapolate)
    !$acc routine seq

    ! Return the array of values Value_V corresponding to arguments
    ! Arg1, Arg2 and Arg3.
    ! If DoExtrapolate is not present, stop with an error if the arguments
    ! are out of range. If it is present and false, return the value of
    ! the closest element in the table. If it is present and true, do a
    ! linear extrapolation.

    integer, intent(in) :: iTable            ! table
    real,    intent(in) :: Arg1, Arg2, Arg3  ! input arguments
    real,    intent(out):: Value_V(:)        ! output values
    logical, optional, intent(in):: DoExtrapolate ! optional extrapolation

    !--------------------------------------------------------------------------
    call interpolate_arg_array(iTable, [Arg1, Arg2, Arg3], Value_V, &
         DoExtrapolate)

  end subroutine interpolate_arg3
  !============================================================================

  subroutine interpolate_arg4(iTable, Arg1, Arg2, Arg3, Arg4, Value_V, &
       DoExtrapolate)
    !$acc routine seq

    ! Return the array of values Value_V corresponding to arguments
    ! Arg1, Arg2, Arg3, and Arg4.
    ! If DoExtrapolate is not present, stop with an error if the arguments
    ! are out of range. If it is present and false, return the value of
    ! the closest element in the table. If it is present and true, do a
    ! linear extrapolation.

    integer, intent(in) :: iTable                 ! table
    real,    intent(in) :: Arg1, Arg2, Arg3, Arg4 ! input arguments
    real,    intent(out):: Value_V(:)             ! output values
    logical, optional, intent(in):: DoExtrapolate ! optional extrapolation

    !--------------------------------------------------------------------------
    call interpolate_arg_array(iTable, [Arg1, Arg2, Arg3, Arg4], Value_V, &
         DoExtrapolate)

  end subroutine interpolate_arg4
  !============================================================================

  subroutine interpolate_arg5(iTable, Arg1, Arg2, Arg3, Arg4, Arg5, Value_V, &
       DoExtrapolate)
    !$acc routine seq

    ! Return the array of values Value_V corresponding to arguments
    ! Arg1, Arg2, Arg3, Arg4, and Arg5.
    ! If DoExtrapolate is not present, stop with an error if the arguments
    ! are out of range. If it is present and false, return the value of
    ! the closest element in the table. If it is present and true, do a
    ! linear extrapolation.

    integer, intent(in) :: iTable                       ! table
    real,    intent(in) :: Arg1, Arg2, Arg3, Arg4, Arg5 ! input arguments
    real,    intent(out):: Value_V(:)                   ! output values
    logical, optional, intent(in):: DoExtrapolate ! optional extrapolation

    !--------------------------------------------------------------------------
    call interpolate_arg_array(iTable, [Arg1, Arg2, Arg3, Arg4, Arg5], &
         Value_V, DoExtrapolate)

  end subroutine interpolate_arg5
  !============================================================================

  subroutine interpolate_arg1_scalar(iTable, Arg1, Value, DoExtrapolate)
    !$acc routine seq

    ! Return the scalar Value corresponding to argument Arg1.
    ! If DoExtrapolate is not present, stop with an error if the arguments
    ! are out of range. If it is present and false, return the value of
    ! the closest element in the table. If it is present and true, do a
    ! linear extrapolation.

    integer, intent(in) :: iTable            ! table
    real,    intent(in) :: Arg1              ! input arguments
    real,    intent(out):: Value             ! output value
    logical, optional, intent(in):: DoExtrapolate ! optional extrapolation

    real:: Value_V(1)

    !--------------------------------------------------------------------------
    call interpolate_arg_array(iTable, [Arg1], Value_V, DoExtrapolate)
    Value = Value_V(1)

  end subroutine interpolate_arg1_scalar
  !============================================================================

  subroutine interpolate_arg2_scalar(iTable, Arg1, Arg2, Value, DoExtrapolate)
    !$acc routine seq

    ! Return the scalar Value corresponding to arguments Arg1 and Arg2.
    ! If DoExtrapolate is not present, stop with an error if the arguments
    ! are out of range. If it is present and false, return the value of
    ! the closest element in the table. If it is present and true, do a
    ! linear extrapolation.

    integer, intent(in) :: iTable            ! table
    real,    intent(in) :: Arg1, Arg2        ! input arguments
    real,    intent(out):: Value             ! output value
    logical, optional, intent(in):: DoExtrapolate ! optional extrapolation

    real:: Value_V(1)

    !--------------------------------------------------------------------------
    call interpolate_arg_array(iTable, [Arg1, Arg2], Value_V, DoExtrapolate)
    Value = Value_V(1)

  end subroutine interpolate_arg2_scalar
  !============================================================================

  subroutine interpolate_arg3_scalar(iTable, Arg1, Arg2, Arg3, Value, &
       DoExtrapolate)
    !$acc routine seq

    ! Return the scalar Value corresponding to arguments Arg1, Arg2 and Arg3.
    ! If DoExtrapolate is not present, stop with an error if the arguments
    ! are out of range. If it is present and false, return the value of
    ! the closest element in the table. If it is present and true, do a
    ! linear extrapolation.

    integer, intent(in) :: iTable            ! table
    real,    intent(in) :: Arg1, Arg2, Arg3  ! input arguments
    real,    intent(out):: Value             ! output value
    logical, optional, intent(in):: DoExtrapolate ! optional extrapolation

    real:: Value_V(1)

    !--------------------------------------------------------------------------
    call interpolate_arg_array(iTable, [Arg1, Arg2, Arg3], Value_V, &
         DoExtrapolate)
    Value = Value_V(1)

  end subroutine interpolate_arg3_scalar
  !============================================================================

  subroutine interpolate_arg4_scalar(iTable, Arg1, Arg2, Arg3, Arg4, Value, &
       DoExtrapolate)
    !$acc routine seq

    ! Return the scalar Value corresponding to arguments Arg1 ... Arg4.
    ! If DoExtrapolate is not present, stop with an error if the arguments
    ! are out of range. If it is present and false, return the value of
    ! the closest element in the table. If it is present and true, do a
    ! linear extrapolation.

    integer, intent(in) :: iTable                 ! table
    real,    intent(in) :: Arg1, Arg2, Arg3, Arg4 ! input arguments
    real,    intent(out):: Value                  ! output value
    logical, optional, intent(in):: DoExtrapolate ! optional extrapolation

    real:: Value_V(1)

    !--------------------------------------------------------------------------
    call interpolate_arg_array(iTable, [Arg1, Arg2, Arg3, Arg4], Value_V, &
         DoExtrapolate)
    Value = Value_V(1)

  end subroutine interpolate_arg4_scalar
  !============================================================================

  subroutine interpolate_arg5_scalar(iTable, Arg1, Arg2, Arg3, Arg4, Arg5, &
       Value, DoExtrapolate)
    !$acc routine seq

    ! Return the scalar Value corresponding to arguments Arg1 ... Arg5.
    ! If DoExtrapolate is not present, stop with an error if the arguments
    ! are out of range. If it is present and false, return the value of
    ! the closest element in the table. If it is present and true, do a
    ! linear extrapolation.

    integer, intent(in) :: iTable                       ! table
    real,    intent(in) :: Arg1, Arg2, Arg3, Arg4, Arg5 ! input arguments
    real,    intent(out):: Value                        ! output value
    logical, optional, intent(in):: DoExtrapolate       ! extrapolation

    real:: Value_V(1)

    !--------------------------------------------------------------------------
    call interpolate_arg_array(iTable, [Arg1, Arg2, Arg3, Arg4, Arg5], &
         Value_V, DoExtrapolate)
    Value = Value_V(1)

  end subroutine interpolate_arg5_scalar
  !============================================================================

  subroutine interpolate_arg_array(iTable, ArgIn_I, Value_V, DoExtrapolate)
    !$acc routine seq

    ! Return the array of values Value_V corresponding to arguments
    ! ArgIn_I in iTable. Use linear interpolation.
    ! If DoExtrapolate is not present, stop with an error if the arguments
    ! are out of range. If it is present and false, return the value of
    ! the closest element in the table. If it is present and true, do a
    ! linear extrapolation.

    integer, intent(in) :: iTable            ! table
    real,    intent(in) :: ArgIn_I(:)        ! input arguments
    real,    intent(out):: Value_V(:)        ! output values

    logical, optional, intent(in):: DoExtrapolate ! optional extrapolation

    real :: Arg_I(5)
    type(TableType), pointer:: Ptr

    integer :: MinIndex_I(5)

#ifndef _OPENACC
    character(len=*), parameter:: NameSub = 'interpolate_arg_array'
    !--------------------------------------------------------------------------
    if(iTable < 1 .or. iTable > nTable) call CON_stop(NameSub// &
         ': incorrect value for iTable=', iTable)
#endif
    MinIndex_I = 1

    Ptr => Table_I(iTable)
    Arg_I(1:Ptr%nIndex) = ArgIn_I(1:Ptr%nIndex)

    ! This line is broken so that emacs does not get confused
    where(Ptr%IsLogIndex_I) &
         Arg_I(1:Ptr%nIndex) = log10(Arg_I(1:Ptr%nIndex))

    ! Normalize Arg_I for interpolation if uniform
    where(Ptr%IsUniform_I) &
         Arg_I(1:Ptr%nIndex) = &
         (Arg_I(1:Ptr%nIndex) - Ptr%IndexMin_I)/Ptr%dIndex_I  + 1

    ! Interpolate values
    if(allocated(Ptr%Value4_VC))then
       Value_V = interpolate_vector4(Ptr%Value4_VC, &
            Ptr%nValue, Ptr%nIndex, MinIndex_I(1:Ptr%nIndex), Ptr%nIndex_I, &
            Arg_I(1:Ptr%nIndex),  &
            Ptr%Index1_I, Ptr%Index2_I, Ptr%Index3_I, Ptr%Index4_I, &
            Ptr%Index5_I, &
            DoExtrapolate = DoExtrapolate)
    else
       Value_V = interpolate_vector(Ptr%Value_VC, &
            Ptr%nValue, Ptr%nIndex, MinIndex_I(1:Ptr%nIndex), Ptr%nIndex_I, &
            Arg_I(1:Ptr%nIndex),  &
            Ptr%Index1_I, Ptr%Index2_I, Ptr%Index3_I, Ptr%Index4_I, &
            Ptr%Index5_I, &
            DoExtrapolate = DoExtrapolate)
    end if
  end subroutine interpolate_arg_array
  !============================================================================
  subroutine interpolate_with_known_val(&
       iTable, iVal, ValIn, Arg2In, Value_V, &
       Arg1Out, Arg1In, Arg2Out, DoExtrapolate)

    ! Return the array of values Value_V corresponding to the argument
    ! Arg2In in iTable, with Arg1Out calculated from the condition, that
    ! Value_V(iVal) equals the given value, ValIn.
    !
    ! Use linear interpolation.
    ! If DoExtrapolate is not present, stop with an error if the arguments
    ! are out of range. If it is present and false, return the value of
    ! the closest element in the table. If it is present and true, do a
    ! linear extrapolation.

    integer, intent(in) :: iTable            ! table
    integer, intent(in) :: iVal              ! which value is known
    real,    intent(in) :: ValIn             ! known table value
    real,    optional,   intent(in) :: Arg1In! first input argument if any
    real,    optional,   intent(in) :: Arg2In! second input argument if any
    real,    intent(out):: Value_V(:)        ! output values

    real, optional, intent(out) :: Arg1Out        ! optional calculated Arg
    real, optional, intent(out) :: Arg2Out        ! optional calculated Arg
    logical, optional, intent(in):: DoExtrapolate ! optional extrapolation

    real :: Arg1, Arg2
    type(TableType), pointer:: Ptr

    real    :: Dx1, Dx2, Dy1, Dy2
    integer :: i1, i2, j1, j2

    character(len=*), parameter:: NameSub = 'interpolate_with_known_val'
    !--------------------------------------------------------------------------
    if(iTable < 1 .or. iTable > nTable) call CON_stop(NameSub// &
         ': incorrect value for iTable=', iTable)

    Ptr => Table_I(iTable)
    if(present(Arg1In))then
       Arg1 = Arg1In

       If(Ptr%IsLogIndex_I(1)) Arg1 = log10(Arg1)

       call find_cell(1, Ptr%nIndex_I(1), &
            (Arg1- Ptr%IndexMin_I(1))/Ptr%dIndex_I(1)+ 1 , &
            i1, Dx1, &
            DoExtrapolate=DoExtrapolate, &
            StringError = 'Called from '//NameSub)

       i2 = i1 + 1; Dx2 = 1.0 - Dx1
       if(allocated(Ptr%Value4_VC))then
          call find_cell(1, Ptr%nIndex_I(2), ValIn, j1, Dy1, &
               Dx2*Ptr%Value4_VC(iVal,i1,:,1,1,1) + &
               Dx1*Ptr%Value4_VC(iVal,i2,:,1,1,1), &
               DoExtrapolate, 'Called from '//NameSub)

          j2 = j1 + 1; Dy2 = 1.0 - Dy1

          ! If value is outside table, use last value (works well for constant)
          Value_V = Dy2*( Dx2*Ptr%Value4_VC(:,i1,j1,1,1,1)   &
               +          Dx1*Ptr%Value4_VC(:,i2,j1,1,1,1))  &
               +    Dy1*( Dx2*Ptr%Value4_VC(:,i1,j2,1,1,1)   &
               +          Dx1*Ptr%Value4_VC(:,i2,j2,1,1,1))
       else
          call find_cell(1, Ptr%nIndex_I(2), ValIn, j1, Dy1, &
               Dx2*Ptr%Value_VC(iVal,i1,:,1,1,1) + &
               Dx1*Ptr%Value_VC(iVal,i2,:,1,1,1), &
               DoExtrapolate, 'Called from '//NameSub)

          j2 = j1 + 1; Dy2 = 1.0 - Dy1

          ! If value is outside table, use last value (works well for constant)
          Value_V = Dy2*( Dx2*Ptr%Value_VC(:,i1,j1,1,1,1)   &
               +          Dx1*Ptr%Value_VC(:,i2,j1,1,1,1))  &
               +    Dy1*( Dx2*Ptr%Value_VC(:,i1,j2,1,1,1)   &
               +          Dx1*Ptr%Value_VC(:,i2,j2,1,1,1))
       end if
       if(present(Arg2Out))then
          Arg2Out = (j1 - 1 + Dy1)*Ptr%dIndex_I(2) + Ptr%IndexMin_I(2)
          if(Ptr%IsLogIndex_I(2)) Arg2Out = 10**Arg2Out
       end if
       RETURN
    elseif(present(Arg2In))then
       Arg2 = Arg2In

       If(Ptr%IsLogIndex_I(2)) Arg2 = log10(Arg2)

       call find_cell(1, Ptr%nIndex_I(2), &
            (Arg2- Ptr%IndexMin_I(2))/Ptr%dIndex_I(2)+ 1 , &
            j1, Dy1, &
            DoExtrapolate=DoExtrapolate, &
            StringError = 'Called from '//NameSub)

       j2 = j1 + 1; Dy2 = 1.0 - Dy1
    else
       j1 = 1; j2 = 1; Dy1 = 1.0; Dy2 = 0.0
    end if

    if(allocated(Ptr%Value4_VC))then
       call find_cell(1, Ptr%nIndex_I(1), ValIn, i1, Dx1, &
            Dy2*Ptr%Value4_VC(iVal,:,j1,1,1,1) + &
            Dy1*Ptr%Value4_VC(iVal,:,j2,1,1,1), &
            DoExtrapolate, 'Called from '//NameSub)

       i2 = i1 + 1; Dx2 = 1.0 - Dx1

       ! If value is outside table, use last value (works well for constant)
       Value_V = Dy2*( Dx2*Ptr%Value4_VC(:,i1,j1,1,1,1)   &
            +          Dx1*Ptr%Value4_VC(:,i2,j1,1,1,1))  &
            +    Dy1*( Dx2*Ptr%Value4_VC(:,i1,j2,1,1,1)   &
            +          Dx1*Ptr%Value4_VC(:,i2,j2,1,1,1))
    else
       call find_cell(1, Ptr%nIndex_I(1), ValIn, i1, Dx1, &
            Dy2*Ptr%Value_VC(iVal,:,j1,1,1,1) + &
            Dy1*Ptr%Value_VC(iVal,:,j2,1,1,1), &
            DoExtrapolate, 'Called from '//NameSub)

       i2 = i1 + 1; Dx2 = 1.0 - Dx1

       ! If value is outside table, use last value (works well for constant)
       Value_V = Dy2*( Dx2*Ptr%Value_VC(:,i1,j1,1,1,1)   &
            +          Dx1*Ptr%Value_VC(:,i2,j1,1,1,1))  &
            +    Dy1*( Dx2*Ptr%Value_VC(:,i1,j2,1,1,1)   &
            +          Dx1*Ptr%Value_VC(:,i2,j2,1,1,1))
    end if
    if(present(Arg1Out))then
       Arg1Out = (i1 - 1 + Dx1)*Ptr%dIndex_I(1) + Ptr%IndexMin_I(1)
       if(Ptr%IsLogIndex_I(1)) Arg1Out = 10**Arg1Out
    end if

  end subroutine interpolate_with_known_val
  !============================================================================
  subroutine get_lookup_table(iTable, &
       iParamIn, Param, nParam, Param_I, nValue, nIndex, nIndex_I, &
       IndexMin_I, IndexMax_I, IsLogIndex_I, &
       NameVar, StringDescription, &
       Index1_I, Index2_I, Index3_I, Index4_I, Index5_I, Time)

    integer,                    intent(in) :: iTable     ! table index
    integer,          optional, intent(in) :: iParamIn   ! index of a parameter
    real,             optional, intent(out):: Param      ! one parameter
    integer,          optional, intent(out):: nParam     ! number of parameters
    real,             optional, intent(out):: Param_I(:) ! array of parameters
    integer,          optional, intent(out):: nValue     ! number of columns
    integer,          optional, intent(out):: nIndex     ! number of indexes
    integer,          optional, intent(out):: nIndex_I(:)! number of points
    real,             optional, intent(out):: IndexMin_I(:)   ! minimum indexes
    real,             optional, intent(out):: IndexMax_I(:)   ! maximum indexes
    logical,          optional, intent(out):: IsLogIndex_I(:) ! is logarithmic
    character(len=*), optional, intent(out):: NameVar         ! variable names
    character(len=*), optional, intent(out):: StringDescription ! description
    real,             optional, intent(out):: Index1_I(:) ! nonuniform indices
    real,             optional, intent(out):: Index2_I(:)
    real,             optional, intent(out):: Index3_I(:)
    real,             optional, intent(out):: Index4_I(:)
    real,             optional, intent(out):: Index5_I(:)
    real,             optional, intent(out):: Time

    ! Get various parameters of table iTable.

    type(TableType), pointer:: Ptr
    integer                 :: iParam
    character(len=*), parameter:: NameSub = 'get_lookup_table'
    !--------------------------------------------------------------------------
    if(iTable < 1 .or. iTable > nTable) call CON_stop(NameSub// &
         ': incorrect value for iTable=', iTable)

    Ptr => Table_I(iTable)
    if(present(StringDescription)) StringDescription = Ptr%StringDescription
    if(present(NameVar))           NameVar           = Ptr%NameVar
    if(present(nValue))            nValue            = Ptr%nValue
    if(present(nIndex))            nIndex            = Ptr%nIndex
    if(present(nIndex_I))          nIndex_I          = Ptr%nIndex_I
    if(present(IndexMin_I))        IndexMin_I        = Ptr%IndexMin_I
    if(present(IndexMax_I))        IndexMax_I        = Ptr%IndexMax_I
    if(present(IsLogIndex_I))      IsLogIndex_I      = Ptr%IsLogIndex_I
    if(present(nParam))            nParam            = Ptr%nParam
    if(present(Index1_I))          Index1_I          = Ptr%Index1_I
    if(present(Index2_I))          Index2_I          = Ptr%Index2_I
    if(present(Index3_I))          Index3_I          = Ptr%Index3_I
    if(present(Index4_I))          Index4_I          = Ptr%Index4_I
    if(present(Index5_I))          Index5_I          = Ptr%Index5_I
    if(present(Time))              Time              = Ptr%Time
    if(present(Param_I))then
       ! return an array of parameters
       if(Ptr%nParam == 0)then
          Param_I = -777.77
       elseif(Ptr%nParam <= size(Param_I))then
          Param_I(1:Ptr%nParam) = Ptr%Param_I
       else
          Param_I = Ptr%Param_I(1:size(Param_I))
       end if
    end if
    if(present(Param))then
       ! return one parameter, the first one by default
       iParam = 1
       if(present(iParamIn)) iParam = iParamIn
       if(iParam < 1 .or. iParam > Ptr%nParam)then
          Param = -777.77
       else
          Param = Ptr%Param_I(iParam)
       end if
    end if

  end subroutine get_lookup_table
  !============================================================================
  subroutine copy_lookup_table_to_gpu

    ! Copy new tables to the GPU

    integer :: iTable
    integer, save:: nTableOnGpu = 0
    !--------------------------------------------------------------------------
    if(nTable == nTableOnGpu) RETURN
    !$acc update device(nTable)
    !$acc update device(Table_I(nTableOnGPU+1:nTable))
    do iTable = nTableOnGpu + 1, nTable
       !$acc enter data attach(Table_I(iTable)%nIndex_I) &
       !$acc copyin(Table_I(iTable)%nIndex_I)

       !$acc enter data attach(Table_I(iTable)%IndexMin_I) &
       !$acc copyin(Table_I(iTable)%IndexMin_I)

       !$acc enter data attach(Table_I(iTable)%IndexMax_I) &
       !$acc copyin(Table_I(iTable)%IndexMax_I)

       !$acc enter data attach(Table_I(iTable)%dIndex_I) &
       !$acc copyin(Table_I(iTable)%dIndex_I)

       !$acc enter data attach(Table_I(iTable)%IsLogIndex_I) &
       !$acc copyin(Table_I(iTable)%IsLogIndex_I)

       !$acc enter data attach(Table_I(iTable)%Value_VC) &
       !$acc copyin(Table_I(iTable)%Value_VC)

       !$acc enter data attach(Table_I(iTable)%Value4_VC) &
       !$acc copyin(Table_I(iTable)%Value4_VC)

       !$acc enter data attach(Table_I(iTable)%Param_I) &
       !$acc copyin(Table_I(iTable)%Param_I)

       !$acc enter data attach(Table_I(iTable)%IsUniform_I) &
       !$acc copyin(Table_I(iTable)%IsUniform_I)

       !$acc enter data attach(Table_I(iTable)%Index1_I) &
       !$acc copyin(Table_I(iTable)%Index1_I)

       !$acc enter data attach(Table_I(iTable)%Index2_I) &
       !$acc copyin(Table_I(iTable)%Index2_I)

       !$acc enter data attach(Table_I(iTable)%Index3_I) &
       !$acc copyin(Table_I(iTable)%Index3_I)

       !$acc enter data attach(Table_I(iTable)%Index4_I) &
       !$acc copyin(Table_I(iTable)%Index4_I)

       !$acc enter data attach(Table_I(iTable)%Index5_I) &
       !$acc copyin(Table_I(iTable)%Index5_I)
    end do
    nTableOnGpu = nTable

  end subroutine copy_lookup_table_to_gpu
  !============================================================================
  subroutine test_lookup_table

    use ModNumConst, ONLY: cPi, cHalfPi

    ! testing the read_lookup_table_param is left for the functionality tests

    type(TableType), pointer :: Ptr, Ptr2
    integer :: iTable, iProc, iError, nParam, nIndex
    real :: Value_I(3), ValueGood_I(3), Arg, Value, ValueGood
    character(len=100):: String

    character(len=*), parameter:: NameSub = 'test_lookup_table'
    !--------------------------------------------------------------------------
    call MPI_comm_rank(MPI_COMM_WORLD,iProc,iError)

    if(iProc==0) write(*,*)'testing 2D lookup table'

    call init_lookup_table(&
         NameTable   = "RhoE",                  &
         NameCommand = "save",                  &
         NameVar     = "logrho e pXe pBe pPl zXe zBe zPl",  &
         NameFile    = "test_lookup_table1.out",&
         TypeFile    = "ascii",                 &
         nIndex_I    = [15, 10],              &
         IndexMin_I  = [0.001,   1.0],        &
         IndexMax_I  = [1000.0, 10.0],        &
         Param_I     = [ 54.0, 4.0, -4.0 ] )

    if(iProc==0) write(*,*)'testing i_lookup_table'

    iTable = i_lookup_table("xxx")
    if(iTable /= -1)then
       write(*,*)'iTable = ',iTable,' should be -1'
       call CON_stop(NameSub)
    end if
    iTable = i_lookup_table("RhoE")
    if(iTable /= 1)then
       write(*,*)'iTable = ',iTable,' should be 1'
       call CON_stop(NameSub)
    end if
    Ptr => Table_I(1)

    if(iProc==0) write(*,*)'testing get_lookup_table'
    call get_lookup_table(1, StringDescription=String, nParam=nParam, &
         Param_I=Value_I, iParamIn = 2, Param = Arg, nIndex=nIndex)
    if(nIndex /= 2)then
       write(*,*)'nIndex=',nIndex,' is different from 2'
       call CON_stop(NameSub)
    end if
    if(nParam /= 3) then
       write(*,*)'nParam=',nParam,' is different from 3'
       call CON_stop(NameSub)
    end if
    ValueGood_I = [ 54.0, 4.0, -4.0 ]
    if(any(Value_I /= ValueGood_I)) then
       write(*,*)'Param_I=', Value_I,' is different from ', ValueGood_I
       call CON_stop(NameSub)
    end if
    if(Arg /= ValueGood_I(2))then
       write(*,*)'Param=', Arg,' is different from ', ValueGood_I(2)
       call CON_stop(NameSub)
    end if

    if(iProc==0) write(*,*)'testing make_lookup_table'
    call make_lookup_table(iTable, eos_rho_e, MPI_COMM_WORLD)

    if(iProc==0) write(*,*)'testing interpolate_lookup_table'
    call interpolate_lookup_table(iTable, 1.0, 2.0, Value_I)
    ! rho=1.0 is exactly in the middle, e=2.0 is also an index, so exact result

    ValueGood_I = [ 4./3., 4./5., 3. ]
    if(any(abs(Value_I - ValueGood_I) > 1e-5))then
       write(*,*)'Value_I=',Value_I,' is different from ValueGood_I=', &
            ValueGood_I
       call CON_stop(NameSub)
    end if

    if(iProc==0) write(*,*)'testing load_lookup_table'

    ! Load the saved file into the second table
    call init_lookup_table(                     &
         NameTable   = "RhoE2"                 ,&
         NameCommand = "load"                  ,&
         NameFile    = "test_lookup_table1.out",&
         TypeFile    = "ascii")

    if(iProc==0) write(*,*)'testing i_lookup_table for table 2'
    iTable = i_lookup_table("RhoE2")
    if(iTable /= 2)then
       write(*,*)'iTable = ',iTable,' should be 2'
       call CON_stop(NameSub)
    end if
    Ptr2 => Table_I(iTable)
    call compare_tables

    if(iProc==0) write(*,*)'testing interpolate_lookup_table for loaded table'
    call interpolate_lookup_table(2, 1.0, 2.0, Value_I)
    if(any(abs(Value_I - ValueGood_I) > 1e-5))then
       write(*,*)'Value_I=',Value_I,' is different from ValueGood_I=', &
            ValueGood_I
       call CON_stop(NameSub)
    end if

    ! Test with the condition that Value_I(3)=3.0 and find first argument Arg
    call interpolate_lookup_table(2, 3, 3.0, 2.0, Value_I, Arg)
    if(any(abs(Value_I - ValueGood_I) > 1e-5) .or. abs(Arg - 1.0) >  1e-5)then
       write(*,*)'Value_I=',Value_I, ' Arg =', Arg,&
            ' are different from ValueGood_I=',ValueGood_I, ' ArgGood = 1.0'
       call CON_stop(NameSub)
    end if

    ! Test with the condition that Value_I(3)=3.0 and find first argument Arg
    call interpolate_lookup_table(2, 3, 3.0, &
         Value_V=Value_I, Arg1In=1.0, Arg2Out=Arg)
    if(any(abs(Value_I - ValueGood_I) > 1e-5) .or. abs(Arg - 2.0) >  1e-5)then
       write(*,*)'Value_I=',Value_I, ' Arg =', Arg,&
            ' are different from ValueGood_I=',ValueGood_I, ' ArgGood = 2.0'
       call CON_stop(NameSub)
    end if

    ! Test 1D lookup table
    if(iProc==0) write(*,*)'testing 1D lookup table'

    call init_lookup_table(&
         NameTable   = "RadCool",               &
         NameCommand = "save",                  &
         NameVar     = "logTe lambdaT ChiantiVersion",  &
         NameFile    = "test_lookup_table_radcool.out",&
         TypeFile    = "ascii",                 &
         nIndex_I    = [100],                 &
         IndexMin_I  = [1e3],                 &
         IndexMax_I  = [1e8],                 &
         Param_I     = [7.0] )

    iTable = i_lookup_table("RadCool")
    if(iTable /= 3)then
       write(*,*)'iTable = ',iTable,' should be 3'
       call CON_stop(NameSub)
    end if
    Ptr => Table_I(iTable)

    if(iProc==0) write(*,*)'testing make_lookup_table_1d'
    call make_lookup_table_1d(iTable, radcool_te, MPI_COMM_WORLD)

    if(iProc==0) write(*,*)'testing interpolate_lookup_table'
    call interpolate_lookup_table(iTable, 1e6, Value)

    ValueGood = 0.6
    if(abs(Value - ValueGood) > 1e-5)then
       write(*,*)'Value=',Value,' is different from ValueGood=',ValueGood
       call CON_stop(NameSub)
    end if

    if(iProc==0) write(*,*)'testing 1D load_lookup_table'

    ! Load the saved file into the fourth table
    call init_lookup_table(                     &
         NameTable   = "RadCool2"              ,&
         NameCommand = "load"                  ,&
         NameFile    = "test_lookup_table_radcool.out",&
         TypeFile    = "ascii")

    if(iProc==0) write(*,*)'testing i_lookup_table for table 4'
    iTable = i_lookup_table("RadCool2")
    if(iTable /= 4)then
       write(*,*)'iTable = ',iTable,' should be 4'
       call CON_stop(NameSub)
    end if
    Ptr2 => Table_I(iTable)
    call compare_tables

    if(iProc==0) write(*,*)'testing interpolate_lookup_table for loaded table'
    call interpolate_lookup_table(iTable, 1e6, Value)
    if(abs(Value - ValueGood) > 1e-5)then
       write(*,*)'Value=', Value,' is different from ValueGood=', ValueGood
       call CON_stop(NameSub)
    end if

    ! Test 3D lookup table
    if(iProc==0) write(*,*)'testing 3D lookup table'

    call init_lookup_table(&
         NameTable   = "B0(r,Lon,Lat)",                &
         NameCommand = "save",                         &
         NameVar     = "logR Lon Lat B0x B0y B0z",     &
         NameFile    = "test_lookup_table_b0.out",     &
         TypeFile    = "ascii",                        &
         nIndex_I    = [11,21,11],                   &
         IndexMin_I  = [ 1.0, 0.0, -cHalfPi],        &
         IndexMax_I  = [36.0, cPi, +cHalfPi]         )

    iTable = i_lookup_table("B0(r,Lon,Lat)")
    if(iTable /= 5)then
       write(*,*)'iTable = ',iTable,' should be 5'
       call CON_stop(NameSub)
    end if
    Ptr => Table_I(iTable)

    if(iProc==0) write(*,*)'testing make_lookup_table_3d'
    call make_lookup_table_3d(iTable, b0_rlonlat, MPI_COMM_WORLD)

    if(iProc==0) write(*,*)'testing interpolate_lookup_table'
    call interpolate_lookup_table(iTable, 6.0, cHalfPi, 0.0, Value_I)

    ValueGood_I = [0.0, 6.0, 0.0]
    if(any(abs(Value_I - ValueGood_I) > 1e-5))then
       write(*,*)'Value_I=',Value_I,' is different from ValueGood_I=', &
            ValueGood_I
       call CON_stop(NameSub)
    end if

    if(iProc==0) write(*,*)'testing 3D load_lookup_table'

    ! Load the saved file into the fourth table
    call init_lookup_table(                        &
         NameTable   = "B02",                      &
         NameCommand = "load",                     &
         NameFile    = "test_lookup_table_b0.out", &
         TypeFile    = "ascii")

    if(iProc==0) write(*,*)'testing i_lookup_table for table 6'
    iTable = i_lookup_table("B02")
    if(iTable /= 6)then
       write(*,*)'iTable = ',iTable,' should be 6'
       call CON_stop(NameSub)
    end if
    Ptr2 => Table_I(iTable)
    call compare_tables

    if(iProc==0) write(*,*)'testing interpolate_lookup_table for loaded table'
    call interpolate_lookup_table(iTable, 6.0, 0.0, cHalfPi, Value_I)
    ValueGood_I = [0.0, 0.0, 6.0]
    if(abs(Value - ValueGood) > 1e-5)then
       write(*,*)'Value=', Value,' is different from ValueGood=', ValueGood
       call CON_stop(NameSub)
    end if

    ! Test non-uniform 1D table
    call init_lookup_table(                               &
         NameTable   = "NonU1D"                          ,&
         NameCommand = "load"                            ,&
         NameFile    = "test_lookup_table_nonuniform.dat",&
         TypeFile    = "ascii")
    iTable = i_lookup_table("NonU1D")
    if(iTable /= 7)then
       write(*,*)'iTable = ',iTable,' should be 7'
       call CON_stop(NameSub)
    end if
    Ptr => Table_I(iTable)
    ! Test interpolate_lookup_table at provided index.
    call interpolate_lookup_table(iTable, 4.0, Value_I)
    ValueGood_I = [8.0, 2.0, 1.0]
    if(any(abs(Value_I - ValueGood_I) > 1e-5))then
       if(iProc==0)write(*,*)'Testing values at provided index:'
       if(iProc==0)write(*,*)'Value_I=', Value_I,' is different from ValueGood_I=', ValueGood_I
       call CON_stop(NameSub)
    end if
    ! Test interpolate_lookup_table at skipped index.
    call interpolate_lookup_table(iTable, 3.0, Value_I)
    ValueGood_I = [7.0, 3.0, 1.0]
    if(any(abs(Value_I - ValueGood_I) > 1e-5))then
       if(iProc==0)write(*,*)'Testing values at interpolated index:'
       if(iProc==0)write(*,*)'Value_I=', Value_I,' is different from ValueGood_I=', ValueGood_I
       call CON_stop(NameSub)
    end if

  contains
    !==========================================================================

    subroutine compare_tables

      ! Compare tables pointed at by Ptr and Ptr2

      !------------------------------------------------------------------------
      if(Ptr2%StringDescription /= Ptr%StringDescription) &
           call CON_stop(NameSub // &
           ' Description='//trim(Ptr2%StringDescription)// &
           ' is different from '// trim(Ptr%StringDescription))
      if(Ptr2%NameVar /= Ptr%NameVar) call CON_stop(NameSub // &
           ' NameVar='//trim(Ptr2%NameVar)//' is different from '// &
           trim(Ptr2%NameVar))
      if(Ptr2%nValue /= Ptr%nValue)then
         write(*,*)'nValue=',Ptr2%nValue,' is different from ',Ptr%nValue
         call CON_stop(NameSub)
      end if
      if(Ptr2%nIndex /= Ptr%nIndex)then
         write(*,*)'nIndex=',Ptr2%nIndex,' is different from ',Ptr%nIndex
         call CON_stop(NameSub)
      end if
      if(any(abs(Ptr2%IndexMin_I - Ptr%IndexMin_I) > 1e-5))then
         write(*,*)'IndexMin_I=',Ptr2%IndexMin_I,' is different from ',&
              Ptr%IndexMin_I
         call CON_stop(NameSub)
      end if
      if(any(abs(Ptr2%IndexMax_I - Ptr%IndexMax_I) > 1e-5))then
         write(*,*)'IndexMax_I=',Ptr2%IndexMax_I,' is different from ',&
              Ptr%IndexMax_I
         call CON_stop(NameSub)
      end if
      if(any(abs(Ptr2%dIndex_I - Ptr%dIndex_I) > 1e-6))then
         write(*,*)'dIndex_I=',Ptr2%dIndex_I,' is different from ',&
              Ptr%dIndex_I
         call CON_stop(NameSub)
      end if
      if(any(Ptr2%IsLogIndex_I .neqv. Ptr%IsLogIndex_I))then
         write(*,*)'IsLogIndex_I=',Ptr2%IsLogIndex_I,' is different from ',&
              Ptr%IsLogIndex_I
         call CON_stop(NameSub)
      end if
      if(Ptr2%nParam /= Ptr%nParam)then
         write(*,*)'nParam=',Ptr2%nParam,' is different from ',Ptr%nParam
         call CON_stop(NameSub)
      end if
      if(Ptr%nParam > 0)then
         if(any(Ptr2%Param_I /= Ptr%Param_I))then
            write(*,*)'Param_I=',Ptr2%Param_I,' is different from ',Ptr%Param_I
            call CON_stop(NameSub)
         end if
      end if

    end subroutine compare_tables
    !==========================================================================

  end subroutine test_lookup_table
  !============================================================================

  subroutine eos_rho_e(iTable, rho, e, p_I)
    ! This is an example for the subroutine passed to make_lookup_table

    integer, intent(in):: iTable
    real, intent(in)   :: rho, e
    real, intent(out)  :: p_I(:)
    !--------------------------------------------------------------------------
    p_I(1) = (2./3.)*e
    p_I(2) = (2./5.)*e
    p_I(3) = e + rho

  end subroutine eos_rho_e
  !============================================================================

  subroutine radcool_te(iTable, Te, Value_I)
    ! This is an example for the subroutine passed to make_lookup_table_1d

    integer, intent(in):: iTable
    real, intent(in)   :: Te
    real, intent(out)  :: Value_I(:)
    !--------------------------------------------------------------------------
    Value_I(1) = 0.1*log10(Te)

  end subroutine radcool_te
  !============================================================================

  subroutine b0_rlonlat(iTable, r, Lon, Lat, B0_D)
    ! This is an example for the subroutine passed to make_lookup_table_3d

    integer, intent(in):: iTable
    real, intent(in)   :: r, Lon, Lat
    real, intent(out)  :: B0_D(:)
    !--------------------------------------------------------------------------
    B0_D(1) = r*cos(Lon)*cos(Lat)
    B0_D(2) = r*sin(Lon)*cos(Lat)
    B0_D(3) = r*sin(Lat)

  end subroutine b0_rlonlat
  !============================================================================

end module ModLookupTable
!==============================================================================
