module io

  !Impose strong typing
  use trace, only : trace_entry, trace_exit,trace_stack,trace_finalise
  use comms, only : rank,nprocs,on_root_node,comms_abort,max_version_length,comms_arch,COMMS_VERSION &
       & ,COMMS_LIBRARY_VERSION,COMMS_FINALISE,comms_barrier
  use iso_fortran_env, only :real64,compiler_version
  implicit none


  logical,           public                :: no_param,read_params=.true.
  character(len=128),public                :: version  = "N_BODY v.1.0, Z. Hawkhead"
  character(len=128),public                :: info     = "Parallel code for n-body gravitational simulations."
  integer,           public,parameter      :: stdout = 984
  integer,           public,parameter      :: dp = real64
  real,              public,parameter      :: pi=3.1415926535_dp
  real(dp),          public,parameter      :: G= 6.67430E-11
  logical,           public                :: file_exists
  logical,           public                :: file_exists_struct

  character(100),dimension(:),allocatable  :: present_array

  character(100),dimension(:),allocatable  :: keys_array
  character(100),dimension(:),allocatable  :: keys_description
  character(100),dimension(:),allocatable  :: keys_default
  character(100),dimension(:),allocatable  :: keys_allowed
  character(100),dimension(:),allocatable  :: keys_type


  integer                                  :: max_params=1

  type  parameters
     ! %Begin: parameters

     !Calculation parameters
     real(dp)         :: calc_len          = 100.0
     real(dp)         :: time_step         = 0.5

     logical          :: debuging         = .false.
     logical          :: dry_run           = .false.

     character(len=30) :: diff_method = 'euler'
     real(dp) :: epsilon =    0.0000000000000000     
     ! %End: parameters
  end type parameters

  ! %Begin: keys
  character(len=30),parameter,public :: key_calc_len         = "calculation_length"
  character(len=30),parameter,public :: key_time_step        = "time_step"
  character(len=30),parameter,public :: key_debug            = "debug"
  character(len=30),parameter,public :: key_dry_run          = "dryrun"
  character(len=30),parameter,public :: key_diff_method      = 'ode_solver'
  character(len=30),parameter,public ::key_epsilon   = 'pot_soften'
  ! %End: keys


  type results
     real(dp)  :: birth_rate

  end type results


  type(parameters),public,save             :: current_params
  integer,parameter::max_keys=           6
  ! %End: max_param



  type structure
     real(dp)         ,allocatable,dimension(:,:) :: positions     ! Grid of positions, all initialisations default to this
     real(dp)         ,allocatable,dimension(:)   :: masses        ! Massess of all the particles
     real(dp)         ,allocatable,dimension(:,:) :: init_velocity ! Initial velocities of all particles
     character(len=15),allocatable,dimension(:)   :: labels        ! Labels, can be custom or automatically defined
     integer                                      :: n_bodies      ! Number of particles, used to define the arrays
     logical                                      :: init_radial   ! Gives the radial initialisation scheme
     logical                                      :: init_cart     ! Gives the cartesian initialisation scheme
     logical                                      :: init_grid     ! Gives the grid initialisation scheme
     integer                                      :: nx, ny, nz    ! Grid dimensions for grid initialisation n_bodies=nx*ny*nz
     real(dp)                                     :: dnx,dny,dnz   ! Grid spacings for grid initalisation

  end type structure

  type(structure),public,save :: current_structure

  !-------------------------------------------------------!
  !              P U B L I C  R O U T I N E S             !
  !-------------------------------------------------------!
  public :: io_initialise
  public :: io_errors

contains



  subroutine io_initialise()
    !==============================================================================!
    !                          I O _ I N I T I A L I S E                           !
    !==============================================================================!
    ! Subroutine for initialising all input/output for the parallel code N_BODY    !
    !------------------------------------------------------------------------------!
    ! Arguments:                                                                   !
    !           None                                                               !
    !------------------------------------------------------------------------------!
    ! Author:   Z. Hawkhead  19/01/2020                                            !
    !==============================================================================!
    ! Should be called in top level file
    implicit none

    integer :: i ! counters
    character(10) :: line
    integer :: stat
    integer :: demo_length ,year, age, age_pop
    logical :: demo_file


    ! Some junk variables


    call trace_entry("io_initialise")
    call io_cl_parser() ! Read the commandline arguments



    ! Get the length of the parameters file
    inquire(file="param.n_body",exist=file_exists)

    if (file_exists)then
       open(unit=1,file="param.n_body",iostat=stat,status="OLD",access="stream",form="formatted")
       do while (stat.eq.0)
          read(1,'(A60)',iostat=stat) line
          max_params=max_params+1
       end do
       close(1)
    end if
    ! Allocate space for the params array
    allocate(present_array(0:max_params))

    ! Fist things first, try to read paramteters
    if (read_params) call io_read_param(current_params)


    ! Try and open the struct file
    inquire(file="struct.n_body",exist=file_exists_struct)
    if (.not.file_exists_struct)call io_errors("Error in I/O: No structure file")
    if (file_exists_struct) call io_read_structure()


    ! Open up the main file for the output
    open(stdout,file="out.n_body",RECL=8192,form="FORMATTED",access="append")


    call io_header()
    call trace_exit("io_initialise")


    return
  end subroutine io_initialise


  subroutine io_read_param(dummy_params)
    !==============================================================================!
    !                          I O _ R E A D _ P A R A M                           !
    !==============================================================================!
    ! Subroutine for reading parameters from the file "param.n_body" to the        !
    ! parameter type current_params                                                !
    !------------------------------------------------------------------------------!
    ! Arguments:                                                                   !
    !           dummy_params,      intent :: inout                                 !
    !------------------------------------------------------------------------------!
    ! Author:   Z. Hawkhead  18/11/2020                                            !
    !==============================================================================!
    implicit none
    !The inout stuff
    type(parameters),intent(inout)  :: dummy_params


    !The boring stuff to make the whole shebang work
    integer           :: stat
    integer           :: read_stat=0
    integer           :: i,j,k           !counter

    character(len=60) :: line        ! charcter string into which each line is read, overwritten in loop
    character(len=30) :: key         ! the keyword used
    character(len=30) :: param       ! the value of the param
    logical           :: comment     ! boolean for comment line, will skip
    real(dp)          :: real_dump   ! a dump for handling scientific
    call trace_entry("io_read_param")

    !Open the parameter file
    if (file_exists) then
       open(unit=1,file="param.n_body",iostat=stat,status="OLD",access="stream",form="formatted")


       if (stat.ne.0) call io_errors("Error in I/O: Open file 'param.n_body'")
       ! now we can do the reading
       k=0
       do i=0,max_params
          !first thing, read new line into 'line' variable
          read(1,'(A)',iostat=read_stat) line

          write(present_array(i),*)i

          ! Check for blank line
          if (trim(line).eq."") cycle


          !Read everying into a thing
          call io_freeform_read(line,key,param,comment)

          if (comment) cycle ! skip if comment

          !Do some trimming
          key=adjustl(trim(io_case(key)))
          param=adjustl(trim(io_case(param)))


          ! %Begin: case_read
          select case(key)
          case(key_calc_len)
             read(param,*,iostat=stat) dummy_params%calc_len
             if (stat.ne.0) call io_errors("Error in I/O: Error parsing value: "//param)
             present_array(i)=key
          case(key_time_step)
             read(param,*,iostat=stat) dummy_params%time_step
             if (stat.ne.0) call io_errors("Error in I/O: Error parsing value: "//param)
             present_array(i)=key
          case(key_debug)
             read(param,*,iostat=stat) dummy_params%debuging
             if (stat.ne.0) call io_errors("Error in I/O: Error parsing value: "//param)
             present_array(i)=key
          case(key_dry_run)
             read(param,*,iostat=stat) dummy_params%dry_run
             if (stat.ne.0) call io_errors("Error in I/O: Error parsing value: "//param)
             present_array(i)=key
          case(key_diff_method                                                 )
             read(param,*,iostat=stat) dummy_params%diff_method
             if (stat.ne.0) call io_errors("Error in I/O: Error parsing value: "//param)
             present_array(i)=key
          case default
             call io_errors("Error in I/O: Error parsing keyword: "//key)
          end select

       end do


    else
       no_param=.true.
       call trace_exit("io_read_param")
       return
    end if
    ! Check for duplicates
    do i=0,max_params
       do j=0,max_params
          if (i.eq.j)cycle
          !print*,present_array
          if (present_array(i).eq.present_array(j))then
             call io_errors("Error in I/O: Duplicate parameter found: "//present_array(i))
          end if
       end do
    end do

    call trace_exit("io_read_param")
    return
  end subroutine io_read_param




  subroutine io_freeform_read(line_unparsed,key,val,com)
    !==============================================================================!
    !                       I O _ F R E E F O R M _ R E A D                        !
    !==============================================================================!
    ! Subroutine for parsing keys and params from general line with delimiter      !
    !------------------------------------------------------------------------------!
    ! Arguments:                                                                   !
    !           line_unparsed,     intent :: in                                    !
    !           key,               intent :: out                                   !
    !           val,               intent :: out                                   !
    !           com,               intent :: inout                                 !
    !------------------------------------------------------------------------------!
    ! Author:   Z. Hawkhead  19/01/2020                                            !
    !==============================================================================!
    !subroutine that parses the lines from param.pop into the key and the the value
    implicit none
    character(*),intent(in)           :: line_unparsed
    character(*),intent(out)          :: key
    character(*),intent(out)          :: val
    logical,     intent(inout)        :: com

    integer                           :: j


    if (line_unparsed(1:1).eq."!" .or. line_unparsed(1:1).eq."#") then
       com=.true.
       return
    else
       com=.false.
    end if

    do j=1,len_trim(line_unparsed)

       if (line_unparsed(j:j).eq.':' .or. line_unparsed(j:j).eq."=")then
          key=line_unparsed(1:j-1)
          val=line_unparsed(j+1:len_trim(line_unparsed))

          exit
       else if (j.eq.len_trim(line_unparsed))then
          call io_errors("Error in I/O: Error parsing line:  "//trim(line_unparsed))
       end if
    end do



    return
  end subroutine io_freeform_read


  subroutine io_errors(message)
    !==============================================================================!
    !                              I O _ E R R O R S                               !
    !==============================================================================!
    ! Subroutine handling all errors writing to the errors file                    !
    !------------------------------------------------------------------------------!
    ! Arguments:                                                                   !
    !           message,           intent :: in                                    !
    !------------------------------------------------------------------------------!
    ! Author:   Z. Hawkhead  19/01/2020                                            !
    !==============================================================================!
    implicit none
    character(*)       :: message

    ! internal variable for rank processing
    character(len=20)  :: file_name

    write(file_name,'("err.",I0.4,".n_body")') rank

    open(2,file=trim(file_name),RECL=8192)
    write(*,*)"Error: called io_abort"
    write(2,*) message

    call trace_stack(2,rank)
    stop
    return
  end subroutine io_errors

  function io_case( string , upper) result (new)
    !==============================================================================!
    !                                I O _ C A S E                                 !
    !==============================================================================!
    ! Low level subroutine for modifying the case of user given arguments to       !
    ! lowercase                                                                    !
    !------------------------------------------------------------------------------!
    ! Arguments:                                                                   !
    !           None                                                               !
    !------------------------------------------------------------------------------!
    ! Result:                                                                      !
    !           strin                                                              !
    !------------------------------------------------------------------------------!
    ! Author:   Z. Hawkhead  19/01/2020                                            !
    !==============================================================================!
    character(len=*)           :: string

    character(len=len(string)) :: new

    integer                    :: i
    integer                    :: k
    integer                    :: length
    logical, optional          :: upper

    length = len(string)
    new    = string
    do i = 1,len(string)
       k = iachar(string(i:i))
       if ( k >= iachar('A') .and. k <= iachar('Z') ) then
          k = k + iachar('a') - iachar('A')
          new(i:i) = achar(k)
       end if
    end do

    if (present(upper))then
       if (upper) then
          do i = 1,len(string)
             k = iachar(string(i:i))
             if ( k >= iachar('a') .and. k <= iachar('z') ) then
                k = k + iachar('A') - iachar('a')
                new(i:i) = achar(k)
             end if
          end do
       end if
    end if
  end function io_case


  subroutine io_cl_parser()
    !==============================================================================!
    !                           I O _ C L _ P A R S E R                            !
    !==============================================================================!
    ! Subroutine for the handling of commandline arguments.                        !
    !------------------------------------------------------------------------------!
    ! Arguments:                                                                   !
    !           None                                                               !
    !------------------------------------------------------------------------------!
    ! Author:   Z. Hawkhead  03/02/2020                                            !
    !==============================================================================!
    implicit none
    integer          ::   nargs     !Number of args
    integer          ::   arg_index !The index
    character(50)    ::   name      !The name of the argument

    logical          ::   search = .false.
    logical          ::   help   = .false.
    nargs=command_argument_count()

    if (nargs.gt.0)then
       do arg_index=1,nargs

          call get_command_argument(arg_index,name)
          select case(adjustl(trim(name)))
          case("-h","--help")
             write(*,*) trim(version)
             write(*,*) trim(info)
             read_params=.false.
             call io_list_params(.false.)
             help=.true.
             if (arg_index.eq.nargs) call io_help
          case("-s","--search")
             write(*,*) trim(version)
             write(*,*) trim(info)
             write(*,*)
             read_params=.false.
             call io_list_params(.false.)
             search=.true.
             if (arg_index.eq.nargs) call io_help

          case("-v")
             write(*,*) trim(version)
             write(*,*) trim(info)
             read_params=.false.
             stop
          case("-d","--dryrun")
             current_params%dry_run=.true.
          case("-l","--list")
             write(*,*) trim(version)
             write(*,*) trim(info)
             read_params=.false.
             call io_list_params(.true.)
             stop
          case default
             if (help)then
                call io_help(name)
                help=.false.
             elseif(search)then
                call io_search(io_case(name))
                search=.false.
                stop
             else
                write(*,*) "Unknown argument: ",adjustl(name)
                write(*,*) trim(version)
                write(*,*) trim(info)
                call io_help()
             end if
             stop
          end select
       end do
    end if
    return
  end subroutine io_cl_parser

  subroutine io_help(string)
    !==============================================================================!
    !                                I O _ H E L P                                 !
    !==============================================================================!
    ! Subroutine for printing the help information to the terminal.                !
    !------------------------------------------------------------------------------!
    ! Arguments:                                                                   !
    !           None                                                               !
    !------------------------------------------------------------------------------!
    ! Author:   Z. Hawkhead  03/02/2020                                            !
    !==============================================================================!
    implicit none
    character(*),optional   :: string
    integer                 :: i !counter
    logical                 :: found=.false.
    if (present(string))then

       do i=1,max_keys
          if (trim(keys_array(i)).eq.io_case(trim(string))) then
             found=.true.
             write(*,*)
             write(*,12) repeat("*",70-len(keys_array(i))/2),&
                  & trim(io_case(keys_array(i),upper=.true.)),repeat("*",70-len(keys_array(i))/2)
             write(*,*) trim(keys_description(i))
             write(*,*)
             write(*,*) "Allowed Values: ",trim(keys_allowed(i))
             write(*,*) "Default:        ",trim(keys_default(i))
             write(*,*)
12           format(1x,a,1x,a,1x,a)

             exit
          end if
       end do
       if (.not.found)then
          write(*,*)
          write(*,*) "************** NO MATCHING PARAMETERS **************"
          write(*,*)
       end if
    else
30     format(4x,A,T40,A)
       write(*,30) "-v","Print version information."
       write(*,30) "-h,--help   <keyword>","Get help and commandline options."
       write(*,30) "-s,--search <keyword>", "Search list of available parameters"
       write(*,30) "-l,--list","Get list of parameters avilable for the user."
       write(*,30) "-d,--dryrun","Run calculation to check input files"
    end if
    return
  end subroutine io_help


  subroutine io_search(string)
    !==============================================================================!
    !                              I O _ S E A R C H                               !
    !==============================================================================!
    ! Subroutine for using the command to search for availible variables that      !
    ! the user can change in a calculation                                         !
    !------------------------------------------------------------------------------!
    ! Arguments:                                                                   !
    !           string,            intent :: in                                    !
    !------------------------------------------------------------------------------!
    ! Author:   Z. Hawkhead  29/02/2020                                            !
    !==============================================================================!
    implicit none
    character(*)     :: string
    logical          :: found
    integer          :: i,scan_res

    do i=1,max_keys
       scan_res=index(trim(keys_array(i)),trim(string))

       if (scan_res.gt.0)then
          found=.true.
100       format(1x,A,T35,A)
          write(*,100)io_case(keys_array(i),.true.),keys_description(i)
       end if
    end do
    if (.not.found)then
       write(*,*)
       write(*,*) "************** NO MATCHING PARAMETERS **************"
       write(*,*)
    end if

    return
  end subroutine io_search

  subroutine io_list_params(print_flag)
    !==============================================================================!
    !                         I O _ L I S T _ P A R A M S                          !
    !==============================================================================!
    ! Subroutine used to print the total list of variables to the terminal and     !
    ! also to allocate all of the global arrays which contain the variable data.   !
    !------------------------------------------------------------------------------!
    ! Arguments:                                                                   !
    !           print_flag,        intent :: in                                    !
    !------------------------------------------------------------------------------!
    ! Author:   Z. Hawkhead  29/02/2020                                            !
    !==============================================================================!
    implicit none
    logical  :: print_flag

    character(60) :: junk
    integer       :: i ! loops
    ! Allocate all the arrays for the parameters
    allocate(keys_array(1:max_keys))
    allocate(keys_default(1:max_keys))
    allocate(keys_description(1:max_keys))
    allocate(keys_allowed(1:max_keys))




    ! assign the keys
    ! %Begin: assign_keys
    keys_array(1)=trim(key_calc_len)
    keys_array(2)=trim(key_time_step)
    keys_array(3)=trim(key_dry_run)
    keys_array(4)=trim(key_debug)
    keys_array(5)=trim(key_diff_method)
    keys_array(6)=trim(key_epsilon)
    ! %End: assign_keys

    ! %Begin: assign_default
    write(junk,*)current_params%calc_len
    keys_default(1)=trim(adjustl(junk))
    write(junk,*)current_params%time_step
    keys_default(2)=trim(adjustl(junk))
    write(junk,*)current_params%dry_run
    keys_default(3)=trim(adjustl(junk))
    write(junk,*)current_params%debuging
    keys_default(4)=trim(adjustl(junk))
    write(junk,*)current_params%diff_method
    keys_default(5)=trim(adjustl(junk))
    write(junk,*)current_params%epsilon                                                     
    keys_default(6)=trim(adjustl(junk))
    ! %End: assign_default

    ! %Begin: assign_description
    keys_description(1)="Length of the calculation in days"
    keys_description(2)="Size of the time step used for integration"
    keys_description(3)="Initialise calculation to check input files"
    keys_description(4)="Turn on profilling and debugging"
    keys_description(5)='Chose method for solving ODEs.'
    keys_description(6)='Amount of softening used in the gravitational potential to avoid singularities'
    ! %End: assign_description

    ! %Begin: assign_allowed
    keys_allowed(1)= "(any real) > 0"
    keys_allowed(2)= "(any real) > 0"
    keys_allowed(3)= "Boolean"
    keys_allowed(4)= "Boolean"
    keys_allowed(5)='euler'
    keys_allowed(6)='(any real) > 0'
    ! %End: assign_allowed

    ! do the loop for printing stuff

    if (print_flag)then
100    format(1x,A,T35,A)
       write(*,*)
       do i=1,max_keys
          write(*,100) io_case(trim(keys_array(i)),.true.),trim(keys_description(i))

       end do
    end if


    return
  end subroutine io_list_params


  subroutine io_header()
    !==============================================================================!
    !                          I O _  H E A D E R                                  !
    !==============================================================================!
    ! Subroutine used to write the io_header of the main Mandelbrot output file    !
    ! "out.mand".                                                                  !
    !------------------------------------------------------------------------------!
    ! Arguments:                                                                   !
    !           None                                                               !
    !------------------------------------------------------------------------------!
    ! Author:   Z. Hawkhead  16/08/2019                                            !
    !==============================================================================!
    implicit none
    integer::file
    integer :: maj_mpi,min_mpi,min_char
    character(len=max_version_length) :: mpi_c_version
    character(len=3) :: MPI_version_num
    character(len=100):: compile_version,cpuinfo
    character(len=5) :: opt

    call trace_entry("io_header")



    if (comms_arch.eq."MPI")then
       call COMMS_LIBRARY_VERSION(mpi_c_version)
       call COMMS_VERSION(min_mpi,maj_mpi)

       write(mpi_version_num,97)min_mpi,maj_mpi
97     format(i1,"."i1)
       min_char=scan(mpi_c_version,".")
       !print*, mpi_c_version,mpi_version_num
    end if

#ifdef __INTEL_COMPILER
#define compiler "Intel Compiler"

#endif
#ifdef __GFORTRAN__
#define compiler "GNU Fortran"
    !#define compile_version __VERSION__
#endif



#define opt opt_strat



    compile_version=compiler_version()
    if (compiler.eq."Intel Compiler")then
       compile_version=compiler_version()

       compile_version=trim(compile_version(87:97))
    end if

    write(stdout,*) '+==============================================================================================+'
    write(stdout,*) '|                       .                        ___                                           |'
    write(stdout,*) '|           .                      *          ,o88888                               *          |'
    write(stdout,*) '|                                          ,o8888888"           .                              |'
    write(stdout,*) '|  o                 ,:o:o:oooo.        ,8O88Pd8888"                                           |'
    write(stdout,*) '|                ,.::.::o:ooooOoOoO. ,oO8O8Pd888""                       *                     |'
    write(stdout,*) '|       *      ,.:.::o:ooOoOoOO8O8OOo.8OOPd8O8O"                                      .        |'
    write(stdout,*) '|             , ..:.::o:ooOoOOOO8OOOOo.FdO8O8"          *                     o                |'
    write(stdout,*) '|            , ..:.::o:ooOoOO8O888O8O,COCOO"                                                   |'
    write(stdout,*) '|           , . ..:.::o:ooOoOOOO8OOOOCOCO"                                                     |'
    write(stdout,*) '|    *      . ..:.::o:ooOoOoOO8O8OCCCC"o          88b 88       88""Yb  dP"Yb  8888b.  Yb  dP   |'
    write(stdout,*) '|              . ..:.::o:ooooOoCoCCC"o:o          88Yb88       88__dP dP   Yb  8I  Yb  YbdP    |'
    write(stdout,*) '|              . ..:.::o:o:,cooooCo"oo:o:         88 Y88       88""Yb Yb   dP  8I  dY   8P     |'
    write(stdout,*) '|           `   . . ..:.:cocoooo""o:o:::"         88  Y8 ooooo 88oodP  YbodP  8888Y"   dP      |'
    write(stdout,*) '|  o        .`   . ..::ccccoc""o:o:o:::"                                                       |'
    write(stdout,*) '|          :.:.    ,c:cccc"":.:.:.:.:."                                   .                    |'
    write(stdout,*) '|        ..:.:""`::::c:""..:.:.:.:.:."        o         .                                *     |'
    write(stdout,*) '|      ...:.".:.::::""    . . . . ."                                                           |'
    write(stdout,*) '|     .. . ....:."" `   .  . . ""                              .                   |           |'
    write(stdout,*) '|   . . . ....""                     o                                            -O-      .   |'
    write(stdout,*) '|   .. . .""      *                              *                                 |           |'
    write(stdout,*) '|  .                                                                    o                      |'
    write(stdout,*) '|                            .                                                                 |'
    write(stdout,*) '+----------------------------------------------------------------------------------------------+'
    write(stdout,'(1x,"|",10x,a,T97,"|")') trim(version)
    write(stdout,*) '+==============================================================================================+'
    write(stdout,*) "Compiled with ",compiler," ",Trim(compile_version), " on ", __DATE__, " at ",__TIME__
    write(stdout,*) "Communications architechture: ",comms_arch
    if (comms_arch.eq."MPI")then
       write(stdout,*) "MPI Version: ",mpi_c_version(1:min_char+1)
    end if
    write(stdout,*) "Optimisation Strategy: ",opt
    write(stdout,*)
    call trace_exit("io_header")
  end subroutine io_header


  subroutine io_dryrun()
    !==============================================================================!
    !                              I O _ D R Y R U N                               !
    !==============================================================================!
    ! Subroutine for handlind the dryrun command which allows for parameter        !
    ! checking                                                                     !
    !------------------------------------------------------------------------------!
    ! Arguments:                                                                   !
    !           None                                                               !
    !------------------------------------------------------------------------------!
    ! Author:   Z. Hawkhead  23/02/2020                                            !
    !==============================================================================!
    implicit none
    call trace_entry("io_dryrun")
    if (on_root_node)then
       write(stdout,*) " "
       write(stdout,'(28x,A)') "****************************************"
       write(stdout,'(28x,A)') "*                                      *"
       write(stdout,'(28x,A)') "*         Dryrun complete....          *"
       write(stdout,'(28x,A)') "*          No errors found             *"
       write(stdout,'(28x,A)') "*                                      *"
       write(stdout,'(28x,A)') "****************************************"
    end if

    call trace_exit("io_dryrun")
    call trace_exit("mc_pop")
    call COMMS_FINALISE()
    call trace_finalise(current_params%debuging,rank)
    stop

  end subroutine io_dryrun



  subroutine io_write_params()
    !==============================================================================!
    !                        I O _ W R I T E _ P A R A M S                         !
    !==============================================================================!
    ! Subroutine that writes the input parameters to the stdout out.pop            !
    !------------------------------------------------------------------------------!
    ! Arguments:                                                                   !
    !           None                                                               !
    !------------------------------------------------------------------------------!
    ! Author:   Z. Hawkhead  23/02/2020                                            !
    !==============================================================================!
    implicit none
    character(50)   :: sec_title
    integer         :: width=97,length
    ! Stuff for getting run time
    character(len=3),dimension(12)  :: months
    integer                         :: d_t(8)
    character*10                    :: b(3)


    call trace_entry("io_write_params")

    call date_and_time(b(1), b(2), b(3), d_t)
    months=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

    write(stdout,*)"+"//repeat("-",(width-15)/2)//" RUN STARTED "//repeat("-",(width-16)/2)//"+"
    write(stdout,1000) d_t(5),d_t(6),d_t(7),trim(months(d_t(2))),d_t(3),d_t(1)
    write(stdout,*)"+"//repeat("-",width-3)//"+"
    write(stdout,*) " "
1000 FORMAT(1x,"|",36x,i2.2,":",i2.2,":",i2.2,",",1x,A,1x,i2.2,1x,i4,37x,"|")
10  format(1x,A,T44,':',5x,I12,1x,A)    !integer
11  format(1x,A,T44,":",5x,f12.2,1x,A)  !real
12  format(1x,A,T44,":",5x,L12,1x,A)    !logical
13  format(1x,A,T44,":",5x,A12,1x,A)    !character
14  format(1x,A,T44,":",5x,ES12.2,1x,A)  !Science




    ! Do some printing about the parameters file

    if (file_exists)then
       write(stdout,*)"Parameters file found, using custom parameters"
    else
       write(stdout,*)"Parameters file not found, using defaults"
    end if


    call trace_exit("io_write_params")
    return
  end subroutine io_write_params




  subroutine io_write_results()


    implicit none
    logical               :: survived

    character(50)   :: sec_title
    integer         :: width=97,length
    ! Stuff for getting run time
    character(len=3),dimension(12)  :: months
    integer                         :: d_t(8)
    character*10                    :: b(3)
    call trace_entry("io_write_results")

10  format(1x,A,T44,':',5x,I9,1x,A)    !integer
11  format(1x,A,T44,":",5x,f9.2,1x,A)  !real
12  format(1x,A,T44,":",5x,L9,1x,A)    !logical



    call date_and_time(b(1), b(2), b(3), d_t)
    months=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

    write(stdout,*)"+"//repeat("-",(width-15)/2)//" RUN FINISHED "//repeat("-",(width-18)/2)//"+"
    write(stdout,1000) d_t(5),d_t(6),d_t(7),trim(months(d_t(2))),d_t(3),d_t(1)
    write(stdout,*)"+"//repeat("-",width-3)//"+"

1000 FORMAT(1x,"|",36x,i2.2,":",i2.2,":",i2.2,",",1x,A,1x,i2.2,1x,i4,37x,"|")

    call trace_exit("io_write_results")
    return
  end subroutine io_write_results

  function io_present(key) result(is_present)
    !==============================================================================!
    !                             I O _ P R E S E N T                              !
    !==============================================================================!
    ! Function used to determine if there is a keyword present in the input file   !
    !                                                                              !
    !------------------------------------------------------------------------------!
    ! Arguments:                                                                   !
    !           key,               intent :: in                                    !
    !------------------------------------------------------------------------------!
    ! Result:                                                                      !
    !           is_present                                                         !
    !------------------------------------------------------------------------------!
    ! Author:   Z. Hawkhead  23/02/2020                                            !
    !==============================================================================!
    implicit none
    logical      :: is_present
    character(*) :: key
    call trace_entry("io_present")
    if (any(present_array.eq.key))then
       is_present=.TRUE.
    else
       is_present=.FALSE.
    end if
    call trace_exit("io_present")
  end function io_present

  subroutine io_flush(unit)
    !==============================================================================!
    !                               I O _ F L U S H                                !
    !==============================================================================!
    ! Subroutine wrapper for the intrinsic function that forces the system to      !
    ! clear the cache so that I/O can be written to file. Mainly used with the     !
    ! GNU fortran compiler, programs compiled with that compiler tend to hold on   !
    ! to the cache much longer.                                                    !
    !------------------------------------------------------------------------------!
    ! Arguments:                                                                   !
    !           unit,              intent :: in                                    !
    !------------------------------------------------------------------------------!
    ! Author:   Z. Hawkhead  29/02/2020                                            !
    !==============================================================================!
    implicit none
    integer  :: unit
    !call trace_entry("io_flush")
    call flush(unit)
    !call trace_exit("io_flush")
    return
  end subroutine io_flush



  subroutine io_read_structure()
    implicit none
    integer :: n
    integer :: stat

    character(len=100) :: line

    call trace_entry("io_read_structure")
    


    
    call trace_exit("io_read_structure")
    return
  end subroutine io_read_structure
  


  subroutine io_cart_to_radial(pos_array)
    implicit none

    real(dp),dimension(:,:), intent(inout) :: pos_array
    real(dp),dimension(:,:),allocatable  :: rad_array 
    integer  :: stat

    call trace_entry("io_cart_to_radial")


    allocate(rad_array(1:current_structure%n_bodies,1:3),stat=stat)
    if (stat.ne.0)call io_errors("Error in I/O: allocate rad_array")

    rad_array(:,1)=sqrt(pos_array(:,1)**2 &
         & + pos_array(:,2)**2 &
         & + pos_array(:,3)**2)

    rad_array(:,2)=atan2(pos_array(:,2),pos_array(:,1))
    rad_array(:,3)=atan(sqrt((pos_array(:,1)**2 &
         & + pos_array(:,2)**2) /&
         & pos_array(:,3)))
    pos_array(:,1)=rad_array(:,1)
    pos_array(:,2)=rad_array(:,2)
    pos_array(:,3)=rad_array(:,3)

    if (allocated(rad_array))deallocate(rad_array)
    call trace_exit("io_cart_to_radial")
    return 
  end subroutine io_cart_to_radial


  subroutine io_radial_to_cart(pos_array)
    implicit none

    real(dp),dimension(:,:), intent(inout) :: pos_array
    real(dp),dimension(:,:),allocatable  :: rad_array 
    integer  :: stat
 
    call trace_entry("io_radial_to_cart")


    allocate(rad_array(1:current_structure%n_bodies,1:3),stat=stat)
    if (stat.ne.0)call io_errors("Error in I/O: allocate rad_array")

    rad_array(:,1)=pos_array(:,1)*sin(pos_array(:,3))*cos(pos_array(:,2))
    rad_array(:,2)=pos_array(:,1)*sin(pos_array(:,3))*cos(pos_array(:,2))
    rad_array(:,3)=pos_array(:,1)*cos(pos_array(:,3))

    pos_array(:,1)=rad_array(:,1)
    pos_array(:,2)=rad_array(:,2)
    pos_array(:,3)=rad_array(:,3)
    
    if (allocated(rad_array))deallocate(rad_array)
    call trace_exit("io_radial_to_cart")
    return
  end subroutine io_radial_to_cart


  
end module io
 
 
