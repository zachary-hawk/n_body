!---- File documented by Fortran Documenter, Z.Hawkhead
module pot
  use trace, only : trace_entry, trace_exit
  use io
  use comms
  implicit none


  type potential
     real(dp), allocatable,dimension(:,:) :: pot_array
     logical                              :: is_allocated = .false.
     !real(dp), allocatable,dimension(:,:) :: pot_dir
  end type potential




contains

  subroutine pot_allocate(dummy_pot,struct)
    !==============================================================================!
    !                           P O T _ A L L O C A T E                            !
    !==============================================================================!
    ! Subroutine for allocating memory for holding potential data type             !
    !------------------------------------------------------------------------------!
    ! Arguments:                                                                   !
    !           dummy_pot,         intent :: inout                                 !
    !           struct,            intent :: in                                    !
    !------------------------------------------------------------------------------!
    ! Author:   Z. Hawkhead  03/12/2020                                            !
    !==============================================================================!
    implicit none
    type(potential),intent(inout) :: dummy_pot
    type(structure),intent(in) :: struct
    integer  :: n
    integer :: stat,ni



    call trace_entry("pot_allocate")
    ! The potential array will be allocated on all processes and then zeroed
    if (allocated(dummy_pot%pot_array)) call io_errors("Error in POT: potential dummy_pot already allocated")
    n=struct%n_bodies

    allocate(dummy_pot%pot_array(1:n,1:3),stat=stat)
    dummy_pot%pot_array(1:n,1:3)=0_dp
    if (stat.ne.0) call io_errors("Error in POT: allocate pot_array")
    !allocate(dummy_pot%pot_dir(1:n,1:3),stat=stat)
    !dummy_pot%pot_dir(1:n,:)=0_dp
    if (stat.ne.0) call io_errors("Error in POT: allocate pot_dir")

    dummy_pot%is_allocated=.true.

    call trace_exit("pot_allocate")
    return
  end subroutine pot_allocate

  subroutine pot_calculate(dummy_pot,struct)
    !==============================================================================!
    !                          P O T _ C A L C U L A T E                           !
    !==============================================================================!
    ! Subroutine for calculating the gravitational potential and total energy of   !
    ! a system of objects                                                          !
    !------------------------------------------------------------------------------!
    ! Arguments:                                                                   !
    !           dummy_pot,         intent :: inout                                 !
    !           struct,            intent :: inout                                 !
    !------------------------------------------------------------------------------!
    ! Author:   Z. Hawkhead  03/12/2020                                            !
    !==============================================================================!
    implicit none
    type(potential),intent(inout) :: dummy_pot
    type(structure),intent(inout) :: struct
    real(dp),dimension(:,:,:),allocatable :: U_sum
    real(dp),dimension(:,:),allocatable :: E_sum
    real(dp),dimension(:,:,:),allocatable :: dir_sum
    integer   :: ni,u,i
    call trace_entry("pot_calculate")


    if (.not.dummy_pot%is_allocated) call io_errors("Error in POT: dummy_pot is not allocated")
    allocate(U_sum(1:struct%n_bodies,1:struct%n_bodies,1:3))
    allocate(E_sum(1:struct%n_bodies,1:struct%n_bodies))
    allocate(dir_sum(1:struct%n_bodies,1:struct%n_bodies,1:3))
    U_sum(:,:,:)=0_dp
    E_sum(:,:)=0_dp
    dir_sum(:,:,:)=0_dp
    !dummy_pot%pot_dir(:,:)=0_dp
    dummy_pot%pot_array(:,:)=0_dp
    !U_sum(:,3)=-4.1154282512119600E-081_dp
    do ni  = comms_scheme_array(rank,1),comms_scheme_array(rank,2)
       do u = comms_scheme_array(rank,3),comms_scheme_array(rank,4)
          !print*,u
          if (ni.eq.u)then
             cycle

          else

             U_sum(ni,u,:)=-G_AU*struct%masses(ni)*struct%masses(u) * &
                  & (struct%positions(ni,:)-struct%positions(u,:)) &
                  & / (current_params%epsilon+sqrt(sum((struct%positions(ni,:)-struct%positions(u,:))**2))**3)
             !print*,ni,u,U_sum(ni,u)


             E_sum(ni,u)=-G*struct%masses(ni)*struct%masses(u)*M_sol**2 &
                  & / (AU*(current_params%epsilon+abs(sum(struct%positions(ni,:)-struct%positions(u,:)))))&
                  & + 0.5_dp * struct%masses(ni)*M_sol &
                  & * (sum((AU*struct%init_velocity(ni,:)/days)**2))
             !dir_sum(ni,u,:)=struct%positions(ni,:)-struct%positions(u,:)
          end if

       end do
    end do
    ! Now we have to reduce and bcast
    !if (nprocs.gt.1)then 
       do ni=1,struct%n_bodies
          call COMMS_REDUCE(sum(U_sum(ni,:,1)),dummy_pot%pot_array(ni,1),3,"MPI_SUM")
          call COMMS_REDUCE(sum(U_sum(ni,:,2)),dummy_pot%pot_array(ni,2),3,"MPI_SUM")
          call COMMS_REDUCE(sum(U_sum(ni,:,3)),dummy_pot%pot_array(ni,3),3,"MPI_SUM")
          call COMMS_REDUCE(sum(E_sum(ni,:)),struct%tot_energy,1,"MPI_SUM")
          !call comms_reduce_double(sum(dir_sum(ni,:,1)),dummy_pot%pot_dir(ni,1),1,"MPI_SUM")
          !call comms_reduce_double(sum(dir_sum(ni,:,2)),dummy_pot%pot_dir(ni,2),1,"MPI_SUM")
          !call comms_reduce_double(sum(dir_sum(ni,:,3)),dummy_pot%pot_dir(ni,3),1,"MPI_SUM")
       end do


       !do  ni=1,struct%n_bodies
       !   dummy_pot%pot_dir(ni,:)=dummy_pot%pot_dir(ni,:)/sqrt(sum(dummy_pot%pot_dir(ni,:)**2))
       !
       !end do


       call comms_bcast(struct%tot_energy,1)
    !end if
    !call comms_bcast_double_array(dummy_pot%pot_array,size(dummy_pot%pot_array))
    !call comms_bcast_double_array(dummy_pot%pot_dir,size(dummy_pot%pot_dir))


    !print*,dummy_pot%pot_array(:)
    !do  ni=1,struct%n_bodies
    !   print*,dummy_pot%pot_dir(ni,:)
    !end do
    call trace_exit("pot_calculate")


    return
  end subroutine pot_calculate





end module pot



