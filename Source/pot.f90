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
    real(dp),dimension(:,:),allocatable :: U_loc
    integer   :: ni,u,i
    call trace_entry("pot_calculate")


    ! First thing is to communicate
    ! send and receive

    !    print*,"Rank:",rank,"Before comm",struct%positions(1,:)
    do ni  = comms_scheme_array(rank,1),comms_scheme_array(rank,2)     
       if (.not.on_root_node)then
          call comms_send(struct%positions(ni,:),size(struct%positions(ni,:)),0,rank)
          !call comms_send(struct%init_velocity(ni,:),size(struct%positions(ni,:)),0,ni+struct%n_bodies)
       end if
    end  do


    ! Now we do the recv                                                                                                                                                             
    if (on_root_node)then
       do i = 1,nprocs-1
          do ni=comms_scheme_array(i,1),comms_scheme_array(i,2)
             !print*,struct%positions(ni,1:3)
             call comms_recv(struct%positions(ni,1:3),3,i,i)
             !call comms_recv(struct%init_velocity(ni,:),3,i,ni+struct%n_bodies)
          end do
       end do
    end if
    ! now bcast                                                                                                                                                                      

    
    do ni=1,struct%n_bodies      
       call comms_bcast(struct%positions(ni,1:3),3)
       !call comms_bcast(struct%positions(ni,:),2)
       !call comms_bcast(struct%positions(ni,:),2)
    end do
    print*,"Rank:",rank,struct%positions(2,1:3) 
    !    print*,"Rank:",rank,"after comm",struct%positions(1,:)

    !print*,"Rank:",rank,struct%labels(2),":",struct%positions(2,:)


    if (.not.dummy_pot%is_allocated) call io_errors("Error in POT: dummy_pot is not allocated")
    allocate(U_sum(1:struct%n_bodies,1:struct%n_bodies,1:3))
    allocate(E_sum(1:struct%n_bodies,1:struct%n_bodies))
    allocate(dir_sum(1:struct%n_bodies,1:struct%n_bodies,1:3))
    allocate(U_loc(1:struct%n_bodies,1:3))
    U_sum(:,:,:)=0_dp
    E_sum(:,:)=0_dp
    dir_sum(:,:,:)=0_dp
    U_loc(:,:)=0_dp
    dummy_pot%pot_array(:,:)=0_dp
    do ni  = comms_scheme_array(rank,1),comms_scheme_array(rank,2)
       do u = comms_scheme_array(rank,3),comms_scheme_array(rank,4)
          if (ni.eq.u)then
             cycle

          else

             U_sum(ni,u,:)=-G_AU*struct%masses(ni)*struct%masses(u) * &
                  & (struct%positions(ni,:)-struct%positions(u,:)) &
                  & / (current_params%epsilon+sqrt(sum((struct%positions(ni,:)-struct%positions(u,:))**2))**3)

             E_sum(ni,u)=-G*struct%masses(ni)*struct%masses(u)*M_sol**2 &
                  & / (AU*(current_params%epsilon+abs(sum(struct%positions(ni,:)-struct%positions(u,:)))))&
                  & + 0.5_dp * struct%masses(ni)*M_sol &
                  & * (sum((AU*struct%init_velocity(ni,:)/days)**2))
             !dir_sum(ni,u,:)=struct%positions(ni,:)-struct%positions(u,:)
          end if

       end do
    end do


    ! Now we have to reduce and bcast

    do ni=1,struct%n_bodies!comms_scheme_array(rank,1),comms_scheme_array(rank,2)

       U_loc(ni,1)=sum(U_sum(ni,:,1))
       U_loc(ni,2)=sum(U_sum(ni,:,2))
       U_loc(ni,3)=sum(U_sum(ni,:,3))
       !print*,"Rank:",rank,"ni:",ni,"SUM",U_loc(:,1)
       call comms_reduce(U_loc(ni,:),dummy_pot%pot_array(ni,:),3,"MPI_SUM")

       call comms_bcast(dummy_pot%pot_array(ni,1),1)
       call comms_bcast(dummy_pot%pot_array(ni,2),1)      
       call comms_bcast(dummy_pot%pot_array(ni,3),1)
       !print*,"Rank:",rank,"ni:",ni,"POT:",dummy_pot%pot_array(ni,1)
    end do
    call comms_bcast(struct%tot_energy,1)




    call trace_exit("pot_calculate")


    return
  end subroutine pot_calculate





end module pot



