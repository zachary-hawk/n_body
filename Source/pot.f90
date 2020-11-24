module pot
  use trace, only : trace_entry, trace_exit
  use io
  use comms
  implicit none
  

  type potential
     real(dp), allocatable,dimension(:) :: pot_array
     logical                            :: is_allocated = .false.
  end type potential

  type(potential),public,save :: G_pot

  
contains
  
  subroutine pot_allocate(dummy_pot,n)
    implicit none
    type(potential),intent(inout) :: dummy_pot
    integer,intent(in) :: n
    integer :: stat,ni
    
    call trace_entry("pot_allocate")
    ! The potential array will be allocated on all processes and then zeroed
    if (allocated(dummy_pot%pot_array)) call io_errors("Error in POT: potential dummy_pot already allocated")


    allocate(dummy_pot%pot_array(1:n),stat=stat)
    dummy_pot%pot_array(1:n)=0_dp
    if (stat.ne.0) call io_errors("Error in POT: allocate pot_array")

    dummy_pot%is_allocated=.true.
    
    call trace_exit("pot_allocate")
    return
  end subroutine pot_allocate

  subroutine pot_calculate(dummy_pot,struct)
    implicit none
    type(potential),intent(inout) :: dummy_pot
    type(structure),intent(inout) :: struct
    real(dp),dimension(:,:),allocatable :: U_sum
    real(dp),dimension(:,:),allocatable :: E_sum
    integer   :: ni,u
    call trace_entry("pot_calculate")


    if (.not.dummy_pot%is_allocated) call io_errors("Error in POT: dummy_pot is not allocated")
    allocate(U_sum(1:struct%n_bodies,1:struct%n_bodies))
    allocate(E_sum(1:struct%n_bodies,1:struct%n_bodies))
    U_sum(:,:)=0_dp
    E_sum(:,:)=0_dp
    do ni  = comms_scheme_array(rank,1),comms_scheme_array(rank,2)
       do u = comms_scheme_array(rank,3),comms_scheme_array(rank,4)
          if (ni.eq.u)then
             cycle
          else
             U_sum(ni,u)=-G_AU*struct%masses(ni)*struct%masses(u) &                  
                  & / (current_params%epsilon+sum(struct%positions(ni,:)-struct%positions(u,:))**2)
             E_sum(ni,u)=-G*struct%masses(ni)*struct%masses(u)*M_sol**2 &
                  & / (AU*(current_params%epsilon+abs(sum(struct%positions(ni,:)-struct%positions(u,:)))))&
                  & + 0.5_dp * struct%masses(ni)*M_sol &
                  & * (sum((AU*struct%init_velocity(ni,:)/days)**2))
          end if
       end do

    end do
    ! Now we have to reduce and bcast
    do ni=1,struct%n_bodies
       call COMMS_REDUCE_DOUBLE(sum(U_sum(ni,:)),dummy_pot%pot_array(ni),1,"MPI_SUM")
       call COMMS_REDUCE_DOUBLE(sum(E_sum(ni,:)),struct%tot_energy,1,"MPI_SUM")
    end do


    call comms_bcast(struct%tot_energy,1)
    call comms_bcast_double_array(dummy_pot%pot_array,size(dummy_pot%pot_array))

    call trace_exit("pot_calculate")
    return
  end subroutine pot_calculate
  


  

end module pot



