module diff
  use comms
  use io, only : current_params,current_structure,io_errors,structure,dp,&
       & io_write_config
  use trace, only : trace_entry,trace_exit
  use pot, only : pot_allocate, pot_calculate,potential
  implicit none

contains

  subroutine diff_solver(struct)
    implicit none
    type(structure), intent(inout) :: struct
    type(potential)   :: G_pot
    integer :: ni,u,i,j,k
    integer :: t
    call trace_entry("diff_solver")

    ! Allocate the potential
    call pot_allocate(G_pot,struct)

    ! Calculate the initial potential

    open(newunit=u,file="earth.dat",form="FORMATTED",status="unknown",access="stream")
    open(newunit=k,file="energy.dat",form="FORMATTED",status="unknown",access="stream")
    do while (struct%sys_time.lt.current_params%calc_len)
       call pot_calculate(G_pot,struct)
       call diff_time_step(G_pot)

       ! Select a method based on whats in the params file

       do ni  = comms_scheme_array(rank,1),comms_scheme_array(rank,2)






          select case(current_params%diff_method)
          case('euler')



             ! update the position
             do i=1,3
                struct%positions(ni,i) = struct%positions(ni,i) + &
                     & current_params%time_step * struct%init_velocity(ni,i) 
             end  do
             ! update the Momentum
             do i=1,3
                struct%momentum(ni,i)=struct%momentum(ni,i) + &
                     & current_params%time_step * G_pot%pot_array(ni,i)! * & ! this is the momentum... 
                     !& G_pot%pot_dir(ni,i) !/ struct%masses(ni)
                struct%init_velocity(ni,i)=struct%momentum(ni,i)/struct%masses(ni)

             end do
          case default
             call io_errors("Error in diff_method: Unknown method "//trim(current_params%diff_method))
          end select





          ! send and receive
          if (.not.on_root_node)then
             call comms_send(struct%positions(ni,:),size(struct%positions(ni,:)),0,ni)
             call comms_send(struct%init_velocity(ni,:),size(struct%positions(ni,:)),0,ni+struct%n_bodies)
          end if
       end  do

       ! Now we do the recv
       if (on_root_node)then 
          do i = 0,nprocs-1
             do ni=comms_scheme_array(i,1),comms_scheme_array(i,2)
                call comms_recv(struct%positions(ni,:),3,i,ni)
                call comms_recv(struct%init_velocity(ni,:),3,i,ni+struct%n_bodies)
             end do
          end do
          ! now bcast
          do ni=1,struct%n_bodies
             call comms_bcast(struct%positions(ni,:),3)
             call comms_bcast(struct%init_velocity(ni,:),3)
          end do
       end if

       ! Advance  the system clock
       struct%sys_time=struct%sys_time+current_params%time_step

       ! If we are writing the file, do it here
       if (current_params%write_config)call io_write_config(struct)

    end do


    call trace_exit("diff_solver")
    return
  end  subroutine diff_solver







  
  subroutine diff_time_step(pot)
    implicit none
    type(potential),intent(in)   :: pot
    real(dp)                     :: dt
    call trace_entry("diff_time_step")


    if (current_params%var_time.eq."fixed")then
       dt=current_params%time_step

    else if (current_params%var_time.eq."variable")then
       dt=1E-3_dp/maxval(abs(pot%pot_array(:,:)))
       if (dt.gt.0.1_dp) dt=0.1_dp

    end if

    ! reduce and bcast
    call comms_reduce(dt,current_params%time_step,1,"MPI_MIN")
    if (on_root_node) call comms_bcast(current_params%time_step,1)
    call trace_exit("diff_time_step")

    return
  end subroutine diff_time_step

  
end module diff
