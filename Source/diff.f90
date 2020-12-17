module diff
  use comms
  use io, only : current_params,current_structure,io_errors,structure,dp,&
       & io_write_config,stdout,io_flush
  use trace, only : trace_entry,trace_exit
  use pot, only : pot_allocate, pot_calculate,potential
  implicit none
  real(dp),private :: perc_prog=10.0
contains

  subroutine diff_solver(struct)
    implicit none
    type(structure), intent(inout) :: struct
    type(potential)   :: G_pot
    integer :: ni,u,i,j,k
    integer :: t
    integer         :: width=97

    call trace_entry("diff_solver")

    ! Write out the time steps

    write(stdout,*)
    write(stdout,*)"+",repeat("=",width-3),"+ <--DIFF"
    write(stdout,'(1x,A,T40,a,T97,a)')'|',"DIFFERENTIAL SOLVER","| <--DIFF"
    write(stdout,*)"+",repeat("=",width-3),"+ <--DIFF"
    write(stdout,'(" |",T10,a,T75,a,T97,"| <--DIFF")')"Percentage Complete (%) ","Time (s)"
    write(stdout,*)"+",repeat("=",width-3),"+ <--DIFF"



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
             !print*,"Positions",rank,ni
             call comms_send(struct%positions(ni,:),size(struct%positions(ni,:)),0,ni)
             !print*,"Velocity",rank,ni
             call comms_send(struct%init_velocity(ni,:),size(struct%positions(ni,:)),0,ni+struct%n_bodies)
          end if
       end  do

       ! Now we do the recv
       if (on_root_node.and.nprocs.gt.1)then 
          do i = 1,nprocs-1
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
       if (on_root_node)call diff_progress(struct,current_params%time_step)
       ! If we are writing the file, do it here
       if (current_params%write_config)call io_write_config(struct)

    end do

    write(stdout,*)"+",repeat("=",94),"+ <--DIFF"

    call trace_exit("diff_solver")
    return
  end  subroutine diff_solver


  subroutine diff_progress(struct,time_step)
    use trace, only :trace_entry, trace_exit,global_start
    type(structure),intent(in) :: struct
    real(dp),intent(in) :: time_step
    real(dp) :: exact,ctime

    exact=current_params%calc_len*perc_prog/100.0_dp
    !print*,exact, struct%sys_time,global_start
    if (struct%sys_time-time_step.le.exact.and.struct%sys_time.gt.exact)then
       ctime=comms_wtime()
       call cpu_time(ctime)
       write(stdout,'(" |",T10,f5.1,1x,"%",T70,f14.6,T97,"| <--DIFF")')perc_prog,ctime-global_start
       perc_prog=perc_prog+10
       call io_flush(stdout)
    end if

    
    return
  end subroutine diff_progress




  
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
