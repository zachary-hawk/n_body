module diff
  use comms
  use io, only : current_params,current_structure,io_errors,structure
  use trace, only : trace_entry,trace_exit
  use pot, only : pot_allocate, pot_calculate,potential
  implicit none

contains

  subroutine diff_euler(struct)
    implicit none
    type(structure), intent(inout) :: struct
    type(potential)   :: G_pot
    integer :: ni,u,i,j,k
    integer :: t
    call trace_entry("diff_euler")

    ! Allocate the potential
    call pot_allocate(G_pot,struct)

    ! Calculate the initial potential

    open(newunit=u,file="earth.dat",form="FORMATTED",status="unknown",access="stream")
    open(newunit=k,file="energy.dat",form="FORMATTED",status="unknown",access="stream")
    do t=1,10000
       call pot_calculate(G_pot,struct)
       !print*,G_pot%pot_dir(2,:), G_pot%pot_array(2)
       do ni  = comms_scheme_array(rank,1),comms_scheme_array(rank,2)
          ! update the position
          do i=1,3
             struct%positions(ni,i) = struct%positions(ni,i) + &
                  & current_params%time_step * struct%init_velocity(ni,i) 
          end  do
          ! update the Momentum
          do i=1,3
             struct%momentum(ni,i)=struct%momentum(ni,i) + &
                  & current_params%time_step * G_pot%pot_array(ni) * & ! this is the momentum... 
                  & G_pot%pot_dir(ni,i) !/ struct%masses(ni)
             struct%init_velocity(ni,i)=struct%momentum(ni,i)/struct%masses(ni)

          end do
       end  do
       write(u,*)(struct%positions(i,:),i=1,struct%n_bodies)
       !write(u,*) (G_pot%pot_array(i),i=1,struct%n_bodies)
       write(k,*) struct%tot_energy
       !print*,struct%tot_energy
       !print*,current_params%time_step * G_pot%pot_array(3) * & ! this is the momentum...                                                                                                                                                                                                                                                                            
       !     & G_pot%pot_dir(3,:)
       !do i = 1,struct%n_bodies
       !   print*,struct%labels(i),G_pot%pot_array(i)
       !end do
    end do
    print*,struct%labels(2),"vel",struct%init_velocity(2,:)
    print*,struct%labels(2),"pos",struct%positions(2,:)
    print*,struct%labels(2),"dir",G_pot%pot_dir(2,:)

    call trace_exit("diff_euler")
    return
  end  subroutine diff_euler


end module diff
