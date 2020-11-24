program n_body
  use iso_fortran_env
  use comms
  use trace , only : trace_init,trace_entry,trace_exit,trace_finalise
  use io , only : io_initialise,current_params,io_write_results,io_write_params,io_dryrun,current_structure,io_cart_to_radial,&
       & io_radial_to_cart,dp
  use pot, only : pot_allocate,G_pot, pot_calculate

  implicit none
  real(dp), dimension(1,1:3) :: pos
  integer ::  i
  call trace_init()
  call trace_entry("n_body")

  call comms_init()

  call io_initialise()

  ! After the io_initialise call, the number of objects has been settled we can now work out the parallelism 
  call comms_scheme(current_structure%n_bodies)


  if (on_root_node)then
     call io_write_params()
     if (current_params%dry_run) call io_dryrun()
     !print*,current_structure%init_velocity
  end if

  call pot_allocate(G_pot,current_structure%n_bodies)
  !initial pot calculate
  call pot_calculate(G_pot,current_structure)



  call trace_exit("n_body")
  call trace_finalise(current_params%debuging,rank)
  call comms_finalise()


end program n_body
