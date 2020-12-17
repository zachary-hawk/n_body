program n_body
  use comms
  use trace , only : trace_init,trace_entry,trace_exit,trace_finalise
  use io , only : io_initialise,current_params,io_write_results,io_write_params,io_dryrun,current_structure,io_cart_to_radial,&
       & io_radial_to_cart,dp,io_write_fmt_structure
  !use pot, only : pot_allocate pot_calculate
  use diff, only : diff_solver
  implicit none
  real(dp), dimension(1,1:3) :: pos
  integer ::  i



  call trace_init()
  call trace_entry("n_body")

  call comms_init()
  
  call io_initialise()

  ! After the io_initialise call, the number of objects has been settled we can now work out the parallelism 

  call comms_scheme(current_structure%n_bodies)
  print*,comms_scheme_array
  if (on_root_node)then
     call io_write_params()
     if (current_params%dry_run) call io_dryrun()
     !print*,current_structure%init_velocity
  end if

  print*,"before diff solver"


  ! do a test  call to diff
  call diff_solver(current_structure)
  print*,"After diff"
  
  

  
  if (on_root_node) then
     call io_write_results()
     if (current_params%write_fmt)&
          & call io_write_fmt_structure(current_structure)
  end if



  call trace_exit("n_body")
  print*,"trace_finalise"
  call trace_finalise(current_params%debuging,rank)
  print*,"comms_finalise"
  call comms_finalise()
  print*,"Fin"

end program n_body
