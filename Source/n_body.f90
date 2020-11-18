program n_body
  use iso_fortran_env
  use comms
  use trace , only : trace_init,trace_entry,trace_exit,trace_finalise
  use io , only : io_initialise,current_params,io_write_results,io_write_params,io_dryrun


  call trace_init()
  call trace_entry("n_body")

  call comms_init()

  call io_initialise()
  if (on_root_node)then
     call io_write_params()



     if (current_params%dry_run) call io_dryrun()

     
     call io_write_results()
  end if
  call trace_exit("n_body")
  call trace_finalise(current_params%debuging,rank)
  call comms_finalise()


end program n_body
