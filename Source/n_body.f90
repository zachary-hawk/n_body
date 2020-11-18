program n_body
  use iso_fortran_env
  use comms
  use trace , only : trace_init,trace_entry,trace_exit,trace_finalise



  call trace_init()
  call trace_entry("n_body")

  call comms_init()


  call trace_exit("n_body")
  call trace_finalise(.true.,rank)
  call comms_finalise()


end program n_body
