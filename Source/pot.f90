module pot
  use trace, only : trace_entry, trace_exit
  use io
  use comms
  implicit none
  

  type potential
     real(dp), dimension(:),allocatable :: pot_array
     logical                            :: is_allocated = .false.
  end type potential

  
  subroutine pot_allocate(dummy_pot,n)
    implicit none
    type(potential),intent(inout) :: dummy_pot

    ! lets take the number of bodies


  end subroutine pot_allocate


end  module pot
