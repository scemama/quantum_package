program mrcc_sto
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  print *, "========================"
  print *, "========================"
  print *, "========================"
  print *, "MRCC_STO not implemented - acts as a unittest for dress_zmq"
  print *, "========================"
  print *, "========================"
  print *, "========================"
  call dress_zmq()
end



!! TESTS MINILIST
subroutine dress_with_alpha_buffer(delta_ij_loc, minilist, n_minilist, alpha)
  use bitmasks
  implicit none
  double precision,intent(inout)  :: delta_ij_loc(N_states,N_det,2) 
  integer, intent(in)            :: n_minilist
  integer(bit_kind),intent(in)  :: alpha(N_int, 2)
  integer, intent(in) :: minilist(n_minilist)
  integer :: a, i, deg
  integer :: refc(N_det), testc(N_det)
  
  refc = 0
  testc = 0
  do i=1,N_det
    call get_excitation_degree(psi_det_sorted(1,1,i), alpha, deg, N_int) 
    if(deg <= 2) refc(i) = refc(i) + 1
  end do
  do i=1,n_minilist
    call get_excitation_degree(psi_det_sorted(1,1,minilist(i)), alpha, deg, N_int) 
    if(deg <= 2) then
      testc(minilist(i)) += 1
    else
      stop "NON LINKED IN MINILIST"
    end if
  end do
  
  do i=1,N_det
    if(refc(i) /= testc(i)) then
      print *, "MINILIST FAIL ", sum(refc), sum(testc), n_minilist
      exit
    end if
  end do

  delta_ij_loc = 0d0
end subroutine


