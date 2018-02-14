
! BEGIN_PROVIDER [ logical, do_dress_with_alpha ]
!&BEGIN_PROVIDER [ logical, do_dress_with_alpha_buffer ]
!&BEGIN_PROVIDER [ logical, do_dress_with_generator ]
!  implicit none
!  do_dress_with_alpha =          .false.
!  do_dress_with_alpha_buffer =   .true.
!  do_dress_with_generator =      .false.
!END_PROVIDER

subroutine dress_with_alpha_buffer(delta_ij_loc, minilist, n_minilist, abuf, n_abuf)
  use bitmasks
  implicit none
  double precision,intent(inout) :: delta_ij_loc(N_states,N_det,2) 
  integer, intent(in)            :: n_minilist, n_abuf
  integer(bit_kind),intent(in)  :: abuf(N_int, 2, n_abuf), minilist(N_int, 2, n_minilist)
  integer :: a, i, nref, nobt, deg
  
  do a=1,n_abuf
    nref=0
    do i=1,N_det
      call get_excitation_degree(psi_det(1,1,i), abuf(1,1,a), deg, N_int) 
      if(deg <= 2) nref = nref + 1
    end do
    nobt=0
    do i=1,n_minilist
      call get_excitation_degree(minilist(1,1,i), abuf(1,1,a), deg, N_int) 
      if(deg <= 2) nobt = nobt + 1
    end do

    if(nref /= nobt) stop "foireous minilist"
  end do

  delta_ij_loc = 1d0
end subroutine


