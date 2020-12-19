program diag_and_save
 implicit none
 read_wf = .True.
 touch read_wf
 call routine
end

subroutine routine
 implicit none
 print*,'N_det = ',N_det
 call diagonalize_CI
 write(*,*)'Which state would you like to save ?'
 integer :: igood_state
 read(5,*)igood_state
 double precision, allocatable :: psi_coef_tmp(:)
 allocate(psi_coef_tmp(n_det))
 integer :: i
 do i = 1, N_det
  psi_coef_tmp(i) = psi_coef(i,igood_state)
 enddo
 call save_wavefunction_general(N_det,1,psi_det,n_det,psi_coef_tmp)
 deallocate(psi_coef_tmp)



end
