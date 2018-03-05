BEGIN_PROVIDER [ double precision, CI_energy_dressed, (N_states_diag) ]
  implicit none
  BEGIN_DOC
  ! N_states lowest eigenvalues of the CI matrix
  END_DOC
  
  integer                        :: j
  character*(8)                  :: st
  call write_time(6)
  do j=1,min(N_det,N_states_diag)
    CI_energy_dressed(j) = CI_electronic_energy_dressed(j) + nuclear_repulsion
  enddo
  do j=1,min(N_det,N_states)
    write(st,'(I4)') j
    call write_double(6,CI_energy_dressed(j),'Energy of state '//trim(st))
    call write_double(6,CI_eigenvectors_s2_dressed(j),'S^2 of state '//trim(st))
  enddo

END_PROVIDER

 BEGIN_PROVIDER [ double precision, CI_electronic_energy_dressed, (N_states_diag) ]
&BEGIN_PROVIDER [ double precision, CI_eigenvectors_dressed, (N_det,N_states_diag) ]
&BEGIN_PROVIDER [ double precision, CI_eigenvectors_s2_dressed, (N_states_diag) ]
  BEGIN_DOC
  ! Eigenvectors/values of the CI matrix
  END_DOC
  implicit none
  double precision               :: ovrlp,u_dot_v
  integer                        :: i_good_state
  integer, allocatable           :: index_good_state_array(:)
  logical, allocatable           :: good_state_array(:)
  double precision, allocatable  :: s2_values_tmp(:)
  integer                        :: i_other_state
  double precision, allocatable  :: eigenvectors(:,:), eigenvectors_s2(:,:), eigenvalues(:)
  integer                        :: i_state
  double precision               :: e_0
  integer                        :: i,j,k,mrcc_state
  double precision, allocatable  :: s2_eigvalues(:)
  double precision, allocatable  :: e_array(:)
  integer, allocatable           :: iorder(:)

  PROVIDE threshold_davidson nthreads_davidson
  ! Guess values for the "N_states" states of the CI_eigenvectors_dressed
  do j=1,min(N_states,N_det)
    do i=1,N_det
      CI_eigenvectors_dressed(i,j) = psi_coef(i,j)
    enddo
  enddo

  do j=min(N_states,N_det)+1,N_states_diag
    do i=1,N_det
      CI_eigenvectors_dressed(i,j) = 0.d0
    enddo
  enddo

  if (diag_algorithm == "Davidson") then
    
    allocate (eigenvectors(size(CI_eigenvectors_dressed,1),size(CI_eigenvectors_dressed,2)),&
        eigenvectors_s2(size(CI_eigenvectors_dressed,1),size(CI_eigenvectors_dressed,2)),&
        eigenvalues(size(CI_electronic_energy_dressed,1)))
    do j=1,min(N_states,N_det)
      do i=1,N_det
        eigenvectors(i,j) = psi_coef(i,j)
      enddo
    enddo
    do mrcc_state=1,N_states
      do j=mrcc_state,min(N_states,N_det)
        do i=1,N_det
          eigenvectors(i,j) = psi_coef(i,j)
        enddo
      enddo
      call davidson_diag_HS2(psi_det,eigenvectors, eigenvectors_s2,    &
          size(eigenvectors,1),                                        &
          eigenvalues,N_det,min(N_det,N_states),min(N_det,N_states_diag),N_int,&
          mrcc_state)
      CI_eigenvectors_dressed(1:N_det,mrcc_state) = eigenvectors(1:N_det,mrcc_state)
      CI_electronic_energy_dressed(mrcc_state) = eigenvalues(mrcc_state)
    enddo
    do k=N_states+1,N_states_diag
      CI_eigenvectors_dressed(1:N_det,k) = eigenvectors(1:N_det,k)
      CI_electronic_energy_dressed(k) = eigenvalues(k)
    enddo
    call u_0_S2_u_0(CI_eigenvectors_s2_dressed,CI_eigenvectors_dressed,N_det,psi_det,N_int,&
        N_states_diag,size(CI_eigenvectors_dressed,1))
    
    deallocate (eigenvectors,eigenvalues)
    
    
  else if (diag_algorithm == "Lapack") then
    
    allocate (eigenvectors(size(H_matrix_dressed,1),N_det))
    allocate (eigenvalues(N_det))

    do j=1,min(N_states,N_det)
      do i=1,N_det
        eigenvectors(i,j) = psi_coef(i,j)
      enddo
    enddo
    do mrcc_state=1,N_states
      do j=mrcc_state,min(N_states,N_det)
        do i=1,N_det
          eigenvectors(i,j) = psi_coef(i,j)
        enddo
      enddo

    call lapack_diag(eigenvalues,eigenvectors,                         &
        H_matrix_dressed(1,1,mrcc_state),size(H_matrix_dressed,1),N_det)
      CI_eigenvectors_dressed(1:N_det,mrcc_state) = eigenvectors(1:N_det,mrcc_state)
      CI_electronic_energy_dressed(mrcc_state) = eigenvalues(mrcc_state)
    enddo
    do k=N_states+1,N_states_diag
      CI_eigenvectors_dressed(1:N_det,k) = eigenvectors(1:N_det,k)
      CI_electronic_energy_dressed(k) = eigenvalues(k)
    enddo
    call u_0_S2_u_0(CI_eigenvectors_s2_dressed,CI_eigenvectors_dressed,N_det,psi_det,N_int,&
        N_states_diag,size(CI_eigenvectors_dressed,1))

    deallocate(eigenvectors,eigenvalues)
  endif

END_PROVIDER

subroutine diagonalize_CI_dressed
  implicit none
  BEGIN_DOC
!  Replace the coefficients of the CI states by the coefficients of the 
!  eigenstates of the CI matrix
  END_DOC
  integer :: i,j
  do j=1,N_states
    do i=1,N_det
      psi_coef(i,j) = CI_eigenvectors_dressed(i,j)
    enddo
  enddo
  SOFT_TOUCH psi_coef 
end



BEGIN_PROVIDER [ double precision, h_matrix_dressed, (N_det,N_det,N_states) ]
 implicit none
 BEGIN_DOC
 ! Dressed H with Delta_ij
 END_DOC
 integer                        :: i, j, ii,jj, dressing_state
 do dressing_state = 1,N_states
   do j=1,N_det
     do i=1,N_det
       h_matrix_dressed(i,j,dressing_state) = h_matrix_all_dets(i,j) 
     enddo
   enddo
   i = dressed_column_idx(dressing_state)
   do j = 1, N_det
     h_matrix_dressed(i,j,dressing_state) += dressing_column_h(j,dressing_state)
     h_matrix_dressed(j,i,dressing_state) += dressing_column_h(j,dressing_state)
   enddo
   h_matrix_dressed(i,i,dressing_state) -= dressing_column_h(i,dressing_state)
 enddo
END_PROVIDER

