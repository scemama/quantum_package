subroutine dress_slave
  implicit none
  BEGIN_DOC
! Helper subroutine to compute the dress in distributed mode.
  END_DOC
  read_wf = .False.
  distributed_davidson = .False.
  SOFT_TOUCH read_wf distributed_davidson
  call provide_everything
  call switch_qp_run_to_master
  call run_wf
end

subroutine provide_everything
  PROVIDE H_apply_buffer_allocated mo_bielec_integrals_in_map psi_det_generators psi_coef_generators psi_det_sorted_bit psi_selectors n_det_generators n_states generators_bitmask zmq_context
end

subroutine run_wf
  use f77_zmq
  implicit none

  integer(ZMQ_PTR), external :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR) :: zmq_to_qp_run_socket
  double precision :: energy(N_states_diag)
  character*(64) :: states(1)
  integer :: rc, i
  
  call provide_everything
  
  zmq_context = f77_zmq_ctx_new ()
  states(1) = 'dress'

  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()

  do

    call wait_for_states(states,zmq_state,1)

    if(trim(zmq_state) == 'Stopped') then

      exit

    else if (trim(zmq_state) == 'dress') then

      ! Selection
      ! ---------

      print *,  'dress'
      call zmq_get_psi(zmq_to_qp_run_socket,1,energy,N_states)

      PROVIDE psi_bilinear_matrix_columns_loc psi_det_alpha_unique psi_det_beta_unique
      PROVIDE psi_bilinear_matrix_rows psi_det_sorted_order psi_bilinear_matrix_order
      PROVIDE psi_bilinear_matrix_transp_rows_loc psi_bilinear_matrix_transp_columns
      PROVIDE psi_bilinear_matrix_transp_order

      !$OMP PARALLEL PRIVATE(i)
      i = omp_get_thread_num()
      call dress_slave_tcp(i, energy)
      !$OMP END PARALLEL
      print *,  'dress done'

    endif

  end do
end

subroutine dress_slave_tcp(i,energy)
  implicit none
  double precision, intent(in) :: energy(N_states_diag)
  integer, intent(in)            :: i
  logical :: lstop
  lstop = .False.
  call run_dress_slave(0,i,energy,lstop)
end

