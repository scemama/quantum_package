
program mrcc_sto
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  call dress_zmq()
end


! BEGIN_PROVIDER [ double precision, hij_cache_, (N_det,Nproc) ]
!&BEGIN_PROVIDER [ double precision, sij_cache_, (N_det,Nproc) ]
 BEGIN_PROVIDER [ double precision, fock_diag_tmp_, (2,mo_tot_num+1,Nproc) ]
&BEGIN_PROVIDER [ integer, current_generator_, (Nproc) ]
  implicit none
!  allocate(fock_diag_tmp(2,mo_tot_num+1))
  current_generator_(:) = 0
 END_PROVIDER



subroutine dress_with_alpha_buffer(delta_ij_loc, i_gen, minilist, det_minilist, n_minilist, alpha, iproc)
 use bitmasks
 implicit none
  BEGIN_DOC
  !delta_ij_loc(:,:,1) : dressing column for H
  !delta_ij_loc(:,:,2) : dressing column for S2
  !minilist : indices of determinants connected to alpha ( in psi_det_sorted )
  !n_minilist : size of minilist
  !alpha : alpha determinant
  END_DOC
  integer(bit_kind), intent(in)   :: alpha(N_int,2), det_minilist(N_int, 2, n_minilist)
  integer,intent(in)              :: minilist(n_minilist), n_minilist, iproc, i_gen
  double precision, intent(inout) :: delta_ij_loc(N_states,N_det,2)
  double precision :: hii, hij, sij, delta_e
  double precision, external :: diag_H_mat_elem_fock
  integer                        :: i,j,k,l,m, l_sd
  

  if(current_generator_(iproc) /= i_gen) then
    current_generator_(iproc) = i_gen
    call build_fock_tmp(fock_diag_tmp_(1,1,iproc),psi_det_generators(1,1,i_gen),N_int)
  end if
  !return
  hii = diag_H_mat_elem_fock(psi_det_generators(1,1,i_gen),alpha,fock_diag_tmp_(1,1,iproc),N_int)
  

  do l_sd=1,n_minilist
    call i_h_j_s2(alpha,det_minilist(1,1,l_sd),N_int,hij, sij)
    do i=1,N_states
      delta_ij_loc(i, minilist(l_sd), 1) += hij / hii * psi_coef(minilist(l_sd), i)
    end do
  end do
end subroutine






