
program mrcc_sto
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  call dress_zmq()
end


 BEGIN_PROVIDER [ double precision, fock_diag_tmp_, (2,mo_tot_num+1,Nproc) ]
&BEGIN_PROVIDER [ integer, current_generator_, (Nproc) ]
  implicit none
  current_generator_(:) = 0
 END_PROVIDER



subroutine dress_with_alpha_buffer(Nstates,Ndet,Nint,delta_ij_loc, i_gen, minilist, det_minilist, n_minilist, alpha, iproc)
 use bitmasks
 implicit none
  BEGIN_DOC
  !delta_ij_loc(:,:,1) : dressing column for H
  !delta_ij_loc(:,:,2) : dressing column for S2
  !i_gen : generator index in psi_det_generators
  !minilist : indices of determinants connected to alpha ( in psi_det_sorted )
  !n_minilist : size of minilist
  !alpha : alpha determinant
  END_DOC
  integer, intent(in)             :: Nint, Ndet, Nstates, n_minilist, iproc, i_gen
  integer(bit_kind), intent(in)   :: alpha(Nint,2), det_minilist(Nint, 2, n_minilist)
  integer,intent(in)              :: minilist(n_minilist)
  double precision, intent(inout) :: delta_ij_loc(Nstates,N_det,2)
  double precision :: hii, hij, sij, delta_e
  double precision, external :: diag_H_mat_elem_fock
  integer                        :: i,j,k,l,m, l_sd
  double precision, save :: tot = 0d0
  double precision :: de(N_states), val, tmp
  

  if(current_generator_(iproc) /= i_gen) then
    current_generator_(iproc) = i_gen
    call build_fock_tmp(fock_diag_tmp_(1,1,iproc),psi_det_generators(1,1,i_gen),N_int)
  end if

  hii = diag_H_mat_elem_fock(psi_det_generators(1,1,i_gen),alpha,fock_diag_tmp_(1,1,iproc),N_int)
  do i=1,N_states
    de(i) = (E0_denominator(i) - hii)
  end do
  
  do i=1,N_states
    val = 0D0
    do l_sd=1,n_minilist
      call i_h_j_s2(alpha,det_minilist(1,1,l_sd),N_int,hij, sij)
      val += hij
    end do
    val = 2d0 * val
    tmp = dsqrt(de(i)**2 + val**2)
    if(de(i) < 0d0) tmp = -tmp
    delta_ij_loc(i, minilist(l_sd), 1) += 0.5d0 * (tmp - de(i)) ! * psi_coef(minilist(l_sd), i)
  end do
end subroutine





BEGIN_PROVIDER [ logical, initialize_E0_denominator ]
    implicit none
    BEGIN_DOC
    ! If true, initialize pt2_E0_denominator
    END_DOC
    initialize_E0_denominator = .True.
END_PROVIDER
   
BEGIN_PROVIDER [ double precision, E0_denominator, (N_states) ]
  implicit none
  BEGIN_DOC
  ! E0 in the denominator of the PT2
  END_DOC
  if (initialize_E0_denominator) then
   E0_denominator(1:N_states) = psi_energy(1:N_states)
 ! call ezfio_get_full_ci_zmq_energy(pt2_E0_denominator(1))
 ! pt2_E0_denominator(1) -= nuclear_repulsion
 ! pt2_E0_denominator(1:N_states) = HF_energy - nuclear_repulsion
 ! pt2_E0_denominator(1:N_states) = barycentric_electronic_energy(1:N_states)
  else
    E0_denominator = -huge(1.d0)
  endif
END_PROVIDER
