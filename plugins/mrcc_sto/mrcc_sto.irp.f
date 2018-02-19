
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

 BEGIN_PROVIDER [ integer, idx_non_ref_from_sorted, (N_det) ]
&BEGIN_PROVIDER [ integer, psi_from_sorted, (N_det) ]
  implicit none
  integer :: i,inpsisor
  
  idx_non_ref_from_sorted = 0
  psi_from_sorted = 0

  do i=1,N_det
    psi_from_sorted(psi_det_sorted_order(i)) = i
    inpsisor = psi_det_sorted_order(i)
    if(inpsisor <= 0) stop "idx_non_ref_from_sorted"
    idx_non_ref_from_sorted(inpsisor) = idx_non_ref_rev(i)
  end do
END_PROVIDER

subroutine dress_with_alpha_buffer(delta_ij_loc, minilist, n_minilist, alpha)
 use bitmasks
 implicit none

  integer(bit_kind), intent(in)   :: alpha(N_int,2)
  integer,intent(in)              :: minilist(n_minilist), n_minilist
  double precision, intent(inout) :: delta_ij_loc(N_states,N_det,2)


  integer                        :: i,j,k,l,m
  integer                        :: degree1, degree2, degree
  integer, allocatable           :: idx_alpha(:)

  double precision               :: hIk, hla, hIl, sla, dIk(N_states), dka(N_states), dIa(N_states), hka
  double precision, allocatable  :: dIa_hla(:,:), dIa_sla(:,:)
  double precision               :: phase, phase2
  double precision               :: ci_inv(N_states)
  integer                        :: exc(0:2,2,2)
  integer                        :: h1,h2,p1,p2,s1,s2
  integer(bit_kind)              :: tmp_det(N_int,2)
  integer                        :: i_state, k_sd, l_sd, m_sd, ll_sd, i_I
  double precision, allocatable  :: hij_cache(:), sij_cache(:)
  double precision :: Delta_E_inv(N_states)
  double precision :: sdress, hdress
  double precision :: c0(N_states)
  logical :: ok
  

  if (perturbative_triples) then
    PROVIDE one_anhil fock_virt_total fock_core_inactive_total one_creat
  endif
  allocate(hij_cache(N_det), sij_cache(N_det)) 
  allocate (dIa_hla(N_states,N_det), dIa_sla(N_states,N_det))
  allocate (idx_alpha(0:n_minilist))
  
  do i_state=1,N_states
    c0(i_state) = 1.d0/psi_coef(dressed_column_idx(i_state),i_state)
  enddo
  ll_sd = 0
  do l_sd=1,n_minilist
    ok = .true.
    k_sd = minilist(l_sd)
    !if(idx_non_ref_rev(k_sd) == 0) cycle
    
    do i_I=1,N_det_ref
      call get_excitation_degree(psi_det_sorted(1,1,k_sd),psi_ref(1,1,i_I),degree1,N_int)
      if(degree1 == 0) then
        ok = .false.
        exit
      end if
    end do
    
    if( xor(ok, idx_non_ref_from_sorted(k_sd) > 0)) stop "BUGUE"
    if(ok) then
      ll_sd += 1
      idx_alpha(ll_sd) = k_sd
!     call i_h_j(alpha,psi_non_ref(1,1,idx_alpha(l_sd)),N_int,hij_cache(k_sd))
!     call get_s2(alpha,psi_non_ref(1,1,idx_alpha(l_sd)),N_int,sij_cache(k_sd))
     call i_h_j(alpha,psi_det_sorted(1,1,k_sd),N_int,hij_cache(k_sd))
     call get_s2(alpha,psi_det_sorted(1,1,k_sd),N_int,sij_cache(k_sd))
    end if
  enddo
  
  idx_alpha(0) = ll_sd
  

  do i_I=1,N_det_ref
    call get_excitation_degree(alpha,psi_ref(1,1,i_I),degree1,N_int)
    if (degree1 > 4) then
      cycle
    endif
    
    do i_state=1,N_states
      dIa(i_state) = 0.d0
    enddo

    do k_sd=1,idx_alpha(0)
      !print *, "idx ", k_sd
      !print *, "idx2", idx_alpha(k_sd)
      !print *, "ref
      !if(idx_non_ref_rev(idx_alpha(k_sd)) == 0) cycle
      call get_excitation_degree(psi_ref(1,1,i_I),psi_det_sorted(1,1,idx_alpha(k_sd)),degree,N_int)
!      print *, "diden"
      if (degree > 2) then
        cycle
      endif
      
      call get_excitation(psi_det_sorted(1,1,idx_alpha(k_sd)),alpha,exc,degree2,phase,N_int)
      !print *, "DEG", degree2
      call decode_exc(exc,degree2,h1,p1,h2,p2,s1,s2)
      do k=1,N_int
        tmp_det(k,1) = psi_ref(k,1,i_I)
        tmp_det(k,2) = psi_ref(k,2,i_I)
      enddo
      call apply_excitation(psi_ref(1,1,i_I), exc, tmp_det, ok, N_int)
      
      do i_state=1,N_states
        dIK(i_state) = dij(i_I, idx_non_ref_from_sorted(idx_alpha(k_sd)), i_state)
      enddo
      
      ! <I| \l/ |alpha>
      do i_state=1,N_states
        dka(i_state) = 0.d0
      enddo

      if (ok) then
        do l_sd=k_sd+1,idx_alpha(0)
          call get_excitation_degree(tmp_det,psi_det_sorted(1,1,idx_alpha(l_sd)),degree,N_int)
          if (degree == 0) then
            call get_excitation(psi_ref(1,1,i_I),psi_det_sorted(1,1,idx_alpha(l_sd)),exc,degree,phase2,N_int)
            do i_state=1,N_states
              dka(i_state) = dij(i_I, idx_non_ref_from_sorted(idx_alpha(l_sd)), i_state) * phase * phase2
            enddo
            exit
          endif
        enddo
      else if (perturbative_triples) then
          hka = hij_cache(idx_alpha(k_sd))
          if (dabs(hka) > 1.d-12) then
            call get_delta_e_dyall_general_mp(psi_ref(1,1,i_I),alpha,Delta_E_inv)

            do i_state=1,N_states
              ASSERT (Delta_E_inv(i_state) < 0.d0)
              dka(i_state) = hka / Delta_E_inv(i_state)
            enddo
          endif
      endif

      if (perturbative_triples.and. (degree2 == 1) ) then
          call i_h_j(psi_ref(1,1,i_I),tmp_det,N_int,hka)
          hka = hij_cache(idx_alpha(k_sd)) - hka
          if (dabs(hka) > 1.d-12) then
            call get_delta_e_dyall_general_mp(psi_ref(1,1,i_I),alpha,Delta_E_inv)
            do i_state=1,N_states
              ASSERT (Delta_E_inv(i_state) < 0.d0)
              dka(i_state) = hka / Delta_E_inv(i_state)
            enddo
          endif
      endif
      do i_state=1,N_states
        dIa(i_state) = dIa(i_state) + dIk(i_state) * dka(i_state)
      enddo
    enddo
    
    do i_state=1,N_states
      ci_inv(i_state) = psi_ref_coef_inv(i_I,i_state)
    enddo
    do l_sd=1,idx_alpha(0)
      k_sd = idx_alpha(l_sd)
      hla = hij_cache(k_sd)
      sla = sij_cache(k_sd)
      do i_state=1,N_states
        dIa_hla(i_state,k_sd) = dIa(i_state) * hla
        dIa_sla(i_state,k_sd) = dIa(i_state) * sla
      enddo
    enddo
    do i_state=1,N_states
      do l_sd=1,idx_alpha(0)
        !print *, "DRES"
        !print *, i_state, idx_alpha(l_sd)
        k_sd = idx_alpha(l_sd)
        m_sd = psi_from_sorted(k_sd)
        if(psi_det(1,1,m_sd) /= psi_det_sorted(1,1,k_sd)) stop "psi_from_sorted foireous"
        hdress = dIa_hla(i_state,k_sd) * psi_ref_coef(i_I,i_state) * c0(i_state)
        sdress = dIa_sla(i_state,k_sd) * psi_ref_coef(i_I,i_state) * c0(i_state)
        !$OMP ATOMIC
        delta_ij_loc(i_state,m_sd,1) += hdress
        !$OMP ATOMIC
        delta_ij_loc(i_state,m_sd,2) += sdress
        !print *, "ENDRES" 
      enddo
    enddo
  enddo
end subroutine





!! TESTS MINILIST
subroutine test_minilist(minilist, n_minilist, alpha)
  use bitmasks
  implicit none
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
end subroutine


