
program mrcc_sto
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  call dress_zmq()
end

 BEGIN_PROVIDER [ integer, psi_from_sorted, (N_det) ]
&BEGIN_PROVIDER [ integer, idx_non_ref_from_sorted, (N_det) ]
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


 BEGIN_PROVIDER [ double precision, hij_cache_, (N_det,Nproc) ]
&BEGIN_PROVIDER [ double precision, sij_cache_, (N_det,Nproc) ]
&BEGIN_PROVIDER [ double precision, dIa_hla_, (N_states,N_det,Nproc) ]
&BEGIN_PROVIDER [ double precision, dIa_sla_, (N_states,N_det,Nproc) ]
&BEGIN_PROVIDER [ integer, idx_alpha_, (0:N_det,Nproc) ]
 BEGIN_DOC
 ! temporay arrays for dress_with_alpha_buffer. Avoids realocation.
END_DOC
END_PROVIDER

subroutine dress_with_alpha_buffer(delta_ij_loc, minilist, n_minilist, alpha, iproc)
 use bitmasks
 implicit none
  BEGIN_DOC
  !delta_ij_loc(:,:,1) : dressing column for H
  !delta_ij_loc(:,:,2) : dressing column for S2
  !minilist : indices of determinants connected to alpha ( in psi_det_sorted )
  !n_minilist : size of minilist
  !alpha : alpha determinant
  END_DOC
  integer(bit_kind), intent(in)   :: alpha(N_int,2)
  integer,intent(in)              :: minilist(n_minilist), n_minilist, iproc
  double precision, intent(inout) :: delta_ij_loc(N_states,N_det,2)


  integer                        :: i,j,k,l,m
  integer                        :: degree1, degree2, degree

  double precision               :: hIk, hla, hIl, sla, dIk(N_states), dka(N_states), dIa(N_states), hka
  double precision               :: phase, phase2
  integer                        :: exc(0:2,2,2)
  integer                        :: h1,h2,p1,p2,s1,s2
  integer(bit_kind)              :: tmp_det(N_int,2), ctrl
  integer                        :: i_state, k_sd, l_sd, m_sd, ll_sd, i_I
  double precision :: Delta_E_inv(N_states)
  double precision :: sdress, hdress
  logical :: ok, ok2
  integer :: old_ninc
  double precision :: shdress
  
  PROVIDE mo_class

  
  
  if(n_minilist == 1) return

  shdress = 0d0
  old_ninc = ninc

  if (perturbative_triples) then
    PROVIDE one_anhil fock_virt_total fock_core_inactive_total one_creat
  endif
  
  do i_I=1,N_det_ref
    call get_excitation_degree(alpha,psi_ref(1,1,i_I),degree1,N_int)
    if(degree1 <= 2) return
  end do

  

  ll_sd = 0
  do l_sd=1,n_minilist
    ok = .true.
    k_sd = minilist(l_sd)
    !if(idx_non_ref_rev(k_sd) == 0) cycle
    
    !do i_I=1,N_det_ref
    !  call get_excitation_degree(psi_det_sorted(1,1,k_sd),psi_ref(1,1,i_I),degree1,N_int)
    !  if(degree1 == 0) then
    !    ok = .false.
    !    exit
    !  end if
    !end do
    if(idx_non_ref_from_sorted(k_sd) == 0) ok = .false.

   !if(ok) then
   !  call get_excitation(psi_det_sorted(1,1,k_sd),alpha,exc,degree1,phase,N_int)
   !  if(degree1 == 0 .or. degree1 > 2) stop "minilist error"
   !  call decode_exc(exc,degree1,h1,p1,h2,p2,s1,s2)
   !  
   !  ok = (mo_class(h1)(1:1) == 'A' .or. mo_class(h1)(1:1) == 'I') .and. &
   !       (mo_class(p1)(1:1) == 'A' .or. mo_class(p1)(1:1) == 'V') 
   !  if(ok .and. degree1 == 2) then
   !          ok = (mo_class(h2)(1:1) == 'A' .or. mo_class(h2)(1:1) == 'I') .and. &
   !               (mo_class(p2)(1:1) == 'A' .or. mo_class(p2)(1:1) == 'V') 
   !  end if
   !end if
    
    if(ok) then
      ll_sd += 1
      idx_alpha_(ll_sd,iproc) = k_sd
      call i_h_j(alpha,psi_det_sorted(1,1,k_sd),N_int,hij_cache_(k_sd,iproc))
      call get_s2(alpha,psi_det_sorted(1,1,k_sd),N_int,sij_cache_(k_sd,iproc))
    end if
  enddo
  if(ll_sd <= 1) return 
  idx_alpha_(0,iproc) = ll_sd
  

  do i_I=1,N_det_ref
    call get_excitation_degree(alpha,psi_ref(1,1,i_I),degree1,N_int)
    if (degree1 > 4) then
      cycle
    endif
    
    do i_state=1,N_states
      dIa(i_state) = 0.d0
    enddo

    do k_sd=1,idx_alpha_(0,iproc)
      call get_excitation_degree(psi_ref(1,1,i_I),psi_det_sorted(1,1,idx_alpha_(k_sd,iproc)),degree,N_int)
!      print *, "diden"
      if (degree > 2) then
        cycle
      endif
      
      call get_excitation(psi_det_sorted(1,1,idx_alpha_(k_sd,iproc)),alpha,exc,degree2,phase,N_int)
      !print *, "DEG", degree2
      call decode_exc(exc,degree2,h1,p1,h2,p2,s1,s2)
      do k=1,N_int
        tmp_det(k,1) = psi_ref(k,1,i_I)
        tmp_det(k,2) = psi_ref(k,2,i_I)
      enddo
      call apply_excitation(psi_ref(1,1,i_I), exc, tmp_det, ok, N_int)
      
      ok2 = .false.
      do i_state=1,N_states
        dIK(i_state) = dij(i_I, idx_non_ref_from_sorted(idx_alpha_(k_sd,iproc)), i_state)
        if(dIK(i_state) /= 0d0) then
          ok2 = .true.
          exit
        endif
      enddo
      if(.not. ok2) cycle
      
      ! <I| \l/ |alpha>
      do i_state=1,N_states
        dka(i_state) = 0.d0
      enddo

      if (ok) then
        do l_sd=k_sd+1,idx_alpha_(0,iproc)
          call get_excitation_degree(tmp_det,psi_det_sorted(1,1,idx_alpha_(l_sd,iproc)),degree,N_int)
          if (degree == 0) then
            call get_excitation(psi_ref(1,1,i_I),psi_det_sorted(1,1,idx_alpha_(l_sd,iproc)),exc,degree,phase2,N_int)
            do i_state=1,N_states
              dka(i_state) = dij(i_I, idx_non_ref_from_sorted(idx_alpha_(l_sd,iproc)), i_state) * phase * phase2
            enddo
            exit
          endif
        enddo
      else if (perturbative_triples) then
          hka = hij_cache_(idx_alpha_(k_sd,iproc),iproc)
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
          hka = hij_cache_(idx_alpha_(k_sd,iproc),iproc) - hka
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
    
    ok2 = .false.
    do i_state=1,N_states
      if(dIa(i_state) /= 0d0) ok2 = .true.
    enddo
    if(.not. ok2) cycle

    do l_sd=1,idx_alpha_(0,iproc)
      k_sd = idx_alpha_(l_sd,iproc)
      hla = hij_cache_(k_sd,iproc)
      sla = sij_cache_(k_sd,iproc)
      do i_state=1,N_states
  !     dIa_hla_(i_state,k_sd,iproc) = dIa(i_state) * hla
  !     dIa_sla_(i_state,k_sd,iproc) = dIa(i_state) * sla
  !   enddo
  ! enddo
  ! do l_sd=1,idx_alpha_(0,iproc)
  !   do i_state=1,N_states
        k_sd = idx_alpha_(l_sd,iproc)
        m_sd = psi_from_sorted(k_sd)
        hdress =  dIa(i_state) * hla * psi_ref_coef(i_I,i_state)
        sdress =  dIa(i_state) * sla * psi_ref_coef(i_I,i_state)
      ! hdress = dIa_hla_(i_state,k_sd,iproc) * psi_ref_coef(i_I,i_state)
      ! sdress = dIa_sla_(i_state,k_sd,iproc) * psi_ref_coef(i_I,i_state)
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


