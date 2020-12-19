use bitmasks


BEGIN_PROVIDER [ integer, N_mrcc_teeth ]
  N_mrcc_teeth = 10
END_PROVIDER

BEGIN_PROVIDER [ double precision, mrcc_norm_acc, (0:N_det_non_ref, N_states) ]
&BEGIN_PROVIDER [ double precision, mrcc_norm, (0:N_det_non_ref, N_states) ]
&BEGIN_PROVIDER [ double precision, mrcc_teeth_size, (0:N_det_non_ref, N_states) ]
&BEGIN_PROVIDER [ integer, mrcc_teeth, (0:N_mrcc_teeth+1, N_states) ]
  implicit none
  integer :: i, j, st, nt
  double precision :: norm_sto, jump, norm_mwen, norm_loc
  
  if(N_states /= 1) stop "mrcc_sto may not work with N_states /= 1"
  
  do st=1,N_states

    !norm_non_ref = 1d0
    !do i=1,N_det_ref
    !  norm_non_ref -= psi_ref_coef(i,st)**2
    !end do
    
    mrcc_teeth(0,st) = 1
    !norm_non_sto = norm_non_ref
    norm_sto = 1d0
    do i=1,N_det
      mrcc_teeth(1,st) = i
      jump = (1d0 / dfloat(N_mrcc_teeth)) * norm_sto
      if(psi_coef_generators(i,1)**2 < jump / 2d0) exit
      !if(i==80) exit
      norm_sto -= psi_coef_generators(i,1)**2
    end do
    
    print *, "FIRST", mrcc_teeth(1,st)

    norm_loc = 0d0
    mrcc_norm_acc(0,st) = 0d0
    nt = 1
    
    do i=1,mrcc_teeth(1,st)-1
      mrcc_norm_acc(i,st) = mrcc_norm_acc(i-1,st) + psi_coef_generators(i,st)**2
    end do

    do i=mrcc_teeth(1,st), N_det_generators!-mrcc_teeth(1,st)+1
      norm_mwen = psi_coef_generators(i,st)**2!-1+mrcc_teeth(1,st),st)**2
      mrcc_norm_acc(i,st) = mrcc_norm_acc(i-1,st) + norm_mwen 
      norm_loc += norm_mwen
      if(norm_loc > (jump*dfloat(nt))) then
        nt = nt + 1
        mrcc_teeth(nt,st) = i
      end if
    end do
    if(nt > N_mrcc_teeth+1) then
      print *, "foireouse mrcc_teeth", nt, mrcc_teeth(nt,st), N_det_non_ref
      stop
    end if

    mrcc_teeth(N_mrcc_teeth+1,st) = N_det_non_ref+1
    !mrcc_norm_acc(:,st) = mrcc_norm_acc(:,st) / norm_non_ref
    norm_loc = 0d0
    do i=N_mrcc_teeth, 0, -1
      mrcc_teeth_size(i,st) = mrcc_norm_acc(mrcc_teeth(i+1,st)-1,st) - mrcc_norm_acc(mrcc_teeth(i,st)-1, st)
      mrcc_norm_acc(mrcc_teeth(i,st):mrcc_teeth(i+1,st)-1,st) -= mrcc_norm_acc(mrcc_teeth(i,st)-1, st)
      mrcc_norm_acc(mrcc_teeth(i,st):mrcc_teeth(i+1,st)-1,st) = &
          mrcc_norm_acc(mrcc_teeth(i,st):mrcc_teeth(i+1,st)-1,st) / mrcc_teeth_size(i,st)
      mrcc_norm(mrcc_teeth(i,st), st) = mrcc_norm_acc(mrcc_teeth(i,st), st)
      do j=mrcc_teeth(i,st)+1, mrcc_teeth(i+1,1)-1
        mrcc_norm(j,1) = mrcc_norm_acc(j, st) - mrcc_norm_acc(j-1, st)
      end do
    end do
  end do
END_PROVIDER
 

! BEGIN_PROVIDER [ double precision, delta_ij_mrcc_sto,(N_states,N_det_non_ref) ]
!&BEGIN_PROVIDER [ double precision, delta_ij_s2_mrcc_sto, (N_states,N_det_non_ref) ]
!  use bitmasks
!  implicit none
!  integer :: gen, h, p, n, t, i, j, h1, h2, p1, p2, s1, s2, iproc
!  integer(bit_kind) :: mask(N_int, 2), omask(N_int, 2)
!  integer(bit_kind),allocatable :: buf(:,:,:)
!  logical :: ok
!  logical, external :: detEq
!  integer, external :: omp_get_thread_num
!  double precision :: coefs(N_det_non_ref), myCoef
!  integer :: n_in_teeth
!  double precision :: contrib(N_states), curn, in_teeth_step, curlim, curnorm
!  
!  contrib = 0d0
!  read(*,*) n_in_teeth
!  !n_in_teeth = 2
!  in_teeth_step = 1d0 / dfloat(n_in_teeth)
!  !double precision :: delta_ij_mrcc_tmp,(N_states,N_det_non_ref) 
!  !double precision :: delta_ij_s2_mrcc_tmp(N_states,N_det_non_ref)
!
!  coefs = 0d0
!  coefs(:mrcc_teeth(1,1)-1) = 1d0
!  
!  do i=1,N_mrcc_teeth
!    print *, "TEETH SIZE", i, mrcc_teeth(i+1,1)-mrcc_teeth(i,1)
!    if(mrcc_teeth(i+1,1) - mrcc_teeth(i,1) <= n_in_teeth) then
!      coefs(mrcc_teeth(i,1):mrcc_teeth(i+1,1)-1) = 1d0
!    else if(.false.) then
!      curnorm = 0d0
!      curn = 0.5d0
!      curlim = curn / dfloat(n_in_teeth)
!      do j=mrcc_teeth(i,1), mrcc_teeth(i+1,1)-1
!        if(mrcc_norm_acc(j,1) >= curlim) then
!          coefs(j) = 1d0
!          curnorm += mrcc_norm(j,1)
!          do while(mrcc_norm_acc(j,1) > curlim)
!            curn += 1d0
!            curlim = curn / dfloat(n_in_teeth)
!          end do
!        end if
!      end do
!      do j=mrcc_teeth(i,1), mrcc_teeth(i+1,1)-1
!        coefs(j) = coefs(j) / curnorm ! 1d0 / norm computed in teeth
!      end do
!    else if(.true.) then
!      coefs(mrcc_teeth(i,1):mrcc_teeth(i,1)+n_in_teeth-1) = 1d0 / mrcc_norm_acc(mrcc_teeth(i,1)+n_in_teeth-1, 1)
!    else
!      curnorm = 0d0
!      n = mrcc_teeth(i+1,1) - mrcc_teeth(i,1)
!      do j=1,n_in_teeth
!        t = int((dfloat(j)-0.5d0) * dfloat(n) / dfloat(n_in_teeth)) + 1 + mrcc_teeth(i,1) - 1
!        curnorm += mrcc_norm(t,1)
!        coefs(t) = 1d0
!      end do
!      do j=mrcc_teeth(i,1), mrcc_teeth(i+1,1)-1
!        coefs(j) = coefs(j) / curnorm ! 1d0 / norm computed in teeth
!      end do
!    end if
!    !coefs(mrcc_teeth(i,1)) = 
!  end do
!
!  !coefs = coefs * dfloat(N_det_generators)
!
!
!  delta_ij_mrcc_sto = 0d0
!  delta_ij_s2_mrcc_sto = 0d0
!  PROVIDE dij
!  provide hh_shortcut psi_det_size! lambda_mrcc
!  !$OMP PARALLEL DO default(none)  schedule(dynamic) &
!  !$OMP shared(psi_ref, psi_non_ref, hh_exists, pp_exists, N_int, hh_shortcut) &
!  !$OMP shared(N_det_generators, coefs,N_det_non_ref, delta_ij_mrcc_sto) &
!  !$OMP shared(contrib,psi_det_generators, delta_ij_s2_mrcc_sto) &
!  !$OMP private(i,j,curnorm,myCoef, h, n, mask, omask, buf, ok, iproc) 
!  do gen= 1,N_det_generators
!    if(coefs(gen) == 0d0) cycle
!    myCoef = coefs(gen)
!    allocate(buf(N_int, 2, N_det_non_ref))
!    iproc = omp_get_thread_num() + 1
!    if(mod(gen, 1000) == 0) print *, "mrcc_sto ", gen, "/", N_det_generators
!    
!    do h=1, hh_shortcut(0)
!      call apply_hole_local(psi_det_generators(1,1,gen), hh_exists(1, h), mask, ok, N_int)
!      if(.not. ok) cycle
!      omask = 0_bit_kind
!      if(hh_exists(1, h) /= 0) omask = mask
!      n = 1
!      do p=hh_shortcut(h), hh_shortcut(h+1)-1
!        call apply_particle_local(mask, pp_exists(1, p), buf(1,1,n), ok, N_int)
!        if(ok) n = n + 1
!        if(n > N_det_non_ref) stop "Buffer too small in MRCC..."
!      end do
!      n = n - 1
!      if(n /= 0) then
!        call mrcc_part_dress(delta_ij_mrcc_sto, delta_ij_s2_mrcc_sto, &
!            gen,n,buf,N_int,omask,myCoef,contrib)
!      endif
!    end do
!    deallocate(buf)
!    end do
!    !$OMP END PARALLEL DO
!
!
!
!    curnorm = 0d0
!    do j=1,N_det_non_ref
!      curnorm += delta_ij_mrcc_sto(1,j)*delta_ij_mrcc_sto(1,j)
!    end do
!    print *, "NORM DELTA ", dsqrt(curnorm)
!
!END_PROVIDER



 BEGIN_PROVIDER [ double precision, delta_ij_cancel, (N_states,N_det_non_ref) ]
&BEGIN_PROVIDER [ double precision, delta_ij_s2_cancel, (N_states,N_det_non_ref) ]
  use bitmasks
  implicit none


  integer                         :: i_state, i, i_I, J, k, k2, k1, kk, ll, m,l, deg, ni, m2
  integer                         :: n(2)
  integer                         :: p1,p2,h1,h2,s1,s2, blok, I_s, J_s, kn
  logical                         :: ok
  double precision                :: phase_ia, phase_Ik, phase_Jl, phase_Ji,phase_la, phase_ka, phase_tmp
  double precision                :: Hka, Hla, Ska, Sla, tmp
  double precision                :: diI, hIi, hJi, delta_JI, dkI, HkI,ci_inv(N_states), cj_inv(N_states)
  double precision                :: contrib, contrib_s2, wall, iwall
  integer, dimension(0:2,2,2)     :: exc_iI, exc_Ik, exc_IJ, exc
  integer(bit_kind)               :: det_tmp(N_int, 2), det_tmp2(N_int, 2),inac, virt
  integer, external               :: get_index_in_psi_det_sorted_bit, searchDet,detCmp
  logical, external               :: is_in_wavefunction
     
  provide dij
  
  delta_ij_cancel = 0d0


  do i=1,N_det_ref
    !$OMP PARALLEL DO default(shared) private(kk, k, blok, exc_Ik,det_tmp2,ok,deg,phase_Ik, l,ll) &
    !$OMP private(contrib, contrib_s2, i_state)
    do kk = 1, nlink(i)
      k = det_cepa0_idx(linked(kk, i))
      blok = blokMwen(kk, i)
      call get_excitation(psi_ref(1,1,i),psi_non_ref(1,1,k),exc_Ik,deg,phase_Ik,N_int)
      
      do j=1,N_det_ref
        if(j == i) cycle
        call apply_excitation(psi_ref(1,1,J),exc_Ik,det_tmp2,ok,N_int)
        if(.not. ok) cycle
          
        l = searchDet(det_cepa0(1,1,cepa0_shortcut(blok)), det_tmp2,cepa0_shortcut(blok+1)-cepa0_shortcut(blok), N_int)
        if(l == -1) cycle
        ll = cepa0_shortcut(blok)-1+l
        l = det_cepa0_idx(ll)
        ll = child_num(ll, J)
     
        do i_state = 1, N_states 
          contrib =  (dij(j, l, i_state) - dij(i, k, i_state)) * delta_cas(i,j,i_state)! * Hla *phase_ia * phase_ik
          contrib_s2 =  dij(j, l, i_state) - dij(i, k, i_state)! * Sla*phase_ia * phase_ik
          !$OMP ATOMIC
          delta_ij_cancel(i_state,l) += contrib * psi_ref_coef(i,i_state) 
          !$OMP ATOMIC
          delta_ij_s2_cancel(i_state,l) += contrib_s2* psi_ref_coef(i,i_state) 
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  enddo

END_PROVIDER


 BEGIN_PROVIDER [ double precision, delta_ij_mrcc, (N_states,N_det_non_ref) ]
&BEGIN_PROVIDER [ double precision, delta_ij_s2_mrcc, (N_states,N_det_non_ref) ]
  use bitmasks
  implicit none
  integer :: gen, h, p, n, t, i, h1, h2, p1, p2, s1, s2, iproc
  integer(bit_kind) :: mask(N_int, 2), omask(N_int, 2)
  integer(bit_kind),allocatable :: buf(:,:,:)
  logical :: ok
  logical, external :: detEq
  integer, external :: omp_get_thread_num
  double precision :: contrib(N_states)
  
  provide dij hh_shortcut psi_det_size  
  
  contrib = 0d0
  delta_ij_mrcc = 0d0
  delta_ij_s2_mrcc = 0d0


  !$OMP PARALLEL DO default(none)  schedule(dynamic) &
  !$OMP shared(contrib,psi_det_generators, N_det_generators, hh_exists, pp_exists, N_int, hh_shortcut) &
  !$OMP shared(N_det_non_ref, N_det_ref, delta_ij_mrcc, delta_ij_s2_mrcc) &
  !$OMP private(h, n, mask, omask, buf, ok, iproc) 
  do gen= 1, N_det_generators
    allocate(buf(N_int, 2, N_det_non_ref))
    iproc = omp_get_thread_num() + 1
    if(mod(gen, 1000) == 0) print *, "mrcc ", gen, "/", N_det_generators
    do h=1, hh_shortcut(0)
      call apply_hole_local(psi_det_generators(1,1,gen), hh_exists(1, h), mask, ok, N_int)
      if(.not. ok) cycle
      omask = 0_bit_kind
      if(hh_exists(1, h) /= 0) omask = mask
      n = 1
      do p=hh_shortcut(h), hh_shortcut(h+1)-1
        call apply_particle_local(mask, pp_exists(1, p), buf(1,1,n), ok, N_int)
        if(ok) n = n + 1
        if(n > N_det_non_ref) stop "Buffer too small in MRCC..."
      end do
      n = n - 1

      if(n /= 0) then
        call mrcc_part_dress(delta_ij_mrcc, delta_ij_s2_mrcc, gen,n,buf,N_int,omask,contrib)
      endif

    end do
    deallocate(buf)
  end do
  !$OMP END PARALLEL DO
END_PROVIDER


! subroutine blit(b1, b2)
!   double precision :: b1(N_states,N_det_non_ref), b2(N_states,N_det_non_ref)
!   b1 = b1 + b2
! end subroutine


subroutine mrcc_part_dress(delta_ij_, delta_ij_s2_, i_generator,n_selected,det_buffer,Nint,key_mask,contrib)
 use bitmasks
 implicit none

  integer, intent(in)            :: i_generator,n_selected, Nint
  double precision, intent(inout) :: delta_ij_(N_states,N_det_non_ref)
  double precision, intent(inout) :: delta_ij_s2_(N_states,N_det_non_ref)

  integer(bit_kind), intent(in)  :: det_buffer(Nint,2,n_selected)
  integer                        :: i,j,k,l,m
  integer,allocatable            :: idx_alpha(:), degree_alpha(:)
  logical                        :: good, fullMatch

  integer(bit_kind),allocatable  :: tq(:,:,:)
  integer                        :: N_tq, c_ref ,degree1, degree2, degree

  double precision               :: hIk, hla, hIl, sla, dIk(N_states), dka(N_states), dIa(N_states), hka
  double precision, allocatable  :: dIa_hla(:,:), dIa_sla(:,:)
  double precision               :: haj, phase, phase2
  double precision               :: f(N_states), ci_inv(N_states)
  integer                        :: exc(0:2,2,2)
  integer                        :: h1,h2,p1,p2,s1,s2
  integer(bit_kind)              :: tmp_det(Nint,2)
  integer                        :: iint, ipos
  integer                        :: i_state, k_sd, l_sd, i_I, i_alpha
  
  integer(bit_kind),allocatable  :: miniList(:,:,:)
  integer(bit_kind),intent(in)   :: key_mask(Nint, 2)
  integer,allocatable            :: idx_miniList(:)
  integer                        :: N_miniList, ni, leng
  double precision, allocatable  :: hij_cache(:), sij_cache(:)
  
  integer(bit_kind), allocatable :: microlist(:,:,:), microlist_zero(:,:,:)
  integer, allocatable           :: idx_microlist(:), N_microlist(:), ptr_microlist(:), idx_microlist_zero(:)
  integer :: mobiles(2), smallerlist
  logical, external :: detEq, is_generable
  !double precision, external :: get_dij, get_dij_index
  double precision :: Delta_E_inv(N_states)
  double precision, intent(inout) :: contrib(N_states)
  double precision :: sdress, hdress


  if (perturbative_triples) then
    PROVIDE one_anhil fock_virt_total fock_core_inactive_total one_creat
  endif


  leng = max(N_det_generators, N_det_non_ref)
  allocate(miniList(Nint, 2, leng), tq(Nint,2,n_selected), idx_minilist(leng), hij_cache(N_det_non_ref), sij_cache(N_det_non_ref))
  allocate(idx_alpha(0:psi_det_size), degree_alpha(psi_det_size))
  call create_minilist_find_previous(key_mask, psi_det_generators, miniList, i_generator-1, N_miniList, fullMatch, Nint)
  
  allocate(ptr_microlist(0:mo_tot_num*2+1),  &
      N_microlist(0:mo_tot_num*2) )
  allocate(   microlist(Nint,2,N_minilist*4),               &
    idx_microlist(N_minilist*4))
  
  if(key_mask(1,1) /= 0) then
    call create_microlist(miniList, N_minilist, key_mask, microlist, idx_microlist, N_microlist, ptr_microlist, Nint)
    call filter_tq_micro(i_generator,n_selected,det_buffer,Nint,tq,N_tq,microlist,ptr_microlist,N_microlist,key_mask)
  else
    call filter_tq(i_generator,n_selected,det_buffer,Nint,tq,N_tq,miniList,N_minilist)
  end if
  
  
  
  deallocate(microlist, idx_microlist)
  
  allocate (dIa_hla(N_states,N_det_non_ref), dIa_sla(N_states,N_det_non_ref))
  
  ! |I>
  
  ! |alpha>
  
  if(N_tq > 0) then
    call create_minilist(key_mask, psi_non_ref, miniList, idx_minilist, N_det_non_ref, N_minilist, Nint)
    if(N_minilist == 0) return
    
    
    if(sum(abs(key_mask(1:N_int,1))) /= 0) then
      allocate(microlist_zero(Nint,2,N_minilist), idx_microlist_zero(N_minilist))
      
      allocate(   microlist(Nint,2,N_minilist*4),               &
        idx_microlist(N_minilist*4))
      call create_microlist(miniList, N_minilist, key_mask, microlist, idx_microlist, N_microlist, ptr_microlist, Nint)
      
      
      do i=0,mo_tot_num*2
        do k=ptr_microlist(i),ptr_microlist(i+1)-1
          idx_microlist(k) = idx_minilist(idx_microlist(k))
        end do
      end do
      
      do l=1,N_microlist(0)
        do k=1,Nint
          microlist_zero(k,1,l) = microlist(k,1,l)
          microlist_zero(k,2,l) = microlist(k,2,l)
        enddo
        idx_microlist_zero(l) = idx_microlist(l)
      enddo
    end if
  end if
      
 
  do i_alpha=1,N_tq
    if(key_mask(1,1) /= 0) then
      call getMobiles(tq(1,1,i_alpha), key_mask, mobiles, Nint)
      
      if(N_microlist(mobiles(1)) < N_microlist(mobiles(2))) then
        smallerlist = mobiles(1)
      else
        smallerlist = mobiles(2)
      end if
      
      
      do l=0,N_microlist(smallerlist)-1
        microlist_zero(:,:,ptr_microlist(1) + l) = microlist(:,:,ptr_microlist(smallerlist) + l)
        idx_microlist_zero(ptr_microlist(1) + l) = idx_microlist(ptr_microlist(smallerlist) + l)
      end do
      
      call get_excitation_degree_vector(microlist_zero,tq(1,1,i_alpha),degree_alpha,Nint,N_microlist(smallerlist)+N_microlist(0),idx_alpha)
      do j=1,idx_alpha(0)
        idx_alpha(j) = idx_microlist_zero(idx_alpha(j))
      end do
      
    else
      call get_excitation_degree_vector(miniList,tq(1,1,i_alpha),degree_alpha,Nint,N_minilist,idx_alpha)
      do j=1,idx_alpha(0)
        idx_alpha(j) = idx_miniList(idx_alpha(j))
      end do
    end if
    
    do l_sd=1,idx_alpha(0)
      k_sd = idx_alpha(l_sd)
      call i_h_j(tq(1,1,i_alpha),psi_non_ref(1,1,idx_alpha(l_sd)),Nint,hij_cache(k_sd))
      call get_s2(tq(1,1,i_alpha),psi_non_ref(1,1,idx_alpha(l_sd)),Nint,sij_cache(k_sd))
      !if(sij_cache(k_sd) /= 0D0) PRINT *, "SIJ ", sij_cache(k_sd)
    enddo
    
    ! |I>
    do i_I=1,N_det_ref
      ! Find triples and quadruple grand parents
      call get_excitation_degree(tq(1,1,i_alpha),psi_ref(1,1,i_I),degree1,Nint)
      if (degree1 > 4) then
        cycle
      endif
      
      do i_state=1,N_states
        dIa(i_state) = 0.d0
      enddo
      
      ! <I|  <>  |alpha>
      do k_sd=1,idx_alpha(0)
        
        call get_excitation_degree(psi_ref(1,1,i_I),psi_non_ref(1,1,idx_alpha(k_sd)),degree,Nint)
        if (degree > 2) then
          cycle
        endif
        
        ! <I| /k\ |alpha>
        
        ! |l> = Exc(k -> alpha) |I>
        call get_excitation(psi_non_ref(1,1,idx_alpha(k_sd)),tq(1,1,i_alpha),exc,degree2,phase,Nint)
        call decode_exc(exc,degree2,h1,p1,h2,p2,s1,s2)
        do k=1,N_int
          tmp_det(k,1) = psi_ref(k,1,i_I)
          tmp_det(k,2) = psi_ref(k,2,i_I)
        enddo
        logical                        :: ok
        call apply_excitation(psi_ref(1,1,i_I), exc, tmp_det, ok, Nint)
        
        do i_state=1,N_states
          dIK(i_state) = dij(i_I, idx_alpha(k_sd), i_state)
        enddo
        
        ! <I| \l/ |alpha>
        do i_state=1,N_states
          dka(i_state) = 0.d0
        enddo
        
        if (ok) then
          do l_sd=k_sd+1,idx_alpha(0)
            call get_excitation_degree(tmp_det,psi_non_ref(1,1,idx_alpha(l_sd)),degree,Nint)
            if (degree == 0) then
              call get_excitation(psi_ref(1,1,i_I),psi_non_ref(1,1,idx_alpha(l_sd)),exc,degree,phase2,Nint)
              do i_state=1,N_states
                dka(i_state) = dij(i_I, idx_alpha(l_sd), i_state) * phase * phase2
              enddo
              exit
            endif
          enddo
          
        else if (perturbative_triples) then
          ! Linked
          
          hka = hij_cache(idx_alpha(k_sd))
          if (dabs(hka) > 1.d-12) then
            call get_delta_e_dyall_general_mp(psi_ref(1,1,i_I),tq(1,1,i_alpha),Delta_E_inv)
            
            do i_state=1,N_states
              ASSERT (Delta_E_inv(i_state) < 0.d0)
              dka(i_state) = hka / Delta_E_inv(i_state)
            enddo
          endif
          
        endif
        
        if (perturbative_triples.and. (degree2 == 1) ) then
          call i_h_j(psi_ref(1,1,i_I),tmp_det,Nint,hka)
          hka = hij_cache(idx_alpha(k_sd)) - hka
          if (dabs(hka) > 1.d-12) then
            call get_delta_e_dyall_general_mp(psi_ref(1,1,i_I),tq(1,1,i_alpha),Delta_E_inv)
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
          k_sd = idx_alpha(l_sd)
          hdress = dIa_hla(i_state,k_sd) * psi_ref_coef(i_I,i_state) 
          sdress = dIa_sla(i_state,k_sd) * psi_ref_coef(i_I,i_state) 
          !!$OMP ATOMIC
          !$OMP ATOMIC
          contrib(i_state) += hdress * psi_coef(dressed_column_idx(i_state), i_state) * psi_non_ref_coef(k_sd, i_state)
          !$OMP ATOMIC
          delta_ij_(i_state,k_sd) += hdress
          !$OMP ATOMIC
          delta_ij_s2_(i_state,k_sd) += sdress
        enddo
      enddo
    enddo
  enddo

  deallocate (dIa_hla,dIa_sla,hij_cache,sij_cache)
  deallocate(miniList, idx_miniList)
end


BEGIN_PROVIDER [ double precision, mrcc_previous_E, (N_states) ]
  implicit none
  BEGIN_DOC
    !energy difference between last two mrcc iterations
  END_DOC
  mrcc_previous_E(:) = ref_bitmask_energy
END_PROVIDER


 BEGIN_PROVIDER [ double precision, delta_ij_mrcc_zmq, (N_states,N_det_non_ref) ]
&BEGIN_PROVIDER [ double precision, delta_ij_s2_mrcc_zmq, (N_states,N_det_non_ref) ]
  use bitmasks
  implicit none

  integer                        :: i,j,k
  
  double precision, allocatable  :: mrcc(:)
  double precision               :: E_CI_before!, relative_error
  double precision, save :: target_error = 0d0

  allocate(mrcc(N_states))


  delta_ij_mrcc_zmq = 0d0
  delta_ij_s2_mrcc_zmq = 0d0

  !call random_seed()
  E_CI_before = mrcc_E0_denominator(1) + nuclear_repulsion
  threshold_selectors = 1.d0
  threshold_generators = 1d0 
  if(target_error /= 0d0) then
    target_error = target_error / 2d0 ! (-mrcc_E0_denominator(1) + mrcc_previous_E(1)) / 1d1
  else
    target_error = 1d-4
  end if
  target_error = 0d0
  call ZMQ_mrcc(E_CI_before, mrcc, delta_ij_mrcc_zmq, delta_ij_s2_mrcc_zmq, abs(target_error))

  mrcc_previous_E(:) = mrcc_E0_denominator(:)
END_PROVIDER


 BEGIN_PROVIDER [ double precision, delta_ij, (N_states,N_det_non_ref) ]
&BEGIN_PROVIDER [ double precision, delta_ij_s2, (N_states,N_det_non_ref) ]
  use bitmasks
  implicit none
  integer                        :: i, j, i_state
  !mrmode : 1=mrcepa0, 2=mrsc2 add, 3=mrcc, 4=stoch
!  if(mrmode == 4) then
!      do j = 1, N_det_non_ref
!        do i_state = 1, N_states
!          delta_ij(i_state,j) = delta_ij_mrcc_sto(i_state,j)
!          delta_ij_s2(i_state,j) = delta_ij_s2_mrcc_sto(i_state,j)
!        enddo
!      end do
!  else if(mrmode == 10) then
!      do j = 1, N_det_non_ref
!        do i_state = 1, N_states
!          delta_ij(i_state,j) = delta_ij_mrsc2(i_state,j)
!          delta_ij_s2(i_state,j) = delta_ij_s2_mrsc2(i_state,j)
!        enddo
!      end do
  if(mrmode == 5) then
      do j = 1, N_det_non_ref
        do i_state = 1, N_states
          delta_ij(i_state,j) = delta_ij_mrcc_zmq(i_state,j)
          delta_ij_s2(i_state,j) = delta_ij_s2_mrcc_zmq(i_state,j)
        enddo
      end do
  else  if(mrmode == 3) then
      do j = 1, N_det_non_ref
        do i_state = 1, N_states
          delta_ij(i_state,j) = delta_ij_mrcc(i_state,j)
          delta_ij_s2(i_state,j) = delta_ij_s2_mrcc(i_state,j)
        enddo
      end do
  else if(mrmode == 2) then
      do j = 1, N_det_non_ref
        do i_state = 1, N_states
          delta_ij(i_state,j) = delta_ij_old(i_state,j)
          delta_ij_s2(i_state,j) = delta_ij_s2_old(i_state,j)
        enddo
      end do
  else if(mrmode == 1) then
      do j = 1, N_det_non_ref
        do i_state = 1, N_states
          delta_ij(i_state,j) = delta_mrcepa0_ij(j,i_state)
          delta_ij_s2(i_state,j) = delta_mrcepa0_ij_s2(j,i_state)
        enddo
      end do
  else
    stop "invalid mrmode"
  end if

  !if(mrmode == 2 .or. mrmode == 3) then
  !    do j = 1, N_det_non_ref
  !      do i_state = 1, N_states
  !        delta_ij(i_state,j) += delta_ij_cancel(i_state,j)
  !      enddo
  !    end do
  !end if
END_PROVIDER


 BEGIN_PROVIDER [ integer, cepa0_shortcut, (0:N_det_non_ref+1) ]
&BEGIN_PROVIDER [ integer, det_cepa0_idx, (N_det_non_ref) ]
&BEGIN_PROVIDER [ integer(bit_kind), det_cepa0_active, (N_int,2,N_det_non_ref) ]
&BEGIN_PROVIDER [ integer(bit_kind), det_ref_active, (N_int,2,N_det_ref) ]
&BEGIN_PROVIDER [ integer(bit_kind), active_sorb, (N_int,2) ]
&BEGIN_PROVIDER [ integer(bit_kind), det_cepa0, (N_int,2,N_det_non_ref) ]
&BEGIN_PROVIDER [ integer, nlink, (N_det_ref) ]
&BEGIN_PROVIDER [ integer, linked, (N_det_non_ref,N_det_ref) ]
&BEGIN_PROVIDER [ integer, blokMwen, (N_det_non_ref,N_det_ref) ]
&BEGIN_PROVIDER [ double precision, searchance, (N_det_ref) ]
&BEGIN_PROVIDER [ integer, child_num, (N_det_non_ref,N_det_ref) ]

  use bitmasks
  implicit none
  
  integer(bit_kind),allocatable :: det_noactive(:,:,:)
  integer, allocatable :: shortcut(:), idx(:)
  integer(bit_kind) :: nonactive_sorb(N_int,2), det(N_int, 2)
  integer i, II, j, k, n, ni, blok, degree
  logical, external :: detEq
  
  allocate(det_noactive(N_int, 2, N_det_non_ref))
  allocate(idx(N_det_non_ref), shortcut(0:N_det_non_ref+1))
  print *, "pre start"
  active_sorb(:,:) = 0_8
  nonactive_sorb(:,:) = not(0_8)
  
  if(N_det_ref > 1) then
    do i=1, N_det_ref
    do k=1, N_int
      active_sorb(k,1) = ior(psi_ref(k,1,i), active_sorb(k,1))
      active_sorb(k,2) = ior(psi_ref(k,2,i), active_sorb(k,2))
      nonactive_sorb(k,1) = iand(psi_ref(k,1,i), nonactive_sorb(k,1))
      nonactive_sorb(k,2) = iand(psi_ref(k,2,i), nonactive_sorb(k,2))
    end do
    end do
    do k=1, N_int
      active_sorb(k,1) = iand(active_sorb(k,1), not(nonactive_sorb(k,1)))
      active_sorb(k,2) = iand(active_sorb(k,2), not(nonactive_sorb(k,2)))
    end do
  end if
  
  
  do i=1, N_det_non_ref
  do k=1, N_int
    det_noactive(k,1,i) = iand(psi_non_ref(k,1,i), not(active_sorb(k,1)))
    det_noactive(k,2,i) = iand(psi_non_ref(k,2,i), not(active_sorb(k,2))) 
  end do
  end do
  
  call sort_dets_ab(det_noactive, det_cepa0_idx, cepa0_shortcut, N_det_non_ref, N_int)
  
  do i=1,N_det_non_ref
    det_cepa0(:,:,i) = psi_non_ref(:,:,det_cepa0_idx(i))
  end do
  
  cepa0_shortcut(0) = 1
  cepa0_shortcut(1) = 1
  do i=2,N_det_non_ref
    if(.not. detEq(det_noactive(1,1,i), det_noactive(1,1,i-1), N_int)) then
      cepa0_shortcut(0) += 1
      cepa0_shortcut(cepa0_shortcut(0)) = i
    end if
  end do
  cepa0_shortcut(cepa0_shortcut(0)+1) = N_det_non_ref+1
  
  if(.true.) then 
  do i=1,cepa0_shortcut(0)
    n = cepa0_shortcut(i+1) - cepa0_shortcut(i)
    call sort_dets_ab(det_cepa0(1,1,cepa0_shortcut(i)), idx, shortcut, n, N_int)
    do k=1,n
      idx(k) = det_cepa0_idx(cepa0_shortcut(i)-1+idx(k))
    end do
    det_cepa0_idx(cepa0_shortcut(i):cepa0_shortcut(i)+n-1) = idx(:n)
  end do
  end if
  
  
  do i=1,N_det_ref
  do k=1, N_int
    det_ref_active(k,1,i) = iand(psi_ref(k,1,i), active_sorb(k,1))
    det_ref_active(k,2,i) = iand(psi_ref(k,2,i), active_sorb(k,2))
  end do
  end do
  
  do i=1,N_det_non_ref
  do k=1, N_int
    det_cepa0_active(k,1,i) = iand(psi_non_ref(k,1,det_cepa0_idx(i)), active_sorb(k,1))
    det_cepa0_active(k,2,i) = iand(psi_non_ref(k,2,det_cepa0_idx(i)), active_sorb(k,2))
  end do
  end do

  do i=1,N_det_non_ref
    if(.not. detEq(psi_non_ref(1,1,det_cepa0_idx(i)), det_cepa0(1,1,i),N_int)) stop "STOOOP"
  end do
  
  searchance = 0d0
  child_num = 0
  do J = 1, N_det_ref
    nlink(J) = 0
    do blok=1,cepa0_shortcut(0)
    do k=cepa0_shortcut(blok), cepa0_shortcut(blok+1)-1
      call get_excitation_degree(psi_ref(1,1,J),det_cepa0(1,1,k),degree,N_int)
      if(degree <= 2) then
        nlink(J) += 1
        linked(nlink(J),J) = k
        child_num(k, J) = nlink(J)
        blokMwen(nlink(J),J) = blok
        searchance(J) += 1d0 + log(dfloat(cepa0_shortcut(blok+1) - cepa0_shortcut(blok)))
      end if
    end do
    end do
  end do
  print *, "pre done"
END_PROVIDER

 
 
 BEGIN_PROVIDER [ double precision, delta_cas, (N_det_ref, N_det_ref, N_states) ]
&BEGIN_PROVIDER [ double precision, delta_cas_s2, (N_det_ref, N_det_ref, N_states) ]
  use bitmasks
  implicit none
  integer :: i,j,k
  double precision :: Sjk,Hjk, Hki, Hij
  !double precision, external :: get_dij
  integer i_state, degree
   
  provide lambda_mrcc dIj
  do i_state = 1, N_states
    !$OMP PARALLEL DO default(none) schedule(dynamic) private(j,k,Sjk,Hjk,Hki,degree) shared(lambda_mrcc,i_state, N_det_non_ref,psi_ref, psi_non_ref,N_int,delta_cas,delta_cas_s2,N_det_ref,dij)
    do i=1,N_det_ref
      do j=1,i
        call get_excitation_degree(psi_ref(1,1,i), psi_ref(1,1,j), degree, N_int)
        delta_cas(i,j,i_state) = 0d0
        delta_cas_s2(i,j,i_state) = 0d0
        do k=1,N_det_non_ref

          call i_h_j(psi_ref(1,1,j), psi_non_ref(1,1,k),N_int,Hjk)
          call get_s2(psi_ref(1,1,j), psi_non_ref(1,1,k),N_int,Sjk)
          
          delta_cas(i,j,i_state) += Hjk * dij(i, k, i_state) ! * Hki * lambda_mrcc(i_state, k)
          delta_cas_s2(i,j,i_state) += Sjk * dij(i, k, i_state) ! * Ski * lambda_mrcc(i_state, k)
        end do
        delta_cas(j,i,i_state) = delta_cas(i,j,i_state)
        delta_cas_s2(j,i,i_state) = delta_cas_s2(i,j,i_state)
      end do
    end do
    !$OMP END PARALLEL DO
  end do
 END_PROVIDER
 
 


logical function isInCassd(a,Nint)
  use bitmasks
  implicit none
  integer, intent(in) :: Nint
  integer(bit_kind), intent(in) :: a(Nint,2)
  integer(bit_kind) :: inac, virt
  integer :: ni, i, deg

  
  isInCassd = .false.
     
  deg = 0
   do i=1,2
    do ni=1,Nint
      virt = iand(not(HF_bitmask(ni,i)), not(active_sorb(ni,i)))
      deg += popcnt(iand(virt, a(ni,i)))
      if(deg > 2) return
    end do
  end do
  
  deg = 0
  do i=1,2
    do ni=1,Nint
      inac = iand(HF_bitmask(ni,i), not(active_sorb(ni,i)))
      deg += popcnt(xor(iand(inac,a(ni,i)), inac))
      if(deg > 2) return
    end do
  end do
  isInCassd = .true.
end function


subroutine getHP(a,h,p,Nint)
  use bitmasks
  implicit none
  integer, intent(in) :: Nint
  integer(bit_kind), intent(in) :: a(Nint,2)
  integer, intent(out) :: h, p
  integer(bit_kind) :: inac, virt
  integer :: ni, i, deg

  
  !isInCassd = .false.
  h = 0
  p = 0
   
  deg = 0
  lp : do i=1,2
    do ni=1,Nint
      virt = iand(not(HF_bitmask(ni,i)), not(active_sorb(ni,i)))
      deg += popcnt(iand(virt, a(ni,i)))
      if(deg > 2) exit lp
    end do
  end do lp
  p = deg

  deg = 0
  lh : do i=1,2
    do ni=1,Nint
      inac = iand(HF_bitmask(ni,i), not(active_sorb(ni,i)))
      deg += popcnt(xor(iand(inac,a(ni,i)), inac))
      if(deg > 2) exit lh
    end do
  end do lh
  h = deg
  !isInCassd = .true.
end subroutine


 BEGIN_PROVIDER [ double precision, delta_mrcepa0_ij, (N_det_non_ref,N_states) ]
&BEGIN_PROVIDER [ double precision, delta_mrcepa0_ij_s2, (N_det_non_ref,N_states) ]
  use bitmasks
  implicit none
  
  integer                         :: i_state, i, i_I, J, k, degree, degree2, m, l, deg, ni
  integer                         :: p1,p2,h1,h2,s1,s2, p1_,p2_,h1_,h2_,s1_,s2_, sortRefIdx(N_det_ref)
  logical                         :: ok
  double precision                :: phase_iI, phase_Ik, phase_Jl, phase_IJ, phase_al, diI, hIi, hJi, delta_JI, dkI(1), HkI, ci_inv(1), dia_hla(1)
  double precision                :: contrib, contrib_s2, HIIi, HJk, wall
  integer, dimension(0:2,2,2)     :: exc_iI, exc_Ik, exc_IJ
  integer(bit_kind)               :: det_tmp(N_int, 2), made_hole(N_int,2), made_particle(N_int,2), myActive(N_int,2)
  integer(bit_kind),allocatable   :: sortRef(:,:,:)
  integer, allocatable            :: idx_sorted_bit(:)
  integer, external               :: get_index_in_psi_det_sorted_bit, searchDet
  logical, external               :: is_in_wavefunction, detEq
  !double precision, external      :: get_dij
  integer :: II, blok
  integer*8, save :: notf = 0

  call wall_time(wall)
  allocate(idx_sorted_bit(N_det), sortRef(N_int,2,N_det_ref))
  
  sortRef(:,:,:) = det_ref_active(:,:,:)
  call sort_det(sortRef, sortRefIdx, N_det_ref, N_int)
  
  idx_sorted_bit(:) = -1
  do i=1,N_det_non_ref
    idx_sorted_bit(get_index_in_psi_det_sorted_bit(psi_non_ref(1,1,i), N_int)) = i
  enddo
    
  ! To provide everything
  contrib = dij(1, 1, 1)
  
  delta_mrcepa0_ij(:,:) = 0d0
  delta_mrcepa0_ij_s2(:,:) = 0d0

  do i_state = 1, N_states
    !$OMP PARALLEL DO default(none) schedule(dynamic) shared(delta_mrcepa0_ij, delta_mrcepa0_ij_s2)       &
    !$OMP private(m,i,II,J,k,degree,myActive,made_hole,made_particle,hjk,contrib,contrib_s2) &
    !$OMP shared(active_sorb, psi_non_ref, psi_non_ref_coef, psi_ref, psi_ref_coef, cepa0_shortcut, det_cepa0_active)     &
    !$OMP shared(N_det_ref, N_det_non_ref,N_int,det_cepa0_idx,lambda_mrcc,det_ref_active, delta_cas, delta_cas_s2) &
    !$OMP shared(notf,i_state, sortRef, sortRefIdx, dij)
    do blok=1,cepa0_shortcut(0)
    do i=cepa0_shortcut(blok), cepa0_shortcut(blok+1)-1
      do II=1,N_det_ref
        call get_excitation_degree(psi_ref(1,1,II),psi_non_ref(1,1,det_cepa0_idx(i)),degree,N_int)
        if (degree > 2 ) cycle
        
        do ni=1,N_int
          made_hole(ni,1) = iand(det_ref_active(ni,1,II), xor(det_cepa0_active(ni,1,i), det_ref_active(ni,1,II)))
          made_hole(ni,2) = iand(det_ref_active(ni,2,II), xor(det_cepa0_active(ni,2,i), det_ref_active(ni,2,II)))
          
          made_particle(ni,1) = iand(det_cepa0_active(ni,1,i), xor(det_cepa0_active(ni,1,i), det_ref_active(ni,1,II)))
          made_particle(ni,2) = iand(det_cepa0_active(ni,2,i), xor(det_cepa0_active(ni,2,i), det_ref_active(ni,2,II)))
        end do
      
        
  kloop: do k=cepa0_shortcut(blok), cepa0_shortcut(blok+1)-1 !i
          !if(lambda_mrcc(i_state, det_cepa0_idx(k)) == 0d0) cycle
          
          do ni=1,N_int
            if(iand(made_hole(ni,1), det_cepa0_active(ni,1,k)) /= 0) cycle kloop
            if(iand(made_particle(ni,1), det_cepa0_active(ni,1,k)) /= made_particle(ni,1)) cycle kloop
            if(iand(made_hole(ni,2), det_cepa0_active(ni,2,k)) /= 0) cycle kloop
            if(iand(made_particle(ni,2), det_cepa0_active(ni,2,k)) /= made_particle(ni,2)) cycle kloop
          end do
          do ni=1,N_int
            myActive(ni,1) = xor(det_cepa0_active(ni,1,k), made_hole(ni,1))
            myActive(ni,1) = xor(myActive(ni,1), made_particle(ni,1))
            myActive(ni,2) = xor(det_cepa0_active(ni,2,k), made_hole(ni,2))
            myActive(ni,2) = xor(myActive(ni,2), made_particle(ni,2))
          end do
          
          j = searchDet(sortRef, myActive, N_det_ref, N_int)
          if(j == -1) then
           cycle
          end if
          j = sortRefIdx(j)
          !$OMP ATOMIC
          notf = notf+1

          contrib = delta_cas(II, J, i_state)* dij(J, det_cepa0_idx(k), i_state)
          contrib_s2 = delta_cas_s2(II, J, i_state)  * dij(J, det_cepa0_idx(k), i_state)
          
          !$OMP ATOMIC
          delta_mrcepa0_ij(det_cepa0_idx(i), i_state) += contrib  * psi_ref_coef(J,i_state)
          delta_mrcepa0_ij_s2(det_cepa0_idx(i), i_state) += contrib_s2  * psi_ref_coef(J,i_state)

        end do kloop
      end do
    end do
    end do
    !$OMP END PARALLEL DO
  end do
  deallocate(idx_sorted_bit)
  call wall_time(wall)
  print *, "cepa0", wall, notf

END_PROVIDER


BEGIN_PROVIDER [ double precision, delta_sub_ij, (N_det_non_ref,N_states) ]
  use bitmasks
  implicit none
  
  integer                         :: i_state, i, i_I, J, k, degree, degree2, l, deg, ni
  integer                         :: p1,p2,h1,h2,s1,s2, p1_,p2_,h1_,h2_,s1_,s2_
  logical                         :: ok
  double precision                :: phase_Ji, phase_Ik, phase_Ii
  double precision                :: contrib, delta_IJk, HJk, HIk, HIl
  integer, dimension(0:2,2,2)     :: exc_Ik, exc_Ji, exc_Ii
  integer(bit_kind)               :: det_tmp(N_int, 2), det_tmp2(N_int, 2)
  integer, allocatable            :: idx_sorted_bit(:)
  integer, external               :: get_index_in_psi_det_sorted_bit
  
  integer :: II, blok
  
  provide delta_cas lambda_mrcc
  allocate(idx_sorted_bit(N_det))
  idx_sorted_bit(:) = -1
  do i=1,N_det_non_ref
    idx_sorted_bit(get_index_in_psi_det_sorted_bit(psi_non_ref(1,1,i), N_int)) = i
  enddo

    
  do i_state = 1, N_states
    delta_sub_ij(:,:) = 0d0
    
    provide mo_bielec_integrals_in_map
    
    
    !$OMP PARALLEL DO default(none) schedule(dynamic,10) shared(delta_sub_ij)       &
    !$OMP private(i, J, k, degree, degree2, l, deg, ni)       &
    !$OMP private(p1,p2,h1,h2,s1,s2, p1_,p2_,h1_,h2_,s1_,s2_)     &
    !$OMP private(ok, phase_Ji, phase_Ik, phase_Ii, contrib, delta_IJk, HJk, HIk, HIl, exc_Ik, exc_Ji, exc_Ii) &
    !$OMP private(det_tmp, det_tmp2, II, blok)    &
    !$OMP shared(idx_sorted_bit, N_det_non_ref, N_det_ref, N_int, psi_non_ref, psi_non_ref_coef, psi_ref, psi_ref_coef)   &
    !$OMP shared(i_state,lambda_mrcc, hf_bitmask, active_sorb)
    do i=1,N_det_non_ref
      if(mod(i,1000) == 0) print *, i, "/", N_det_non_ref
      do J=1,N_det_ref
        call get_excitation(psi_ref(1,1,J),psi_non_ref(1,1,i),exc_Ji,degree,phase_Ji,N_int)
        if(degree == -1) cycle
        
        
        do II=1,N_det_ref
          call apply_excitation(psi_ref(1,1,II),exc_Ji,det_tmp,ok,N_int)

          if(.not. ok) cycle
          l = get_index_in_psi_det_sorted_bit(det_tmp, N_int)
          if(l == 0) cycle
          l = idx_sorted_bit(l)
          
          call i_h_j(psi_ref(1,1,II), det_tmp, N_int, HIl)
          
          do k=1,N_det_non_ref
            if(lambda_mrcc(i_state, k) == 0d0) cycle
            call get_excitation(psi_ref(1,1,II),psi_non_ref(1,1,k),exc_Ik,degree2,phase_Ik,N_int)
            
            det_tmp(:,:) = 0_bit_kind
            det_tmp2(:,:) = 0_bit_kind
            
            ok = .true.
            do ni=1,N_int
              det_tmp(ni,1) = iand(xor(HF_bitmask(ni,1), psi_non_ref(ni,1,k)), not(active_sorb(ni,1)))
              det_tmp(ni,2) = iand(xor(HF_bitmask(ni,1), psi_non_ref(ni,1,i)), not(active_sorb(ni,1)))
              ok = ok .and. (popcnt(det_tmp(ni,1)) + popcnt(det_tmp(ni,2)) == popcnt(xor(det_tmp(ni,1), det_tmp(ni,2))))

              det_tmp(ni,1) = iand(xor(HF_bitmask(ni,2), psi_non_ref(ni,2,k)), not(active_sorb(ni,2)))
              det_tmp(ni,2) = iand(xor(HF_bitmask(ni,2), psi_non_ref(ni,2,i)), not(active_sorb(ni,2)))
              ok = ok .and. (popcnt(det_tmp(ni,1)) + popcnt(det_tmp(ni,2)) == popcnt(xor(det_tmp(ni,1), det_tmp(ni,2))))
            end do
            
            if(ok) cycle
            
            
            call i_h_j(psi_ref(1,1,J), psi_non_ref(1,1,k), N_int, HJk)
            call i_h_j(psi_ref(1,1,II), psi_non_ref(1,1,k), N_int, HIk)
            if(HJk == 0) cycle
            !assert HIk == 0
            delta_IJk = HJk * HIk * lambda_mrcc(i_state, k)
            call apply_excitation(psi_non_ref(1,1,i),exc_Ik,det_tmp,ok,N_int)
            if(ok) cycle
            contrib = delta_IJk * HIl * lambda_mrcc(i_state,l)   
            !$OMP ATOMIC
            delta_sub_ij(i, i_state) += contrib* psi_ref_coef(II,i_state)
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end do
  deallocate(idx_sorted_bit)
END_PROVIDER


subroutine set_det_bit(det, p, s)
   implicit none
   integer(bit_kind),intent(inout) :: det(N_int, 2)
   integer, intent(in) :: p, s
   integer :: ni, pos
 
   ni = (p-1)/bit_kind_size + 1
   pos = mod(p-1, bit_kind_size)
   det(ni,s) = ibset(det(ni,s), pos)
end subroutine


 BEGIN_PROVIDER [ double precision, h_cache, (N_det_ref,N_det_non_ref) ]
&BEGIN_PROVIDER [ double precision, s2_cache, (N_det_ref,N_det_non_ref) ]
  implicit none
  integer :: i,j
  do i=1,N_det_ref
  do j=1,N_det_non_ref
    call i_h_j(psi_ref(1,1,i), psi_non_ref(1,1,j), N_int, h_cache(i,j))
    call get_s2(psi_ref(1,1,i), psi_non_ref(1,1,j), N_int, s2_cache(i,j))
  end do
  end do
END_PROVIDER



subroutine filter_tq(i_generator,n_selected,det_buffer,Nint,tq,N_tq,miniList,N_miniList)

 use bitmasks
 implicit none

  integer, intent(in)            :: i_generator,n_selected, Nint

  integer(bit_kind), intent(in)  :: det_buffer(Nint,2,n_selected)
  integer                        :: i,j,k,m
  logical                        :: is_in_wavefunction
  integer,allocatable            :: degree(:)
  integer,allocatable            :: idx(:)
  logical                        :: good

  integer(bit_kind), intent(inout) :: tq(Nint,2,n_selected) !! intent(out)
  integer, intent(out)           :: N_tq
  
  integer                        :: nt,ni
  logical, external              :: is_connected_to, is_generable
  
  integer(bit_kind),intent(in)  :: miniList(Nint,2,N_det_generators)
  integer,intent(in)            :: N_miniList
  
  allocate(degree(psi_det_size)) 
  allocate(idx(0:psi_det_size))
  N_tq = 0

  i_loop : do i=1,N_selected
    do k=1, N_minilist
      if(is_generable(miniList(1,1,k), det_buffer(1,1,i), Nint)) cycle i_loop
    end do

    ! Select determinants that are triple or quadruple excitations
    ! from the ref
    good = .True.
    call get_excitation_degree_vector(psi_ref,det_buffer(1,1,i),degree,Nint,N_det_ref,idx) 
    !good=(idx(0) == 0) tant que degree > 2 pas retourné par get_excitation_degree_vector
    do k=1,idx(0)
      if (degree(k) < 3) then
        good = .False.
        exit
      endif
    enddo
    if (good) then
      if (.not. is_in_wavefunction(det_buffer(1,1,i),Nint)) then
        N_tq += 1
        do k=1,N_int
          tq(k,1,N_tq) = det_buffer(k,1,i)
          tq(k,2,N_tq) = det_buffer(k,2,i)
        enddo
      endif
    endif
  enddo i_loop
end


subroutine filter_tq_micro(i_generator,n_selected,det_buffer,Nint,tq,N_tq,microlist,ptr_microlist,N_microlist,key_mask)

 use bitmasks
 implicit none

  integer, intent(in)            :: i_generator,n_selected, Nint

  integer(bit_kind), intent(in)  :: det_buffer(Nint,2,n_selected)
  integer                        :: i,j,k,m
  logical                        :: is_in_wavefunction
  integer,allocatable            :: degree(:)
  integer,allocatable            :: idx(:)
  logical                        :: good

  integer(bit_kind), intent(inout) :: tq(Nint,2,n_selected) !! intent(out)
  integer, intent(out)           :: N_tq
  
  integer                        :: nt,ni
  logical, external              :: is_connected_to, is_generable

  integer(bit_kind),intent(in)  :: microlist(Nint,2,*)
  integer,intent(in)  :: ptr_microlist(0:*)
  integer,intent(in)            :: N_microlist(0:*)
  integer(bit_kind),intent(in)   :: key_mask(Nint, 2)
  
  integer :: mobiles(2), smallerlist
  
  
  allocate(degree(psi_det_size)) 
  allocate(idx(0:psi_det_size))
  N_tq = 0
  
  i_loop : do i=1,N_selected
    call getMobiles(det_buffer(1,1,i), key_mask, mobiles, Nint) 
    if(N_microlist(mobiles(1)) < N_microlist(mobiles(2))) then
      smallerlist = mobiles(1)
    else
      smallerlist = mobiles(2)
    end if
    
    if(N_microlist(smallerlist) > 0) then
      do k=ptr_microlist(smallerlist), ptr_microlist(smallerlist)+N_microlist(smallerlist)-1
        if(is_generable(microlist(1,1,k), det_buffer(1,1,i), Nint)) cycle i_loop
      end do
    end if
    
    if(N_microlist(0) > 0) then
      do k=1, N_microlist(0)
        if(is_generable(microlist(1,1,k), det_buffer(1,1,i), Nint)) cycle i_loop
      end do
    end if

    ! Select determinants that are triple or quadruple excitations
    ! from the ref
    good = .True.
    call get_excitation_degree_vector(psi_ref,det_buffer(1,1,i),degree,Nint,N_det_ref,idx) 
    !good=(idx(0) == 0) tant que degree > 2 pas retourné par get_excitation_degree_vector
    do k=1,idx(0)
      if (degree(k) < 3) then
        good = .False.
        exit
      endif
    enddo
    if (good) then
      if (.not. is_in_wavefunction(det_buffer(1,1,i),Nint)) then
        N_tq += 1
        do k=1,N_int
          tq(k,1,N_tq) = det_buffer(k,1,i)
          tq(k,2,N_tq) = det_buffer(k,2,i)
        enddo
      endif
    endif
  enddo i_loop
end




