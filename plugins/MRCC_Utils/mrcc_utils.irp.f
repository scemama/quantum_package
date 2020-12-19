use bitmasks

 BEGIN_PROVIDER [ integer, mrmode ]
  mrmode = 0
END_PROVIDER
 
 

 BEGIN_PROVIDER [ double precision, lambda_mrcc, (N_states, N_det_non_ref) ]
&BEGIN_PROVIDER [ integer, lambda_mrcc_pt2, (0:psi_det_size) ]
&BEGIN_PROVIDER [ integer, lambda_mrcc_kept, (0:psi_det_size) ]
  implicit none
  BEGIN_DOC
  ! cm/<Psi_0|H|D_m> or perturbative 1/Delta_E(m)
  END_DOC
  integer :: i,k
  double precision               :: ihpsi_current(N_states)
  integer                        :: i_pert_count
  double precision               :: hii, lambda_pert
  integer                        :: N_lambda_mrcc_pt2, N_lambda_mrcc_pt3
  
  i_pert_count = 0
  lambda_mrcc = 0.d0
  N_lambda_mrcc_pt2 = 0
  N_lambda_mrcc_pt3 = 0
  lambda_mrcc_pt2(0) = 0
  lambda_mrcc_kept(0) = 0

  do i=1,N_det_non_ref
    call i_h_psi(psi_non_ref(1,1,i), psi_ref, psi_ref_coef, N_int, N_det_ref,&
        size(psi_ref_coef,1), N_states,ihpsi_current)
    call i_H_j(psi_non_ref(1,1,i),psi_non_ref(1,1,i),N_int,hii)
    do k=1,N_states
      if (ihpsi_current(k) == 0.d0) then
        ihpsi_current(k) = 1.d-32
      endif
!      lambda_mrcc(k,i) = psi_non_ref_coef(i,k)/ihpsi_current(k) 
      lambda_mrcc(k,i) = min(-1.d-32,psi_non_ref_coef(i,k)/ihpsi_current(k) )
      lambda_pert = 1.d0 / (psi_ref_energy_diagonalized(k)-hii)
      if (lambda_pert / lambda_mrcc(k,i)  < 0.5d0) then
        ! Ignore lamdba
        i_pert_count += 1
        lambda_mrcc(k,i) = 0.d0
        if (lambda_mrcc_pt2(N_lambda_mrcc_pt2) /= i) then
          N_lambda_mrcc_pt2 += 1
          lambda_mrcc_pt2(N_lambda_mrcc_pt2) = i
        endif
      else
        ! Keep lamdba
        if (lambda_mrcc_kept(N_lambda_mrcc_pt3) /= i) then
          N_lambda_mrcc_pt3 += 1
          lambda_mrcc_kept(N_lambda_mrcc_pt3) = i
        endif
      endif
    enddo
  enddo
  lambda_mrcc_pt2(0) = N_lambda_mrcc_pt2
  lambda_mrcc_kept(0) = N_lambda_mrcc_pt3
  print*,'N_det_non_ref = ',N_det_non_ref
  print*,'psi_coef_ref_ratio = ',psi_ref_coef(2,1)/psi_ref_coef(1,1)
  print*,'lambda max = ',maxval(dabs(lambda_mrcc))
  print*,'Number of ignored determinants = ',i_pert_count  

END_PROVIDER

! BEGIN_PROVIDER [ double precision, lambda_mrcc, (N_states, N_det_non_ref) ]
!&BEGIN_PROVIDER [ integer, lambda_mrcc_pt2, (0:psi_det_size) ]
!&BEGIN_PROVIDER [ integer, lambda_mrcc_kept, (0:psi_det_size) ]
!&BEGIN_PROVIDER [ double precision, lambda_pert, (N_states, N_det_non_ref) ] 
!  implicit none
!  BEGIN_DOC
!  ! cm/<Psi_0|H|D_m> or perturbative 1/Delta_E(m)
!  END_DOC
!  integer :: i,k
!  double precision               :: ihpsi_current(N_states)
!  integer                        :: i_pert_count
!  double precision               :: hii, E2(N_states), E2var(N_states)
!  integer                        :: N_lambda_mrcc_pt2, N_lambda_mrcc_pt3
!  
!  i_pert_count = 0
!  lambda_mrcc = 0.d0
!  N_lambda_mrcc_pt2 = 0
!  N_lambda_mrcc_pt3 = 0
!  lambda_mrcc_pt2(0) = 0
!  lambda_mrcc_kept(0) = 0
!
!  E2 = 0.d0
!  E2var = 0.d0
!  do i=1,N_det_non_ref
!    call i_h_psi(psi_non_ref(1,1,i), psi_ref, psi_ref_coef, N_int, N_det_ref,&
!        size(psi_ref_coef,1), N_states,ihpsi_current)
!    call i_H_j(psi_non_ref(1,1,i),psi_non_ref(1,1,i),N_int,hii)
!    do k=1,N_states
!      if (ihpsi_current(k) == 0.d0) then
!        ihpsi_current(k) = 1.d-32
!      endif
!      lambda_mrcc(k,i) = psi_non_ref_coef(i,k)/ihpsi_current(k) 
!      lambda_pert(k,i) = 1.d0 / (psi_ref_energy_diagonalized(k)-hii)
!      E2(k) += ihpsi_current(k)*ihpsi_current(k) / (psi_ref_energy_diagonalized(k)-hii)
!      E2var(k) += ihpsi_current(k) * psi_non_ref_coef(i,k)
!    enddo
!  enddo
!
!  do i=1,N_det_non_ref
!    call i_h_psi(psi_non_ref(1,1,i), psi_ref, psi_ref_coef, N_int, N_det_ref,&
!        size(psi_ref_coef,1), N_states,ihpsi_current)
!    call i_H_j(psi_non_ref(1,1,i),psi_non_ref(1,1,i),N_int,hii)
!    do k=1,N_states
!      if (ihpsi_current(k) == 0.d0) then
!        ihpsi_current(k) = 1.d-32
!      endif
!      lambda_mrcc(k,i) = psi_non_ref_coef(i,k)/ihpsi_current(k) 
!      lambda_pert(k,i) = 1.d0 / (psi_ref_energy_diagonalized(k)-hii) * E2var(k)/E2(k)
!    enddo
!  enddo
!  lambda_mrcc_pt2(0) = N_lambda_mrcc_pt2
!  lambda_mrcc_kept(0) = N_lambda_mrcc_pt3
!  print*,'N_det_non_ref = ',N_det_non_ref
!  print*,'psi_coef_ref_ratio = ',psi_ref_coef(2,1)/psi_ref_coef(1,1)
!  print*,'lambda max = ',maxval(dabs(lambda_mrcc))
!  print*,'Number of ignored determinants = ',i_pert_count  
!
!END_PROVIDER



BEGIN_PROVIDER [ double precision, hij_mrcc, (N_det_non_ref,N_det_ref) ]
 implicit none
 BEGIN_DOC
 ! < ref | H | Non-ref > matrix
 END_DOC
 integer :: i_I, k_sd
  do i_I=1,N_det_ref
    do k_sd=1,N_det_non_ref
      call i_h_j(psi_ref(1,1,i_I),psi_non_ref(1,1,k_sd),N_int,hij_mrcc(k_sd,i_I))
    enddo
  enddo

END_PROVIDER




logical function is_generable(det1, det2, Nint)
  use bitmasks
  implicit none
  integer, intent(in) :: Nint
  integer(bit_kind) :: det1(Nint, 2), det2(Nint, 2)
  integer :: degree, f, exc(0:2, 2, 2), t
  integer :: h1, h2, p1, p2, s1, s2
  integer, external :: searchExc
  logical, external :: excEq
  double precision :: phase
  integer :: tmp_array(4)
  
  is_generable = .false.
  call get_excitation(det1, det2, exc, degree, phase, Nint)
  if(degree == -1) return
  if(degree == 0) then
    is_generable = .true.
    return
  end if
  if(degree > 2) stop "?22??"
  
  call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
  
  if(degree == 1) then
    h2 = h1
    p2 = p1
    s2 = s1
    h1 = 0
    p1 = 0
    s1 = 0
  end if
  
  if(h1 + (s1-1)*mo_tot_num < h2 + (s2-1)*mo_tot_num) then
    tmp_array = (/s1, h1, s2, h2/)
  else
    tmp_array = (/s2, h2, s1, h1/)
  end if
  f = searchExc(hh_exists(1,1), tmp_array,  hh_shortcut(0))

  if(p1 + (s1-1)*mo_tot_num < p2 + (s2-1)*mo_tot_num) then
    tmp_array = (/s1, p1, s2, p2/)
  else
    tmp_array = (/s2, p2, s1, p1/)
  end if
  if (f /= -1) then
    f = searchExc(pp_exists(1,hh_shortcut(f)), tmp_array, hh_shortcut(f+1)-hh_shortcut(f))
  endif

  is_generable = (f /= -1) 
end function



integer function searchDet(dets, det, n, Nint)
  implicit none
  use bitmasks
  
  integer(bit_kind),intent(in) :: dets(Nint,2,n), det(Nint,2)
  integer, intent(in) :: nint, n
  integer :: l, h, c
  integer, external :: detCmp
  logical, external :: detEq

  l = 1
  h = n
  do while(.true.)
    searchDet = (l+h)/2
    c = detCmp(dets(1,1,searchDet), det(1,1), Nint)
    if(c == 0) then
      return
    else if(c == 1) then
      h = searchDet-1
    else
      l = searchDet+1
    end if
    if(l > h) then
      searchDet = -1
      return
    end if
    
  end do
end function


integer function unsortedSearchDet(dets, det, n, Nint)
  implicit none
  use bitmasks
  
  integer(bit_kind),intent(in) :: dets(Nint,2,n), det(Nint,2)
  integer, intent(in) :: nint, n
  integer :: l, h, c
  integer, external :: detCmp
  logical, external :: detEq

  do l=1, n
    if(detEq(det, dets(1,1,l), N_int)) then
      unsortedSearchDet = l
      return
    end if
  end do
  unsortedSearchDet = -1
end function


integer function searchExc(excs, exc, n)
  implicit none
  use bitmasks
  
  integer, intent(in) :: n
  integer,intent(in) :: excs(4,n), exc(4)
  integer :: l, h, c
  integer, external :: excCmp
  logical, external :: excEq

  l = 1
  h = n
  do
    searchExc = (l+h)/2
    c = excCmp(excs(1,searchExc), exc(1))
    if(c == 0) return
    if(c == 1) then
      h = searchExc-1
    else
      l = searchExc+1
    end if
    if(l > h) then
      searchExc = -1
      return
    end if
  end do
end function


subroutine sort_det(key, idx, N_key, Nint)
  implicit none
  

  integer, intent(in)                   :: Nint, N_key
  integer(8),intent(inout)       :: key(Nint,2,N_key)
  integer,intent(inout)                   :: idx(N_key)
  integer(8)                     :: tmp(Nint, 2)
  integer                               :: tmpidx,i,ni
  
  do i=1,N_key
    idx(i) = i
  end do
  
  do i=N_key/2,1,-1
    call tamiser(key, idx, i, N_key, Nint, N_key)
  end do
  
  do i=N_key,2,-1
    do ni=1,Nint
      tmp(ni,1) = key(ni,1,i)
      tmp(ni,2) = key(ni,2,i)
      key(ni,1,i) = key(ni,1,1)
      key(ni,2,i) = key(ni,2,1)
      key(ni,1,1) = tmp(ni,1)
      key(ni,2,1) = tmp(ni,2)
    enddo

    tmpidx = idx(i)
    idx(i) = idx(1)
    idx(1) = tmpidx
    call tamiser(key, idx, 1, i-1, Nint, N_key)
  end do
end subroutine 


subroutine sort_exc(key, N_key)
  implicit none
  

  integer, intent(in)                   :: N_key
  integer,intent(inout)       :: key(4,N_key)
  integer                     :: tmp(4)
  integer                               :: i,ni
  
  
  do i=N_key/2,1,-1
    call tamise_exc(key, i, N_key, N_key)
  end do
  
  do i=N_key,2,-1
    do ni=1,4
      tmp(ni) = key(ni,i)
      key(ni,i) = key(ni,1)
      key(ni,1) = tmp(ni)
    enddo

    call tamise_exc(key, 1, i-1, N_key)
  end do
end subroutine 


logical function exc_inf(exc1, exc2)
  implicit none
  integer,intent(in) :: exc1(4), exc2(4)
  integer :: i
  exc_inf = .false.
  do i=1,4
    if(exc1(i) < exc2(i)) then
      exc_inf = .true.
      return
    else if(exc1(i) > exc2(i)) then
      return
    end if
  end do
end function


subroutine tamise_exc(key, no, n, N_key)
  use bitmasks
  implicit none
  
  BEGIN_DOC
! Uncodumented : TODO
  END_DOC
  integer,intent(in)            :: no, n, N_key
  integer,intent(inout)       :: key(4, N_key)
  integer                       :: k,j
  integer                     :: tmp(4)
  logical                       :: exc_inf
  integer                       :: ni
  
  k = no
  j = 2*k
  do while(j <= n)
    if(j < n) then
      if (exc_inf(key(1,j), key(1,j+1))) then
        j = j+1
      endif
    endif
    if(exc_inf(key(1,k), key(1,j))) then
      do ni=1,4
        tmp(ni)   = key(ni,k)
        key(ni,k) = key(ni,j)
        key(ni,j) = tmp(ni)
      enddo
      k = j
      j = k+k
    else
      return
    endif
  enddo
end subroutine


subroutine dec_exc(exc, h1, h2, p1, p2)
  implicit none
  integer, intent(in)            :: exc(0:2,2,2)
  integer, intent(out)           :: h1, h2, p1, p2
  integer                        :: degree, s1, s2
  
  degree = exc(0,1,1) + exc(0,1,2)
  
  h1 = 0
  h2 = 0
  p1 = 0
  p2 = 0
    
  if(degree == 0) return
  
  call decode_exc(exc, degree, h1, p1, h2, p2, s1, s2)
  
  h1 += mo_tot_num * (s1-1)
  p1 += mo_tot_num * (s1-1)
  
  if(degree == 2) then
    h2 += mo_tot_num * (s2-1)
    p2 += mo_tot_num * (s2-1)
    if(h1 > h2) then
      s1 = h1
      h1 = h2
      h2 = s1
    end if
    if(p1 > p2) then
      s1 = p1
      p1 = p2
      p2 = s1
    end if
  else
    h2 = h1
    p2 = p1
    p1 = 0
    h1 = 0
  end if
end subroutine


 BEGIN_PROVIDER [ integer, N_hh_exists ]
&BEGIN_PROVIDER [ integer, N_pp_exists ]
&BEGIN_PROVIDER [ integer, N_ex_exists ]
  implicit none
  integer :: exc(0:2, 2, 2), degree, n, on, s, l, i
  integer :: h1, h2, p1, p2
  double precision :: phase
  logical,allocatable :: hh(:,:) , pp(:,:)
  
  allocate(hh(0:mo_tot_num*2, 0:mo_tot_num*2))
  allocate(pp(0:mo_tot_num*2, 0:mo_tot_num*2))
  hh = .false.
  pp = .false.
  N_hh_exists = 0
  N_pp_exists = 0
  N_ex_exists = 0

  n = 0
  !TODO Openmp
  do i=1, N_det_ref
    do l=1, N_det_non_ref
      call get_excitation(psi_ref(1,1,i), psi_non_ref(1,1,l), exc, degree, phase, N_int)
      if(degree == -1) cycle
      call dec_exc(exc, h1, h2, p1, p2)
      N_ex_exists += 1
      if(.not. hh(h1,h2)) N_hh_exists = N_hh_exists + 1
      if(.not. pp(p1,p2)) N_pp_exists = N_pp_exists + 1
      hh(h1,h2) = .true.
      pp(p1,p2) = .true.
    end do
  end do
  N_pp_exists = min(N_ex_exists, N_pp_exists * N_hh_exists)
END_PROVIDER



 BEGIN_PROVIDER [ integer(bit_kind), psi_non_ref_sorted, (N_int, 2, N_det_non_ref) ]
&BEGIN_PROVIDER [ integer, psi_non_ref_sorted_idx, (N_det_non_ref) ]
  implicit none
  psi_non_ref_sorted = psi_non_ref
  call sort_det(psi_non_ref_sorted, psi_non_ref_sorted_idx, N_det_non_ref, N_int)
END_PROVIDER


 BEGIN_PROVIDER [ double precision, dIj_unique, (hh_nex, N_states) ]
&BEGIN_PROVIDER [ double precision, rho_mrcc, (N_det_non_ref, N_states) ]
  implicit none
  logical                        :: ok
  integer                        :: i, j, k, s, II, pp, ppp, hh, ind, wk, a_col, at_row
  integer, external              :: searchDet, unsortedSearchDet
  integer(bit_kind)              :: myDet(N_int, 2), myMask(N_int, 2)
  integer                        :: N, INFO, r1, r2
  double precision , allocatable :: AtB(:), x(:), x_new(:), A_val_mwen(:,:), t(:)
  double precision               :: norm, cx, res
  integer, allocatable           :: lref(:), A_ind_mwen(:)
  double precision               :: phase
  
  
  double precision, allocatable :: rho_mrcc_inact(:)
  integer                        :: a_coll, at_roww
  
  print *, "TI", hh_nex, N_det_non_ref

  allocate(rho_mrcc_inact(N_det_non_ref))
  allocate(x_new(hh_nex))
  allocate(x(hh_nex), AtB(hh_nex))

  do s=1,N_states

    AtB(:) = 0.d0
    !$OMP PARALLEL default(none) shared(k, psi_non_ref_coef, active_excitation_to_determinants_idx,&
        !$OMP   active_excitation_to_determinants_val, N_det_ref, hh_nex, N_det_non_ref)          &
        !$OMP private(at_row, a_col, i, j, r1, r2, wk, A_ind_mwen, A_val_mwen, a_coll, at_roww)&
        !$OMP shared(N_states,mrcc_col_shortcut, mrcc_N_col, AtB, mrcc_AtA_val, mrcc_AtA_ind, s, n_exc_active, active_pp_idx)
    
    !$OMP DO schedule(static, 100)
    do at_roww = 1, n_exc_active ! hh_nex
      at_row = active_pp_idx(at_roww)
      do i=1,active_excitation_to_determinants_idx(0,at_roww)
          AtB(at_row) = AtB(at_row) + psi_non_ref_coef(active_excitation_to_determinants_idx(i, at_roww), s) * active_excitation_to_determinants_val(s,i, at_roww)
      end do
    end do
    !$OMP END DO
   
    !$OMP END PARALLEL

    X(:) = 0d0
    
    
    do a_coll = 1, n_exc_active
      a_col = active_pp_idx(a_coll)
      X(a_col) = AtB(a_col)
    end do
    
    rho_mrcc_inact(:) = 0d0
    
    allocate(lref(N_det_ref))
    do hh = 1, hh_shortcut(0)
      do pp = hh_shortcut(hh), hh_shortcut(hh+1)-1
        if(is_active_exc(pp)) cycle
        lref = 0
        AtB(pp) = 0.d0
        do II=1,N_det_ref
          call apply_hole_local(psi_ref(1,1,II), hh_exists(1, hh), myMask, ok, N_int)
          if(.not. ok) cycle
          call apply_particle_local(myMask, pp_exists(1, pp), myDet, ok, N_int)
          if(.not. ok) cycle
          ind = searchDet(psi_non_ref_sorted(1,1,1), myDet(1,1), N_det_non_ref, N_int)
          if(ind == -1) cycle
          ind = psi_non_ref_sorted_idx(ind)
          call get_phase(myDet(1,1), psi_ref(1,1,II), phase, N_int)
          AtB(pp) += psi_non_ref_coef(ind, s) * psi_ref_coef(II, s) * phase
          lref(II) = ind
          if(phase < 0.d0) lref(II) = -ind
        end do
        X(pp) =  AtB(pp) 
        do II=1,N_det_ref
          if(lref(II) > 0) then
            rho_mrcc_inact(lref(II)) = psi_ref_coef(II,s) * X(pp)
          else if(lref(II) < 0) then
            rho_mrcc_inact(-lref(II)) = -psi_ref_coef(II,s) * X(pp)
          end if
        end do
      end do
    end do
    deallocate(lref)

    x_new = x
    
    double precision               :: factor, resold
    factor = 1.d0
    resold = huge(1.d0)

    do k=0,hh_nex/4
      res = 0.d0
      do a_coll = 1, n_exc_active
        a_col = active_pp_idx(a_coll)
        cx = 0.d0
        do i=mrcc_col_shortcut(a_coll), mrcc_col_shortcut(a_coll) + mrcc_N_col(a_coll) - 1
          cx = cx + x(mrcc_AtA_ind(i)) * mrcc_AtA_val(s,i)
        end do
        x_new(a_col) = AtB(a_col) + cx * factor
        res = res + (X_new(a_col) - X(a_col))*(X_new(a_col) - X(a_col))
        X(a_col) = X_new(a_col)
      end do
      
      if (res > resold) then
        factor = factor * 0.5d0
      endif
      
      if(iand(k, 127) == 0) then
        print *, k, res, 1.d0 - res/resold
      endif
      
      if ( res < 1d-10 ) then
        exit
      endif
      if ( (res/resold > 0.99d0) ) then
        exit
      endif
      resold = res

    end do
    dIj_unique(1:size(X), s) = X(1:size(X))
    print *, k, res, 1.d0 - res/resold


    do i=1,N_det_non_ref
      rho_mrcc(i,s) = 0.d0
    enddo

    do a_coll=1,n_exc_active
      a_col = active_pp_idx(a_coll)
      do j=1,N_det_non_ref
        i = active_excitation_to_determinants_idx(j,a_coll)
        if (i==0) exit
        if (rho_mrcc_inact(i) /= 0.d0) then
          call debug_det(psi_non_ref(1,1,i),N_int)
          stop
        endif
        rho_mrcc(i,s) = rho_mrcc(i,s) + active_excitation_to_determinants_val(s,j,a_coll) * dIj_unique(a_col,s)
      enddo
    end do

    double precision :: norm2_ref, norm2_inact, a, b, c, Delta
    ! Psi = Psi_ref + Psi_inactive + f*Psi_active
    ! Find f to normalize Psi

    norm2_ref = 0.d0
    do i=1,N_det_ref
      norm2_ref = norm2_ref + psi_ref_coef(i,s)*psi_ref_coef(i,s)
    enddo

    a = 0.d0
    do i=1,N_det_non_ref
      a = a + rho_mrcc(i,s)*rho_mrcc(i,s)
    enddo

    norm = a + norm2_ref
    print *, "norm : ", sqrt(norm)

    norm = sqrt((1.d0-norm2_ref)/a)

    ! Renormalize Psi+A.X
    do i=1,N_det_non_ref
      rho_mrcc(i,s) = rho_mrcc(i,s) * norm 
    enddo

!norm = norm2_ref
!do i=1,N_det_non_ref
!  norm = norm + rho_mrcc(i,s)**2
!enddo
!print *,  'check', norm
!stop

     
        
     norm = 0.d0
     double precision               :: f, g, gmax
     gmax = maxval(dabs(psi_non_ref_coef(:,s)))
     do i=1,N_det_non_ref
       if (lambda_type == 2) then
         f = 1.d0
       else
        if (rho_mrcc(i,s) == 0.d0) then
          cycle
        endif
        ! f is such that f.\tilde{c_i} = c_i
        f = psi_non_ref_coef(i,s) / rho_mrcc(i,s)

        ! Avoid numerical instabilities
        g = 2.d0+100.d0*exp(-20.d0*dabs(psi_non_ref_coef(i,s)/gmax))
        f = min(f, g)
        f = max(f,-g)

      endif

       norm = norm + (rho_mrcc(i,s)*f)**2
       rho_mrcc(i,s) = f
     enddo
     ! rho_mrcc now contains the mu_i factors

     print *,  'norm of |T Psi_0> = ', dsqrt(norm)
     if (norm > 1.d0) then
       stop 'Error : Norm of the SD larger than the norm of the reference.'
     endif

  end do

END_PROVIDER




BEGIN_PROVIDER [ double precision, dij, (N_det_ref, N_det_non_ref, N_states) ]
  integer :: s,i,j
  double precision, external :: get_dij_index
  print *, "computing amplitudes..."
  do s=1, N_states
    do i=1, N_det_non_ref
      do j=1, N_det_ref
        !DIR$ FORCEINLINE
        dij(j, i, s) = get_dij_index(j, i, s, N_int)
      end do
    end do
  end do
  print *, "done computing amplitudes"
END_PROVIDER



!double precision function f_fit(x)
!  implicit none
!  double precision :: x
!  f_fit = 0.d0
!  return
!  if (x < 0.d0) then
!    f_fit = 0.d0
!  else if (x < 1.d0) then
!    f_fit = 1.d0/0.367879441171442 * ( x**2 * exp(-x**2))
!  else
!    f_fit = 1.d0
!  endif
!end
!
!double precision function get_dij_index(II, i, s, Nint)
!  integer, intent(in) :: II, i, s, Nint
!  double precision, external :: get_dij
!  double precision :: HIi, phase, c, a, b, d
!
!  call i_h_j(psi_ref(1,1,II), psi_non_ref(1,1,i), Nint, HIi)
!  call get_phase(psi_ref(1,1,II), psi_non_ref(1,1,i), phase, N_int)
!
!  a = lambda_pert(s,i)
!  b = lambda_mrcc(s,i)
!  c = f_fit(a/b)
!
!  d = get_dij(psi_ref(1,1,II), psi_non_ref(1,1,i), s, Nint) * phase* rho_mrcc(i,s)
!
!  c = f_fit(a*HIi/d)
!
!  get_dij_index = HIi * a * c + (1.d0 - c) * d
!  get_dij_index = d
!  return
!
!  if(lambda_type == 0) then
!    call get_phase(psi_ref(1,1,II), psi_non_ref(1,1,i), phase, N_int)
!    get_dij_index = get_dij(psi_ref(1,1,II), psi_non_ref(1,1,i), s, Nint) * phase
!    get_dij_index = get_dij_index * rho_mrcc(i,s) 
!  else if(lambda_type == 1) then
!    call i_h_j(psi_ref(1,1,II), psi_non_ref(1,1,i), Nint, HIi)
!    get_dij_index = HIi * lambda_mrcc(s, i)
!  else if(lambda_type == 2) then
!    call get_phase(psi_ref(1,1,II), psi_non_ref(1,1,i), phase, N_int)
!    get_dij_index = get_dij(psi_ref(1,1,II), psi_non_ref(1,1,i), s, Nint) * phase
!    get_dij_index = get_dij_index * rho_mrcc(i,s) 
!  end if
!end function

double precision function get_dij_index(II, i, s, Nint)
  integer, intent(in) :: II, i, s, Nint
  double precision, external :: get_dij
  double precision :: HIi, phase

  if(lambda_type == 0) then
    call get_phase(psi_ref(1,1,II), psi_non_ref(1,1,i), phase, N_int)
    get_dij_index = get_dij(psi_ref(1,1,II), psi_non_ref(1,1,i), s, Nint) * phase
    get_dij_index = get_dij_index * rho_mrcc(i,s) 
  else if(lambda_type == 1) then
    call i_h_j(psi_ref(1,1,II), psi_non_ref(1,1,i), Nint, HIi)
    get_dij_index = HIi * lambda_mrcc(s, i)
  else if(lambda_type == 2) then
    call get_phase(psi_ref(1,1,II), psi_non_ref(1,1,i), phase, N_int)
    get_dij_index = get_dij(psi_ref(1,1,II), psi_non_ref(1,1,i), s, Nint) * phase
    get_dij_index = get_dij_index * rho_mrcc(i,s) 
  end if
end function


double precision function get_dij(det1, det2, s, Nint)
  use bitmasks
  implicit none
  integer, intent(in) :: s, Nint
  integer(bit_kind) :: det1(Nint, 2), det2(Nint, 2)
  integer :: degree, f, exc(0:2, 2, 2), t
  integer :: h1, h2, p1, p2, s1, s2
  integer, external :: searchExc
  logical, external :: excEq
  double precision :: phase
  integer :: tmp_array(4)
  
  get_dij = 0d0
  call get_excitation(det1, det2, exc, degree, phase, Nint)
  if(degree == -1) return
  if(degree == 0) then
    stop "get_dij"
  end if
  
  call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
  
  if(degree == 1) then
    h2 = h1
    p2 = p1
    s2 = s1
    h1 = 0
    p1 = 0
    s1 = 0
  end if
  
  if(h1 + (s1-1)*mo_tot_num < h2 + (s2-1)*mo_tot_num) then
    tmp_array = (/s1, h1, s2, h2/)
  else
    tmp_array = (/s2, h2, s1, h1/)
  end if
  f = searchExc(hh_exists(1,1), tmp_array,  hh_shortcut(0))

  if(f == -1) return
  
  if(p1 + (s1-1)*mo_tot_num < p2 + (s2-1)*mo_tot_num) then
    tmp_array = (/s1, p1, s2, p2/)
  else
    tmp_array = (/s2, p2, s1, p1/)
  end if
  t = searchExc(pp_exists(1,hh_shortcut(f)), tmp_array, hh_shortcut(f+1)-hh_shortcut(f))

  if(t /= -1) then
    get_dij = dIj_unique(t - 1 + hh_shortcut(f), s)
  end if
end function


 BEGIN_PROVIDER [ integer, hh_exists, (4, N_hh_exists) ]
&BEGIN_PROVIDER [ integer, pp_exists, (4, N_pp_exists) ]
&BEGIN_PROVIDER [ integer, hh_shortcut, (0:N_hh_exists + 1) ]
&BEGIN_PROVIDER [ integer, hh_nex ]
  implicit none
  BEGIN_DOC
  !
  ! hh_exists : 
  !
  ! pp_exists : 
  !
  ! hh_shortcut : 
  !
  ! hh_nex : Total number of excitation operators
  !
  END_DOC
  integer,allocatable :: num(:,:)
  integer :: exc(0:2, 2, 2), degree, n, on, s, l, i
  integer :: h1, h2, p1, p2
  double precision :: phase
  logical, external :: excEq
  
  allocate(num(4, N_ex_exists+1))
  
  hh_shortcut = 0
  hh_exists = 0
  pp_exists = 0
  num = 0
  
  n = 0
  do i=1, N_det_ref
    do l=1, N_det_non_ref
      call get_excitation(psi_ref(1,1,i), psi_non_ref(1,1,l), exc, degree, phase, N_int)
      if(degree == -1) cycle
      call dec_exc(exc, h1, h2, p1, p2)
      n += 1
      num(:, n) = (/h1, h2, p1, p2/)
    end do
  end do
  
  call sort_exc(num, n)
  
  hh_shortcut(0) = 1
  hh_shortcut(1) = 1
  hh_exists(:,1) = (/1, num(1,1), 1, num(2,1)/)
  pp_exists(:,1) = (/1, num(3,1), 1, num(4,1)/)
  s = 1
  do i=2,n
    if(.not. excEq(num(1,i), num(1,s))) then
      s += 1
      num(:, s) = num(:, i)
      pp_exists(:,s) = (/1, num(3,s), 1, num(4,s)/)
      if(hh_exists(2, hh_shortcut(0)) /= num(1,s) .or. &
            hh_exists(4, hh_shortcut(0)) /= num(2,s)) then
        hh_shortcut(0) += 1
        hh_shortcut(hh_shortcut(0)) = s
        hh_exists(:,hh_shortcut(0)) = (/1, num(1,s), 1, num(2,s)/)
      end if
    end if
  end do
  hh_shortcut(hh_shortcut(0)+1) = s+1
  
  if (hh_shortcut(0) > N_hh_exists) then
    print *,  'Error in ', irp_here
    print *,  'hh_shortcut(0) :', hh_shortcut(0)
    print *,  'N_hh_exists : ', N_hh_exists
    print *,  'Is your active space defined?'
    stop
  endif

  if (hh_shortcut(hh_shortcut(0)+1)-1 > N_pp_exists) then
    print *,  'Error 1 in ', irp_here
    print *,  'hh_shortcut(hh_shortcut(0)+1)-1 :', hh_shortcut(hh_shortcut(0)+1)-1 
    print *,  'N_pp_exists : ', N_pp_exists
    print *,  'Is your active space defined?'
    stop
  endif

  do s=2,4,2
    do i=1,hh_shortcut(0)
      if(hh_exists(s, i) == 0) then
        hh_exists(s-1, i) = 0
      else if(hh_exists(s, i) > mo_tot_num) then
        hh_exists(s, i) -= mo_tot_num
        hh_exists(s-1, i) = 2
      end if
    end do
    

    do i=1,hh_shortcut(hh_shortcut(0)+1)-1
      if(pp_exists(s, i) == 0) then
        pp_exists(s-1, i) = 0
      else if(pp_exists(s, i) > mo_tot_num) then
        pp_exists(s, i) -= mo_tot_num
        pp_exists(s-1, i) = 2
      end if
    end do
  end do
  hh_nex = hh_shortcut(hh_shortcut(0)+1)-1
END_PROVIDER


logical function excEq(exc1, exc2)
  implicit none
  integer, intent(in) :: exc1(4), exc2(4)
  integer :: i
  excEq = .false.
  do i=1, 4
    if(exc1(i) /= exc2(i)) return
  end do
  excEq = .true.
end function


integer function excCmp(exc1, exc2)
  implicit none
  integer, intent(in) :: exc1(4), exc2(4)
  integer :: i
  excCmp = 0
  do i=1, 4
    if(exc1(i) > exc2(i)) then
      excCmp = 1
      return
    else if(exc1(i) < exc2(i)) then
      excCmp = -1
      return
    end if
  end do
end function


subroutine apply_hole_local(det, exc, res, ok, Nint)
  use bitmasks
  implicit none
  integer, intent(in) :: Nint
  integer, intent(in) :: exc(4)
  integer :: s1, s2, h1, h2
  integer(bit_kind),intent(in) :: det(Nint, 2)
  integer(bit_kind),intent(out) :: res(Nint, 2)
  logical, intent(out) :: ok
  integer :: ii, pos
  
  ok = .false.
  s1 = exc(1)
  h1 = exc(2)
  s2 = exc(3)
  h2 = exc(4)
  res = det
  
  if(h1 /= 0) then
    ii = (h1-1)/bit_kind_size + 1 
    pos = iand(h1-1,bit_kind_size-1) ! mod 64
    if(iand(det(ii, s1), ishft(1_bit_kind, pos)) == 0_8) then
      return
    endif
    res(ii, s1) = ibclr(res(ii, s1), pos)
  end if
  
  ii = (h2-1)/bit_kind_size + 1 
  pos = iand(h2-1,bit_kind_size-1) ! mod 64
  if(iand(det(ii, s2), ishft(1_bit_kind, pos)) == 0_8) then
    return
  endif
  res(ii, s2) = ibclr(res(ii, s2), pos)
  ok = .true.
end subroutine


subroutine apply_particle_local(det, exc, res, ok, Nint)
  use bitmasks
  implicit none
  integer, intent(in) :: Nint
  integer, intent(in) :: exc(4)
  integer :: s1, s2, p1, p2
  integer(bit_kind),intent(in) :: det(Nint, 2)
  integer(bit_kind),intent(out) :: res(Nint, 2)
  logical, intent(out) :: ok
  integer :: ii, pos 
  
  ok = .false.
  s1 = exc(1)
  p1 = exc(2)
  s2 = exc(3)
  p2 = exc(4)
  res = det 
  
  if(p1 /= 0) then
    ii = (p1-1)/bit_kind_size + 1 
    pos = iand(p1-1,bit_kind_size-1)
    if(iand(det(ii, s1), ishft(1_bit_kind, pos)) /= 0_8) then
      return
    endif
    res(ii, s1) = ibset(res(ii, s1), pos)
  end if

  ii = (p2-1)/bit_kind_size + 1 
  pos = iand(p2-1,bit_kind_size-1)
  if(iand(det(ii, s2), ishft(1_bit_kind, pos)) /= 0_8) then
    return
  endif
  res(ii, s2) = ibset(res(ii, s2), pos)

  
  ok = .true.
end subroutine




