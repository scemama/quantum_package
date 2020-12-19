BEGIN_PROVIDER [ integer, fragment_first ]
  implicit none
  fragment_first = first_det_of_teeth(1)
END_PROVIDER

subroutine ZMQ_pt2(E, pt2,relative_error, absolute_error, error)
  use f77_zmq
  use selection_types
  
  implicit none
  
  character(len=64000)           :: task
  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket, zmq_socket_pull
  type(selection_buffer)         :: b
  integer, external              :: omp_get_thread_num
  double precision, intent(in)   :: relative_error, absolute_error, E(N_states)
  double precision, intent(out)  :: pt2(N_states),error(N_states)
  
  
  double precision, allocatable  :: pt2_detail(:,:), comb(:)
  logical, allocatable           :: computed(:)
  integer, allocatable           :: tbc(:)
  integer                        :: i, j, k, Ncomb, i_generator_end
  integer, external              :: pt2_find
  
  double precision               :: sumabove(comb_teeth), sum2above(comb_teeth), Nabove(comb_teeth)
  double precision, external     :: omp_get_wtime
  double precision               :: state_average_weight_save(N_states), w(N_states)
  double precision               :: time
  integer(ZMQ_PTR), external     :: new_zmq_to_qp_run_socket
  
  if (N_det < max(10,N_states)) then
    pt2=0.d0
    call ZMQ_selection(0, pt2)
    error(:) = 0.d0
  else
    
    state_average_weight_save(:) = state_average_weight(:)
    do pt2_stoch_istate=1,N_states
      SOFT_TOUCH pt2_stoch_istate
      state_average_weight(:) = 0.d0
      state_average_weight(pt2_stoch_istate) = 1.d0
      TOUCH state_average_weight 
      
      allocate(pt2_detail(N_states,N_det_generators+1), comb(N_det_generators), computed(N_det_generators), tbc(0:size_tbc))
      sumabove = 0d0
      sum2above = 0d0
      Nabove = 0d0
      
      provide nproc fragment_first fragment_count mo_bielec_integrals_in_map mo_mono_elec_integral pt2_weight psi_selectors 
      
      computed = .false.
      
      tbc(0) = first_det_of_comb - 1
      do i=1, tbc(0)
        tbc(i) = i
        computed(i) = .true.
      end do
      
      Ncomb=size(comb)
      call get_carlo_workbatch(computed, comb, Ncomb, tbc)

      pt2_detail = 0d0
      print *, '========== ================= ================= ================='
      print *, ' Samples        Energy         Stat. Error         Seconds      '
      print *, '========== ================= ================= ================='
      
      call new_parallel_job(zmq_to_qp_run_socket, zmq_socket_pull, 'pt2')

      integer, external              :: zmq_put_psi
      integer, external              :: zmq_put_N_det_generators
      integer, external              :: zmq_put_N_det_selectors
      integer, external              :: zmq_put_dvector
      if (zmq_put_psi(zmq_to_qp_run_socket,1) == -1) then
        stop 'Unable to put psi on ZMQ server'
      endif
      if (zmq_put_N_det_generators(zmq_to_qp_run_socket, 1) == -1) then
        stop 'Unable to put N_det_generators on ZMQ server'
      endif
      if (zmq_put_N_det_selectors(zmq_to_qp_run_socket, 1) == -1) then
        stop 'Unable to put N_det_selectors on ZMQ server'
      endif
      if (zmq_put_dvector(zmq_to_qp_run_socket,1,'energy',pt2_e0_denominator,size(pt2_e0_denominator)) == -1) then
        stop 'Unable to put energy on ZMQ server'
      endif
      call create_selection_buffer(1, 1*2, b)
      
      integer                        :: ipos
      ipos=1

      integer, external :: add_task_to_taskserver
      
      do i=1,tbc(0)
        if(tbc(i) > fragment_first) then
          write(task(ipos:ipos+20),'(I9,1X,I9,''|'')') 0, tbc(i)
          ipos += 20
          if (ipos > 63980) then
            if (add_task_to_taskserver(zmq_to_qp_run_socket,trim(task(1:ipos))) == -1) then
              stop 'Unable to add task to task server'
            endif
            ipos=1
          endif
        else
          do j=1,fragment_count
            write(task(ipos:ipos+20),'(I9,1X,I9,''|'')') j, tbc(i)
            ipos += 20
            if (ipos > 63980) then
              if (add_task_to_taskserver(zmq_to_qp_run_socket,trim(task(1:ipos))) == -1) then
                stop 'Unable to add task to task server'
              endif
              ipos=1
            endif
          end do
        end if
      end do
      if (ipos > 1) then
        if (add_task_to_taskserver(zmq_to_qp_run_socket,trim(task(1:ipos))) == -1) then
          stop 'Unable to add task to task server'
        endif
      endif
      
      integer, external :: zmq_set_running
      if (zmq_set_running(zmq_to_qp_run_socket) == -1) then
        print *,  irp_here, ': Failed in zmq_set_running'
      endif

      
      !$OMP PARALLEL DEFAULT(shared) NUM_THREADS(nproc+1)            &
          !$OMP  PRIVATE(i)
      i = omp_get_thread_num()
      if (i==0) then
        call pt2_collector(zmq_socket_pull,E(pt2_stoch_istate), b, tbc, comb, Ncomb, computed, pt2_detail, sumabove, sum2above, Nabove, relative_error, absolute_error, w, error)
        pt2(pt2_stoch_istate) = w(pt2_stoch_istate)
      else
        call pt2_slave_inproc(i)
      endif
      !$OMP END PARALLEL
      call end_parallel_job(zmq_to_qp_run_socket, zmq_socket_pull, 'pt2')
      call delete_selection_buffer(b)
      
      print *, '========== ================= ================= ================='
      
      deallocate(pt2_detail, comb, computed, tbc)
    enddo
    FREE pt2_stoch_istate 
    state_average_weight(:) = state_average_weight_save(:)
    TOUCH state_average_weight
  endif
  do k=N_det+1,N_states
    pt2(k) = 0.d0
  enddo

end subroutine


subroutine do_carlo(tbc, Ncomb, comb, pt2_detail, computed, sumabove, sum2above, Nabove)
  integer, intent(in) :: tbc(0:size_tbc), Ncomb
  logical, intent(in) :: computed(N_det_generators)
  double precision, intent(in) :: comb(Ncomb), pt2_detail(N_states,N_det_generators)
  double precision, intent(inout) :: sumabove(comb_teeth), sum2above(comb_teeth), Nabove(comb_teeth)
  integer :: i, dets(comb_teeth)
  double precision :: myVal, myVal2

  mainLoop : do i=1,Ncomb
    call get_comb(comb(i), dets, comb_teeth)
    do j=1,comb_teeth
      if(.not.(computed(dets(j)))) then
        exit mainLoop
      end if
    end do
    
    myVal = 0d0
    myVal2 = 0d0
    do j=comb_teeth,1,-1
      myVal += pt2_detail(pt2_stoch_istate,dets(j)) * pt2_weight_inv(dets(j)) * comb_step
      sumabove(j) += myVal
      sum2above(j) += myVal*myVal
      Nabove(j) += 1
    end do
  end do mainLoop
end subroutine


subroutine pt2_slave_inproc(i)
  implicit none
  integer, intent(in)            :: i

  call run_pt2_slave(1,i,pt2_e0_denominator)
end

subroutine pt2_collector(zmq_socket_pull, E, b, tbc, comb, Ncomb, computed, pt2_detail, sumabove, sum2above, Nabove, relative_error, absolute_error, pt2,error)
  use f77_zmq
  use selection_types
  use bitmasks
  implicit none

  
  integer, intent(in) :: Ncomb
  integer(ZMQ_PTR), intent(in)   :: zmq_socket_pull
  double precision, intent(inout) :: pt2_detail(N_states, N_det_generators)
  double precision, intent(in) :: comb(Ncomb), relative_error, absolute_error, E
  logical, intent(inout) :: computed(N_det_generators)
  integer, intent(in) :: tbc(0:size_tbc)
  double precision, intent(inout) :: sumabove(comb_teeth), sum2above(comb_teeth), Nabove(comb_teeth)
  double precision, intent(out)  :: pt2(N_states),error(N_states)


  type(selection_buffer), intent(inout) :: b
  double precision, allocatable      :: pt2_mwen(:,:)
  integer(ZMQ_PTR),external      :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket


  integer :: msg_size, rc, more
  integer :: acc, i, j, robin, N, n_tasks
  double precision, allocatable :: val(:)
  integer(bit_kind), allocatable :: det(:,:,:)
  integer, allocatable :: task_id(:)
  integer, allocatable :: index(:)
  double precision :: time0
  double precision :: time, timeLast, Nabove_old
  double precision, external :: omp_get_wtime
  integer :: tooth, firstTBDcomb, orgTBDcomb, n_tasks_max
  integer, allocatable :: parts_to_get(:)
  logical, allocatable :: actually_computed(:)
  double precision :: eqt
  character*(512) :: task
  Nabove_old = -1.d0
  n_tasks_max = N_det_generators/100+1
  
  allocate(actually_computed(N_det_generators), parts_to_get(N_det_generators), &
    pt2_mwen(N_states, n_tasks_max) )

  pt2_mwen(1:N_states, 1:n_tasks_max) = 0.d0
  do i=1,N_det_generators
    actually_computed(i) = computed(i)
  enddo
  
  parts_to_get(:) = 1
  if(fragment_first > 0) then
    do i=1,fragment_first
      parts_to_get(i) = fragment_count
    enddo
  endif

  do i=1,tbc(0)
    actually_computed(tbc(i)) = .false.
  end do
  
  orgTBDcomb = int(Nabove(1))
  firstTBDcomb = 1

  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()
  allocate(val(b%N), det(N_int, 2, b%N), task_id(n_tasks_max), index(n_tasks_max))
  more = 1
  call wall_time(time0)
  timeLast = time0

  call get_first_tooth(actually_computed, tooth)
  Nabove_old = Nabove(tooth)
  
  logical :: loop
  loop = .True.
  pullLoop : do while (loop)

    call pull_pt2_results(zmq_socket_pull, index, pt2_mwen, task_id, n_tasks)
    do i=1,n_tasks
      pt2_detail(1:N_states, index(i)) += pt2_mwen(1:N_states,i)
      parts_to_get(index(i)) -= 1
      if(parts_to_get(index(i)) < 0) then 
        print *, i, index(i), parts_to_get(index(i))
        print *, "PARTS ??"
        print *, parts_to_get
        stop "PARTS ??"
      end if
      if(parts_to_get(index(i)) == 0) actually_computed(index(i)) = .true.
    enddo

    integer, external :: zmq_delete_tasks
    if (zmq_delete_tasks(zmq_to_qp_run_socket,zmq_socket_pull,task_id,n_tasks,more) == -1) then
        stop 'Unable to delete tasks'
    endif
    if (more == 0) then
      loop = .False.
    endif

    time = omp_get_wtime()
  
    if(time - timeLast > 10d0 .or. (.not.loop)) then
      timeLast = time
      do i=1, first_det_of_teeth(1)-1
        if(.not.(actually_computed(i))) then
          cycle pullLoop
        end if
      end do
      
      integer, external :: zmq_abort

      if (firstTBDcomb > Ncomb) then
        if (zmq_abort(zmq_to_qp_run_socket) == -1) then
          call sleep(1)
          if (zmq_abort(zmq_to_qp_run_socket) == -1) then
            print *, irp_here, ': Error in sending abort signal (1)'
          endif
        endif
        exit pullLoop
      endif

      double precision :: E0, avg, prop
      call do_carlo(tbc, Ncomb+1-firstTBDcomb, comb(firstTBDcomb), pt2_detail, actually_computed, sumabove, sum2above, Nabove)
      firstTBDcomb = int(Nabove(1)) - orgTBDcomb + 1
      if(Nabove(1) < 5d0) cycle
      call get_first_tooth(actually_computed, tooth)
     
      E0 = sum(pt2_detail(pt2_stoch_istate,:first_det_of_teeth(tooth)-1))
      if (tooth <= comb_teeth) then
        prop = ((1d0 - dfloat(comb_teeth - tooth + 1) * comb_step) - pt2_cweight(first_det_of_teeth(tooth)-1))
        prop = prop * pt2_weight_inv(first_det_of_teeth(tooth))
        E0 += pt2_detail(pt2_stoch_istate,first_det_of_teeth(tooth)) * prop
        avg = E0 + (sumabove(tooth) / Nabove(tooth))
        eqt = sqrt(1d0 / (Nabove(tooth)-1) * abs(sum2above(tooth) / Nabove(tooth) - (sumabove(tooth)/Nabove(tooth))**2))
      else
        eqt = 0.d0
      endif
      call wall_time(time)
      if ( (dabs(eqt/avg) < relative_error) .or. (dabs(eqt) < absolute_error) ) then
        ! Termination
        pt2(pt2_stoch_istate) = avg
        error(pt2_stoch_istate) = eqt
        print '(G10.3, 2X, F16.10, 2X, G16.3, 2X, F16.4, A20)', Nabove(tooth), avg+E, eqt, time-time0, ''
        if (zmq_abort(zmq_to_qp_run_socket) == -1) then
          call sleep(1)
          if (zmq_abort(zmq_to_qp_run_socket) == -1) then
            print *, irp_here, ': Error in sending abort signal (2)'
          endif
        endif
      else
        if (Nabove(tooth) > Nabove_old) then
          print '(G10.3, 2X, F16.10, 2X, G16.3, 2X, F16.4, A20)', Nabove(tooth), avg+E, eqt, time-time0, ''
          Nabove_old = Nabove(tooth)
        endif
      endif
    end if
  end do pullLoop
  
  E0 = sum(pt2_detail(pt2_stoch_istate,:first_det_of_teeth(tooth)-1))
  prop = ((1d0 - dfloat(comb_teeth - tooth + 1) * comb_step) - pt2_cweight(first_det_of_teeth(tooth)-1))
  prop = prop * pt2_weight_inv(first_det_of_teeth(tooth))
  E0 += pt2_detail(pt2_stoch_istate,first_det_of_teeth(tooth)) * prop
  pt2(pt2_stoch_istate) = E0 + (sumabove(tooth) / Nabove(tooth))

  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
  call sort_selection_buffer(b)
end subroutine

integer function pt2_find(v, w, sze, imin, imax)
  implicit none
  integer, intent(in) :: sze, imin, imax
  double precision, intent(in) :: v, w(sze)
  integer :: i,l,h
  integer, parameter :: block=64

  l = imin
  h = imax-1

  do while(h-l >= block)
    i = ishft(h+l,-1)
    if(w(i+1) > v) then
      h = i-1
    else
      l = i+1
    end if
  end do
  !DIR$ LOOP COUNT (64)
  do pt2_find=l,h
    if(w(pt2_find) >= v) then
      exit
    end if
  end do
end function


BEGIN_PROVIDER [ integer, comb_teeth ]
  implicit none
  comb_teeth = 100
END_PROVIDER



subroutine get_first_tooth(computed, first_teeth)
  implicit none
  logical, intent(in) :: computed(N_det_generators)
  integer, intent(out) :: first_teeth
  integer :: i, first_det

  first_det = 1
  first_teeth = 1
  do i=first_det_of_comb, N_det_generators
    if(.not.(computed(i))) then
      first_det = i
      exit
    end if
  end do
  
  do i=comb_teeth, 1, -1
    if(first_det_of_teeth(i) < first_det) then
      first_teeth = i
      exit
    end if
  end do

end subroutine


BEGIN_PROVIDER [ integer*8, size_tbc ]
  implicit none
  BEGIN_DOC
! Size of the tbc array
  END_DOC
  size_tbc = int((comb_teeth+1),8)*int(N_det_generators,8) + fragment_count*fragment_first
END_PROVIDER

subroutine get_carlo_workbatch(computed, comb, Ncomb, tbc)
  implicit none
  integer, intent(inout)         :: Ncomb
  double precision, intent(out)  :: comb(Ncomb)
  integer, intent(inout)         :: tbc(0:size_tbc)
  logical, intent(inout)         :: computed(N_det_generators)
  integer                        :: i, j, last_full, dets(comb_teeth)
  integer                        :: icount, n
  integer                        :: k, l
  l=first_det_of_comb
  call RANDOM_NUMBER(comb)
  do i=1,size(comb)
    comb(i) = comb(i) * comb_step
    !DIR$ FORCEINLINE
    call add_comb(comb(i), computed, tbc, size_tbc, comb_teeth)
    Ncomb = i
    if (tbc(0) == N_det_generators) return
    do while (computed(l))
      l=l+1
    enddo
    k=tbc(0)+1
    tbc(k) = l
    computed(l) = .True.
    tbc(0) = k
  enddo
  
end subroutine



subroutine get_comb(stato, dets, ct)
  implicit none
  integer, intent(in) :: ct
  double precision, intent(in) :: stato
  integer, intent(out) :: dets(ct)
  double precision :: curs
  integer :: j
  integer, external :: pt2_find

  curs = 1d0 - stato
  do j = comb_teeth, 1, -1
    !DIR$ FORCEINLINE
    dets(j) = pt2_find(curs, pt2_cweight,size(pt2_cweight), first_det_of_teeth(j), first_det_of_teeth(j+1))
    curs -= comb_step
  end do
end subroutine


subroutine add_comb(comb, computed, tbc, stbc, ct)
  implicit none
  integer*8, intent(in) :: stbc
  integer, intent(in) :: ct
  double precision, intent(in) :: comb
  logical, intent(inout) :: computed(N_det_generators)
  integer, intent(inout) :: tbc(0:stbc)
  integer :: i, k, l, dets(ct)
  
  !DIR$ FORCEINLINE
  call get_comb(comb, dets, ct)
  
  k=tbc(0)+1
  do i = 1, ct
    l = dets(i)
    if(.not.(computed(l))) then
      tbc(k) = l
      k = k+1
      computed(l) = .true.
    end if
  end do
  tbc(0) = k-1
end subroutine


BEGIN_PROVIDER [ integer, pt2_stoch_istate ]
 implicit none
 BEGIN_DOC
 ! State for stochatsic PT2 
 END_DOC
 pt2_stoch_istate = 1
END_PROVIDER


 BEGIN_PROVIDER [ double precision, pt2_weight, (N_det_generators) ]
&BEGIN_PROVIDER [ double precision, pt2_cweight, (N_det_generators) ]
&BEGIN_PROVIDER [ double precision, pt2_cweight_cache, (N_det_generators) ]
&BEGIN_PROVIDER [ double precision, comb_step ]
&BEGIN_PROVIDER [ integer, first_det_of_teeth, (comb_teeth+1) ]
&BEGIN_PROVIDER [ integer, first_det_of_comb ]
  implicit none
  integer :: i
  double precision :: norm_left, stato
  integer, external :: pt2_find  

  pt2_weight(1)  = psi_coef_generators(1,pt2_stoch_istate)**2
  pt2_cweight(1) = psi_coef_generators(1,pt2_stoch_istate)**2
  
  do i=1,N_det_generators
    pt2_weight(i) = psi_coef_generators(i,pt2_stoch_istate)**2
  enddo

  ! Important to loop backwards for numerical precision
  pt2_cweight(N_det_generators) = pt2_weight(N_det_generators)
  do i=N_det_generators-1,1,-1
    pt2_cweight(i) = pt2_weight(i) + pt2_cweight(i+1) 
  end do
  
  do i=1,N_det_generators
    pt2_weight(i)  = pt2_weight(i) / pt2_cweight(1)
    pt2_cweight(i) = pt2_cweight(i) / pt2_cweight(1)
  enddo

  do i=1,N_det_generators-1
    pt2_cweight(i) = 1.d0 - pt2_cweight(i+1) 
  end do
  pt2_cweight(N_det_generators) = 1.d0
  
  norm_left = 1d0
  
  comb_step = 1d0/dfloat(comb_teeth)
  first_det_of_comb = 1
  do i=1,N_det_generators
    if(pt2_weight(i)/norm_left < .5d0*comb_step) then
      first_det_of_comb = i
      exit
    end if
    norm_left -= pt2_weight(i)
  end do
  first_det_of_comb = max(2,first_det_of_comb)
  call write_int(6, first_det_of_comb-1, 'Size of deterministic set')
  
  comb_step =  (1d0 - pt2_cweight(first_det_of_comb-1)) * comb_step
  
  stato = 1d0 - comb_step
  iloc = N_det_generators
  do i=comb_teeth, 1, -1
    integer :: iloc
    iloc = pt2_find(stato, pt2_cweight, N_det_generators, 1, iloc)
    first_det_of_teeth(i) = iloc
    stato -= comb_step
  end do
  first_det_of_teeth(comb_teeth+1) = N_det_generators + 1
  first_det_of_teeth(1) = first_det_of_comb
  if(first_det_of_teeth(1) /= first_det_of_comb) then
     print *, 'Error in ', irp_here
     stop "comb provider"
  endif
  
END_PROVIDER

BEGIN_PROVIDER [ double precision, pt2_weight_inv, (N_det_generators) ]
  implicit none
  BEGIN_DOC
!  Inverse of pt2_weight array
  END_DOC
  integer :: i
  do i=1,N_det_generators
    pt2_weight_inv(i) = 1.d0/pt2_weight(i)
  enddo

END_PROVIDER






