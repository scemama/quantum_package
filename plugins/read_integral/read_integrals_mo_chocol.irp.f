program read_integrals
  BEGIN_DOC
! Reads the integrals from the following files:
! - kinetic_mo
! - nuclear_mo
! - bielec_mo
  END_DOC

  integer :: iunit
  integer :: getunitandopen
  integer :: i,j,n

  PROVIDE ezfio_filename
  call ezfio_set_integrals_monoelec_disk_access_mo_one_integrals("None")

  logical :: has
  call ezfio_has_mo_basis_mo_tot_num(has)
  if (.not.has) then

    iunit = getunitandopen('nuclear_mo','r')
    n=0
    do 
      read (iunit,*,end=12) i
      n = max(n,i)
    enddo
    12 continue
    close(iunit)
    call ezfio_set_mo_basis_mo_tot_num(n)

    call ezfio_has_ao_basis_ao_num(has)
    mo_label = "None"
    if (has) then
      call huckel_guess
    else
      call ezfio_set_ao_basis_ao_num(n)
    endif
  endif
  call run
end

subroutine run
  use map_module
  implicit none
  
  integer :: iunit
  integer :: getunitandopen

  integer ::i,j,k,l
  double precision :: integral
  double precision, allocatable :: A(:,:)

  integer             :: n_integrals 
  integer(key_kind), allocatable   :: buffer_i(:) 
  real(integral_kind), allocatable :: buffer_values(:)
  integer(key_kind)  :: key
   
  call ezfio_get_mo_basis_mo_tot_num(mo_tot_num)

  allocate (A(mo_tot_num,mo_tot_num))
  A = 0.d0
  
  iunit = getunitandopen('kinetic_mo','r')
  do 
    read (iunit,*,end=10) i,j, integral
    A(i,j) = integral
  enddo
  10 continue
  close(iunit)
  call write_one_e_integrals('mo_kinetic_integral', A, size(A,1), size(A,2))


  iunit = getunitandopen('nuclear_mo','r')
  do 
    read (iunit,*,end=12) i,j, integral
    A(i,j) = integral
  enddo
  12 continue
  close(iunit)
  call write_one_e_integrals('mo_ne_integral', A, size(A,1), size(A,2))

  call write_one_e_integrals('mo_pseudo_integral', mo_pseudo_integral,&
        size(mo_pseudo_integral,1), size(mo_pseudo_integral,2))


  call ezfio_set_integrals_monoelec_disk_access_mo_one_integrals("Read")
end
