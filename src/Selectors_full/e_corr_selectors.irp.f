
use bitmasks
 BEGIN_PROVIDER [integer, exc_degree_per_selectors, (N_det_selectors)]
&BEGIN_PROVIDER [integer, double_index_selectors, (N_det_selectors)]
&BEGIN_PROVIDER [integer, n_double_selectors]
 implicit none
 BEGIN_DOC
 ! degree of excitation respect to Hartree Fock for the wave function
 !
 ! for the all the selectors determinants
 !
 ! double_index_selectors = list of the index of the double excitations
 ! 
 ! n_double_selectors = number of double excitations in the selectors determinants
 END_DOC
 integer :: i,degree
 n_double_selectors = 0
 do i = 1, N_det_selectors
  call get_excitation_degree(psi_selectors(1,1,i),ref_bitmask,degree,N_int)
  exc_degree_per_selectors(i) = degree
  if(degree==2)then
   n_double_selectors += 1
   double_index_selectors(n_double_selectors) =i
  endif
 enddo
END_PROVIDER
  BEGIN_PROVIDER[double precision, coef_hf_selector]
 &BEGIN_PROVIDER[double precision, E_corr_per_selectors, (N_det_selectors)]
 implicit none
 BEGIN_DOC
 ! energy of correlation per determinant respect to the Hartree Fock determinant
 !
 ! for the all the double excitations in the selectors determinants
 !
 ! E_corr_per_selectors(i) = <D_i|H|HF> * c(D_i)/c(HF) if |D_i> is a double excitation
 ! 
 ! E_corr_per_selectors(i) = -1000.d0 if it is not a double excitation 
 !
 ! coef_hf_selector = coefficient of the Hartree Fock determinant in the selectors determinants
 END_DOC
 integer :: i,degree
 double precision :: hij,inv_selectors_coef_hf
 do i = 1, N_det_selectors
  if(exc_degree_per_selectors(i)==2)then
   call i_H_j(ref_bitmask,psi_selectors(1,1,i),N_int,hij)
   E_corr_per_selectors(i) = psi_selectors_coef(i,1) * hij
  elseif(exc_degree_per_selectors(i) == 0)then
   coef_hf_selector = psi_selectors_coef(i,1)
   E_corr_per_selectors(i) = -1000.d0
  else
   E_corr_per_selectors(i) = -1000.d0
  endif
 enddo
 if (dabs(coef_hf_selector) > 1.d-8) then
   inv_selectors_coef_hf = 1.d0/coef_hf_selector
 else
   inv_selectors_coef_hf = 0.d0
 endif
 do i = 1,n_double_selectors
  E_corr_per_selectors(double_index_selectors(i)) *=inv_selectors_coef_hf
 enddo

 END_PROVIDER
