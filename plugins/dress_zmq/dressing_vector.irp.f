
 BEGIN_PROVIDER [ double precision, dressing_column_h, (N_det,N_states) ]
&BEGIN_PROVIDER [ double precision, dressing_column_s, (N_det,N_states) ]
 implicit none
 BEGIN_DOC
 ! Null dressing vectors
 END_DOC
 dressing_column_h(:,:) = 0.d0
 dressing_column_s(:,:) = 0.d0

 integer :: i,ii,k,j,jj, l
 double precision :: f, tmp
 double precision, external :: u_dot_v
 
 do k=1,N_states
   l = dressed_column_idx(k)
   f = 1.d0/psi_coef(l,k)
   do jj = 1, n_det
     j = jj !idx_non_ref(jj)
     dressing_column_h(j,k) = delta_ij   (k,jj,1) * f
     dressing_column_s(j,k) = delta_ij   (k,jj,2) * f!delta_ij_s2(k,jj)
   enddo
   tmp = u_dot_v(dressing_column_h(1,k), psi_coef(1,k), N_det)
   dressing_column_h(l,k) -= tmp * f
   tmp = u_dot_v(dressing_column_s(1,k), psi_coef(1,k), N_det)
   dressing_column_s(l,k) -= tmp * f
 enddo
END_PROVIDER

