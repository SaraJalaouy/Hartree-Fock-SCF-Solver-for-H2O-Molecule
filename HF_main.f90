program hf
use hfmod
implicit none
real(kind=8), dimension(:,:), allocatable :: overlap,oneint,SD,U,sqrt_inv_SD,sqrt_inv_S, D,C,F,UT,USD,FS,F_prime,c_prime, F_diag
real(kind=8),dimension(:,:,:,:),allocatable :: twoint,P,MO
real(kind=8) :: repuls,E_HF,old_EHF,conv,E_MP2
character(len=15) :: filename
integer :: n, nbf,nocc,max_i,i,j, k, l
logical :: cont


!n=1
nbf=7
nocc=5


allocate(overlap(nbf,nbf))
allocate(oneint(nbf,nbf))
allocate(twoint(nbf,nbf,nbf,nbf))
allocate(SD(nbf,nbf))
allocate(U(nbf,nbf))
allocate(sqrt_inv_SD(nbf,nbf))
allocate(sqrt_inv_S(nbf,nbf))
allocate(D(nbf,nbf))
allocate(C(nbf,nbf))
allocate(P(nbf,nbf,nbf,nbf))
allocate(F(nbf,nbf))
allocate(UT(nbf,nbf))
allocate(USD(nbf,nbf))
allocate(FS(nbf,nbf))
allocate(F_prime(nbf,nbf))
allocate(c_prime(nbf,nbf))
allocate(F_diag(nbf,nbf))
allocate(MO(nbf,nbf,nbf,nbf))


write(*,*) 'Overlap Matrix:'
call read_overlap(overlap)

write(*,*) 'One electron Integral :'
call read_oneint(oneint)

write(*,*) 'Two electron Integral :'
call read_twoint(twoint)

write(*,*) 'Repuls Matrix:'
call read_repuls(repuls)
write(*,*) repuls

write(*,*) 'Integrals read successfully!'

write(*,*) 'finde SD = U†SU?'
!call EIG(overlap,SD,U)
SD=overlap
call eig(SD,U,0,nbf,0)
write(*,*)'SD=', SD


call SD_inv_sqrt(nbf,SD,sqrt_inv_SD)
write(*,*) 'SD^-1/2=', sqrt_inv_SD



write(*,*) 'S−1/2 = USD−1/2U†'
call dgemm('N','N',nbf,nbf,nbf,1.0d0,U,nbf,sqrt_inv_SD,nbf,0.0d0,USD,nbf)
!call matrixmulti(U,sqrt_inv_SD,USD,nbf,nbf,nbf)
write(*,*) 'USD−1/2', USD


UT=transpose(U)
call dgemm('N','N',nbf,nbf,nbf,1.0d0,USD,nbf,UT,nbf,0.0d0,sqrt_inv_S,nbf)
!call matrixmulti(USD,UT,sqrt_inv_S,nbf,nbf,nbf)
write(*,*) 'S^-1/2=', sqrt_inv_S


write(*,*) 'Dichtematrix..'
C=0.d0
call density_matrix(D,C,nocc,nbf)
write(*,*) D

call supermatrix(twoint,nbf,P)
write(*,*) 'P='
DO i =1, nbf
    DO j =1, nbf
        DO k =1, nbf
            DO l =1, nbf
write(*,*) i, j, k, l, P(i,j,k,l)
            End do
        end do
    end do
end do



write(*,*) 'Setting up the Fock matrix..'
call fockmatrix(oneint,D,P,F,nbf)
write(*,*) 'F=',F


write(*,*) 'Transformation of the Fock matrix into the orthogonal basis...(F′=S−1/2FS−1/2)'
call Fprime(F,sqrt_inv_S,FS,F_prime,nbf)

!(F′C′ = C′ε)

write(*,*) 'Diagonalization of F_prime...'
call diag_fockmatrix(F_prime,F_diag,c_prime,nbf)
write(*,*) 'Orbital enr: F_diag=', F_diag


write(*,*) 'Back transformation of the C′ into the original AO basis..'
call back_trans_C(C, c_prime,sqrt_inv_S ,nbf)
write(*,*) 'C = S−1/2C =', C


call density_matrix(D,C,nocc,nbf)
write(*,*) 'Dichtematrix : D=',D

call fockmatrix(oneint,D,P,F,nbf) 
write(*,*) 'New F: F=',F !, 'P=',P, 'ONEINT:',oneint

!call Fprime(F,sqrt_inv_S,FS,F_prime)
!call diag_fockmatrix(F_prime,F_diag,c_prime,nbf)


write(*,*) 'Hartree-Fock Energie:'
call HF_Energie(E_HF,D,oneint,F,repuls,nbf)
write(*,*) 'E_HF=',E_HF

WRITE(*,*) 'F_DIAG=', F_diag



! Start the SCF loop
  i = 0
	max_i=100
  old_EHF=0.d0
	conv=1.d-8
	do !while (cont .and. i < max_i)
     !i = i + 1

    ! call density_matrix(D, C, nocc, nbf)
     call Fprime(F, sqrt_inv_S, FS, F_prime,nbf)
     call diag_fockmatrix(F_prime, F_diag, c_prime, nbf)
     call back_trans_C(C, c_prime,sqrt_inv_S ,nbf) 
		 call density_matrix(D, C, nocc, nbf)
     call fockmatrix(oneint, D, P, F, nbf)
		 call HF_Energie(E_HF, D, oneint, F, repuls, nbf)
     call check_conv(E_HF, old_EHF, conv, max_i, i, cont)
old_EHF = E_HF
i=i+1
     old_EHF = E_HF

     ! If converged, exit the loop
     if (.not. cont) then  ! Exit loop if converged
      write(*,*) 'Converged at iteration:', i
      exit
   end if 
!i=i+1 
  end do

if (i == max_i) then
   write(*,*) 'Maximum iterations reached. Final energy:', E_HF
else
   write(*,*) 'Final Hartree-Fock energy:', E_HF
end if
!write(*,*)'E_HF:',E_HF 


!filename ='wasser.out'
!Call the update_hazwei subroutine to update hazwei.out with new values
!call update_molden_output (nbf,nocc,F_diag,C, filename)


!call atomic_charge(nbf,nocc,D,overlap)


!MP2:


 
!call mp2_mo(nbf,C,twoint,MO)
call MOs(nbf,C,twoint,MO)
call mp2_energy(nbf,nocc,MO,F_diag,E_MP2)
print*, 'MP2 Energy Correction: ', E_MP2


deallocate(overlap)
deallocate(oneint)
deallocate(twoint)
deallocate(SD)
deallocate(U)
deallocate(sqrt_inv_SD)
deallocate(sqrt_inv_S)
deallocate(D)
deallocate(C)
deallocate(P)
deallocate(F)
deallocate(USD)
deallocate(UT)
deallocate(FS)
deallocate(F_prime)
deallocate(c_prime)
deallocate(F_diag)
deallocate(MO)

end program

