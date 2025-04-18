module hfmod
implicit none

contains


subroutine read_matrix_overlap(mtrx,n,m,filename)
implicit none
character(len=*) :: filename
real(kind=8), dimension(n,n),intent(out) :: mtrx
integer :: i, j,row, io
integer, intent(in) :: n,m
real(kind=8) :: val


open(unit=1,file=filename)

!Matrix X:
    do row=1,n*m
        read(1,*,iostat=io) val,i,j
         !mtrx(i,j)=val
            if(io/=0) exit
           ! end if
                mtrx(i,j)=val
                mtrx(j,i)=val
           ! end if
    end do

close(unit=1)

end subroutine

subroutine read_overlap(overlap)
implicit none
real(kind=8), dimension(:,:), intent(inout) :: overlap
character(len=15) :: filename


filename='overlap'
call read_matrix_overlap(overlap,7,7,filename)
!write(*,*)'Overlap:', overlap

end subroutine


subroutine read_matrix_oneint(mtrx,n,m,filename)
implicit none
character(len=*) :: filename
real(kind=8), dimension(n,m),intent(out) :: mtrx
integer :: i, j,row, io
integer :: n, m
real(kind=8) :: val


open(unit=1,file=filename)

!Matrix X:
    do row=1,n*m
        read(1,*,iostat=io) val,i,j
         !mtrx(i,j)=val
            if(io/=0) exit
                mtrx(i,j)=val
                mtrx(j,i)=val
            !end if
    end do

close(unit=1)
end subroutine read_matrix_oneint


subroutine read_oneint(oneint)
implicit none 
real(kind=8), dimension(:,:), intent(inout) :: oneint
character(len=10) :: filename

filename='oneint'
call read_matrix_oneint(oneint,7,7,filename) !oben
!write(*,*)'Oneint:', oneint

end subroutine



SUBROUTINE read_2int(m,twoint,filename)
implicit none
integer :: u,o,v,p,row !Indice of twoint
INTEGER, INTENT(IN) :: m 
CHARACTER (len=10), INTENT(IN) :: filename 
REAL(kind=8),DIMENSION (m,m,m,m), INTENT (OUT) :: twoint
real(kind=8) :: val
integer :: io


OPEN (UNIT=1, FILE=filename)
	DO row = 1,m*m*m*m
			!write(*,*) "Rows twoint: ", row
			read (1,*,iostat=io) val, u, v, o, P
            if(io/=0) exit
			write(*,*) u,o,v,p, val
     		twoint(u,o,v,p)=val
				twoint(v,o,u,p)=val
				twoint(u,p,v,o)=val
				twoint(v,p,u,o)=val
				twoint(o,u,p,v)=val
				twoint(o,v,p,u)=val
				twoint(p,u,o,v)=val
				twoint(p,v,o,u)=val
end do
CLOSE (UNIT=1)
END SUBROUTINE

subroutine read_twoint(twoint)
implicit none 
real(kind=8), dimension(:,:,:,:), intent(inout) :: twoint
character(len=10) :: filename

filename='twoint'
call read_2int(7,twoint,filename)
!write(*,*)'Twoint:', twoint

end subroutine









subroutine read_repuls(repuls)
implicit none
character(len=10) :: filename
real(kind=8), intent(out) :: repuls

filename='repuls'

open(unit=1,file=filename)

		read(1,*) repuls

close(unit=1)

end subroutine read_repuls



!Λ = U†AU (subroutine for this is in another place with the name diag_sub.f90)

subroutine SD_inv_sqrt(n,X,sqrt_inv_SD)
	implicit none
	integer, intent(in) :: n
	integer :: i
	real(kind=8),dimension(n,n), intent(inout) :: X
	real(kind=8), dimension(n,n), intent(out) ::sqrt_inv_SD

sqrt_inv_SD = 0.0d0
!n=2
	do i=1,n
     if (X(i,i) > 0.0d0) then
 		sqrt_inv_SD(i,i) = 1.0d0/SQRT(X(i,i))
     else
        sqrt_inv_SD(i,i) = 0.0d0  !avoid division by zero
     end if
	end do
write(*,*) 'sqrt_inv_SD=',sqrt_inv_SD
end subroutine


subroutine density_matrix(D,C,nocc,n)
	implicit none
  integer, intent(in) :: nocc, n
  real(kind=8), dimension(n,n), intent(in) :: C
  real(kind=8), dimension(n,n), intent(out) :: D
  integer :: i, j,k
 
D= 0.0d0

  do i = 1, n
     do j = 1,n
				do k=1,nocc
			D(i,j)= D(i,j)+ 2.0d0*(C(i,k)*C(j,k))
			end do      
     end do
  end do

end subroutine



subroutine supermatrix(twoint,nmax,P)
	implicit none
	integer, intent(in) :: nmax
	integer :: mu,rho,nu,sigma  !(μ,ρ,ν und σ): mu,rho,nu,sigma
	real(kind=8), dimension(nmax,nmax,nmax,nmax), intent(in) :: twoint
	real(kind=8), dimension(nmax,nmax,nmax,nmax), intent(out) :: P

P=0.0d0 
	do mu=1,nmax
		do rho = 1,nmax
			do nu = 1,nmax
				do sigma = 1,nmax
P(mu,nu,sigma,rho) = twoint(mu,sigma,nu,rho) - 0.25 * twoint(mu,sigma,rho,nu) - 0.25 * twoint(mu,rho,sigma,nu)
				end do
			end do
		end do
	end do



end subroutine



subroutine fockmatrix(oneint,D,P,F,n)

  implicit none
  integer, intent(in) :: n
  real(kind=8), dimension(n,n), intent(in) :: oneint, D
  real(kind=8), dimension(n,n,n,n), intent(in) :: P
  real(kind=8), dimension(n,n), intent(out) :: F
  integer :: mu, nu, sigma, rho

  do mu = 1, n
     do nu = 1, n
  F(mu,nu) = oneint(mu,nu)
				do sigma=1,n
					do rho=1,n
        F(mu,nu) = F(mu,nu) + (P(mu,nu,sigma,rho) * D(sigma,rho))
					end do
				end do
     end do
  end do

end subroutine

subroutine Fprime(F,sqrt_inv_S,FS,F_prime,nbf)
implicit none
integer, intent(in) :: nbf
real(kind=8), dimension(nbf,nbf), intent(in) :: F,sqrt_inv_S
real(kind=8), dimension(nbf,nbf), intent(out) :: FS,F_prime


!F′ = S−1/2FS−1/2
FS=0.d0

!call DGEMM('N','N',M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
call dgemm('N','N',nbf,nbf,nbf,1.0d0,F,nbf,sqrt_inv_S,nbf,0.0d0,FS,nbf)
!call matrixmulti(F,sqrt_inv_S,FS,2,2,2)
write(*,*) 'FS=',FS

F_prime=0.d0
call dgemm('N','N',nbf,nbf,nbf,1.0d0,sqrt_inv_S,nbf,FS,nbf,0.0d0,F_prime,nbf)
!call matrixmulti(sqrt_inv_S,FS,F_prime,2,2,2)
write(*,*) 'F′=',F_prime


end subroutine

subroutine diag_fockmatrix(F_prime,F_diag,c_prime,n)
	implicit none
	integer,intent(in) :: n
	real(kind=8), dimension(n,n), intent(in) :: F_prime
	real(kind=8), dimension(n,n), intent(out) :: c_prime, F_diag

 
!call EIG(F_prime,F_diag,c_prime)
F_diag=F_prime
call eig(F_diag,c_prime,0,n,0)
end subroutine




subroutine back_trans_C(C, c_prime,sqrt_inv_S ,n)
	implicit none 
	integer, intent(in):: n
	real(kind=8), dimension (n,n), intent(in) :: c_prime,sqrt_inv_S
	real(kind=8), dimension(n,n), intent(out) :: C

!call matrixmulti(sqrt_inv_S,c_prime,C,n,n,n)
call dgemm('N','N',n,n,n,1.0d0,sqrt_inv_S,n,c_prime,n,0.0d0,C,n)

end subroutine 




subroutine HF_Energie(E_HF,D,oneint,F,repuls,n)
implicit none
integer, intent(in) :: n
real(kind=8), dimension(n,n), intent(in):: D,oneint,F
real(kind=8),intent(in):: repuls
real(kind=8), intent(out) :: E_HF
integer :: mu, nu

E_HF=0.d0

	do mu=1,n
		do nu=1, n
			E_HF = E_HF + 0.5d0*(D(mu, nu)*(oneint (mu, nu) +F (mu, nu)))
		end do
	end do

E_HF = E_HF + repuls

end subroutine


subroutine check_conv(E_HF, old_EHF, conv, max_i, i, cont)
  implicit none
  real(kind=8), intent(in) :: E_HF, old_EHF, conv 
  integer, intent(in) :: max_i, i
  logical, intent(out) :: cont

  cont = .true.

  ! Check if convergence criterion is met
  if (abs(E_HF - old_EHF) < conv) then
     cont = .false.
     write(*,*) 'Convergence reached in', i, "iterations"
  elseif (i == max_i) then
     write(*,*) 'Maximum number of iterations reached'
     cont = .false.
  end if

end subroutine



SUBROUTINE update_molden_output (nbf, nocc,A,B,filename )
IMPLICIT NONE
INTEGER, INTENT (IN) :: nbf, nocc
INTEGER :: i, j
CHARACTER (*), INTENT (IN) :: filename
REAL (kind=8), DIMENSION(nbf,nbf), INTENT (IN) :: A, B

OPEN(UNIT=1, FILE=filename, ACCESS='append', FORM='formatted', STATUS='old' )

    DO i = 1, nbf
        WRITE (1, ' (a, f15.10) ') "Ene= ", A(i,i)
        WRITE (1, ' (a, f7.5) ') "Spin= Alpha"
            if (i<=nocc) then
        WRITE (1, ' (a, f7.5) ') "Occup=", 2.d0
            else
        WRITE (1, ' (a, f7.5) ') "Occup=", 0.d0
            end if
    DO j= 1, nbf
        WRITE (1, ' (i0, f15.10) ' ) j, B(j,i)
    END DO
    END DO
        WRITE (1,*) "end"

CLOSE (UNIT=1)

END SUBROUTINE


subroutine atomic_charge(n,nocc,D,overlap)
IMPLICIT NONE 
INTEGER :: i
INTEGER,INTENT(IN):: n, nocc
REAL(Kind=8), DIMENSION(n,n), INTENT(IN) ::D,overlap
REAL(Kind=8), DIMENSION(n,n) :: n_A, DS
REAL(Kind=8) :: Z_O, Z_H1, Z_H2, n_H1, n_H2,q_H1,q_H2, q_O, val

DS= 0.0
Z_O = 8 !Kernladungszahl of oxygen
Z_H1 = 1 !Kernladungszahl of Hydrogen
Z_H2 = 1
n_A = 0

CALL matrixmulti(D,overlap,DS,n,n,n)
	do i = 1, n
n_A(i,i) = DS(i,i)
	end do

!Wasserstoff(Hydrogen):
n_H1 = n_A (1,1)
q_H1 = Z_H1 - n_H1
WRITE(*,*) "The partial charge of the first hydrogen is: "
WRITE (*, ' (f15.10)') q_H1

n_H2 = n_A (7,7)
q_H2 = Z_H2 - n_H2
WRITE(*,*) "The partial charge of the second hydrogen is: "
WRITE (*,'(f15.10)') q_H2

!Sauerstoff (oxygen)
	do i = 2,6
		val = val + n_A(i,i)
	end do
q_O = Z_O - val
WRITE (*,*) "Die Partialladung von Sauerstoff ist: "
WRITE (*, ' (f15.10) ') q_O

end subroutine

subroutine mp2_mo(n,C,twoint,MO)
implicit none
real(kind=8), dimension (n,n,n,n),intent(in) :: twoint
real(kind=8), dimension (n,n,n,n),intent(out):: MO
real(kind=8), dimension (n,n),intent(in):: C
integer:: p,q,r,s,a,b,i,j
integer:: mu,nu,sigma,rho
integer,intent(in):: n

do p=1,n
	do q=1,n
		do r=1,n
			do s=1,n
				do mu= 1,n
					do nu=1,n
						do sigma=1,n
							do rho=1,n
					MO(p,q,r,s)=MO(p,q,r,s)+((twoint(mu,nu,sigma,rho) *C(mu, p)*C(nu,q)*C(sigma,r)*C(rho,s)))
							end do
						end do
					end do
				end do
			end do
		end do
	end do
end do


end subroutine


subroutine mp2_energy(n,nocc,MO,e,E_MP2)
implicit none
real(kind=8), dimension (n,n,n,n) :: MO
real(kind=8), dimension (n,n),intent(in)::e
real(kind=8),intent(out) :: E_MP2
integer:: i,j,a,b
integer,intent(in):: n, nocc

E_MP2=0.0d0
do i=1,nocc
	do j=1,nocc
		do a=nocc+1,n
			do b=nocc+1,n
				E_MP2=E_MP2+(((2*MO(i,j,a,b)-MO(i,j,b,a)) *MO(i,j,a,b))/ (e(i,i)+e(j,j)-e(a,a)-e(b,b)))
			end do
		end do
	end do
end do


end subroutine

subroutine MOs(n,C,twoint,MO)
implicit none
integer,intent(in)::n
real(kind=8), dimension (n,n,n,n),intent(in):: twoint
real(kind=8), dimension (n,n,n,n),intent(out) ::MO
real(kind=8), dimension (n,n,n,n)::MO1,MO2,MO3
real(kind=8), dimension (n,n),intent(in):: C
integer:: p,q,r,s
integer:: mu,nu,sigma,rho


MO1=0.0d0
MO2=0.0d0
MO3=0.0d0
MO=0.0d0

do mu=1,n
	do nu=1, n
		do sigma=1,n
			do s=1,n
				do rho=1,n
						MO1(mu,nu,sigma,s) =MO1(mu,nu,sigma,s)+(twoint(mu,nu,sigma,rho) *C(rho, s))
				end do
			end do
		end do
	end do
end do

do mu=1,n
	do nu=1, n
		do r=1,n
			do s=1, n
				do sigma=1,n
					MO2(mu,nu,r,s) =MO2(mu, nu,r,s)+(MO1(mu,nu,sigma,s) *C(sigma, r))
				end do
			end do
		end do
	end do
end do

do mu=1,n
	do q=1,n
		do r=1,n
			do s=1,n
				do nu=1,n
			MO3(mu,q,r,s)=MO3(mu,q,r,s)+(MO2(mu, nu,r,s) *C(nu, q))
				end do
			end do
		end do
	end do
end do


do p=1,n
	do q=1,n
		do r=1,n
			do s=1,n
				do mu= 1,n
					MO(p,q,r,s)=MO(p,q,r,s)+(MO3(mu,q,r,s) *C(mu, p))
				end do
			end do
		end do
	end do
end do

!write(*,*) 'MO1 :', MO1
!write(*,*) 'MO2 :', MO2
!write(*,*) 'MO3 :', MO3
!write(*,*) 'MO :', MO
end subroutine

end
