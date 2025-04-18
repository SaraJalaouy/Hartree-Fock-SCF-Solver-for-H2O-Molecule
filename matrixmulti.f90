subroutine matrixmulti(A, B, C, n, m, p)
  implicit none
  real(kind=8), dimension(n, m), intent(in) :: A   ! Matrix A (n x m)
  real(kind=8), dimension(m, p), intent(in) :: B   ! Matrix B (m x p)
  real(kind=8), dimension(n, p), intent(out) :: C  ! Matrix C (n x p)
  integer :: i, j, k,n,m,p

  C = 0.0d0

  ! Matrix multiplication: C = A * B
  do i = 1, n
    do j = 1, p
      do k = 1, m
        C(i, j) = C(i, j) + A(i, k) * B(k, j)
      end do
    end do
  end do

!  write (*,*) 'The product of the given matrices:',C
end subroutine





