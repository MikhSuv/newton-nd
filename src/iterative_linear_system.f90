module iterative_linear_system
   use precision_mod
   implicit none
    private
   integer, PARAMETER :: max_iter = 100000000
   real(dp), PARAMETER :: eps = 1.0e2_dp * epsilon(1.0_dp)
    public :: read_linear_system, write_result, iterative_solve, residual
contains
   subroutine read_linear_system(filename, A, B, n)
      character(len=*), intent(in) :: filename
      real(dp), allocatable, intent(out) :: A(:, :), B(:) ! СЛУ (A|B)
      integer, intent(out), optional :: n

      integer :: iunit, iostatus, i
      character(len=256) :: line ! для чтения первой строки

      open (newunit=iunit, file=filename, status='old', &
            action='read', iostat=iostatus)
      if (iostatus /= 0) then
         error stop 'Error occured while opening file'
      end if

      read (iunit, '(a)', iostat=iostatus) line

      if (iostatus /= 0) then
         error stop 'Error occured while reading line'
      end if
      read (line(2:), *) n

      ! Выделение памяти под матрицу A
      allocate (A(n, n))
      ! Выделение памяти под столбец B
      allocate (B(n))
      do i = 1, n
         read (iunit, *) A(i, :)
      end do

      do i = 1, n
         read (iunit, *) B(i)
      end do

      close (iunit)

   end subroutine read_linear_system

   subroutine write_result(filename, X)
      character(len=*), intent(in) :: filename
      real(dp), intent(in) ::  X(:) ! Столбец результата

      integer :: n
      integer :: ounit, i, iostatus

      n = size(X)
      open (newunit=ounit, file=filename, action='write', iostat=iostatus)
      if (iostatus /= 0) then
         error stop 'Error occured while opening file'
      end if

      write (ounit, '("# ", i0)') n

      do i = 1, n
         write (ounit, '(*(e23.15, 1x))') X(i)
      end do

      close (ounit)

   end subroutine write_result
   function is_diagonally_dominant(A) result(check)
      ! Проверяет, обладает ли матрица A диагоняльным преобладанием
      real(dp), intent(in) :: A(:, :) ! Матрица
      real(dp), allocatable :: diag_vals(:), row_summs(:)
      logical :: check
      integer :: n, i
      check = .false.
      n = size(A(:, 1))

      allocate (diag_vals(n), row_summs(n))
      diag_vals = [(abs(A(i, i)), i=1, n)]
      row_summs = sum(abs(A), dim=2) - diag_vals
      check = all(diag_vals > row_summs)
   end function is_diagonally_dominant

   function iterative_solve(A, B, method) result(X)
      character(len=*), intent(in) :: method
      real(dp), intent(in) :: A(:, :), B(:)
      real(dp), allocatable :: X(:)

      select case (method)
      case ("jacobi")
         X = jacobi_solve(A, B)
      case ("seidel")
         X = seidel_solve(A, B)
      case ("relaxation")
         X = relaxation_solve(A, B)
      case default
         error stop "Invalid method parameter"
      end select

   end function iterative_solve

   function jacobi_solve(A, B) result(X)
      real(dp), intent(in) :: A(:, :), B(:)
      real(dp), allocatable :: X(:), X_new(:)
      real(dp), allocatable :: invD(:), Z(:, :), G(:)
      integer :: n, i

      if (.not. is_diagonally_dominant(A)) then
         print *, "WARDING: Not diagonally dominant matrix!"
      end if
      n = size(B)
      allocate (invD(n), G(n), Z(n, n))
      ! Обратная диагональ
      invD = 1.0_dp/[(A(i, i), i=1, n)]
      G = B*invD
      ! Начальное приближение
      X = G

      Z = A
      !$omp parallel do private(i) shared(A, invD, Z)
      do i = 1, n
         Z(i, :) = -invD(i)*A(i, :)
      end do
      !$omp end parallel do
      forall (i=1:n) Z(i, i) = 0.0_dp

      do i = 1, max_iter
         X_new = matmul(Z, X) + G
         if (norm2(X_new - X) < eps) then
            return
         end if
         X = X_new
      end do
      print *, "Ineration limit has been reached"
   end function jacobi_solve

   function seidel_solve(A, B) result(X)
      real(dp), intent(in) :: A(:, :), B(:)
      real(dp), allocatable :: X(:), X_old(:)
      real(dp), allocatable :: invD(:), Q(:), P(:, :)
      integer :: n, i, j

      if (.not. is_diagonally_dominant(A)) then
         print *, "WARDING: Not diagonally dominant matrix!"
      end if

      n = size(B)
      allocate (invD(n), Q(n), P(n, n), X_old(n))
      ! Обратная диагональ
      invD = 1.0_dp/[(A(i, i), i=1, n)]
      !$omp parallel do private(i) shared(A, invD, P)
      do i = 1, n
         P(i, :) = -invD(i)*A(i, :)
      end do
      !$omp end parallel do
      forall (i=1:n) P(i, i) = 0.0_dp

      Q = B * invD
      X = Q

      do j = 1, max_iter
         X_old = X 
         do i = 1, n
         X(i) = dot_product(P(i, 1:i-1), X(1:i-1)) &
           + dot_product(P(i, i+1:n), X_old(i+1:n)) &
           + Q(i)
         end do
         if (norm2(X_old - X) < eps) return
      end do
      print *, "Ineration limit has been reached"
   end function seidel_solve

   function relaxation_solve(A, B) result(X)
      real(dp), intent(in) :: A(:, :), B(:)
      real(dp), allocatable :: X(:)
      real(dp), allocatable :: invD(:), Q(:), P(:, :)
      integer :: n, i, j

      if (.not. is_diagonally_dominant(A)) then
         print *, "WARDING: Not diagonally dominant matrix!"
      end if

      n = size(B)
      allocate (invD(n), Q(n), P(n, n))

      ! Обратная диагональ
      invD = 1.0_dp/[(A(i, i), i=1, n)]
      P = A
      !$omp parallel do private(i) shared(A, invD, P)
      do i = 1, n
         P(i, :) = -invD(i)*A(i, :)
      end do
      !$omp end parallel do
      forall (i=1:n) P(i, i) = 0.0_dp

      Q = B*invD
      allocate (X(n))
      ! Начальное приближение
      X = 0.0_dp

      do i = 1, max_iter
         j = maxloc(abs(Q), dim = 1)

         if (abs(Q(j)) < eps) then
            return
         end if

         X(j) = X(j) + Q(j)
         Q = Q + P(:, j)*Q(j)
         Q(j) = 0.0_dp
      end do

      print *, "Ineration limit has been reached"
   end function relaxation_solve

   function residual(A, B, X) result(r)
      real(dp), intent(in) :: A(:, :), B(:), X(:)
      real(dp) :: r ! модуль вектора невязки

      r = norm2(matmul(A, X) - B)
   end function residual
end module iterative_linear_system
