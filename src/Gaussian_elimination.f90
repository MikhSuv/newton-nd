module Gaussian_elimination
   use precision_mod
   implicit none
   private

   public :: read_linear_system, write_result, solve_linear_system, residual
   real(dp), PARAMETER :: eps = 1.0e2_dp * epsilon(1.0_dp)

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

   function solve_linear_system(A, B, method) result(X)
     real(dp), intent(in) :: A(:,:), B(:)
     real(dp), allocatable :: X(:)
     character(len=*), intent(in) :: method

     select case (method)
      case("gauss")
        X = gaussian_solve(A, B)
      case("jordan")
        X = jordan_solve(A, B)
      case("choice")
        X = gaussian_choice(A, B)
      case default
        error stop "Invalid method parameter"
     end select

   end function solve_linear_system

   function gaussian_solve(AA, B) result(X)
      real(dp), intent(in) :: AA(:, :), B(:)
      real(dp), allocatable :: X(:) ! столбец результата
      real(dp), allocatable :: A(:, :) ! расширенная матрица системы

      real(dp) :: leading_entry ! ведуший элемент
      integer :: n
      integer :: i, k
      n = size(B)
      allocate (A(n, n + 1))
      A(:, 1:n) = AA
      A(:, n + 1) = B
      ! Прамой ход
      do k = 1, n
         leading_entry = a(k, k)
         if (abs(leading_entry) <= eps) then
            print *, "WARNING: the leading_entry in step", k, " is too small: ", leading_entry
         end if

         a(k, k:n + 1) = a(k, k:n + 1)/leading_entry

         ! do concurrent(i=k + 1:n)
         !$OMP PARALLEL DO PRIVATE(i) SHARED(a, k, n)
         do i=k + 1,n
            a(i, k:n + 1) = a(i, k:n + 1) - a(i, k)*a(k, k:n + 1)
         end do
         !$OMP END PARALLEL DO
      end do

      ! Обратный ход
      allocate (X(n))
      do i = n, 1, -1
         x(i) = a(i, n + 1) - dot_product(a(i, i + 1:n), x(i + 1:n))
      end do
   end function gaussian_solve

   function jordan_solve(AA, B) result(X)
      real(dp), intent(in) :: AA(:, :), B(:)
      real(dp), allocatable :: X(:) ! столбец результата
      real(dp), allocatable :: A(:, :) ! расширенная матрица системы

      real(dp) :: leading_entry
      integer :: n
      integer :: i, k
      n = size(B)
      allocate (A(n, n + 1))
      A(:, 1:n) = AA
      A(:, n + 1) = B
      ! Прамой ход
      do k = 1, n
         leading_entry = a(k, k)
         if (abs(leading_entry) <= eps) then
            print *, "WARNING: the leading_entry in step", k, " is too small: ", leading_entry
         end if

         a(k, k:n + 1) = a(k, k:n + 1)/leading_entry

         ! do concurrent(i=1:k-1)
         !$OMP PARALLEL DO PRIVATE(i) SHARED(a, k, n)
         do i=1,k-1
            a(i, k:n + 1) = a(i, k:n + 1) - a(i, k)*a(k, k:n + 1)
         end do
         !$OMP END PARALLEL DO

         ! do concurrent(i=k+1:n)
         !$OMP PARALLEL DO PRIVATE(i) SHARED(a, k, n)
         do i=k+1,n
            a(i, k:n + 1) = A(i, k:n + 1) - A(i, k)*A(k, k:n + 1)
         end do
         !$OMP END PARALLEL DO
      end do

      allocate (X(n))
      X = A(:, n + 1)
   end function jordan_solve

   function gaussian_choice(AA, B) result(X)
      real(dp), intent(in) :: AA(:, :), B(:)
      real(dp), allocatable :: X(:) ! столбец результата
      real(dp), allocatable :: A(:, :) ! расширенная матрица системы
      real(dp) ::factor
      integer :: n
      integer :: i, k, loc, lead_idx
      integer, allocatable :: perm(:)

      n = size(B)
      allocate(A(n, n+1))
      A(:, 1:n) = AA
      A(:, n + 1) = B
      perm = [(k, k=1,n)]      
      do k = 1, n-1
        lead_idx = maxloc(abs(A(perm(k:n), k)), dim = 1) + k - 1
        if (lead_idx /= k) then 
          perm([k, lead_idx]) =  perm([lead_idx, k])
        end if 

        if (abs(A(perm(k), k)) < eps) then 
          print *, "WARNING: the leading_entry in step", k, " is too small: ", A(perm(k), k)
        end if 

        ! do concurrent (i = k+1:n) local(factor)
        !$OMP PARALLEL DO PRIVATE(i, factor) SHARED(a, k, n, perm)
        do i = k+1,n
          factor = A(perm(i), k) / A(perm(k), k)
          A(perm(i),  k+1:n+1) = A(perm(i), k+1:n+1) - factor * A(perm(k), k+1:n+1)
        end do
         !$OMP END PARALLEL DO
      end do
      if (abs(A(perm(n), n)) < eps) then 
        print *, "WARNING: the leading_entry in step", n, " is too small: ", A(perm(n), n)
      end if 

      allocate(x(n))
      x(n) = A(perm(n), n+1)/A(perm(n),n)
      do k = n-1, 1, -1
        x(k) = (A(perm(k),n+1) - dot_product(A(perm(k), k+1:n), x(k+1:n)))/A(perm(k), k)
      end do
   end function gaussian_choice
   function residual(A, B, X) result(r)
      real(dp), intent(in) :: A(:, :), B(:), X(:)
      real(dp) :: r ! модуль вектора невязки

      r = sqrt(sum((matmul(A, X) - B)**2))
   end function residual
end module Gaussian_elimination
