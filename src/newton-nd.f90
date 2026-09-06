module newton_nd
   use precision_mod
   use derivative_matrix
   use function_interfaces
   use Gaussian_elimination, only : solve_linear_system
   implicit none
   private

   real(dp), parameter :: tol = 1.0e2_dp * epsilon(1.0_dp)

   public :: newton_solve, newton

contains

  function newton_solve(F, X0, max_iter) result(X)
    procedure(vector_func) :: F
    real(dp), intent(in) :: X0(:)
    integer, intent(in) :: max_iter
    real(dp) :: X(size(X0)), X_(size(X0))

    real(dp) :: H(size(X0), size(X0)), B(size(X0))
    integer :: i

    X_ = X0

    do i = 1, max_iter
      H = calc_derivatives(X_, F)
      B = matmul(H, X_) - F(X_)
       X = solve_linear_system(H, B, "choice")
      if (maxval(abs(X-X_)) < tol) then
        return 
      else 
        X_ = X
      end if
    end do
    print *, "Достигнуто максимальное число итераций"
  end function newton_solve

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

      do i = 1, n
         write (ounit, '(*(e23.15, 1x))') X(i)
      end do

      close (ounit)

   end subroutine write_result

  subroutine newton(F, X0, max_iter, filename) 
    procedure(vector_func) :: F
    real(dp), intent(in) :: X0(:)
    integer, intent(in) :: max_iter
    character(len=*), intent(in), optional :: filename

    real(dp) :: X(size(X0))

    X = newton_solve(F, X0, max_iter)

    if (present(filename)) then 
      call write_result(filename, X)
    else 
      call write_result("result.dat", X)
    end if

    print *, "||F(X)|| = ", NORM2(F(X))

  end subroutine newton

end module newton_nd
