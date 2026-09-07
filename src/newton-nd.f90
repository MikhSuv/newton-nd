! Module implementing the multi-dimensional Newton method
! for solving a system of nonlinear equations F(X) = 0.
!
! Algorithm at each iteration:
!   1. Compute the Jacobian matrix H = dF/dX at current guess X_
!   2. Solve the linear system H*X = H*X_ - F(X_) for the new X
!   3. Check convergence: ||X - X_|| < tol
module newton_nd
   use precision_mod
   use derivative_matrix
   use function_interfaces
   use Gaussian_elimination, only: solve_linear_system
   ! use iterative_linear_system, only: iterative_solve
   implicit none
   private

   real(dp), parameter :: tol = 1.0e2_dp*epsilon(1.0_dp)

   public :: newton_solve, newton

contains

   ! Solve F(X) = 0 using Newton iteration.
   !   F      — vector function F: R^n -> R^n
   !   X0     — initial guess
   !   max_iter — maximum number of iterations
   ! Returns the approximate root X such that ||F(X)|| ~ 0.
   function newton_solve(F, X0, max_iter) result(X)
      procedure(vector_func) :: F
      real(dp), intent(in) :: X0(:)
      integer, intent(in) :: max_iter
      real(dp) :: X(size(X0)), X_(size(X0))

      real(dp) :: H(size(X0), size(X0))
      integer :: i

      X_ = X0

      do i = 1, max_iter
         H = calc_derivatives(X_, F)
         ! Right-hand side: H*X_ - F(X_)
         X = solve_linear_system(H, matmul(H, X_) - F(X_), "choice")
         ! X = iterative_solve(H, matmul(H, X_) - F(X_), "relaxation")
         if (maxval(abs(X - X_)) < tol) then
            return
         else
            X_ = X
         end if
      end do
      print *, "Maximum number of iterations reached"
   end function newton_solve

   ! Write solution vector X to a file
   subroutine write_result(filename, X)
      character(len=*), intent(in) :: filename
      real(dp), intent(in) ::  X(:) ! Result column vector

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

   ! High-level driver: solve F(X)=0 via newton_solve, write result to file,
   ! and print the residual norm ||F(X)||.
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
