! Test program for verifying the Newton solver (newton_solve).
! Checks that the solver converges to a root of F or G from various
! initial guesses by verifying ||F(X)|| < tol at the solution.
!
! Known roots:
!  F(x1,x2) = [x1^2 + x2^2 - 2, x1^2 - x2]  has roots (1,1) and (-1,1).
!  G(X) (50D cyclic system) has trivial roots (1, 1, ..., 1) and (-2, -2, ..., -2).
program test_newton

   use precision_mod
   use function_interfaces
   use test_functions
   use newton_nd, only: newton_solve

   implicit none

   real(dp), parameter :: tol = 1.0e-6_dp
   integer, parameter :: max_iter = 1000000
   integer :: failures, total
   real(dp), allocatable :: x(:)

   failures = 0
   total = 0

   print *, "=== Newton solver tests ==="
   print *

   ! Tests for F: R^2 -> R^2
   call test_case(F, [100.0_dp, 100.0_dp], "F from (100, 100)")
   call test_case(F, [-2.0_dp, 2.0_dp], "F from (-2, 2)")

   ! Tests for G: R^50 -> R^50
   allocate (x(g_dim))
   x = 100.0_dp
   call test_case(G, x, "G from (100, ..., 100)")
   x = -100.0_dp
   call test_case(G, x, "G from (-100, ..., -100)")
   deallocate (x)

   print *
   write (*, '(A, I0, A, I0, A)') "Results: ", total - failures, " / ", total, " passed"

   if (failures > 0) then
      print *, "FAILED"
      stop 1
   else
      print *, "ALL TESTS PASSED"
   end if

contains

   ! Run a single test case: solve F(X)=0 via newton_solve starting from X0,
   ! then check that ||F(X)|| < tol.
   subroutine test_case(F, X0, label)
      procedure(vector_func) :: F
      real(dp), intent(in) :: X0(:)
      character(len=*), intent(in) :: label

      real(dp) :: X(size(X0))
      real(dp) :: residual

      X = newton_solve(F, X0, max_iter)
      residual = norm2(F(X))

      total = total + 1
      if (residual < tol) then
         write (*, '(A, A, A, ES12.5)') "  PASS: ", trim(label), "  ||F(X)|| = ", residual
      else
         write (*, '(A, A, A, ES12.5)') "  FAIL: ", trim(label), "  ||F(X)|| = ", residual
         failures = failures + 1
      end if
   end subroutine test_case

end program test_newton
