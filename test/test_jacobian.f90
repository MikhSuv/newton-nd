! Test program for verifying the numerical Jacobian matrix computation.
! Compares the forward-difference Jacobian against the analytic Jacobian
! for the test function F(x1,x2) = [x1^2 + x2^2 - 2, x1^2 - x2].
!
! Analytic Jacobian:
!   J = [[2*x1, 2*x2],
!        [2*x1,   -1]]
program test_jacobian

   use precision_mod
   use function_interfaces
   use test_functions
   use derivative_matrix

   implicit none

   real(dp), parameter :: tol = 1.0e-6_dp
   real(dp), allocatable :: x(:)
   integer :: i
   integer :: failures, total

   failures = 0
   total = 0

   print *, "=== Jacobian matrix tests ==="
   print *

   call test_case([1.0_dp, 1.0_dp], F, F_jacobian, "point (1, 1)")
   call test_case([-1.0_dp, 1.0_dp], F, F_jacobian, "point (-1, 1)")
   call test_case([0.0_dp, 0.0_dp], F, F_jacobian, "point (0, 0)")
   call test_case([-1.0_dp, -1.0_dp], F, F_jacobian, "point (-1, -1)")
   call test_case([0.5_dp, 2.0_dp], F, F_jacobian, "point (0.5, 2)")

   allocate (x(g_dim), source=1.0_dp)
   call test_case(x, G, G_jacobian, "point (x_i = 1), i = 1,...,5")
   deallocate (x)
   allocate (x(g_dim), source=-2.0_dp)
   call test_case(x, G, G_jacobian, "point (x_i = -2), i = 1,...,5")
   deallocate (x)
   allocate (x(g_dim), source=0.0_dp)
   call test_case(x, G, G_jacobian, "point (x_i = 0), i = 1,...,5")

   do i = 1, g_dim
      x(i) = real(i, dp)
   end do
   call test_case(x, G, G_jacobian, "point (x_i = i), i = 1,...,5")
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

   ! Run a single test case: compute numerical Jacobian at point X via
   ! calc_derivatives, compare each entry against the analytic Jacobian from D.
   subroutine test_case(X, Func, D, label)
      real(dp), intent(in) :: X(:) ! test point in R^n
      procedure(vector_func) :: Func ! vector function F: R^n -> R^n
      procedure(jacobian_func) :: D ! analytic Jacobian of F
      character(len=*), intent(in) :: label ! test case name for output

      real(dp) :: J_exact(size(x), size(x))
      real(dp) :: H(size(X), size(X))
      integer :: i, j, n

      n = size(X)
      H = calc_derivatives(X, Func)
      J_exact = D(X)

      print *, "--- Test: "//trim(label)//" ---"
      do i = 1, n
         do j = 1, n
            total = total + 1
            if (abs(H(i, j) - J_exact(i, j)) > tol) then
               write (*, '(A,I0,A,I0,A,ES12.5,A,ES12.5,A,ES12.5)') &
                  "  FAIL: H(", i, ",", j, ") = ", H(i, j), &
                  "  expected ", J_exact(i, j), &
                  "  diff ", abs(H(i, j) - J_exact(i, j))
               failures = failures + 1
            else
               write (*, '(A,I0,A,I0,A,ES12.5)') &
                  "  PASS: H(", i, ",", j, ") = ", H(i, j)
            end if
         end do
      end do
      print *
   end subroutine test_case

   function F_jacobian(X) result(J)
      real(dp), intent(in) :: X(:)
      real(dp) :: J(size(X), size(X))

      J(:, 1) = 2.0_dp*x(1)
      J(1, 2) = 2.0_dp*x(2)
      J(2, 2) = -1.0_dp

   end function F_jacobian

   function G_jacobian(X) result(J)
      real(dp), intent(in) :: X(:)
      real(dp) :: J(size(X), size(X))

      integer :: n, i
      n = size(X)

      do i = 1, n - 1
         J(i, i) = 2.0_dp*x(i)
         J(i, i + 1) = 1.0_dp
      end do

      j(n, n) = 2.0_dp*x(n)
      j(n, 1) = 1.0_dp

   end function G_jacobian

end program test_jacobian
