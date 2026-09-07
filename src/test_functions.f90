! Test functions for the Newton solver and Jacobian tests.
! Defines two vector-valued functions with known analytic Jacobians.
module test_functions

   use precision_mod

   implicit none

   ! Dimension of the test system G: R^g_dim -> R^g_dim
   integer, parameter :: g_dim = 50

contains

   ! Function F maps R^2 -> R^2.
   ! F(X) = [x1^2 + x2^2 - 2, x1^2 - x2]^T
   ! Roots: (1, 1) and (-1, 1).
   function F(X) result(Y)

      real(dp), intent(in) :: X(:)
      real(dp)             :: Y(size(x))

      Y(1) = X(1)**2 + X(2)**2 - 2.0_dp
      Y(2) = X(1)**2 - X(2)
   end function F

   ! Function G maps R^g_dim -> R^g_dim (cyclic coupled system).
   ! G(i) = X(i)^2 + X(i+1) - 2  for i = 1..g_dim-1
   ! G(g_dim) = X(g_dim)^2 + X(1) - 2  (cyclic wrap-around)
   ! Trivial roots: X = (1, 1, ..., 1) and (-2, -2, ..., -2).
   function G(X) result(Y)

      real(dp), intent(in) :: X(:)
      real(dp)             :: Y(size(x))
      integer :: i

      do i = 1, g_dim - 1
         Y(i) = X(i)**2 + X(i + 1) - 2.0_dp
      end do
      Y(g_dim) = X(g_dim)**2 + X(1) - 2.0_dp

   end function G

end module test_functions
