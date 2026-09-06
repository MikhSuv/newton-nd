! Module containing function F that computes a vector-valued transformation
module test_functions

   use precision_mod

   implicit none

   ! Dimension of the test system G: R^g_dim -> R^g_dim
   integer, parameter :: g_dim = 5

contains

   ! Function F maps R^2 -> R^2
   ! Input X = (x1, x2) is an input vector of length 2
   ! Output Y returns the computed transformation
   function F(X) result(Y)

      real(dp), intent(in) :: X(:)
      real(dp)             :: Y(size(x))

      Y(1) = X(1)**2 + X(2)**2 - 2.0_dp
      Y(2) = X(1)**2 - X(2)
   end function F

   ! Function G maps R^g_dim -> R^g_dim
   ! Input X is an input vector of length g_dim
   ! Output Y returns the computed transformation
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
