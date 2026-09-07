! Module for computing the numerical matrix of partial derivatives
! of a vector-valued function F: R^n -> R^n using forward finite differences.
module derivative_matrix

   use precision_mod
   use function_interfaces
   implicit none
   PRIVATE

   PUBLIC :: calc_derivatives

contains

   ! Compute the n x n Jacobian matrix of function F at point X.
   ! Uses forward finite differences with step size delta_x = sqrt(machine epsilon).
   ! Arguments:
   !   X(:) — point in R^n at which to evaluate the Jacobian
   !   F    — vector-valued function conforming to the vector_func interface
   ! Returns:
   !   H(n,n) — numerical approximation of the Jacobian matrix
   function calc_derivatives(X, F) result(H)
      procedure(vector_func) :: F
      real(dp), intent(in) :: X(:)
      real(dp) :: H(size(X), size(X))
      real(dp) :: X1(size(X))
      real(dp) :: delta_x = sqrt(epsilon(1.0_dp))
      integer :: i, j, n

      n = size(X)

      ! Forward difference: dF_j/dX_i ≈ (F(X + delta_x*e_i) - F(X)) / delta_x
      do i = 1, n
         X1 = X
         X1(i) = X1(i) + delta_x
         H(:, i) = F(X1) - F(X)
      end do

      H = H/delta_x
   end function calc_derivatives

end module derivative_matrix
