module derivative_matrix

     use precision_mod
     use function_interfaces
     implicit none

     contains

       function calc_derivatives(X, F) result(H)
         procedure(vector_func) :: F
         real(dp), intent(in) :: X(:) ! точка, в которой считается матрица
         real(dp) :: H(size(X), size(X))
         real(dp) :: X1(size(X))
         real(dp) :: delta_x = sqrt(epsilon(1.0_dp))
         integer :: i, j, n

         n = size(X)

         do i = 1, n
          X1 = X
          X1(i) = X1(i) + delta_x
          H(:, i) = F(X1) - F(X)
         end do

         H = H / delta_x
       end function calc_derivatives

   end module derivative_matrix
