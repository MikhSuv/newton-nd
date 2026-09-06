! Module defining abstract interfaces for vector-valued functions.
! Used to declare functions F: R^n -> R^n for use in numerical methods.
module function_interfaces
   use precision_mod
   implicit none

   abstract interface
      ! Interface for a vector function F mapping R^n -> R^n.
      ! Input:  x(:)   — input vector of length n
      ! Output: y(:)   — result vector of length n
      function vector_func(x) result(y)
         import :: dp
         real(dp), intent(in) :: x(:)
         real(dp) :: y(size(x))
      end function vector_func
   end interface

   abstract interface
      ! Interface for computing the Jacobian matrix of a function F: R^n -> R^n.
      ! Input:  x(:)— point in R^n at which to evaluate the Jacobian
      ! Output: j(n,n) - n x n matrix of partial derivatives
      function jacobian_func(x) result(j)
         import :: dp
         real(dp), intent(in) :: x(:)
         real(dp) :: j(size(x), size(x))
      end function jacobian_func
   end interface

end module function_interfaces
