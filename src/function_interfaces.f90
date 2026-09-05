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

end module function_interfaces
