! Module containing function F that computes a vector-valued transformation
module test_functions

  use precision_mod
  
  implicit none
  
contains
    
    ! Function F maps R^2 -> R^2
    ! Input X = (x1, x2) is an input vector of length 2
    ! Output Y returns the computed transformation where:
    !   y(1) = x1**2 + x2**2 - 2.0_dp  
    !   y(2) = x1**2 - x2
    function F(X) result(Y)

      real(dp), intent(in) :: X(:)
      real(dp)             :: Y(size(x))

      Y(1) = X(1)**2 + X(2)**2 - 2.0_dp
      Y(2) = X(1)**2 - X(2)
    end function F


end module test_functions
