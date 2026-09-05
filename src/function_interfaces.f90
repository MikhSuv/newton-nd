module function_interfaces
  use precision_mod
  implicit none

  abstract interface
    function vector_func(x) result(y)
      import :: dp
      real(dp), intent(in) :: x(:)
      real(dp) :: y(size(x))
    end function vector_func
  end interface

end module function_interfaces
