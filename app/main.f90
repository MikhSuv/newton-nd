program main
   use precision_mod
   use newton_nd
   use test_functions

   implicit none
   real(dp) :: x(50) = -12.0_dp

   call newton(G, x, 1000000)

end program main
