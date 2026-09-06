! Main application entry point for the newton-nd project.
program main
  use precision_mod
  use newton_nd
  use test_functions

  real(dp) :: x1(50) = 100.0_dp
  real(dp) :: x2(50) = -100.0_dp

  real(dp) :: x4(2) = 100.0_dp
  real(dp) :: x5(2) = [-2.0_dp, 2.0_dp]
  integer :: max_iter = 1000000

  call newton(G, x1, max_iter)
  call newton(G, x2, max_iter)
  call newton(F, x4, max_iter)
  call newton(F, x5, max_iter)

  
end program main
