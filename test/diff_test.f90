! Program to test computation of the Jacobian matrix.
! Prints the numerical Jacobian at several test points for visual inspection.
program test_diff

   use precision_mod
   use test_functions
   use derivative_matrix

   implicit none
   real(dp) :: H(2, 2)
   real(dp) :: m1(2) = [1, 1]
   real(dp) :: m2(2) = [-1, 1]
   real(dp) :: m3(2) = [0, 0]
   real(dp) :: m4(2) = [-1, -1]
   integer :: i

   ! Evaluate and print Jacobian at (1, 1)
   H = calc_derivatives(m1, F)
   do i = 1, size(m1)
      write (*, '(*(F15.8))') H(i, :)
   end do
   ! Evaluate and print Jacobian at (-1, 1)
   H = calc_derivatives(m2, F)
   do i = 1, size(m1)
      write (*, '(*(F15.8))') H(i, :)
   end do
   ! Evaluate and print Jacobian at (0, 0)
   H = calc_derivatives(m3, F)
   do i = 1, size(m1)
      write (*, '(*(F15.8))') H(i, :)
   end do
   ! Evaluate and print Jacobian at (-1, -1)
   H = calc_derivatives(m4, F)
   do i = 1, size(m1)
      write (*, '(*(F15.8))') H(i, :)
   end do
end program test_diff
