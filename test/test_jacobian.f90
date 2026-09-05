! Test program for verifying the numerical Jacobian matrix computation.
! Compares the forward-difference Jacobian against the analytic Jacobian
! for the test function F(x1,x2) = [x1^2 + x2^2 - 2, x1^2 - x2].
!
! Analytic Jacobian:
!   J = [[2*x1, 2*x2],
!        [2*x1,   -1]]
program test_jacobian

  use precision_mod
  use test_functions
  use derivative_matrix

  implicit none

  real(dp), parameter :: tol = 1.0e-6_dp
  integer :: failures, total

  failures = 0
  total = 0

  print *, "=== Jacobian matrix tests ==="
  print *

  call test_case([1.0_dp, 1.0_dp],  &
       reshape([2.0_dp, 2.0_dp, 2.0_dp, -1.0_dp], [2,2]), "point (1, 1)")

  call test_case([-1.0_dp, 1.0_dp], &
       reshape([-2.0_dp, -2.0_dp, 2.0_dp, -1.0_dp], [2,2]), "point (-1, 1)")

  call test_case([0.0_dp, 0.0_dp],  &
       reshape([0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp], [2,2]), "point (0, 0)")

  call test_case([-1.0_dp, -1.0_dp], &
       reshape([-2.0_dp, -2.0_dp, -2.0_dp, -1.0_dp], [2,2]), "point (-1, -1)")

  call test_case([0.5_dp, 2.0_dp],  &
       reshape([1.0_dp, 1.0_dp, 4.0_dp, -1.0_dp], [2,2]), "point (0.5, 2)")

  print *
  write(*, '(A, I0, A, I0, A)') "Results: ", total - failures, " / ", total, " passed"

  if (failures > 0) then
    print *, "FAILED"
    stop 1
  else
    print *, "ALL TESTS PASSED"
  end if

contains

  ! Run a single test case: compute numerical Jacobian at point X,
  ! compare each entry against the expected analytic Jacobian J_exact.
  subroutine test_case(X, J_exact, label)
    real(dp), intent(in) :: X(:)
    real(dp), intent(in) :: J_exact(:, :)
    character(len=*), intent(in) :: label

    real(dp) :: H(size(X), size(X))
    integer :: i, j, n

    n = size(X)
    H = calc_derivatives(X, F)

    print *, "--- Test: " // trim(label) // " ---"
    do i = 1, n
      do j = 1, n
        total = total + 1
        if (abs(H(i,j) - J_exact(i,j)) > tol) then
          write(*, '(A,I0,A,I0,A,ES12.5,A,ES12.5,A,ES12.5)') &
            "  FAIL: H(", i, ",", j, ") = ", H(i,j), &
            "  expected ", J_exact(i,j), &
            "  diff ", abs(H(i,j) - J_exact(i,j))
          failures = failures + 1
        else
          write(*, '(A,I0,A,I0,A,ES12.5)') &
            "  PASS: H(", i, ",", j, ") = ", H(i,j)
        end if
      end do
    end do
    print *
  end subroutine test_case

end program test_jacobian
