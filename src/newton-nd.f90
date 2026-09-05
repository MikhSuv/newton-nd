module newton_nd
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, newton-nd!"
  end subroutine say_hello
end module newton_nd
