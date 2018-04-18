subroutine testMatrix(x)
  implicit none

  integer :: x(4, 3), i

  print *, ""

  do i = 1, 4
    print *, x(i, :)
  end do
end subroutine testMatrix
