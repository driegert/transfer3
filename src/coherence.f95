subroutine cohMsc(coh, yk1, yk2, nfreq1, nfreq2, k, nOffsets)
  implicit none

  integer :: nfreq1, nfreq2, k, nOffsets, i
  real*8 :: s1(nfreq1), s2(nfreq2), coh(nOffsets, nfreq1)
  complex*16 :: yk1(nfreq1, k), yk2(nfreq2, k), cs(nfreq1)

  s1 = sum(abs(yk1)**2, 2)
  s2 = sum(abs(yk2)**2, 2)

  do i = 1, nOffsets
    cs(:) = sum( yk1(:, :) * conjg(yk2(i:(i+nfreq1-1), :)), 2 )
    coh(i, :) = realpart(cs(:))**2 + imagpart(cs(:))**2
    coh(i, :) = coh(i, :) / ( s1(:) * s2(i:(i+nfreq1-1)) )
  end do

end subroutine cohMsc

subroutine msc2indicator(msc, nrow, ncol, ind, level, nOff)
! determines if the mcs's are local maxes (by column, i.e., central freq)
! and if the value exceeds the level provided (calculated based on the
! normal distribution usually, but will depend what values are in the
! msc matrix)
  real*8 :: msc(nrow, ncol), cohLmax(nrow, ncol), level
  integer :: i, j, ind(nrow, ncol), nrow, ncol, nOff, cmax(ncol)

  ind(:, :) = 0
  cohLmax(:, :) = -1

! only keep the local maxes above the cutoff.
  do i = 2, nrow-1
    do j = 1, ncol
      if (msc(i,j) >= level .and. msc(i, j) > msc(i-1, j) .and. &
        msc(i, j) > msc(i+1, j)) then
        cohLmax(i, j) = msc(i, j)
      end if
    end do
  end do

  do i = 1, nOff
    cmax = maxloc(msc, dim = 1, mask = (ind .eq. 0 .and. cohLmax > 0))
    do j = 1, ncol
      if (cmax(j) .eq. 0) cycle
      ind(cmax(j), j) = i
    end do
  end do
end subroutine msc2indicator
