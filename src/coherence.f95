subroutine cohMsc(coh, yk1, yk2, nfreq1, nfreq2, k, nOffsets)
  implicit none

  integer :: nfreq1, nfreq2, k, nOffsets, i
  real*8 :: s1(nfreq1), s2(nfreq2), coh(nOffsets, nfreq1)
  complex*16 :: yk1(nfreq1, k), yk2(nfreq2, k), cs(nfreq1)

  s1 = sum(abs(yk1)**2, 2)
  s2 = sum(abs(yk2)**2, 2)

  do i = 1, nOffsets
    cs = sum( yk1 * conjg(yk2(i:(i+nfreq1-1), :)), 2 )
    coh(i, :) = realpart(cs)**2 + imagpart(cs)**2
    coh(i, :) = coh(i, :) / ( s1 * s2(i:(i+nfreq1-1)) )
  end do

end subroutine cohMsc
