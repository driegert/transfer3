
!!!! This does NOT include blocking.
subroutine tf(H, yk1, yk2, cohInd, nrow1, nrow2, npred, nBlocks &
  , nOffsets, nUniqOffsets, k)
  implicit none

  integer :: cohInd(nOffsets, nrow1, npred), nrow1, nrow2, npred &
    , nUniqOffsets, k, nDesCol(npred+1), f, p, oFreqIdx(nOffsets) &
    , nOffsets, nBlocks, j, b, nHcoef, offsetIndAll(nOffsets, npred) &
    , HnColByPred(npred + 1)
  integer, allocatable :: Hind(:)
  real*8, allocatable :: stdErr(:), svd_ev(:)
  complex*16 :: H(nrow1, nUniqOffsets), yk1(nrow1, k, nBlocks) &
    , yk2(nrow2, k, nBlocks, npred), Y(nBlocks*k)
  complex*16, allocatable :: design(:, :), Hwrk(:)

  interface
    subroutine zSvdRegression(Y, X, m, n, beta, stdErr, svd_ev)
      implicit none
      integer :: i, m, n, lda, ldu, ldvt, lwork, info, lrwork
      real*8 :: svd_ev(n), stdErr(n), lworkopt
      real*8, allocatable :: rwork(:), s(:)
      complex*16 :: Y(m), X(m, n), Xcopy(m,n), beta(n)
      complex*16, allocatable :: svd_work(:), u(:, :), vt(:, :)
      character(1) :: jobu, jobvt
    end subroutine zSvdRegression

    subroutine designMatrix(Y, design, yk1, yk2, f, nrow1, nrow2 &
      , npred, k, nBlocks, nDesCol, nHcoef, cohInd, oFreqIdx, nOffsets)
      implicit none
      integer :: f, nrow1, nrow2, npred, k, nBlocks &
        , cohInd(nOffsets, nrow1, npred) &
        , oFreqIdx(nOffsets), p, b, nHcoef, nOffsets, nDesCol(npred+1)
      complex*16 :: Y(nBlocks*k), design(nBlocks*k, nHcoef) &
        , yk1(nrow1, k, nBlocks), yk2(nrow2, k, nBlocks, npred)
      end subroutine designMatrix

      subroutine determineHcol(cohInd, offsetIndAll, HnColByPred, npred &
          , oFreqIdx, nOffsets)
        implicit none
        integer :: cohInd(:, :, :), offsetIndAll(nOffsets, npred), npred &
          , HnColByPred(npred + 1), oFreqIdx(nOffsets), p, ones(nOffsets) &
          , nOffsets, i, cur
      end subroutine determineHcol

      subroutine findHcol(cohIndCol, offsetIndAll, Hind &
        , nOffsets, npred, nHcoef)
        implicit none
        integer :: cohIndCol(nOffsets, npred), offsetIndAll(nOffsets, npred) &
          , Hind(nHcoef), nOffsets, npred, nHcoef
      end subroutine findHcol
  end interface

  nDesCol(1) = 1
  ! used to get rows of yk2 to use - row indices in cohInd basically...
  oFreqIdx = (/ (f, f = 1, nOffsets) /)

  call determineHcol(cohInd, offsetIndAll, HnColByPred, npred &
    , oFreqIdx, nOffsets) ! set this up

  do f = 1, nrow1 ! frequency / yk1 row index / cohInd col index
    do p = 1, npred ! predictor
      nDesCol(p+1) = nDesCol(p) + sum(cohInd(:, f, p))
      ! print *, "sum(cohInd(:, f, p)) = ", sum(cohInd(:, f, p))
    end do

    nHcoef = nDesCol(npred+1)-1 ! number of columns in design mat at this f

    allocate(design(nBlocks*k, nHcoef), Hwrk(nHcoef), Hind(nHcoef))
    allocate(stdErr(nHcoef), svd_ev(nHcoef))

    call designMatrix(Y = Y, design = design, yk1 = yk1, yk2 = yk2 &
      , f = f, nrow1 = nrow1, nrow2 = nrow2, npred = npred, k = k &
      , nBlocks = nBlocks, nDesCol = nDesCol, nHcoef = nHcoef &
      , cohInd = cohInd, oFreqIdx = oFreqIdx, nOffsets = nOffsets)

    call zSvdRegression(Y = Y, X = design, m = nBlocks*k, n = nHcoef &
      , beta = Hwrk, stdErr = stdErr, svd_ev = svd_ev)

    call findHcol(cohInd(:, f, :), offsetIndAll, Hind, nOffsets &
      , npred, nHcoef)

    H(f, Hind) = Hwrk

    deallocate(design, Hwrk, stdErr, svd_ev, Hind)
  end do
end subroutine tf


subroutine designMatrix(Y, design, yk1, yk2, f, nrow1, nrow2 &
  , npred, k, nBlocks, nDesCol, nHcoef, cohInd, oFreqIdx, nOffsets)
  implicit none

  integer :: f, nrow1, nrow2, npred, k, nBlocks &
    , cohInd(nOffsets, nrow1, npred) &
    , oFreqIdx(nOffsets), p, b, nHcoef, nOffsets, nDesCol(npred+1)
  complex*16 :: Y(nBlocks*k), design(nBlocks*k, nHcoef) &
    , yk1(nrow1, k, nBlocks), yk2(nrow2, k, nBlocks, npred)

  do p = 1, npred
    do b = 0, nBlocks - 1
      design( (b*k+1):((b+1)*k), nDesCol(p):(nDesCol(p+1)-1)) = &
        transpose(yk2( (f-1)+pack(oFreqIdx, cohInd(:, f, p) > 0), :, b+1, p) )

      if (p .eq. 1) then
        Y( (b*k+1):((b+1)*k) ) = yk1(f, :, b+1)
      end if
    end do
  end do
end subroutine designMatrix

subroutine determineHcol(cohInd, offsetIndAll, HnColByPred, npred &
    , oFreqIdx, nOffsets)
  implicit none

  integer :: cohInd(:, :, :), offsetIndAll(nOffsets, npred), npred &
    , HnColByPred(npred + 1), oFreqIdx(nOffsets), p, ones(nOffsets) &
    , nOffsets, i, cur

  ones(:) = 1
  offsetIndAll(:, :) = 0
  HnColByPred(1) = 1
  cur = 0

  print *, "npred = ", npred
  do p = 1, npred
    offsetIndAll(pack(oFreqIdx(:), sum(cohInd(:, :, p), 2) > 0), p) = 1
    HnColByPred(p+1) = HnColByPred(p) + sum(offsetIndAll(:, p))

    do i = 1, nOffsets
      if (offsetIndAll(i,p) > 0) then
        cur = cur + offsetIndAll(i, p)
        offsetIndAll(i,p) = cur
      end if
    end do
  end do
end subroutine determineHcol

subroutine findHcol(cohIndCol, offsetIndAll, Hind &
  , nOffsets, npred, nHcoef)
  implicit none

  integer :: cohIndCol(nOffsets, npred), offsetIndAll(nOffsets, npred) &
    , Hind(nHcoef), nOffsets, npred, nHcoef

  Hind = pack(offsetIndAll, cohIndCol > 0)
end subroutine findHcol

subroutine zSvdRegression(Y, X, m, n, beta, stdErr, svd_ev)
! uses the SVD of the matrix X to calculate the regression coefficients for
! the model: Y = XB (where B is a vector of beta coefficients)
!
! as a side note: this mxn notation is bad... as a result of lapack docs
!
! [in] Y      - complex*16(m)  - response vector
! [in] X      - complex*16(m,n)  - design matrix
! [in] m      - integer - number of rows of both Y and X
! [in] n      - integer - number of columns of X
! [out] beta  - complex*16(n) - the regression coefficients
! [out] stdErr  - real*8(n) - the standard error on the beta estimates
! [out] ev    - real*8(n) - the square of the singular values
  implicit none

  integer :: i, m, n, lda, ldu, ldvt, lwork, info, lrwork
  real*8 :: svd_ev(n), stdErr(n), lworkopt
  real*8, allocatable :: rwork(:), s(:)
  complex*16 :: Y(m), X(m, n), Xcopy(m,n), beta(n)
  complex*16, allocatable :: svd_work(:), u(:, :), vt(:, :)
  character(1) :: jobu, jobvt

  ! 'S' says that we want the left and right singular vectors
  jobu = 'S'
  jobvt = 'S'
  lda = m
  ldu = m
  ldvt = n
  Xcopy = X
  ! set values according to:
  ! http://www.netlib.org/lapack/explore-html/index.html
  ! search for zgesvd
  lrwork = 5*min(m, n)
  allocate(rwork(lrwork))
  allocate(s(min(m,n)))
  allocate(u(ldu, min(m,n)))
  allocate(vt(ldvt, n))

  ! obtain optimal size for lwork
  lwork = -1
  call zgesvd(jobu, jobvt, m, n, Xcopy, lda, s, u, ldu, vt, ldvt &
    , lworkopt, lwork, rwork, info)

  ! allocate the work array
  lwork = nint(lworkopt)
  allocate(svd_work(lwork))

  ! perform the svd
  call zgesvd(jobu, jobvt, m, n, Xcopy, lda, s, u, ldu, vt, ldvt &
    , svd_work, lwork, rwork, info)

  ! calculate them betas '*' is matrix mult here
  ! beta = V * s^(-1) * [ t(u) Y ] (mandel - eqn's (17) and (12))
  beta = matmul( transpose(conjg(vt)), matmul( transpose(conjg(u)), Y ) / s )

  ! calculate the standard error estiamtes on the betas
  ! mandel eqn (20) NOT YET IMPLEMENTED
  do i = 1, n
    stdErr(i) = -1.0D0 !sum(vt(:, i)x)
  end do

  do i = 1, n
    svd_ev(i) =  -1.0D0 ! NOT YET IMPLEMENTED
  end do
end subroutine zSvdRegression
