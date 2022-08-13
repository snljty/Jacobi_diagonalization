subroutine reference_diagonalization(n, a, w, v, uplo, status)
    implicit none
    integer, intent(in) :: n
    double precision, dimension(n, n), intent(in) :: a
    double precision, dimension(n), intent(out) :: w
    double precision, dimension(n, n), intent(out) :: v
    character(len=1), intent(in) :: uplo
    integer, intent(out) :: status

    integer :: lwork
    double precision, dimension(:), allocatable :: work

    call dlacpy("A", n, n, a, n, v, n)
    lwork = -1
    allocate(work(- lwork))
    call dsyev("V", uplo, n, v, n, w, work, lwork, status)
    lwork = nint(work(1))
    deallocate(work)
    allocate(work(lwork))
    call dsyev("V", uplo, n, v, n, w, work, lwork, status)
    deallocate(work)

    return
end subroutine reference_diagonalization

