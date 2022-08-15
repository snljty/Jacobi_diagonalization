! Jacobi diagonalization algorithm to get the eigenvalues and eigenvectors
! currently there is no method for handling degenerate situation

program main
    implicit none

    interface
        subroutine Jacobi_diagonalization(n, a, w, v, vecp, vecq, uplo, status, destroy)
            implicit none
            integer, intent(in) :: n
            double precision, dimension(n, n), intent(inout), target :: a
            double precision, dimension(n), intent(out) :: w
            double precision, dimension(n, n), intent(out) :: v
            double precision, dimension(n), intent(inout) :: vecp, vecq
            character(len=1), intent(in) :: uplo
            integer, intent(out) :: status
            logical, optional :: destroy
        end subroutine Jacobi_diagonalization
    end interface

    integer :: argc, iarg, argl, argl_max
    character(len=:), allocatable :: argv0
    character(len=:), dimension(:), allocatable :: argv

    integer :: n
    integer, parameter :: n_def = 5
    double precision, dimension(:, :), pointer :: a
    double precision, dimension(:, :), allocatable :: v
    double precision, dimension(:), allocatable :: w
    double precision, dimension(:), allocatable :: vecp, vecq
    integer :: status
    integer :: i, j

    double precision, dimension(:, :), allocatable :: v_cmp
    double precision, dimension(:), allocatable :: w_cmp
    logical, parameter :: is_compare = .true.

    ! get command arguments
    argc = command_argument_count()
    call get_command_argument(0, length = argl)
    allocate(character(len=argl) :: argv0)
    if (argc > 0) then
        argl_max = 0
        do iarg = 1, argc
            call get_command_argument(iarg, length = argl)
            if (argl > argl_max) then
                argl_max = argl
            end if
        end do
        argl = argl_max
        allocate(character(len=argl_max) :: argv(argc))
        do iarg = 1, argc
            call get_command_argument(iarg, argv(iarg))
        end do
    end if

    ! pharse command arguments
    if (argc == 1) then
        read(argv(1), *) n
    else
        n = n_def
    end if

    ! release memory of command arguments
    deallocate(argv0)
    if (argc > 0) then
        deallocate(argv)
    end if

    ! allocate memory for arrays
    allocate(a(n, n), v(n, n), w(n), vecp(n), vecq(n))

    ! initialize a
    a = 0.d0
    do i = 1, n
        j = 1 + mod(i, n)
        a(i, j) = 1.d0
        a(j, i) = 1.d0
    end do
    a(1, n) = 0.d0
    a(n, 1) = 0.d0

    ! or
    ! call random_number(a)
    ! a = matmul(transpose(a), a)

    write(*, '(a)') "Matrix:"
    call print_square_matrix(n, a)
    write(*, '()')

    if (is_compare) then
        allocate(v_cmp(n, n), w_cmp(n))
        call reference_diagonalization(n, a, w_cmp, v_cmp, "L", status)
        if (status /= 0) then
            write(0, '(a, i0, a)') "Error! LAPACK dsyev returns ", status, " ."
        end if
    end if

    call Jacobi_diagonalization(n, a, w, v, vecp, vecq, "L", status)
    if (status /= 0) then
        write(0, '(a)') repeat("#", 64)
        write(0, '(a, a, i2, a, a)') repeat("#", 8), " Error! The algorithm returns ", status, " instead of 0 . ", repeat("#", 8)
        write(0, '(a)') repeat("#", 64)
    end if
    call sort_eigen(n, w, v)

    write(*, '(a)') "Eigenvalues:"
    call print_vector(n, w)
    write(*, '()')
    write(*, '(a)') "Eigenvectors:"
    call print_square_matrix(n, v)
    write(*, '()')

    deallocate(a, v, w, vecp, vecq)
    nullify(a)

    if (is_compare) then
        write(*, '(a)') "Reference: "
        write(*, '()')
        write(*, '(a)') "Eigenvalues:"
        call print_vector(n, w_cmp)
        write(*, '()')
        write(*, '(a)') "Eigenvectors:"
        call print_square_matrix(n, v_cmp)
        write(*, '()')
        deallocate(v_cmp, w_cmp)
    end if

    stop
end program main

subroutine print_vector(n, a)
    implicit none
    integer, intent(in) :: n
    double precision, dimension(n), intent(in) :: a

    integer :: j

    do j = 1, n
        if (j /= 1) then
            write(*, "(2x)", advance = "no")
        end if
        write(*, "(f10.7)", advance = "no") a(j)
    end do
    write(*, "()")

    return
end subroutine print_vector

subroutine print_matrix(m, n, a)
    implicit none
    integer, intent(in) :: m, n
    double precision, dimension(m, n), intent(in) :: a

    integer :: i

    do i = 1, m
        call print_vector(n, a(i, :))
    end do

    return
end subroutine print_matrix

subroutine print_square_matrix(n, a)
    implicit none
    integer, intent(in) :: n
    double precision, dimension(n, n), intent(in) :: a

    call print_matrix(n, n, a)

    return    
end subroutine print_square_matrix

subroutine Jacobi_diagonalization(n, a, w, v, vecp, vecq, uplo, status, destroy)
    implicit none
    integer, intent(in) :: n
    double precision, dimension(n, n), intent(inout), target :: a
    double precision, dimension(n), intent(out) :: w
    double precision, dimension(n, n), intent(out) :: v
    double precision, dimension(n), intent(inout) :: vecp, vecq
    character(len=1), intent(in) :: uplo
    integer, intent(out) :: status
    logical, optional :: destroy

    integer :: i, j, p, q
    double precision, dimension(:, :), allocatable, target :: tmp
    double precision, dimension(:, :), pointer :: mat
    double precision :: abs_pq

    integer :: iloop, nloops_max
    double precision, parameter :: tol = epsilon(1.d0)

    ! check uplo
    if (uplo /= "U" .and. uplo /= "u" .and. uplo /= "L" .and. uplo /= "l") then
        status = -2
        return
    end if

    ! check n and if we can quick return
    status = 0
    if (n < 1) then
        status = -1
        return
    else if (n == 1) then
        v(1, 1) = 1.d0
        w(1) = a(1, 1)
        return
    else
        continue
    end if

    ! allocate temporary matrix if needed, assign pointer to the actual matrix
    if (present(destroy)) then
        if (destroy) then
            mat => a
            goto 100
        end if
    end if
    allocate(tmp(n, n))
    tmp(1:n, 1:n) = a(1:n, 1:n)
    mat => tmp
100 continue

    ! initialize eigenvectors as unit matrix
    v = 0.d0
    forall (i = 1:n)
        v(i, i) = 1.d0
    end forall

    ! maximum number of loops 
    nloops_max = 30 * n ** 2 ! do not set less than 2 * n ** 2
    ! 30 is from OpenCV

    iloop = 0
    ! the upper triangle will be used
    if (uplo == "U" .or. uplo == "u") then
        do while (iloop < nloops_max)
            p = 2
            q = 1
            abs_pq = abs(mat(p, q))
            do j = 2, n
                do i = 1, j - 1
                    if (abs(mat(i, j)) > abs_pq) then
                        p = i
                        q = j
                        abs_pq = abs(mat(p, q))
                    end if
                end do
            end do
            if (abs_pq < tol) then
                exit
            end if
            call Jacobi_rotate(n, p, q, mat, v, vecp, vecq, status)
            iloop = iloop + 1
        end do
    ! the lower triangle will be used
    else if (uplo == "L" .or. uplo == "l") then
        do while (iloop < nloops_max)
            p = 1
            q = 2
            abs_pq = abs(mat(p, q))
            do j = 1, n - 1
                do i = j + 1, n
                    if (abs(mat(i, j)) > abs_pq) then
                        p = i
                        q = j
                        abs_pq = abs(mat(p, q))
                    end if
                end do
            end do
            if (abs_pq < tol) then
                exit
            end if
            ! rotate once
            call Jacobi_rotate(n, p, q, mat, v, vecp, vecq, status)
            iloop = iloop + 1
        end do
    else
        continue ! never happen
    end if

    ! get the eigenvalues from the diagonal of the target matrix
    forall (i = 1:n)
        w(i) = mat(i, i)
    end forall

    ! if not converged, returns 1
    if (iloop == nloops_max) then
        status = 1
    end if

    ! nullify pointer, release memory if allocated in this subroutine
    nullify(mat)
    if (allocated(tmp)) then
        deallocate(tmp)
    end if

    return
end subroutine Jacobi_diagonalization

subroutine Jacobi_rotate(n, p, q, mat, v, vecp, vecq, status)
    implicit none
    integer, intent(in) :: n, p, q
    double precision, dimension(n, n), intent(inout) :: mat, v
    double precision, dimension(n), intent(inout) :: vecp, vecq
    integer, intent(out) :: status

    double precision :: theta
    double precision :: c, s
    double precision :: pp, pq, qq

    ! quick return if wrong in input
    status = 0
    if (p < 1 .or. p > n .or. q < 1 .or. q > n) then
        status = -3
        return
    end if

    ! get rotation angle and sin as well as cos values
    pp = mat(p, p)
    pq = mat(p, q)
    qq = mat(q, q)
    if (abs(pp - qq) > 0.d0) then
        theta = atan(2.d0 * pq / (pp - qq)) / 2.d0
    else
        theta = atan(1.d0)
    end if
    c = cos(theta)
    s = sin(theta)

    ! update eigenvectors
    vecp = v(:, p)
    vecq = v(:, q)
    v(:, p) = vecq * s + vecp * c
    v(:, q) = vecq * c - vecp * s

    ! update target matrix
    vecp = mat(p, :)
    vecq = mat(q, :)
    mat(p, :) = vecq * s + vecp * c
    mat(q, :) = vecq * c - vecp * s
    vecp = mat(:, p)
    vecq = mat(:, q)
    mat(:, p) = vecq * s + vecp * c
    mat(:, q) = vecq * c - vecp * s
    mat(p, p) = pp * c ** 2 + qq * s ** 2 + 2.d0 * pq * s * c
    mat(q, q) = pp * s ** 2 + qq * c ** 2 - 2.d0 * pq * s * c
    mat(p, q) = (qq - pp) * s * c + pq * (c ** 2 - s ** 2)
    ! mat(p, q) = 0.d0
    ! mat(q, p) = mat(p, q)

    return
end subroutine Jacobi_rotate

subroutine sort_eigen(n, w, v)
    implicit none
    integer, intent(in) :: n
    double precision, dimension(n), intent(inout) :: w
    double precision, dimension(n, n), intent(inout) :: v

    integer :: i
    integer :: dad, son
    integer, dimension(:), allocatable :: ind
    double precision, dimension(:), allocatable :: swp
    integer :: pos

    ! initialize the sort sequence, allocate memory for swap
    allocate(ind(n))
    forall (i = 1:n)
        ind(i) = i
    end forall
    allocate(swp(n))

    ! heap sort of w, and get the sort sequence
    do i = n / 2 - 1, 0, -1
        dad = i + 1
        son = dad * 2
        do while (son <= n)
            if (son <= n - 1) then
                if (w(son) < w(son + 1)) then
                    son = son + 1
                end if
            end if
            if (w(son) <= w(dad)) then
                exit
            else
                call swap_double(w(dad), w(son))
                call swap_int(ind(dad), ind(son))
                dad = son
                son = dad * 2
            end if
        end do
    end do
    do i = n - 1, 1, -1
        call swap_double(w(1), w(i + 1))
        call swap_int(ind(1), ind(i + 1))
        dad = 1
        son = dad * 2
        do while (son <= i)
            if (son <= i - 1) then
                if (w(son) < w(son + 1)) then
                    son = son + 1
                end if
            end if
            if (w(son) <= w(dad)) then
                exit
            else
                call swap_double(w(dad), w(son))
                call swap_int(ind(dad), ind(son))
                dad = son
                son = dad * 2
            end if
        end do
    end do

    ! reorder v according to ind
    i = 1
    do while (i < n)
        if (ind(i) == i) then
            i = i + 1
            cycle
        end if
        pos = i
        swp(:) = v(:, i)
        do while (ind(i) /= pos)
            v(:, i) = v(:, ind(i))
            call swap_int(i, ind(i))
        end do
        v(:, i) = swp(:)
        ind(i) = i
        i = pos + 1
    end do

    ! release memory
    deallocate(ind)
    deallocate(swp)

    return
end subroutine sort_eigen

subroutine swap_double(a, b)
    double precision, intent(inout) :: a, b

    double precision :: t

    t = a
    a = b
    b = t

    return
end subroutine swap_double

subroutine swap_int(a, b)
    integer, intent(inout) :: a, b

    integer :: t

    t = a
    a = b
    b = t

    return
end subroutine swap_int

