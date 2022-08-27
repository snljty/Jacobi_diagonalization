program main
    implicit none

    integer :: argc, iarg, argl, argl_max
    character(len=:), allocatable :: argv0
    character(len=:), dimension(:), allocatable :: argv

    integer :: n
    integer, parameter :: n_def = 5

    integer :: ofl
    character(len=7), parameter :: ofl_name = "mat.dat"

    double precision, dimension(:, :), allocatable :: a

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

    ! allocate memory
    allocate(a(n, n))

    ! generate random matrix a
    call random_seed()
    call random_number(a)

    ! make a symmetric
    a = matmul(transpose(a), a)

    ! output matrix
    open(newunit = ofl, file = ofl_name, status = "replace", action = "write", &
        access = "direct", form = "unformatted", recl = 8 * n ** 2)
    write(ofl, rec = 1) a
    close(ofl)

    ! release memory
    deallocate(a)

    stop
end program main

