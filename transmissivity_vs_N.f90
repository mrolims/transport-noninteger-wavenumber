program main
    use params
    use functions
    implicit none
    ! --- Parameters --- !
    real(wp), parameter :: a = 0.805_wp
    real(wp), parameter :: b = 0.597_wp
    real(wp), parameter :: c = 0.005_wp
    ! --- Variables --- !
    integer(8) :: i, j, N
    integer :: num_ic, num_remaining
    real(wp), allocatable :: x(:)
    real(wp), allocatable :: y(:)
    real(wp), allocatable :: x_temp(:)
    real(wp), allocatable :: y_temp(:)
    real(wp) :: m
    real(wp) :: y0
    real(wp) :: esc_y
    real(wp) :: T
    real(wp) :: sali
    logical, allocatable :: mask(:)
    character :: path*100, datafile*200, arg*20
    ! --- Seed variables --- !
    integer :: ns ! seed variables
    integer, allocatable :: seed(:) ! seed array

    path = "Data/"

    call getarg(1, arg)
    read(arg, *) num_ic
    call getarg(2, arg)
    read(arg, *) m
    call getarg(3, arg)
    read(arg, *) y0
    call getarg(4, arg)
    read(arg, *) esc_y
    call getarg(5, arg)
    read(arg, *) N

    write(*,*)""
    write(*, "(a, i0)")"num_ic = ", num_ic
    write(*, "(a, f0.5)")"m = ", m
    write(*, "(a, f0.5)")"y0 = ", y0
    write(*, "(a, i0)")"N = ", N
    write(*, "(a, f0.5)")"esc_y = ", esc_y

    1000 format(a, "transmissivity_vs_N_m=", f0.5, "_y0=", f0.5, "_nic=", i0, "_N=", i0, "_escy=", f0.5, ".dat")
    write(unit=datafile, fmt=1000) trim(path), m, y0, num_ic, N, esc_y
    open(10, file=trim(datafile))

    ! Initialize arrays
    allocate(x(num_ic), y(num_ic), mask(num_ic))

    ! Set the seed
    call random_seed(size=ns)
    allocate(seed(ns))
    seed = 13
    call random_seed(put=seed)
    deallocate(seed)

    ! Generate random chaotic initial conditions
    call random_number(x)
    y = y0
    do j = 1, num_ic
        do while (.true.)
            call final_SALI(x(j), y(j), a, b, c, m, 1000, sali)
            if (sali < 1e-12) then
                exit
            else
                call random_number(x(j))
            end if
        end do
    end do

    T = 0.0_wp
    do i = 1, N
        call mapping_ensemble(x, y, a, b, c, m)
        ! Check for escaped particles
        mask = (esc_y > 0.0_wp .and. y > esc_y) .or. (esc_y < 0.0_wp .and. y < esc_y)
        
        ! Update the number of escaped particles
        T = T + real(sum(merge(1, 0, mask)), wp)

        ! Count the number of remaining particles
        num_remaining = count(.not. mask)

        ! Pack the arrays to remove escaped particles
        if (allocated(x_temp)) deallocate(x_temp)
        if (allocated(y_temp)) deallocate(y_temp)
        allocate(x_temp(num_remaining), y_temp(num_remaining))

        x_temp = pack(x, .not. mask)
        y_temp = pack(y, .not. mask)

        ! Resize the original arrays (optional, if you want to keep them)
        if (allocated(x)) deallocate(x)
        if (allocated(y)) deallocate(y)
        allocate(x(num_remaining), y(num_remaining))
        x = x_temp
        y = y_temp

        deallocate(x_temp, y_temp)  ! Clean up temporary arrays

        if (should_process(i)) then
            write(10, "(i10, f30.16)") i, T / real(num_ic, wp)
        end if
    end do
    
    close(10)
    
end program main
