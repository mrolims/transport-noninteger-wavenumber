program search_for_fixed_points
    use params
    use functions
    implicit none
    ! --- Parameters --- !
    integer, parameter :: period = 11
    integer, parameter :: grid_size = 1000
    real(wp), parameter :: a = 0.805_wp
    real(wp), parameter :: b = 0.597_wp
    real(wp), parameter :: c = 0.005_wp
    ! --- Variables --- !
    integer :: i, j, num_fixed_points, check_sali
    real(wp) :: m
    real(wp) :: fixed_point(2)
    real(wp) :: x_range(2), y_range(2)
    real(wp) :: fixed_points(grid_size ** 2, 2)
    real(wp) :: eps
    real(wp) :: sali
    character :: path*100, datafile*200, arg*33

    path = "Data/"

    call getarg(1, arg)
    read(arg, *) m
    call getarg(2, arg)
    read(arg, *) x_range(1)
    call getarg(3, arg)
    read(arg, *) x_range(2)
    call getarg(4, arg)
    read(arg, *) y_range(1)
    call getarg(5, arg)
    read(arg, *) y_range(2)
    call getarg(6, arg)
    read(arg, *) check_sali
    call getarg(7, arg)

    1000 format(a , a, "_elliptic_point_a=", f0.5, "_b=", f0.5, "_c=", f0.5, "_m=", f0.5, "_period=", i0, ".dat")
    write(unit=datafile, fmt=1000)trim(path), trim(arg), a, b, c, m, period

    write(*, "(a, f0.5)")"m = ", m
    write(*, "(a, f0.5)")"x_ini = ", x_range(1)
    write(*, "(a, f0.5)")"x_end = ", x_range(2)
    write(*, "(a, f0.5)")"y_ini = ", y_range(1)
    write(*, "(a, f0.5)")"y_end = ", y_range(2)
    write(*, "(a, a)")"region = ", trim(arg)
    write(*, "(a)")trim(datafile)

    eps = 2.0_wp / grid_size
    fixed_point = 0.0_wp
    do j = 1, 200
        call find_fixed_points(x_range, y_range, a, b, c, m, period, eps, grid_size, fixed_points, num_fixed_points)
        
        if (num_fixed_points == 0) then
            
            if (j == 1) fixed_point = -1.0_wp
            exit
        else
            fixed_point = 0.0_wp
        end if

        x_range(1) = huge(x_range(1))
        y_range(1) = huge(y_range(1))
        x_range(2) = -huge(x_range(2))
        y_range(2) = -huge(y_range(2))
        do i = 1, num_fixed_points
            if (fixed_points(i, 1) /= 0.0_wp .and. fixed_points(i, 2) /= 0.0_wp) then

                if (check_sali == 1) then
                    call final_SALI(fixed_points(i, 1), fixed_points(i, 2), a, b, c, m, 500, sali)
                    if (sali > 1e-4) then
                        fixed_point(1) = fixed_point(1) + fixed_points(i, 1)
                        fixed_point(2) = fixed_point(2) + fixed_points(i, 2)
                        if (fixed_points(i, 1) < x_range(1)) x_range(1) = fixed_points(i, 1)
                        if (fixed_points(i, 1) > x_range(2)) x_range(2) = fixed_points(i, 1)

                        if (fixed_points(i, 2) < y_range(1)) y_range(1) = fixed_points(i, 2)
                        if (fixed_points(i, 2) > y_range(2)) y_range(2) = fixed_points(i, 2)
                    end if
                else
                    fixed_point(1) = fixed_point(1) + fixed_points(i, 1)
                    fixed_point(2) = fixed_point(2) + fixed_points(i, 2)
                    if (fixed_points(i, 1) < x_range(1)) x_range(1) = fixed_points(i, 1)
                    if (fixed_points(i, 1) > x_range(2)) x_range(2) = fixed_points(i, 1)

                    if (fixed_points(i, 2) < y_range(1)) y_range(1) = fixed_points(i, 2)
                    if (fixed_points(i, 2) > y_range(2)) y_range(2) = fixed_points(i, 2)
                end if
            end if
        end do
        
        fixed_point = fixed_point / num_fixed_points

        if (abs(x_range(1) - x_range(2)) < TOL .and. abs(y_range(1) - y_range(2)) < TOL) then
            exit
        end if

        eps = eps / 2

    end do

    open(10, file=trim(datafile))
    write(10, "(*(f40.32))") fixed_point(1), fixed_point(2)
    close(10)

end program search_for_fixed_points
