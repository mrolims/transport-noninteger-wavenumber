program manifolds
    use functions
    implicit none
    ! --- Parameters --- !
    integer, parameter :: period = 11
    real(wp), parameter :: a = 0.805_wp
    real(wp), parameter :: b = 0.597_wp
    real(wp), parameter :: c = 0.005_wp
    real(wp), parameter :: eps = 1e-6
    ! --- Variables --- !
    integer :: i, j, k, l, ii
    integer :: num_ic
    integer :: times(2, 2, 2)
    real(wp) :: m
    real(wp) :: x, y
    real(wp) :: xy_coordinates(2, 2)
    real(wp) :: eigvals(2), eigvectors(2, 2)
    real(wp), allocatable, dimension(:, :, :) :: w
    character(len=10), dimension(2) :: chains = ["upper", "lower"]
    character(len=10), dimension(2) :: stabilities = ["stable  ", "unstable"]
    character(len=10), dimension(2) :: branches = ["upper", "lower"]
    character :: path*100, datafile*200, user*20, home*100, arg*16

    call get_environment_variable("USER", user)
    call get_environment_variable("HOME", home)

    if (trim(user) == "jdanilo") then
        path = trim(home) // "/Matheus/Pesquisa/extended-standard-nontwist-map/"
    else
        path = trim(home) // "/Pesquisa/extended-standard-nontwist-map/"
    end if

    ! --- Get the arguments --- !
    call getarg(1, arg)
    read(arg, *) m
    call getarg(2, arg)
    read(arg, *) num_ic

    ! --- Formats for the files used in the program --- !
    1000 format(a, "Data/", a, "_hyperbolic_point_a=", f0.5, "_b=", f0.5, "_c=", f0.5, "_m=", f0.5, "_period=", i0, ".dat")
    1100 format(a, "Data/", a , "_", a, "_manifold_a=", f0.5, "_b=", f0.5, "_c=", f0.5, "_m=", f0.5, "_period=", i0, "_nic=", i0, ".dat")

    ! --- Get the iteration times for each manifold --- !
    call get_manifolds_times(m, times)

    ! --- Extract the hyperbolic points --- !
    ! Upper chain
    write(unit=datafile, fmt=1000)trim(path), "upper", a, b, c, m, period
    open(10, file=trim(datafile), status="old", action="read")
    read(10, *) (xy_coordinates(1, j), j = 1, 2)
    close(10)

    ! Lower chain
    write(unit=datafile, fmt=1000)trim(path), "lower", a, b, c, m, period
    open(10, file=trim(datafile), status="old", action="read")
    read(10, *) (xy_coordinates(2, j), j = 1, 2)
    close(10)

    do i = 1, size(chains)
        x = xy_coordinates(i, 1)
        y = xy_coordinates(i, 2)
        call eigenvalues_and_eigenvectors(x, y, a, b, c, m, period, eigvals, eigvectors)
        do j = 1, size(stabilities)
            write(unit=datafile, fmt=1100) trim(path), trim(chains(i)), trim(stabilities(j)), a, b, c, m, period, num_ic
            open(11, file=trim(datafile))
            do k = 1, size(branches)
                allocate(w(num_ic, times(k, j, i), 2))
                print *, "Calculating the ", trim(branches(k)), " branch of the " , trim(stabilities(j)), " manifold of the ", trim(chains(i)), " chain..."
                call calculate_manifold(x, y, a, b, c, m, eigvectors, branches(k), trim(stabilities(j)), times(k, j, i), num_ic, eps, w)
                do l = 1, num_ic
                    do ii = 1, times(k, j, i)
                        if (w(l, ii, 1) >= 0.0_wp .and. w(l, ii, 2) <= 1.0_wp .and. w(l, ii, 2) >= -1.0_wp .and. w(l, ii, 2) <= 1.0_wp) then
                            write(11, "(*(f30.16))") w(l, ii, 1), w(l, ii, 2)
                        end if
                    end do
                end do
                deallocate(w)
            end do
            close(11)
        end do
    end do

    
    write(*, "(a)")"Finished."


end program manifolds
