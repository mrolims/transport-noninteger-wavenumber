program main
    use params
    use functions
    implicit none
    ! --- Parameters --- !
    real(wp), parameter :: a = 0.805_wp
    real(wp), parameter :: b = 0.597_wp
    real(wp), parameter :: c = 0.005_wp
    real(wp), parameter :: m = -9.0_wp
    integer, parameter :: num_y = 1000
    integer, parameter :: N = 1e8
    ! --- Variables --- !
    integer :: i 
    real(wp) :: x
    real(wp) :: y, y_ini, y_end
    real(wp) :: rn
    character :: path*100, datafile*200

    path = "Data/"

    1000 format(a, "rotation_number_m=", f0.5, "_a=", f0.5, "_b=", f0.5, "_c=", f0.5, "_N=", i0, ".dat")
    write(unit=datafile, fmt=1000)trim(path), m, a, b, c, N
    open(10, file=trim(datafile))

    x = 0.5_wp
    y_ini = 0.21_wp
    y_end = 0.225_wp

    !$omp parallel do ordered schedule(dynamic) private(i, y, rn)
    do i = 1, num_y
        y = y_ini + (y_end - y_ini) * i / num_y
        call rotation_number(x, y, a, b, c, m, N, rn)
        !$omp ordered
        write(10, *) y, rn
        !$omp end ordered
    end do
    !$omp end parallel do
    close(10)
    
end program main