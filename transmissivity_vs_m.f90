program main
    use params
    use functions
    implicit none
    ! --- External --- !
    external :: getenv
    ! --- Parameters --- !
    real(wp), parameter :: a = 0.805_wp
    real(wp), parameter :: b = 0.597_wp
    real(wp), parameter :: c = 0.005_wp
    ! --- Variables --- !
    integer :: j, N
    integer :: num_ic
    integer, allocatable :: has_escaped(:)
    real(wp), allocatable :: x(:)
    real(wp), allocatable :: y(:)
    real(wp) :: m, m_ini, m_end, dm
    real(wp) :: y0
    real(wp) :: esc_y
    real(wp) :: T
    real(wp) :: sali
    logical, allocatable :: mask(:)
    character :: path*100, datafile*200, arg*20
    ! --- Seed variables --- !
    integer :: ns, clock ! seed variables
    integer, allocatable :: seed(:) ! seed array

    path = "Data/"

    call getarg(1, arg)
    read(arg, *) num_ic
    call getarg(2, arg)
    read(arg, *) m_ini
    call getarg(3, arg)
    read(arg, *) m_end
    call getarg(4, arg)
    read(arg, *) dm
    call getarg(5, arg)
    read(arg, *) y0
    call getarg(6, arg)
    read(arg, *) esc_y
    call getarg(7, arg)
    read(arg, *) N

    write(*,*)""
    write(*, "(a, i0)")"num_ic = ", num_ic
    write(*, "(a, f0.5)")"m_ini = ", m_ini
    write(*, "(a, f0.5)")"m_end = ", m_end
    write(*, "(a, f0.5)")"dm = ", dm
    write(*, "(a, f0.5)")"y0 = ", y0
    write(*, "(a, f0.5)")"esc_y = ", esc_y
    write(*, "(a, i0)")"N = ", N

    1000 format(a, "transmissivity_vs_m_y0=", f0.5, "_m0=", f0.5, "_m1=", f0.5, "_dm=", f0.5, "_nic=", i0, "_N=", i0, "_escy=", f0.5, ".dat")
    write(unit=datafile, fmt=1000) trim(path), y0, m_ini, m_end, dm, num_ic, N, esc_y
    open(10, file=trim(datafile))

    ! Allocate arrays
    allocate(x(num_ic), y(num_ic), has_escaped(num_ic), mask(num_ic))

    m = m_ini
    do while (m < m_end + dm)

        ! Set the seed
        call random_seed(size=ns)
        allocate(seed(ns))
        
        ! Get the system clock
        call system_clock(count=clock)
        seed = clock + 37 * (/ (j - 1, j = 1, ns) /)
        
        ! Change the seed according to the system clock
        call random_seed(put=seed)
        deallocate(seed)
        
        ! Generate random initial conditions
        call random_number(x)
        
        ! Set the initial y values
        y = y0
        
        ! Check if all initial conditions are chaotic
        !$omp parallel do schedule(dynamic) private(sali)
        do j = 1, num_ic
            do while (.true.)
                call final_SALI(x(j), y(j), a, b, c, m, 1000, sali)
                if (sali <  1e-12) then
                    exit
                else
                    call random_number(x(j))
                end if
            end do
        end do
        !$omp end parallel do

        ! Check for escaped particles
        !$omp parallel do schedule(dynamic)
        do j = 1, num_ic
            call check_for_escape(x(j), y(j), a, b, c, m, N, esc_y, has_escaped(j))
        end do
        !$omp end parallel do
        mask = (has_escaped == 1)
        T = real(sum(merge(1, 0, mask)), kind=wp) / real(num_ic, kind=wp)

        ! Write the transmissivity to the file
        write(10, "(*(f30.16))") m, T
        
        ! Update m
        m = m + dm
    end do
    
    close(10)
    
end program main