module functions
    use params
    implicit none
    
    contains

    subroutine mapping(x, y, a, b, c, m)
        implicit none
        real(wp), intent(inout) :: x, y
        real(wp), intent(in) :: a, b, c, m

        y = y - b * sin(2.0_wp * PI * x) - c * sin(2.0_wp * PI * m * x)
        x = modulo(x + a * (1.0_wp - y ** 2), 1.0_wp)

        return

    end subroutine mapping
    
    subroutine mapping_ensemble(x, y, a, b, c, m)
        implicit none
        real(wp), intent(inout) :: x(:), y(:)
        real(wp), intent(in) :: a, b, c, m
        integer :: i, num_ic

        num_ic = size(x)
        !$omp parallel do private(i)
        do i = 1, num_ic
            call mapping(x(i), y(i), a, b, c, m)
        end do
        !$omp end parallel do

        return
    end subroutine mapping_ensemble

    subroutine mapping_backward(x, y, a, b, c, m)
        implicit none
        real(wp), intent(inout) :: x, y
        real(wp), intent(in) :: a, b, c, m

        x = modulo(x - a * (1 - y ** 2), 1.0_wp)
        y = y + b * sin(2.0_wp * PI * x) + c * sin(2.0_wp * PI * m * x)

        return

    end subroutine mapping_backward

    subroutine rotation_number(x0, y0, a, b, c, m, N, rn)
        use, intrinsic :: ieee_arithmetic
        implicit none
        real(wp), intent(in) :: x0, y0, a, b, c, m
        integer, intent(in) :: N
        real(wp), intent(out) :: rn
        real(wp) :: x_old, y_old, x_new, y_new, S, w(N), u, aux_rn(2)
        integer :: i
    
        ! Initialize variables
        x_old = x0
        y_old = y0
        aux_rn = 0.0_wp
    
        ! Compute u and w
        do i = 1, N
            u = real(i, wp) / real(N, wp)
            if (u == 0.0_wp .or. u == 1.0_wp) then
                w(i) = 0.0_wp
            else
                w(i) = exp(-1.0_wp / (u * (1.0_wp - u)))
            end if
        end do
    
        ! Normalize w
        S = sum(w)
        w = w / S
    
        ! Compute the rotation number
        do i = 1, N
            x_new = x_old
            y_new = y_old
            call mapping(x_new, y_new, a, b, c, m)  ! Assuming mapping is a subroutine
            aux_rn(1) = aux_rn(1) + w(i) * modulo(x_new - x_old, 1.0_wp)
            aux_rn(2) = aux_rn(2) + w(i) * ( a * (1 - y_new ** 2))
            x_old = x_new  ! Update x_old for the next iteration
            y_old = y_new  ! Update y_old for the next iteration
        end do

        if (abs(aux_rn(1) - aux_rn(2)) < 1e-15) then
            rn = aux_rn(1)
        else
            rn = ieee_value(rn, ieee_quiet_nan)
        end if
    
    end subroutine rotation_number

    subroutine is_fixed_point(x0, y0, a, b, c, m, period, eps, ifp)
        implicit none
        real(wp), intent(in) :: x0, y0, a, b, c, m, eps
        integer, intent(in) :: period
        logical, intent(out) :: ifp
        integer :: i
        real(wp) :: x, y

        ! Initialize initial condition
        x = x0
        y = y0

        ! Iterate the initial condition
        do i = 1, period
            call mapping(x, y, a, b, c, m)
        end do

        ifp = .false.

        if (abs(x - x0) < eps .and. abs(y - y0) < eps) then
            ifp = .true.
        end if

        return

    end subroutine is_fixed_point

    subroutine find_fixed_points(x_range, y_range, a, b, c, m, period, eps, grid_size, fixed_points, num_fixed_points)
        implicit none
        real(wp), intent(in) :: x_range(:), y_range(:), a, b, c, m, eps
        integer, intent(in) :: period, grid_size
        real(wp), intent(out) :: fixed_points(:, :)
        integer, intent(out) :: num_fixed_points
        integer :: i, j, count
        real(wp) :: x, y
        logical :: ifp

        fixed_points = 0.0_wp
        count = 0
        do i = 0, grid_size
            x = x_range(1) + (x_range(2) - x_range(1)) * i / grid_size
            do j = 0, grid_size
                y = y_range(1) + (y_range(2) - y_range(1)) * j / grid_size
                call is_fixed_point(x, y, a, b, c, m, period, eps, ifp)
                if (ifp .eqv. .true.) then
                    fixed_points(count + 1, 1) = x
                    fixed_points(count + 1, 2) = y
                    count = count + 1
                end if
            end do
        end do

        num_fixed_points = count

        return

    end subroutine find_fixed_points

    subroutine time_series(x_ini, y_ini, a, b, c, m, N, u)
        implicit none
        real(wp), intent(in) :: x_ini, y_ini, a, b, c, m
        integer, intent(in) :: N
        real(wp), intent(out) :: u(:, :)
        integer :: i
        real(wp) :: x, y

        x = x_ini
        y = y_ini

        u(1, 1) = x
        u(1, 2) = y
        do i = 2, N
            call mapping(x, y, a, b, c, m)
            u(i, 1) = x
            u(i, 2) = y
        end do

        return
        
    end subroutine time_series

    subroutine time_series_backward(x_ini, y_ini, a, b, c, m, N, u)
        implicit none
        real(wp), intent(in) :: x_ini, y_ini, a, b, c, m
        integer, intent(in) :: N
        real(wp), intent(out) :: u(:, :)
        integer :: i
        real(wp) :: x, y

        x = x_ini
        y = y_ini

        u(1, 1) = x
        u(1, 2) = y
        do i = 2, N
            call mapping_backward(x, y, a, b, c, m)
            
            u(i, 1) = x
            u(i, 2) = y
        end do

        return
        
    end subroutine time_series_backward

    subroutine jacobian(x, y, a, b, c, m, J)
        implicit none
        real(wp), intent(in) :: x, y, a, b, c, m
        real(wp), intent(out) :: J(:, :)
        real(wp) :: y_new

        y_new = y - b * sin(2.0_wp * PI * x) - c * sin(2.0_wp * PI * m * x)

        J(2, 1) = -2.0_wp * PI * b * cos(2.0_wp * PI * x) - 2.0_wp * PI * c * m * cos(2.0_wp * PI * m * x)
        J(2, 2) = 1.0_wp
        J(1, 1) = 1.0_wp - 2.0_wp * a * y_new * J(2, 1)
        J(1, 2) = - 2.0_wp * a * y_new * J(2, 2)

        return

    end subroutine jacobian

    subroutine eigenvalues_and_eigenvectors(x, y, a, b, c, m, period, eigvals, eigvectors)
        implicit none
        real(wp), intent(in) :: x, y, a, b, c, m
        integer, intent(in) :: period
        real(wp), intent(out) :: eigvals(2)
        real(wp), intent(out) :: eigvectors(2, 2)
        integer :: i
        real(wp) :: ts(period, 2)
        real(wp) :: J(2, 2), aux_J(2, 2)
        real(wp) :: trace, determinant, discriminant, norm

        ! Initialize J as the identity matrix
        J(1, :) = [1.0_wp, 0.0_wp]
        J(2, :) = [0.0_wp, 1.0_wp]
        call time_series(x, y, a, b, c, m, period, ts)

        do i = 1, period
            call jacobian(ts(i, 1), ts(i, 2), a, b, c, m, aux_J)
            J = matmul(aux_J, J)
        end do
        ! Compute trace and determinant
        trace = J(1, 1) + J(2, 2)
        determinant = J(1, 1) * J(2, 2) - J(1, 2) * J(2, 1)
        ! Compute discriminant
        discriminant = trace**2 - 4.0_wp * determinant

        ! Check for complex eigenvalues
        if (discriminant < 0.0_wp) then
            print *, "Negative discriminant. Complex eigenvalues."
            stop
        end if

        ! Compute the eigenvalues
        eigvals(1) = (trace + sqrt(discriminant)) / 2.0_wp
        eigvals(2) = (trace - sqrt(discriminant)) / 2.0_wp
        
        ! Compute the eigenvectors        
        do i = 1, 2
            if (J(1, 2) /= 0.0_wp) then
                eigvectors(:, i) = [J(1, 2), eigvals(i) - J(1, 1)]
            else if (J(2, 1) /= 0.0_wp) then
                eigvectors(:, i) = [eigvals(i) - J(2, 2), J(2, 1)]
            else
                eigvectors(:, i) = [real(i - 1, wp), real(2 - i, wp)]  ! Trivial eigenvectors for diagonal matrix
            end if

            ! Normalize the eigenvector
            norm = sqrt(eigvectors(1, i)**2 + eigvectors(2, i)**2)
            eigvectors(:, i) = eigvectors(:, i) / norm
        end do

        return

    end subroutine eigenvalues_and_eigenvectors

    subroutine calculate_manifold(x_PO, y_PO, a, b, c, m, eigenvectors, branch, stability, N, num_ic, eps, w)
        implicit none
        real(wp), intent(in) :: x_PO, y_PO, a, b, c, m, eps
        real(wp), intent(in) :: eigenvectors(2, 2)
        character(len=*), intent(in) :: branch, stability
        integer, intent(in) :: N, num_ic
        real(wp), intent(out) :: w(:, :, :)
        integer :: i
        real(wp) :: y_ini, y_end, x(num_ic), y(num_ic)
        real(wp) :: theta
        real(wp) :: ts(N, 2)

        ! Determine theta based on stability
        if (stability == "stable") then
            theta = atan2(eigenvectors(2, 2), eigenvectors(1, 2))
            theta = modulo(theta, PI)
        else if (stability == "unstable") then
            theta = atan2(eigenvectors(2, 1), eigenvectors(1, 1))
            theta = modulo(theta, PI)
        else
            print *, "Invalid input: ", stability, ". Expected 'stable' or 'unstable'."
            stop
        end if

        ! Determine y_ini and y_end based on branch
        if (branch == "upper") then
            y_ini = y_PO
            y_end = y_PO + sin(theta) * eps
        else if (branch == "lower") then
            y_ini = y_PO - sin(theta) * eps
            y_end = y_PO
        else
            print *, "Invalid input: ", branch, ". Expected 'upper' or 'lower'."
            stop
        end if

        call linspace(y, y_ini, y_end, num_ic)
        ! call logspace(y, log10(y_ini), log10(y_end), num_ic)
        x = (y - y_PO) / tan(theta) + x_PO

        ! Compute the manifold
        do i = 1, num_ic
            if (stability == "stable") then
                call time_series_backward(x(i), y(i), a, b, c, m, N, ts)
            else
                call time_series(x(i), y(i), a, b, c, m, N, ts)
            end if
            w(i, :, 1) = ts(:, 1)
            w(i, :, 2) = ts(:, 2)
        end do

        return

    end subroutine calculate_manifold

    subroutine check_for_escape(x0, y0, a, b, c, m, N, esc_pos, has_escaped)
        implicit none
        real(wp), intent(in) :: x0, y0, a, b, c, m, esc_pos
        integer, intent(in) :: N
        integer, intent(out) :: has_escaped
        integer :: i
        real(wp) :: x, y

        x = x0
        y = y0

        has_escaped = 0
        do i = 1, N
            call mapping(x, y, a, b, c, m)
            if ((esc_pos > 0 .and. y > esc_pos) .or. (esc_pos < 0 .and. y < esc_pos)) then
                has_escaped = 1
                return
            end if
        end do

        return

    end subroutine check_for_escape

    subroutine final_SALI(x0, y0, a, b, c, m, N, sali)
        implicit none
        real(wp), intent(in) :: x0, y0, a, b, c, m
        integer, intent(in) :: N
        real(wp), intent(out) :: sali
        real(wp) :: x, y, J(2, 2), v1(2), v2(2), PAI, AAI
        integer :: i
        ! Initialize main trajectory
        x = x0
        y = y0

        ! Initialize two deviation vectors
        v1 = (/ 1.0_wp, 0.0_wp /)
        v2 = (/ 0.0_wp, 1.0_wp /)

        do i = 1, N
            ! Evolve the main trajectory
            call mapping(x, y, a, b, c, m)

            ! Jacobian
            call jacobian(x, y, a, b, c, m, J)

            ! Evolve the deviation vectors
            v1 = matmul(J, v1)
            v2 = matmul(J, v2)

            ! Normalize the deviation vectors
            v1 = v1 / norm2(v1)
            v2 = v2 / norm2(v2)

            ! Compute the aligment indices
            PAI = norm2(v1 + v2)
            AAI = norm2(v1 - v2)

            ! Compute SALI
            sali = min(PAI, AAI)

            if (sali < TOL) then
                return
            end if

        end do

        return

    end subroutine final_SALI

    subroutine get_manifolds_times(m, times)
        implicit none
        real(wp), intent(in) :: m
        integer, intent(out) :: times(:, :, :)
        ! Chain - Stability - Branch
        !                    usu  usl  uuu  uul | lsu  lsl  luu  lul
        if (m == 1.0_wp) then
            times = reshape([180, 1, 1, 200, 200, 1, 180, 1], shape(times))
            !times = reshape([1, 100, 1, 1, 1, 1, 1, 100], shape(times))
            !times = 110
        elseif (m == -1.0_wp) then
            times = 80
        elseif (m == 2.0_wp) then
            times = reshape([120, 100, 110, 120, 130, 90, 140, 85], shape(times))
        elseif (m == -2.0_wp) then
            times = reshape([120, 100, 120, 90, 125, 100, 120, 100], shape(times))
        elseif (m == 3.0_wp) then
            times = reshape([120, 100, 120, 90, 125, 100, 120, 100], shape(times))
        elseif (m == -3.0_wp) then
            times = reshape([150, 120, 120, 150, 150, 120, 150, 120], shape(times))
        elseif (m == 4.0_wp) then
            times = reshape([120, 100, 110, 110, 125, 110, 120, 110], shape(times))
        elseif (m == -4.0_wp) then
            times = reshape([120, 110, 110, 120, 120, 105, 120, 100], shape(times))
        elseif (m == 0.5_wp) then
            times = reshape([150, 1, 170, 1, 1, 170, 150, 1], shape(times))
        elseif (m == -0.5_wp) then
            times = reshape([110, 1, 170, 1, 1, 170, 150, 1], shape(times))
            times = reshape([80, 150, 80, 170, 170, 80, 110, 150], shape(times))
        elseif (m == 1.5_wp) then
            times = reshape([150, 90, 95, 120, 120, 110, 150, 90], shape(times))
        elseif (m == -1.5_wp) then
            times = reshape([110, 110, 110, 130, 140, 100, 110, 110], shape(times))
        elseif (m == 2.5_wp) then
            times = reshape([160, 100, 100, 120, 120, 110, 130, 100], shape(times))
        elseif (m == -2.5_wp) then
            times = reshape([120, 100, 100, 140, 150, 100, 120, 110], shape(times))
        elseif (m == 3.4_wp) then
            times = reshape([120, 100, 100, 140, 150, 100, 120, 110], shape(times))
        elseif (m == 3.5_wp) then
            times = reshape([120, 100, 100, 140, 150, 100, 120, 110], shape(times))
        elseif (m == 3.8_wp) then
            times = reshape([120, 100, 100, 140, 150, 100, 120, 110], shape(times))    
        elseif (m == -3.6_wp) then
            times = reshape([120, 100, 100, 140, 150, 100, 120, 110], shape(times))
        elseif (m == -9.0_wp) then
            times = reshape([120, 100, 100, 140, 150, 100, 120, 110], shape(times))
        else
            times = reshape([150, 110, 120, 140, 145, 130, 140, 120], shape(times))
        end if
    end subroutine get_manifolds_times

    subroutine linspace(result, start, stop, num)
        real(wp), intent(in) :: start      ! Starting value
        real(wp), intent(in) :: stop       ! Ending value
        integer, intent(in) :: num     ! Number of points
        real(wp), dimension(:), intent(out) :: result  ! Output array

        ! Calculate the spacing between points
        real(wp) :: step
        integer :: i

        if (num < 2) then
            print *, "Number of points must be at least 2."
            stop
        endif

        step = (stop - start) / real(num - 1)

        ! Fill the array with linearly spaced values
        do i = 1, num
            result(i) = start + (i - 1) * step
        end do
    end subroutine linspace

    subroutine logspace(result, start_exp, stop_exp, num)
        ! Generate logarithmically spaced values between 10**start_exp and 10**stop_exp
        real(wp), intent(in) :: start_exp      ! Exponent of starting value (10^start_exp)
        real(wp), intent(in) :: stop_exp       ! Exponent of ending value (10^stop_exp)
        integer, intent(in) :: num             ! Number of points
        real(wp), dimension(:), intent(out) :: result  ! Output array
    
        real(wp) :: step, exponent
        integer :: i
    
        if (num < 2) then
            print *, "Number of points must be at least 2."
            stop
        endif
    
        step = (stop_exp - start_exp) / real(num - 1)
    
        do i = 1, num
            exponent = start_exp + (i - 1) * step
            result(i) = 10.0_wp ** exponent
        end do
    end subroutine logspace

    function should_process(i) result(process)
        implicit none
        integer(kind=8), intent(in) :: i
        logical :: process
    
        process = .false.
    
        if (i <= 1e4) then
            if (mod(i - 1, int(2)) == 0) process = .true.
        else if (i > 1e4 .and. i <= 1e5) then
            if (mod(i, int(5e2)) == 0) process = .true.
        else if (i > 1e5 .and. i <= 1e6) then
            if (mod(i, int(5e3)) == 0) process = .true.
        else if (i > 1e6 .and. i <= 1e7) then
            if (mod(i, int(5e4)) == 0) process = .true.
        else if (i > 1e7 .and. i <= 1e8) then
            if (mod(i, int(5e5)) == 0) process = .true.
        else if (i > 1e8 .and. i <= 1e9) then
            if (mod(i, int(5e6)) == 0) process = .true.
        else if (i > 1e9 .and. i <= 1e10) then
            if (mod(i, int(5e7)) == 0) process = .true.
        end if
    end function should_process

end module functions