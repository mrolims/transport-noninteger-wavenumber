module params
    implicit none
    
    ! --- Quad precision --- !
    integer, parameter :: wp = selected_real_kind(33, 4931)
    real(wp), parameter :: TOL = 1e-31

    ! --- Pi constant --- !
    real(wp), parameter :: PI = acos(-1.0_wp)
    
end module params