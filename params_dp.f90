module params
    implicit none
    
    ! --- Double precision --- !
    integer, parameter :: wp = selected_real_kind(15, 307)
    real(wp), parameter :: TOL = 1e-15

    ! --- Pi constant --- !
    real(wp), parameter :: PI = acos(-1.0_wp)
    
end module params
