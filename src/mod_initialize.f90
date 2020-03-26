module mod_initialize

implicit none

integer :: n_rotors,&    ! Number of rotors
           n_blades, &    ! Number of blades
           n_elements, &  ! Number of elements in each blade
           n_panels,&       ! Number of panels in the wake in the downstream direction
           n_boundaries  ! Number of boundaries in the domain


real(kind=8) :: rho, & ! Air density
                sound_speed, &
                tf, & ! Final time
                dt, & ! Time interval
                toll_influ, & ! toll distance
                toll_gamma ! toll on gamma

logical :: dyn_stall


real(kind=8), parameter :: pi = 3.14159265358979323846264338327950288419716939937510

real(kind=8), dimension(3) :: Vinf ! Air velocity at the infinite upstream
real(kind=8), dimension(:), allocatable :: center ! Coordinates of the centers of the rotors
real(kind=8), dimension(:), allocatable :: omega,&  ! Angular velocity
                                           R, &             ! Radius
                                           c, &             ! Chord
                                           bs, &            ! Adimensional length at which the blade start
                                           e, &             ! Dimensional hinge position
                                           width, &
                                           offset_hub, &
                                           fill_percent,&
                                           Vtip, &          ! Tip velocity
                                           alpha_sh_deg,  & ! AoA in degrees of the shaft with respect to the TPP
                                           theta_deg, &     ! Pitch angle of the blades (no twist)
                                           theta_s_deg, &
                                           alpha_sh, &
                                           theta, &
                                           theta_s, &
                                           beta_deg,&
                                           sigma, &         ! Solidity
                                           A,&             ! Area of the rotor
                                           initial_radius,&
                                           Ibeta, & ! Moment of inertia
                                           Sbeta, & ! Satic moment
                                           boundary_first_vert,&
                                           boundary_second_vert,&
                                           boundary_third_vert,&
                                           boundary_forth_vert


! Global declarations
real(kind=8), dimension(:, :), allocatable :: y,& ! vettore dimensionale delle cooridnate
                                                  ! dei punti di collocazione sulla pala
                                              beta
real(kind=8), dimension(:, :, :), allocatable :: boundaries

type wake_panel
   real(kind=8)                 :: gamma  ! Circulation in the panel
   real(kind=8), dimension(4,3) :: vertex   ! Coordinates of the 4 vertexes of the wake panel
   real(kind=8), dimension(2,3) :: vp ! Velocity associeted to the 3rd and 4th vertex of the panel
   real(kind=8), dimension(2,3) :: vp_old ! Velocity at the previous time step
   real(kind=8)                 :: radius ! Radius of the "vortex"
end type wake_panel

type segment
   real(kind=8)                 :: gamma  ! Circulation in the panel
   real(kind=8), dimension(2,3) :: vertex   ! Coordinates of the 4 vertexes of the wake panel
   real(kind=8)                 :: radius ! Radius of the "vortex"
 end type segment

type particle
   real(kind=8), dimension(3)   :: alpha_p ! Circulation associeted to the particle
   real(kind=8), dimension(3)   :: position
   real(kind=8), dimension(3)   :: vp ! Velocity associeted to the particle
   real(kind=8), dimension(3)   :: vp_old ! Velocity at the previous time step
   real(kind=8)                 :: radius,& ! Radius of the particle
                                   radius_old ! Old radius of the particle
   logical ::                      on_ground
end type particle

type blade_panel
   real(kind=8)                 :: gamma, alpha_old  ! Circulation in the blade
   real(kind=8), dimension(3)   :: collocation  ! Coordinates of the collocation point
   real(kind=8), dimension(4,3) :: vertex  ! Coordinates of the 4 vertexes of the blade panel
   real(kind=8)                 :: radius ! Radius of the vortex
   real(kind=8)                 :: dFz, dQ ! z-Force and torque of the element
end type blade_panel

contains

subroutine init()


  namelist /input/ &
    n_rotors,&    ! Number of rotors
    n_blades, &    ! Number of blades
    n_elements, &  ! Number of elements in each blade
    n_panels, &    ! Number of panels in the wake in the downstream direction

    n_boundaries

    open(8,file='vp_setUp.in', status='old',action='read',form='formatted')
    read(8,input)
    close(8)


endsubroutine init

subroutine init2()

  integer :: i

  namelist /input/ &
    center,&          ! Coordinates of the centers of the rotors
    R, &              ! Radius
    c, &              ! Chord
    bs, &             ! Adimensional length at which the blade start
    e, &              ! Dimensional position oh the flapping hinge
    width, &
    offset_hub, &
    fill_percent,&
    Vtip, &           ! Tip velocity
    rho , &           ! Air density
    sound_speed, &
    alpha_sh_deg, &   ! AoA in degrees of the shaft with respect to the TPP
    theta_deg, &
    theta_s_deg, &
    beta_deg,&
    Vinf, &
    Ibeta, &
    Sbeta, &

    tf, &
    dt, &
    toll_influ,&
    toll_gamma, &

    dyn_stall


    if(.not.allocated(center)) allocate(center(3*n_rotors))
    if(.not.allocated(omega)) allocate(omega(n_rotors))
    if(.not.allocated(R)) allocate(R(n_rotors))
    if(.not.allocated(c)) allocate(c(n_rotors))
    if(.not.allocated(bs)) allocate(bs(n_rotors))
    if(.not.allocated(e)) allocate(e(n_rotors))
    if(.not.allocated(width)) allocate(width(n_rotors))
    if(.not.allocated(offset_hub)) allocate(offset_hub(n_rotors))
    if(.not.allocated(fill_percent)) allocate(fill_percent(n_rotors))
    if(.not.allocated(Vtip)) allocate(Vtip(n_rotors))
    if(.not.allocated(alpha_sh_deg)) allocate(alpha_sh_deg(n_rotors))
    if(.not.allocated(theta_deg)) allocate(theta_deg(n_rotors))
    if(.not.allocated(theta_s_deg)) allocate(theta_s_deg(n_rotors))
    if(.not.allocated(alpha_sh)) allocate(alpha_sh(n_rotors))
    if(.not.allocated(theta)) allocate(theta(n_rotors))
    if(.not.allocated(theta_s)) allocate(theta_s(n_rotors))
    if(.not.allocated(beta_deg)) allocate(beta_deg(n_rotors))
    if(.not.allocated(sigma)) allocate(sigma(n_rotors))
    if(.not.allocated(Ibeta)) allocate(Ibeta(n_rotors))
    if(.not.allocated(Sbeta)) allocate(Sbeta(n_rotors))
    if(.not.allocated(A)) allocate(A(n_rotors))
    if(.not.allocated(initial_radius)) allocate(initial_radius(n_rotors))
    if(.not.allocated(beta))  allocate(beta(n_rotors, n_blades))

    open(8,file='vortex_particle.in', status='old',action='read',form='formatted')
    read(8,input)
    close(8)

    !omega(1) = Vtip/R             ! Angular velocity
    omega = Vtip/R             ! Angular velocity
    alpha_sh = alpha_sh_deg*pi/180
    theta = theta_deg*pi/180
    theta_s = theta_s_deg*pi/180
    sigma = n_blades*c/(pi*R) ! Solidity
    A = pi*R**2             ! Area of the rotor
    initial_radius = 0.05 * c

    do i = 1, n_rotors
      beta(i, :) = beta_deg(i) * pi / 180
    enddo

endsubroutine init2

subroutine init3()

  integer :: i

  namelist /input/ &
    boundary_first_vert,&
    boundary_second_vert,&
    boundary_third_vert,&
    boundary_forth_vert

  if(.not.allocated(boundaries)) allocate(boundaries(n_boundaries, 4, 3))
  if(.not.allocated(boundary_first_vert)) allocate(boundary_first_vert(3*n_boundaries))
  if(.not.allocated(boundary_second_vert)) allocate(boundary_second_vert(3*n_boundaries))
  if(.not.allocated(boundary_third_vert)) allocate(boundary_third_vert(3*n_boundaries))
  if(.not.allocated(boundary_forth_vert)) allocate(boundary_forth_vert(3*n_boundaries))

    open(9,file='boundary.in', status='old',action='read',form='formatted')
    read(9,input)
    close(9)


    do i = 1, n_boundaries

      boundaries(i, 1, :) = boundary_first_vert((i-1) * 3 + 1: i * 3)
      boundaries(i, 2, :) = boundary_second_vert((i-1) * 3 + 1: i * 3)
      boundaries(i, 3, :) = boundary_third_vert((i-1) * 3 + 1: i * 3)
      boundaries(i, 4, :) = boundary_forth_vert((i-1) * 3 + 1: i * 3)

    enddo



endsubroutine init3

end module mod_initialize
