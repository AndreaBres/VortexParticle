 module mod_advectVortices

use mod_types
use mod_check

implicit none

real(kind=8), parameter :: k_radius = 0.095

interface translations_particles
  module procedure translations_particles_s, translations_particles_v
end interface translations_particles

interface translations_wake_panels
  module procedure translations_wake_panels_s, translations_wake_panels_v
end interface translations_wake_panels

interface radius_particle
  module procedure radius_particle_s, radius_particle_v
end interface radius_particle

interface radius_panel
  module procedure radius_panel_s, radius_panel_v
end interface radius_panel

contains

subroutine panel2particles(wake_panels, rotor,   new_particles, wake_segment)
!  wake_panels:   tabella contenente tutti i pannelli di scia
!  wake_segment:  tabella contenente la CIRCOLAZIONE di tutti i segmenti rimasti
!  new_particles: array contenente le nuove particelle
   type(wake_panel), dimension(:, :), intent(in)     ::  wake_panels
   type(segment),    dimension(:, :), intent(inout)  ::  wake_segment
   type(particle),   dimension(:),    intent(out)    ::  new_particles
   integer, intent(in) :: rotor

   real(kind=8) :: gamma
   real(kind=8), dimension(4, 3) :: vertex
   real(kind=8) :: radius
   real(kind=8), dimension(2, 3) :: vp_panel

   real(kind=8) :: gamma_p
   real(kind=8), dimension(4, 3) :: vertex_p
   real(kind=8) :: radius_p

   real(kind=8), dimension(3) :: alpha_p
   real(kind=8), dimension(3) :: position
   real(kind=8), dimension(3) :: vp
   real(kind=8), dimension(3) :: vp_old
                   !                1           2          3         4
                   ! in ordine : top-left   top-right, bot-right, bot-left
   real(kind=8), dimension(3) :: vertex_tl, vertex_tr, vertex_br, vertex_bl, &
                                 vertex_bl_p ! serve per l'elemento precedente
   real(kind=8), dimension(3) :: vertical_segment, right_semisegment, left_semisegment
   logical :: on_ground

   integer :: i, k


   do i = 1, size(wake_panels, 1) ! n_blades
      do k = 1, size(wake_panels, 2) + 1 ! n_elements + 1
         if (k==1) then
            call wake_panel_deconstructor(wake_panels(i, k),  gamma, vertex, radius, vp_panel)

            vertex_tl = vertex(1, :)
            vertex_bl = vertex(4, :)
            vertex_br = vertex(3, :)
            vertical_segment = vertex_bl - vertex_tl
            right_semisegment = (vertex_br - vertex_bl)/2

            alpha_p = gamma * (vertical_segment + right_semisegment) &
                    - wake_segment(i, k)%gamma * right_semisegment
            position = vertex_bl
            vp_old = vp_panel(2, :) ! bot-left (4)

         else if (k == size(wake_panels, 2) + 1) then
            call wake_panel_deconstructor(wake_panels(i, k-1),  gamma, vertex, radius, vp_panel)

            vertex_tr = vertex(2,:)
            vertex_br = vertex(3,:)
            vertex_bl = vertex(4,:)
            vertical_segment = vertex_tr - vertex_br
            left_semisegment = (vertex_br - vertex_bl)/2

            alpha_p = gamma * (vertical_segment + left_semisegment) &
                    - wake_segment(i, k-1)%gamma * left_semisegment
            position = vertex_br
            vp_old = vp_panel(1, :) ! bot-right (3)

         else
            call wake_panel_deconstructor(wake_panels(i, k-1),   gamma_p, vertex_p, radius_p)
            call wake_panel_deconstructor(wake_panels(i, k),   gamma, vertex, radius, vp_panel)

            vertex_tl = vertex(1,:)
            vertex_bl = vertex(4,:)
            vertex_br = vertex(3,:)
            vertex_bl_p = vertex_p(4,:)
            vertical_segment = vertex_bl - vertex_tl
            left_semisegment = (vertex_bl - vertex_bl_p)/2
            right_semisegment = (vertex_br - vertex_bl)/2

            alpha_p = gamma * (right_semisegment + vertical_segment) &
                    + gamma_p * (left_semisegment - vertical_segment) &
                    - wake_segment(i, k)%gamma * right_semisegment &
                    - wake_segment(i, k-1)%gamma * left_semisegment
            position = vertex_bl
            vp_old = vp_panel(2, :) ! bot-left (4)
         endif

         vp = 0
         on_ground = .false.
         call particle_constructor(alpha_p, position, vp, vp_old, initial_radius(rotor), initial_radius(rotor), on_ground,  &
              new_particles((i-1) * (size(wake_panels, 2) + 1) + k))
      enddo
   enddo

   ! Update of the segments
   do i = 1,size(wake_panels, 1)
      do k=1, size(wake_panels, 2)

        call wake_segment_constructor(wake_panels(i, k)%gamma, wake_panels(i, k)%vertex(3:4, :),&
                                      wake_panels(i, k)%radius,  wake_segment(i, k))

      enddo
   enddo

end subroutine panel2particles


  !!!!! INUTILIZZABILE VETTORIALMENTE PERCHè SEMPLICEMENTE DO SCALARE !!!!!!!!!!!!!
    subroutine translations_particles_v(particles, blade_extremes, distance, check)

        real(kind=8), dimension(:, :, :), intent(in)    :: blade_extremes
        type(particle), dimension(:), intent(inout) :: particles
        real(kind=8), dimension(:), intent(out) :: distance
        logical, intent(out) :: check


        real(kind=8), dimension(size(particles, 1), 3) :: position, position_n
        real(kind=8), dimension(size(particles, 1), 3) :: vp, vp_old

        integer :: k

        do k=1, size(particles,1)

           position(k, :)=particles(k)%position
           vp(k, :)=particles(k)%vp
           vp_old(k, :)=particles(k)%vp_old

         enddo

           position_n=position+(3./2*vp-1./2*vp_old)*dt

          !!!!!!!!!! OCCCHIO, MANCA L'IMPATTO CON PALA!!!!!!!
         !do k = 1, size(particles, 1)
        !   call check_blade(position(k, :), position_n(k, :),    blade_extremes)
         !enddo

        do k=1, size(particles,1)

           particles(k)%position=position_n(k, :)
           particles(k)%vp_old=vp(k, :)

        enddo


    end subroutine translations_particles_v


    subroutine translations_particles_s(particles, blade_extremes, rot, rig, col, t, distance, check)

           real(kind=8), dimension(:, :, :, :), intent(in)    :: blade_extremes
           real(kind=8), intent(in) :: t
           type(particle), intent(inout) :: particles
           real(kind=8), intent(out) :: distance
           logical, intent(out) :: check
           integer, intent(in) :: rot, rig, col

           real(kind=8), dimension(3) :: position, position_n
           real(kind=8), dimension(3) :: vp, vp_old
           logical :: chec_blade, chec_hub
           integer :: i, cont = 0, cont_1 = 0, cont_2

           position=particles%position
           vp=particles%vp
           vp_old=particles%vp_old
           position_n=position+(3./2*vp-1./2*vp_old)*dt

           ! Calculating distance travelled
           distance = norm2(position - position_n)

           ! cont_2 = 0
           ! if(.not.particles%on_ground) cont_2 = 1;

           check = .false.

           call check_blade(position, blade_extremes, particles%radius, particles%radius_old,&
                            position_n, chec_blade, distance)

           call check_hub(position, particles%radius, particles%radius_old,&
                          position_n, chec_hub, distance)

           ! Continue cycling if there's an intersection
           !!!!! IT HAS TO BE CORRECTED IF THE ENDING POINT IS OUTSIDE THE VOLUME !!!!
           !do while(check)
           !cont = 0;  ! Taking track of an intersection

           !do i = 1, size(boundaries, 1)
          !   call check_boundary(position, boundaries(i, :, :), particles%radius,&
          !                particles%radius_old,   position_n, check, distance)
           !enddo

           if(chec_hub .or. chec_blade) then
             particles%on_ground = .true.
             check = .true.
           !   if(cont_2 == 1)   particles%alpha_p = -particles%alpha_p;
           ! else
           !   if(cont_2 == 0)   particles%alpha_p = -particles%alpha_p;
           endif


           particles%position=position_n
           particles%vp_old=vp


    end subroutine translations_particles_s


    subroutine translations_wake_panels_s(panels)

           type(wake_panel) , intent(inout) :: panels

           real(kind=8), dimension(4,3) :: vertex
           real(kind=8), dimension(2,3) :: vp, vp_old
           real(kind=8), dimension(3) :: vp3, vp3_old, vp4, vp4_old

           real(kind=8), dimension(3) :: position_bl, position_br, position_bl_n, position_br_n

           integer :: k

           vertex=panels%vertex

           vp3=panels%vp(1,:)
           vp3_old=panels%vp_old(1,:)
           vp4=panels%vp(2,:)
           vp4_old=panels%vp_old(2,:)

           position_bl=vertex(4,:)
           position_br=vertex(3,:)

           position_bl_n=position_bl+(3./2*vp4-1./2*vp4_old)*dt
           position_br_n=position_br+(3./2*vp3-1./2*vp3_old)*dt

           panels%vertex(4,:)=position_bl_n
           panels%vertex(3,:)=position_br_n
           panels%vp_old(1,:)=vp3
           panels%vp_old(2,:)=vp4

    end subroutine translations_wake_panels_s



         subroutine translations_wake_panels_v(panels)

                type(wake_panel), dimension(:), intent(inout) :: panels

                real(kind=8), dimension(size(panels,1),4,3) :: vertex
                real(kind=8), dimension(size(panels,1),2,3) :: vp, vp_old
                real(kind=8), dimension(size(panels,1),3) :: vp3, vp3_old, vp4, vp4_old
                real(kind=8), dimension(size(panels,1),3) :: position_bl, position_br, position_bl_n, position_br_n

                integer :: i

                do i = 1,size(panels, 1)
                      vertex(i,:,:)=panels(i)%vertex
                      vp3(i,:)=panels(i)%vp(1,:)
                      vp3_old(i,:)=panels(i)%vp_old(1,:)
                      vp4(i,:)=panels(i)%vp(2,:)
                      vp4_old(i,:)=panels(i)%vp_old(2,:)

                      position_bl(i,:)=vertex(i, 4,:)
                      position_br(i,:)=vertex(i, 3,:)

                enddo

                      position_bl_n=position_bl+(3./2*vp4-1./2*vp4_old)*dt
                      position_br_n=position_br+(3./2*vp3-1./2*vp3_old)*dt

               do i = 1,size(panels, 1)
                      panels(i)%vertex(4,:)=position_bl_n(i, :)
                      panels(i)%vertex(3,:)=position_br_n(i, :)
                      panels(i)%vp_old(1,:)=vp3(i, :)
                      panels(i)%vp_old(2,:)=vp4(i, :)
               enddo


          end subroutine translations_wake_panels_v

subroutine radius_particle_s(particle_2mod, rotor,   time_particle)

  type(particle), intent(inout) ::  particle_2mod
  real(kind=8),   intent(in)    ::  time_particle
  integer, intent(in) :: rotor

  real(kind=8), dimension(3) :: alpha_p
  real(kind=8) :: new_radius

  ! Extracting the vector of alpha_p
  alpha_p = particle_2mod%alpha_p

  ! Computing new radius
  new_radius = ((3./4) * k_radius * time_particle * abs(sum(alpha_p)) + initial_radius(rotor)**3)**(1./3)

  ! Updating radius of the particle
  particle_2mod%radius = new_radius

endsubroutine radius_particle_s

subroutine radius_particle_v(particles_2mod, rotor,   time_particle)

  ! Assuming all the particles given have the same time (given by columns)
  type(particle), dimension(:), intent(inout)  ::  particles_2mod
  real(kind=8),                 intent(in)     ::  time_particle
  integer, intent(in) :: rotor

  real(kind=8), dimension(size(particles_2mod), 3) :: alpha_p
  real(kind=8), dimension(size(particles_2mod)) :: new_radius, time_particles, initial_radius_multi

  integer :: i

  ! Extracting the vector of alpha_p
  do i = 1, size(particles_2mod)
    alpha_p(i, :) = particles_2mod(i)%alpha_p
  enddo

  ! Creating a vector of times
  time_particles = spread(time_particle, dim = 1, ncopies = size(particles_2mod))

  ! Creating a vector of initial radius
  initial_radius_multi = spread(initial_radius(rotor), dim = 1, ncopies = size(particles_2mod))

  ! Computing new radius
  !!!!!!!!!!!!!!!!!!!!!! DA AGGIUSTARE CON ABS(SUM)) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  new_radius = ((3./4) * k_radius * time_particles * norm2(alpha_p, 2) + initial_radius_multi**3)**(1./3)

  ! Updating radius of the particle
  do i = 1, size(particles_2mod)
    particles_2mod(i)%radius = new_radius(i)
  enddo

endsubroutine radius_particle_v

subroutine radius_panel_s(wake_panel_2mod, rotor,  time_panel)

  type(wake_panel), intent(inout) ::  wake_panel_2mod
  real(kind=8),     intent(in)    ::  time_panel
  integer, intent(in) :: rotor

  real(kind=8) :: new_radius, gamma

  ! Extracting gamma from the panel
  gamma = wake_panel_2mod%gamma

  ! Computing new radius
  new_radius = sqrt(k_radius * time_panel * abs(gamma) / pi + initial_radius(rotor)**2)

  ! Updating radius of the particle
  wake_panel_2mod%radius = new_radius

endsubroutine radius_panel_s

subroutine radius_panel_v(wake_panels_2mod, rotor,   time_panel)

  type(wake_panel), dimension(:),  intent(inout) ::  wake_panels_2mod
  real(kind=8),                    intent(in)    ::  time_panel
  integer, intent(in) :: rotor

  real(kind=8), dimension(size(wake_panels_2mod)) :: new_radius, gamma,&
                                                     time_panels, initial_radius_multi

  ! Extracting gamma from the panel
  gamma = wake_panels_2mod%gamma

  ! Creating vector of time
  time_panels = spread(time_panel, dim=1, ncopies=size(wake_panels_2mod))

  ! Creating a vector of initial radius
  initial_radius_multi = spread(initial_radius(rotor), dim = 1, ncopies = size(wake_panels_2mod))

  ! Computing new radius
  new_radius = sqrt(k_radius * time_panels * abs(gamma) / pi + initial_radius_multi**2)

  ! Updating radius of the particle
  wake_panels_2mod%radius = new_radius

endsubroutine radius_panel_v

end module mod_advectVortices



! left-most wake segment
!
! || tl
! ||
! ||
! ||
! ||=====!=====
!  bl         br

!
! O <- new particle


! right-most wake segment
!
!         tr ||
!            ||
!            ||
! bl         ||
! =====!=====|| br
!

! new      -> O
! particle
