module mod_velocity

   use mod_geometry
   use mod_initialize
   use mod_aero_coeff
   use mod_types

   implicit none

contains


subroutine double_time (psi, beta_here, beta_dot, blade_panels,   wake_panels, particles, wake_segment)
! Versione scalare
! Ciclo interno, all'inizio di ogni passo temporale:
! calcolo la circolazione sui pannelli di pala considerando la scia congelata
! INPUT:    - psi           ->  vettore con la posizione di tutte le pale
!           - blade_panels  ->  struttura con tutti i pannelli di pala
!           - wake_panels   ->  struttura con tutti i pannelli di scia (opzionale)
!           - particles     ->  struttura con tutte le particelle (opzionale)
!           - wake_segment  ->  struttura con le gamma residue da particelle (opzionale)
! OUTPUT:   - blade_panels  ->  per ogni pannello viene cambiata solo la gamma

   implicit none

   real(kind=8),      dimension(:, :),       intent(in)           :: psi
   real(kind=8),      dimension(:, :),       intent(in)           :: beta_here
   real(kind=8),      dimension(:, :, :),    intent(in)           :: beta_dot
   type(blade_panel), dimension(:, :, :),    intent(inout)        :: blade_panels
   type(segment),     dimension(:, :, :),    intent(in), optional :: wake_segment
   type(wake_panel),  dimension(:, :, :, :), intent(in), optional :: wake_panels
   type(particle),    dimension(:, :, :),    intent(in), optional :: particles


   real(kind=8), dimension(n_rotors,&
                           n_blades, &
                           n_elements, &
                           n_rotors,&
                           n_blades, &
                           n_elements, 3)  :: influences

   real(kind=8), dimension(n_rotors, &
                           n_blades, &
                           n_elements, 3)  :: vi_wake

   real(kind=8), dimension(n_rotors, &
                           n_blades, &
                           n_elements)     :: gamma_old

   real(kind=8), dimension(n_rotors, &
                           n_blades, &
                           n_elements)     :: gamma_blade

   real(kind=8), dimension(n_rotors, &
                           n_blades, &
                           n_elements)     :: difference

   real(kind=8) :: alpha, mach, Cl, U, phi, q, Cd, S
   integer :: i, j, k, l, m, n, iter
   integer :: nmax = 300

   ! Initializing influences matrix
   influences = 0.

   ! calcolo inziale dei coefficienti di influenza dei pannelli e delle
   ! velocità indotte da scia e particelle (sono conglelati quindi posso
   ! accumulare tutti i contributi sul singolo punto di collocazione in un
   ! unico vettore!)
   do m = 1, n_rotors
     do j = 1, n_blades            ! | Ciclo sui punti di collocazione in cui
        do i = 1, n_elements       ! | verranno calcolate le influenze/velocità

           ! Ciclo su tutti i pannelli "sorgenti"
           do n = 1, n_rotors
             do k = 1, n_blades
                do l = 1, n_elements
                ! influences: struttura contenente l'influenza di TUTTI i pannelli
                ! su TUTTI i pannelli: (j, i) -> punto di collocazione
                !                      (k, l) -> pannello sorgente
                !                      ultimo indice -> componenti (x, y, z)
                   call influence_blade_panel(blade_panels(m, j, i)%collocation, &
                                        blade_panels(n, k,l),  influences(m, j, i, n, k, l, :))
                enddo
             enddo
           enddo

           if(present(wake_panels)) then
              if(.not.present(particles)) then
              ! vi_wake: vettore contenente la velocità indotta da ogni elemento
              ! di scia su TUTTI i pannelli di pala:
              ! (j, i) -> punto di collocazione
              ! ultimo indice -> componenti (x, y, z)
                 call vi_calc_wake(vi_wake(m, j, i, :), blade_panels(m, j, i)%collocation, &
                                   wake_panels)
              else
                 call vi_calc_wake(vi_wake(m, j, i, :), blade_panels(m, j, i)%collocation, &
                                   wake_panels, particles, wake_segment)
              end if
           end if

           gamma_old(m, j, i) = blade_panels(m, j, i)%gamma

        end do
     end do
   enddo

   ! Ciclo interno: 1) la circolazione in scia genera una v_i sulla pala
   !                2) dalla v_i calcolo l'angolo di incidenza della pala e
   !                   quindi il Cl
   !                3) Kutta-Joukowski: Cl -> circolazione che cambia la v_i
   !                   sulla pala => ** itero fino a convergenza sulla gamma **
   difference = toll_gamma + 1
   iter = 1
   do while (maxval(difference) > toll_gamma .and. iter < nmax)
      iter = iter + 1
      do m = 1, n_rotors
        do j = 1, n_blades
           do i = 1, n_elements
           ! Calcolo la gamma del singolo pannello: passo ad alpha_mach le influenze
           ! di TUTTI i pannelli (k=:, l=:) sul pannello (j, i) -> influence ha ora
           ! 3 indici; e di tutta la scia (vi_wake contiene già tutti i contributi)
           ! sul pannello (j, i) -> vi_wake ha 1 indice
              if(present(wake_panels)) then
                 call alpha_mach(alpha, mach, U,  y(m, i), psi(m, j), beta_here(m, j), beta_dot(m, j, 1), m, blade_panels, &
                              influences(m, j, i, :, :, :, :), vi_wake(m, j, i, :))
              else
                 call alpha_mach(alpha, mach, U,  y(m, i), psi(m, j), beta_here(m, j), beta_dot(m, j, 1), m, blade_panels, &
                              influences(m, j, i, :, :, :, :))
              endif

              if(dyn_stall) then
                 call dyn_stall_boeing(alpha, mach, blade_panels(m,j,i), m,   Cl)
              else
                 call cpcrcm(alpha, mach,  Cl)
              endif

              call gamma_calc(Cl, gamma_blade(m, j, i), U, m)


              ! calcolo errore (uso come criterio il massimo)
              difference(m, j, i) = abs(gamma_old(m, j, i) - gamma_blade(m, j, i))
              gamma_old(m, j, i) = gamma_blade(m, j, i)

           end do
        end do
      enddo

      !write(*,*) "maxval(difference) = ", maxval(difference)
      !write(*,*) "iter = ", iter

      do m = 1, n_rotors
        do j = 1, n_blades
          do i = 1, n_elements
            blade_panels(m, j, i)%gamma = gamma_blade(m, j, i)
          enddo
        enddo
      enddo
   end do !while
   if (iter == nmax) write(*,*) 'WARNING: Numero massimo di iterazioni raggiunte '&
                           'nella subroutine double_time'

   ! Calcolo prestazioni: dFz e dQ per ogni pannello
   do m = 1, n_rotors
     !!! valido solo per pannelli uniformi e rettangolari !!!
     S = c(m)*(y(m, 2) - y(m, 1)) ! Area del pannello
     do j = 1, n_blades
        do i = 1, n_elements
           if(present(wake_panels)) then
              call alpha_mach(alpha, mach, U,   y(m, i), psi(m, j), beta_here(m, j), beta_dot(m, j, 1), m, blade_panels, &
                           influences(m, j, i, :, :, :, :), vi_wake(m, j, i, :))
           else
              call alpha_mach(alpha, mach, U,  y(m, i), psi(m, j), beta_here(m, j), beta_dot(m, j, 1), m, blade_panels, &
                           influences(m, j, i, :, :, :, :))
           endif

           call cpcrcm(alpha, mach,  Cl, Cd) ! Calcola il Cd in ogni caso

           if(dyn_stall) & ! Ricalcola il Cl
                call dyn_stall_boeing(alpha, mach, blade_panels(m,j,i), m,   Cl)

           phi = theta(m) - theta_s(m) * sin(psi(m, j)) - alpha
           q = 0.5d0*rho*(U**2)
           blade_panels(m,j,i)%dFz = q*S*cos(beta_here(m, j))*(Cl*cos(phi) - Cd*sin(phi))
           blade_panels(m,j,i)%dQ = q*S*y(m, i)*(sin(phi)*Cl + cos(phi)*Cd)
           blade_panels(m,j,i)%alpha_old = alpha
        end do
     end do
   enddo

end subroutine double_time


subroutine alpha_mach(alpha, mach, U,    y_value, psi, beta_here, beta_dot, rotor, blade_panels, influences, vi_wake)

   implicit none

   real(kind=8),                             intent(in)           :: y_value
   real(kind=8),                             intent(in)           :: psi
   real(kind=8),                             intent(in)           :: beta_here
   real(kind=8),                             intent(in)           :: beta_dot
   real(kind=8),      dimension(:, :, :, :), intent(in)           :: influences
   type(blade_panel), dimension(:, :, :),    intent(in)           :: blade_panels
   real(kind=8),      dimension(:),          intent(in), optional :: vi_wake
   integer,                                  intent(in)           :: rotor

   real(kind=8), intent(out) :: alpha, mach, U

   real(kind=8), dimension(3) :: vi, vi_blade, Vinf_blade
   real(kind=8) :: Up, Ut
   integer :: j, k

   vi = 0.

   ! calcolo la velocita' indotta
   do k = 1, n_rotors
     do j= 1, n_blades
       vi = vi + matmul(transpose(influences(k, j, :, :)), blade_panels(k, j, :)%gamma)
     end do
   enddo

   if(present(vi_wake)) then
      vi = vi + vi_wake
   end if


   call xyz2alpha(psi, beta_here, vi, rotor,   vi_blade) ! Induced velocity in blade coordinate
   call xyz2alpha(psi, beta_here, Vinf, rotor,   Vinf_blade) ! Asimptotic velocity in blade coordinate

! note tutte le velocità indotte possiamo calcolare Up,Ut, e quindi alpha e Mach
   Up = - vi_blade(3) - Vinf_blade(3) + (y_value - e(rotor)) * beta_dot
   Ut = vi_blade(1) + Vinf_blade(1) + sign(1.d0, omega(rotor)) * omega(rotor) * y_value

   alpha = theta(rotor) - theta_s(rotor) * sin(psi) - atan(Up/Ut)
   U = sqrt(Up**2 + Ut**2)
   mach = U / sound_speed

end subroutine alpha_mach


subroutine influence_blade_panel(point, panel,     influence, j, k)
! influence_blade_panel: Subroutine used to calc the influence of a blade panel vortex on a point

! Input:  point: vector containing the 3 coords of the induced point
!         panel:  structure containing the information about the panel vortex considered

! Output: influence: Vector containing the 3 components of the influence of a blade panel vortex on the point

  real(kind=8), dimension(:), intent(in) :: point
  type(blade_panel), intent(in) :: panel
  integer, intent(in), optional :: j, k

  real(kind=8), dimension(3), intent(inout) :: influence

  real(kind=8), dimension(2, 3) :: vertex
  real(kind=8), dimension(3) :: vertex_left, vertex_right, influence_seg
  type(segment) :: segmentino
  real(kind=8) :: radius
  integer :: i, i_l, i_r

  influence = 0  ! Initialization of the influence of the whole panel

  radius = panel%radius

  !write(*,*) "Influenza del pannello di pala"

  ! Cycle on all of the segment which belongs to the panel
  do i = 1, 4

    i_l = i      ! selecting the left corner of the segment
    i_r = i + 1  ! selecting the right corner of the segment
    if(i_r > 4)  i_r = 1;  ! back to vertex 1 if i = 4

    vertex_left   =  panel%vertex(i_l, :)  ! Extracting the coordinates of the left corner
    vertex_right  =  panel%vertex(i_r, :)  ! Extracting the coordinates of the right corner

    vertex(1, :) = vertex_left;  vertex(2, :) = vertex_right
    radius = panel%radius  ! Extracting information about the radius of the vortex

    call wake_segment_constructor(panel%gamma, vertex, radius,   segmentino)

    if(present(k)) then
      call influence_segment(point, segmentino,     influence_seg, j, k)
    else
      call influence_segment(point, segmentino,     influence_seg)
    endif


    influence = influence + influence_seg  ! Adding the influence of each vortex of the panel

  enddo

  ! Extracting informations about the vortex
  !gamma = panel%gamma

  !vi = influence * gamma  ! Computing the vector of the induced velocity

endsubroutine  influence_blade_panel

subroutine influence_wake_panel(point, panel,     influence,  k, l)
! influence_wake_panel: Subroutine used to calc the influence of the wake panel vortex on a point

! Input:  point: vector containing the 3 coords of the induced point
!         panel:  structure containing the information about the panel vortex considered

! Output: influence: Vector containing the 3 components of the influence of the wake panel vortex on the point


  real(kind=8), dimension(:), intent(in) :: point
  type(wake_panel), intent(in) :: panel
  integer, intent(in), optional :: k, l

  real(kind=8), dimension(3), intent(inout) :: influence

  real(kind=8), dimension(2, 3) :: vertex
  real(kind=8), dimension(3) :: vertex_left, vertex_right, influence_seg
  type(segment) :: segmentino
  real(kind=8) :: radius
  integer :: i, i_l, i_r

  influence = 0  ! Initialization of the influence of the whole panel
  !write(*,*) " Influenza del pannello di scia"

  ! Cycle on all of the segment which belongs to the panel
  do i = 1, 4

    i_l = i      ! selecting the left corner of the segment
    i_r = i + 1  ! selecting the right corner of the segment
    if(i_r > 4)  i_r = 1;  ! back to vertex 1 if i = 4

    vertex_left   =  panel%vertex(i_l, :)  ! Extracting the coordinates of the left corner
    vertex_right  =  panel%vertex(i_r, :)  ! Extracting the coordinates of the right corner

    vertex(1, :) = vertex_left;  vertex(2, :) = vertex_right
    radius = panel%radius  ! Extracting information about the radius of the vortex

    call wake_segment_constructor(panel%gamma, vertex, radius,   segmentino)

    call influence_segment(point, segmentino,     influence_seg)

    influence = influence + influence_seg  ! Adding the influence of each vortex of the panel

  enddo
  !write(*,*) "influence = ", influence

  ! Extracting informations about the vortex
  !gamma = panel%gamma

  !vi = influence * gamma  ! Computing the vector of the induced velocity

endsubroutine  influence_wake_panel

subroutine influence_segment(point, segmentino,     influence, k, l)
! influence_wake_panel: Subroutine used to calc the influence of the wake panel vortex on a point

! Input:  point: vector containing the 3 coords of the induced point
!         panel:  structure containing the information about the panel vortex considered

! Output: influence: Vector containing the 3 components of the influence of the wake panel vortex on the point

  type(segment), intent(in) :: segmentino
  real(kind=8), dimension(:), intent(in) :: point
  integer, intent(in), optional :: k, l

  real(kind=8), dimension(3), intent(inout) :: influence

  real(kind=8), dimension(3) :: vertex_left, vertex_right,&
                                r1, r2,&
                                l12, e_vec, cross
  real(kind=8) :: h, radius, costeta1, costeta2, factor

  influence = 0  ! Initialization of the influence of the segment

  vertex_left   =  segmentino%vertex(1, :)  ! Extracting the coordinates of the left corner
  vertex_right  =  segmentino%vertex(2, :)  ! Extracting the coordinates of the right corner

  l12 = vertex_left - vertex_right   ! Computing the components of the congiungente leftcorner ---> rightcorner
  r2 = point - vertex_right          ! Computing the components of the congiungente rightcorner ---> point
  r1 = point - vertex_left           ! Computing the components of the congiungente leftcorner ---> point

  costeta1=DOT_PRODUCT(l12,r1)/(norm2(l12)*norm2(r1))
  costeta2=DOT_PRODUCT(l12,r2)/(norm2(l12)*norm2(r2))

  h = norm2(r1) * sqrt(1 - costeta1**2)   ! Computing the distance point-vortex_segment

  !if(h > toll_influ) then   ! Check if the point is one of the segment.
                    ! In that case there's no contribution from them

  radius = segmentino%radius  ! Extracting information about the radius of the vortex
  factor=costeta1-costeta2

  call cross_product(l12, r1,   cross)  ! Computing cross product

  if(norm2(cross) > toll_influ) then

    e_vec = cross / norm2(cross)  ! Computing the unit vector indicating the direction of the induced velocity

    influence = influence + e_vec * h * factor/ (4 * pi * (h**2 + radius**2))  ! Adding the influence of each vortex of the panel
  else
    !write(*,*) "è successo, cross = ", cross
    !write(*,*) "point = ", point
    !write(*,*) "r1 = ", r1
    !write(*,*) "l12 = ", l12
    !write(*,*) "vertex_left = ", vertex_left
    !write(*,*) "vertex_right = ", vertex_right
    !write(*,*) " "
  endif

endsubroutine  influence_segment

subroutine influence_particle(point, vortex_particle,     influence)
! influence_particle: Subroutine used to calc the influence of a particle on a point

! Input:  point: vector containing the 3 coords of the induced point
!         vortex_particle:  structure containing the information about the particle considered

! Output: influence: Vector containing the 3 components of the influence of the particle on the point

  real(kind=8), dimension(:), intent(in) :: point
  type(particle), intent(in) :: vortex_particle

  real(kind=8), dimension(3), intent(inout) :: influence

  real(kind=8), dimension(3) :: center_coord, distance
  real(kind=8) :: radius, rho_particle, troublesome_component

  !vi = 0  ! Initialization of the influence of the particle

  center_coord =  vortex_particle%position    ! Extracting the coordinates of the center of the particle
  !alpha_p  =  vortex_particle%alpha_p       ! Extracting informations about the vortex

  !if(point(1) /= center(1)) then   ! Check if the point is the center of the particle.
                                   ! In that case there's no contribution from it

    distance = point - center_coord       ! Computing the components of the congiungente center ---> point

    if(norm2(distance) < toll_influ) then   ! if r is too small it will converge to 1

      troublesome_component = 1
      !distance = [0, 0, 0]     ! Make the particle not to contribute to the vi since the point
                                ! is near the center

    else

      radius = vortex_particle%radius  ! Extracting information about the radius of the particle
      rho_particle = norm2(distance) / radius ! Computing the exponent rho

      troublesome_component = (1 - exp(-rho_particle**3)) / (norm2(distance)**3)

    endif

    !vi = - 1/ (4 * pi) * troublesome_component * cross_product(r, alpha_p)

    influence = - 1/ (4 * pi) * troublesome_component * distance

endsubroutine influence_particle

subroutine vi_calc(vi,   point,  blade_panels, wake_panels, particles, wake_segment)  ! Make things such that wake panels and particles can be not given
! vi_calc: subroutine that computes the induced velocity in a point

! Input:  point: vector of coordinates of the point where we want to calculate the induced velocity at
!         blade_panels: structure containing all of the panels on the blade
!         wake_panels: (optional) structure containing all of the wake panels that are present
!         particles: (optional) structure containing all of the particles that are present
!         wake_segment: (optional) structure containing the residual segment on the wake

! Output: vi: vector that refers to the 3 components of the total induced velocity calculated on the input point

  type(wake_panel),   dimension(:, :, :, :), intent(in), optional    :: wake_panels
  type(segment),      dimension(:, :, :),    intent(in), optional    :: wake_segment
  type(particle),     dimension(:, :, :),    intent(in), optional    :: particles
  type(blade_panel),  dimension(:, :, :),    intent(in)              :: blade_panels
  real(kind=8),       dimension(:),          intent(in)              :: point

  real(kind=8), dimension(3),                intent(inout)           :: vi

  real(kind=8), dimension(3) :: vi_toadd, influence, vi_blade, vi_wake,&
                                vi_wake_segment, vi_particles
  integer :: i, j, k, m

  ! Initialization of the induced velocity
  vi = 0;   vi_toadd = 0;   vi_blade = 0;   vi_wake = 0
  vi_wake_segment = 0;   vi_particles = 0; influence = 0.

  do m = 1, n_rotors
    do i = 1, n_blades ! Moving on the n_blades
      do j = 1, n_elements  ! Moving on the n_elements

        ! Computing all of the influences induced by all of the elements of the blades
        call influence_blade_panel(point, blade_panels(m, i, j),  influence)

        ! Computing induced velocity
        vi_toadd = blade_panels(m, i, j)%gamma * influence

        ! Update of the velocity induced by the blade
        vi_blade = vi_blade + vi_toadd

        ! Reset vi_toadd to 0
        vi_toadd = 0

        influence = 0.

        if(present(wake_panels)) then  ! Check if there are panels in the wake

          do k = 1, size(wake_panels, 4)

            ! Computing all of the influences induced by all of the panels in the wake
            call influence_wake_panel(point, wake_panels(m, i, j, k), influence)

            ! Computing induced velocity
            vi_toadd = wake_panels(m, i, j, k)%gamma * influence

            ! Update of the velocity induced by the wake
            vi_wake = vi_wake + vi_toadd

            ! Reset vi_toadd to 0
            vi_toadd = 0

            influence = 0.

          enddo

        endif

      enddo
    enddo
  enddo


  if(present(particles)) then ! Check if there are particles in the wake
    do m = 1, n_rotors
      do i = 1, size(particles, 2)
        do j = 1, size(particles, 3)

          ! Computing all of the influences induced by all of the particles in the wake
          call influence_particle(point, particles(m, i, j), influence)

          ! Computing induced velocity
          call  cross_product(influence, particles(m, i, j)%alpha_p,   vi_toadd)

          ! Update of the velocity induced by the particles in the wake
          vi_particles = vi_particles + vi_toadd

          ! Reset vi_toadd to 0
          vi_toadd = 0

          influence = 0.

        enddo
      enddo
    enddo

    ! If particles are present then also the wake segment is present
    do m = 1, n_rotors
      do i = 1, size(wake_segment, 2)
        do j = 1, size(wake_segment, 3)

          ! Computing the influence of the wake segment
          call influence_segment(point, wake_segment(m, i, j), influence)

          ! Computing induced velocity
          vi_toadd = wake_segment(m, i, j)%gamma * influence

          ! Update of the velocity induced by the wake
          vi_wake_segment = vi_wake_segment + vi_toadd

          ! Reset vi_toadd to 0
          vi_toadd = 0

          influence = 0.

        enddo
      enddo
    enddo

  endif

  !write(*,*) "vi_blade = ", vi_blade
  !write(*,*) "vi_wake = ", vi_wake

  ! Update of the total induced velocity
  vi = vi_blade + vi_wake + vi_wake_segment + vi_particles

endsubroutine vi_calc

subroutine vi_calc_wake(vi,   point, wake_panels, particles, wake_segment)  ! Make things such that wake panels and particles can be not given

  type(wake_panel),  dimension(:, :, :, :), intent(in)            :: wake_panels
  type(segment),     dimension(:, :, :),    intent(in), optional  :: wake_segment
  type(particle),    dimension(:, :, :),    intent(in), optional  :: particles
  real(kind=8),      dimension(:),          intent(in)            :: point

  real(kind=8),      dimension(3),          intent(inout)         :: vi

  real(kind=8), dimension(3) :: vi_toadd, influence, vi_wake, vi_wake_segment, vi_particles
  integer :: i, j, k, m

  ! Initialization of the induced velocity
  vi = 0;   vi_toadd = 0;   vi_wake = 0;  vi_wake_segment = 0;   vi_particles = 0;

  do m = 1, n_rotors
    do i = 1, n_blades ! Moving on the n_blades
      do j = 1, n_elements  ! Moving on the n_elements

        ! Computing all of the influences induced by all of the panels in the wake
        do k = 1, size(wake_panels, 4)

           call influence_wake_panel(point, wake_panels(m, i, j, k), influence)

           ! Computing induced velocity
           vi_toadd = wake_panels(m, i, j, k)%gamma * influence

           ! Update of the velocity induced by the wake
           vi_wake = vi_wake + vi_toadd

           ! Reset vi_toadd to 0
           vi_toadd = 0

        enddo
      enddo
    enddo
  enddo


  if(present(particles)) then ! Check if there are particles in the wake

    ! Computing all of the influences induced by all of the particles in the wake
    do m = 1, n_rotors
      do i = 1, size(particles, 2)
         do j = 1, size(particles, 3)

            call influence_particle(point, particles(m, i, j), influence)

            ! Computing induced velocity
            call  cross_product(influence, particles(m, i, j)%alpha_p,   vi_toadd)

            ! Update of the velocity induced by the particles in the wake
            vi_particles = vi_particles + vi_toadd

            ! Reset vi_toadd to 0
            vi_toadd = 0

         enddo
      enddo
    enddo

    ! If particles are present then also the wake segment is present
    do m = 1, n_rotors
      do i = 1, size(wake_segment, 2)
        do j = 1, size(wake_segment, 3)

          ! Computing the influence of the wake segment
          call influence_segment(point, wake_segment(m, i, j), influence)

          ! Computing induced velocity
          vi_toadd = wake_segment(m, i, j)%gamma * influence

          ! Update of the velocity induced by the wake
          vi_wake_segment = vi_wake_segment + vi_toadd

          ! Reset vi_toadd to 0
          vi_toadd = 0

        enddo
      enddo
    enddo

  endif

  ! Update of the total induced velocity
  vi = vi_wake + vi_particles + vi_wake_segment

endsubroutine vi_calc_wake

subroutine gamma_calc(Cl, gamma, U, rotor)

   implicit none

   integer, intent(in) :: rotor
   real(kind=8), intent(in)    :: Cl, U
   real(kind=8), intent(inout) :: gamma

   gamma = - 0.5 * U * c(rotor) * Cl

end subroutine gamma_calc

subroutine dyn_stall_boeing(alpha, mach, this_blade_panel, rotor_index,  Cl)

   real(kind=8), intent(in) :: alpha, mach
   type(blade_panel), intent(in) :: this_blade_panel
   integer :: rotor_index

   real(kind=8), intent(out) :: Cl

   real(kind=8) :: alpha_dot, MN1, MN0, max_gamma2, thickness_ratio, gamma2, gamma1, &
        break_value, K1, alpha0, alpha_ref, valore, one

   one = 1.

   !!! Warning vale solo per il NACA0012 !!!
   thickness_ratio = 0.12
   alpha0 = 0.
   !!! end warning !!!

   alpha_dot = (alpha - this_blade_panel%alpha_old)/dt
call cpcrcm(alpha, mach,  Cl)
   if (alpha_dot == 0.) then ! Non serve modello dinamico
      call cpcrcm(alpha, mach,  Cl)

   else ! serve modello dinamico
      MN1 = 0.4 + 5*(0.06 - thickness_ratio)
      MN0 = 0.9 + 2.5*(0.06 - thickness_ratio)
      max_gamma2 = 1.4 - 6*(0.06 - thickness_ratio)
      break_value = 0.6 + 1.5*(0.06 - thickness_ratio)

      if (mach < MN1) then
         gamma2 = max_gamma2
      elseif (mach > MN0) then
         gamma2 = 0
      else
         gamma2 = max_gamma2*(mach - MN0) / (MN1 - MN0)
      endif

      gamma1 = gamma2/2.

      K1 = 1
      if (sign(one, alpha_dot) < 0) K1 = 0.5

      valore = sqrt(abs(c(rotor_index)*alpha_dot/(2*mach*sound_speed)))

      if (valore < break_value) then
         alpha_ref = alpha - K1*gamma1*valore*sign(one, alpha_dot)
      else
         alpha_ref = alpha - K1*(gamma1*break_value + gamma2*(valore - break_value))*sign(one, alpha_dot)
      endif

      call cpcrcm(alpha_ref, mach,  Cl)

      Cl = Cl*alpha/(alpha_ref - alpha0)

   endif


end subroutine dyn_stall_boeing

end module mod_velocity
