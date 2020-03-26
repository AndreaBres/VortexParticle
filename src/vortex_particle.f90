program vortex_particle

use mod_further_beyond
use mod_types
use mod_velocity
use mod_geometry
use mod_aero_coeff
use mod_advectVortices
use mod_flapping
use mod_output
use omp_lib

implicit none


! Varibales definition

type(wake_panel),  dimension(:, :, :, :),    allocatable :: wake_panels
type(segment),     dimension(:, :,:),        allocatable :: wake_segment ! Residual gamma from panel2particles transofrmations
type(particle),    dimension(:, :, :),       allocatable :: particles, particles_copy
type(particle),    dimension(:),             allocatable :: new_particles
type(blade_panel), dimension(:, :, :),       allocatable :: blade_panels
real(kind=8)                                             :: time   ! Advancing time variable
real(kind=8),      dimension(:, :),          allocatable :: psi    ! Blade rotation
real(kind=8),      dimension(:, :, :),       allocatable :: beta_dot
real(kind=8),      dimension(:, :, :),       allocatable :: beta_ddot
real(kind=8),      dimension(:, :),          allocatable :: Ma    ! Blade rotation
real(kind=8),      dimension(:,:),           allocatable :: ct    ! Thrust coefficient
real(kind=8),      dimension(:,:),           allocatable :: cp    ! Power coefficient
real(kind=8),      dimension(:,:,:),         allocatable :: bp_vertex ! matrice dei vertici dei pannelli
real(kind=8),      dimension(:,:,:,:),       allocatable :: blade_extremes ! Matrice dei vertici delle pale
real(kind=8),      dimension(:),             allocatable :: comp_time_vi, num_particles,&
                                                            comp_time_gamma, integral, torque
real(kind=8),      dimension(:, :, :, :),    allocatable :: vi_on_blades
real(kind=8),      dimension(:, :),          allocatable :: Re


integer :: n_particles,&  ! Maximum number of particles
           max_temp_iter, &  ! Maximum number of temporal iterations
           j, i, k, l, &
           n_particles_toadd, iter_temp = 1

real(kind=8), dimension(3)    :: alpha_p, position, collocation, vp3, vp4, vi, &
                                 new_vertex1,new_vertex2, velocity_vertex3, velocity_vertex4
real(kind=8), dimension(2, 3) :: vp, vp_old
real(kind=8), dimension(:,:, :,:,:), allocatable :: new_vertex12, velocity_vertex34
real(kind=8)                  :: radius, up, ut, gamma, dFz, dQ, alpha
real(kind=8), dimension(4,3)  :: vertex
real(kind=8), dimension(2,3)  :: vertex_segment
integer :: intersection_particles
logical :: check_intersection
!real(kind=8) :: initial_radius

! Computational time variables
real(kind = 8) :: start, finish, start_int, finish_int, start_iter, finish_iter
character (len = 50) :: filename, format_string

integer :: f, u

real(kind=8) :: distance_max = 0., distance

! Directories
character(len = 10) :: dir
logical :: ex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialization of variables from input file
call CPU_TIME(start)

call init()
call init2()
if(n_boundaries > 0)  call init3();

!call OMP_set_nested(.true.)

! check for directories
dir = './solution'
INQUIRE (FILE=dir, EXIST=ex)
if(ex) then
else
  call system('mkdir -p ./solution' )
  call system('mkdir -p ./solution/blade_coord' )
  call system('mkdir -p ./solution/blade_coord_Fz' )
  call system('mkdir -p ./solution/blade_Fz' )
  call system('mkdir -p ./solution/wake_panel_coord' )
  call system('mkdir -p ./solution/particle_coord' )
  call system('mkdir -p ./solution/particle_coords' )
  call system('mkdir -p ./solution/particle_intersect' )
  call system('mkdir -p ./solution/vi_on_blade' )
endif

call CPU_TIME(finish)
write(*,*) " "
write(*,'(A, F6.4, A)') "Elapsed time for Data acquisition = ", finish - start, " s"

!!! Variable allocation !!!
call initial_allocation()

! Pre process: blade geometry
do i = 1, n_rotors
  call collocation_points(i, y(i, :))
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Inizializzazione segmenti sfigati di scia
!gamma_segment = 0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initial blade position
do l = 1, n_rotors
  do j = 1, n_blades
     psi(l, j) = (j - 1)*2*pi/n_blades
  enddo
enddo

! Initialization of the radius for both panels and particles
!initial_radius = 0.05 * c

! Blade panel inizialization
gamma = 0.
dFz = 0.
dQ = 0.
alpha = 0.
do k = 1, n_rotors
  call blade_matrix(k,    bp_vertex)
  do j = 1, n_blades
     do i = 1, n_elements
        vertex = blade_panel_vertex(bp_vertex(1,i,:), bp_vertex(1,i+1,:), &
                    bp_vertex(2,i+1,:), bp_vertex(2,i,:), psi(k, j), beta(k, j), k) + &
                    spread(center((1 + (k-1)*3 ): (k * 3)), dim=1, ncopies=4)
        collocation = global_blade_vector(psi(k, j), beta(k, j), [0.d0, y(k, i)], k) + center((1 + (k-1)*3 ): (k * 3))
        call blade_panel_constructor(gamma, initial_radius(k), collocation, vertex, dFz, dQ, alpha,  blade_panels(k, j, i))
     end do
  end do
enddo

! Wake panel inizialization
gamma = 0.
vertex = 0.
vp = 0.
vp_old = 0.
do l = 1, n_rotors
  radius = initial_radius(l)
  do j = 1, n_blades
     do i = 1, n_elements
        do k = 1, n_panels
           call wake_panel_constructor(gamma, vertex, vp, vp_old, radius,   wake_panels(l, j, i, k))
        enddo
     enddo
  enddo
enddo

! Wake segment inizialization
gamma = 0.
vertex_segment = 0.
do l = 1, n_rotors
  radius = initial_radius(l)
  do j = 1, n_blades
     do i = 1, n_elements
        call wake_segment_constructor(gamma, vertex_segment, radius,   wake_segment(l, j, i))
     enddo
  enddo
enddo

! beta Initialization
beta_dot = 0.
beta_ddot = 0.
integral = 0.
torque = 0.
Ma = 0.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! STARTING TEMPORAL ITERATION !

call CPU_TIME(start_int)

OPEN ( UNIT = 20, FILE = "solution/flapping.dat", ACTION= 'READWRITE', STATUS = 'REPLACE')

do while(time < tf)

  !if (time > 0.25)   Vinf = [0., 0., 1.]

   call CPU_TIME(start_iter)

   !call CPU_TIME(start)

   if(allocated(particles)) then
     num_particles(iter_temp) = size(particles, 1) * size(particles, 2) * size(particles, 3)
   else
     num_particles(iter_temp) = 0
   endif


  if(allocated(particles)) then
  ! Cycling on all the particles
  do i = 1, size(particles, 1)
    do j = 1, size(particles, 2)
      do k = 1, size(particles, 3)
        particles(i, j, k)%radius_old = particles(i, j, k)%radius
           ! if(.not.particles(i, j, k)%on_ground) call radius_particle(particles(i, j, k), i,  k * dt)
          call radius_particle(particles(i, j, k), i,  k * dt)
      enddo
    enddo
  enddo
  !!! BISOGNA CAMBIARLO ANCHE PER IL SEGMENTINO??? !!!
  endif

  if(iter_temp <= n_panels + 1 .and. iter_temp > 1) then
    ! Cycling on all the wake panels
    do l = 1, n_rotors
      do i = 1, n_blades
        do j = 1, n_elements
          do k = 1, iter_temp-1
            ! Passing columns of panels: this way every panel has the same time
            call radius_panel(wake_panels(l, i, j, k), l,  k * dt)
          enddo
        enddo
      enddo
    enddo
  elseif(iter_temp > n_panels + 1) then
      ! Cycling on all the wake panels
      do l = 1, n_rotors
        do i = 1, n_blades
          do j = 1, n_elements
            do k = 1, n_panels
              ! Passing columns of panels: this way every panel has the same time
              call radius_panel(wake_panels(l, i, j, k), l,  k * dt)
            enddo
          enddo
        enddo
      enddo
  endif

   !call CPU_TIME(finish)
   !write(*,*) " "
   !write(*,'(A, F6.4, A)') "Elapsed time for Enlarging of vortices = ", finish - start, " s"

   call CPU_TIME(start)

   !write(*,*) " Sto calcolando gamma blade"

!  1) Inner time loop: calcolo la circolazione sugli elementi di pala
   if (iter_temp == 1) then ! solo pala
      call double_time (psi, beta, beta_dot, blade_panels)
   elseif (iter_temp <= n_panels + 1 ) then ! solo pala e (parte) della scia
      call double_time (psi, beta, beta_dot, blade_panels, wake_panels(:, :, :, 1:iter_temp - 1))
   else ! ho tutte le strutture!
      call double_time (psi, beta, beta_dot, blade_panels, wake_panels(:, :, :, 1:n_panels), particles, wake_segment)
   end if
   !write(*,*) " Finito gamma blade"
   !write(*,*) " "
   !write(*,*) " "

   call CPU_TIME(finish)
   write(*,*) " "
   write(*,'(A, F8.4, A)') "Elapsed time for Computation of gamma blade = ", finish - start, " s"

   comp_time_gamma(iter_temp) = finish - start

!  2) Calcolo il campo di moto (nei vertici dei pannelli di scia e nei centri
!     delle particelle

call CPU_TIME(start)

! Computing induced velocities on all the vertexes of the wake panels


  ! Cycling on all the blade panels
  if(iter_temp == 1) then
    do l = 1, n_rotors
     do j = 1, n_blades
        do i = 1, n_elements

            !!! TOP-RIGHT corner of the blade panel !!!
            ! Computing induced velocity
            call vi_calc(vi,   blade_panels(l, j, i)%vertex(3,:),  blade_panels)
            ! Storing induced velocity
            vi_on_blades(l, j, i, :) = vi
            ! Computing total velocity
            velocity_vertex3 = vi + Vinf
            ! Computing position of the vertex after translation
            new_vertex1 = blade_panels(l, j,i)%vertex(3,:) + velocity_vertex3*dt

            !vi_blade(index_vi_blade) = norm2(vi)

            !!! BOTTOM-RIGHT corner of the blade panel !!!
            ! Computing induced velocity
            call vi_calc(vi,   blade_panels(l, j,i)%vertex(4,:),  blade_panels)
            ! Storing induced velocity
            if (i == n_elements)  vi_on_blades(l, j, i+1, :) = vi;
            ! Computing total velocity
            velocity_vertex4 = vi + Vinf
            ! Computing position of the vertex after translation
            new_vertex2 = blade_panels(l, j,i)%vertex(4,:) + velocity_vertex4*dt

            !vi_blade(index_vi_blade+1) = norm2(vi)

            ! Update the big structure containing all the vertexes of blade panels translated
            new_vertex12(1,:, l,j,i) = new_vertex1
            new_vertex12(2,:, l,j,i) = new_vertex2

            ! Update the big structure containing all the velocities of vertexes of blade panels
            velocity_vertex34(1,:, l,j,i) = velocity_vertex3
            velocity_vertex34(2,:, l,j,i) = velocity_vertex4

            !index_vi_blade = index_vi_blade + 2
          enddo
        enddo
      enddo

  elseif (iter_temp <= n_panels + 1 .and. iter_temp > 1) then  ! wake_panels are present
    do l = 1, n_rotors
      do j = 1, n_blades
         do i = 1, n_elements
            !!! TOP-RIGHT corner of the blade panel !!!
            ! Computing induced velocity
            call vi_calc(vi,   blade_panels(l,j,i)%vertex(3,:),  blade_panels, wake_panels(:, :, :, 1:iter_temp-1))
            ! Storing induced velocity
            vi_on_blades(l, j, i, :) = vi
            ! Computing total velocity
            velocity_vertex3 = vi + Vinf
            ! Computing position of the vertex after translation
            new_vertex1 = blade_panels(l,j,i)%vertex(3,:) + velocity_vertex3*dt

            !vi_blade(index_vi_blade) = norm2(vi)

            !!! BOTTOM-RIGHT corner of the blade panel !!!
            ! Computing induced velocity
            call vi_calc(vi,   blade_panels(l,j,i)%vertex(4,:),  blade_panels, wake_panels(:, :, :, 1:iter_temp-1))
            ! Storing induced velocity
            if (i == n_elements)  vi_on_blades(l, j, i+1, :) = vi;
            ! Computing total velocity
            velocity_vertex4 = vi + Vinf
            ! Computing position of the vertex after translation
            new_vertex2 = blade_panels(l,j,i)%vertex(4,:) + velocity_vertex4*dt

            !vi_blade(index_vi_blade+1) = norm2(vi)

            ! Update the big structure containing all the vertexes of blade panels translated
            new_vertex12(1,:, l,j,i) = new_vertex1
            new_vertex12(2,:, l,j,i) = new_vertex2

            ! Update the big structure containing all the velocities of vertexes of blade panels
            velocity_vertex34(1,:, l,j,i) = velocity_vertex3
            velocity_vertex34(2,:, l,j,i) = velocity_vertex4

            !index_vi_blade = index_vi_blade + 2

          enddo
        enddo
      enddo
    else
      do l = 1, n_rotors
        do j = 1, n_blades
           do i = 1, n_elements
            !!! TOP-RIGHT corner of the blade panel !!!
            ! Computing induced velocity
            call vi_calc(vi,   blade_panels(l,j,i)%vertex(3,:),  blade_panels, &
                         wake_panels(:, :, :, 1:n_panels), particles, wake_segment)
            ! Storing induced velocity
            vi_on_blades(l, j, i, :) = vi
            ! Computing total velocity
            velocity_vertex3 = vi + Vinf

            ! Computing position of the vertex after translation
            new_vertex1 = blade_panels(l,j,i)%vertex(3,:) + velocity_vertex3*dt

            !vi_blade(index_vi_blade) = norm2(vi)

            !!! BOTTOM-RIGHT corner of the blade panel !!!
            ! Computing induced velocity
            call vi_calc(vi,   blade_panels(l,j,i)%vertex(4,:),  blade_panels, &
                         wake_panels(:, :, :, 1:n_panels), particles, wake_segment)
            ! Storing induced velocity
            if (i == n_elements)  vi_on_blades(l, j, i+1, :) = vi;
            ! Computing total velocity
            velocity_vertex4 = vi + Vinf
            ! Computing position of the vertex after translation
            new_vertex2 = blade_panels(l,j,i)%vertex(4,:) + velocity_vertex4*dt

            !vi_blade(index_vi_blade+1) = norm2(vi)

            ! Update the big structure containing all the vertexes of blade panels translated
            new_vertex12(1,:, l,j,i) = new_vertex1
            new_vertex12(2,:, l,j,i) = new_vertex2

            ! Update the big structure containing all the velocities of vertexes of blade panels
            velocity_vertex34(1,:, l,j,i) = velocity_vertex3
            velocity_vertex34(2,:, l,j,i) = velocity_vertex4

            !index_vi_blade = index_vi_blade + 2

          enddo
        enddo
      enddo
    endif

  if (iter_temp <= n_panels + 1 .and. iter_temp > 1) then  ! wake_panels are present


      ! Cycling on all wake_panels present
      do l = 1, n_rotors
        do j = 1, n_blades
           do i = 1, n_elements
              do k = 1,iter_temp-1

                 !!! TOP-RIGHT corner of the blade panel !!!
                 ! Computing induced velocity
                 call vi_calc(vi,   wake_panels(l,j,i,k)%vertex(3,:),  blade_panels, wake_panels(:,:,:,1:iter_temp-1))

                 ! Updating vp
                 wake_panels(l,j,i,k)%vp(1,:) = vi + Vinf

                 !vi_wake(index_vi_wake) = norm2(vi)

                 !!! BOTTOM-RIGHT corner of the blade panel !!!
                 ! Computing induced velocity
                 call vi_calc(vi,   wake_panels(l,j,i,k)%vertex(4,:),  blade_panels, wake_panels(:,:,:,1:iter_temp-1))
                 ! Updating vp
                 wake_panels(l,j,i,k)%vp(2,:) = vi + Vinf

                 !vi_wake(index_vi_wake+1) = norm2(vi)

                 !index_vi_wake = index_vi_wake + 2

              enddo
            enddo
          enddo
        enddo

    elseif(iter_temp > n_panels + 1) then ! Every structure is present

      ! Cycling on all the wake_panels
      do l = 1, n_rotors
        do j = 1, n_blades
           do i = 1, n_elements
              do k = 1, n_panels

                 !!! TOP-RIGHT corner of the blade panel !!!

                 ! Computing induced velocity
                 call vi_calc(vi,   wake_panels(l,j,i,k)%vertex(3,:),  blade_panels,&
                              wake_panels(:, :, :, 1:n_panels), particles, wake_segment)
                 !write(*,*) "vi = ", vi

                 ! Updating vp
                 wake_panels(l,j,i,k)%vp(1,:) = vi + Vinf

                 !vi_wake(index_vi_wake) = norm2(vi)

                 !!! BOTTOM-RIGHT corner of the blade panel !!!
                 ! Computing induced velocity
                 call vi_calc(vi,   wake_panels(l,j,i,k)%vertex(4,:),  blade_panels,&
                              wake_panels(:, :, :, 1:n_panels), particles, wake_segment)
                 ! Updating vp
                 wake_panels(l,j,i,k)%vp(2,:) = vi + Vinf

                 !vi_wake(index_vi_wake+1) = norm2(vi)

                 !index_vi_wake = index_vi_wake + 2

              enddo
           enddo
        enddo
      enddo
   endif


   ! Checking if particles are present
   if(allocated(particles))then
     ! Cycling on all the particles
     !$OMP PARALLEL
     do l = 1, size(particles, 1)

       !$OMP DO PRIVATE (vi)
       do j = 1, size(particles, 2)

         do i = 1, size(particles, 3)

           ! Computing induced velocity on the particle
           call vi_calc(vi,   particles(l, j, i)%position,  blade_panels, wake_panels(:, :, :, 1:n_panels), particles, wake_segment)
           ! Update of the particle velocity
           particles(l, j, i)%vp = vi + Vinf

           !vi_particle(index_vi_particle+1) = norm2(vi)

           !index_vi_particle = index_vi_particle + 1

         enddo


       enddo
       !$OMP END DO
     enddo

     !$OMP END PARALLEL

   endif

   call CPU_TIME(finish)
   write(*,*) " "
   write(*,'(A, F12.4, A)') "Elapsed time for Computation of all the induced velocities = ", finish - start, " s"

   comp_time_vi(iter_temp) = finish - start

   ! 2 bis) FLAPPING DYNAMICS

   do l = 1, n_rotors
     call flapping(beta(l, :), beta_dot(l, :, :), beta_ddot(l, :, :), Ma(l, :), l)
     write (20,*) beta(l, :)
     write (20,*) beta_dot(l, :, 1)
     write (20,*) Ma(l, :)
   enddo

   call CPU_TIME(start)

!  ** Output to file: tutte le quantità sono aggiornate allo stesso istante
!     temporale (ora o mai più!!)

   if(allocated(particles)) then
     call write2file(iter_temp, psi, blade_panels, vi_on_blades, wake_panels, particles)
     call write2screen(iter_temp, blade_panels, wake_panels, particles, distance_max)
   else
     call write2file(iter_temp, psi, blade_panels, vi_on_blades, wake_panels)
     call write2screen(iter_temp, blade_panels, wake_panels)
   endif

   call CPU_TIME(finish)
   write(*,*) " "
   write(*,'(A, F10.4, A)') "Elapsed time for storing of data = ", finish - start, " s"

!  3) AVANZO AL TIME STEP SUCCESSIVO: sposto scia, particelle e pala
   time = time + dt
!   write(*,*)  'Time = ', time

   do j = 1, n_rotors
     psi(j, :) = psi(j, :) + omega(j) * dt;
     where(psi(j, :) > 2 * pi)
       psi(j, :) = psi(j, :) - 2 * pi
     endwhere
     where(psi(j, :) <  -2 * pi)
      psi(j, :) = psi(j, :) + 2 * pi
     endwhere
   enddo

!  a) Sposto particelle e scia

     if ( iter_temp > n_panels + 1 ) then
       do l = 1, n_rotors
        do j = 1, n_blades
           do i = 1, n_elements
              do k = 1, n_panels
                 call translations_wake_panels_s(wake_panels(l, j, i, k))
              end do
           end do
        end do
      enddo
   elseif ( iter_temp > 1 ) then
      do l = 1, n_rotors
        do j = 1, n_blades
           do i = 1, n_elements
              do k = 1, iter_temp - 1
                 call translations_wake_panels_s(wake_panels(l, j, i, k))
              end do
           end do
        end do
      enddo
   end if

   distance_max = 0.

   intersection_particles = 0

   if(allocated(particles)) then

     if(iter_temp<10) then
       format_string = "(A33,I1,A4)"
     elseif (iter_temp<100 .and. iter_temp>=10) then
       format_string = "(A33,I2,A4)"
     elseif (iter_temp<1000 .and. iter_temp>=100) then
       format_string = "(A33,I3,A4)"
     elseif (iter_temp<10000 .and. iter_temp>=1000) then
       format_string = "(A33,I4,A4)"
     else
       format_string = "(A33,I5,A4)"
     endif

     write(filename, format_string) "solution/particle_intersect/coord", iter_temp,".dat"
     OPEN ( UNIT = 40, FILE = filename, ACTION= 'READWRITE', STATUS= 'REPLACE')

     do j = 1, n_rotors
      call blade_vertex(psi(j, :), beta(j, :),  blade_extremes(j, :, :, :), center((1 + (j-1) * 3):(3*j)), j)
    enddo
       do l = 1, n_rotors
        do j = 1, size(particles, 2)
           do i = 1, size(particles, 3)
              call translations_particles(particles(l, j, i), blade_extremes, l , j, i, time, distance, check_intersection)
              if(check_intersection) then
                 intersection_particles = intersection_particles + 1
                 write(40,*) particles(l, j, i)%position
              endif
              distance_max = max(distance, distance_max)
           end do
        end do
      enddo
      close(40)
!      write(*,*) "Number of particles intersected = ", intersection_particles
      do l = 1, n_rotors
        do j = 1, n_blades
           do i = 1, n_elements
              wake_segment(l, j, i)%vertex(1, :) = wake_panels(l, j, i, 2)%vertex(4, :)
              wake_segment(l, j, i)%vertex(2, :) = wake_panels(l, j, i, 2)%vertex(3, :)
           end do
        end do
      enddo
   end if

!  b) Sposto pala
 do l = 1, n_rotors
   do j = 1, n_blades
      do i = 1, n_elements
         blade_panels(l,j,i)%vertex = blade_panel_vertex(bp_vertex(1,i,:), bp_vertex(1,i+1,:), &
                                        bp_vertex(2,i+1,:), bp_vertex(2,i,:), psi(l, j), beta(l, j), l) + &
                                        spread(center((1 + (l-1)*3 ): (l * 3)), dim=1, ncopies=4)
         blade_panels(l,j,i)%collocation = global_blade_vector(psi(l, j), beta(l, j), [0.d0, y(l, i)], l) + &
                                           center((1 + (l-1)*3 ): (l * 3))
      end do
   end do
 enddo

   call CPU_TIME(finish)
   write(*,*) " "
   write(*,'(A, F6.4, A)') "Elapsed time for moving all elements = ", finish - start, " s"


   call CPU_TIME(start)

!  4) RIORDINO le strutture:

!  a) Sposto la posizione della memoria della scia
if(iter_temp > 1 .and. iter_temp <= n_panels + 1) then
  do l = 1, n_rotors
    do i = iter_temp, 2, -1
       wake_panels(l,:,:,i) = wake_panels(l,:,:,i-1)
    enddo
  enddo
else
  do l = 1, n_rotors
    do i = n_panels + 1, 2, -1
       wake_panels(l,:,:,i) = wake_panels(l,:,:,i-1)
    enddo
  enddo
endif

!  b) Cucio il primo fronte della pala alla scia generando il nuovo fronte di
!     vortici di scia

  do l = 1, n_rotors
   do j = 1, n_blades
      do i = 1, n_elements

         wake_panels(l,j,i,1)%gamma = blade_panels(l,j,i)%gamma
         wake_panels(l,j,i,1)%vertex(1,:) = blade_panels(l,j,i)%vertex(4,:)
         wake_panels(l,j,i,1)%vertex(2,:) = blade_panels(l,j,i)%vertex(3,:)

         wake_panels(l, j,i,1)%vp_old = velocity_vertex34(:,:, l, j, i)
         wake_panels(l, j,i,1)%vertex(3:4,:) = new_vertex12(:, :, l, j, i)

         if(iter_temp <= n_panels + 1 .and. iter_temp > 1) then
            do k = iter_temp, 2, -1
               wake_panels(l,j,i,k)%vertex(2,:) = wake_panels(l,j,i,k-1)%vertex(3,:)
               wake_panels(l,j,i,k)%vertex(1,:) = wake_panels(l,j,i,k-1)%vertex(4,:)
            enddo

         else
            do  k = n_panels + 1, 2, -1
               wake_panels(l,j,i,k)%vertex(2,:) = wake_panels(l,j,i,k-1)%vertex(3,:)
               wake_panels(l,j,i,k)%vertex(1,:) = wake_panels(l,j,i,k-1)%vertex(4,:)
            enddo
         endif
      enddo
   enddo
 enddo

  !  c) Creo nuove particelle dall'ultimo fronte della scia
    if(iter_temp >= n_panels+1) call particle_generation() ! It's time to add particles because there are 3 wake panels

   call CPU_TIME(finish)
   write(*,*) " "
   write(*,'(A, F6.4, A)') "Elapsed time for organizing structures = ", finish - start, " s"



   !!!!!! OUTPUT SECTION !!!!!!

   ! Output: time_iteration, Max gamma on blade, max gamma on blade panels, max gamma on wake panel
   !         max position of particles, max vp of particles, max vp of wake panels, max vp of blade panels
   !         file with position of particles, blades and wake panels each timestep, ct at the end


   call CPU_TIME(finish_iter)
   write(*,*) " "
   write(*,'(A, F12.4, A)') "Elapsed time for Computation in this time-step = ", finish_iter - start_iter, " s"

   ! Calcolo prestazioni
   integral = 0.
   Ma = 0.
   torque = 0.


   do l = 1, n_rotors
      do i = 1, n_blades
       do j = 1, n_elements

         integral(l) = integral(l) + blade_panels(l, i, j)%dFz
         Ma(l, i) = Ma(l, i) + (y(l, j) - e(l)) * blade_panels(l, i, j)%dFz
         torque(l) = torque(l) + blade_panels(l, i, j)%dQ

       enddo
      enddo
      ct(l, iter_temp) = integral(l) / (rho * A(l) * (Vtip(l)**2))
      cp(l, iter_temp) = torque(l) / (rho * A(l) * Vtip(l)**2 * R(l))
      write(*,*) "ct rotor n°",l,"= ", ct(l, iter_temp)
      write(*,*) "cp rotor n°",l,"= ", cp(l, iter_temp)
   enddo



   iter_temp = iter_temp + 1   ! Update of temporal iteration

enddo ! Time loop

close(20)

filename = "solution/CT.dat"
OPEN ( UNIT = 1, FILE = filename, ACTION= 'READWRITE', STATUS= 'REPLACE')

filename = "solution/comp_time_vi.dat"
OPEN ( UNIT = 2, FILE = filename, ACTION= 'READWRITE', STATUS= 'REPLACE')

filename = "solution/comp_time_gamma.dat"
OPEN ( UNIT = 3, FILE = filename, ACTION= 'READWRITE', STATUS= 'REPLACE')

filename = "solution/num_particles.dat"
OPEN ( UNIT = 4, FILE = filename, ACTION= 'READWRITE', STATUS= 'REPLACE')

filename = "solution/CP.dat"
OPEN ( UNIT = 5, FILE = filename, ACTION= 'READWRITE', STATUS= 'REPLACE')

do l = 1, n_rotors
  do i = 1, size(comp_time_vi)
    write(1, *) ct(l, i)
    write(5, *) cp(l, i)
  enddo
enddo
do i = 1, size(comp_time_vi)
  write(2, *) comp_time_vi(i)
  write(3, *) comp_time_gamma(i)
  write(4, *) num_particles(i)
enddo

close(1)
close(2)
close(3)
close(4)
close(5)

call CPU_TIME(finish_int)
write(*,*) " "
write(*,'(A, F12.4, A)') "Elapsed time for Total Computation = ", finish_int - start_int, " s"

contains


! Useful subroutines -----------------------------------------------------------
subroutine initial_allocation()

implicit none

   call CPU_TIME(start)

   ! Computing maximum number of particles
   max_temp_iter = ceiling(tf / dt)  ! Computing max number of temporal iterations
   !n_particles = (max_temp_iter - n_panels - 1) * (n_elements + 1)

   ! Allocation of the psi and y vectors
   if(.not.allocated(psi)) allocate(psi(n_rotors, n_blades))
   if(.not.allocated(Ma)) allocate(Ma(n_rotors, n_blades))
   if(.not.allocated(beta_dot)) allocate(beta_dot(n_rotors, n_blades, 2))
   if(.not.allocated(beta_ddot)) allocate(beta_ddot(n_rotors, n_blades, 2))
   if(.not.allocated(y)) allocate(y(n_rotors, n_elements))
   if(.not.allocated(bp_vertex)) allocate(bp_vertex(2, n_elements+1, 2))
   if(.not.allocated(blade_extremes)) allocate(blade_extremes(n_rotors, n_blades, 4, 3))

   ! Allocation of the particles vector
   !if(.not.allocated(particles)) allocate(particles(n_particles))

   n_particles_toadd = n_blades * (n_elements + 1)  ! Computing the number of particles that has to be added each time step

   ! Allocation of the new_particles vector
   if(.not.allocated(new_particles)) allocate(new_particles(n_particles_toadd))

   ! Allocation of the wake panel matrix
   if(.not.allocated(wake_panels)) allocate(wake_panels(n_rotors, n_blades, n_elements, n_panels + 1))

   ! Allocation of the vector for the residual gamma
   if(.not.allocated(wake_segment)) allocate(wake_segment(n_rotors, n_blades, n_elements))

   ! Allocation of the blade panel matrix
   if(.not.allocated(blade_panels)) allocate(blade_panels(n_rotors, n_blades, n_elements))
   ! Allocating blade panels induced velocity
   if(.not.allocated(vi_on_blades)) allocate(vi_on_blades(n_rotors, n_blades, n_elements+1, 3))

   if(.not.allocated(new_vertex12)) allocate(new_vertex12(2, 3,  n_rotors, n_blades, n_elements))
   if(.not.allocated(velocity_vertex34)) allocate(velocity_vertex34(2, 3,  n_rotors, n_blades, n_elements))

   if(.not.allocated(ct)) allocate(ct(n_rotors, ceiling(tf/dt)))
   if(.not.allocated(cp)) allocate(cp(n_rotors, ceiling(tf/dt)))
   if(.not.allocated(Re)) allocate(Re(n_rotors, ceiling(tf/dt)))
   if(.not.allocated(comp_time_vi)) allocate(comp_time_vi(ceiling(tf/dt)))
   if(.not.allocated(comp_time_gamma)) allocate(comp_time_gamma(ceiling(tf/dt)))
   if(.not.allocated(num_particles)) allocate(num_particles(ceiling(tf/dt)))
   if(.not.allocated(integral)) allocate(integral(n_rotors))
   if(.not.allocated(torque)) allocate(torque(n_rotors))

   call CPU_TIME(finish)
   write(*,*) " "
   write(*,'(A, F6.4, A)') "Elapsed time for Data Allocation = ", finish - start, " s"

end subroutine initial_allocation

subroutine particle_generation()

  !!!! Generation of the new particles
!  write(*,*) " Creation of particles"

  if(allocated(particles))  then ! Check if it's the first trasformation eva

    ! Allocation of the matrix that will contain the particles already existing
    allocate(particles_copy(n_rotors, n_particles_toadd, size(particles, 3)))

    ! Storing already existing particles
    particles_copy = particles

    ! Deallocating the matrix particles
    deallocate(particles)

    ! Reallocating particles matrix with a new column
    allocate(particles(n_rotors, n_particles_toadd, size(particles_copy, 3) + 1))

    ! Cycling on all the rotors
    do i = 1, n_rotors

      call panel2particles(wake_panels(i, :, :, n_panels + 1), i,   new_particles, wake_segment(i, :, :))

      ! Re-storing the old particles in the tail of the particles matrix
      particles(i, :, 2:size(particles, 3)) = particles_copy(i, :, :)

      ! Add new particles in the first column
      particles(i, :, 1) = new_particles

    enddo

    ! Deallocating the matrix for the particles copy
    deallocate(particles_copy)

  else

    ! Allocation of the first column of particles
    allocate(particles(n_rotors, n_particles_toadd, 1))

    ! Cycling on all the rotors
    do i = 1, n_rotors

      call panel2particles(wake_panels(i, :, :, n_panels + 1), i,   new_particles, wake_segment(i, :, :))

      ! Add new particles in the first column
      particles(i, :, 1) = new_particles

    enddo


  endif


end subroutine particle_generation

endprogram vortex_particle
