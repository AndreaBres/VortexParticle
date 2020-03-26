module mod_output

  use mod_initialize
  use mod_geometry

  implicit none

contains

  subroutine write2file(time_step, psi, blade_panels, vi_on_blades, wake_panels, particles)

    type(blade_panel), dimension(:, :, :), intent(in) :: blade_panels
    type(wake_panel), dimension(:, :, :, :), intent(in), optional :: wake_panels
    type(particle), dimension(:, :, :), intent(in), optional :: particles
    real(kind=8), dimension(:, :, :, :), intent(in) :: vi_on_blades
    real(kind=8), dimension(:, :), intent(in) :: psi
    integer, intent(in) :: time_step

    real(kind=8), dimension(size(vi_on_blades, 1),&
                            size(vi_on_blades, 2),&
                            size(vi_on_blades, 3),&
                            size(vi_on_blades, 4))  :: vi_2write
    character (len = 50) :: filename, format_string
    integer :: i, j, k, l, m

!    write(*,*) "----------------------------------------------------"
!    write(*,*) " SALVATAGGIO DATI H, U E V"
!    write(*,*) "----------------------------------------------------"


    if(time_step<10) then
      format_string = "(A26,I1,A4)"
    elseif (time_step<100 .and. time_step>=10) then
      format_string = "(A26,I2,A4)"
    elseif (time_step<1000 .and. time_step>=100) then
      format_string = "(A26,I3,A4)"
    elseif (time_step<10000 .and. time_step>=1000) then
      format_string = "(A26,I4,A4)"
    else
      format_string = "(A26,I5,A4)"
    endif

    write(filename, format_string) "solution/blade_coord/coord", time_step,".dat"
    OPEN ( UNIT = 1, FILE = filename, ACTION= 'READWRITE', STATUS= 'REPLACE')

    if(time_step<10) then
      format_string = "(A32,I1,A4)"
    elseif (time_step<100 .and. time_step>=10) then
      format_string = "(A32,I2,A4)"
    elseif (time_step<1000 .and. time_step>=100) then
      format_string = "(A32,I3,A4)"
    elseif (time_step<10000 .and. time_step>=1000) then
      format_string = "(A32,I4,A4)"
    else
      format_string = "(A32,I5,A4)"
    endif

    write(filename, format_string) "solution/vi_on_blade/vi_on_blade", time_step,".dat"
    OPEN ( UNIT = 10, FILE = filename, ACTION= 'READWRITE', STATUS= 'REPLACE')


    if(time_step<10) then
      format_string = "(A20,I1,A4)"
    elseif (time_step<100 .and. time_step>=10) then
      format_string = "(A20,I2,A4)"
    elseif (time_step<1000 .and. time_step>=100) then
      format_string = "(A20,I3,A4)"
    elseif (time_step<10000 .and. time_step>=1000) then
      format_string = "(A20,I4,A4)"
    else
      format_string = "(A20,I5,A4)"
    endif

    write(filename, format_string) "solution/blade_Fz/Fz", time_step,".dat"
    OPEN ( UNIT = 5, FILE = filename, ACTION= 'READWRITE', STATUS= 'REPLACE')

    if(time_step<10) then
      format_string = "(A31,I1,A4)"
    elseif (time_step<100 .and. time_step>=10) then
      format_string = "(A31,I2,A4)"
    elseif (time_step<1000 .and. time_step>=100) then
      format_string = "(A31,I3,A4)"
    elseif (time_step<10000 .and. time_step>=1000) then
      format_string = "(A31,I4,A4)"
    else
      format_string = "(A31,I5,A4)"
    endif

      write(filename, format_string) "solution/wake_panel_coord/coord", time_step,".dat"
      OPEN ( UNIT = 2, FILE = filename, ACTION= 'READWRITE', STATUS= 'REPLACE')

      if(time_step<10) then
        format_string = "(A29,I1,A4)"
      elseif (time_step<100 .and. time_step>=10) then
        format_string = "(A29,I2,A4)"
      elseif (time_step<1000 .and. time_step>=100) then
        format_string = "(A29,I3,A4)"
      elseif (time_step<10000 .and. time_step>=1000) then
        format_string = "(A29,I4,A4)"
      else
        format_string = "(A29,I5,A4)"
      endif

      write(filename, format_string) "solution/blade_coord_Fz/coord", time_step,".dat"
      OPEN ( UNIT = 8, FILE = filename, ACTION= 'READWRITE', STATUS= 'REPLACE')


    do m = 1, n_rotors
      do i = 1, n_blades
        write(5, *) blade_panels(m, i, :)%dFz

        do j = 1, n_elements
          write (8, *) (blade_panels(m, i, j)%vertex(1, 1:2) + blade_panels(m, i, j)%vertex(2, 1:2)) / 2
          do k = 1, 4
            write (1, *) blade_panels(m, i, j)%vertex(k, :)
          enddo

          if(present(wake_panels)) then

            !!!! COORREGGERE IN BASE ALL'INIZIO DEL TIME_STEP !!!!
            if(time_step <= n_panels+1 .and. time_step > 1) then
              do k = 1, time_step
                do l = 1, 4
                  write (2, *) wake_panels(m, i, j, k)%vertex(l, :)
                enddo

              enddo
            else

              do k = 1, n_panels+1

                do l = 1, 4
                  write (2, *) wake_panels(m, i, j, k)%vertex(l, :)
                enddo

              enddo

            endif
          endif
        enddo
        do j = 1, n_elements + 1
          call xyz2blade(psi(m, i), beta(m, i), vi_on_blades(m, i, j, :), m, vi_2write(m, i, j, :)) ! Induced velocity in blade coordinate
          write(10, *) vi_2write(m, i, j, :)
        enddo
      enddo
    enddo


    close(1)
    close(2)
    close(5)
    close(8)
    close(10)

    if (present(particles)) then

      if(time_step<10) then
        format_string = "(A29,I1,A4)"
      elseif (time_step<100 .and. time_step>=10) then
        format_string = "(A29,I2,A4)"
      elseif (time_step<1000 .and. time_step>=100) then
        format_string = "(A29,I3,A4)"
      elseif (time_step<10000 .and. time_step>=1000) then
        format_string = "(A29,I4,A4)"
      else
        format_string = "(A29,I5,A4)"
      endif

      write(filename, format_string) "solution/particle_coord/coord", time_step,".dat"
      OPEN ( UNIT = 3, FILE = filename, ACTION= 'READWRITE', STATUS= 'REPLACE')

      do m = 1, n_rotors
        do i = 1, size(particles, 2)
          do j = 1, size(particles, 3)

            write(3, *) particles(m, i, j)%position

          enddo
        enddo
      enddo

      close(3)

!      if(time_step<10) then
!        format_string = "(A31,I1,A4)"
!      elseif (time_step<100 .and. time_step>=10) then
!        format_string = "(A31,I2,A4)"
!      elseif (time_step<1000 .and. time_step>=100) then
!        format_string = "(A31,I3,A4)"
!      elseif (time_step<10000 .and. time_step>=1000) then
!        format_string = "(A31,I4,A4)"
!      else
!        format_string = "(A31,I5,A4)"
!      endif
!
!      write(filename, format_string) "solution/particle_coords/coords", time_step,".dat"
!      OPEN ( UNIT = 4, FILE = filename, ACTION= 'READWRITE', STATUS= 'REPLACE')
!
!      do m = 1, n_rotors
!        do j = 1, size(particles, 3)
!          do i = 1, size(particles, 2)

!            write(4, *) particles(m, i, j)%position

!          enddo
!        enddo
!      enddo

!      close(4)

    endif

  endsubroutine  write2file


  subroutine write2screen(time_step, blade_panels, wake_panels, particles, distance_max)

    type(blade_panel), dimension(:, :, :), intent(in) :: blade_panels
    type(wake_panel), dimension(:, :, :, :), intent(in), optional :: wake_panels
    type(particle), dimension(:, :, :), intent(in), optional :: particles
    real(kind=8), intent(in), optional :: distance_max
    integer, intent(in) :: time_step

    real(kind=8) :: max_gamma_blade_panel, max_gamma_wake_panel, max_gamma_particle,&
                    max_distance_blade_panel, max_distance_wake_panel, max_distance_particle,&
                    min_gamma_blade_panel

    integer :: i, j, k, m

    ! Computing maximum value of gamma for all of the blade panels
    max_gamma_blade_panel = maxval(blade_panels%gamma)
    min_gamma_blade_panel = minval(blade_panels%gamma)

    max_distance_wake_panel = 0
    max_distance_blade_panel = 0
    do m = 1, n_rotors
      do i = 1, n_blades
        do j = 1, n_elements

          ! Computing maximum distance of a point on the blade panels
          max_distance_blade_panel = max(max_distance_blade_panel, maxval(norm2(blade_panels(m, i, j)%vertex, 2)))

          !!!! COORREGGERE IN BASE ALL'INIZIO DEL TIME_STEP !!!!
          if(time_step <= n_panels+1 .and. time_step > 1) then
            do k = 1, time_step - 1

              ! Computing maximum distance of a point on the blade panels
              max_distance_wake_panel = max(max_distance_wake_panel, maxval(norm2(wake_panels(m, i, j, k)%vertex, 2)))

            enddo

          else

            do k = 1, n_panels

              ! Computing maximum distance of a point on the blade panels
              max_distance_wake_panel = max(max_distance_wake_panel, maxval(norm2(wake_panels(m, i, j, k)%vertex, 2)))

            enddo

          endif

        enddo
      enddo
    enddo

    ! Computing maximum value of gamma for all of the wake panels
    max_gamma_wake_panel = maxval(wake_panels%gamma)

    if(present(particles)) then
      ! Computing maximum value of gamma and of distance for all of the particles
      max_gamma_particle = 0
      max_distance_particle = 0
      do m = 1, n_rotors
        do i = 1, size(particles, 2)
          do j = 1, size(particles, 3)

            max_gamma_particle = max(max_gamma_particle, norm2(particles(m, i, j)%alpha_p))
            max_distance_particle = max(max_distance_particle, norm2(particles(m, i, j)%position))

          enddo
        enddo
      enddo
    endif

    write(*,*) "--------------------------------------------------------"
    write(*,*) "Computing the solution at time = ", time_step * dt
    write(*,*) "Maximum value of gamma on the blade = ", max_gamma_blade_panel
    write(*,*) "Minimum value of gamma on the blade = ", min_gamma_blade_panel
    if(time_step > 1) then
      write(*,*) "Maximum value of gamma on the wake panel = ", max_gamma_wake_panel
      if (present(particles)) then
        write(*,*) "Maximum value of gamma on the particle = ", max_gamma_particle
      endif
    endif
    write(*,*) "--------------------------------------------------------"
    write(*,*) "Maximum distance of a point on the blade = ", max_distance_blade_panel
    if(time_step > 1) then
      write(*,*) "Maximum distance of a point on the wake panel = ", max_distance_wake_panel
      if (present(particles)) then
        write(*,*) "Maximum distance of the particles = ", max_distance_particle
        write(*,*) "Max Distance travelled by a particle = ", distance_max
      endif
    endif
    write(*,*) "--------------------------------------------------------"



  endsubroutine write2screen

endmodule mod_output
