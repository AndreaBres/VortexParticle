module mod_check

  ! Module containing all of the subroutine that are used in order to check something

  use mod_geometry

  implicit none

contains

  ! which_side: gives which side a point is wrt to a plane

  ! Input:   point: vector of coordinates of the point
  !          plane_vertexes: matrix of 4x3 elements containings the coordinates of the 4 vertexes of the blade

  ! Output:  side: it is either -1, 0 or +1, depending if the point is respectively
  !                below, on or above the plane
  subroutine which_side(point, plane_vertexes,   side)

    real(kind=8), dimension(:, :), intent(in)  ::  plane_vertexes
    real(kind=8), dimension(:),    intent(in)  ::  point

    real(kind=8),                  intent(out) ::  side

    real(kind=8), dimension(3) :: normal, vec1, vec2, congiungente
    real(kind=8) :: scalar_product

    ! Initialization of the side variable
    side = 1.

    ! Computing the 2 vectors that belongs to the plane of the i-th blade
    vec1 = plane_vertexes(2, :) - plane_vertexes(1, :)
    vec2 = plane_vertexes(3, :) - plane_vertexes(2, :)

    ! Computing the normal vector to the plane of the i-th blade
    call cross_product(vec1, vec2,   normal)

    normal = normal / norm2(normal)

    ! Computing the congiungente between the starting position and a vertex of the plane of the i-th blade
    congiungente = point - plane_vertexes(1, :)

    ! Checking if the starting position is below the plane of the i-th blade
    scalar_product = dot_product(congiungente, normal)
    side = sign(side, scalar_product)

  endsubroutine which_side


  subroutine intersection_line_plane(starting_position, ending_position, plane_vertexes,   point_onplane)

    real(kind=8), dimension(:),    intent(in)  ::  starting_position, ending_position
    real(kind=8), dimension(:, :), intent(in)  ::  plane_vertexes

    real(kind=8), dimension(:),    intent(out) ::  point_onplane

    real(kind=8), dimension(3) :: normal, vec1, vec2, congiungente
    real(kind=8) :: scalar_product, distance

    ! Computing the 2 vectors that belongs to the plane
    vec1 = plane_vertexes(1, :) - plane_vertexes(2, :)
    vec2 = plane_vertexes(3, :) - plane_vertexes(2, :)

    ! Computing the normal vector to the plane of the i-th blade
    call cross_product(vec1, vec2,   normal)

    normal = normal / norm2(normal)

    ! Computing the congiungente between the starting position and the ending position
    congiungente = ending_position - starting_position

    ! Computing the distance on the congiungente
    distance = dot_product(plane_vertexes(1, :) - starting_position, normal) / &
               dot_product(congiungente, normal)

    ! Computing the coordinates of the intersection point
    point_onplane = distance * congiungente + starting_position


  endsubroutine intersection_line_plane



  subroutine projection_point_plane(ending_position, plane_vertexes)

    real(kind=8), dimension(:),    intent(inout)  ::   ending_position
    real(kind=8), dimension(:, :), intent(in)     ::   plane_vertexes

    real(kind=8), dimension(3) ::  point_onplane
    real(kind=8), dimension(3) :: normal, vec1, vec2, congiungente
    real(kind=8) :: distance

    ! Computing the 2 vectors that belongs to the plane
    vec1 = plane_vertexes(2, :) - plane_vertexes(1, :)
    vec2 = plane_vertexes(3, :) - plane_vertexes(2, :)

    ! Computing the normal vector to the plane of the i-th blade
    call cross_product(vec1, vec2,   normal)

    normal = normal / norm2(normal)

    ! Computing the congiungente between the ending position and the plane
    congiungente = -normal

    ! Computing the distance on the congiungente
    distance = dot_product(plane_vertexes(1, :) - ending_position, normal) / &
               dot_product(congiungente, normal)

    ! Computing the coordinates of the projection point
    point_onplane = distance * congiungente + ending_position

    ! Update of the ending position
    ending_position = point_onplane


  endsubroutine projection_point_plane

  subroutine expand_blade(blade_vertexes, radius, width_here,   blade_vertexes_expanded)

    ! expand_blade: subroutine that modigies the coordinates of the vertexes of the blades in order
    !               to take into account of the radius of the particle.
    !
    ! INPUT:  blade_vertexes: coords of the 4 vertexes of the blade
    !         radius:         radius of the particle
    !
    ! OUTPUT: blade_vertexes_expanded: coords of the 8 new vertexes of the blade
    !                                  (it will now have 2 parallel planes due to the width)

    real(kind=8), dimension(:, :), intent(in)  :: blade_vertexes
    real(kind=8),                  intent(in)  :: radius, width_here
    real(kind=8), dimension(:, :, :), intent(out) :: blade_vertexes_expanded

    real(kind=8), dimension(size(blade_vertexes, 1), size(blade_vertexes, 2)) :: blade_vertex_enlarged
    real(kind=8), dimension(3) :: direction, normal, vec1, vec2

    !                             TOP VIEW (enlargement)
    !
    !   Blade before                                  Blade after
    !                                                  r         r
    !                                                 |-|       |-|
    !   2.-------.3                                  2.-----------.3
    !    |       |                                    |           | ) r
    !    |       |                                    |           |
    !    |       |                                    |           |
    !    |       |                                    |           |
    !    |       |                                    |           |
    !   1!_______!4                                   |           |
    !       HUB                                       |           |
    !                                                1!___________!4) r
    !                                                      HUB
    !
    !                           TIP VIEW (increase in width)
    !
    !   Blade before                                  Blade after
    !   2.-------.3                                  2.___________.3
    !                                                 |           | ) r
    !                                                6!___________!7) r

    !! Enlarging of the segment 1-2
    direction = blade_vertexes(2, :) - blade_vertexes(1, :)  ! Computing the direction of the segment
    direction = direction / norm2(direction)
    ! Computing new position of the vertex 2
    blade_vertex_enlarged(2, :) = blade_vertexes(2, :) + direction * radius
    ! Computing new position of the vertex 1
    blade_vertex_enlarged(1, :) = blade_vertexes(1, :) - direction * radius

    !!! HP : BLADE IS A PLANE ---> SEGMENT 1-2 // SEGMENT 3-4
    ! Computing new position of the vertex 3
    blade_vertex_enlarged(3, :) = blade_vertexes(3, :) + direction * radius
    ! Computing new position of the vertex 4
    blade_vertex_enlarged(4, :) = blade_vertexes(4, :) - direction * radius

    direction = blade_vertexes(3, :) - blade_vertexes(2, :)  ! Computing the direction of the segment
    direction = direction / norm2(direction)
    ! Computing new position of the vertex 2
    blade_vertex_enlarged(2, :) = blade_vertexes(2, :) - direction * radius
    ! Computing new position of the vertex 1
    blade_vertex_enlarged(1, :) = blade_vertexes(1, :) - direction * radius

    !!! HP : BLADE IS A PLANE ---> SEGMENT 2-3 // SEGMENT 1-4
    ! Computing new position of the vertex 3
    blade_vertex_enlarged(3, :) = blade_vertexes(3, :) + direction * radius
    ! Computing new position of the vertex 4
    blade_vertex_enlarged(4, :) = blade_vertexes(4, :) + direction * radius

    !! Introducing width
    ! Computing the 2 vectors that belongs to the plane
    vec1 = blade_vertexes(2, :) - blade_vertexes(1, :)
    vec2 = blade_vertexes(3, :) - blade_vertexes(2, :)
    !write(*,*) "vec1 = ", vec1
    !write(*,*) "vec2 = ", vec2
    !write(*,*) "blade_vertexes(1, :) = ", blade_vertexes(1, :)
    !write(*,*) "blade_vertexes(2, :) = ", blade_vertexes(2, :)
    !write(*,*) "blade_vertexes(3, :) = ", blade_vertexes(3, :)

    ! Computing the normal vector to the plane of the i-th blade
    call cross_product(vec1, vec2,   normal)
      !write(*,*) "normal = ", normal
    normal = normal / norm2(normal)
    !write(*,*) "normal = ", normal

    ! Computing the top plane vertexes
    blade_vertexes_expanded(2, 1, :) = blade_vertex_enlarged(1, :) + normal * (radius + width_here)
    blade_vertexes_expanded(2, 2, :) = blade_vertex_enlarged(2, :) + normal * (radius + width_here)
    blade_vertexes_expanded(2, 3, :) = blade_vertex_enlarged(3, :) + normal * (radius + width_here)
    blade_vertexes_expanded(2, 4, :) = blade_vertex_enlarged(4, :) + normal * (radius + width_here)

    ! Computing the bottom plane vertexes
    blade_vertexes_expanded(1, 1, :) = blade_vertex_enlarged(1, :) - normal * (radius + width_here)
    blade_vertexes_expanded(1, 2, :) = blade_vertex_enlarged(4, :) - normal * (radius + width_here)
    blade_vertexes_expanded(1, 3, :) = blade_vertex_enlarged(3, :) - normal * (radius + width_here)
    blade_vertexes_expanded(1, 4, :) = blade_vertex_enlarged(2, :) - normal * (radius + width_here)

  endsubroutine expand_blade


  subroutine check_distance (vertexes, distance, starting_position,    check_dist)

    ! check_distance: subroutine that checks if the particle is near enough to the blade
    !                 in order to make it possible to go through it
    !
    ! INPUT:   vertexes:  matrix 2x4x3 containings coords of the blade extended by the radius
    !          distance: distance traveled by the particle during its motion
    !          starting_position:  vector containings coords of the initial position of the particle
    !
    ! OUTPUT:  check_dist: true if distance > min(distance_starting_point - blade_spigoli_extended), false otherwise

    real(kind=8), dimension(:, :, :), intent(in)  ::  vertexes
    real(kind=8), dimension(:),       intent(in)  ::  starting_position
    real(kind=8),                     intent(in)  ::  distance
    logical,                          intent(out) ::  check_dist

    real(kind=8), dimension(6, 4, 3) :: planes
    real(kind=8), dimension(3) :: spigolo, spigolo_1, congiungente, normal
    real(kind=8) :: min_dist, dist
    integer :: i, k, i_l, i_r

    ! Initialization of the min_dist with the value of the distance of the particle
    min_dist = distance
    ! Initialization of the check variable
    check_dist = .false.

    call planes_gen (vertexes, planes)

    ! Computing minimum distance between the particle and all the planes of the blade
    do i = 1, size(planes, 1)


        call distance_point2plane(planes(i, :, :), starting_position, dist)
        !WRITE(*,*) "dist = ", dist
        !writE(*,*) "min_dist = ", min_dist

        ! Check if the minimum distance between the particle and the blade is lower than the distance traveled by the particle
        if (dist < min_dist) then
          ! It can already exit the cycle
          check_dist = .true.
          !write(*,*) "the plane with probable intersection is plane n°", i
          !write(*,*) "planes(i, i_r, :)  = ", planes(i, i_r, :)
          !write(*,*) "planes(i, i_l, :) = ", planes(i, i_l, :)
          !write(*,*) "planes(i, i_r+1, :) = ", planes(i, i_r+1, :)

          go to 200
        endif
    enddo

    200  continue

  endsubroutine check_distance

  subroutine distance_point2plane(plane, starting_position, dist)

    real(kind=8), dimension(:, :), intent(in)    :: plane
    real(kind=8), dimension(:),    intent(in)    :: starting_position
    real(kind=8),                  intent(inout) :: dist

    real(kind=8), dimension(3) :: spigolo, spigolo_1, congiungente, normal



    spigolo = plane(2, :) - plane(1, :)
    spigolo_1 = plane(3, :) - plane(2, :)
    spigolo = spigolo / norm2(spigolo)
    spigolo_1 = spigolo_1 / norm2(spigolo_1)
    congiungente = starting_position - plane(1, :)

    call cross_product(spigolo , spigolo_1, normal)
    !write(*,*) "panels = ", i
    !write(*,*) "normal = ", normal
    !write(*,*) "spigolo = ", spigolo
    !write(*,*) "spigolo_1 =", spigolo_1
    !write(*,*) "congiungente = ", congiungente

    dist = abs(DOT_PRODUCT(normal, congiungente)) / norm2(normal)

  endsubroutine distance_point2plane

  ! check_intersecion_blade: subroutine that checks if the intersection point is in the blade or outside it
  !
  ! INPUT:   point_onplane: vector containing the coordinates of the intersection point
  !          plane_vertex:  vector 4x3 containing the coordinates of the 4 vertices
  !
  ! OUTPUT:  check_inblade; logical variable. true if the intersection poin is in the blade, false if not

  subroutine check_intersecion_blade(point_onplane, plane_vertex,   check_inblade)

    real(kind=8), dimension(:, :), intent(in)  ::  plane_vertex
    real(kind=8), dimension(:),    intent(in)  ::  point_onplane
    logical,                       intent(out) ::  check_inblade

    ! real(kind=8), dimension(3) :: vec1, vec2
    real(kind=8), dimension(3) :: vertex1, vertex2

    ! How to find out if the intersection is in the rectangular that is defined by the blade?
    ! By checking the dot-product of the congiungente point-diagonal vertices 2 and 4:
    ! if the dot product is negative then it's inside

    ! Initialization of the logical variable
    check_inblade = .false.

    vertex1 = plane_vertex(1, :)
    vertex2 = plane_vertex(3, :)

    if((point_onplane(1) - vertex1(1))*(point_onplane(1) - vertex2(1)) <= 0) then
      if((point_onplane(2) - vertex1(2))*(point_onplane(2) - vertex2(2)) <= 0)  then
        if((point_onplane(3) - vertex1(3))*(point_onplane(3) - vertex2(3)) <= 0)  then
          check_inblade = .true.
        endif
      endif
    endif



  endsubroutine check_intersecion_blade



  ! check_blade: checks if the particle, after translation, goes through the blade

  ! Input:  starting_position: initial position of the particle
  !         ending_position: final position of the particle after translation and
  !                          before executing the check
  !         blades_extremes: 3D matrix n_blades x n_vertexes(4) x n_coords (3).
  !                          It contains the positions of the extremities of the blades
  !                          used in order to compute blades plane

  ! Output: ending_position: untouched if the check is false, otherwise it is modified

  subroutine check_blade(starting_position, blade_extremes, radius, radius_old,&
                         ending_position, check, distance)

    real(kind=8), dimension(:, :, :, :), intent(in)    ::  blade_extremes
    real(kind=8), dimension(:),          intent(inout)    ::  starting_position
    real(kind=8),                        intent(in)    ::  radius, radius_old
    real(kind=8), dimension(:),          intent(inout) ::  ending_position
    real(kind=8),                        intent(inout) ::  distance
    logical,                             intent(out)   ::  check

    real(kind=8), dimension(6, 4, 3) :: planes
    real(kind=8), dimension(2, size(blade_extremes, 3), size(blade_extremes, 4)) :: blade_extremes_extended
    real(kind=8), dimension(3) :: point_onplane, position_s
    real(kind=8) :: side_check_start, side_check_end, dist, dist_1, width_here = 0
    logical :: check_dist, check_inblade, check_inplane
    integer :: i, j, k, cont, l, t

    ! Starting iteration on the number of blades
    do j = 1, n_rotors
      do i = 1, size(blade_extremes, 2)

        ! Initialize the check variable
        check = .false.;  side_check_start = 1.;  side_check_end = 1.;   check_dist = .false.;  check_inblade = .false.

        ! Expand the vertexes of the i-th blade in order to take into account the radius of the particle
        call expand_blade(blade_extremes(j, i, :, :), radius, width_here,       blade_extremes_extended)

        ! Check if it's inside the box due to previous iterations
        call planes_gen (blade_extremes_extended, planes)

        cont = 0

        do l = 1, size(planes, 1)
          position_s = starting_position
          call check_intersecion_blade(position_s, planes(l, :, :),   check_inplane)
          if(check_inplane) then
            call which_side(starting_position, planes(l, :, :),   side_check_start)
            if (side_check_start<0) cont = cont + 1;
          endif
        enddo

        ! Check if the dot product is less than 0
        if(cont == size(planes, 1)) then ! it's inside
          ! Check if the poin is più vicino al top plane o al bottom plane del vecchio parallelepipedo
          ! Expand the vertexes of the i-th blade in order to take into account the old radius of the particle
          call expand_blade(blade_extremes(j, i, :, :), radius_old, width_here,       blade_extremes_extended)

          ! Check if it's inside the box due to previous iterations
          call planes_gen (blade_extremes_extended, planes)

          dist = 100
          do k = 1, size(planes, 1)
            call distance_point2plane(planes(k, :, :), starting_position, dist_1)
            dist = min(dist, dist_1)
            if(dist == dist_1)  t = k;
          enddo
          ! Expand the vertexes of the i-th blade in order to take into account the new radius
          call expand_blade(blade_extremes(j, i, :, :), radius, width_here,       blade_extremes_extended)
          call planes_gen (blade_extremes_extended, planes)

          call projection_point_plane(starting_position, planes(t, :, :))
        endif

        ! Check if it's possible for the particle to cover the distance between it and the i-th blade during its motion
        !call check_distance(blade_extremes_extended, distance, starting_position,    check_dist)
        check_dist = .true.

        ! If check_distance == .true. then check if it goes through the blade_plane
        if(check_dist) then

          call expand_blade(blade_extremes(j, i, :, :), radius, width_here,       blade_extremes_extended)
          call planes_gen (blade_extremes_extended, planes)

          do k = 1, size(planes, 1)

            ! Checking the side of the starting position of the particle wrt the bottom plane of the i-th blade
            call which_side(starting_position, planes(k, :, :),   side_check_start)
            call which_side(ending_position, planes(k, :, :),   side_check_end)
            if (side_check_start * side_check_end < 0) then

              !write(*,*) "Attraversa pannello n° ",k," della pala ", i

              ! Find intersection with the plane
              call intersection_line_plane(starting_position, ending_position, planes(k, :, :),   point_onplane)

              ! Checking if it's inside the rectangle
              call check_intersecion_blade(point_onplane, planes(k, :, :),   check_inblade)

              ! If the point is inside
              if(check_inblade) then

                ! Finding point on the plane where to stop the particle (projection of the point on the plane)
                call projection_point_plane(ending_position, planes(k, :, :))

                ! Changing logical variable to true
                check = .true.

                ! Updating distance travelled
                distance = norm2(ending_position - starting_position)

                ! Exit from the cycle on all the blades
                go to 600

              endif
            endif

          enddo
        endif
      enddo
    enddo

    600 continue

  endsubroutine check_blade


  subroutine check_boundary(starting_position, boundaries_here, radius, radius_old,&
                            ending_position, check, distance)

    real(kind=8), dimension(:, :),       intent(in)    ::  boundaries_here
    real(kind=8), dimension(:),          intent(inout)    ::  starting_position
    real(kind=8),                        intent(in)    ::  radius, radius_old
    real(kind=8), dimension(:),          intent(inout) ::  ending_position
    real(kind=8),                        intent(inout) ::  distance
    logical,                             intent(out)   ::  check

    real(kind=8), dimension(6, 4, 3) :: planes
    real(kind=8), dimension(2, size(boundaries_here, 1), size(boundaries_here, 2)) :: boundaries_extremes_extended
    real(kind=8), dimension(3) :: point_onplane, cong1, cong2, position_s
    real(kind=8) :: side_check_start, side_check_end, dot, width_here = 0, dist, dist_1
    logical :: check_dist, check_inboundary, check_inplane
    integer :: i, cont, t, k
    character (len = 50) :: filename, format_string

    ! Initialize the check variable
    check = .false.;  side_check_start = 1.;  side_check_end = 1.;   check_dist = .false.;  check_inboundary = .false.

    ! Expand the vertexes of the i-th blade in order to take into account the radius of the particle
    call expand_blade(boundaries_here, radius, width_here,       boundaries_extremes_extended)

    ! Check if it's inside the box due to previous iterations
    call planes_gen (boundaries_extremes_extended, planes)

    cont = 0
    check_inplane = .false.

    do i = 1, size(planes, 1)
      position_s = starting_position
      call check_intersecion_blade(position_s, planes(i, :, :),   check_inplane)
      if(check_inplane) then
        call which_side(starting_position, planes(i, :, :),   side_check_start)
        if (side_check_start<0) cont = cont + 1;
      endif
    enddo

    ! Check if the dot product is less than 0
    if(cont == size(planes, 1)) then ! it's inside

      ! Check if the poin is più vicino al top plane o al bottom plane del vecchio parallelepipedo
      ! Expand the vertexes of the i-th blade in order to take into account the old radius of the particle
      call expand_blade(boundaries_here, radius_old, width_here,       boundaries_extremes_extended)

      ! Check if it's inside the box due to previous iterations
      dist = 100
      call planes_gen (boundaries_extremes_extended, planes)

        do k = 1, size(planes, 1)

          call distance_point2plane(planes(k, :, :), starting_position, dist_1)
          dist = min(dist, dist_1)
          if (dist == dist_1) t = k;

        enddo

        ! Expand the vertexes of the i-th blade in order to take into account the new radius
        call expand_blade(boundaries_here, radius, width_here,       boundaries_extremes_extended)
        call planes_gen (boundaries_extremes_extended, planes)

        call projection_point_plane(starting_position, planes(t, :, :))

    endif

    ! Check if it's possible for the particle to cover the distance between it and the i-th blade during its motion
    !call check_distance(boundaries_extremes_extended, distance, starting_position,    check_dist)

    ! Expand the vertexes of the i-th blade in order to take into account the radius of the particle
    call expand_blade(boundaries_here, radius, width_here,       boundaries_extremes_extended)
    ! Creation of the planes
    call planes_gen (boundaries_extremes_extended, planes)

    check_dist = .true.

    ! If check_distance == .true. then check if it goes through the blade_plane
    if(check_dist) then

      do k = 1, size(planes, 1)

        ! Checking the side of the starting position of the particle wrt the bottom plane of the i-th blade
        call which_side(starting_position, planes(k, :, :),   side_check_start)
        call which_side(ending_position, planes(k, :, :),   side_check_end)
        if (side_check_start * side_check_end < 0) then

          ! Find intersection with the plane
          call intersection_line_plane(starting_position, ending_position, planes(k, :, :),   point_onplane)

          ! Checking if it's inside the rectangle
          call check_intersecion_blade(point_onplane, planes(k, :, :),   check_inboundary)

          ! If the point is inside
          if(check_inboundary) then

            ! Finding point on the plane where to stop the particle (projection of the point on the plane)
            call projection_point_plane(ending_position, planes(k, :, :))

            ! Changing logical variable to true
            check = .true.

            ! Updating distance travelled
            distance = norm2(ending_position - starting_position)

            ! Exit from the cycle on all the rotors
            go to 999

          endif
        endif
      enddo
      endif

    999 continue

  endsubroutine check_boundary

  subroutine check_hub(starting_position, radius, radius_old,&
                       ending_position, check, distance)


  real(kind=8), dimension(:),          intent(inout) ::  starting_position
  real(kind=8),                        intent(in)    ::  radius, radius_old
  real(kind=8), dimension(:),          intent(inout) ::  ending_position
  real(kind=8),                        intent(inout) ::  distance
  logical,                             intent(out)   ::  check

  real(kind=8), dimension(6, 4, 3) :: planes
  real(kind=8), dimension(2, 4, 3) :: hub_extremes_extended
  real(kind=8), dimension(n_rotors) :: hub_radius
  real(kind=8), dimension(100, 3) :: circumference
  real(kind=8), dimension(4, 3) :: vertexes
  real(kind=8), dimension(2, 3) :: diameter
  real(kind=8), dimension(3) :: point_onplane, vec1, vec2, center_here, position_s
  real(kind=8) :: dist_2hub, dx, dist_1, dist
  real(kind=8) :: side_check_start, side_check_end, dot
  integer :: i, k, l, cont, t
  logical :: check_dist, check_inhub, check_inplane

  ! The hub is a circle of radius equal to 80% of e*R
  hub_radius = fill_percent * bs * R

  check = .false.

  ! Cycling on all rotors
  do i = 1, n_rotors

    center_here = center(3* (i-1) + 1: 3 * i)

    check_dist = .true.

    if (check_dist) then

      vertexes(1, 1) = center_here(1) + hub_radius(i) ;  vertexes(3, 1) = center_here(1) - hub_radius(i) ;
      vertexes(2, 1) = center_here(1);                   vertexes(4, 1) = center_here(1);

      vertexes(1, 2) = center_here(2);                   vertexes(3, 2) = center_here(2);
      vertexes(2, 2) = center_here(2) - hub_radius(i);   vertexes(4, 2) = center_here(2) + hub_radius(i);

      vertexes(:, 3) = center_here(3) - offset_hub(i);

      ! Expand the vertexes of the i-th blade in order to take into account the radius of the particle
      call expand_blade(vertexes, radius, width(i),       hub_extremes_extended)

      ! Check if it's inside the box due to previous iterations
      call planes_gen (hub_extremes_extended, planes)

      cont = 0

      do l = 1, size(planes, 1)
        position_s = starting_position
        call check_intersecion_blade(position_s, planes(l, :, :),   check_inplane)
        if(check_inplane) then
          call which_side(starting_position, planes(l, :, :),   side_check_start)
          if (side_check_start<0) cont = cont + 1;
        endif
      enddo

      ! Check if the dot product is less than 0
      if(cont == size(planes, 1)) then ! it's inside
        ! Check if the poin is più vicino al top plane o al bottom plane del vecchio parallelepipedo
        ! Expand the vertexes of the i-th blade in order to take into account the old radius of the particle
        call expand_blade(vertexes, radius_old, width(i),       hub_extremes_extended)

        ! Check if it's inside the box due to previous iterations
        dist = 100
        call planes_gen (hub_extremes_extended, planes)

          do k = 1, size(planes, 1)

            call distance_point2plane(planes(k, :, :), starting_position, dist_1)
            dist = min(dist, dist_1)
            if (dist == dist_1) t = k;

          enddo

          ! Expand the vertexes of the i-th blade in order to take into account the new radius
          call expand_blade(vertexes, radius, width(i),       hub_extremes_extended)
          call planes_gen (hub_extremes_extended, planes)

          call projection_point_plane(starting_position, planes(t, :, :))

      endif

      ! Expand the vertexes of the i-th blade in order to take into account the radius of the particle
      call expand_blade(vertexes, radius, width(i),       hub_extremes_extended)
      ! Creation of the planes
      call planes_gen (hub_extremes_extended, planes)

      do k = 1, size(planes, 1)

        ! Checking the side of the starting position of the particle wrt the bottom plane of the i-th blade
        call which_side(starting_position, planes(k, :, :),   side_check_start)
        call which_side(ending_position, planes(k, :, :),   side_check_end)
        if (side_check_start * side_check_end < 0) then

          ! Find intersection with the plane
          call intersection_line_plane(starting_position, ending_position, planes(k, :, :),   point_onplane)

          if (k == 1 .or. k == 2) then
            ! Checking if it's inside the circle
            vec1 = planes(k, 1, :) - point_onplane
            vec2 = planes(k, 3, :) - point_onplane
            dot = dot_product(vec1, vec2)
            if(dot <= 0) check_inhub = .true.
          else
            ! Checking if it's inside the rectangle
            call check_intersecion_blade(point_onplane, planes(k, :, :),   check_inhub)
          endif


          ! If the point is inside
          if(check_inhub) then

            ! Finding point on the plane where to stop the particle (projection of the point on the plane)
            call projection_point_plane(ending_position, planes(k, :, :))

            ! Changing logical variable to true
            check = .true.

            ! Updating distance travelled
            distance = norm2(ending_position - starting_position)

            ! Exit from the cycle on all the rotors
            go to 500

          endif
        endif

      enddo

    endif
  enddo

  500 continue


  endsubroutine check_hub

  subroutine planes_gen (vertexes, planes)

    real(kind=8), dimension(:, :, :), intent(in)  ::  vertexes
    real(kind=8), dimension(:, :, :), intent(out) ::  planes

    planes(1, :, :) = vertexes(1, :, :)
    planes(2, :, :) = vertexes(2, :, :)

    planes(3, 1, :) = vertexes(1, 3, :);   planes(3, 2, :) = vertexes(1, 2, :);
    planes(3, 3, :) = vertexes(2, 4, :);   planes(3, 4, :) = vertexes(2, 3, :);

    planes(4, 1, :) = vertexes(1, 4, :);   planes(4, 2, :) = vertexes(1, 3, :);
    planes(4, 3, :) = vertexes(2, 3, :);   planes(4, 4, :) = vertexes(2, 2, :);

    planes(5, 1, :) = vertexes(1, 1, :);   planes(5, 2, :) = vertexes(1, 4, :);
    planes(5, 3, :) = vertexes(2, 2, :);   planes(5, 4, :) = vertexes(2, 1, :);

    planes(6, 1, :) = vertexes(1, 2, :);   planes(6, 2, :) = vertexes(1, 1, :);
    planes(6, 3, :) = vertexes(2, 1, :);   planes(6, 4, :) = vertexes(2, 4, :);

  endsubroutine planes_gen

endmodule mod_check
