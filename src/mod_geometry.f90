module mod_geometry

   use mod_initialize

   implicit none

   interface cross_product

     module procedure cross_product_s, cross_product_v

   endinterface cross_product

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine collocation_points(rotor, y_vec)
! SI POTREBBE CHIAMARE d LA y PER NON CONFONDERLA CON UNA COMPONENTE DEL VETTORE
! POSIZIONE
! bb sta per blade begin, il punto ADIMENSIONALE in cui comincia la pala
! costruisce il vettore dimensionale delle cooridnate dei punti
! di collocazione sulla pala
   implicit none

    integer, intent(in) :: rotor
    real(kind=8), dimension(:), intent(out) :: y_vec

   integer :: i
   real(kind=8) :: dy


! definisco il passo di y
   dy = R(rotor)*(1 - bs(rotor))/n_elements
! costruisco il vettore dimensionale delle y
   do  i = 1, n_elements
      y_vec(i) =  R(rotor)*bs(rotor) + dy/2 + (i - 1)*dy
   end do

end subroutine collocation_points


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine blade_matrix (rotor,    bp_vertex)
! Costruisce la matrice V che raccoglie le coordinate bidimensionali x e y (terzo indice)
! dei vertici dei pannelli. I primi due indici numerano i vertici come:
!  ----->y   1,1-----1,2-----1,3-----...-----1,n+1
!  '          '       ,       ,       ,       '
!  '          '       ,       ,       ,       '
!  V          '       ,       ,       ,       '
!  x         2,1-----2,2-----2,3-----...-----2,n+1

  integer, intent(in) :: rotor
  real(kind=8), dimension(2, n_elements+1, 2), intent(out) :: bp_vertex

  real(kind=8) :: dy
  integer :: j

  dy = R(rotor)*(1 - bs(rotor))/n_elements


     do  j = 1, n_elements + 1
        bp_vertex(:,j,2) = R(rotor)*bs(rotor)  + (j - 1)*dy
        bp_vertex(1,j,1) = 0.d0
        bp_vertex(2,j,1) = c(rotor)
     end do


endsubroutine blade_matrix

subroutine blade_vertex(psi, beta_here,  blades_vertexes, center_coord, rotor)

  integer, intent(in) :: rotor
  real(kind=8), dimension(:),       intent(in)  ::  psi, beta_here, center_coord
  real(kind=8), dimension(:, :, :), intent(out) ::  blades_vertexes

  real(kind=8), dimension(2, n_elements+1, 2) :: bp_vertex
  real(kind=8), dimension(2, 2, 2) :: b_vertex
  integer :: i

  ! Construction of the coordinates of the vertexes of the blade panels
  call blade_matrix (rotor,    bp_vertex)

  ! Tranlsation of the vertex from the panel one to the blade one
  do i = 1, size(bp_vertex, 2)
    bp_vertex(1, i, 1) = bp_vertex(1, i, 1) - c(rotor)/4
    bp_vertex(2, i, 1) = bp_vertex(2, i, 1) - c(rotor)/4
  enddo

  ! Acquisition of the 4 vertexes of the blade
  b_vertex(1, 1, :) = bp_vertex(1, 1, :);    b_vertex(1, 2, :) = bp_vertex(1, n_elements+1, :);
  b_vertex(2, 1, :) = bp_vertex(2, 1, :);    b_vertex(2, 2, :) = bp_vertex(2, n_elements+1, :);

  ! Rotation of the blade vertexes around psi
  do i = 1, size(psi, 1)

    blades_vertexes(i, 1, :) = global_blade_vector(psi(i), beta_here(i), b_vertex(1, 1, :), rotor) + center_coord
    blades_vertexes(i, 2, :) = global_blade_vector(psi(i), beta_here(i), b_vertex(1, 2, :), rotor) + center_coord
    blades_vertexes(i, 3, :) = global_blade_vector(psi(i), beta_here(i), b_vertex(2, 2, :), rotor) + center_coord
    blades_vertexes(i, 4, :) = global_blade_vector(psi(i), beta_here(i), b_vertex(2, 1, :), rotor) + center_coord

  enddo

endsubroutine blade_vertex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function global_blade_vector(psi, beta_here, p_blade, rotor) result(p_global)
! esprime in coordinate globali i punti p_blade in coodinate pala, espressi tramite
! le corrispettive coordinate x e y
!
!____________________________________________________
!------>X
!|
!|
!V___________________________________________________
!Y

   real(kind=8), intent(in) :: psi
   real(kind=8), intent(in) :: beta_here
   real(kind=8), dimension(:), intent(in)  :: p_blade
   integer,                    intent(in)  :: rotor
   real(kind=8)                            :: th
   real(kind=8), dimension(3)              :: p_global

th = theta(rotor) - theta_s(rotor)  * sin(psi)

p_global(1) = cos(alpha_sh(rotor))*cos(psi)*(e(rotor) + p_blade(2)*cos(beta_here) - &
              sin(beta_here)*sin(th)*(e(rotor) - p_blade(1))) &
            - sin(alpha_sh(rotor))*(p_blade(2)*sin(beta_here) + &
              cos(beta_here)*sin(th)*(e(rotor) - p_blade(1)))  &
            - cos(alpha_sh(rotor))*cos(th)*sin(psi)*(e(rotor) - p_blade(1))

p_global(2) = sin(psi)*(e(rotor) + p_blade(2)*cos(beta_here) - &
              sin(beta_here)*sin(th)*(e(rotor) - p_blade(1))) &
            + cos(psi)*cos(th)*(e(rotor) - p_blade(1))

p_global(3) = cos(alpha_sh(rotor))*(p_blade(2)*sin(beta_here) + &
              cos(beta_here)*sin(th)*(e(rotor) - p_blade(1))) &
            + cos(psi)*sin(alpha_sh(rotor))*(e(rotor) + p_blade(2)*cos(beta_here) - &
              sin(beta_here)*sin(th)*(e(rotor) - p_blade(1))) &
            - sin(alpha_sh(rotor))*cos(th)*sin(psi)*(e(rotor) - p_blade(1))

end function global_blade_vector


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine xyz2blade(psi, beta_here, p, rotor,   q)
! come da titolo cambia la base con cui viene espresso un vettore portandolo in
! sistema di riferimento blade

! Vectorial function

   implicit none

   integer, intent(in) :: rotor
   real(kind=8), dimension(:), intent(in) :: p
   real(kind=8),               intent(in) :: psi
   real(kind=8),               intent(in) :: beta_here
   real(kind=8)                           :: th
   real(kind=8), dimension(size(p)), intent(out) :: q

th = theta(rotor) - theta_s(rotor) * sin(psi)

! script di verifica su Matlab\Octave
q(1) = p(1)*(cos(alpha_sh(rotor))*(cos(th)*sin(psi) + cos(psi)*sin(beta_here)*sin(th)) +&
       cos(beta_here)*sin(alpha_sh(rotor))*sin(th)) - &
       p(2)*(cos(psi)*cos(th) - sin(beta_here)*sin(psi)*sin(th)) + &
       p(3)*(sin(alpha_sh(rotor))*(cos(th)*sin(psi) + &
       cos(psi)*sin(beta_here)*sin(th)) - cos(alpha_sh(rotor))*cos(beta_here)*sin(th))

q(2) = p(3)*(cos(alpha_sh(rotor))*sin(beta_here) + &
       cos(beta_here)*cos(psi)*sin(alpha_sh(rotor))) - &
       p(1)*(sin(alpha_sh(rotor))*sin(beta_here) - &
       cos(alpha_sh(rotor))*cos(beta_here)*cos(psi)) + p(2)*cos(beta_here)*sin(psi)

q(3) = p(1)*(cos(alpha_sh(rotor))*(sin(psi)*sin(th) - cos(psi)*sin(beta_here)*cos(th)) -&
       cos(beta_here)*sin(alpha_sh(rotor))*cos(th)) - &
       p(2)*(cos(psi)*sin(th) + sin(beta_here)*cos(th)*sin(psi)) + &
       p(3)*(sin(alpha_sh(rotor))*(sin(psi)*sin(th) - &
       cos(psi)*sin(beta_here)*cos(th)) + cos(alpha_sh(rotor))*cos(beta_here)*cos(th))


endsubroutine xyz2blade


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine xyz2alpha(psi, beta_here, p, rotor,   q)
! come da titolo cambia la base con cui viene espresso un vettore portandolo in
! sistema di riferimento blade

! Vectorial function

   implicit none

   integer, intent(in) :: rotor
   real(kind=8), dimension(:), intent(in) :: p
   real(kind=8),               intent(in) :: psi
   real(kind=8),               intent(in) :: beta_here
   real(kind=8)                           :: th
   real(kind=8), dimension(size(p)), intent(out) :: q

th = 0

! script di verifica su Matlab\Octave
q(1) = p(1)*(cos(alpha_sh(rotor))*(cos(th)*sin(psi) + cos(psi)*sin(beta_here)*sin(th)) +&
       cos(beta_here)*sin(alpha_sh(rotor))*sin(th)) - &
       p(2)*(cos(psi)*cos(th) - sin(beta_here)*sin(psi)*sin(th)) + &
       p(3)*(sin(alpha_sh(rotor))*(cos(th)*sin(psi) + &
       cos(psi)*sin(beta_here)*sin(th)) - cos(alpha_sh(rotor))*cos(beta_here)*sin(th))

q(2) = p(3)*(cos(alpha_sh(rotor))*sin(beta_here) + &
       cos(beta_here)*cos(psi)*sin(alpha_sh(rotor))) - &
       p(1)*(sin(alpha_sh(rotor))*sin(beta_here) - &
       cos(alpha_sh(rotor))*cos(beta_here)*cos(psi)) + p(2)*cos(beta_here)*sin(psi)

q(3) = p(1)*(cos(alpha_sh(rotor))*(sin(psi)*sin(th) - cos(psi)*sin(beta_here)*cos(th)) -&
       cos(beta_here)*sin(alpha_sh(rotor))*cos(th)) - &
       p(2)*(cos(psi)*sin(th) + sin(beta_here)*cos(th)*sin(psi)) + &
       p(3)*(sin(alpha_sh(rotor))*(sin(psi)*sin(th) - &
       cos(psi)*sin(beta_here)*cos(th)) + cos(alpha_sh(rotor))*cos(beta_here)*cos(th))


endsubroutine xyz2alpha


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cross_product_s(vec1, vec2,    res)
  real(kind=8), dimension(:), intent(in) :: vec1, vec2
  real(kind=8), dimension(:), intent(inout) :: res

  res(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
  res(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
  res(3) = vec1(1)*vec2(2) - vec2(1)*vec1(2)

end subroutine cross_product_s

subroutine cross_product_v(vec1, vec2,    res)
  real(kind=8), dimension(:, :), intent(in)     ::  vec1, vec2
  real(kind=8), dimension(:, :), intent(inout)  ::  res

  res(:, 1) = vec1(:, 2)*vec2(:, 3) - vec1(:, 3)*vec2(:, 2)
  res(:, 2) = vec1(:, 3)*vec2(:, 1) - vec1(:, 1)*vec2(:, 3)
  res(:, 3) = vec1(:, 1)*vec2(:, 2) - vec2(:, 1)*vec1(:, 2)

end subroutine cross_product_v

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function blade_panel_vertex (v1, v2, v3, v4, psi, beta_here, rotor) result(vertex)

   real(kind=8), dimension(:)   :: v1, v2, v3, v4
   integer                      :: rotor
   real(kind=8)                 :: psi
   real(kind=8)                 :: beta_here
   real(kind=8), dimension(4,3) :: vertex

   vertex(1,:) = global_blade_vector(psi, beta_here, v1, rotor)
   vertex(2,:) = global_blade_vector(psi, beta_here, v2, rotor)
   vertex(3,:) = global_blade_vector(psi, beta_here, v3, rotor)
   vertex(4,:) = global_blade_vector(psi, beta_here, v4, rotor)

end function blade_panel_vertex

function area(panel) result(S)

   implicit none

   type(blade_panel), intent(in) :: panel
   real(kind=8) :: S

   real(kind=8), dimension(3) :: v1, v2, v3

   v1 = panel%vertex(2,:) - panel%vertex(1,:)
   v2 = panel%vertex(3,:) - panel%vertex(2,:)

   call cross_product(v1, v2,  v3)

   S = norm2(v3)

end function area

end module mod_geometry
