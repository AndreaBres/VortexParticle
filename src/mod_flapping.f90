MODULE mod_flapping

use mod_initialize

contains

subroutine flapping(beta_here, beta_dot, beta_ddot, Ma, rotor)

   real(kind=8), dimension(:)    :: beta_here
   real(kind=8), dimension(:, :) :: beta_dot
   real(kind=8), dimension(:, :) :: beta_ddot
   real(kind=8), dimension(:)    :: Ma
   integer, intent(in) :: rotor

   integer :: i

   ! calcolo derivate
   beta_ddot(:, 1) = - omega(rotor)**2 * (1 + e(rotor)*Sbeta(rotor)/Ibeta(rotor)) * beta_here + Ma/Ibeta(rotor)

   ! propagazione alla Adam-Bashforth
   beta_dot(:, 1) = beta_dot(:, 2) + (1.5 * beta_ddot(:, 1) - 0.5 * beta_ddot(:, 2)) * dt
   beta_here = beta_here + (1.5 * beta_dot(:, 1) - 0.5 * beta_dot(:, 2)) * dt

   do i = 1, size(beta_here)
     do while(beta_here(i) > 2*pi)
       beta_here(i) = beta_here(i) - 2 * pi
     enddo
   enddo

   beta_dot(:, 2) = beta_dot(:, 1)
   beta_ddot(:, 2) = beta_ddot(:, 1)

end subroutine flapping

END MODULE mod_flapping
