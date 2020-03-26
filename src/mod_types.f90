module mod_types
   use mod_initialize

   implicit none

contains

   subroutine blade_panel_constructor(gamma, radius, collocation, vertex, dFz, dQ, alpha_old,   this_blade_panel)

      real(kind=8), dimension(:,:), intent(in) :: vertex
      real(kind=8), dimension(:), intent(in) :: collocation
      real(kind=8), intent(in) :: gamma, radius, dFz, dQ, alpha_old
      type(blade_panel), intent(out) :: this_blade_panel

      this_blade_panel%gamma = gamma
      this_blade_panel%radius = radius
      this_blade_panel%collocation = collocation
      this_blade_panel%vertex = vertex
      this_blade_panel%dFz = dFz
      this_blade_panel%dQ = dQ
      this_blade_panel%alpha_old = alpha_old

   end subroutine blade_panel_constructor

   subroutine blade_panel_deconstructor(this_blade_panel,  gamma, radius, collocation, vertex, dFz, dQ, alpha_old)

      type(blade_panel), intent(in):: this_blade_panel
      real(kind=8), dimension(:,:), intent(out) :: vertex
      real(kind=8), dimension(:), intent(out) :: collocation
      real(kind=8), intent(out) :: gamma, radius, dFz, dQ, alpha_old

      gamma = this_blade_panel%gamma
      radius = this_blade_panel%radius
      collocation = this_blade_panel%collocation
      vertex = this_blade_panel%vertex
      dFz = this_blade_panel%dFz
      dQ = this_blade_panel%dQ
      alpha_old = this_blade_panel%alpha_old

   end subroutine blade_panel_deconstructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine particle_constructor(alpha_p, position, vp, vp_old, radius, radius_old, on_ground,   this_particle)

      real(kind=8), dimension(:), intent(in) :: alpha_p, position, vp, vp_old
      real(kind=8), intent(in) :: radius, radius_old
      logical, intent(in) :: on_ground
      type(particle), intent(out) :: this_particle

      this_particle%alpha_p = alpha_p
      this_particle%position = position
      this_particle%vp = vp
      this_particle%vp_old = vp_old
      this_particle%radius = radius
      this_particle%radius_old = radius_old
      this_particle%on_ground = on_ground

   end subroutine particle_constructor

   subroutine particle_deconstructor(this_particle,   alpha_p, position, radius, radius_old, on_ground, vp, vp_old)

      type(particle), intent(in) :: this_particle

      real(kind=8), dimension(:), intent(out) :: alpha_p, position
      real(kind=8), dimension(:), intent(out), optional ::  vp, vp_old
      logical, intent(out) :: on_ground
      real(kind=8), intent(out) :: radius, radius_old

      alpha_p = this_particle%alpha_p
      position = this_particle%position
      if(present(vp)) then
        vp = this_particle%vp
        vp_old = this_particle%vp_old
      endif
      radius = this_particle%radius
      radius_old = this_particle%radius_old
      on_ground = this_particle%on_ground

   end subroutine particle_deconstructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine wake_panel_constructor(gamma, vertex, vp, vp_old, radius,   this_wake_panel)

      real(kind=8), intent(in) :: gamma  ! Circulation in the panel
      real(kind=8), dimension(:, :), intent(in) :: vertex
      real(kind=8), dimension(:, :), intent(in) :: vp, vp_old
      real(kind=8), intent(in) :: radius ! Radius of the "vortex"

      type(wake_panel), intent(out) :: this_wake_panel

      this_wake_panel%gamma = gamma
      this_wake_panel%vertex = vertex
      this_wake_panel%vp = vp
      this_wake_panel%vp_old = vp_old
      this_wake_panel%radius = radius

   end subroutine wake_panel_constructor

   subroutine wake_panel_deconstructor(this_wake_panel,   gamma, vertex, radius, vp, vp_old)
      type(wake_panel), intent(in) :: this_wake_panel

      real(kind=8), intent(out) :: gamma  ! Circulation in the panel
      real(kind=8), dimension(:,:), intent(out) :: vertex
      real(kind=8), dimension(:,:), intent(out), optional :: vp, vp_old
      real(kind=8), intent(out) :: radius ! Radius of the "vortex"

      gamma = this_wake_panel%gamma
      vertex = this_wake_panel%vertex
      if(present(vp)) then
        vp = this_wake_panel%vp
        if(present(vp_old)) then
          vp_old = this_wake_panel%vp_old
        endif
      endif
      radius = this_wake_panel%radius

   end subroutine wake_panel_deconstructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wake_segment_constructor(gamma, vertex, radius,   this_wake_segment)

   real(kind=8), intent(in) :: gamma  ! Circulation in the panel
   real(kind=8), dimension(:, :), intent(in) :: vertex
   real(kind=8), intent(in) :: radius ! Radius of the "vortex"

   type(segment), intent(out) :: this_wake_segment

   this_wake_segment%gamma = gamma
   this_wake_segment%vertex = vertex
   this_wake_segment%radius = radius

endsubroutine wake_segment_constructor


subroutine wake_segment_deconstructor(this_wake_segment,     gamma, vertex, radius)

   type(segment), intent(in) :: this_wake_segment

   real(kind=8), intent(out) :: gamma  ! Circulation in the panel
   real(kind=8), dimension(2, 3), intent(out) :: vertex
   real(kind=8), intent(out) :: radius ! Radius of the "vortex"

   gamma = this_wake_segment%gamma
   vertex = this_wake_segment%vertex
   radius = this_wake_segment%radius

endsubroutine wake_segment_deconstructor

endmodule mod_types
