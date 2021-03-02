!
module pusher_tetra_func_mod
!
contains
!
    subroutine pusher_handover2neighbour(ind_tetr,ind_tetr_out,iface_inout,x,iper_phi)
!    
      use constants, only: pi
      use tetra_physics_mod, only: coord_system, tetra_skew_coord
      use tetra_grid_mod, only: tetra_grid
      use tetra_grid_settings_mod, only: n_field_periods
      use gorilla_settings_mod, only: handover_processing_kind
!      
      implicit none
!
      integer, intent(in)       :: ind_tetr
      integer, intent(inout)    :: iface_inout
      integer, intent(out)      :: ind_tetr_out
      double precision, dimension(3), intent(inout) :: x
      double precision, dimension(3) :: x_lin_cart, b
      integer, intent(out) :: iper_phi
      integer :: iface,iper_theta
!
      !Save input variables
      iface = iface_inout
!     
      !Set output quantities through logic of the mesh
      ind_tetr_out=tetra_grid(ind_tetr)%neighbour_tetr(iface)
      iface_inout=tetra_grid(ind_tetr)%neighbour_face(iface)
      iper_phi=tetra_grid(ind_tetr)%neighbour_perbou_phi(iface)
      iper_theta=tetra_grid(ind_tetr)%neighbour_perbou_theta(iface)
!
      !Position exchange in between tetrahedra
      !Distinguish in between different kinds of processing handling (handover_processing_kind)
      select case(handover_processing_kind)
        case(1) !Normal processing with treatment of periodic boundaries and manipulation of periodic values
!
          select case(coord_system)
            case(1) !R,phi,Z --> Cylindrical coordinate system
              if(iper_phi.eq.1) then
                x(2) = x(2) - 2.d0*pi/n_field_periods
              elseif(iper_phi.eq.-1) then
                x(2) = x(2) + 2.d0*pi/n_field_periods
              endif
!
            case(2) !s,theta,phi --> Symmetry flux coordinates
              if(iper_phi.eq.1) then
                x(3) = x(3) - 2.d0*pi/n_field_periods
              elseif(iper_phi.eq.-1) then
                x(3) = x(3) + 2.d0*pi/n_field_periods
              endif
              if(iper_theta.eq.1) then
                x(2) = x(2) - 2.d0*pi
              elseif(iper_theta.eq.-1) then
                x(2) = x(2) + 2.d0*pi
              endif
          end select
!        
        case(2) !Position exchange via Cartesian variables (skew_coordinates)
!          
          !Calculation with already pre-computed matrices and inverse matrices
!         
            !---------------------Compute position in 'leaving'-tetrahedron --------------------------------!
!
            !Vector with right handside of linear Equation set A*X=B, where A ... skew coordinates mat
            b(1) = x(1) - tetra_skew_coord(ind_tetr)%skew_ref_x1x2x3(1,iface)
            b(2) = x(2) - tetra_skew_coord(ind_tetr)%skew_ref_x1x2x3(2,iface)
            b(3) = x(3) - tetra_skew_coord(ind_tetr)%skew_ref_x1x2x3(3,iface)
!          
            !Linear transformation to Cartesian coordinates (!Not exact transformation!)
            x_lin_cart = matmul(tetra_skew_coord(ind_tetr)%skew_coord_xyz(:,:,iface), &
                              & matmul(tetra_skew_coord(ind_tetr)%inv_skew_coord_x1x2x3(:,:,iface),b) )
            x_lin_cart(1) = x_lin_cart(1) + tetra_skew_coord(ind_tetr)%skew_ref_xyz(1,iface)
            x_lin_cart(2) = x_lin_cart(2) + tetra_skew_coord(ind_tetr)%skew_ref_xyz(2,iface)
            x_lin_cart(3) = x_lin_cart(3) + tetra_skew_coord(ind_tetr)%skew_ref_xyz(3,iface)
!
            !-----------------Compute position in neighbouring 'entering'-tetrahedron ----------------------!
!
            !Matrix with right handside of linear Equation set A*X=B, where A ... skew coordinates mat
            b(1) = x_lin_cart(1) - tetra_skew_coord(ind_tetr_out)%skew_ref_xyz(1,iface_inout)
            b(2) = x_lin_cart(2) - tetra_skew_coord(ind_tetr_out)%skew_ref_xyz(2,iface_inout)
            b(3) = x_lin_cart(3) - tetra_skew_coord(ind_tetr_out)%skew_ref_xyz(3,iface_inout)
!
            !Linear transformation to calc-(coord_system) coordinates (!Not exact transformation!)
            x = matmul(tetra_skew_coord(ind_tetr_out)%skew_coord_x1x2x3(:,:,iface_inout), &
                              & matmul(tetra_skew_coord(ind_tetr_out)%inv_skew_coord_xyz(:,:,iface_inout),b) )
            x(1) = x(1) + tetra_skew_coord(ind_tetr_out)%skew_ref_x1x2x3(1,iface_inout)
            x(2) = x(2) + tetra_skew_coord(ind_tetr_out)%skew_ref_x1x2x3(2,iface_inout)
            x(3) = x(3) + tetra_skew_coord(ind_tetr_out)%skew_ref_x1x2x3(3,iface_inout)
!          
      end select
!
    end subroutine pusher_handover2neighbour
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module pusher_tetra_func_mod
