program test_jorek_find_tetra_continuity
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use find_tetra_mod, only: find_tetra
    use gorilla_settings_mod, only: load_gorilla_inp
    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use tetra_grid_mod, only: ntetr, tetra_grid, verts_rphiz
    use tetra_grid_settings_mod, only: load_tetra_grid_inp
    use tetra_physics_mod, only: tetra_physics

    implicit none

    integer, parameter :: target_faces = 512
    real(dp), parameter :: inward_fraction = 1.0e-6_dp
    type :: continuity_metrics_t
        integer :: samples = 0
        integer :: missing = 0
        integer :: wrong_owner = 0
        integer :: constant_phi_faces = 0
        integer :: spanning_phi_faces = 0
        real(dp) :: max_face_jump = 0.0_dp
        real(dp) :: max_near_jump = 0.0_dp
    end type continuity_metrics_t

    call load_tetra_grid_inp()
    call load_gorilla_inp()
    call initialize_gorilla()
    call verify_continuity()

contains

    subroutine verify_continuity()
        type(continuity_metrics_t) :: metrics
        integer :: candidate_faces, stride

        candidate_faces = count_candidate_faces()
        if (candidate_faces < target_faces) &
            error stop 'too few production JOREK faces for continuity sample'
        stride = max(1, candidate_faces/target_faces)
        call sample_shared_faces(stride, metrics)
        print '(A, I0, A, I0, A, I0)', 'sampled=', metrics%samples, &
            ' missing=', metrics%missing, ' wrong owner=', metrics%wrong_owner
        print '(A, I0, A, I0)', &
            'constant-phi faces=', metrics%constant_phi_faces, &
            ' phi-spanning faces=', metrics%spanning_phi_faces
        print '(A, ES14.6)', 'maximum face-limit relative B jump=', &
            metrics%max_face_jump
        print '(A, ES14.6)', 'maximum near-face relative B jump=', &
            metrics%max_near_jump
        if (metrics%samples /= target_faces .or. metrics%missing /= 0 &
                .or. metrics%wrong_owner /= 0) &
            error stop 'production JOREK face ownership is discontinuous'
        if (metrics%constant_phi_faces == 0 &
                .or. metrics%spanning_phi_faces == 0) &
            error stop 'production JOREK face sample lost an orientation class'
        if (.not. ieee_is_finite(metrics%max_face_jump) &
                .or. .not. ieee_is_finite(metrics%max_near_jump)) &
            error stop 'production JOREK continuity metric is non-finite'
        if (metrics%max_face_jump > 1.0e-12_dp &
                .or. metrics%max_near_jump > 1.0e-7_dp) &
            error stop 'production JOREK field is discontinuous across a shared face'
        print '(A)', 'PASS: production JOREK find_tetra face continuity'
    end subroutine verify_continuity

    integer function count_candidate_faces()
        integer :: face, tetra

        count_candidate_faces = 0
        do tetra = 1, ntetr
            do face = 1, 4
                if (is_candidate_face(tetra, face)) &
                    count_candidate_faces = count_candidate_faces + 1
            end do
        end do
    end function count_candidate_faces

    subroutine sample_shared_faces(stride, metrics)
        integer, intent(in) :: stride
        type(continuity_metrics_t), intent(inout) :: metrics

        integer :: candidate, face, neighbour, tetra

        candidate = 0
        sample_loop: do tetra = 1, ntetr
            do face = 1, 4
                if (.not. is_candidate_face(tetra, face)) cycle
                candidate = candidate + 1
                if (mod(candidate - 1, stride) /= 0) cycle
                neighbour = tetra_grid(tetra)%neighbour_tetr(face)
                call sample_face(tetra, face, neighbour, metrics)
                if (metrics%samples == target_faces) exit sample_loop
            end do
        end do sample_loop
    end subroutine sample_shared_faces

    logical function is_candidate_face(tetra_index, face_index)
        integer, intent(in) :: tetra_index, face_index

        is_candidate_face = &
            tetra_grid(tetra_index)%neighbour_tetr(face_index) > tetra_index &
            .and. tetra_grid(tetra_index)%neighbour_face(face_index) > 0 &
            .and. tetra_grid(tetra_index)%neighbour_perbou_phi(face_index) == 0 &
            .and. tetra_phi_span(tetra_index) <= acos(-1.0_dp)
    end function is_candidate_face

    subroutine sample_face(tetra_index, face_index, neighbour, metrics)
        integer, intent(in) :: tetra_index, face_index, neighbour
        type(continuity_metrics_t), intent(inout) :: metrics

        real(dp) :: face_field_left(3), face_field_right(3), face_point(3)
        real(dp) :: left_center(3), left_field(3), left_point(3)
        real(dp) :: right_center(3), right_field(3), right_point(3)
        integer :: face_nodes(3), iface, left_found, right_found

        face_nodes = pack([1, 2, 3, 4], [1, 2, 3, 4] /= face_index)
        face_point = sum(verts_rphiz(:, &
            tetra_grid(tetra_index)%ind_knot(face_nodes)), dim=2)/3.0_dp
        call count_face_orientation(tetra_index, face_nodes, metrics)
        call tetra_center(tetra_index, left_center)
        call tetra_center(neighbour, right_center)
        left_point = (1.0_dp - inward_fraction)*face_point &
            + inward_fraction*left_center
        right_point = (1.0_dp - inward_fraction)*face_point &
            + inward_fraction*right_center
        call find_tetra(left_point, 1.0e6_dp, 1.0e6_dp, left_found, iface)
        call find_tetra(right_point, 1.0e6_dp, 1.0e6_dp, right_found, iface)
        if (left_found < 1 .or. right_found < 1) then
            metrics%missing = metrics%missing + 1
            return
        end if
        if (left_found /= tetra_index .or. right_found /= neighbour) &
            metrics%wrong_owner = metrics%wrong_owner + 1
        call field_at(tetra_index, face_point, face_field_left)
        call field_at(neighbour, face_point, face_field_right)
        call field_at(left_found, left_point, left_field)
        call field_at(right_found, right_point, right_field)
        metrics%max_face_jump = max(metrics%max_face_jump, &
            relative_difference(face_field_left, face_field_right))
        metrics%max_near_jump = max(metrics%max_near_jump, &
            relative_difference(left_field, right_field))
        metrics%samples = metrics%samples + 1
    end subroutine sample_face

    subroutine count_face_orientation(tetra_index, face_nodes, metrics)
        integer, intent(in) :: tetra_index, face_nodes(3)
        type(continuity_metrics_t), intent(inout) :: metrics

        real(dp) :: phi(3)

        phi = verts_rphiz(2, tetra_grid(tetra_index)%ind_knot(face_nodes))
        if (maxval(phi) - minval(phi) < 1.0e-12_dp) then
            metrics%constant_phi_faces = metrics%constant_phi_faces + 1
        else
            metrics%spanning_phi_faces = metrics%spanning_phi_faces + 1
        end if
    end subroutine count_face_orientation

    real(dp) function tetra_phi_span(tetra_index)
        integer, intent(in) :: tetra_index

        real(dp) :: values(4)

        values = verts_rphiz(2, tetra_grid(tetra_index)%ind_knot)
        tetra_phi_span = maxval(values) - minval(values)
    end function tetra_phi_span

    subroutine tetra_center(tetra_index, center)
        integer, intent(in) :: tetra_index
        real(dp), intent(out) :: center(3)

        real(dp) :: coordinates(3, 4)

        coordinates = verts_rphiz(:, tetra_grid(tetra_index)%ind_knot)
        if (maxval(coordinates(2, :)) - minval(coordinates(2, :)) &
                > acos(-1.0_dp)) then
            where (coordinates(2, :) < acos(-1.0_dp))
                coordinates(2, :) = coordinates(2, :) + 2.0_dp*acos(-1.0_dp)
            end where
        end if
        center = sum(coordinates, dim=2)/4.0_dp
    end subroutine tetra_center

    subroutine field_at(tetra_index, point, field)
        integer, intent(in) :: tetra_index
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: field(3)

        real(dp) :: bmod, h(3), offset(3)

        offset = point - tetra_physics(tetra_index)%x1
        bmod = tetra_physics(tetra_index)%bmod1 &
            + sum(tetra_physics(tetra_index)%gB*offset)
        h(1) = tetra_physics(tetra_index)%h1_1 &
            + sum(tetra_physics(tetra_index)%gh1*offset)
        h(2) = tetra_physics(tetra_index)%h2_1 &
            + sum(tetra_physics(tetra_index)%gh2*offset)
        h(3) = tetra_physics(tetra_index)%h3_1 &
            + sum(tetra_physics(tetra_index)%gh3*offset)
        field = [h(1)*bmod, h(2)*bmod/point(1), h(3)*bmod]
    end subroutine field_at

    real(dp) function relative_difference(left, right)
        real(dp), intent(in) :: left(3), right(3)

        relative_difference = sqrt(sum((left - right)**2)) &
            /max(sqrt(sum(left**2)), sqrt(sum(right**2)), tiny(1.0_dp))
    end function relative_difference

end program test_jorek_find_tetra_continuity
