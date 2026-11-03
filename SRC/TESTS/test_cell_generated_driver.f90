program test_cell_generated_driver
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use gorilla_cell_generated, only: gorilla_affine_coefficients2, &
                                     gorilla_poly_evaluate2, gorilla_poly_step1, &
                                     gorilla_poly_step2, &
                                     gorilla_poly_step3, gorilla_poly_step4
    implicit none
    real(dp) :: amat(4,4), b(4), z(4), coef(4,3), expected(4), actual(4)
    real(dp) :: jets(4,7), tau, tau_values(5), derivative(4), weight
    integer :: i, j, k, order

    tau_values = [-1.25_dp, -0.01_dp, 0.0_dp, 0.03125_dp, 2.0_dp]
    ! A nilpotent acceleration system has an exact quadratic trajectory.
    amat = 0.0_dp
    amat(1,3) = 1.0_dp
    amat(2,4) = 1.0_dp
    z = [1.0_dp, -2.0_dp, 0.5_dp, -0.75_dp]
    b = [0.2_dp, -0.3_dp, 1.5_dp, -2.0_dp]
    call gorilla_affine_coefficients2(amat,b,z,coef)
    do k = 1, size(tau_values)
        tau = tau_values(k)
        expected(1:2) = z(1:2) + tau*(z(3:4)+b(1:2)) + 0.5_dp*tau**2*b(3:4)
        expected(3:4) = z(3:4) + tau*b(3:4)
        call gorilla_poly_evaluate2(coef,tau,actual)
        call check('nilpotent trajectory',actual,expected)
    end do

    ! Dense systems check matrix orientation independently of scalar emission.
    do j = 1, 4
        do i = 1, 4
            amat(i,j) = real(3*i-2*j,dp)/7.0_dp
        end do
    end do
    call gorilla_affine_coefficients2(amat,b,z,coef)
    call check('constant coefficient',coef(:,1),z)
    call check('linear coefficient',coef(:,2),matmul(amat,z)+b)
    expected = 0.5_dp*(matmul(matmul(amat,amat),z)+matmul(amat,b))
    call check('quadratic coefficient',coef(:,3),expected)

    ! Cached jets need not equal fresh products after pusher corrections.
    do j = 1, 7
        do i = 1, 4
            jets(i,j) = sin(real(2*i+3*j,dp))
        end do
    end do
    do k = 1, size(tau_values)
        tau = tau_values(k)
        do order = 1, 4
            expected = z + tau*(b+jets(:,1))
            if (order >= 2) expected = expected + tau**2*0.5_dp*(jets(:,2)+jets(:,3))
            if (order >= 3) expected = expected + tau**3/6.0_dp*(jets(:,4)+jets(:,5))
            if (order >= 4) expected = expected + tau**4/24.0_dp*(jets(:,6)+jets(:,7))
            select case(order)
            case(1)
                call gorilla_poly_step1(z,b,jets(:,1),tau,actual)
            case(2)
                call gorilla_poly_step2(z,b,jets(:,1),jets(:,2),jets(:,3),tau,actual)
            case(3)
                call gorilla_poly_step3(z,b,jets(:,1),jets(:,2),jets(:,3), &
                                       jets(:,4),jets(:,5),tau,actual)
            case(4)
                call gorilla_poly_step4(z,b,jets(:,1),jets(:,2),jets(:,3), &
                                       jets(:,4),jets(:,5),jets(:,6),jets(:,7), &
                                       tau,actual)
            end select
            call check('cached Taylor jet',actual,expected)
        end do
    end do

    jets(:,1) = matmul(amat,z)
    jets(:,2) = matmul(amat,b)
    jets(:,3) = matmul(amat,jets(:,1))
    jets(:,4) = matmul(amat,jets(:,2))
    jets(:,5) = matmul(amat,jets(:,3))
    jets(:,6) = matmul(amat,jets(:,4))
    jets(:,7) = matmul(amat,jets(:,5))
    do k = 1, size(tau_values)
        tau = tau_values(k)
        expected = z
        derivative = matmul(amat,z) + b
        weight = 1.0_dp
        do order = 1, 4
            weight = weight*tau/real(order,dp)
            expected = expected + weight*derivative
            select case(order)
            case(1)
                call gorilla_poly_step1(z,b,jets(:,1),tau,actual)
            case(2)
                call gorilla_poly_step2(z,b,jets(:,1),jets(:,2),jets(:,3),tau,actual)
            case(3)
                call gorilla_poly_step3(z,b,jets(:,1),jets(:,2),jets(:,3), &
                                       jets(:,4),jets(:,5),tau,actual)
            case(4)
                call gorilla_poly_step4(z,b,jets(:,1),jets(:,2),jets(:,3), &
                                       jets(:,4),jets(:,5),jets(:,6),jets(:,7), &
                                       tau,actual)
            end select
            call check('matrix Taylor series',actual,expected)
            derivative = matmul(amat,derivative)
        end do
    end do
    print *, 'generated cell arithmetic oracle: PASS'

contains

    subroutine check(label,actual,expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual(4), expected(4)
        if (any(abs(actual-expected) > &
                128.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(expected)))) then
            print *, label, actual, expected
            error stop 'generated arithmetic mismatch'
        end if
    end subroutine check

end program test_cell_generated_driver
