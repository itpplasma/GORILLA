module test_supporting_functions_mod
   use supporting_functions_mod
   use funit
   implicit none
!
! Background:
! The module supporting_functions_mod is a collection of various functions/subroutines that
! have no specific association with other modules and are intended to be used from several
! places of the program. Two examples are the boole-to-integer converter logical2integer and 
! the eigenvalue routine calc_eigval(N,A,WR,VR). The latter provides not only the eigenvalues 
! as an array WR, but also the eigenvectors in form of a matrix VR.
!
! Georg Graßler (25.01.2023)
!      
contains

	@test
    ! Description of test_logical2integer_true:
    ! The tests test_logical2integer_true() checks if the boolean .TRUE. 
    ! is correctly converted to the integer 1.
    !
    ! Georg Graßler (25.01.2023)
    !
	subroutine test_logical2integer_true()
		logical :: boole_in
		integer :: value
		
		boole_in = .TRUE.
		value = logical2integer(boole_in)
		
		@assertEqual(value,1,message = "Testing logical2integer with true value")
		
	end subroutine test_logical2integer_true
	
	@test
    ! Description of test_logical2integer_false:
    ! The tests test_logical2integer_false() checks if the boolean .FALSE. 
    ! is correctly converted to the integer 0. .FALSE. -> 0). 
    !
    ! Georg Graßler (25.01.2023)
    !
	subroutine test_logical2integer_false()
		logical :: boole_in
		integer :: value
		
		boole_in = .FALSE.
		value = logical2integer(boole_in)
		
		@assertEqual(value,0,message = "Testing logical2integer with false value")
		
	end subroutine test_logical2integer_false
	
	@test
    ! Description of test_calc_eigval:
    ! In test_calc_eigval() a trial (3 x 3) matrix A is provided, of which the eigenvalues/eigenvectors 
    ! are known analytically. Not only are these exact results compared to the output of calc_eigval, 
    ! but it is also checked if the numerical found eigenvectors/eigenvalues selfconsistently fullfil 
    ! the eigenvalue equation
    !
    ! matmul(A,VR(:,i)) = WR(i) * VR(:,i) // Matrix x eigenvector = eigenvector scaled by eigenvalue
    !
    ! Georg Graßler (25.01.2023)
    !
	subroutine test_calc_eigval()
		
		integer :: N = 3
		double precision,dimension(3,3) :: A, VR
		double precision,dimension(3) :: WR
		double precision,dimension(3) :: eigvals
		double precision,dimension(3) :: eigvector1, eigvector2, eigvector3
				
		A = reshape((/3.d0, 2.d0, -2.d0, -1.d0, 0.d0, 2.d0, 0.d0, 0.d0, -1.d0/), (/N, N/))
		
		call calc_eigval(N,A,WR,VR)
	  	
	  	eigvals = [-1.d0,1.d0,2.d0]
	  	@assertEqual(WR, eigvals, tolerance=1e-13, message="Test calc_eigval, compare eigenvalues")
	  	
	  	eigvector1 = [0.d0,0.d0,1.d0]
	  	eigvector2 = [-1.d0,-2.d0,-1.d0] / SQRT(6.d0)
	  	eigvector3 = [-1.d0,-1.d0,0.d0] / SQRT(2.d0)	  	
	  	
	  	@assertEqual(VR(:,1), eigvector1, tolerance=1e-13, message="Test calc_eigval, compare first eigenvector") 
	  	@assertEqual(VR(:,2), eigvector2, tolerance=1e-13, message="Test calc_eigval, compare second eigenvector")
	  	@assertEqual(VR(:,3), eigvector3, tolerance=1e-13, message="Test calc_eigval, compare third eigenvector")
	  	
	  	!A needs to get declared again, because it gets overwritten in the calc_eigval subroutine
	  	A = reshape((/3.d0, 2.d0, -2.d0, -1.d0, 0.d0, 2.d0, 0.d0, 0.d0, -1.d0/), (/N, N/))
	  	
	  	@assertEqual(matmul(A, VR(:,1)), VR(:,1)*WR(1), tolerance=1e-13, message="Test calc_eigval, use eigenvalue equation")
	  	@assertEqual(matmul(A, VR(:,2)), VR(:,2)*WR(2), tolerance=1e-13, message="Test calc_eigval, use eigenvalue equation")
	  	@assertEqual(matmul(A, VR(:,3)), VR(:,3)*WR(3), tolerance=1e-13, message="Test calc_eigval, use eigenvalue equation")
		
	end subroutine test_calc_eigval

end module test_supporting_functions_mod
