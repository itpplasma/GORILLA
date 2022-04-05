module test_supporting_functions_mod
   use supporting_functions_mod
   use funit
   implicit none
   
contains

	@test
	subroutine test_logical2integer_true()
		logical :: boole_in
		integer :: value
		
		boole_in = .TRUE.
		value = logical2integer(boole_in)
		
		@assertEqual(value,1,message = "Testing logical2integer with true value")
		
	end subroutine test_logical2integer_true
	
	@test
	subroutine test_logical2integer_false()
		logical :: boole_in
		integer :: value
		
		boole_in = .FALSE.
		value = logical2integer(boole_in)
		
		@assertEqual(value,0,message = "Testing logical2integer with false value")
		
	end subroutine test_logical2integer_false
	
	@test
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
