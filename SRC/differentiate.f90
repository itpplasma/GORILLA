!
  subroutine differentiate(x,y,z,n,f,fx,fy,fz)
!
! Computes derivatives of linear functions specified at tetrahedron vertices
!
! Input parameters:
!           Formal: x(4),y(4),z(4)    - coordinates (x,y,z) of 4 vertices
!                   n                 - number of functions to differentiate
!                   f(4,n)            - values of functions at the vertices
! Output parameters:
!           Formal: fx(n),fy(n),fz(n) - deivatives of functions over x,y and z, resp.
!
! Called routines:  dgesv             - from the Lapack library
!
  implicit none
!
  integer                          :: n
  double precision, dimension(n)   :: fx,fy,fz
  double precision, dimension(4)   :: x,y,z
  double precision, dimension(4,n) :: f
!
  integer                          :: ierr
  integer, dimension(3)            :: ipiv
  double precision, dimension(3,3) :: a,b
  double precision, dimension(:,:), allocatable :: df
!
  allocate(df(3,n))
!
  df(1,:)=f(2,:)-f(1,:)
  df(2,:)=f(3,:)-f(1,:)
  df(3,:)=f(4,:)-f(1,:)
!
  a(1,1)=x(2)-x(1)
  a(1,2)=x(3)-x(1)
  a(1,3)=x(4)-x(1)
  a(2,1)=y(2)-y(1)
  a(2,2)=y(3)-y(1)
  a(2,3)=y(4)-y(1)
  a(3,1)=z(2)-z(1)
  a(3,2)=z(3)-z(1)
  a(3,3)=z(4)-z(1)
!
  b=0.d0
  b(1,1)=1.d0
  b(2,2)=1.d0
  b(3,3)=1.d0
!
  call dgesv(3,3,a,3,ipiv,b,3,ierr)
!
  df=matmul(transpose(b),df)
!
  fx=df(1,:)
  fy=df(2,:)
  fz=df(3,:)
!
  deallocate(df)
!
  end subroutine differentiate
