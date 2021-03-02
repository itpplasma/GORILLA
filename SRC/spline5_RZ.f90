!
  subroutine spl_five_reg(n,h,a,b,c,d,e,f)
!
  implicit none
!
  integer :: n,i,ip1,ip2
  double precision :: h,rhop,rhom,fac,fpl31,fpl40,fmn31,fmn40          ,x
  double precision :: a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,det
  double precision :: abeg,bbeg,cbeg,dbeg,ebeg,fbeg
  double precision :: aend,bend,cend,dend,eend,fend
  double precision, dimension(n) :: a,b,c,d,e,f
  double precision, dimension(:), allocatable :: alp,bet,gam
!
  rhop=13.d0+sqrt(105.d0)
  rhom=13.d0-sqrt(105.d0)
!
  a11=1.d0
  a12=1.d0/4.d0
  a13=1.d0/16.d0
  a21=3.d0
  a22=27.d0/4.d0
  a23=9.d0*27.d0/16.d0
  a31=5.d0
  a32=125.d0/4.d0
  a33=5.d0**5/16.d0
  det=a11*a22*a33+a12*a23*a31+a13*a21*a32-a12*a21*a33-a13*a22*a31-a11*a23*a32
  b1=a(4)-a(3)
  b2=a(5)-a(2)
  b3=a(6)-a(1)
  bbeg=b1*a22*a33+a12*a23*b3+a13*b2*a32-a12*b2*a33-a13*a22*b3-b1*a23*a32
  bbeg=bbeg/det
  dbeg=a11*b2*a33+b1*a23*a31+a13*a21*b3-b1*a21*a33-a13*b2*a31-a11*a23*b3
  dbeg=dbeg/det
  fbeg=a11*a22*b3+a12*b2*a31+b1*a21*a32-a12*a21*b3-b1*a22*a31-a11*b2*a32
  fbeg=fbeg/det
  b1=a(n-2)-a(n-3)
  b2=a(n-1)-a(n-4)
  b3=a(n)-a(n-5)
  bend=b1*a22*a33+a12*a23*b3+a13*b2*a32-a12*b2*a33-a13*a22*b3-b1*a23*a32
  bend=bend/det
  dend=a11*b2*a33+b1*a23*a31+a13*a21*b3-b1*a21*a33-a13*b2*a31-a11*a23*b3
  dend=dend/det
  fend=a11*a22*b3+a12*b2*a31+b1*a21*a32-a12*a21*b3-b1*a22*a31-a11*b2*a32
  fend=fend/det
  a11=2.d0
  a12=1.d0/2.d0
  a13=1.d0/8.d0
  a21=2.d0
  a22=9.d0/2.d0
  a23=81.d0/8.d0
  a31=2.d0
  a32=25.d0/2.d0
  a33=625.d0/8.d0
  det=a11*a22*a33+a12*a23*a31+a13*a21*a32-a12*a21*a33-a13*a22*a31-a11*a23*a32
  b1=a(4)+a(3)
  b2=a(5)+a(2)
  b3=a(6)+a(1)
  abeg=b1*a22*a33+a12*a23*b3+a13*b2*a32-a12*b2*a33-a13*a22*b3-b1*a23*a32
  abeg=abeg/det
  cbeg=a11*b2*a33+b1*a23*a31+a13*a21*b3-b1*a21*a33-a13*b2*a31-a11*a23*b3
  cbeg=cbeg/det
  ebeg=a11*a22*b3+a12*b2*a31+b1*a21*a32-a12*a21*b3-b1*a22*a31-a11*b2*a32
  ebeg=ebeg/det
  b1=a(n-2)+a(n-3)
  b2=a(n-1)+a(n-4)
  b3=a(n)+a(n-5)
  aend=b1*a22*a33+a12*a23*b3+a13*b2*a32-a12*b2*a33-a13*a22*b3-b1*a23*a32
  aend=aend/det
  cend=a11*b2*a33+b1*a23*a31+a13*a21*b3-b1*a21*a33-a13*b2*a31-a11*a23*b3
  cend=cend/det
  eend=a11*a22*b3+a12*b2*a31+b1*a21*a32-a12*a21*b3-b1*a22*a31-a11*b2*a32
  eend=eend/det
!
  allocate(alp(n),bet(n),gam(n))
!
  alp(1)=0.0d0
  bet(1)=ebeg*(2.d0+rhom)-5.d0*fbeg*(3.d0+1.5d0*rhom) !gamma1
!
  do i=1,n-4
    ip1=i+1
    alp(ip1)=-1.d0/(rhop+alp(i))
    bet(ip1)=alp(ip1)*(bet(i)- &
             5.d0*(a(i+4)-4.d0*a(i+3)+6.d0*a(i+2)-4.d0*a(ip1)+a(i)))
  enddo
!
  gam(n-2)=eend*(2.d0+rhom)+5.d0*fend*(3.d0+1.5d0*rhom) !gamma
  do i=n-3,1,-1
    gam(i)=gam(i+1)*alp(i)+bet(i)
  enddo
!
  alp(1)=0.0d0
  bet(1)=ebeg-2.5d0*5.d0*fbeg !e1
!
  do i=1,n-2
    ip1=i+1
    alp(ip1)=-1.d0/(rhom+alp(i))
    bet(ip1)=alp(ip1)*(bet(i)-gam(i))
  enddo
!
  e(n)=eend+2.5d0*5.d0*fend
  e(n-1)=e(n)*alp(n-1)+bet(n-1)
  f(n-1)=(e(n)-e(n-1))/5.d0
  e(n-2)=e(n-1)*alp(n-2)+bet(n-2)
  f(n-2)=(e(n-1)-e(n-2))/5.d0
  d(n-2)=dend+1.5d0*4.d0*eend+1.5d0**2*10.d0*fend
!
  do i=n-3,1,-1
    e(i)=e(i+1)*alp(i)+bet(i)
    f(i)=(e(i+1)-e(i))/5.d0
    d(i)=(a(i+3)-3.d0*a(i+2)+3.d0*a(i+1)-a(i))/6.d0 &
        -(e(i+3)+27.d0*e(i+2)+93.d0*e(i+1)+59.d0*e(i))/30.d0
    c(i)=0.5d0*(a(i+2)+a(i))-a(i+1)-0.5d0*d(i+1)-2.5d0*d(i) &
        -0.1d0*(e(i+2)+18.d0*e(i+1)+31.d0*e(i))
    b(i)=a(i+1)-a(i)-c(i)-d(i)-0.2d0*(4.d0*e(i)+e(i+1))
  enddo
!
  do i=n-3,n
    b(i)=b(i-1)+2.d0*c(i-1)+3.d0*d(i-1)+4.d0*e(i-1)+5.d0*f(i-1)
    c(i)=c(i-1)+3.d0*d(i-1)+6.d0*e(i-1)+10.d0*f(i-1)
    d(i)=d(i-1)+4.d0*e(i-1)+10.d0*f(i-1)
    if(i.ne.n) f(i)= a(i+1)-a(i)-b(i)-c(i)-d(i)-e(i)
  enddo
  f(n)=f(n-1)
!
  fac=1.d0/h
  b=b*fac
  fac=fac/h
  c=c*fac
  fac=fac/h
  d=d*fac
  fac=fac/h
  e=e*fac
  fac=fac/h
  f=f*fac
!
  deallocate(alp,bet,gam)
!
  return
  end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine s2dcut(nx,ny,hx,hy,f,imi,ima,jmi,jma,icount,spl,ipoint)
!
! Calculates coefficients of a 2D spline for a convex domain
! (for a non-rectangular domain the interpolation is not continious)
! equidistant mesh, but hx must not be = hy
! 
!  Input parameters:
!                    nx           - horizontal size of the mesh (over x)
!                    ny           - vertical size of the mesh (over y)
!                    hx           - step of the mesh over x
!                    hy           - step of the mesh over y
!                    f(i,j)       - array of values to be interpolated
!                                   (i = 1, ..., nx; j = 1, ..., ny).
!
!                                   For the case of non-rectangular domain:
!                                   numbers of the mesh points which
!                                   correspond to the boundaries of the
!                                   interpolation region:
!                    imi(j)       - left boundary of the row j (j=1,...,ny)
!                    ima(j)       - right boundary of the row j (j=1,...,ny)
!                    jmi(i)       - lower boundary of the column i (i=1,...,nx)
!                    jma(i)       - upper boundary of the column i (i=1,...,nx)
!                                   in a rectangle should be:
!                                   imi(:) = 1
!                                   ima(:) = nx
!                                   jmi(:) = 1
!                                   jma(:) = ny
!
!                    icount       - maximum number of entries in spl
!
! Output parameters:
!                    spl(l,m,k)   - spline coefficients (i,j = 1, ... , n;
!                    ipoint(i,j)    l,m = 1, ..., 4; i,j - numbers of the
!                                   mesh point in horizontal and vertical
!                                   direction (over x and over y), l,m -
!                                   the numbers of expansion power over x
!                                   and y (~ dx**(l-1)*dy**(m-1) ))
!                                   ipoint(i,j) contains the pointer to k
!
  implicit double precision (a-h,o-z)
!
  dimension f(nx,ny),spl(6,6,icount),ipoint(nx,ny)
  dimension imi(ny),ima(ny),jmi(nx),jma(nx)
!
  double precision, dimension(:), allocatable :: ai,bi,ci,di,ei,fi
!
  nmax=max(nx,ny)
!
  allocate( ai(nmax),bi(nmax),ci(nmax),di(nmax),ei(nmax),fi(nmax) )
!
  spl=0.d0
  ipoint=-1
!
!  spline along Y-axis
!
  ic = 0
  do i=1,nx
    if(jmi(i).gt.0) then
      nsi=jma(i)-jmi(i)+1
      do j=jmi(i),jma(i)
        ai(j-jmi(i)+1)=f(i,j)
      enddo
      call spl_five_reg(nsi,hy,ai,bi,ci,di,ei,fi)
      do j=jmi(i),jma(i)
        jj=j-jmi(i)+1
        ic = ic+1
        ipoint(i,j)=ic
        spl(1,1,ic)=ai(jj)
        spl(1,2,ic)=bi(jj)
        spl(1,3,ic)=ci(jj)
        spl(1,4,ic)=di(jj)
        spl(1,5,ic)=ei(jj)
        spl(1,6,ic)=fi(jj)
      enddo
    endif
  enddo
!
  if (ic .ne. icount) then
    write (6,*) 'Warning, ic, icount:  ',ic,icount
  endif
!
!  spline along X-axis
!
  do j=1,ny
    if(imi(j).gt.0) then
      nsi=ima(j)-imi(j)+1
      do l=1,6
        do i=imi(j),ima(j)
          ai(i-imi(j)+1)=spl(1,l,ipoint(i,j))
        enddo   
        call spl_five_reg(nsi,hx,ai,bi,ci,di,ei,fi)
        do i=imi(j),ima(j)
          ii=i-imi(j)+1
          spl(2,l,ipoint(i,j))=bi(ii)
          spl(3,l,ipoint(i,j))=ci(ii)
          spl(4,l,ipoint(i,j))=di(ii)
          spl(5,l,ipoint(i,j))=ei(ii)
          spl(6,l,ipoint(i,j))=fi(ii)
        enddo   
      enddo   
    endif
  enddo   
!
  deallocate( ai,bi,ci,di,ei,fi )
!
  return
  end
!
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine spline(nx,ny,x,y,hx,hy,icount,spl,ipoint,xb,yb,u,ux,uy, &
    uxx,uxy,uyy,ierr)
!
! Evaluates interpolated value u(x,y) and its derivatives ux,uy,uxx,uyy,uxy 
! using arrays calculated by  s2dcut
! see also comments in subroutine s2dcut
!  ierr = 1 - point out of the domain

  implicit double precision (a-h,o-z)
!
  dimension spl(6,6,icount),x(nx),y(ny),ipoint(nx,ny)
  dimension a(6),ax(6),axx(6)
!
  xk=(xb-x(1))/hx
  kx=int(xk)+1
  kx = min(nx,max(1,kx))
  yk=(yb-y(1))/hy
  ky=int(yk)+1
  ky = min(ny,max(1,ky))
!  if( kx.lt.1 .or. kx.gt.nx .or. ky.lt.1 .or. ky.gt.ny &
!     .or. ipoint(kx,ky).lt.0 ) then
!   print *,'spline: out of range', xb,yb
!    ierr=1
!    return
!  endif
  ierr=0
  dx=xb-x(kx)
  dy=yb-y(ky)
  do l=1,6
    a(l) =     spl(1,l,ipoint(kx,ky)) &
         + dx*(spl(2,l,ipoint(kx,ky)) &
         + dx*(spl(3,l,ipoint(kx,ky)) &
         + dx*(spl(4,l,ipoint(kx,ky)) &
         + dx*(spl(5,l,ipoint(kx,ky)) &
         + dx* spl(6,l,ipoint(kx,ky))))))
    ax(l) =          spl(2,l,ipoint(kx,ky)) &
          + dx*(2.d0*spl(3,l,ipoint(kx,ky)) &
          + dx*(3.d0*spl(4,l,ipoint(kx,ky)) &
          + dx*(4.d0*spl(5,l,ipoint(kx,ky)) &
          + dx* 5.d0*spl(6,l,ipoint(kx,ky)))))
    axx(l) = 2.d0*spl(3,l,ipoint(kx,ky)) &
           + dx*(6.d0*spl(4,l,ipoint(kx,ky)) &
           + dx*(12.d0*spl(5,l,ipoint(kx,ky)) &
           + dx*(20.d0*spl(6,l,ipoint(kx,ky)))))
  enddo

  u = a(1) + dy*(a(2) + dy*(a(3) + dy*(a(4) + dy*(a(5) + dy*a(6)))))
  ux = ax(1) + dy*(ax(2) + dy*(ax(3) + dy*(ax(4) + dy*(ax(5) + dy*ax(6)))))
  uy = a(2) + dy*(2.d0*a(3) + dy*(3.d0*a(4) + dy*(4.d0*a(5) + dy*5.d0*a(6))))
  uxx =  axx(1) + dy*(axx(2) + dy*(axx(3) + dy*(axx(4) + dy*(axx(5) + dy*axx(6)))))
  uxy= ax(2) + dy*(2.d0*ax(3) + dy*(3.d0*ax(4) + dy*(4.d0*ax(5) + dy*5.d0*ax(6))))
  uyy = 2.d0*a(3) + dy*(6.d0*a(4) + dy*(12.d0*a(5) + dy*20.d0*a(6)))

  return
  end
