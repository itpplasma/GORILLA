!
  subroutine make_grid_rect(tetra,verts_rphiz,grid_size,Rmin,Rmax,Zmin,Zmax)
!
  use tetra_grid_mod, only: ntetr,nvert,tetrahedron_grid
  use constants, only: pi
!   
  implicit none
!
  type(tetrahedron_grid), dimension(ntetr), intent(out) :: tetra
  double precision,  dimension(3,nvert), intent(out) :: verts_rphiz
  integer, dimension(3), intent(in) :: grid_size
  double precision, intent(in) :: Rmin,Rmax,Zmin,Zmax
  integer :: ir,iphi,iz,iv,itetr,i,j,k,l
  integer :: ind_tetr,ind_tetr1,ind_tetr2,iface1,iface2
  integer :: nvertinner,nper
  integer, dimension(4) :: knots1,knots2  
  integer, dimension(3,4,3) :: ivert_prism
  double precision :: hr,hphi,hz,r,phi,z,x,y,s
  double precision :: R_c,Z_c
  double precision, dimension(4)   :: sp,pp,zp
  integer,          dimension(:,:,:),       allocatable :: inodes,itetrbeg
  integer :: nr,nphi,nz
!
  nr = grid_size(1)
  nphi = grid_size(2)
  nz = grid_size(3)
!
! vertices of the first tetrahedron:
! point 1:
  ivert_prism(1,1,1)=0
  ivert_prism(2,1,1)=0
  ivert_prism(3,1,1)=0
! point 2:
  ivert_prism(1,2,1)=1
  ivert_prism(2,2,1)=0
  ivert_prism(3,2,1)=0
! point 3:
  ivert_prism(1,3,1)=0
  ivert_prism(2,3,1)=1
  ivert_prism(3,3,1)=0
! point 4:
  ivert_prism(1,4,1)=0
  ivert_prism(2,4,1)=0
  ivert_prism(3,4,1)=1
!
! vertices of the second tetrahedron:
! point 1:
  ivert_prism(:,1,2)=ivert_prism(:,2,1)
! point 2:
  ivert_prism(:,2,2)=ivert_prism(:,3,1)
! point 3:
  ivert_prism(:,3,2)=ivert_prism(:,4,1)
! point 4:
  ivert_prism(1,4,2)=1
  ivert_prism(2,4,2)=0
  ivert_prism(3,4,2)=1
!
! vertices of the third tetrahedron:
! point 1:
  ivert_prism(:,1,3)=ivert_prism(:,3,1)
! point 2:
  ivert_prism(:,2,3)=ivert_prism(:,4,1)
! point 3:
  ivert_prism(:,3,3)=ivert_prism(:,4,2)
! point 4:
  ivert_prism(1,4,3)=0
  ivert_prism(2,4,3)=1
  ivert_prism(3,4,3)=1
!
  hr=(Rmax-Rmin)/nr
  hphi=2.d0*pi/nphi
  hz=(Zmax-Zmin)/nz
!
  R_c=0.5d0*(Rmax+Rmin)
  Z_c=0.5d0*(Zmax+Zmin)
  nper=0
!
!
  allocate(inodes(0:nr,0:nz,0:nphi))
!
  iv=0
!
  do iphi=0,nphi
    phi=hphi*iphi
    do ir=0,nr
      do iz=0,nz
        r=rmin+hr*ir
        z=zmin+hz*iz
!
        x=(r-R_c)
        y=(z-Z_c)
        r=R_c+x*cos(nper*phi)+y*sin(nper*phi)
        z=Z_c-x*sin(nper*phi)+y*cos(nper*phi)! +hz/3000
        s=r
!
        iv=iv+1
        inodes(ir,iz,iphi)=iv
!
        verts_rphiz(:,iv) = [s,phi,z]
!
      enddo
    enddo
  enddo
!
  allocate(itetrbeg(nr,nz,nphi))
!
  ind_tetr=0
!
  do iphi=1,nphi
    do ir=1,nr
      do iz=1,nz
!
        itetrbeg(ir,iz,iphi)=ind_tetr
!
        do itetr=1,3
          ind_tetr=ind_tetr+1
          do i=1,4
            iv=inodes(ir-1+ivert_prism(1,i,itetr),iz-1+ivert_prism(2,i,itetr),iphi-1+ivert_prism(3,i,itetr))
            tetra(ind_tetr)%ind_knot(i)=iv
            tetra(ind_tetr)%neighbour_tetr(i)=-1
            tetra(ind_tetr)%neighbour_face(i)=-1
            tetra(ind_tetr)%neighbour_perbou_phi(i)=0
          enddo
        enddo
!
        do itetr=1,3
          ind_tetr=ind_tetr+1
          do i=1,4
            iv=inodes(ir-ivert_prism(1,i,itetr),iz-ivert_prism(2,i,itetr),iphi-ivert_prism(3,i,itetr))
            tetra(ind_tetr)%ind_knot(i)=iv
            tetra(ind_tetr)%neighbour_tetr(i)=-1
            tetra(ind_tetr)%neighbour_face(i)=-1
            tetra(ind_tetr)%neighbour_perbou_phi(i)=0
          enddo
        enddo
!
        do i=1,6
          ind_tetr1=itetrbeg(ir,iz,iphi)+i
          do j=1,6
            if(j.eq.i) cycle
            ind_tetr2=itetrbeg(ir,iz,iphi)+j
            call check_neighbour(tetra(ind_tetr1)%ind_knot(:),              &
                                 tetra(ind_tetr2)%ind_knot(:),iface1,iface2)
            if(iface1.ne.-1) then
              tetra(ind_tetr1)%neighbour_tetr(iface1)=ind_tetr2
              tetra(ind_tetr1)%neighbour_face(iface1)=iface2
            endif
          enddo
        enddo
      enddo
    enddo
  enddo
!
!
  do iphi=1,nphi
    do ir=1,nr
      do iz=1,nz
!
        if(ir.gt.1) then
          do i=1,3
            ind_tetr1=itetrbeg(ir,iz,iphi)+i
            do j=4,6
              ind_tetr2=itetrbeg(ir-1,iz,iphi)+j
              call check_neighbour(tetra(ind_tetr1)%ind_knot(:),              &
                                   tetra(ind_tetr2)%ind_knot(:),iface1,iface2)
              if(iface1.ne.-1) then
                tetra(ind_tetr1)%neighbour_tetr(iface1)=ind_tetr2
                tetra(ind_tetr1)%neighbour_face(iface1)=iface2
              endif
            enddo
          enddo
        endif
!
        if(ir.lt.nr) then
          do i=4,6
            ind_tetr1=itetrbeg(ir,iz,iphi)+i
            do j=1,3
              ind_tetr2=itetrbeg(ir+1,iz,iphi)+j
              call check_neighbour(tetra(ind_tetr1)%ind_knot(:),              &
                                   tetra(ind_tetr2)%ind_knot(:),iface1,iface2)
              if(iface1.ne.-1) then
                tetra(ind_tetr1)%neighbour_tetr(iface1)=ind_tetr2
                tetra(ind_tetr1)%neighbour_face(iface1)=iface2
              endif
            enddo
          enddo
        endif
!
        if(iz.gt.1) then
          do i=1,3
            ind_tetr1=itetrbeg(ir,iz,iphi)+i
            do j=4,6
              ind_tetr2=itetrbeg(ir,iz-1,iphi)+j
              call check_neighbour(tetra(ind_tetr1)%ind_knot(:),              &
                                   tetra(ind_tetr2)%ind_knot(:),iface1,iface2)
              if(iface1.ne.-1) then
                tetra(ind_tetr1)%neighbour_tetr(iface1)=ind_tetr2
                tetra(ind_tetr1)%neighbour_face(iface1)=iface2
              endif
            enddo
          enddo
        endif
!
        if(iz.lt.nz) then
          do i=4,6
            ind_tetr1=itetrbeg(ir,iz,iphi)+i
            do j=1,3
              ind_tetr2=itetrbeg(ir,iz+1,iphi)+j
              call check_neighbour(tetra(ind_tetr1)%ind_knot(:),              &
                                   tetra(ind_tetr2)%ind_knot(:),iface1,iface2)
              if(iface1.ne.-1) then
                tetra(ind_tetr1)%neighbour_tetr(iface1)=ind_tetr2
                tetra(ind_tetr1)%neighbour_face(iface1)=iface2
              endif
            enddo
          enddo
        endif
!
        if(iphi.gt.1) then
          do i=1,3
            ind_tetr1=itetrbeg(ir,iz,iphi)+i
            do j=1,3
              ind_tetr2=itetrbeg(ir,iz,iphi-1)+j
              call check_neighbour(tetra(ind_tetr1)%ind_knot(:),              &
                                   tetra(ind_tetr2)%ind_knot(:),iface1,iface2)
              if(iface1.ne.-1) then
                tetra(ind_tetr1)%neighbour_tetr(iface1)=ind_tetr2
                tetra(ind_tetr1)%neighbour_face(iface1)=iface2
              endif
            enddo
          enddo
          do i=4,6
            ind_tetr1=itetrbeg(ir,iz,iphi)+i
            do j=4,6
              ind_tetr2=itetrbeg(ir,iz,iphi-1)+j
              call check_neighbour(tetra(ind_tetr1)%ind_knot(:),              &
                                   tetra(ind_tetr2)%ind_knot(:),iface1,iface2)
              if(iface1.ne.-1) then
                tetra(ind_tetr1)%neighbour_tetr(iface1)=ind_tetr2
                tetra(ind_tetr1)%neighbour_face(iface1)=iface2
              endif
            enddo
          enddo
        endif
!
        if(iphi.lt.nphi) then
          do i=1,3
            ind_tetr1=itetrbeg(ir,iz,iphi)+i
            do j=1,3
              ind_tetr2=itetrbeg(ir,iz,iphi+1)+j
              call check_neighbour(tetra(ind_tetr1)%ind_knot(:),              &
                                   tetra(ind_tetr2)%ind_knot(:),iface1,iface2)
              if(iface1.ne.-1) then
                tetra(ind_tetr1)%neighbour_tetr(iface1)=ind_tetr2
                tetra(ind_tetr1)%neighbour_face(iface1)=iface2
              endif
            enddo
          enddo
          do i=4,6
            ind_tetr1=itetrbeg(ir,iz,iphi)+i
            do j=4,6
              ind_tetr2=itetrbeg(ir,iz,iphi+1)+j
              call check_neighbour(tetra(ind_tetr1)%ind_knot(:),              &
                                   tetra(ind_tetr2)%ind_knot(:),iface1,iface2)
              if(iface1.ne.-1) then
                tetra(ind_tetr1)%neighbour_tetr(iface1)=ind_tetr2
                tetra(ind_tetr1)%neighbour_face(iface1)=iface2
              endif
            enddo
          enddo
        endif
!
      enddo
    enddo
  enddo
!
! Neighbours through the periodic boundary
!
  nvertinner=(nr+1)*(nz+1)*nphi
!
  do ir=1,nr
    do iz=1,nz
!
      do i=1,3
        ind_tetr1=itetrbeg(ir,iz,1)+i
        knots1(:)=tetra(ind_tetr1)%ind_knot(:)
        do j=1,3
          ind_tetr2=itetrbeg(ir,iz,nphi)+j
          knots2(:)=tetra(ind_tetr2)%ind_knot(:)
          knots2(:)=modulo(knots2(:),nvertinner)
          call check_neighbour(knots1,knots2,iface1,iface2)
          if(iface1.ne.-1) then
            tetra(ind_tetr1)%neighbour_tetr(iface1)=ind_tetr2
            tetra(ind_tetr1)%neighbour_face(iface1)=iface2
            tetra(ind_tetr1)%neighbour_perbou_phi(iface1)=-1
            tetra(ind_tetr2)%neighbour_tetr(iface2)=ind_tetr1
            tetra(ind_tetr2)%neighbour_face(iface2)=iface1
            tetra(ind_tetr2)%neighbour_perbou_phi(iface2)=1
          endif
        enddo
      enddo
!
      do i=4,6
        ind_tetr1=itetrbeg(ir,iz,1)+i
        knots1(:)=tetra(ind_tetr1)%ind_knot(:)
        do j=4,6
          ind_tetr2=itetrbeg(ir,iz,nphi)+j
          knots2(:)=tetra(ind_tetr2)%ind_knot(:)
          knots2(:)=modulo(knots2(:),nvertinner)
          call check_neighbour(knots1,knots2,iface1,iface2)
          if(iface1.ne.-1) then
            tetra(ind_tetr1)%neighbour_tetr(iface1)=ind_tetr2
            tetra(ind_tetr1)%neighbour_face(iface1)=iface2
            tetra(ind_tetr1)%neighbour_perbou_phi(iface1)=-1
            tetra(ind_tetr2)%neighbour_tetr(iface2)=ind_tetr1
            tetra(ind_tetr2)%neighbour_face(iface2)=iface1
            tetra(ind_tetr2)%neighbour_perbou_phi(iface2)=1
          endif
        enddo
      enddo
!
    enddo
  enddo
!
! Check the number of tetrahedrons at the XY boundary
  j=0
!
b: do ind_tetr=1,ntetr
  do i=1,4
     if(tetra(ind_tetr)%neighbour_tetr(i).lt.1) then
       j=j+1
       cycle b
     endif
   enddo
  enddo b
!
  deallocate(inodes,itetrbeg)
!
  end subroutine make_grid_rect
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine check_neighbour(knots1,knots2,iface1,iface2)
!
  implicit none
!
  logical :: result
  integer :: iface1,iface2,i,j,matches
  integer, dimension(4) :: knots1,knots2
  integer, dimension(3) :: match1,match2
  logical, dimension(4) :: unused
!
  iface1=-1
  iface2=-1
  matches=0
!
  a: do i=1,4
    do j=1,4
      if(knots1(i).eq.knots2(j)) then
        matches=matches+1
        match1(matches)=i
        match2(matches)=j
        if(matches.eq.3) exit a
      endif
    enddo
    if(i-matches.gt.1) return
  enddo a
!
  unused= .true.
  do i=1,3
    unused(match1(i))= .false.
  enddo
  do i=1,4
    if(unused(i)) then
      iface1=i
      exit
    endif
  enddo
!
  unused= .true.
  do i=1,3
    unused(match2(i))= .false.
  enddo
  do i=1,4
    if(unused(i)) then
      iface2=i
      exit
    endif
  enddo
!
  end subroutine check_neighbour
!
!ccccccccccccccccccccccccccccccccccccccccccc
!
