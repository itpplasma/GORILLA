  subroutine preload_for_SYNCH
!
  implicit none
!
  integer :: nstep,nsurfmax,nlabel,ntheta,i
!
  double precision :: rmn,rmx,zmn,zmx,raxis,zaxis
  double precision, dimension(:),   allocatable :: rbeg,rsmall,qsaf,psisurf,phitor
  double precision, dimension(:,:), allocatable :: R_st,Z_st,bmod_st,sqgnorm_st
!
  open(1,file='preload_for_SYNCH.inp')
  read (1,*) nstep    !number of integration steps
  read (1,*) nlabel   !grid size over radial variabl
  read (1,*) ntheta   !grid size over poloidal angle
  read (1,*) nsurfmax !number of starting points between the
                      !magnetic axis and right box boundary
                      !when searching for the separatrix
  close(1)
!
  allocate(rbeg(nlabel),rsmall(nlabel),qsaf(nlabel),psisurf(nlabel),phitor(nlabel))
  allocate(R_st(nlabel,ntheta),Z_st(nlabel,ntheta),bmod_st(nlabel,ntheta),sqgnorm_st(nlabel,ntheta))
!
  call field_line_integration_for_SYNCH(nstep,nsurfmax,nlabel,ntheta,    &
                                        rmn,rmx,zmn,zmx,raxis,zaxis,     &
                                        rbeg,rsmall,qsaf,psisurf,phitor, &
                                        R_st,Z_st,bmod_st,sqgnorm_st)
!
  open(1,form='formatted',file='box_size_axis.dat')
  write (1,*) rmn,rmx, '<= rmn, rmx (cm)'
  write (1,*) zmn,zmx, '<= zmn, zmx (cm)'
  write (1,*) raxis,zaxis, '<= raxis, zaxis (cm)'
  close(1)
!
  open(1,form='formatted',file='flux_functions.dat')
  write (1,*) '# R_beg, r,  q, psi_pol, psi_tor'
! 2020-02-20: Lukas Bauer: rbeg is uninitialized and never used for computation. Still 'random' values are written in file.  
  do i=1,nlabel
    write (1,*) rbeg(i),rsmall(i),qsaf(i),psisurf(i),phitor(i)
  enddo
  close(1)
!
  open(1,form='formatted',file='twodim_functions.dat')
  write (1,*) nlabel, ntheta, '<= nlabel, ntheta'
  write (1,*) 'R(label,theta)'
  do i=1,nlabel
    write (1,*) R_st(i,:)
  enddo
  write (1,*) 'Z(label,theta)'
  do i=1,nlabel
    write (1,*) Z_st(i,:)
  enddo
  write (1,*) 'B(label,theta)'
  do i=1,nlabel
    write (1,*) bmod_st(i,:)
  enddo
  write (1,*) 'sqrtg_norm(label,theta)'
  do i=1,nlabel
    write (1,*) sqgnorm_st(i,:)
  enddo
  close(1)
!
  deallocate(rbeg,rsmall,qsaf,psisurf,phitor)
  deallocate(R_st,Z_st,bmod_st,sqgnorm_st)
!
  end subroutine preload_for_SYNCH
