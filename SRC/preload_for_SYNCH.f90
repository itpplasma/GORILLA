module preload_for_SYNCH_mod

  implicit none

  contains

  subroutine preload_for_SYNCH

    use field_line_integration_for_SYNCH_mod, only: field_line_integration_for_SYNCH
!
  implicit none
!
  integer :: nstep,nsurfmax,nlabel,ntheta,i,iunit
!
  double precision :: rmn,rmx,zmn,zmx,raxis,zaxis
  double precision, dimension(:),   allocatable :: rbeg,rsmall,qsaf,psisurf,phitor
  double precision, dimension(:,:), allocatable :: R_st,Z_st,bmod_st,sqgnorm_st
!
  open(newunit=iunit,file='preload_for_SYNCH.inp')
  read (iunit,*) nstep    !number of integration steps
  read (iunit,*) nlabel   !grid size over radial variabl
  read (iunit,*) ntheta   !grid size over poloidal angle
  read (iunit,*) nsurfmax !number of starting points between the
                          !magnetic axis and right box boundary
                          !when searching for the separatrix
  close(iunit)
!
  allocate(rbeg(nlabel),rsmall(nlabel),qsaf(nlabel),psisurf(nlabel),phitor(nlabel))
  allocate(R_st(nlabel,ntheta),Z_st(nlabel,ntheta),bmod_st(nlabel,ntheta),sqgnorm_st(nlabel,ntheta))
!
  call field_line_integration_for_SYNCH(nstep,nsurfmax,nlabel,ntheta,    &
                                        rmn,rmx,zmn,zmx,raxis,zaxis,     &
                                        rbeg,rsmall,qsaf,psisurf,phitor, &
                                        R_st,Z_st,bmod_st,sqgnorm_st)
!
  open(newunit=iunit,form='formatted',file='box_size_axis.dat')
  write (iunit,*) rmn,rmx, '<= rmn, rmx (cm)'
  write (iunit,*) zmn,zmx, '<= zmn, zmx (cm)'
  write (iunit,*) raxis,zaxis, '<= raxis, zaxis (cm)'
  close(iunit)
!
  open(newunit=iunit,form='formatted',file='flux_functions.dat')
  write (iunit,*) '# R_beg, r,  q, psi_pol, psi_tor'
! 2020-02-20: Lukas Bauer: rbeg is uninitialized and never used for computation. Still 'random' values are written in file.
  do i=1,nlabel
    write (iunit,*) rbeg(i),rsmall(i),qsaf(i),psisurf(i),phitor(i)
  enddo
  close(iunit)
!
  open(newunit=iunit,form='formatted',file='twodim_functions.dat')
  write (iunit,*) nlabel, ntheta, '<= nlabel, ntheta'
  write (iunit,*) 'R(label,theta)'
  do i=1,nlabel
    write (iunit,*) R_st(i,:)
  enddo
  write (iunit,*) 'Z(label,theta)'
  do i=1,nlabel
    write (iunit,*) Z_st(i,:)
  enddo
  write (iunit,*) 'B(label,theta)'
  do i=1,nlabel
    write (iunit,*) bmod_st(i,:)
  enddo
  write (iunit,*) 'sqrtg_norm(label,theta)'
  do i=1,nlabel
    write (iunit,*) sqgnorm_st(i,:)
  enddo
  close(iunit)
!
  deallocate(rbeg,rsmall,qsaf,psisurf,phitor)
  deallocate(R_st,Z_st,bmod_st,sqgnorm_st)
!
  end subroutine preload_for_SYNCH

end module preload_for_SYNCH_mod
