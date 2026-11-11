! Build the GORILLA grid from tetra_grid.inp and dump the first toroidal slice:
!   verts_2d.dat : i, R, Z, (s, theta if available)
!   tets_slice1.dat : 4 knot indices of every tetrahedron in the first slice
!   core_settings.nml : grid_kind 2 settings and g-file signature (edge mesh settings record)
program dump_mesh_2d
  use tetra_grid_settings_mod, only: load_tetra_grid_inp, n2, grid_kind, write_edge_mesh_core_record
  use tetra_grid_mod, only: make_tetra_grid, verts_rphiz, verts_sthetaphi, tetra_grid, ntetr, nvert
  implicit none
  integer :: i, iu, nvs, nts
  call load_tetra_grid_inp()
  call make_tetra_grid()
  nvs = nvert / n2
  nts = ntetr / n2
  open(newunit=iu, file='verts_2d.dat')
  write(iu,*) nvs, nvert, ntetr, n2
  do i = 1, nvs
    if (allocated(verts_sthetaphi)) then
      write(iu,'(i8,4es24.16)') i, verts_rphiz(1,i), verts_rphiz(3,i), verts_sthetaphi(1,i), verts_sthetaphi(2,i)
    else
      write(iu,'(i8,4es24.16)') i, verts_rphiz(1,i), verts_rphiz(3,i), 0d0, 0d0
    end if
  end do
  close(iu)
  open(newunit=iu, file='tets_slice1.dat')
  do i = 1, nts
    write(iu,'(4i10)') tetra_grid(i)%ind_knot
  end do
  close(iu)
  call write_edge_mesh_core_record('core_settings.nml')
  print *, 'dumped', nvs, 'verts/slice,', nts, 'tets/slice, grid_kind', grid_kind
end program
