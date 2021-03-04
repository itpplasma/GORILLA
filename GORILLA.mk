#FC       = f95f
#FC       = f95-lah
#FC       = ifort
FC       = gfortran

#OPTS= -M OBJS --chk a,e,s,u,x --trace --trap -g
#OPTS= -M OBJS -O 
#OPTS= -module OBJS
#OPTS= -J OBJS -O -fopenmp 
#OPTS= -J OBJS -O
#vOPTS= -J OBJS -O0 -g -fbounds-check
OPTS= -J OBJS  -g -fbacktrace -ffpe-trap=zero,overflow,invalid  -fbounds-check -fopenmp
#OPTS= -J OBJS  -g -fbacktrace -ffpe-trap=zero,overflow,invalid  -fbounds-check
#OPTS = -J OBJS -Wall -pedantic
#OPTS = -J OBJS -Wuninitialized
#OPTS = -J OBJS -g
# The option below shows in which line the error occurs
#OPTS= -J OBJS  -ffpe-trap=invalid,zero,overflow -g -fopenmp

NCINC ?= -I/usr/include
#NCINC = -I/proj/plasma/Libs/gcc-9/NetCDF/include

NCLIB ?= -L/usr/lib -lnetcdff -lnetcdf -llapack
#NCLIB = -L/proj/plasma/Libs/gcc-9/NetCDF/lib -lnetcdff -lnetcdf

OBJS = OBJS/SetWorkingPrecision.o \
	OBJS/Polynomial234RootSolvers.o \
	OBJS/constants_mod.o \
	OBJS/tetra_grid_settings_mod.o \
	OBJS/gorilla_settings_mod.o \
	OBJS/various_functions_mod.o \
	OBJS/gorilla_diag_mod.o \
	OBJS/canonical_coordinates_mod.o \
	OBJS/nctools_module.o \
	OBJS/runge_kutta_mod.o \
	OBJS/magfie.o \
	OBJS/chamb_m.o \
	OBJS/vmecinm_m.o \
	OBJS/spl_three_to_five_mod.o \
	OBJS/spline_vmec_data.o \
	OBJS/new_vmec_allocation_stuff.o \
	OBJS/binsrc.o \
	OBJS/field_divB0.o \
	OBJS/scaling_r_theta.o\
	OBJS/field_line_integration_for_SYNCH.o \
	OBJS/preload_for_SYNCH.o \
	OBJS/plag_coeff.o \
	OBJS/magdata_in_symfluxcoord.o \
	OBJS/points_2d.o\
	OBJS/circular_mesh.o\
	OBJS/tetra_grid_mod.o \
	OBJS/make_grid_rect.o \
	OBJS/bdivfree.o \
	OBJS/tetra_physics_mod.o \
	OBJS/tetra_physics_poly_precomp_mod.o \
	OBJS/differentiate.o \
	OBJS/spline5_RZ.o \
	OBJS/supporting_functions_mod.o \
	OBJS/pusher_tetra_func_mod.o \
	OBJS/pusher_tetra_poly.o \
	OBJS/pusher_tetra_rk.o \
	OBJS/get_canonical_coordinates.o \
	OBJS/orbit_timestep_gorilla.o \
	OBJS/gorilla_plot_mod.o \
	OBJS/test_gorilla_main.o


test_gorilla_main.x: $(OBJS) Gorilla.mk
	$(FC) $(OPTS) -o test_gorilla_main.x $(OBJS) $(NCLIB)
OBJS/SetWorkingPrecision.o: SRC/SetWorkingPrecision.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/SetWorkingPrecision.f90
	mv SetWorkingPrecision.o OBJS/
OBJS/Polynomial234RootSolvers.o: SRC/Polynomial234RootSolvers.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/Polynomial234RootSolvers.f90
	mv Polynomial234RootSolvers.o OBJS/
OBJS/constants_mod.o: SRC/constants_mod.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/constants_mod.f90
	mv constants_mod.o OBJS/
OBJS/tetra_grid_settings_mod.o: SRC/tetra_grid_settings_mod.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/tetra_grid_settings_mod.f90
	mv tetra_grid_settings_mod.o OBJS/
OBJS/gorilla_settings_mod.o: SRC/gorilla_settings_mod.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/gorilla_settings_mod.f90
	mv gorilla_settings_mod.o OBJS/
OBJS/various_functions_mod.o: SRC/various_functions_mod.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/various_functions_mod.f90
	mv various_functions_mod.o OBJS/
OBJS/gorilla_diag_mod.o: SRC/gorilla_diag_mod.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/gorilla_diag_mod.f90
	mv gorilla_diag_mod.o OBJS/
OBJS/canonical_coordinates_mod.o: SRC/canonical_coordinates_mod.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/canonical_coordinates_mod.f90
	mv canonical_coordinates_mod.o OBJS/
OBJS/nctools_module.o: SRC/nctools_module.f90 Gorilla.mk SRC/canonical_coordinates_mod.f90
	$(FC) $(OPTS) $(NCINC) -c SRC/nctools_module.f90 
	mv nctools_module.o OBJS/	
OBJS/runge_kutta_mod.o: SRC/runge_kutta_mod.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/runge_kutta_mod.f90
	mv runge_kutta_mod.o OBJS/	
OBJS/vmecinm_m.o: SRC/vmecinm_m.f90 Gorilla.mk SRC/canonical_coordinates_mod.f90
	$(FC) $(OPTS) -c SRC/vmecinm_m.f90
	mv vmecinm_m.o OBJS/
OBJS/magfie.o: SRC/magfie.f90 Gorilla.mk SRC/magfie.f90
	$(FC) $(OPTS) -c SRC/magfie.f90
	mv magfie.o OBJS/
OBJS/chamb_m.o: SRC/chamb_m.f90 Gorilla.mk SRC/chamb_m.f90
	$(FC) $(OPTS) -c SRC/chamb_m.f90
	mv chamb_m.o OBJS/
OBJS/spline_vmec_data.o: SRC/spline_vmec_data.f90 Gorilla.mk SRC/canonical_coordinates_mod.f90
	$(FC) $(OPTS) -c SRC/spline_vmec_data.f90
	mv spline_vmec_data.o OBJS/	
OBJS/new_vmec_allocation_stuff.o: SRC/new_vmec_allocation_stuff.f90 Gorilla.mk SRC/canonical_coordinates_mod.f90
	$(FC) $(OPTS) -c SRC/new_vmec_allocation_stuff.f90
	mv new_vmec_allocation_stuff.o OBJS/	
OBJS/make_grid_rect.o: SRC/make_grid_rect.f90 Gorilla.mk SRC/various_functions_mod.f90
	$(FC) $(OPTS) -c SRC/make_grid_rect.f90
	mv make_grid_rect.o OBJS/
OBJS/tetra_grid_mod.o: SRC/tetra_grid_mod.f90 Gorilla.mk SRC/make_grid_rect.f90
	$(FC) $(OPTS) -c SRC/tetra_grid_mod.f90
	mv tetra_grid_mod.o OBJS/
OBJS/tetra_physics_mod.o: SRC/tetra_physics_mod.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/tetra_physics_mod.f90
	mv tetra_physics_mod.o OBJS/
OBJS/tetra_physics_poly_precomp_mod.o: SRC/tetra_physics_poly_precomp_mod.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/tetra_physics_poly_precomp_mod.f90
	mv tetra_physics_poly_precomp_mod.o OBJS/
OBJS/pusher_tetra_func_mod.o: SRC/pusher_tetra_func_mod.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/pusher_tetra_func_mod.f90
	mv pusher_tetra_func_mod.o OBJS/
OBJS/pusher_tetra_poly.o: SRC/pusher_tetra_poly.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/pusher_tetra_poly.f90
	mv pusher_tetra_poly.o OBJS/
OBJS/orbit_timestep_gorilla.o: SRC/orbit_timestep_gorilla.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/orbit_timestep_gorilla.f90
	mv orbit_timestep_gorilla.o OBJS/
OBJS/gorilla_plot_mod.o: SRC/gorilla_plot_mod.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/gorilla_plot_mod.f90
	mv gorilla_plot_mod.o OBJS/
OBJS/pusher_tetra_rk.o: SRC/pusher_tetra_rk.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/pusher_tetra_rk.f90
	mv pusher_tetra_rk.o OBJS/
OBJS/binsrc.o: SRC/binsrc.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/binsrc.f90
	mv binsrc.o OBJS/
OBJS/field_divB0.o: SRC/field_divB0.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/field_divB0.f90
	mv field_divB0.o OBJS/
OBJS/bdivfree.o: SRC/bdivfree.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/bdivfree.f90
	mv bdivfree.o OBJS/
OBJS/spline5_RZ.o: SRC/spline5_RZ.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/spline5_RZ.f90
	mv spline5_RZ.o OBJS/
OBJS/supporting_functions_mod.o: SRC/supporting_functions_mod.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/supporting_functions_mod.f90
	mv supporting_functions_mod.o OBJS/
OBJS/differentiate.o: SRC/differentiate.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/differentiate.f90
	mv differentiate.o OBJS/
OBJS/test_gorilla_main.o: SRC/test_gorilla_main.f90 Gorilla.mk SRC/various_functions_mod.f90
	$(FC) $(OPTS) -c SRC/test_gorilla_main.f90
	mv test_gorilla_main.o OBJS/
OBJS/points_2d.o: SRC/points_2d.f90 Gorilla.mk SRC/field_divB0.f90
	$(FC) $(OPTS) -c SRC/points_2d.f90
	mv points_2d.o OBJS/
OBJS/circular_mesh.o: SRC/circular_mesh.f90 Gorilla.mk SRC/points_2d.f90
	$(FC) $(OPTS) -c SRC/circular_mesh.f90
	mv circular_mesh.o OBJS/
OBJS/scaling_r_theta.o: SRC/scaling_r_theta.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/scaling_r_theta.f90
	mv scaling_r_theta.o OBJS/
OBJS/spl_three_to_five_mod.o: SRC/spl_three_to_five_mod.f90 Gorilla.mk SRC/canonical_coordinates_mod.f90
	$(FC) $(OPTS) -c SRC/spl_three_to_five_mod.f90
	mv spl_three_to_five_mod.o OBJS/
OBJS/plag_coeff.o: SRC/plag_coeff.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/plag_coeff.f90
	mv plag_coeff.o OBJS/
OBJS/field_line_integration_for_SYNCH.o: SRC/field_line_integration_for_SYNCH.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/field_line_integration_for_SYNCH.f90
	mv field_line_integration_for_SYNCH.o OBJS/
OBJS/preload_for_SYNCH.o: SRC/preload_for_SYNCH.f90 SRC/field_line_integration_for_SYNCH.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/preload_for_SYNCH.f90
	mv preload_for_SYNCH.o OBJS/
OBJS/magdata_in_symfluxcoord.o: SRC/magdata_in_symfluxcoord.f90 SRC/spl_three_to_five_mod.f90 SRC/plag_coeff.f90 Gorilla.mk
	$(FC) $(OPTS) -c SRC/magdata_in_symfluxcoord.f90
	mv magdata_in_symfluxcoord.o OBJS/
OBJS/get_canonical_coordinates.o: SRC/get_canonical_coordinates.f90 Gorilla.mk SRC/canonical_coordinates_mod.f90
	$(FC) $(OPTS) -c SRC/get_canonical_coordinates.f90
	mv get_canonical_coordinates.o OBJS/
OBJS/.o: .f90 Gorilla.mk
	$(FC) $(OPTS) -c .f90
	mv .o OBJS/


.PHONY: clean
clean:
	rm -f OBJS/*
	rm -f SRC/*.mod
	rm -f test_gorilla_main.x
