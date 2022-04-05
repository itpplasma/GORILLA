%#######################################################################################################################
% description:
%-----------------------------------------------------------------------------------------------------------------------
% * Compute a collisionless guiding-center orbit with GORILLA for a trapped Deuterium particle.
% * Use a field-aligned grid for a non-axisymmetric VMEC MHD equilibrium.
% * Compare the results of GORILLA with different polynominal orders and Runge-Kutta 4.
% * Create a figure with the Poincaré sections (v_\parallel = 0) in cylindrical and symmetry flux coordinates.
% * Compute the normalized parallel adiabatic invariant as a function of banana bounces.
%
%#######################################################################################################################
% used functions:
%-----------------------------------------------------------------------------------------------------------------------
% InputFile
% NameList
% read_in
% load_copy
%#######################################################################################################################
% authors: Michael Scheidt, Michael Eder
% created: 19.02.2021


%Initialize used paths and files
%-----------------------------------------------------------------------------------------------------------------------

%Name of the current calculation to create folder
name_test_case='example_1';

%path of Matlab script
path_script=pwd;

%main path of GORILLA
c=strsplit(path_script,'/');
path_main=strjoin(c(1:end-1),'/');

%path to run GORILLA code
mkdir(path_main,['EXAMPLES/MATLAB_RUN/',name_test_case]);
path_RUN=[path_main,'/EXAMPLES/MATLAB_RUN/',name_test_case];

%path of input files (blueprints)
path_inp_files=[path_main,'/INPUT'];

%path of the used functions
path_functions=[path_main,'/MATLAB/functions'];

%define path for data and plots and create a new folder
mkdir([path_script,'/data_plots'],name_test_case);
path_data_plots=[path_script,'/data_plots/',name_test_case];

%Change to RUN folder and clean it
cd(path_RUN);

%add path
addpath(path_inp_files);
addpath(path_functions);

%initialize class for inputfiles
gorilla = InputFile([path_inp_files,'/gorilla.inp']);
gorilla_plot = InputFile([path_inp_files,'/gorilla_plot.inp']);
tetra_grid = InputFile([path_inp_files,'/tetra_grid.inp']);

%read default inputfiles
gorilla.read();
gorilla_plot.read();
tetra_grid.read();


%Set input variables for GORILLA
%-----------------------------------------------------------------------------------------------------------------------

%Input file gorilla.inp
%All necessary variables for current calculation

    %Change in the electrostatic potential within the plasma volume in Gaussian units
        gorilla.GORILLANML.eps_Phi=0;

    %Coordinate system
        %1 ... (R,phi,Z) cylindrical coordinate system
        %2 ... (s,theta,phi) symmetry flux coordinate system
        gorilla.GORILLANML.coord_system = 2;

    %particle species
        %1 ... electron, 2 ... deuterium ion, 3 ... alpha particle
        gorilla.GORILLANML.ispecies = 2;

    %Switch for initial periodic coordinate particle re-location (modulo operation)
        % true ... Particles are re-located at initialization in the case of a periodic coordinate, if they are outside the computation domain.
        % false ... Particles are not re-located at initialization (This might lead to error if particles are outside the computation domain)
        gorilla.GORILLANML.boole_periodic_relocation = false;

    %1 ... numerical RK pusher, 2 ... polynomial pusher
        gorilla.GORILLANML.ipusher = 2;

    %Polynomial order for orbit pusher (from 2 to 4)
        gorilla.GORILLANML.poly_order = 2;


%Input file tetra_grid.inp
%All necessary variables for current calculation

    %Grid Size
        %Rectangular: nR, Field-aligned: ns
        tetra_grid.TETRA_GRID_NML.n1 = 100;
        %Rectangular: nphi, Field-aligned: nphi
        tetra_grid.TETRA_GRID_NML.n2 = 30;
        %Rectangular: nZ, Field-aligned: ntheta
        tetra_grid.TETRA_GRID_NML.n3 = 30;

    %Grid kind
        %1 ... rectangular grid for axisymmetric EFIT data
        %2 ... field-aligned grid for axisymmetric EFIT data
        %3 ... field-aligned grid for non-axisymmetric VMEC
        tetra_grid.TETRA_GRID_NML.grid_kind = 3;

    %Switch for selecting number of field periods automatically or manually
        %.true. ... number of field periods is selected automatically (Tokamak = 1, Stellarator depending on VMEC equilibrium)
        %.false. ... number of field periods is selected manually (see below)
        tetra_grid.TETRA_GRID_NML.boole_n_field_periods = true;


%Input file gorilla_plot.inp
%All necessary variables for current calculation

    %Switch for options
        % 1 ... Single orbit - Starting positions and pitch for the orbit are taken from file (see below) [First Line]
        % 2 ... Single orbit - Starting positions and pitch for the orbit are taken from starting drift surfaces (see below)
        % 3 ... Multiple orbits - Starting positions and pitch for orbits are taken from file (see below) [Every Line New Starting position]
        % 4 ... Multiple orbits - Starting positions and pitch for orbits are taken from drift surfaces with regular spacing (see below)
        gorilla_plot.GORILLA_PLOT_NML.i_orbit_options = 2;

    %Total individual orbit flight time for plotting
        gorilla_plot.GORILLA_PLOT_NML.total_orbit_time = 2;

    %Total Energy of particle in eV
        gorilla_plot.GORILLA_PLOT_NML.energy_eV_start = 3000;

    %Switch for plotting Poincar� cuts at toroidal variable $\varphi$ = 0
        gorilla_plot.GORILLA_PLOT_NML.boole_poincare_phi_0 = false;

    %Switch for plotting Poincar� cuts at parallel velocity $v_\parallel$ = 0
        gorilla_plot.GORILLA_PLOT_NML.boole_poincare_vpar_0 = true;

        %Number of skipped (non-printed) Poincar� cuts at parallel velocity $v_\parallel$ = 0
            gorilla_plot.GORILLA_PLOT_NML.n_skip_vpar_0 = 1;

        %Filename for Poincar� cuts at parallel velocity $v_\parallel$ = 0 in cylindrical coordinates (R,$\varphi$,Z)
            gorilla_plot.GORILLA_PLOT_NML.filename_poincare_vpar_0_rphiz = 'poincare_plot_vpar_0_rphiz_order2.dat';

        %Filename for Poincar� cuts at parallel velocity $v_\parallel$ = 0 in symmetry flux coordinates (s,$\vartheta$,$\varphi$)
            gorilla_plot.GORILLA_PLOT_NML.filename_poincare_vpar_0_sthetaphi = 'poincare_plot_vpar_0_sthetaphi_order2.dat';

    %Switch for plotting full orbit
        gorilla_plot.GORILLA_PLOT_NML.boole_full_orbit = false;

    %Plot invariances of motion (ONLY for single orbits)

        %Switch for plotting total particle energy
            gorilla_plot.GORILLA_PLOT_NML.boole_e_tot = false;

        %Switch for plotting canoncial (toroidal) angular momentum $p_\varphi$
            gorilla_plot.GORILLA_PLOT_NML.boole_p_phi = false;

        %Switch for parallel adiabatic invariant $J_\parallel$
            gorilla_plot.GORILLA_PLOT_NML.boole_J_par = true;

        %Filename for parallel adiabatic invariant $J_\parallel$
            gorilla_plot.GORILLA_PLOT_NML.filename_J_par = 'J_par_order2.dat';

    %Single orbit from starting drift surface (i_orbit_options = 2)

        %Starting drift surface
            % = s for (s,$\vartheta$,$\varphi$)
            % = R for (R,$\varphi$,Z)
            gorilla_plot.GORILLA_PLOT_NML.start_pos_x1_beg = 0.5;

        %Starting value for toroidal variable
            % = $\vartheta$ for (s,$\vartheta$,$\varphi$)
            % = $\varphi$ for (R,$\varphi$,Z)
            gorilla_plot.GORILLA_PLOT_NML.start_pos_x2  = 0;

        %Starting value for poloidal variable $\vartheta$
            % = $\varphi$ for (s,$\vartheta$,$\varphi$)
            % = Z for (R,$\varphi$,Z)
            gorilla_plot.GORILLA_PLOT_NML.start_pos_x3  = 0.63;

        %Pitch parameter $\lambda$ = $v_\parallel$ / vmod
            gorilla_plot.GORILLA_PLOT_NML.start_pitch_parameter = 0.3;


%Run GORILLA
%-----------------------------------------------------------------------------------------------------------------------

%write Input files for GORILLA
gorilla.write([path_RUN,'/gorilla.inp']);
gorilla_plot.write([path_RUN,'/gorilla_plot.inp']);
tetra_grid.write([path_RUN,'/tetra_grid.inp']);

%Create softlinks for used files
if isfile('../../../test_gorilla_main.x')
    ! ln -s ../../../test_gorilla_main.x .
elseif isfile('../../../BUILD/SRC/test_gorilla_main.x')
    ! ln -s ../../../BUILD/SRC/test_gorilla_main.x .
else
    disp('GORILLA not built, exiting the MatLab script')
    return
end
! ln -s ../../../MHD_EQUILIBRIA .

%Run GORILLA code
! ./test_gorilla_main.x

%Load Data-files and copy them into data_plots folder
data=struct();
data=load_copy(data,path_RUN,path_data_plots,gorilla_plot);

%Repeat orbit calculation for different polynominal orders and Runge-Kutta 4
%polynominal order = 3
gorilla.GORILLANML.poly_order = 3;
gorilla_plot.GORILLA_PLOT_NML.filename_J_par = 'J_par_order3.dat';
gorilla_plot.GORILLA_PLOT_NML.filename_poincare_vpar_0_rphiz = 'poincare_plot_vpar_0_rphiz_order3.dat';
gorilla_plot.GORILLA_PLOT_NML.filename_poincare_vpar_0_sthetaphi = 'poincare_plot_vpar_0_sthetaphi_order3.dat';
gorilla.write([path_RUN,'/gorilla.inp']);
gorilla_plot.write([path_RUN,'/gorilla_plot.inp']);
! ./test_gorilla_main.x
data=load_copy(data,path_RUN,path_data_plots,gorilla_plot);

%polynominal order = 4
gorilla.GORILLANML.poly_order = 4;
gorilla_plot.GORILLA_PLOT_NML.filename_J_par = 'J_par_order4.dat';
gorilla_plot.GORILLA_PLOT_NML.filename_poincare_vpar_0_rphiz = 'poincare_plot_vpar_0_rphiz_order4.dat';
gorilla_plot.GORILLA_PLOT_NML.filename_poincare_vpar_0_sthetaphi = 'poincare_plot_vpar_0_sthetaphi_order4.dat';
gorilla.write([path_RUN,'/gorilla.inp']);
gorilla_plot.write([path_RUN,'/gorilla_plot.inp']);
! ./test_gorilla_main.x
data=load_copy(data,path_RUN,path_data_plots,gorilla_plot);

%numerical RK4
gorilla.GORILLANML.ipusher = 1;
gorilla.GORILLANML.boole_pusher_ode45 = false;
gorilla_plot.GORILLA_PLOT_NML.filename_J_par = 'J_par_rk4.dat';
gorilla_plot.GORILLA_PLOT_NML.filename_poincare_vpar_0_rphiz = 'poincare_plot_vpar_0_rphiz_rk4.dat';
gorilla_plot.GORILLA_PLOT_NML.filename_poincare_vpar_0_sthetaphi = 'poincare_plot_vpar_0_sthetaphi_rk4.dat';
gorilla.write([path_RUN,'/gorilla.inp']);
gorilla_plot.write([path_RUN,'/gorilla_plot.inp']);
! ./test_gorilla_main.x
data=load_copy(data,path_RUN,path_data_plots,gorilla_plot);


%Create plots of generated data
%-----------------------------------------------------------------------------------------------------------------------

%Normalize toroidal angular momentum
data.J_par_order2(:,2)=data.J_par_order2(:,2)./data.J_par_order2(1,2);
data.J_par_order3(:,2)=data.J_par_order3(:,2)./data.J_par_order3(1,2);
data.J_par_order4(:,2)=data.J_par_order4(:,2)./data.J_par_order4(1,2);
data.J_par_rk4(:,2)=data.J_par_rk4(:,2)./data.J_par_rk4(1,2);

%Poincare plot of a trapped particle in VMEC equilbirum field of a stellarator
figure
subplot(1,3,1)
plot(data.poincare_plot_vpar_0_rphiz_order2(:,1),data.poincare_plot_vpar_0_rphiz_order2(:,3),'s','MarkerSize',2)
hold on
plot(data.poincare_plot_vpar_0_rphiz_order3(:,1),data.poincare_plot_vpar_0_rphiz_order3(:,3),'d','MarkerSize',2)
plot(data.poincare_plot_vpar_0_rphiz_order4(:,1),data.poincare_plot_vpar_0_rphiz_order4(:,3),'v','MarkerSize',2)
plot(data.poincare_plot_vpar_0_rphiz_rk4(:,1),data.poincare_plot_vpar_0_rphiz_rk4(:,3),'^','MarkerSize',2)
legend('GORILLA Poly2','GORILLA Poly3','GORILLA Poly4','GORILLA RK4')
xlabel('R [cm]')
ylabel('Z [cm]')
title('Stellarator: Poincar\.{e} section ($v_{\parallel}=0$)','Interpreter','latex');
hold off

subplot(1,3,2)
plot(data.poincare_plot_vpar_0_sthetaphi_order2(:,2),data.poincare_plot_vpar_0_sthetaphi_order2(:,1),'s','MarkerSize',2)
hold on
plot(data.poincare_plot_vpar_0_sthetaphi_order3(:,2),data.poincare_plot_vpar_0_sthetaphi_order3(:,1),'d','MarkerSize',2)
plot(data.poincare_plot_vpar_0_sthetaphi_order4(:,2),data.poincare_plot_vpar_0_sthetaphi_order4(:,1),'v','MarkerSize',2)
plot(data.poincare_plot_vpar_0_sthetaphi_rk4(:,2),data.poincare_plot_vpar_0_sthetaphi_rk4(:,1),'^','MarkerSize',2)
legend('GORILLA Poly2','GORILLA Poly3','GORILLA Poly4','GORILLA RK4')
xlabel('\vartheta')
ylabel('s')
title('Stellarator: Poincar\.{e} section ($v_{\parallel}=0$)','Interpreter','latex');
hold off

subplot(1,3,3)
plot(data.J_par_order2(:,1),data.J_par_order2(:,2),'s','MarkerSize',2)
hold on
plot(data.J_par_order3(:,1),data.J_par_order3(:,2),'d','MarkerSize',2)
plot(data.J_par_order4(:,1),data.J_par_order4(:,2),'v','MarkerSize',2)
plot(data.J_par_rk4(:,1),data.J_par_rk4(:,2),'^','MarkerSize',2)
legend('GORILLA Poly2','GORILLA Poly3','GORILLA Poly4','GORILLA RK4')
xlabel('N_{mappings}')
ylabel('$J_{\parallel}/J_{\parallel}(t=0)$','Interpreter','latex')
title('Stellarator: Normalized parallel adiabatic invariant','Interpreter','latex');
hold off

%Save figure as Matlab-file and as png
% savefig([path_data_plots,'/poincare_section.fig']);
% saveas(gcf,[path_data_plots,'/poincare_section.png']);

%Go back to path of Matlab script
cd(path_script);
