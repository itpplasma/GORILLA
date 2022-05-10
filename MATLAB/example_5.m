%#######################################################################################################################
% description:
%-----------------------------------------------------------------------------------------------------------------------
% * Compute collisionless guiding-center orbits with GORILLA for a passing and a trapped Deuterium particle.
% * Use a field-aligned grid for an axisymmetric tokamak equilibrium (g-file).
% * Plot the plasma boundary, the guiding-center orbits, and the resulting Poincare plot (\varphi = 0) for both orbits.
%
%#######################################################################################################################
% used functions:
%-----------------------------------------------------------------------------------------------------------------------
% InputFile
% NameList
% read_in
% load_copy
% efit
%#######################################################################################################################
% authors: Michael Scheidt, Michael Eder
% created: 19.02.2021


%Initialize used paths and files
%-----------------------------------------------------------------------------------------------------------------------

%Name of the current calculation to create folder
name_test_case='example_5';

%path of Matlab script
path_script=pwd;

%main path of GORILLA
c=strsplit(path_script,'/');
path_main=strjoin(c(1:end-1),'/');

%path to run GORILLA
mkdir(path_main,['EXAMPLES/MATLAB_RUN/',name_test_case]);
path_RUN=[path_main,'/EXAMPLES/MATLAB_RUN/',name_test_case];

%path of input files (blueprints)
path_inp_files=[path_main,'/INPUT'];

%path of the used functions
path_functions=[path_main,'/MATLAB/functions'];

%define path for data and plots and create a new folder their
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
%All variables get a default value from the blueprints
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
        gorilla.GORILLANML.boole_periodic_relocation = true;

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
        tetra_grid.TETRA_GRID_NML.n2 = 40;
        %Rectangular: nZ, Field-aligned: ntheta
        tetra_grid.TETRA_GRID_NML.n3 = 40;

    %Grid kind
        %1 ... rectangular grid for axisymmetric EFIT data
        %2 ... field-aligned grid for axisymmetric EFIT data
        %3 ... field-aligned grid for non-axisymmetric VMEC
        tetra_grid.TETRA_GRID_NML.grid_kind = 2;

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
        gorilla_plot.GORILLA_PLOT_NML.total_orbit_time = 0.004;

    %Total Energy of particle in eV
        gorilla_plot.GORILLA_PLOT_NML.energy_eV_start = 3000;

    %Switch for plotting Poincar� cuts at toroidal variable $\varphi$ = 0
        gorilla_plot.GORILLA_PLOT_NML.boole_poincare_phi_0 = true;

        %Number of skipped (non-printed) Poincar� cuts at toroidal variable $\varphi$ = 0
            gorilla_plot.GORILLA_PLOT_NML.n_skip_phi_0 = 1;

        %Filename for Poincar� cuts at toroidal variable $\varphi$ = 0 in cylindrical coordinates (R,$\varphi$,Z)
            gorilla_plot.GORILLA_PLOT_NML.filename_poincare_phi_0_rphiz = 'poincare_plot_phi_0_rphiz_passing.dat';

        %Filename for Poincar� cuts at toroidal variable $\varphi$ = 0 in symmetry flux coordinates (s,$\vartheta$,$\varphi$)
            gorilla_plot.GORILLA_PLOT_NML.filename_poincare_phi_0_sthetaphi = 'poincare_plot_phi_0_sthetaphi_passing.dat';

    %Switch for plotting Poincar� cuts at parallel velocity $v_\parallel$ = 0
        gorilla_plot.GORILLA_PLOT_NML.boole_poincare_vpar_0 = false;

    %Switch for plotting full orbit
        gorilla_plot.GORILLA_PLOT_NML.boole_full_orbit = true;

        %Number of skipped (non-printed tetrahedra passings) full orbit
            gorilla_plot.GORILLA_PLOT_NML.n_skip_full_orbit = 1;

        %Filename for full orbit in cylindrical coordinates (R,$\varphi$,Z)
            gorilla_plot.GORILLA_PLOT_NML.filename_full_orbit_rphiz = 'full_orbit_plot_rphiz_passing.dat';

        %Filename for full orbit in symmetry flux coordinates (s,$\vartheta$,$\varphi$)
            gorilla_plot.GORILLA_PLOT_NML.filename_full_orbit_sthetaphi = 'full_orbit_plot_sthetaphi_passing.dat';

    %Plot invariances of motion (ONLY for single orbits)

        %Switch for plotting total particle energy
            gorilla_plot.GORILLA_PLOT_NML.boole_e_tot = false;

        %Switch for plotting canoncial (toroidal) angular momentum $p_\varphi$
            gorilla_plot.GORILLA_PLOT_NML.boole_p_phi = false;

        %Switch for parallel adiabatic invariant $J_\parallel$
            gorilla_plot.GORILLA_PLOT_NML.boole_J_par = false;

    %Single orbit from starting drift surface (i_orbit_options = 2)

        %Starting drift surface
            % = s for (s,$\vartheta$,$\varphi$)
            % = R for (R,$\varphi$,Z)
            gorilla_plot.GORILLA_PLOT_NML.start_pos_x1_beg = 0.5;

        %End drift surface
            % = s for (s,$\vartheta$,$\varphi$)
            % = R for (R,$\varphi$,Z)
            gorilla_plot.GORILLA_PLOT_NML.start_pos_x1_end = 0.9;

        %Number of drift surfaces in between start and end
            gorilla_plot.GORILLA_PLOT_NML.n_surfaces = 30;

        %Starting value for toroidal variable
            % = $\vartheta$ for (s,$\vartheta$,$\varphi$)
            % = $\varphi$ for (R,$\varphi$,Z)
            gorilla_plot.GORILLA_PLOT_NML.start_pos_x2  = 0.1;

        %Starting value for poloidal variable $\vartheta$
            % = $\varphi$ for (s,$\vartheta$,$\varphi$)
            % = Z for (R,$\varphi$,Z)
            gorilla_plot.GORILLA_PLOT_NML.start_pos_x3  = 3.63;

        %Pitch parameter $\lambda$ = $v_\parallel$ / vmod
            gorilla_plot.GORILLA_PLOT_NML.start_pitch_parameter = 0.9;


%Run GORILLA
%-----------------------------------------------------------------------------------------------------------------------

%write Input files for Gorilla
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
! ln -s ../../../INPUT/field_divB0.inp .
! ln -s ../../../INPUT/preload_for_SYNCH.inp .

%Run GORILLA Code
! ./test_gorilla_main.x

%Load Data-files and copy them into data_plots folder
data=struct();
data=load_copy(data,path_RUN,path_data_plots,gorilla_plot);


%Second run for trapped particle
%-----------------------------------------------------------------------------------------------------------------------

%Changes in Input files
gorilla_plot.GORILLA_PLOT_NML.filename_poincare_phi_0_rphiz = 'poincare_plot_phi_0_rphiz_trapped.dat';
gorilla_plot.GORILLA_PLOT_NML.filename_poincare_phi_0_sthetaphi = 'poincare_plot_phi_0_sthetaphi_trapped.dat';
gorilla_plot.GORILLA_PLOT_NML.filename_full_orbit_rphiz = 'full_orbit_plot_rphiz_trapped.dat';
gorilla_plot.GORILLA_PLOT_NML.filename_full_orbit_sthetaphi = 'full_orbit_plot_sthetaphi_trapped.dat';
gorilla_plot.GORILLA_PLOT_NML.start_pitch_parameter = 0.4;

gorilla_plot.write([path_RUN,'/gorilla_plot.inp']);

%Run GORILLA code
! ./test_gorilla_main.x
data=load_copy(data,path_RUN,path_data_plots,gorilla_plot);


%Create plots of generated data
%-----------------------------------------------------------------------------------------------------------------------

%Names for particle calculation in plot
particle_type = {'Passing','Trapped'};

%Number of points shown from full_orbit_files [passing, trapped]
n_show_points = [290,300];

%Loop over GORILLA calculations
for i=1:2
    figure
    hold on

    %Data for passing or trapped paricle
    switch i
        case 1
            full_orbit_rphiz = data.full_orbit_plot_rphiz_passing;
            poincare_rphiz = data.poincare_plot_phi_0_rphiz_passing;
        case 2
            full_orbit_rphiz = data.full_orbit_plot_rphiz_trapped;
            poincare_rphiz = data.poincare_plot_phi_0_rphiz_trapped;
    end

    %XYZ Data of orbit for plot
    full_orbit_xyz = full_orbit_rphiz;
    full_orbit_xyz(:,1) = full_orbit_rphiz(:,1).*cos(full_orbit_rphiz(:,2));
    full_orbit_xyz(:,2) = full_orbit_rphiz(:,1).*sin(full_orbit_rphiz(:,2));

    poincare_xyz = poincare_rphiz;
    poincare_xyz(:,1) = poincare_rphiz(:,1).*cos(poincare_rphiz(:,2));
    poincare_xyz(:,2) = poincare_rphiz(:,1).*sin(poincare_rphiz(:,2));

    %Used points of full orbit
    full_orbit_xyz=full_orbit_xyz(end-n_show_points(i):end,:);

    %Limits of plot
    xlim([-max(full_orbit_rphiz(:,1))-30,max(full_orbit_rphiz(:,1))+30]);
    ylim([-max(full_orbit_rphiz(:,1))-30,max(full_orbit_rphiz(:,1))+30]);
    zlim([min(full_orbit_rphiz(:,3))-48,max(full_orbit_rphiz(:,3))+30]);

    %Full orbit posiitons
    p0 = plot3(full_orbit_xyz(:,1),full_orbit_xyz(:,2),full_orbit_xyz(:,3),'r-');
    p0.LineWidth = 2;

    %Current particle position
    p1 = scatter3(full_orbit_xyz(end,1),full_orbit_xyz(end,2),full_orbit_xyz(end,3),'r');
    p1.MarkerFaceColor = 'r';
    p1.DisplayName = [particle_type{i},' particle'];

    %Poincare plot
    p2 = scatter3(poincare_xyz(:,1),poincare_xyz(:,2),poincare_xyz(:,3),'b');
    p2.MarkerFaceColor = 'b';
    p2.SizeData = 6.0;
    p2.DisplayName = 'Poincar\.{e} plot';


    %Plasma boundaries
    e = efit([path_RUN,'/MHD_EQUILIBRIA/g_file_for_test'], [], []);
    e.read();

    [R, PHI] = meshgrid(e.rbbbs .* 100, linspace(0,1.5*pi,100));

    sb = surf(R.*cos(PHI), R.*sin(PHI), repmat(e.zbbbs.*100, numel(linspace(0,1.5*pi,100)),1));
    sb.FaceColor = [0,0.7,0.7];
    sb.EdgeColor = 'none';
    sb.FaceAlpha = 0.4;
    sb.DisplayName = 'Plasma boundary';
    sb.FaceLighting = 'gouraud';
    light

    %Plane of Poincare plot
    xp = [min(poincare_rphiz(:,1))-30,max(poincare_rphiz(:,1))+30];
    zp = [min(full_orbit_rphiz(:,3))-48,max(full_orbit_rphiz(:,3))+30];
    x1 = [ xp(1) xp(2) xp(2) xp(1)];
    y1 = zeros(size(x1));
    z1 = [ zp(1) zp(1) zp(2) zp(2)];

    v = patch(x1,y1,z1, 'b');
    v.FaceAlpha = 0.2;
    v.EdgeAlpha = 0.2;

    t = text(xp(2)-70,0,zp(1)+15,'$\varphi=0$','Interpreter','latex');
    t.FontSize = 15;
    t.Rotation = -14;

    %Legend
    lh = legend([sb,p1,p2],'Interpreter','latex');
    lh.Position = [0.2279 0.7798 0.2009 0.1452];
    lh.FontSize = 12;

    view(24.2848,60.7185);
    axis off
    hold off

    %Save figure as Matlab-file and as png
%     savefig([path_data_plots,'/poincare_plot_',particle_type{i},'.fig']);
%     saveas(gcf,[path_data_plots,'/poincare_plot_',particle_type{i},'.png']);
end

%Go back to path of Matlab script
cd(path_script);





