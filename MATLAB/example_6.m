%#######################################################################################################################
% description:
%-----------------------------------------------------------------------------------------------------------------------
% * Compute collisionless guiding-center orbits with GORILLA for passing and trapped Deuterium particles (WEST geometry).
% * Construct a 3D extension of the Soledge3x-EIRENE 2D-mesh for an axisymmetric tokamak equilibrium (g-file)
% * Plot the 2D projection of the guiding-center orbits on the original Soledge3x-EIRENE grid for both particle types.
%
%#######################################################################################################################
% used functions:
%-----------------------------------------------------------------------------------------------------------------------
% InputFile
% NameList
% read_in
% load_copy
%#######################################################################################################################
% authors: Georg Graßler
% created: 13.05.2022


%Initialize used paths and files
%-----------------------------------------------------------------------------------------------------------------------

%Name of the current calculation to create folder
name_test_case='example_6';

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
%All variables get a default value from the blueprints
gorilla.read();
gorilla_plot.read();
tetra_grid.read();


%Set input variables for GORILLA
%-----------------------------------------------------------------------------------------------------------------------

%Input file gorilla.inp
%All necessary variables for current calculation

    %Change in the electrostatic potential within the plasma volume in Gaussian units
        gorilla.GORILLANML.eps_Phi=-1e-7;

    %Coordinate system
        %1 ... (R,phi,Z) cylindrical coordinate system
        %2 ... (s,theta,phi) symmetry flux coordinate system
        gorilla.GORILLANML.coord_system = 1;

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
        gorilla.GORILLANML.poly_order = 4;


%Input file tetra_grid.inp
%All necessary variables for current calculation

    %Grid Size
        %Rectangular: nR, Field-aligned: ns
        tetra_grid.TETRA_GRID_NML.n1 = 100;
        %Rectangular: nphi, Field-aligned: nphi
        tetra_grid.TETRA_GRID_NML.n2 = 60;
        %Rectangular: nZ, Field-aligned: ntheta
        tetra_grid.TETRA_GRID_NML.n3 = 60;

    %Grid kind
        %1 ... rectangular grid for axisymmetric EFIT data
        %2 ... field-aligned grid for axisymmetric EFIT data
        %3 ... field-aligned grid for non-axisymmetric VMEC
        %4 ... SOLEDGE3X_EIRENE grid
        tetra_grid.TETRA_GRID_NML.grid_kind = 4;
        
    %MHD equilibrium filename
        tetra_grid.TETRA_GRID_NML.g_file_filename = 'MHD_EQUILIBRIA/g_file_for_test_WEST';
        tetra_grid.TETRA_GRID_NML.convex_wall_filename = 'MHD_EQUILIBRIA/convex_wall_for_test_WEST.dat';

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
        gorilla_plot.GORILLA_PLOT_NML.i_orbit_options = 3;

    %Total individual orbit flight time for plotting
        gorilla_plot.GORILLA_PLOT_NML.total_orbit_time = 0.004;

    %Total Energy of particle in eV
        gorilla_plot.GORILLA_PLOT_NML.energy_eV_start = 3000;

    %Switch for plotting Poincar� cuts at toroidal variable $\varphi$ = 0
        gorilla_plot.GORILLA_PLOT_NML.boole_poincare_phi_0 = false;

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
            
    %Filename for list of starting position(s) of particle(s) in cylindrical coordinates (R,$\varphi$,Z) and pitch ($\lambda$)
            gorilla_plot.GORILLA_PLOT_NML.filename_orbit_start_pos_rphiz = 'orbit_start_rphizlambda_passing.dat';

%Create inputfile for starting positions
%R,phi,Z,lambda
R_left_values = linspace(192.245,205.714,7);
Z_left_values = linspace(57.143,40.000,7);
R_inner_left_values = [207.959];
Z_inner_left_values = [37.143];
R_inner_right_values = [266.327,268.571];
Z_inner_right_values = [-37.143,-40.000];
R_right_values = linspace(270.816,291.020,10);
Z_right_values = linspace(-42.857,-68.571,10);
starting_2D = [R_left_values(:),Z_left_values(:); ...
                R_inner_left_values(:),Z_inner_left_values(:); ...
                R_inner_right_values(:),Z_inner_right_values(:); ...
                R_right_values(:),Z_right_values(:)];
lambda(1) = 0.9;
start_pos = zeros(length(starting_2D),4);
start_pos(:,[1,3]) = starting_2D;
start_pos(:,2) = 0;
start_pos(:,4) = lambda(1);

start_pos_print = start_pos.';
start_pos_print = start_pos_print(:).';

fileID = fopen([path_RUN,'/',gorilla_plot.GORILLA_PLOT_NML.filename_orbit_start_pos_rphiz],'w');
fprintf(fileID,'%.8f %.8f %.8f %.8f\n',start_pos_print);
fclose(fileID);


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
! ln -s ../../../INPUT/field_divB0.inp .
! ln -s ../../../INPUT/preload_for_SYNCH.inp .

%Run GORILLA Code
! ./test_gorilla_main.x

%Load Data-files and copy them into data_plots folder
data=struct();
data=load_copy(data,path_RUN,path_data_plots,gorilla_plot);

%Repeat orbit calculation for different pitch parameter (passing particles)
gorilla_plot.GORILLA_PLOT_NML.filename_full_orbit_rphiz = 'full_orbit_plot_rphiz_trapped.dat';
gorilla_plot.GORILLA_PLOT_NML.filename_full_orbit_sthetaphi = 'full_orbit_plot_sthetaphi_trapped.dat';
gorilla_plot.GORILLA_PLOT_NML.filename_orbit_start_pos_rphiz = 'orbit_start_rphizlambda_trapped.dat';

lambda(2) = 0.3;
start_pos(:,4) = lambda(2);
start_pos(end-10,[1,3]) = start_pos(end-10,[1,3]) + [0.5,-0.5];
start_pos(end-9,:) = [];
start_pos(8,:) = [];
start_pos_print = start_pos.';
start_pos_print = start_pos_print(:).';
fileID = fopen([path_RUN,'/',gorilla_plot.GORILLA_PLOT_NML.filename_orbit_start_pos_rphiz],'w');
fprintf(fileID,'%.8f %.8f %.8f %.8f\n',start_pos_print);
fclose(fileID);

gorilla.write([path_RUN,'/gorilla.inp']);
gorilla_plot.write([path_RUN,'/gorilla_plot.inp']);
! ./test_gorilla_main.x
data=load_copy(data,path_RUN,path_data_plots,gorilla_plot);


%Create plots of generated data
%-----------------------------------------------------------------------------------------------------------------------
%Names for guiding-center calculation in plot
guiding_center_type = {'Passing','Trapped'};

%Color scheme
grid_color = [204,204,204]/256;
grid_thickness = 0.5;
orbit_color = [4,90,141]/256;

%Read in Soledge3x-EIRENE mesh data and prepare 2D grid
fileID = fopen('MHD_EQUILIBRIA/MESH_SOLEDGE3X_EIRENE/knots_for_test.dat');
coordinates = textscan(fileID,'%f %f','HeaderLines',1);
fclose(fileID);
coordinates = cell2mat(coordinates);
n_vertex = length(coordinates);

fileID = fopen('MHD_EQUILIBRIA/MESH_SOLEDGE3X_EIRENE/triangles_for_test.dat');
triangles = textscan(fileID,'%f %f %f','HeaderLines',1);
fclose(fileID);
triangles = cell2mat(triangles);
n_triangles = length(triangles);

grid = nan*ones(n_triangles*5,2);
for t = [1:n_triangles]
    for k = [1:3]
        grid((t-1)*5 + k,:) = coordinates(triangles(t,k),:);
    end
    grid((t-1)*5 + 4,:) = coordinates(triangles(t,1),:);
    grid((t-1)*5 + 5,:) = nan;
end

%Loop over GORILLA calculations
for i=1:2
    figure('Renderer', 'painters', 'Position', [100 100 1300 1000])
    grid_plot = plot(grid(:,1),grid(:,2),'Color',grid_color);
    grid_plot.LineWidth = grid_thickness;
    grid_plot.DisplayName = 'SOLEDGE3X-EIRENE mesh';
    hold on

    %Data for passing or trapped paricle
    switch i
        case 1
            full_orbit_rphiz = data.full_orbit_plot_rphiz_passing;
            start_pos = load('orbit_start_rphizlambda_passing.dat');
        case 2
            full_orbit_rphiz = data.full_orbit_plot_rphiz_trapped;
            start_pos = load('orbit_start_rphizlambda_trapped.dat');
    end

    orbit_plot = plot(full_orbit_rphiz(:,1),full_orbit_rphiz(:,3),'.','Color',orbit_color);
    orbit_plot.DisplayName = 'Guiding-center orbit-projection';
    start_pos_plot = plot(start_pos(:,1),start_pos(:,3),'o','MarkerFaceColor',orbit_color,'MarkerEdgeColor',orbit_color);
    start_pos_plot.DisplayName = 'Guiding-center starting postions';
    
    %Legend
    lh = legend([grid_plot,orbit_plot,start_pos_plot],'Interpreter','latex');
    lh.FontSize = 12;
    lh.Location = 'southwest';

    hold off
    xlabel('R [cm]')
    ylabel('Z [cm]')
    title(['SOLEDGE3X-EIRENE: Toroidal projection ($\varphi=0$) of guiding-center orbits ($\lambda =$ ',num2str(lambda(i)),')'],'Interpreter','latex');

    %Save figure as Matlab-file and as png
    %savefig([path_data_plots,'/orbit_projections_',guiding_center_type{i},'.fig']);
    saveas(gcf,[path_data_plots,'/orbit_projections_',guiding_center_type{i},'.png']);
    
end


%Go back to path of Matlab script
cd(path_script);