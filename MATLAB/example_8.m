%#######################################################################################################################
% description:
%-----------------------------------------------------------------------------------------------------------------------
% * Using a potential that is a scaled poloidal flux (ASDEX-mesh) calculating the corresponding electric field via central differences.
% * Perform each a run with and without the inclusion of the additional terms in the strong electric field Lagrangian using fully ionised Tungsten.
% * Plot the poincare cross-section of both runs, as well as their respective fluctuations of total energy and toroidal momentum.
% * The definition of total energy and toroidal momentum is different in the two runs. Each definition should yield conserved quantities in their respective case.
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
% created: 04.03.2023


%Initialize used paths and files
%-----------------------------------------------------------------------------------------------------------------------

%Name of the current calculation to create folder
name_test_case='example_8';

%Control quantitites
%close all

lambda(1) = -0.4; %0.3 for WEST / -0.9 for ASDEX
eps_Phi = 2*1e-7; %-4*1e-4 for WEST / +1e-4 for ASDEX / +1e-5 ASDEX Tungsten
%eps_Phi = 1e-7;
energy_eV_start = 3e3; %3e3
time_multiplier = 100;
boole_strong_electric_field = true;
boole_calculate_conserved_quantities = false;
chosen_orbit_index =10; %8/11 for strong_electric else 8/10
grid_kind = 2;
ipusher = 2;
ispecies = 4;
start_ASDEX = [210,0];

standart = 1e-7;
standart_energy = 3e3;
standart_lambda = 0.3;
%total_orbit_time = 0.002*standart/max(abs(eps_Phi),standart)*sqrt(standart_energy/energy_eV_start);
total_orbit_time = 2e-5;

if (boole_calculate_conserved_quantities)
    total_orbit_time = time_multiplier*total_orbit_time;
    n_skip_full_orbit = time_multiplier;
    i_orbit_options = 1;
else
    n_skip_full_orbit = 1;
    i_orbit_options = 3;
end

if (grid_kind == 2)
    energy_eV_start = 3e5;
    g_file_filename = 'MHD_EQUILIBRIA/g_file_for_test';
    convex_wall_filename = 'MHD_EQUILIBRIA/convex_wall_for_test.dat';
else
    g_file_filename = 'MHD_EQUILIBRIA/g_file_for_test_WEST';
    convex_wall_filename = 'MHD_EQUILIBRIA/convex_wall_for_test_WEST.dat';
end

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
        gorilla.GORILLANML.eps_Phi= eps_Phi; %-1e-7

    %Coordinate system
        %1 ... (R,phi,Z) cylindrical coordinate system
        %2 ... (s,theta,phi) symmetry flux coordinate system
        gorilla.GORILLANML.coord_system = 1;

    %particle species
        %1 ... electron, 2 ... deuterium ion, 3 ... alpha particle, 4 ... ionised tungsten
        gorilla.GORILLANML.ispecies = ispecies;

    %Switch for initial periodic coordinate particle re-location (modulo operation)
        % true ... Particles are re-located at initialization in the case of a periodic coordinate, if they are outside the computation domain.
        % false ... Particles are not re-located at initialization (This might lead to error if particles are outside the computation domain)
        gorilla.GORILLANML.boole_periodic_relocation = true;

    %1 ... numerical RK pusher, 2 ... polynomial pusher
        gorilla.GORILLANML.ipusher = ipusher;

    %Polynomial order for orbit pusher (from 2 to 4)
        gorilla.GORILLANML.poly_order = 3;
        
    %Adaptive step size scheme to enforce energy conservation
        gorilla.GORILLANML.boole_adaptive_time_steps = false;
        gorilla.GORILLANML.desired_delta_energy = 1.E-14;
       
    %In the guiding center theory of the default implementation the electric field is assumed to be weak 
    %compared to the magnetic field. In case of the investigation of impurities in WEST geometry
    %(Soledge3X-EIRENE, cylindrical coordinates only) additional terms have to be included in the dynamics
    
        %Switch for including strong electric field terms
            %false ... no additional terms
            %true ... includes polarization drift terms in dynamics
        gorilla.GORILLANML.boole_strong_electric_field = boole_strong_electric_field;
        
        %Switch for saving electric field and drift to .dat-files
        gorilla.GORILLANML.boole_save_electric = boole_strong_electric_field;

        %Filenames for storing electric field and drift velocity
        gorilla.GORILLANML.filename_electric_field = 'electric_field.dat';
        gorilla.GORILLANML.filename_electric_drift = 'electric_drift.dat';


%Input file tetra_grid.inp
%All necessary variables for current calculation

    %Grid Size
        %Rectangular: nR, Field-aligned: ns
        tetra_grid.TETRA_GRID_NML.n1 = 10;
        %Rectangular: nphi, Field-aligned: nphi
        tetra_grid.TETRA_GRID_NML.n2 = 10;
        %Rectangular: nZ, Field-aligned: ntheta
        tetra_grid.TETRA_GRID_NML.n3 = 120;

    %Grid kind
        %1 ... rectangular grid for axisymmetric EFIT data
        %2 ... field-aligned grid for axisymmetric EFIT data
        %3 ... field-aligned grid for non-axisymmetric VMEC
        %4 ... SOLEDGE3X_EIRENE grid
        tetra_grid.TETRA_GRID_NML.grid_kind = grid_kind;
        
    %Option for $\theta$-variable being used in grid
        % 1 ... theta scaling in symmetry flux coordinates
        % 2 ... theta scaling in geometrical theta angle
        tetra_grid.TETRA_GRID_NML.theta_geom_flux = 2;
        
    %MHD equilibrium filename
        tetra_grid.TETRA_GRID_NML.g_file_filename = g_file_filename;
        tetra_grid.TETRA_GRID_NML.convex_wall_filename = convex_wall_filename;

    %Switch for selecting number of field periods automatically or manually
        %.true. ... number of field periods is selected automatically (Tokamak = 1, Stellarator depending on VMEC equilibrium)
        %.false. ... number of field periods is selected manually (see below)
        tetra_grid.TETRA_GRID_NML.boole_n_field_periods = true;
        
    %Switch for writing object file with mesh data
        tetra_grid.TETRA_GRID_NML.boole_write_mesh_obj = true;


%Input file gorilla_plot.inp
%All necessary variables for current calculation

    %Switch for options
        % 1 ... Single orbit - Starting positions and pitch for the orbit are taken from file (see below) [First Line]
        % 2 ... Single orbit - Starting positions and pitch for the orbit are taken from starting drift surfaces (see below)
        % 3 ... Multiple orbits - Starting positions and pitch for orbits are taken from file (see below) [Every Line New Starting position]
        % 4 ... Multiple orbits - Starting positions and pitch for orbits are taken from drift surfaces with regular spacing (see below)
        gorilla_plot.GORILLA_PLOT_NML.i_orbit_options = i_orbit_options;

    %Total individual orbit flight time for plotting
        gorilla_plot.GORILLA_PLOT_NML.total_orbit_time = total_orbit_time;

    %Total Energy of particle in eV
        gorilla_plot.GORILLA_PLOT_NML.energy_eV_start = energy_eV_start; %3000

    %Switch for plotting Poincar� cuts at toroidal variable $\varphi$ = 0
        gorilla_plot.GORILLA_PLOT_NML.boole_poincare_phi_0 = false;

    %Switch for plotting Poincar� cuts at parallel velocity $v_\parallel$ = 0
        gorilla_plot.GORILLA_PLOT_NML.boole_poincare_vpar_0 = false;

    %Switch for plotting full orbit
        gorilla_plot.GORILLA_PLOT_NML.boole_full_orbit = true;

        %Number of skipped (non-printed tetrahedra passings) full orbit
            gorilla_plot.GORILLA_PLOT_NML.n_skip_full_orbit = n_skip_full_orbit;

        %Filename for full orbit in cylindrical coordinates (R,$\varphi$,Z)
            gorilla_plot.GORILLA_PLOT_NML.filename_full_orbit_rphiz = 'full_orbit_plot_rphiz_electric.dat';

        %Filename for full orbit in symmetry flux coordinates (s,$\vartheta$,$\varphi$)
            gorilla_plot.GORILLA_PLOT_NML.filename_full_orbit_sthetaphi = 'full_orbit_plot_sthetaphi_electric.dat';

    %Plot invariances of motion (ONLY for single orbits)

        %Switch for plotting total particle energy
            gorilla_plot.GORILLA_PLOT_NML.boole_e_tot = true;
            
        %Filename for total energy
            gorilla_plot.GORILLA_PLOT_NML.filename_e_tot = 'e_tot_electric.dat';

        %Switch for plotting canoncial (toroidal) angular momentum $p_\varphi$
            gorilla_plot.GORILLA_PLOT_NML.boole_p_phi = true;
            
        %Filename for canoncial (toroidal) angular momentum $p_\varphi$
            gorilla_plot.GORILLA_PLOT_NML.filename_p_phi = 'p_phi_electric.dat';

        %Switch for parallel adiabatic invariant $J_\parallel$
            gorilla_plot.GORILLA_PLOT_NML.boole_J_par = false;
            
    %Filename for list of starting position(s) of particle(s) in cylindrical coordinates (R,$\varphi$,Z) and pitch ($\lambda$)
            gorilla_plot.GORILLA_PLOT_NML.filename_orbit_start_pos_rphiz = 'orbit_start_rphizlambda_electric.dat';

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
% R_left_values = [];
% Z_left_values = [];
% R_inner_left_values = [];
% Z_inner_left_values = [];
% R_inner_right_values(1) = [];
% Z_inner_right_values(1) = [];
% R_right_values = [];
% Z_right_values = [];
starting_2D = [R_left_values(:),Z_left_values(:); ...
                R_inner_left_values(:),Z_inner_left_values(:); ...
                R_inner_right_values(:),Z_inner_right_values(:); ...
                R_right_values(:),Z_right_values(:)];
if (boole_calculate_conserved_quantities)
    starting_2D = starting_2D(chosen_orbit_index,:);
end
if (grid_kind == 2)
    starting_2D = start_ASDEX;
end
start_pos = zeros(size(starting_2D,1),4);
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


%% Create plots of generated data
%-----------------------------------------------------------------------------------------------------------------------

%Color scheme
grid_thickness = 0.2;
electric_thickness = 0.2;
drift_thickness = 0.2;

grid_color = [204,204,204]/256;
electric_color = [256,0,0]/256;
drift_color = [0,0,256]/256;

orbit_color = [4,90,141]/256;
FontSize = 12;

%% Setup mesh for plotting
cd(path_RUN);

if (grid_kind == 2)
    
    %Read in ASDEX 2D mesh data and prepare 2D grid
    fileID = fopen('mesh_rphiz.obj');
    FILE = textscan(fileID,'%s %f %f %f');
    fclose(fileID);
    TEXT = FILE{1};
    DATA = cell2mat(FILE(2:4));
    n_vertex = sum(strcmp(TEXT,'v'));

    %Get front faces
    faces = DATA(n_vertex + 1:end,:);
    faces = faces(4:12:end,:);
    n_faces = size(faces,1);

    %Get front points
    n_vertex = n_vertex/2;
    vertex = DATA(1:n_vertex,[1,3]);

    grid = nan*ones(n_faces*5,2);
    for t = [1:n_faces]
        for k = 1:4
            grid((t-1)*5 + k,:) = vertex(faces(t,mod(k-1,3) + 1),:);
        end
        grid((t-1)*5 + 5,:) = nan;
    end

    n_s = tetra_grid.TETRA_GRID_NML.n1 + 1;
    n_theta = tetra_grid.TETRA_GRID_NML.n3;

    skip_theta = 1;
    skip_s = 1;
    
elseif (grid_kind == 4)
    
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
    
else
    error('Invalid grid_kind! Choose among 2 and 4!');
end

%Go back to path of Matlab script
cd(path_script);

%% Draw Poincare plot
cd(path_RUN);

if (~boole_calculate_conserved_quantities)
    
    %Draw grid
    figure('Renderer', 'painters', 'Position', [100 100 1300 1000])
    grid_plot = plot(grid(:,1),grid(:,2),'Color',grid_color);
    grid_plot.LineWidth = grid_thickness;
    grid_plot.DisplayName = 'SOLEDGE3X-EIRENE mesh';
    hold on

    %Data for orbits
    full_orbit_rphiz = data.full_orbit_plot_rphiz_electric;
    start_pos = load('orbit_start_rphizlambda_electric.dat');

    orbit_plot = plot(full_orbit_rphiz(:,1),full_orbit_rphiz(:,3),'.','Color',orbit_color);
    orbit_plot.DisplayName = 'Guiding-center orbit projection';
    start_pos_plot = plot(start_pos(:,1),start_pos(:,3),'o','MarkerFaceColor',orbit_color,'MarkerEdgeColor',orbit_color);
    start_pos_plot.DisplayName = 'Guiding-center starting postions';

    hold off

    %Legend
    lh = legend([grid_plot,orbit_plot,start_pos_plot],'Interpreter','latex');
    lh.FontSize = FontSize;
    lh.Location = 'southwest';

    xlabel('R [cm]')
    ylabel('Z [cm]')
    title(['SOLEDGE3X-EIRENE: Toroidal projection ($\varphi=0$) of electric field'],'Interpreter','latex');

    %Save figure as Matlab-file and as png
    %savefig([path_data_plots,'/orbit_projections_',guiding_center_type{i},'.fig']);
    saveas(gcf,[path_data_plots,'/orbit_projections_','conservation_electric','.png']);
    
else

    %Quantity data after each n_skip-th push
    e_tot_data = data.e_tot_electric;
    p_phi_data = data.p_phi_electric;

    %Energy / toroidal momentum evolution over time
    xlabel_txt = '$t$';
    ylabel_txt = '$E/E(t = 0) - 1$';

    figure('Renderer', 'painters', 'Position', [24 37 1214 580])
    
    subplot(1,2,1)
    plot(e_tot_data(1:end,1),e_tot_data(1:end,2)/e_tot_data(1,2)-1,'Color',orbit_color)
    tit = title(['Energy evolution'],'Interpreter','latex');
    xlab = xlabel(xlabel_txt,'Interpreter','latex');
    ylab = ylabel(ylabel_txt,'Interpreter','latex');
    tit.FontSize = FontSize;
    xlab.FontSize = FontSize;
    ylab.FontSize = FontSize;
    
    ylabel_txt = '$p_{\varphi}/p_{\varphi}(t = 0) - 1$';
    
    subplot(1,2,2)
    plot(p_phi_data(1:end,1),p_phi_data(1:end,2)/p_phi_data(1,2)-1,'Color',orbit_color)
    tit = title(['Toroidal momentum evolution'],'Interpreter','latex');
    xlab = xlabel(xlabel_txt,'Interpreter','latex');
    ylab = ylabel(ylabel_txt,'Interpreter','latex');
    tit.FontSize = FontSize;
    xlab.FontSize = FontSize;
    ylab.FontSize = FontSize;

end

axis equal

%Go back to path of Matlab script
cd(path_script);

%% Electric field
cd(path_RUN);

%Extracting an easy to plot projection of the electric field
arrow_scale_field = (2*1e-5/eps_Phi)*1e+0;
electric_field_data = load('electric_field.dat');

switch(grid_kind)
    case(2)
        electric_field = nan*ones(ceil(n_s/skip_s)*ceil((n_theta/skip_theta)),2);
        counter = 0;
        for i_s = [1:skip_s:n_s]
            for i_theta = [1:skip_theta:n_theta]
                t = (i_s-1)*n_theta + i_theta;
                counter = counter + 1;
                electric_field((counter-1)*3 + 1,:) = electric_field_data(t,[1,3]);
                electric_field((counter-1)*3 + 2,:) = electric_field_data(t,[1,3]) + arrow_scale_field*electric_field_data(t,[4,6]);
                electric_field((counter-1)*3 + 3,:) = nan;
            end
        end
    case(4)
        electric_field = nan*ones(n_vertex*3,2);
        for t = [1:n_vertex]
            electric_field((t-1)*3 + 1,:) = electric_field_data(t,[1,3]);
            electric_field((t-1)*3 + 2,:) = electric_field_data(t,[1,3]) + arrow_scale_field*electric_field_data(t,[4,6]);
            electric_field((t-1)*3 + 3,:) = nan;
        end
end

electric_mod = sqrt(electric_field_data(:,4).^2 + (electric_field_data(:,1).*electric_field_data(:,5)).^2 + electric_field_data(:,6).^2);
average_electric_mod = sum(electric_mod,1)/size(electric_mod,1);

%Extracting an easy to plot ExB drift on first plane
arrow_scale_drift = (2*1e-5/eps_Phi)*1e-6;
electric_drift_data = load('electric_drift.dat');
R = electric_drift_data(:,1);
phi = electric_drift_data(:,2);
Z = electric_drift_data(:,3);
v_E = electric_drift_data(:,4:6);

switch(grid_kind)
    case(2)
        electric_drift_field = nan*ones(ceil(n_s/skip_s)*ceil((n_theta/skip_theta)),3);
        counter = 0;
        for i_s = [1:skip_s:n_s]
            for i_theta = [1:skip_theta:n_theta]
                t = (i_s-1)*n_theta + i_theta;
                counter = counter + 1;
                electric_drift_field((counter-1)*3 + 1,:) = [R(t)*cos(phi(t)),R(t)*sin(phi(t)),Z(t)];
                electric_drift_field((counter-1)*3 + 2,:) = [R(t)*cos(phi(t)),R(t)*sin(phi(t)),Z(t)] + arrow_scale_drift*cylinder2cart(v_E(t,:),R(t),phi(t));
                electric_drift_field((counter-1)*3 + 3,:) = nan;
            end
        end
    case(4)
        electric_drift_field = nan*ones(n_vertex*3,3);
        for t = [1:n_vertex]
            electric_drift_field((t-1)*3 + 1,:) = [R(t)*cos(phi(t)),R(t)*sin(phi(t)),Z(t)];
            electric_drift_field((t-1)*3 + 2,:) = [R(t)*cos(phi(t)),R(t)*sin(phi(t)),Z(t)] + arrow_scale_drift*cylinder2cart(v_E(t,:),R(t),phi(t));
            electric_drift_field((t-1)*3 + 3,:) = nan;
        end
end

v2_E_mod = electric_drift_data(:,7);
v2_E_mod_check = v_E(:,1).^2 + v_E(:,2).^2./R.^2 + v_E(:,3).^2;

% Plotting the field plus drift
% figure('Renderer', 'painters', 'Position', [100 100 1300 1000])
% 
% electric_drift_plot = plot3(electric_drift_field(:,1),electric_drift_field(:,2),electric_drift_field(:,3),'-','Color',drift_color);
% electric_drift_plot.LineWidth = drift_thickness;
% electric_drift_plot.DisplayName = ['Electric drift (ratio 1:1e',num2str(round(log10(1/arrow_scale_drift))),')'];
% 
% hold on
% 
% grid_drift_plot = plot3(grid(:,1),zeros(size(grid(:,1))),grid(:,2),'Color',grid_color);
% grid_drift_plot.LineWidth = grid_thickness;
% grid_drift_plot.DisplayName = 'ASDEX 2D mesh';
% 
% electric_plot = plot3(electric_field(:,1),zeros(size(electric_field(:,1))),electric_field(:,2),'-','Color',electric_color);
% electric_plot.LineWidth = electric_thickness;
% electric_plot.DisplayName = 'Electric field projection';
% 
% hold off
% 
% axis equal
% 
% %Legend
% lh_drift = legend([grid_drift_plot,electric_plot,electric_drift_plot],'Interpreter','latex');
% lh_drift.FontSize = FontSize;
% lh_drift.Location = 'southwest';
    

figure('Renderer', 'painters', 'Position', [100 100 1000 1000])

electric_drift_plot = plot(electric_drift_field(:,1),electric_drift_field(:,3),'-','Color',drift_color);
electric_drift_plot.LineWidth = drift_thickness;
electric_drift_plot.DisplayName = ['Electric drift (ratio 1:1e',num2str(round(log10(1/arrow_scale_drift))),')'];

hold on

grid_drift_plot = plot(grid(:,1),grid(:,2),'Color',grid_color);
grid_drift_plot.LineWidth = grid_thickness;
grid_drift_plot.DisplayName = 'ASDEX 2D mesh';

electric_plot = plot(electric_field(:,1),electric_field(:,2),'-','Color',electric_color);
electric_plot.LineWidth = electric_thickness;
electric_plot.DisplayName = 'Electric field projection';

hold off

axis equal

E_vec = electric_field_data(:,[4:6]);
E_mod = sqrt(sum(E_vec.^2,2));
v_E_mod = sqrt(v2_E_mod);
dot_product = sum(E_vec.*v_E,2)./(E_mod.*v_E_mod);

%Go back to path of Matlab script
cd(path_script);

function [v_cart] = cylinder2cart(v_cylinder,r,phi)
    v_r = v_cylinder(1);    
    v_phi = v_cylinder(2);
    v_z = v_cylinder(3);
    
    v_cart = nan*ones(1,3);
    v_cart(1) = v_r*cos(phi) - v_phi/r*sin(phi);
    v_cart(2) = v_r*sin(phi) + v_phi/r*cos(phi);
    v_cart(3) = v_z;
end