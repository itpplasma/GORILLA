%#######################################################################################################################
% description:
%-----------------------------------------------------------------------------------------------------------------------
% * Using a potential that is a scaled poloidal flux (ASDEX/SOLEDGE3X-mesh) calculating the corresponding electric field via central differences.
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


%Initialize test_case, used paths and files
%-----------------------------------------------------------------------------------------------------------------------

%Name of the current calculation to create folder
name_test_case='example_8';

%Control quantitites
close all

%lambda(1) = -0.9; %0.3 for WEST / -0.9 for ASDEX -0.4
%eps_Phi = 2*1e-5; %-4*1e-4 for WEST / +1e-4 for ASDEX / +1e-5 ASDEX Tungsten
%eps_Phi = 1e-7;

grid_kind = 2;

%Setting configurations depending on test case
switch(grid_kind)
    case(2)
        eps_Phi = 2*1e-5;
        starting_2D = [210,8];
        lambda = -0.9;
        n2 = 40;
        energy_eV_start = 3e5;
        total_orbit_time = 2e-2;
        g_file_filename = 'MHD_EQUILIBRIA/g_file_for_test';
        convex_wall_filename = 'MHD_EQUILIBRIA/convex_wall_for_test.dat';
    case(4)
        eps_Phi = -4*1e-4;
        starting_2D = [210,38];
        lambda = 0.3;
        n2 = 10;
        energy_eV_start = 3e5;
        total_orbit_time = 2e-2;
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
        gorilla.GORILLANML.eps_Phi= eps_Phi;

    %Coordinate system
        %1 ... (R,phi,Z) cylindrical coordinate system
        %2 ... (s,theta,phi) symmetry flux coordinate system
        gorilla.GORILLANML.coord_system = 1;

    %particle species
        %1 ... electron, 2 ... deuterium ion, 3 ... alpha particle, 4 ... ionised tungsten
        gorilla.GORILLANML.ispecies = 4;

    %Switch for initial periodic coordinate particle re-location (modulo operation)
        % true ... Particles are re-located at initialization in the case of a periodic coordinate, if they are outside the computation domain.
        % false ... Particles are not re-located at initialization (This might lead to error if particles are outside the computation domain)
        gorilla.GORILLANML.boole_periodic_relocation = true;

    %1 ... numerical RK pusher, 2 ... polynomial pusher
        gorilla.GORILLANML.ipusher = 2;

    %Polynomial order for orbit pusher (from 2 to 4)
        gorilla.GORILLANML.poly_order = 4;
        
    %Adaptive step size scheme to enforce energy conservation
        gorilla.GORILLANML.boole_adaptive_time_steps = false;
        
    %Allowed relative fluctuation of energy between entry and exit of a tetrahedron
    %Must not be smaller or equal zero!
        gorilla.GORILLANML.desired_delta_energy = 1.E-10;
       
    %In the guiding center theory of the default implementation the electric field is assumed to be weak 
    %compared to the magnetic field. In case of the investigation of impurities in WEST geometry
    %(Soledge3X-EIRENE, cylindrical coordinates only) additional terms have to be included in the dynamics
    
        %Switch for including strong electric field terms
            %false ... no additional terms
            %true ... includes polarization drift terms in dynamics
        gorilla.GORILLANML.boole_strong_electric_field = true;
        
        %Switch for saving electric field and drift to .dat-files
        gorilla.GORILLANML.boole_save_electric = true;

        %Filenames for storing electric field and drift velocity
        gorilla.GORILLANML.filename_electric_field = 'electric_field.dat';
        gorilla.GORILLANML.filename_electric_drift = 'electric_drift.dat';


%Input file tetra_grid.inp
%All necessary variables for current calculation

    %Grid Size
        %Rectangular: nR, Field-aligned: ns
        tetra_grid.TETRA_GRID_NML.n1 = 40;
        %Rectangular: nphi, Field-aligned: nphi
        tetra_grid.TETRA_GRID_NML.n2 = n2;
        %Rectangular: nZ, Field-aligned: ntheta
        tetra_grid.TETRA_GRID_NML.n3 = 40;

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
        
    %Filename for mesh object file in cylindrical coordintes
        tetra_grid.TETRA_GRID_NML.filename_mesh_rphiz = 'mesh_rphiz.obj';


%Input file gorilla_plot.inp
%All necessary variables for current calculation

    %Switch for options
        % 1 ... Single orbit - Starting positions and pitch for the orbit are taken from file (see below) [First Line]
        % 2 ... Single orbit - Starting positions and pitch for the orbit are taken from starting drift surfaces (see below)
        % 3 ... Multiple orbits - Starting positions and pitch for orbits are taken from file (see below) [Every Line New Starting position]
        % 4 ... Multiple orbits - Starting positions and pitch for orbits are taken from drift surfaces with regular spacing (see below)
        gorilla_plot.GORILLA_PLOT_NML.i_orbit_options = 1;

    %Total individual orbit flight time for plotting
        gorilla_plot.GORILLA_PLOT_NML.total_orbit_time = total_orbit_time;

    %Total Energy of particle in eV
        gorilla_plot.GORILLA_PLOT_NML.energy_eV_start = energy_eV_start;

    %Plot Poincaré cuts ($\varphi$ = 0)
        %Switch for plotting Poincaré cuts at toroidal variable $\varphi$ = 0
        gorilla_plot.GORILLA_PLOT_NML.boole_poincare_phi_0 = true;

        %Number of skipped (non-printed) Poincaré cuts at toroidal variable $\varphi$ = 0
        gorilla_plot.GORILLA_PLOT_NML.n_skip_phi_0 = 1;

        %Filename for Poincaré cuts at toroidal variable $\varphi$ = 0 in cylindrical coordinates (R,$\varphi$,Z)
        gorilla_plot.GORILLA_PLOT_NML.filename_poincare_phi_0_rphiz = 'poincare_plot_phi_0_rphiz_electric.dat';

    %Switch for plotting Poincar� cuts at parallel velocity $v_\parallel$ = 0
        gorilla_plot.GORILLA_PLOT_NML.boole_poincare_vpar_0 = false;

    %Switch for plotting full orbit
        gorilla_plot.GORILLA_PLOT_NML.boole_full_orbit = false;

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

start_pos_print = [starting_2D(1),0,starting_2D(2),lambda];
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

%Changes for RUN without considering strong electric field effects
gorilla.GORILLANML.boole_strong_electric_field = false;
gorilla.GORILLANML.boole_save_electric = false;
tetra_grid.TETRA_GRID_NML.boole_write_mesh_obj = false;
gorilla_plot.GORILLA_PLOT_NML.filename_poincare_phi_0_rphiz = 'poincare_plot_phi_0_rphiz.dat';
gorilla_plot.GORILLA_PLOT_NML.filename_e_tot = 'e_tot.dat';
gorilla_plot.GORILLA_PLOT_NML.filename_p_phi = 'p_phi.dat';
gorilla.write([path_RUN,'/gorilla.inp']);
gorilla_plot.write([path_RUN,'/gorilla_plot.inp']);

%Run without considering strong electric field effects
! ./test_gorilla_main.x

data=load_copy(data,path_RUN,path_data_plots,gorilla_plot);


%% Create plots of generated data
%-----------------------------------------------------------------------------------------------------------------------
cd(path_RUN);

%Color scheme and plotting settings
grid_thickness = 0.2;
electric_thickness = 0.2;
drift_thickness = 0.2;

grid_color = [204,204,204]/256;
electric_field_color = [256,0,0]/256;
electric_drift_color = [0,0,256]/256;

color_electric = [217,95,2]/256;
color = [27,158,119]/256;
color_start = [0,0,0]/256;

FontSize = 12;


%Setup of data for plotting
%-----------------------------------------------------------------------------------------------------------------------

%Setup mesh for plotting
fileID = fopen('mesh_rphiz.obj');
FILE = textscan(fileID,'%s %f %f %f');
fclose(fileID);
TEXT = FILE{1};
DATA = cell2mat(FILE(2:4));
n_vertex = sum(strcmp(TEXT,'v'));
faces = DATA(n_vertex + 1:end,:);
faces = faces(4:12:end,:);
n_faces = size(faces,1);
n_vertex = n_vertex/2;
vertex = DATA(1:n_vertex,[1,3]);

grid = nan*ones(n_faces*5,2);
for t = [1:n_faces]
    for k = 1:4
        grid((t-1)*5 + k,:) = vertex(faces(t,mod(k-1,3) + 1),:);
    end
    grid((t-1)*5 + 5,:) = nan;
end

switch(grid_kind)
    case(2)
        grid_name = 'Field-aligned ASDEX mesh';
        n_s = tetra_grid.TETRA_GRID_NML.n1 + 1;
        n_theta = tetra_grid.TETRA_GRID_NML.n3;
        skip_s = 2;
        skip_theta = 1;
    case(4)
        grid_name = 'SOLEDGE3X-EIRENE mesh';
        skip = 1;
end

%Setup poincare data for plotting
poincare_rphiz_electric = data.poincare_plot_phi_0_rphiz_electric;
poincare_rphiz = data.poincare_plot_phi_0_rphiz;
start_pos = load('orbit_start_rphizlambda_electric.dat');

%Setup energy/toroidal momentum for plotting
e_tot_electric = data.e_tot_electric;
p_phi_electric = data.p_phi_electric;
e_tot = data.e_tot;
p_phi = data.p_phi;
e_tot_electric(:,2) = e_tot_electric(:,2)/e_tot_electric(1,2) - 1;
p_phi_electric(:,2) = p_phi_electric(:,2)/p_phi_electric(1,2) - 1;
e_tot(:,2) = e_tot(:,2)/e_tot(1,2) - 1;
p_phi(:,2) = p_phi(:,2)/p_phi(1,2) - 1;


%Setup electric field and ExB drift for plotting
electric_field_data = load('electric_field.dat');
electric_drift_data = load('electric_drift.dat');
electric_field_data = electric_field_data(1:n_vertex,:);
electric_drift_data = electric_drift_data(1:n_vertex,:);

R = electric_field_data(:,1);
phi = electric_field_data(:,2);
Z = electric_field_data(:,3);
E = electric_field_data(:,4:6);
v_E = electric_drift_data(:,4:6);
v2_E_mod = electric_drift_data(:,7);

E_mod = sqrt(E(:,1).^2 + (R.*E(:,2)).^2 + E(:,3).^2);

switch(grid_kind)
    case(2)
        if (eps_Phi==0)
            scale_field = 1;
            scale_drift = 1;
        else
            scale_field = abs(2*1e-5/eps_Phi)*0.5;
            scale_drift = abs(2*1e-5/eps_Phi)*1e-6*1e-1;
        end
        electric_field = nan*ones(ceil(n_s/skip_s)*ceil((n_theta/skip_theta)),2);
        electric_drift = nan*ones(ceil(n_s/skip_s)*ceil((n_theta/skip_theta)),3);
        counter = 0;
        for i_s = [1:skip_s:n_s]
            for i_theta = [1:skip_theta:n_theta]
                t = (i_s-1)*n_theta + i_theta;
                counter = counter + 1;
                electric_field((counter-1)*3 + 1,:) = electric_field_data(t,[1,3]);
                electric_field((counter-1)*3 + 2,:) = electric_field_data(t,[1,3]) + scale_field*electric_field_data(t,[4,6]);
                electric_field((counter-1)*3 + 3,:) = nan;
                electric_drift((counter-1)*3 + 1,:) = [R(t)*cos(phi(t)),R(t)*sin(phi(t)),Z(t)];
                electric_drift((counter-1)*3 + 2,:) = [R(t)*cos(phi(t)),R(t)*sin(phi(t)),Z(t)] + scale_drift*cylinder2cart(v_E(t,:),R(t),phi(t));
                electric_drift((counter-1)*3 + 3,:) = nan;
            end
        end
    case(4)
        if (eps_Phi==0)
            scale_field = 1;
            scale_drift = 1;
        else
            scale_field = abs(1e-4/eps_Phi)*1e-1*0.3;
            scale_drift = 1e-6*abs(1e-4/eps_Phi)*1e-1*0.3;
        end
        electric_field = nan*ones(ceil(n_vertex/skip)*3,2);
        electric_drift = nan*ones(ceil(n_vertex/skip)*3,3);
        counter = 0;
        for t = [1:skip:n_vertex]
            counter = counter + 1;
            electric_field((counter-1)*3 + 1,:) = electric_field_data(t,[1,3]);
            electric_field((counter-1)*3 + 2,:) = electric_field_data(t,[1,3]) + scale_field*electric_field_data(t,[4,6]);
            electric_field((counter-1)*3 + 3,:) = nan;
            electric_drift((counter-1)*3 + 1,:) = [R(t)*cos(phi(t)),R(t)*sin(phi(t)),Z(t)];
            electric_drift((counter-1)*3 + 2,:) = [R(t)*cos(phi(t)),R(t)*sin(phi(t)),Z(t)] + scale_drift*cylinder2cart(v_E(t,:),R(t),phi(t));
            electric_drift((counter-1)*3 + 3,:) = nan;
        end
end

%Setup modulo v_E and (v_E)_y for plotting
v_E_mod = sqrt(v2_E_mod);
modulo_grid = nan*ones(n_faces*5,3);
for t = [1:n_faces]
    for k = 1:4
        modulo_grid((t-1)*5 + k,[1,2]) = vertex(faces(t,mod(k-1,3) + 1),:);
        modulo_grid((t-1)*5 + k,3) = v_E_mod(faces(t,mod(k-1,3) + 1));
    end
    modulo_grid((t-1)*5 + 5,:) = nan;
end

v_E_cart = cylinder2cart(v_E,R,phi);
v_E_tor_grid = nan*ones(n_faces*5,3);
for t = [1:n_faces]
    for k = 1:4
        v_E_tor_grid((t-1)*5 + k,[1,2]) = vertex(faces(t,mod(k-1,3) + 1),:);
        if(v_E_mod(faces(t,mod(k-1,3) + 1))==0)
            v_E_tor_grid((t-1)*5 + k,3) = 0;
        else
            v_E_tor_grid((t-1)*5 + k,3) = v_E_cart(faces(t,mod(k-1,3) + 1),2)/v_E_mod(faces(t,mod(k-1,3) + 1));
        end
    end
    v_E_tor_grid((t-1)*5 + 5,:) = nan;
end

%Calculate Mach number
clight = 3*1e10;
mass_tungsten = 184 * 1.675 * 1e-24;
eV2erg = 1.60218*1e-12;

v_E_mod_average = sum(v_E_mod,1)/size(v_E_mod,1);
v_thermal = sqrt(2*energy_eV_start*eV2erg/mass_tungsten);
mach_number_average = v_E_mod_average/v_thermal;
mach_number_max = max(v_E_mod)/v_thermal;


%%Create Plots
%-----------------------------------------------------------------------------------------------------------------------
close all
%Draw Poincare plot
xlabel_txt = '$R$ [cm]';
ylabel_txt = '$Z$ [cm]';
figure('Renderer', 'painters', 'Position', [24 730 1214 580])

grid_plot = plot(grid(:,1),grid(:,2),'Color',grid_color);
grid_plot.LineWidth = grid_thickness;
grid_plot.DisplayName = grid_name;
hold on
poincare_plot_electric = plot(poincare_rphiz_electric(:,1),poincare_rphiz_electric(:,3),'.','Color',color_electric);
poincare_plot_electric.DisplayName = 'including strong $\bf{E}$ effects';
poincare_plot = plot(poincare_rphiz(:,1),poincare_rphiz(:,3),'.','Color',color);
poincare_plot.DisplayName = 'default mode';
start_pos_plot = plot(start_pos(:,1),start_pos(:,3),'o','MarkerFaceColor',color_start,'MarkerEdgeColor',color_start);
start_pos_plot.DisplayName = 'Guiding-center starting postions';
hold off
axis equal
lh_poincare = legend([grid_plot,poincare_plot_electric,poincare_plot,start_pos_plot],'Interpreter','latex');
lh_poincare.FontSize = FontSize;
lh_poincare.Location = 'best';
tit = title(['Comparison: Poincare plot ($\varphi=0$)'],'Interpreter','latex');
xlab = xlabel(xlabel_txt,'Interpreter','latex');
ylab = ylabel(ylabel_txt,'Interpreter','latex');
tit.FontSize = FontSize;
xlab.FontSize = FontSize;
ylab.FontSize = FontSize;

saveas(gcf,[path_data_plots,'/poincare_',name_test_case,'.png']);

%Energy / toroidal momentum evolution over time
figure('Position', [24 37 1214 580])
xlabel_txt = '$N_\mathrm{Mappings}$';

ylabel_txt = '$E/E(t = 0) - 1$';
subplot(1,2,1)
e_tot_plot_electric = plot(e_tot_electric(:,1),e_tot_electric(:,2),'Color',color_electric);
e_tot_plot_electric.DisplayName = 'including strong $\bf{E}$ effects';
hold on
e_tot_plot = plot(e_tot(:,1),e_tot(:,2),'Color',color);
e_tot_plot.DisplayName = 'default mode';
hold off
lh_energy = legend([e_tot_plot_electric,e_tot_plot],'Interpreter','latex');
lh_energy.FontSize = FontSize;
lh_energy.Location = 'best';
tit = title(['Comparison: Total energy $E$ evolution'],'Interpreter','latex');
xlab = xlabel(xlabel_txt,'Interpreter','latex');
ylab = ylabel(ylabel_txt,'Interpreter','latex');
tit.FontSize = FontSize;
xlab.FontSize = FontSize;
ylab.FontSize = FontSize;

ylabel_txt = '$p_{\varphi}/p_{\varphi}(t = 0) - 1$';
subplot(1,2,2)
p_phi_plot_electric = plot(p_phi_electric(:,1),p_phi_electric(:,2),'Color',color_electric);
p_phi_plot_electric.DisplayName = 'including strong $\bf{E}$ effects';
hold on
p_phi_plot = plot(p_phi(:,1),p_phi(:,2),'Color',color);
p_phi_plot.DisplayName = 'default mode';
hold off
lh_momentum = legend([p_phi_plot_electric,p_phi_plot],'Interpreter','latex');
lh_momentum.FontSize = FontSize;
lh_momentum.Location = 'best';
tit = title(['Comparison: Toroidal momentum $p_{\varphi}$ evolution'],'Interpreter','latex');
xlab = xlabel(xlabel_txt,'Interpreter','latex');
ylab = ylabel(ylabel_txt,'Interpreter','latex');
tit.FontSize = FontSize;
xlab.FontSize = FontSize;
ylab.FontSize = FontSize;

saveas(gcf,[path_data_plots,'/conservation_',name_test_case,'.png']);

%2D Projection of electric field and drift
xlabel_txt = '$R$ [cm]';
ylabel_txt = '$Z$ [cm]';
figure('Renderer', 'painters','Position', [1300 730 1214 580])

projection_grid_plot = plot(grid(:,1),grid(:,2),'Color',grid_color);
projection_grid_plot.LineWidth = grid_thickness;
projection_grid_plot.DisplayName = grid_name;
hold on
projection_drift_plot = plot(electric_drift(:,1),electric_drift(:,3),'-','Color',electric_drift_color);
projection_drift_plot.LineWidth = drift_thickness;
projection_drift_plot.DisplayName = ['$\bf{E}\times \bf{B}$ drift (scaled 1:1e',num2str(round(log10(1/scale_drift))),')'];
projection_field_plot = plot(electric_field(:,1),electric_field(:,2),'-','Color',electric_field_color);
projection_field_plot.LineWidth = electric_thickness;
projection_field_plot.DisplayName = ['$\bf{E}$-field (scaled 1:1e',num2str(round(log10(1/scale_field))),')'];
hold off

axis equal
lh_projection = legend([projection_grid_plot,projection_field_plot,projection_drift_plot],'Interpreter','latex');
lh_projection.FontSize = FontSize;
lh_projection.Location = 'best';
tit = title(['Toroidal projection ($\varphi=0$) of electric field and drift'],'Interpreter','latex');
xlab = xlabel(xlabel_txt,'Interpreter','latex');
ylab = ylabel(ylabel_txt,'Interpreter','latex');
tit.FontSize = FontSize;
xlab.FontSize = FontSize;
ylab.FontSize = FontSize;

saveas(gcf,[path_data_plots,'/projection_',name_test_case,'.png']);

% Modulo of ExB drift
xlabel_txt = '$R$ [cm]';
ylabel_txt = '$Z$ [cm]';
figure('Position', [1300 37 1214 580])

zlabel_txt = '$v_E$ [cm/s]';
tile_set = tiledlayout(1,2);
tile1 = nexttile;
cla
electric_drift_plot = patch(modulo_grid(:,1),modulo_grid(:,2),modulo_grid(:,3),modulo_grid(:,3),'EdgeColor','interp','FaceColor','none');
electric_drift_plot.LineWidth = drift_thickness;
colormap(tile1,parula)
colorbar
h = get(gca,'DataAspectRatio');
if h(3)==1
      set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))])
else
      set(gca,'DataAspectRatio',[1 1 h(3)])
end
tit = title(['Modulo of $\bf{E}\times \bf{B}$ drift'],'Interpreter','latex');
xlab = xlabel(xlabel_txt,'Interpreter','latex');
ylab = ylabel(ylabel_txt,'Interpreter','latex');
zlab = zlabel(zlabel_txt,'Interpreter','latex');
tit.FontSize = FontSize;
xlab.FontSize = FontSize;
ylab.FontSize = FontSize;
zlab.FontSize = FontSize;

zlabel_txt = '$(v_E)_\mathrm{tor}/v_E$ [cm/s]';
tile2 = nexttile;
cla
electric_drift_plot = patch(v_E_tor_grid(:,1),v_E_tor_grid(:,2),v_E_tor_grid(:,3),v_E_tor_grid(:,3),'EdgeColor','interp','FaceColor','none');
electric_drift_plot.LineWidth = drift_thickness;
%colormap(sub2,cool)
colormap(tile2,cool)
colorbar
h = get(gca,'DataAspectRatio');
if h(3)==1
      set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))])
else
      set(gca,'DataAspectRatio',[1 1 h(3)])
end
tit = title(['Relative toroidale component of $\bf{E}\times \bf{B}$ drift'],'Interpreter','latex');
xlab = xlabel(xlabel_txt,'Interpreter','latex');
ylab = ylabel(ylabel_txt,'Interpreter','latex');
zlab = zlabel(zlabel_txt,'Interpreter','latex');
tit.FontSize = FontSize;
xlab.FontSize = FontSize;
ylab.FontSize = FontSize;
zlab.FontSize = FontSize;

saveas(gcf,[path_data_plots,'/modulo_drift_',name_test_case,'.png']);

%Go back to path of Matlab script
cd(path_script);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%AUXILIARY%FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v_cart] = cylinder2cart(v_cylinder,r,phi)
    v_r = v_cylinder(:,1);    
    v_phi = v_cylinder(:,2);
    v_z = v_cylinder(:,3);
    
    v_cart = nan*ones(size(v_cylinder));
    v_cart(:,1) = v_r.*cos(phi) - v_phi./r.*sin(phi);
    v_cart(:,2) = v_r.*sin(phi) + v_phi./r.*cos(phi);
    v_cart(:,3) = v_z;
end

%% REST
%     grid_plot = plot3(grid(:,1),zeros(size(grid(:,1))),grid(:,2),'Color',grid_color);
%     grid_plot.LineWidth = grid_thickness;
%     grid_plot.DisplayName = 'SOLEDGE3X-EIRENE mesh';
%     hold on
% 
%     %Data for orbits
%     full_orbit_rphiz = data.full_orbit_plot_rphiz_electric;
%     start_pos = load('orbit_start_rphizlambda_electric.dat');
% 
%     phi = full_orbit_rphiz(:,2);
%     full_orbit_xyz = full_orbit_rphiz;
%     full_orbit_xyz(:,2) = 0;
%     full_orbit_xyz = cylinder2cart(full_orbit_xyz,1,phi);
%     orbit_plot = plot3(full_orbit_xyz(:,1),full_orbit_xyz(:,2),full_orbit_xyz(:,3),'.','Color',orbit_color);
%     orbit_plot.DisplayName = 'Guiding-center orbit projection';
%     start_pos_plot = plot3(start_pos(:,1),start_pos(:,2),start_pos(:,3),'o','MarkerFaceColor',orbit_color,'MarkerEdgeColor',orbit_color);
%     start_pos_plot.DisplayName = 'Guiding-center starting postions';
% Plotting 3D drift field of first plane

% xlabel_txt = '$R$ [cm]';
% ylabel_txt = '$Y$ [cm]';
% zlabel_txt = '$Z$ [cm]';
% %figure('Renderer', 'painters', 'Position', [1300 37 1214 580])
% figure('Renderer', 'painters','units','normalized','outerposition',[0 0 1 1])
% 
% %electric_drift_plot = plot3(electric_drift(:,1),electric_drift(:,2),electric_drift(:,3),'-','Color',electric_drift_color);
% drift_position = electric_drift(1:3:end,:);
% drift_direction = electric_drift(2:3:end,:) - drift_position;
% electric_drift_plot = quiver3(drift_position(:,1),drift_position(:,2),drift_position(:,3),drift_direction(:,1),drift_direction(:,2),drift_direction(:,3),'Color',electric_drift_color);
% electric_drift_plot.LineWidth = drift_thickness;
% electric_drift_plot.DisplayName = ['$\bf{E}\times \bf{B}$ drift (arbitrary scale)'];
% hold on
% % grid_drift_plot = plot3(grid(:,1),zeros(size(grid(:,1))),grid(:,2),'Color',grid_color*1.2);
% % grid_drift_plot.LineWidth = grid_thickness;
% % grid_drift_plot.DisplayName = grid_name;
% hold off
% 
% axis equal
% %lh_3D = legend([grid_drift_plot,electric_drift_plot],'Interpreter','latex');
% lh_3D = legend([electric_drift_plot],'Interpreter','latex');
% lh_3D.FontSize = FontSize;
% lh_3D.Location = 'best';
% tit = title(['3D visualization of $\bf{E}\times \bf{B}$ drift'],'Interpreter','latex');
% xlab = xlabel(xlabel_txt,'Interpreter','latex');
% ylab = ylabel(ylabel_txt,'Interpreter','latex');
% zlab = zlabel(zlabel_txt,'Interpreter','latex');
% tit.FontSize = FontSize;
% xlab.FontSize = FontSize;
% ylab.FontSize = FontSize;
% 
% saveas(gcf,[path_data_plots,'/3D_drift_',name_test_case,'.png']);