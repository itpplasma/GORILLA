%#######################################################################################################################
% description:
%-----------------------------------------------------------------------------------------------------------------------
% * Compute collisionless guiding-center orbit with GORILLA for a trapped Deuterium particle with/without adaptive scheme.
% * Use a field-aligned grid for a non-axisymmetric VMEC MHD equilibrium.
% * Create a figure with the Poincar� plots (\varphi = 0) in cylindrical and symmetry flux coordinates.
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
% authors: Georg Graßler
% created: 16.09.2022

%For Construction of Demonstration
%-----------------------------------------------------------------------------------------------------------------------
%control of script
close all
boole_recalculate_nonadaptive = false;
boole_recalculate_adaptive = false;
boole_recalculate_control_order4 = false;
boole_recalculate_full_orbit = false;
boole_recalculate_full_orbit_adaptive = false;
boole_load_data = false;
poly_order = 2;
poly_order_adaptive = 2;
total_orbit_time = 1d1;
desired_delta_energy = 1.0d-10;
max_n_intermediate_steps = 10000;

boole_full_orbit = true;
boole_e_tot = true;
boole_poincare_vpar_0 = true;
boole_J_par = true;

boole_show_e_tot = true;
boole_show_e_tot_fluctuation = false;
boole_show_e_tot_fluctuation_full_orbit = true;
boole_show_full_orbit = false;
boole_show_J_par = false;
boole_show_J_par_combined = true;
boole_show_poincare_vpar_0_rphiz = false;
boole_show_poincare_vpar_0_rphiz_combined = true;
boole_show_poincare_vpar_0_sthetaphi_combined = true;

n_skip_poincare = 5;
n_skip_poincare_rphiz = 5;
n_skip_J_par = 20;
n_skip_e_tot_fluctuation = 10;
n_skip_e_tot_fluctuation_full_orbit = 1;
n_skip_e_tot = 1;

% order2Color = [102,194,165]/256;
% adaptiveColor = [141,160,203]/256;
% order4Color = [252,141,98]/256;
order2Color = [27,158,119]/256;
% adaptiveColor = [117,112,179]/256;
% order4Color = [217,95,2]/256;
adaptiveColor = [217,95,2]/256;
order4Color = [117,112,179]/256;

MarkerPoincare = '.';
MarkerSizePoincare = 15;
CentralLineWidth = 2;
MarkerEnergyFluctuation = '.';
MarkerSizeEnergyFluctuation = 15;
MarkerSize = 5;
FontSize = 14;

postprocessing_functions_path = 'postprocessing_functions';
addpath(postprocessing_functions_path);

%tests on lamor radius
start_pitch_parameter = 0.3;
n_skip_full_orbit = 1;
energy_eV_start = 3d3;
scaling_factor = 1/10;
ispecies = 2;
boole_guess = true;
m_ratio = 1/(4*1800);
%changed_start_pitch_parameter = 1 - m_ratio * (1-start_pitch_parameter);
% changed_start_pitch_parameter = sqrt(1-(1-start_pitch_parameter^2)*scaling_factor^(1/3));
% start_pitch_parameter = changed_start_pitch_parameter;
zeta = 1 - start_pitch_parameter^2;
changed_energy = energy_eV_start * m_ratio;

if (boole_recalculate_nonadaptive || boole_recalculate_adaptive || boole_recalculate_control_order4 || boole_load_data || boole_recalculate_full_orbit || boole_recalculate_full_orbit_adaptive)
    %filenames for output
    filename_poincare_vpar_0_rphiz = 'poincare_plot_v_par_0_rphiz.dat';
    filename_poincare_vpar_0_rphiz_adaptive = 'poincare_plot_v_par_0_rphiz_adaptive.dat';
    filename_poincare_vpar_0_rphiz_control_order4 = 'poincare_plot_v_par_0_rphiz_control_order4.dat';
    filename_poincare_vpar_0_sthetaphi = 'poincare_plot_v_par_0_sthetaphi.dat';
    filename_poincare_vpar_0_sthetaphi_adaptive = 'poincare_plot_v_par_0_sthetaphi_adaptive.dat';
    filename_poincare_vpar_0_sthetaphi_control_order4 = 'poincare_plot_v_par_0_sthetaphi_control_order4.dat';
    filename_full_orbit_rphiz = 'full_orbit_plot_rphiz.dat';
    filename_full_orbit_rphiz_adaptive = 'full_orbit_plot_rphiz_adaptive.dat';
    filename_full_orbit_rphiz_control_order4 = 'full_orbit_plot_rphiz_control_order4.dat';
    filename_e_tot = 'e_tot.dat';
    filename_e_tot_adaptive = 'e_tot_adaptive.dat';
    filename_e_tot_control_order4 = 'e_tot_control_order4.dat';
    filename_e_tot_full_orbit = 'e_tot_full_orbit.dat';
    filename_e_tot_full_orbit_adaptive = 'e_tot_full_orbit_adaptive.dat';
    filename_J_par = 'J_par.dat';
    filename_J_par_adaptive = 'J_par_adaptive.dat';
    filename_J_par_control_order4 = 'J_par_control_order4.dat';

    %Check of request
    % if (boole_poincare_vpar_0 && boole_full_orbit)
    %     disp('Error: Requested Poincare Plot but have full orbit calculation on instead!');
    %     return
    % end
    if (boole_show_full_orbit && ~boole_full_orbit)
        disp('Error: Requested Full Orbit Plot but have full orbit calculation turned off!');
        return
    end

    %Initialize used paths and files
    %-----------------------------------------------------------------------------------------------------------------------

    %Name of the current calculation to create folder
    name_test_case='example_7';

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
            gorilla.GORILLANML.eps_Phi=0;

        %Coordinate system
            %1 ... (R,phi,Z) cylindrical coordinate system
            %2 ... (s,theta,phi) symmetry flux coordinate system
            gorilla.GORILLANML.coord_system = 2;

        %particle species
            %1 ... electron, 2 ... deuterium ion, 3 ... alpha particle
            gorilla.GORILLANML.ispecies = ispecies;

        %Switch for initial periodic coordinate particle re-location (modulo operation)
            % true ... Particles are re-located at initialization in the case of a periodic coordinate, if they are outside the computation domain.
            % false ... Particles are not re-located at initialization (This might lead to error if particles are outside the computation domain)
            gorilla.GORILLANML.boole_periodic_relocation = false;

        %1 ... numerical RK pusher, 2 ... polynomial pusher
            gorilla.GORILLANML.ipusher = 2;

        %true ... Face guessing algorithm is used, false ... NOT used
            gorilla.GORILLANML.boole_guess = boole_guess;

        %Polynomial order for orbit pusher (from 2 to 4)
            gorilla.GORILLANML.poly_order = poly_order;

        % Switches to calculate optional quantities
            gorilla.GORILLANML.boole_time_Hamiltonian = false;
            gorilla.GORILLANML.boole_gyrophase = false;

        % Switch for adaptive time step scheme
            % false ... no adaptive scheme
            % true ... adaptive scheme to ensure energy conservation up to specified fluctuation
            gorilla.GORILLANML.boole_adaptive_time_steps = false;

        % Allowed relative fluctuation of energy between entry and exit of a tetrahedron
            gorilla.GORILLANML.desired_delta_energy  = desired_delta_energy;

        % Maximum number of intermediate steps allowed for splitting up an orbit section
            gorilla.GORILLANML.max_n_intermediate_steps = max_n_intermediate_steps;


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
            %4 ... SOLEDGE3X_EIRENE grid
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
            gorilla_plot.GORILLA_PLOT_NML.total_orbit_time = total_orbit_time;

        %Total Energy of particle in eV
            gorilla_plot.GORILLA_PLOT_NML.energy_eV_start = energy_eV_start;

        %Switch for plotting Poincar� cuts at toroidal variable $\varphi$ = 0
            gorilla_plot.GORILLA_PLOT_NML.boole_poincare_phi_0 = false;

        %Switch for plotting Poincar� cuts at parallel velocity $v_\parallel$ = 0
            gorilla_plot.GORILLA_PLOT_NML.boole_poincare_vpar_0 = boole_poincare_vpar_0;

            %Number of skipped (non-printed) Poincar� cuts at parallel velocity $v_\parallel$ = 0
                gorilla_plot.GORILLA_PLOT_NML.n_skip_vpar_0 = 1;

            %Filename for Poincar� cuts at parallel velocity $v_\parallel$ = 0 in cylindrical coordinates (R,$\varphi$,Z)
                gorilla_plot.GORILLA_PLOT_NML.filename_poincare_vpar_0_rphiz = filename_poincare_vpar_0_rphiz;

            %Filename for Poincar� cuts at parallel velocity $v_\parallel$ = 0 in symmetry flux coordinates (s,$\vartheta$,$\varphi$)
                gorilla_plot.GORILLA_PLOT_NML.filename_poincare_vpar_0_sthetaphi = filename_poincare_vpar_0_sthetaphi;

        %Switch for plotting full orbit
            gorilla_plot.GORILLA_PLOT_NML.boole_full_orbit = false;

        %Number of skipped (non-printed tetrahedra passings) full orbit
            gorilla_plot.GORILLA_PLOT_NML.n_skip_full_orbit = n_skip_full_orbit;

        %Filename for full orbit in cylindrical coordinates (R,$\varphi$,Z)
            gorilla_plot.GORILLA_PLOT_NML.filename_full_orbit_rphiz = filename_full_orbit_rphiz;

        %Plot invariances of motion (ONLY for single orbits)

            %Switch for plotting total particle energy
                gorilla_plot.GORILLA_PLOT_NML.boole_e_tot = boole_e_tot;

            %Switch for plotting canoncial (toroidal) angular momentum $p_\varphi$
                gorilla_plot.GORILLA_PLOT_NML.boole_p_phi = false;

            %Switch for parallel adiabatic invariant $J_\parallel$
                gorilla_plot.GORILLA_PLOT_NML.boole_J_par = boole_J_par;

            %Filename for parallel adiabatic invariant $J_\parallel$
                gorilla_plot.GORILLA_PLOT_NML.filename_J_par = filename_J_par;

            %Filename for total particle energy
                gorilla_plot.GORILLA_PLOT_NML.filename_e_tot = filename_e_tot;

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
                gorilla_plot.GORILLA_PLOT_NML.start_pitch_parameter = start_pitch_parameter;


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

    %Run without energy control
    if (boole_recalculate_nonadaptive)
        ! ./test_gorilla_main.x
    end

    %Load Data-files and copy them into data_plots folder
    data=struct();
    data=load_copy(data,path_RUN,path_data_plots,gorilla_plot);
    
    %polynominal order = 4
    gorilla.GORILLANML.poly_order = 4;
    gorilla_plot.GORILLA_PLOT_NML.filename_J_par = filename_J_par_control_order4;
    gorilla_plot.GORILLA_PLOT_NML.filename_e_tot = filename_e_tot_control_order4;
    gorilla_plot.GORILLA_PLOT_NML.filename_full_orbit_rphiz = filename_full_orbit_rphiz_control_order4;
    gorilla_plot.GORILLA_PLOT_NML.filename_poincare_vpar_0_rphiz = filename_poincare_vpar_0_rphiz_control_order4;
    gorilla_plot.GORILLA_PLOT_NML.filename_poincare_vpar_0_sthetaphi = filename_poincare_vpar_0_sthetaphi_control_order4;
    gorilla.write([path_RUN,'/gorilla.inp']);
    gorilla_plot.write([path_RUN,'/gorilla_plot.inp']);
    if(boole_recalculate_control_order4)
        ! ./test_gorilla_main.x
    end
    data=load_copy(data,path_RUN,path_data_plots,gorilla_plot);
    
    %Changes for adaptive RUN
    gorilla.GORILLANML.poly_order = poly_order_adaptive;
    gorilla.GORILLANML.boole_adaptive_time_steps = true;
    gorilla_plot.GORILLA_PLOT_NML.filename_poincare_vpar_0_rphiz = filename_poincare_vpar_0_rphiz_adaptive;
    gorilla_plot.GORILLA_PLOT_NML.filename_full_orbit_rphiz = filename_full_orbit_rphiz_adaptive;
    gorilla_plot.GORILLA_PLOT_NML.filename_poincare_vpar_0_sthetaphi = filename_poincare_vpar_0_sthetaphi_adaptive;
    gorilla_plot.GORILLA_PLOT_NML.filename_e_tot = filename_e_tot_adaptive;
    gorilla_plot.GORILLA_PLOT_NML.filename_J_par = filename_J_par_adaptive;
    gorilla.write([path_RUN,'/gorilla.inp']);
    gorilla_plot.write([path_RUN,'/gorilla_plot.inp']);

    %Run with energy control
    if (boole_recalculate_adaptive)
        ! ./test_gorilla_main.x
    end

    data=load_copy(data,path_RUN,path_data_plots,gorilla_plot);

    %observe energy fluctuation in detail
    gorilla.GORILLANML.poly_order = poly_order;
    gorilla.GORILLANML.boole_adaptive_time_steps = false;
    gorilla_plot.GORILLA_PLOT_NML.boole_poincare_vpar_0 = false;
    gorilla_plot.GORILLA_PLOT_NML.boole_full_orbit = true;
    gorilla_plot.GORILLA_PLOT_NML.boole_J_par = false;
    gorilla_plot.GORILLA_PLOT_NML.total_orbit_time = total_orbit_time/100;
    gorilla_plot.GORILLA_PLOT_NML.filename_e_tot = filename_e_tot_full_orbit;
    gorilla.write([path_RUN,'/gorilla.inp']);
    gorilla_plot.write([path_RUN,'/gorilla_plot.inp']);
    
    %Run with energy control
    if (boole_recalculate_full_orbit)
        ! ./test_gorilla_main.x
    end

    data=load_copy(data,path_RUN,path_data_plots,gorilla_plot);
    
    %observe energy fluctuation in detail
    gorilla.GORILLANML.poly_order = poly_order_adaptive;
    gorilla.GORILLANML.boole_adaptive_time_steps = true;
    gorilla_plot.GORILLA_PLOT_NML.filename_e_tot = filename_e_tot_full_orbit_adaptive;
    gorilla.write([path_RUN,'/gorilla.inp']);
    gorilla_plot.write([path_RUN,'/gorilla_plot.inp']);
    
    %Run with energy control
    if (boole_recalculate_full_orbit_adaptive)
        ! ./test_gorilla_main.x
    end
    
    data=load_copy(data,path_RUN,path_data_plots,gorilla_plot);
    
    cd(path_script)
    
end

%Create plots of generated data
%-----------------------------------------------------------------------------------------------------------------------

%Load results of both runs
if (boole_poincare_vpar_0)
    poincare_rphiz_data = data.poincare_plot_v_par_0_rphiz;
    poincare_rphiz_adaptive_data = data.poincare_plot_v_par_0_rphiz_adaptive;
    poincare_rphiz_control_order4_data = data.poincare_plot_v_par_0_rphiz_control_order4;
    poincare_sthetaphi_data = data.poincare_plot_v_par_0_sthetaphi;
    poincare_sthetaphi_adaptive_data = data.poincare_plot_v_par_0_sthetaphi_adaptive;
    poincare_sthetaphi_control_order4_data = data.poincare_plot_v_par_0_sthetaphi_control_order4;
end
if (boole_show_full_orbit)
    full_orbit_data = data.full_orbit_plot_rphiz;
    full_orbit_adaptive_data = data.full_orbit_plot_rphiz_adaptive;
end
if (boole_show_e_tot_fluctuation_full_orbit && boole_full_orbit)
    e_tot_full_orbit_data = data.e_tot_full_orbit;
    e_tot_full_orbit_adaptive_data = data.e_tot_full_orbit_adaptive;
end
if (boole_e_tot)
    e_tot_data = data.e_tot;
    e_tot_adaptive_data = data.e_tot_adaptive;
end
if (boole_J_par)
    J_par_data = data.J_par;
    J_par_adaptive_data = data.J_par_adaptive;
    J_par_control_order4_data = data.J_par_control_order4;
end

%Energy evolution comparison
if (boole_show_e_tot)
    xlabel_txt = '$N_{Bounces}$';
    ylabel_txt = '$E/E(\tau = 0)$';
    
    figure('Renderer', 'painters', 'Position', [24 37 1214 580])
    subplot(1,2,1)
    plot(e_tot_data(1:n_skip_e_tot:end,1),e_tot_data(1:n_skip_e_tot:end,2)/e_tot_data(1,2)-1,'Color',order2Color)
    tit = title(['Energy evolution (order ',num2str(poly_order),')'],'Interpreter','latex');
    xlab = xlabel(xlabel_txt,'Interpreter','latex');
    ylab = ylabel(ylabel_txt,'Interpreter','latex');
    tit.FontSize = FontSize;
    xlab.FontSize = FontSize;
    ylab.FontSize = FontSize;
    
    subplot(1,2,2)
    plot(e_tot_adaptive_data(1:n_skip_e_tot:end,1),e_tot_adaptive_data(1:n_skip_e_tot:end,2)/e_tot_adaptive_data(1,2)-1,'Color',adaptiveColor)
    xlab = xlabel(xlabel_txt,'Interpreter','latex');
    ylab = ylabel(ylabel_txt,'Interpreter','latex');
    tit = title(['Energy evolution (adaptive/order ',num2str(poly_order_adaptive),')'],'Interpreter','latex');
    xlab.FontSize = FontSize;
    ylab.FontSize = FontSize;
    tit.FontSize = FontSize;
end

%Energy fluctuation comparison
if (boole_show_energy_fluctuation && ~boole_full_orbit)
    xlabel_txt = 'N_{Bounces}';
    ylabel_txt = '$\delta$';
    energy_fluctuation_data = nan*ones([length(e_tot_data)-1,2]);
    %energy_fluctuation_data(:,2) = (e_tot_data(2:end,2) - e_tot_data(1:end-1,2))./(abs(e_tot_data(1:end-1,2))+abs(e_tot_data(1:end-1,2)));
    energy_fluctuation_data(:,2) = e_tot_data(2:end,2)./e_tot_data(1:end-1,2) -1;
    energy_fluctuation_adaptive_data = nan*ones([length(e_tot_adaptive_data)-1,2]);
    %energy_fluctuation_adaptive_data(:,2) = (e_tot_adaptive_data(2:end,2) - e_tot_adaptive_data(1:end-1,2))./(abs(e_tot_adaptive_data(1:end-1,2))+abs(e_tot_adaptive_data(1:end-1,2)));
    energy_fluctuation_adaptive_data(:,2) = e_tot_adaptive_data(2:end,2)./e_tot_adaptive_data(1:end-1,2) -1;
    
    figure('Renderer', 'painters', 'Position', [1288 37 1214 580])
    subplot(1,2,1)
    energy_fluctuation_data(:,1) = e_tot_data(2:end,1);
    plot(energy_fluctuation_data(:,1),energy_fluctuation_data(:,2),'.','MarkerSize',10)
    tit = title(['Energy fluctuation (order ',num2str(poly_order),')'],'Interpreter','latex');
    xlab = xlabel(xlabel_txt,'Interpreter','latex');
    ylab = ylabel(ylabel_txt,'Interpreter','latex');
    tit.FontSize = FontSize;
    xlab.FontSize = FontSize;
    ylab.FontSize = FontSize;
    
    subplot(1,2,2)
    energy_fluctuation_adaptive_data(:,1) = e_tot_adaptive_data(2:end,1);
    plot(energy_fluctuation_adaptive_data(:,1),energy_fluctuation_adaptive_data(:,2),'.','MarkerSize',10)
    tit = title(['Energy fluctuation (adaptive/order ',num2str(poly_order_adaptive),')'],'Interpreter','latex');
    xlab = xlabel(xlabel_txt,'Interpreter','latex');
    ylab = ylabel(ylabel_txt,'Interpreter','latex');
    tit.FontSize = FontSize;
    xlab.FontSize = FontSize;
    ylab.FontSize = FontSize;
end

%Energy fluctuation comparison full orbit
if (boole_show_e_tot_fluctuation_full_orbit && boole_full_orbit)
    xlabel_txt = '$\Delta{}\tau$';
    ylabel_txt = '$\delta$';
    energy_fluctuation_full_orbit_data = nan*ones([length(e_tot_full_orbit_data)-1,2]);
    energy_fluctuation_full_orbit_data(:,2) = abs(e_tot_full_orbit_data(2:end,2)./e_tot_full_orbit_data(1:end-1,2) -1);
    energy_fluctuation_full_orbit_adaptive_data = nan*ones([length(e_tot_full_orbit_adaptive_data)-1,2]);
    energy_fluctuation_full_orbit_adaptive_data(:,2) = abs(e_tot_full_orbit_adaptive_data(2:end,2)./e_tot_full_orbit_adaptive_data(1:end-1,2) -1);
    
    figure('Renderer', 'painters', 'Position', [1288 37 1214 580])
    subplot(1,2,1)
    energy_fluctuation_full_orbit_data(:,1) = e_tot_full_orbit_data(2:end,1)-e_tot_full_orbit_data(1:end-1,1);
    plot(energy_fluctuation_full_orbit_data(1:n_skip_e_tot_fluctuation_full_orbit:end,1),energy_fluctuation_full_orbit_data(1:n_skip_e_tot_fluctuation_full_orbit:end,2),MarkerEnergyFluctuation,'MarkerSize',MarkerSizeEnergyFluctuation,'MarkerEdgeColor',order2Color,'MarkerFaceColor',order2Color)
    tit = title(['Energy fluctuation (order ',num2str(poly_order),')'],'Interpreter','latex');
    xlab = xlabel(xlabel_txt,'Interpreter','latex');
    ylab = ylabel(ylabel_txt,'Interpreter','latex');
    tit.FontSize = FontSize;
    xlab.FontSize = FontSize;
    ylab.FontSize = FontSize;
    
    set(gcf, 'units', 'normalized');
    box_x1 = 0;
    box_y1 = 0;
    box_x2 = 4d-5;
    box_x3 = box_x2;
    box_x4 = box_x1;
    box_y2 = box_y1;
    box_y3 = desired_delta_energy;
    box_y4 = box_y3;
    
    ax1 = gca;
    [xfig(1), yfig(1)] = axescoord2figurecoord((box_x3-box_x4)*1/3+box_x4, box_y3);
    ax2 = axes('position',[.16 .72 .10 .15]); % [left bottom width height];
    box on % put box around new pair of axes
    plot(energy_fluctuation_full_orbit_data(1:n_skip_e_tot_fluctuation_full_orbit:end,1),energy_fluctuation_full_orbit_data(1:n_skip_e_tot_fluctuation_full_orbit:end,2),MarkerEnergyFluctuation,'MarkerSize',MarkerSizeEnergyFluctuation,'MarkerEdgeColor',order2Color,'MarkerFaceColor',order2Color)
    xlim([box_x1,box_x3]);
    ylim([box_y1,box_y3]);
    [xfig(2), yfig(2)] = axescoord2figurecoord((box_x2-box_x1)*1/3+box_x1, box_y1);
    annotation('line', xfig, yfig);
    
    subplot(1,2,2)
    energy_fluctuation_full_orbit_adaptive_data(:,1) = e_tot_full_orbit_adaptive_data(2:end,1)-e_tot_full_orbit_adaptive_data(1:end-1,1);
    plot(energy_fluctuation_full_orbit_adaptive_data(1:n_skip_e_tot_fluctuation_full_orbit:end,1),energy_fluctuation_full_orbit_adaptive_data(1:n_skip_e_tot_fluctuation_full_orbit:end,2),MarkerEnergyFluctuation,'MarkerSize',MarkerSizeEnergyFluctuation,'MarkerEdgeColor',adaptiveColor,'MarkerFaceColor',adaptiveColor);
    tit = title(['Energy fluctuation (adaptive/order ',num2str(poly_order_adaptive),')'],'Interpreter','latex');
    xlab = xlabel(xlabel_txt,'Interpreter','latex');
    ylab = ylabel(ylabel_txt,'Interpreter','latex');
    tit.FontSize = FontSize;
    xlab.FontSize = FontSize;
    ylab.FontSize = FontSize;
    
end

%Poincare rphiz comparison
if (boole_show_poincare_vpar_0_rphiz)
    figure('Renderer', 'painters', 'Position', [100 100 1300 1000])
    subplot(1,2,1)
    plot(poincare_rphiz_data(:,1),poincare_rphiz_data(:,3),'.','MarkerSize',10)
    title(['Poincare cut (order ',num2str(poly_order),')'])
    subplot(1,2,2)
    plot(poincare_rphiz_adaptive_data(:,1),poincare_rphiz_adaptive_data(:,3),'.','MarkerSize',10)
    title(['Poincare cut (adaptive/order ',num2str(poly_order_adaptive),')'])
end

%Poincare rphiz comparison direct
if (boole_show_poincare_vpar_0_rphiz_combined)
    xlabel_txt = '$R$ [cm]';
    ylabel_txt = '$Z$ [cm]';
    if(boole_show_poincare_vpar_0_sthetaphi_combined)
        figure('Renderer', 'painters', 'Position', [30 734 1214 580])
        subplot(1,2,1)
    else
        figure('Renderer', 'painters', 'Position', [100 100 1300 1000])
    end
    hold on
    h1 = plot(poincare_rphiz_data(1:n_skip_poincare_rphiz:end,2),poincare_rphiz_data(1:n_skip_poincare_rphiz:end,1),MarkerPoincare,'MarkerSize',MarkerSizePoincare,'MarkerEdgeColor',order2Color,'MarkerFaceColor',order2Color);
    h3 = plot(poincare_rphiz_control_order4_data(1:n_skip_poincare_rphiz:end,2),poincare_rphiz_control_order4_data(1:n_skip_poincare_rphiz:end,1),MarkerPoincare,'MarkerSize',MarkerSizePoincare,'MarkerEdgeColor',order4Color,'MarkerFaceColor',order4Color);
    h2 = plot(poincare_rphiz_adaptive_data(1:n_skip_poincare_rphiz:end,2),poincare_rphiz_adaptive_data(1:n_skip_poincare_rphiz:end,1),MarkerPoincare,'MarkerSize',MarkerSizePoincare,'MarkerEdgeColor',adaptiveColor,'MarkerFaceColor',adaptiveColor);
    hold off
    tit = title('Poincare cut at $v_{\parallel} = 0$','Interpreter','latex');
    leg = legend([h1 h2 h3],{['order ',num2str(poly_order)],['order ',num2str(poly_order_adaptive),' (adaptive)'],'order 4'},'Interpreter','latex');
    set(leg,'Location','Southeast')
    leg.FontSize = FontSize;
    %set(leg,'Location','Southeast')
    xlab = xlabel(xlabel_txt,'interpreter','latex');
    ylab = ylabel(ylabel_txt,'interpreter','latex');
    tit.FontSize = FontSize;
    xlab.FontSize = FontSize;
    ylab.FontSize = FontSize;
end

%Poincare sthetaphi comparison direct
if (boole_show_poincare_vpar_0_sthetaphi_combined)
    xlabel_txt = '$\vartheta$';
    ylabel_txt = '$s$';
    if (boole_show_poincare_vpar_0_rphiz_combined)
        subplot(1,2,2)
    else
        figure('Renderer', 'painters', 'Position', [100 100 1300 1000])
    end
    hold on
    h1 = plot(poincare_sthetaphi_data(1:n_skip_poincare:end,2),poincare_sthetaphi_data(1:n_skip_poincare:end,1),MarkerPoincare,'MarkerSize',MarkerSizePoincare,'MarkerEdgeColor',order2Color,'MarkerFaceColor',order2Color);
    h3 = plot(poincare_sthetaphi_control_order4_data(1:n_skip_poincare:end,2),poincare_sthetaphi_control_order4_data(1:n_skip_poincare:end,1),MarkerPoincare,'MarkerSize',MarkerSizePoincare,'MarkerEdgeColor',order4Color,'MarkerFaceColor',order4Color);
    h2 = plot(poincare_sthetaphi_adaptive_data(1:n_skip_poincare:end,2),poincare_sthetaphi_adaptive_data(1:n_skip_poincare:end,1),MarkerPoincare,'MarkerSize',MarkerSizePoincare,'MarkerEdgeColor',adaptiveColor,'MarkerFaceColor',adaptiveColor);
    hold off
    tit = title('Poincare cut at $v_{\parallel} = 0$','Interpreter','latex');
    leg = legend([h1 h2 h3],{['order ',num2str(poly_order)],['order ',num2str(poly_order_adaptive),' (adaptive)'],'order 4'},'Interpreter','latex');
    set(leg,'Location','Southeast')
    leg.FontSize = FontSize;
    %set(leg,'Location','Southeast')
    xlab = xlabel(xlabel_txt,'interpreter','latex');
    ylab = ylabel(ylabel_txt,'interpreter','latex');
    tit.FontSize = FontSize;
    xlab.FontSize = FontSize;
    ylab.FontSize = FontSize;
    xlim([0 2*pi])
end

%J_par comparison
if (boole_show_J_par)
    xlabel_txt = '$N_{mappings}$';
    ylabel_txt = '$J_{\parallel}$';
    figure('Renderer', 'painters', 'Position', [1300 100 1200 1000])
    subplot(1,2,1)
    plot(J_par_data(:,1),J_par_data(:,2)/J_par_data(1,2),'.','MarkerSize',10)
    title(['$J_{\parallel}$ over number of bounces (order ',num2str(poly_order),')'],'interpreter','latex')
    xlabel(xlabel_txt,'interpreter','latex')
    ylabel(ylabel_txt,'interpreter','latex')
    subplot(1,2,2)
    plot(J_par_adaptive_data(:,1),J_par_adaptive_data(:,2)/J_par_adaptive_data(1,2),'.','MarkerSize',10)
    title(['$J_{\parallel}$ over number of bounces (adaptive/order ',num2str(poly_order_adaptive),')'],'interpreter','latex')
    xlabel(xlabel_txt,'interpreter','latex')
    ylabel(ylabel_txt,'interpreter','latex')
end

%J_par comparison direct
if (boole_show_J_par_combined)
    xlabel_txt = '$N_{Bounces}$';
    ylabel_txt = '$J_{\parallel}/J_{\parallel}(\tau = 0)$';
    figure('Renderer', 'painters', 'Position', [1288 734 1214 580])
    hold on
    h1 = plot(J_par_data(1:n_skip_J_par:end,1),J_par_data(1:n_skip_J_par:end,2)/J_par_data(1,2),'o','MarkerSize',MarkerSize,'MarkerEdgeColor',order2Color,'MarkerFaceColor',order2Color);
    h3 = plot(J_par_control_order4_data(1:n_skip_J_par:end,1),J_par_control_order4_data(1:n_skip_J_par:end,2)/J_par_control_order4_data(1,2),'o','MarkerSize',MarkerSize,'MarkerEdgeColor',order4Color,'MarkerFaceColor',order4Color);
    h2 = plot(J_par_adaptive_data(1:n_skip_J_par:end,1),J_par_adaptive_data(1:n_skip_J_par:end,2)/J_par_adaptive_data(1,2),'o','MarkerSize',MarkerSize,'MarkerEdgeColor',adaptiveColor,'MarkerFaceColor',adaptiveColor);
    plot(J_par_data(1:n_skip_J_par:end,1),ones(length(J_par_data(1:n_skip_J_par:end,1)),1),'k-','LineWidth',CentralLineWidth)
    hold off
    tit = title(['$J_{\parallel}$ over number of bounces'],'interpreter','latex');
    leg = legend([h1 h2 h3],{['order ',num2str(poly_order)],['order ',num2str(poly_order_adaptive),' (adaptive)'],'order 4'}, 'Interpreter','latex');
    leg.FontSize = FontSize;
    %set(leg,'Location','Southeast')
    xlab = xlabel(xlabel_txt,'interpreter','latex');
    ylab = ylabel(ylabel_txt,'interpreter','latex');
    tit.FontSize = FontSize;
    xlab.FontSize = FontSize;
    ylab.FontSize = FontSize;
    xlim([0 max(J_par_data(:,1))])
end

%Full orbit comparison (2D projection)
if (boole_show_full_orbit)
    figure('Renderer', 'painters', 'Position', [100 100 1300 1000])
    subplot(1,2,1)
    plot(full_orbit_data(:,1),full_orbit_data(:,3))
    hold on
    plot(full_orbit_data(end-1:end,1),full_orbit_data(end-1:end,3),'r')
    plot(full_orbit_data(end,1),full_orbit_data(end,3),'r.','MarkerSize',10)
    hold off
    title(['Full orbit projection (order ',num2str(poly_order),')'])
    subplot(1,2,2)
    plot(full_orbit_adaptive_data(:,1),full_orbit_adaptive_data(:,3))
    hold on
    plot(full_orbit_adaptive_data(end-1:end,1),full_orbit_adaptive_data(end-1:end,3),'r')
    plot(full_orbit_adaptive_data(end,1),full_orbit_adaptive_data(end,3),'r.','MarkerSize',10)
    hold off
    title(['Full orbit projection (adaptive/order ',num2str(poly_order_adaptive),')'])
end