function data=load_copy(data,path_RUN,path_data_plots,gorilla_plot)
if gorilla_plot.GORILLA_PLOT_NML.boole_J_par == true
    name=strsplit(gorilla_plot.GORILLA_PLOT_NML.filename_J_par,'.');
    data.(name{1})=load(gorilla_plot.GORILLA_PLOT_NML.filename_J_par);
    copyfile([path_RUN,'/',gorilla_plot.GORILLA_PLOT_NML.filename_J_par],path_data_plots);
end

if gorilla_plot.GORILLA_PLOT_NML.boole_e_tot == true
    name=strsplit(gorilla_plot.GORILLA_PLOT_NML.filename_e_tot,'.');
    data.(name{1})=load(gorilla_plot.GORILLA_PLOT_NML.filename_e_tot);
    copyfile([path_RUN,'/',gorilla_plot.GORILLA_PLOT_NML.filename_e_tot],path_data_plots);
end

if gorilla_plot.GORILLA_PLOT_NML.boole_p_phi == true
    name=strsplit(gorilla_plot.GORILLA_PLOT_NML.filename_p_phi,'.');
    data.(name{1})=load(gorilla_plot.GORILLA_PLOT_NML.filename_p_phi);
    copyfile([path_RUN,'/',gorilla_plot.GORILLA_PLOT_NML.filename_p_phi],path_data_plots);
end

if gorilla_plot.GORILLA_PLOT_NML.boole_full_orbit == true
    name=strsplit(gorilla_plot.GORILLA_PLOT_NML.filename_full_orbit_rphiz,'.');
    data.(name{1})=load(gorilla_plot.GORILLA_PLOT_NML.filename_full_orbit_rphiz);
    copyfile([path_RUN,'/',gorilla_plot.GORILLA_PLOT_NML.filename_full_orbit_rphiz],path_data_plots);
    name=strsplit(gorilla_plot.GORILLA_PLOT_NML.filename_full_orbit_sthetaphi,'.');
    data.(name{1})=load(gorilla_plot.GORILLA_PLOT_NML.filename_full_orbit_sthetaphi);
    copyfile([path_RUN,'/',gorilla_plot.GORILLA_PLOT_NML.filename_full_orbit_sthetaphi],path_data_plots);
end

if gorilla_plot.GORILLA_PLOT_NML.boole_poincare_vpar_0 == true
    name=strsplit(gorilla_plot.GORILLA_PLOT_NML.filename_poincare_vpar_0_rphiz,'.');
    data.(name{1})=load(gorilla_plot.GORILLA_PLOT_NML.filename_poincare_vpar_0_rphiz);
    copyfile([path_RUN,'/',gorilla_plot.GORILLA_PLOT_NML.filename_poincare_vpar_0_rphiz],path_data_plots);
    name=strsplit(gorilla_plot.GORILLA_PLOT_NML.filename_poincare_vpar_0_sthetaphi,'.');
    data.(name{1})=load(gorilla_plot.GORILLA_PLOT_NML.filename_poincare_vpar_0_sthetaphi);
    copyfile([path_RUN,'/',gorilla_plot.GORILLA_PLOT_NML.filename_poincare_vpar_0_sthetaphi],path_data_plots);
end


if gorilla_plot.GORILLA_PLOT_NML.boole_poincare_phi_0 == true
    name=strsplit(gorilla_plot.GORILLA_PLOT_NML.filename_poincare_phi_0_rphiz,'.');
    data.(name{1})=load(gorilla_plot.GORILLA_PLOT_NML.filename_poincare_phi_0_rphiz);
    copyfile([path_RUN,'/',gorilla_plot.GORILLA_PLOT_NML.filename_poincare_phi_0_rphiz],path_data_plots);
    name=strsplit(gorilla_plot.GORILLA_PLOT_NML.filename_poincare_phi_0_sthetaphi,'.');
    data.(name{1})=load(gorilla_plot.GORILLA_PLOT_NML.filename_poincare_phi_0_sthetaphi);
    copyfile([path_RUN,'/',gorilla_plot.GORILLA_PLOT_NML.filename_poincare_phi_0_sthetaphi],path_data_plots);
end

end
