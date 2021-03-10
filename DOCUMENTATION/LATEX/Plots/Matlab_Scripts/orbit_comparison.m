%%%  Plot orbit comparison
close all

% FontSize
fsize = 18;
fsize_axis = 16;

%
    path_raw_data = '../RAW_Data/orbit_comparison/';
    path_exported_figures = '../Exported_Figures/';
%     
    postprocessing_functions_path = './postprocessing_functions';
    addpath(postprocessing_functions_path);
%
    vertices_rphiz = load([path_raw_data,'vertices_rphiz_20_20_20.dat']);
    ntetra = numel(vertices_rphiz(:,1))/3;
    vec_R = nan(ntetra*5,1);
    vec_Z = nan(ntetra*5,1);
 %
    color_green = myColorFun('green');
 %
     for i = 1:ntetra
           vec_R(i*5-4,1) = vertices_rphiz(i*3-2,1);
           vec_R(i*5-3,1) = vertices_rphiz(i*3-1,1);
           vec_R(i*5-2,1) = vertices_rphiz(i*3,1);
           vec_R(i*5-1,1) = vertices_rphiz(i*3-2,1);
           
           vec_Z(i*5-4,1) = vertices_rphiz(i*3-2,3);
           vec_Z(i*5-3,1) = vertices_rphiz(i*3-1,3);
           vec_Z(i*5-2,1) = vertices_rphiz(i*3,3);
           vec_Z(i*5-1,1) = vertices_rphiz(i*3-2,3);
     end
 %
    orbit_cylindrical = load([path_raw_data,'orbit_comparison_cylindrical.dat']);
    orbit_symflux = load([path_raw_data,'orbit_comparison_symflux.dat']);
    orbit_ode45 = load([path_raw_data,'orbit_comparison_ode45_3.dat']);
 %
     box_a_xlim = [144.5,149];
     box_a_ylim = [19,31];
     box_b_xlim = [192,196.5];
     box_b_ylim = [-10,31];
 
 start_point = [196.06889699753378,8.9223047279345824];
 
 num_points_cyl = 4;
 num_points_symflux = 4;
 num_points_ode = 4;
 
 fh = figure;
 fh.Units = 'normalized';
 fh.Position = [0.3542 0.1380 0.3370 0.7657];
 
 %subplot(2,2,[1,3])
 plot(vec_R,vec_Z,'Color',[0.7,0.7,0.7]);
 hold on
%  plot(orbit_cylindrical([1:num_points_cyl:end],1),orbit_cylindrical([1:num_points_cyl:end],3),'r.');
 plot(orbit_symflux([1:num_points_symflux:end],1),orbit_symflux([1:num_points_symflux:end],3),'.','Color','r');
%  hp1 = plot(230,0,'r.','MarkerSize',5,'MarkerFaceColor','r');
%  hp2 = plot(230,0,'.','Color',color_green,'MarkerSize',5,'MarkerFaceColor',color_green);
%  hp3 = plot(230,0,'b.','MarkerSize',5,'MarkerFaceColor','b');
%  hp4 = plot(orbit_ode45(:,1),orbit_ode45(:,3),'b-','LineWidth',0.7);
%  plot(box_a_xlim([1,2,2,1,1]),box_a_ylim([1,1,2,2,1]),'k-','LineWidth',1);
%  plot(box_b_xlim([1,2,2,1,1]),box_b_ylim([1,1,2,2,1]),'k-','LineWidth',1);
%  text(box_a_xlim(1),box_a_ylim(1) + (box_a_ylim(2)-box_a_ylim(1))/2,'b','HorizontalAlignment','right','FontSize',18,'Color','k')
%  text(box_b_xlim(2),box_b_ylim(1) + (box_b_ylim(2)-box_b_ylim(1))/2,'c','HorizontalAlignment','left','FontSize',18,'Color','k')
 
 %Hidden Plots outside the visible range
%  mh2 = plot(300,5.36,'d','MarkerSize',8,'MarkerFaceColor',color_green,'Color',color_green);
%  mh3 = plot(300,14.41,'^','MarkerSize',8,'MarkerFaceColor','b','Color','b');
%  mh1 = plot(300,12.14,'s','MarkerSize',8,'MarkerFaceColor','r','Color','r');
 
%  xlabel('R [cm]','Interpreter','latex','FontSize',fsize)
%  ylabel('Z [cm]','Interpreter','latex','FontSize',fsize)
%  leg = legend([mh1,mh2,mh3],{'GORILLA (R,$\varphi$,Z)','GORILLA (s,$\vartheta$,$\varphi$)','Exact orbit'},'Interpreter','latex','Location','southeast')
 xlim([105,220]);
 ylim([-95,90]);
 %text(210,-85,'a','HorizontalAlignment','right','FontSize',14,'Color','k')
%  title('$\mathbf{a}$','Interpreter','latex','FontSize',fsize)
%  set(leg,'FontSize',fsize+2)
%  set(gca,'FontSize',fsize_axis)
 
 hold off
 axis tight
 axis off
 
 
 nskip_cylindrical = 500;
 nskip_symflux = 500;
 nskip_ode45 = 70;
 nskip_ode45_2 = 15;
 
 
 
 L1 = orbit_cylindrical(:,1)>=box_a_xlim(1) & orbit_cylindrical(:,1)<=box_a_xlim(2) & orbit_cylindrical(:,3)>=box_a_ylim(1) & orbit_cylindrical(:,3)<=box_a_ylim(2);
 L2 = orbit_symflux(:,1)>=box_a_xlim(1) & orbit_symflux(:,1)<=box_a_xlim(2) & orbit_symflux(:,3)>=box_a_ylim(1) & orbit_symflux(:,3)<=box_a_ylim(2);
 L3 = orbit_ode45(:,1)>=box_a_xlim(1) & orbit_ode45(:,1)<=box_a_xlim(2) & orbit_ode45(:,3)>=box_a_ylim(1) & orbit_ode45(:,3)<=box_a_ylim(2);
 
 %
%  subplot(2,2,[2])
%  plot(vec_R,vec_Z,'Color',[0.7,0.7,0.7]);
%  hold on
%  plot(orbit_cylindrical(L1,1),orbit_cylindrical(L1,3),'r.');
%  
%  plot(orbit_symflux(L2,1),orbit_symflux(L2,3),'Marker','.','Color',color_green,'LineStyle','none');
%  plot(orbit_ode45(:,1),orbit_ode45(:,3),'b','LineWidth',2);
%  plot(orbit_cylindrical([1:nskip_cylindrical:end],1),orbit_cylindrical([1:nskip_cylindrical:end],3),'r','LineStyle','none');
%  plot(orbit_symflux([1:nskip_symflux:end],1),orbit_symflux([1:nskip_symflux:end],3),'Marker','None','Color',color_green,'LineStyle','none');
%  plot(orbit_ode45([1:nskip_ode45:end],1),orbit_ode45([1:nskip_ode45:end],3),'Marker','None','Color','b','LineStyle','none');
%  
%  % Plot additional Markers
%  plot(145.25,21.5,'d','MarkerSize',8,'MarkerFaceColor',color_green,'Color',color_green);
%  plot(148.4,30.15,'^','MarkerSize',8,'MarkerFaceColor','b','Color','b');
%  plot(145.8,21.7,'s','MarkerSize',8,'MarkerFaceColor','r','Color','r');
%  
%  xlim(box_a_xlim);
%  ylim(box_a_ylim);
%  xlabel('R [cm]','Interpreter','latex','FontSize',fsize)
%  ylabel('Z [cm]','Interpreter','latex','FontSize',fsize)
%  hold off
%  title('$\mathbf{b}$','Interpreter','latex','FontSize',fsize)
%  set(gca,'FontSize',fsize_axis)
%  
%  L1 = orbit_cylindrical(:,1)>=box_b_xlim(1) & orbit_cylindrical(:,1)<=box_b_xlim(2) & orbit_cylindrical(:,3)>=box_b_ylim(1) & orbit_cylindrical(:,3)<=box_b_ylim(2);
%  L2 = orbit_symflux(:,1)>=box_b_xlim(1) & orbit_symflux(:,1)<=box_b_xlim(2) & orbit_symflux(:,3)>=box_b_ylim(1) & orbit_symflux(:,3)<=box_b_ylim(2);
%  L3 = orbit_ode45(:,1)>=box_a_xlim(1) & orbit_ode45(:,1)<=box_a_xlim(2) & orbit_ode45(:,3)>=box_a_ylim(1) & orbit_ode45(:,3)<=box_a_ylim(2);
%  
%  %
%  subplot(2,2,[4])
%  plot(vec_R,vec_Z,'Color',[0.7,0.7,0.7]);
%  hold on
%  plot(orbit_cylindrical(L1,1),orbit_cylindrical(L1,3),'r.');
%  
%  plot(orbit_symflux(L2,1),orbit_symflux(L2,3),'.','Color',color_green);
%  plot(orbit_ode45(:,1),orbit_ode45(:,3),'b-','LineWidth',1);
%  plot(orbit_cylindrical([1:nskip_cylindrical:end],1),orbit_cylindrical([1:nskip_cylindrical:end],3),'r.','MarkerSize',5,'MarkerFaceColor','r');
%  plot(orbit_symflux([1:nskip_symflux:end],1),orbit_symflux([1:nskip_symflux:end],3),'Marker','.','Color',color_green,'LineStyle','none','MarkerSize',5,'MarkerFaceColor',color_green);
%  plot(orbit_ode45([1:nskip_ode45_2:end],1),orbit_ode45([1:nskip_ode45_2:end],3),'Marker','.','Color','b','LineStyle','none','MarkerSize',5,'MarkerFaceColor','b');
%  % arrows(start_point(1)+0.5,start_point(2),0.5,270,'Ref',1.5,'FaceColor','k','EdgeColor','k')
%  
%  % Plot additional Markers
%  plot(193.3,5.36,'d','MarkerSize',8,'MarkerFaceColor',color_green,'Color',color_green);
%  plot(193.5,14.41,'^','MarkerSize',8,'MarkerFaceColor','b','Color','b');
%  plot(193.4,12.14,'s','MarkerSize',8,'MarkerFaceColor','r','Color','r');
%  
%  
%  xlim(box_b_xlim);
%  ylim(box_b_ylim);
%  xlabel('R [cm]','Interpreter','latex','FontSize',fsize)
%  ylabel('Z [cm]','Interpreter','latex','FontSize',fsize)
%  hold off
%  title('$\mathbf{c}$','Interpreter','latex','FontSize',fsize)
%  set(gca,'FontSize',fsize_axis)
 
%  print([path_exported_figures,'orbit_comparison'],'-depsc')
