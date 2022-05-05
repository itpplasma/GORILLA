%#######################################################################################################################
% description:
%-----------------------------------------------------------------------------------------------------------------------
% In order to run GORILLA in the WEST geometry with the mesh of
% SOLEDGE3X-EIRENE, two input quantities are needed for GORILLA:
%
% (1) the magnetic field data of an MHD equilibrium in the G-EQDSK format
% (2) the knots and triangles of the corresponding SOLEDGE3X-EIRENE mesh
%
% This script pre-processes the mesh data of SOLEDGE3X-EIRENE
% that is originally provided as a .h5 file. The output of this script is 
% an ASCI file in the .dat format that provides solely the knots and triangles
% information.
%#######################################################################################################################
% author: Michael Eder
% created: 09.05.2022

%% Load data from h5 file

% Relative path and filename of the provided of meshEIRENE file in the .h5 format:
filename_meshEIRENE = 'mesh_for_test.h5';

% Coordinates of the vertices (knots) of the EIRENE mesh
R=hdf5read(filename_meshEIRENE,'/knots/R');
Z=hdf5read(filename_meshEIRENE,'/knots/Z');

% Connectivity matrix: vertex indices for each triangle
TK=hdf5read(filename_meshEIRENE,'/triangles/tri_knots');


%% Write data as .dat file

% knots
fid=fopen('knots.dat','w');
fprintf(fid,'%d %d \n',numel(R),2);
fclose(fid);
writematrix([R,Z],'knots.dat','Delimiter','space','WriteMode','append');


% triangles
fid=fopen('triangles.dat','w');
fprintf(fid,'%d %d \n',numel(TK(:,1)),3);
fclose(fid);
writematrix(TK,'triangles.dat','Delimiter','space','WriteMode','append');
