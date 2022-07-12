%#######################################################################################################################
% description:
%-----------------------------------------------------------------------------------------------------------------------
% In order to run GORILLA in the WEST geometry with the mesh of
% SOLEDGE3X-EIRENE, two input quantities are needed for GORILLA:
%
% (1) the magnetic field data of an MHD equilibrium in the G-EQDSK format
% (2) the knots and triangles of the corresponding SOLEDGE3X-EIRENE mesh
%
% This script pre-processes the magnetic field data of an MHD equilibrium
% that is originally computed by the FEEQS code and is provided as a .mat
% file. The output of this script is in the G-EQDSK format.
%#######################################################################################################################
% author: Michael Eder
% created: 06.05.2022

% Relative path and filename of the provided MHD equilibrium in the .mat format
input_file = 'mag_shot_for_test.mat';

% Name of the output of the MHD equilibrium in the G-EQDSK format
output_name = 'g_file_for_test_WEST';

% Load input
load(input_file);

% Size of rectangular regular grid on which magnetic field data is provided
[nz,nr] = size(Br2D);

% Define format specs for output
form_spec = '%f ';

% Write output in a new file in the G-EQDSK format
fid=fopen(output_name,'w');

% 1) nr nz (two integers separated by space - size of the grid)
fprintf(fid,'%d %d \n',[nr,nz]);

% 2) Bphi2D(1,1)*r2D(1,1) (scalar - a double precision number)
fprintf(fid,[form_spec,'\n'], Bphi2D(1,1)*r2D(1,1));

% 3) r2D(1,:) (column vector, written as one line)
for k = 1:numel(r2D(1,:))
    fprintf(fid,form_spec, r2D(1,k));
end
fprintf(fid,'\n');

% 4) z2D(:,1) (line vector, written as one line)
for k = 1:numel(z2D(:,1))
    fprintf(fid,form_spec, z2D(k,1));
end
fprintf(fid,'\n');

% 5 -- 5+nr) flux2D(:,:) (matrix of flux2D)
for kr = 1:nr
    for kz = 1:nz
        fprintf(fid,form_spec, flux2D(kz,kr));
    end
    fprintf(fid,'\n');
end

% Close file
fclose(fid)
