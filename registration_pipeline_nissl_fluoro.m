function registration_pipeline_nissl_fluoro(atlas_file,input_dir,config_file,detailed_output_dir,output_dir,seg_file)
% run the rna seq mapping pipeline
% 1. specify filename for vtk atlas
% 2. specify directory for input images at 15 um
% they should have a csv file in this directory that indicates their
% geometry, if not, consider running create_location_csv.m first
% this file should have the form
% filename, nx, ny, nz, dx, dy, dz, x0, y0, z0
% 3. specify a config file, see the example
% 4. specify output directory for ALL transformation parameters, they are
% saved as matlab .mat files
% 5. specify output directory for vtk format transformations
% 6. (optional) specify a file containing segmentations, for QC figures only
%
% NOTE all directories should be strings with a trailing slash

close all;
%% 
if nargin < 5
    seg_file = '';
end
%%
% testing on daniels computer
if nargin == 0
    seg_file = '/cis/home/dtward/Documents/ARA/Mouse_CCF/vtk/annotation_25.vtk';
    atlas_file = '/cis/home/dtward/Documents/ARA/Mouse_CCF/vtk/ara_nissl_50.vtk';    
    config_file = '/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/nissl_fluoro_config_v03_rigid.ini';
    
    input_dir = '/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/Xu2Daniel/PMD1238/';
    detailed_output_dir = '/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/PMD1238_rigid_detailed/';
    output_dir ='/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/PMD1238_rigid/';
    
%     keyboard
end

% create geometry
create_location_csv_nissl_fluoro(input_dir, 14.72, 14.72, 20)

%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 1 is to create an initial slice alignment
% we use a simple initialization where we locate the center of mass of each
% slice and align by translation
find_centers_for_initialization_v01_nissl_fluoro(input_dir, detailed_output_dir)
close all;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 2 is to run mapping
% note this currently assumes 0 centered atlas and target
warning('off','MATLAB:griddedInterpolant:MeshgridEval2DWarnId')
warning('off','MATLAB:griddedInterpolant:MeshgridEval3DWarnId')
mouse_map_v24_weightepsilon(atlas_file, input_dir, config_file, detailed_output_dir)
close all;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 3 is to generate standard vtk outputs
apply_deformation_v_0_12_nissl({seg_file,atlas_file}, input_dir, detailed_output_dir, output_dir);
close all;

