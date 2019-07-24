function registration_pipeline_rnaseq_atlas_free_align(atlas_file,input_dir,output_dir,seg_file)
% run the rna seq mapping pipeline
% 1. specify filename for vtk atlas
% 2. specify directory for input images at 15 um
% they should have a csv file in this directory that indicates their
% geometry, if not, consider running create_location_csv.m first
% this file should have the form
% filename, nx, ny, nz, dx, dy, dz, x0, y0, z0
% 3. specify output directory for ALL transformation parameters, they are
% saved as matlab .mat files
% 4. specify output directory for vtk format transformations
% 5. (optional) specify a file containing segmentations, for QC figures only
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
    seg_file = '/cis/home/dtward/Documents/ARA/Mouse_CCF/vtk/annotation_50.vtk';
    atlas_file = '/cis/home/dtward/Documents/ARA/Mouse_CCF/vtk/ara_nissl_50.vtk';
    input_dir = '/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/toDaniel/MD710/';
    output_dir = '710_test/';
    config_file = 'rnaseq_config.ini';
    
    output_dir = '710_test2/';
    
    % ideally this step should be done before hand
%     create_location_csv_rnaseq(input_dir, 14.72, 14.72, 20, 200,[40,41,42,121:126,133:138,214:219])
    
    
    input_dir = '/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/toDaniel/MD705/';
    output_dir = '705_test0/';
%     create_location_csv_rnaseq(input_dir, 14.72, 14.72, 20, 250,[17,107,425])

    input_dir = '/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/718/';
    output_dir = '718_test0/';
    
    input_dir = '/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/719/';
    output_dir = '719_test0/';


    keyboard
end
detailed_output_dir = [output_dir(1:end-1),'_detailed/'];



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 1 is to create an initial slice alignment
% we use a simple initialization where we locate the center of mass of each
% slice and align by translation
% then we update with slice to slice alignment, exclude thick
% then we will get an initial affine alignment

pattern = '*.tif';
find_centers_for_initialization_nissl(input_dir, pattern, detailed_output_dir); 

r = 25;
% this downsampling and iterations is enough for a good initial guess
% not enough for a full accurate reconstruction
downs = [32,16];
niter = 40;
et_factor = 1e-4;
etheta_factor = 1e-11;
skip_thick = 25; % skip thick
load_initializer = 1;
atlas_free_rigid_alignment(input_dir, pattern, detailed_output_dir, r, downs, niter, et_factor, etheta_factor, skip_thick, load_initializer)
close all;

% initial affine
downs = [8,4];
niter = 30;
et_factor = 2e-7;
el_factor = 1e-14;
affine_for_initial_alignment(atlas_file, input_dir, pattern, detailed_output_dir, downs, niter, et_factor, el_factor)
close all;



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 2 is to run mapping
warning('off','MATLAB:griddedInterpolant:MeshgridEval2DWarnId')
warning('off','MATLAB:griddedInterpolant:MeshgridEval3DWarnId')
nonrigid_thick_only = 1;
ThreeD_to_2D_registration(atlas_file, input_dir, pattern, config_file, detailed_output_dir,nonrigid_thick_only)
close all;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 3 is to generate standard vtk outputs
apply_deformation({seg_file,atlas_file}, input_dir, detailed_output_dir, output_dir);
close all;
