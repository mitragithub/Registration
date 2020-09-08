function registration_pipeline_DK39(atlas_file,input_dir,output_dir,seg_file)
% run the mapping pipeline
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
%


close all;
%% 
if nargin < 5
    seg_file = '';
end
%%
% testing on daniels computer
config_file = 'DK39_config.ini';
if nargin == 0
    seg_file = 'atlas_50_vtk/annotation_50.vtk';
    atlas_file = 'atlas_50_vtk/ara_nissl_50.vtk';

  

    input_dir = 'DK39/DK39RGB/';
    output_dir = 'DK39_output/';
    
    
    % metadata
%     Name: DK39
%     Species: mouse
%     Sex: male
%     Genotype: C57/B6J
%     Anesthesia: pentobarbital
%     Perfusing Age in days: 53
%     Exangination method: ACSF
%     Fixative method: Para
%     Side sectioned first: Left
%     Section thickness in microns: 20
%     Orientation: sagittal
%     Mounting: every section
%     Counterstain: NtB
%     Objective: 20x
%     Resolution: 0.325    

    
    
    
    % ideally this step should be done before hand
%     create_location_csv_DK39(input_dir, 0.325*64/4, 0.325*64/4, 20, input_dir, 'png')
    % if I put 10 it will go to 20!
    keyboard
end
% detailed_output_dir = [output_dir(1:end-1),'_detailed/'];
% Xu is doing this
detailed_output_dir = output_dir;



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 1 is to create an initial slice alignment
% we use a simple initialization where we locate the center of mass of each
% slice and align by translation
pattern = '*.png';
find_centers_for_initialization_nissl(input_dir, pattern, detailed_output_dir);
close all;

r = 25;
% this downsampling and iterations is enough for a good initial guess
% not enough for a full accurate reconstruction
downs = [32,16];
niter = 40;
skip_thick = -1; % no thick slices to be skipped
load_initializer = 1;
e = 0.1;
atlas_free_rigid_alignment(input_dir, pattern, detailed_output_dir, r, downs, niter, e, skip_thick, load_initializer)
close all;

% initial affine
downs = [8,4];
niter = 30;
affine_for_initial_alignment(atlas_file, input_dir, pattern, detailed_output_dir, downs, niter)
close all;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 2 is to run mapping
% 3D to 2D transformations for nissl
warning('off','MATLAB:griddedInterpolant:MeshgridEval2DWarnId')
warning('off','MATLAB:griddedInterpolant:MeshgridEval3DWarnId')
pattern = '*DK39*.png';
ThreeD_to_2D_registration(atlas_file, input_dir, pattern, config_file, detailed_output_dir)
close all;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 3 is to generate standard vtk outputs
tic
% apply_deformation({seg_file,atlas_file}, input_dir, detailed_output_dir, output_dir);
apply_deformation_contours({seg_file,atlas_file}, input_dir, detailed_output_dir, output_dir);
% apply_deformation({atlas_file}, input_dir, detailed_output_dir, output_dir);
close all;
toc



return
%%
% step 4, if desired, is to run in edit mode
% edit mode is enabled by including segmentations in the ThreeD to 2D, 
% it will read segmentations in the input directory, these will be binary
% tifs
% then it will update
ThreeD_to_2D_registration({atlas_file,seg_file}, input_dir, pattern, config_file, detailed_output_dir)
combine_nissl_and_fluoro_transforms(detailed_output_dir)
apply_deformation({seg_file,atlas_file}, input_dir, detailed_output_dir, output_dir);
ara_nissl_50