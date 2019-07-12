function registration_pipeline_rnaseq(atlas_file,input_dir,output_dir,seg_file)
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
    detailed_output_dir = '710_test_detailed/';
    output_dir = '710_test/';
    config_file = 'rnaseq_config.ini';
    
    % my test for 710 was a total failure!
    % the affine and deformation kept mving in opposite directions, and the
    % result was terrible
    % also the energy did not converge (reg energy kept growing)
    % let's start by checking 720, which worked previously
    detailed_output_dir = '720_test_detailed/';
    output_dir = '720_test/';
    
    % go back to 710, with the dramatic failure
    input_dir = '/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/toDaniel/MD710/';
    output_dir = '710_test_01_sooner_diff/';
    config_file = 'rnaseq_config_test_sooner_diffeo.ini';
    
    % next one, set post affine reduce to 0
    output_dir = '710_test_02_postaffine0/';
    config_file = 'rnaseq_config_test_postaffine0.ini';

    % third one, affine only for a long time
    output_dir = '710_test_03_affine_only/';
    config_file = 'rnaseq_config_test_03_affine_only.ini';

    % note before I ran this I founda  bug in my 3D z gradient
    % it only affected first and last slices, so I doubt it was a problem
    % fourth one, freeze diffeo after 1000
%     output_dir = '710_test_04_diffeo_freeze/';
%     config_file = 'rnaseq_config_test_sooner_diffeo.ini';

    % the above never finished because I got a nan, I will repeat it
    
    % now let's change the way I'm computing the gradient of the atlas for
    % the affine update, this is because even for affine only it seemed
    % like there were errors
    output_dir = '710_test_05_affine_fix/';
    config_file = 'rnaseq_config_test_sooner_diffeo.ini';
    
    
    output_dir = '710_test_06_final/';
    config_file = 'rnaseq_config.ini';

    % ideally this step should be done before hand
    create_location_csv_rnaseq(input_dir, 14.72, 14.72, 20, 200,[40,41,42,121:126,133:138,214:219])
    
    keyboard
end

detailed_output_dir = [output_dir(1:end-1),'_detailed/'];


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 1 is to create an initial slice alignment
% we use a simple initialization where we locate the center of mass of each
% slice and align by translation
find_centers_for_initialization_nissl(input_dir, detailed_output_dir)
close all;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 2 is to run mapping
warning('off','MATLAB:griddedInterpolant:MeshgridEval2DWarnId')
warning('off','MATLAB:griddedInterpolant:MeshgridEval3DWarnId')
nonrigid_thick_only = 1;
ThreeD_to_2D_registration(atlas_file, input_dir, config_file, detailed_output_dir,nonrigid_thick_only)
% ThreeD_to_2D_registration_diffeo_freeze(atlas_file, input_dir, config_file, detailed_output_dir,nonrigid_thick_only)
close all;



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 3 is to generate standard vtk outputs
apply_deformation({seg_file,atlas_file}, input_dir, detailed_output_dir, output_dir);
close all;
