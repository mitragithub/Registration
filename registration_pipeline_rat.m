function registration_pipeline_rat(atlas_file,input_dir,output_dir,seg_file)
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
config_file = 'rat_nissl_ihc_config.ini';
if nargin == 0
    seg_file = '/cis/home/dtward/Documents/waxholm_rat_atlas/WHS_SD_rat_atlas_v2_pack/WHS_SD_rat_atlas_v2_down2.vtk';
    seg_file = '/cis/home/dtward/Documents/waxholm_rat_atlas/WHS_SD_rat_atlas_v2_pack/WHS_SD_rat_atlas_v2.vtk';
    atlas_file = '/cis/home/dtward/Documents/waxholm_rat_atlas/WHS_SD_rat_atlas_v2_pack/WHS_SD_rat_T2star_v1.01_down2_ss.vtk';

    
    input_dir = '/cis/home/dtward/Documents/csh_data/rat/PTM800/';
    output_dir = 'PTM800/OUTPUT/';
    
    % second try
    input_dir = 'PTM800_BA/PTM800_BA/';
    output_dir = 'PTM800_BA/OUTPUT/';
    
    seg_file = '/mnt/data/waxholm/WHS_SD_rat_atlas_v2_pack/WHS_SD_rat_atlas_v2.vtk';
    atlas_file = '/mnt/data/waxholm/WHS_SD_rat_atlas_v2_pack/WHS_SD_rat_T2star_v1.01_down2_ss.vtk';
    
    seg_file = '/mnt/data/waxholm/WHS_SD_rat_atlas_v2_pack/WHS_SD_rat_atlas_v2_crop.vtk';
    atlas_file = '/mnt/data/waxholm/WHS_SD_rat_atlas_v2_pack/WHS_SD_rat_T2star_v1.01_crop_down2_ss.vtk';
    
    % ideally this step should be done before hand
%     create_location_csv_MBA(input_dir, 14.72, 14.72, 20)
%    keyboard
end
% detailed_output_dir = [output_dir(1:end-1),'_detailed/'];
% Xu is doing this
detailed_output_dir = output_dir;



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 1 is to create an initial slice alignment
% we use a simple initialization where we locate the center of mass of each
% slice and align by translation
pattern = '*-N*.tif';
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
% let's start by matching each fluoro slice to its nearest nissl slice
% note that this is named fluoro, but it supports fluoro or ihc
nissl_pattern = pattern;
files = dir([input_dir '*.tif']);
files = {files.name};
for f = 1 : length(files)
    if strfind(files{f}, '-F')
        fluoro_pattern = '*-F*.tif';
    elseif strfind(files{f},'-IHC')
        fluoro_pattern = '*-IHC*.tif';
    end
end

align_fluoro_to_nissl(input_dir, nissl_pattern, fluoro_pattern, detailed_output_dir);
close all;


% now 3D to 2D transformations for nissl
warning('off','MATLAB:griddedInterpolant:MeshgridEval2DWarnId')
warning('off','MATLAB:griddedInterpolant:MeshgridEval3DWarnId')

ThreeD_to_2D_registration(atlas_file, input_dir, pattern, config_file, detailed_output_dir)
close all;

% interleave these two transformations in appropriate format
combine_nissl_and_fluoro_transforms(detailed_output_dir)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 3 is to generate standard vtk outputs
% something wrong with segmentation file I think
apply_deformation({seg_file,atlas_file}, input_dir, detailed_output_dir, output_dir);
close all;



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
