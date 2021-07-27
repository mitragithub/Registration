% Standard registration pipeline with settings for:
%     marmoset m6328
% All options are set in the first cell and this file is run as a script
% All 3D images are expected to be in vtk format
% 2D images are expected to lie in a single directory and described with a
% geometry csv file.
%
% Overview: 
% 1. All 3D data is collected into a common space.
% 2. An initial alignment of 2D data is performed, followed by a 3D to 2D
% reconstruction
% 3. Outputs are generated including all images mapped to all spaces.
%
%
% arguments:
% output_prefix: where to put outputs, this should be a directory ending in a slash
%
% data_2D_directory: where to find slice images and geometry
%
% data_2D: a cell array.  Each cell contains a cell array with two
% elements, the first containing a dataset name, and the second containing
% a glob file matching pattern
%
% data_2D_register: an index for which of the above datasets should be used
% for registration to 3D data
%
% config_file: a .ini config file for 3D to 2D mapping
%
% A0: an initial guess of affine alignment.  Should account at least for 90
% degree rotatoins.
%
% data_3D_files: A cell of cells of cells.  data_3D_files{1} describes data
% in space 1.  data_3D_files{1}{1} describes image 1 in space 1.  It should
% be a cell array pair with data_3D_files{1}{1}{1} the image name and
% data_3D_files{1}{1}{2} the file location.
% 
% data_3D_spaces: A cell array containing names for all of the above
% spaces (not images)
%
% data_3D_common: An index to say which is the "common space" that all 3D
% data will be collected in.
%
% mapping_3D_pairs: an array of triples of indices.  The first and second
% says which space and image we should map from, the third says which image 
% in the common space to map to.
%
% mapping_3D_2D: a pair saying which dataset should be mapped onto the 2D
% data, the first element is the space and the second is the image.  Note
% it is a good idea for this space to be the "common space".
%
% mapping_3D_options: a cell array of structs for options to pass to the 3D
% mapping algorithm.
%
% preprocessing_3D_options: a cell of cells containing options
% corresponding to preprocessing for each 3D mapping we compute
%
% dx_contour: a size to output coordinate grids
%
%
% in v02b I will try to improve the 3D mappings particularly on the
% anterior surface
% 1229v02b work started on July 13 2021
%
% TODO: try to use the atlas nissl as the mapping one instead of the atlas
% mri
% the labels don't align that well.

addpath Functions/plotting
addpath Functions/downsample
addpath Functions/vtk



clear all;
close all;
fclose all;

% stop for qc?
stop_for_qc = true;

% where to put outputs
output_prefix = 'm1229_v02b_output_v04/';

% 2D data, a directory
data_2D_directory = '/home/dtward/data/csh_data/marmoset/m1229/slices/';
data_2D = {
    {'nissl','*-N*'},
    {'fluoro','*-F*'},
    {'myelin','*-M*'}
}; % how to find the different modalities
data_2D_register = 1; % index for which of the above is for registration
config_file = 'marmoset_nissl_config_test3.ini'; %

% for marmoset initial affine (this is default)
A0 = eye(4);
A0 = [0,1,0,0;
    1,0,0,0;
    0,0,-1,0;
    0,0,0,1]*A0;


% for 3D data, separate into spaces
data_3D_files = {
    {
        {'atlas-mri','/home/dtward/data/csh_data/marmoset/Woodward_2018/bma-1-mri-reorient.vtk'},
        {'atlas-nissl','/home/dtward/data/csh_data/marmoset/Woodward_2018/bma-1-nissl-reorient.vtk'},
        {'atlas-seg','/home/dtward/data/csh_data/marmoset/Woodward_2018/bma-1-region_seg-reorient.vtk'}
    },
    {
        {'exvivo-mri','/home/dtward/data/csh_data/marmoset/m1229/M1229MRI/MRI/exvivo/HR_T2/HR_T2_CM1229F-reorient.vtk'}
    },
    {
        {'invivo-mri','/home/dtward/data/csh_data/marmoset/m1229/M1229MRI/MRI/invivo/T2/T2_CM1129F-reorient.vtk'}
    },
    {
        {'invivo-mri-quant','/home/dtward/data/csh_data/marmoset/m1229/m1229_MRImaps/invivo/t1/T1mapping_T1-reorient.vtk'}
    },
};
% I also need to name the spaces
data_3D_spaces = {'atlas','exvivo','invivo','invivo-quant'};
data_3D_common = 2; % all 3D data gets registered to this space, here ex vivo
% mapping pairs
mapping_3D_pairs = [
    [1,2,1], % 2 for nissl
    [3,1,1],
    [4,1,1]
];
mapping_3D_2D = [2,1]; % which dataset to map onto 2D

% map mri to mri
% input options as a struct
% use common_to_target = True or False, mandatory
% to compute a map in one direction versus the other
% in this example I want to compute the map atlas to ex vivo, but common is
% ex vivo, so I put false
% for now just do true
opts = struct('common_to_target',true,...
        'downs',[4,2,1],'niters',round([100,50,25]*4),'niters0',[1,1,1],...
        'eA',0.5,'eV',2e4,...
        'nt',5,'a',400,'apfactor',2,'p',2,...
        'sigmaR',2e3,'sigmaM',0.5,'sigmaA',2.5,'sigmaB',1.0,...
        'order',3,'prior',[0.9,0.05,0.05]);
opts = struct('common_to_target',true,...
        'downs',[4,2,1],'niters',round([100,50,25]*4),'niters0',[1,1,1],...
        'eA',0.5,'eV',2e4,...
        'nt',5,'a',400,'apfactor',2,'p',2,...
        'sigmaR',1e3,'sigmaM',0.2,'sigmaA',1.0,'sigmaB',0.4,...
        'order',1,'prior',[0.9,0.05,0.05]); % try order 2, maybe avoid crushing the contrast so much
opts = struct('common_to_target',true,...
        'downs',[4,2,1],'niters',round([100,50,25]*4),'niters0',[1,1,1],...
        'eA',0.5,'eV',1e4,...
        'nt',5,'a',400,'apfactor',2,'p',2,...
        'sigmaR',1e3,'sigmaM',0.2,'sigmaA',1.0,'sigmaB',0.4,...
        'order',3,'prior',[0.9,0.05,0.05]); % try nissl, back to order 3
    
opts = struct('common_to_target',true,...
        'downs',[2,1],'niters',round([50,25]*4),'niters0',[1,1,1],...
        'eA',0.5,'eV',1e4,...
        'nt',5,'a',400,'apfactor',2,'p',2,...
        'sigmaR',1e3,'sigmaM',0.2,'sigmaA',1.0,'sigmaB',0.4,...
        'order',1,'prior',[0.9,0.05,0.05]); % try nissl, back to order 3, try 1
opts = struct('common_to_target',true,...
        'downs',[1],'niters',round([50]*4),'niters0',[1,1,1],...
        'eA',0.5,'eV',1e4,...
        'nt',5,'a',400,'apfactor',2,'p',2,...
        'sigmaR',1e3,'sigmaM',0.5,'sigmaA',2.5,'sigmaB',1.0,...
        'order',1,'prior',[0.9,0.05,0.05]); % try nissl, back to order 3, try 1

% trying with more iterations    
opts = struct('common_to_target',true,...
        'downs',[1],'niters',round([500]),'niters0',[1,1,1],...
        'eA',0.5,'eV',1e4,...
        'nt',5,'a',400,'apfactor',2,'p',2,...
        'sigmaR',1e3,'sigmaM',0.5,'sigmaA',2.5,'sigmaB',1.0,...
        'order',1,'prior',[0.9,0.05,0.05]); % try nissl, back to order 3, try 1
% next try with order 2 with a new name
% also try a bit of bias correction on the target (i.e. the nissl)    
% next test is a little less reg, v03
opts = struct('common_to_target',true,...
        'downs',[1],'niters',round([500]),'niters0',[1,1,1],...
        'eA',0.5,'eV',1e4,...
        'nt',5,'a',400,'apfactor',2,'p',2,...
        'sigmaR',2e3,'sigmaM',0.5,'sigmaA',2.5,'sigmaB',1.0,...
        'order',1,'prior',[0.9,0.05,0.05]); % try nissl, back to order 3, try 1
% this time I want to give other compartments a chance
opts = struct('common_to_target',true,...
        'downs',[1],'niters',round([500]),'niters0',[1,1,1],...
        'eA',0.5,'eV',1e4,...
        'nt',5,'a',400,'apfactor',2,'p',2,...
        'sigmaR',2e3,'sigmaM',0.4,'sigmaA',2.0,'sigmaB',0.8,...
        'order',1,'prior',[0.9,0.05,0.05]); % try nissl, back to order 3, try 1
% maybe go back to more reg, yes this one was not so good, too wiggly,
% going back to more reg
opts = struct('common_to_target',true,...
        'downs',[1],'niters',round([500]),'niters0',[1,1,1],...
        'eA',0.5,'eV',1e4,...
        'nt',5,'a',400,'apfactor',1,'p',2,...
        'sigmaR',1e3,'sigmaM',0.4,'sigmaA',2.0,'sigmaB',0.8,...
        'order',1,'prior',[0.9,0.05,0.05]); % try nissl, back to order 3, try 1
% okay I think the above is quite good
% I'd like to move along
mapping_3D_options = {
    opts,
    opts,
    opts
    };

preprocessing_3D_options = {
    {
        struct('name','flip_contrast','which','target','lower_threshold',20,'upper_threshold','none'),
%         struct('name','power','which','target','power',0.5),
        struct('name','resample_isotropic'),
%         struct('name','bias_correct','which','atlas','scale',1500,'n',4), % was 2 by default, was 2000
%         struct('name','mask','which','atlas','n',4),
        struct('name','pad','which','atlas','n',10),
    }, % pretty happy with mapping 1, ex vivo to atlas
    {
        struct('name','resample_isotropic'),
%         struct('name','bias_correct','which','atlas','scale',2000),
%         struct('name','mask','which','atlas'),
        struct('name','pad','which','atlas','n',10),
    },  % mapping 2 is ex vivo to in vivo
    {
        struct('name','resample_isotropic'),
%         struct('name','mask','which','atlas'),
        struct('name','pad','which','atlas','n',10),
    },
    
};
% for outputting coordinate grids
dx_contour = 1000;

% make an output dir
out_dir = './';
[a,b,c] = fileparts(output_prefix);
if ~isempty(a)
    out_dir = [a '/'];
end
if ~isempty(a) && ~exist(a,'dir')
    mkdir(a)
end


%%
% start my 3D mapping, begin loop
close all;
for m = 1 : size(mapping_3D_pairs,1)
    
    pairs = mapping_3D_pairs(m,:);
    
    
    % get names for pairs
    atlas_name = data_3D_files{data_3D_common}{pairs(3)}{1};
    atlas_fname = data_3D_files{data_3D_common}{pairs(3)}{2};
    target_name = data_3D_files{pairs(1)}{pairs(2)}{1};
    target_fname = data_3D_files{pairs(1)}{pairs(2)}{2};
    for i = 1 : 5
        disp(' ')
    end
    disp(['Starting mapping pair ' atlas_name ' to ' target_name])
    
    % if order calls for it, swap names
    if ~mapping_3D_options{m}.common_to_target
        atlas_name_ = atlas_name;
        atlas_fname_ = atlas_fname;
        atlas_name = target_name;
        atlas_fname = target_fname;
        target_name = atlas_name_;
        target_fname = atlas_fname_;
    end
    
    % load images
    [xI,yI,zI,I,titleI,namesI] = read_vtk_image(atlas_fname);
    nxI = [size(I,2),size(I,1),size(I,3)];
    dxI = [xI(2)-xI(1),yI(2)-yI(1),zI(2)-zI(1)];
    [xJ,yJ,zJ,J,titleJ,namesJ] = read_vtk_image(target_fname);
    nxJ = [size(J,2),size(J,1),size(J,3)];
    dxJ = [xJ(2)-xJ(1),yJ(2)-yJ(1),zJ(2)-zJ(1)];
    
    % start preprocessing, each time I'll update one or both of I and J
    for p = 1 : length(preprocessing_3D_options{m})
        close all;
        opts = preprocessing_3D_options{m}{p};
        disp(['Starting preprocessing step ' num2str(p)])
        % 
        name = opts.name;
        if strcmpi(name,'resample_isotropic')
            prefix_outdir = [out_dir 'preprocessing_' num2str(p,'%03d') '_resample_isotropic/'];
            mkdir([prefix_outdir])
            disp(['isotropic resampling'])
            % isotropic resampling
            % we downsample both images
            % we try to get them about the same resolution
            % and about isotropic
            % we will choose the largest dimension
            desired_dx = max([dxI,dxJ]);
            % downsample them by the closest integer factor
            factors_I = round(desired_dx ./ dxI);
            factors_J = round(desired_dx ./ dxJ);
            disp(['Downsampling atlas by factor ' num2str(factors_I)])
            disp(['Downsampling target by factor ' num2str(factors_J)])
            
            [xI,yI,zI,I] = downsample(xI,yI,zI,I,factors_I);
            nxI = [size(I,2),size(I,1),size(I,3)];
            dxI = [xI(2)-xI(1),yI(2)-yI(1),zI(2)-zI(1)];
            
            [xJ,yJ,zJ,J] = downsample(xJ,yJ,zJ,J,factors_J);
            nxJ = [size(J,2),size(J,1),size(J,3)];
            dxJ = [xJ(2)-xJ(1),yJ(2)-yJ(1),zJ(2)-zJ(1)];
            
            danfigure(1);
            sliceView(xI,yI,zI,I)
            title('Resampled Atlas')
            saveas(gcf,[prefix_outdir 'resampled_I.png'])
            danfigure(2);
            sliceView(xJ,yJ,zJ,J)
            title('Resampled Target')
            saveas(gcf,[prefix_outdir 'resampled_J.png'])
        elseif strcmpi(name,'power')    
            for a = [0,1]
                if a == 0 && (strcmpi(opts.which,'atlas') || strcmpi(opts.which,'both'))
                    I = I.^opts.power;
                elseif a == 1 && (strcmpi(opts.which,'target') || strcmpi(opts.which,'both'))
                    J = J.^opts.power;
                end
            end
        elseif strcmpi(name,'flip_contrast')
            % linearly maps lower threshold to upper and vice versa
            % outside of threshold leave as is
            lt = opts.lower_threshold;
            ut = opts.upper_threshold;
            for a = [0,1]
                if a == 0 && (strcmpi(opts.which,'atlas') || strcmpi(opts.which,'both'))
                    if strcmp(lt,'none')
                        lt = min(I(:));
                    end
                    if strcmp(ut,'none')
                        ut = max(I(:));
                    end                    
                    ind = (I>lt & I<ut);                    
                    I(ind) = lt + ut - I(ind);
                elseif a == 1 && (strcmpi(opts.which,'target') || strcmpi(opts.which,'both'))
                    if strcmp(lt,'none')
                        lt = min(J(:));
                    end
                    if strcmp(ut,'none')
                        ut = max(J(:));
                    end                    
                    ind = (J>lt & J<ut);                    
                    J(ind) = lt + ut - J(ind);
                end
                
                
            end
        elseif strcmpi(name,'bias_correct')
            disp('Bias correct')
            prefix_outdir = [out_dir 'preprocessing_' num2str(p,'%03d') '_bias_correct/'];
            mkdir([prefix_outdir])
            
            for a = [0,1] % 0 atlas, 1 target
                if a == 0 && (strcmpi(opts.which,'atlas') || strcmpi(opts.which,'both'))
                    im = I;
                    x = xI;
                    y = yI;
                    z = zI;                    

                elseif a == 1 && (strcmpi(opts.which,'target') || strcmpi(opts.which,'both'))
                    im = J;
                    x = xJ;
                    y = yJ;
                    z = zJ;
                end
                opt_ = struct;
                if ~isfield(opts,'n')
                    opts.n = 2;
                end
                opt_.M = opts.n;
                mask = kmeans(im(:),opt_); % note this is my kmeans, not builtin
                mask = reshape(sum(mask(:,1,2:end),3),size(im)); % fg
                r = opts.scale;
               
                % pad the image
                dx = [(x(2)-x(1)), (y(2)-y(1)), (z(2)-z(1))];
                topad = round(r./dx*2);
                Ip = padarray(im,topad([2,1,3]),0,'both');
                maskp = padarray(mask,topad([2,1,3]),0,'both');

                x_ = 0 : (x(2)-x(1)) : 3*r;
                y_ = 0 : (y(2)-y(1)) : 3*r;
                z_ = 0 : (z(2)-z(1)) : 3*r;
                x_ = [-x_(end:-1:2),x_];
                y_ = [-y_(end:-1:2),y_];
                z_ = [-z_(end:-1:2),z_];                
                [X_,Y_,Z_] = meshgrid(x_,y_,z_);
                niter = 10;
                K = exp(-(X_.^2 + Y_.^2 + Z_.^2)/2.0/(r/niter).^2);
                K = K / sum(K(:));                

                
                Kp = padarray(K,size(Ip)-size(K),'post');
                Kp = circshift(Kp,-[(length(y_)-1),(length(x_)-1),(length(z_)-1)]/2);
                Kphat = fftn(Kp);

                logIp = log(Ip);
                bc = 0;
                logIp(~maskp) = bc;
                
                
                % blur it with boundary conditions
                maskphat = fftn(maskp);
                maskblur = ifftn(maskphat .* Kphat,'symmetric');
                logIpblur = logIp;
                for it = 1 : niter
                    logIpblurhat = fftn(logIpblur);
                    logIpblur = ifftn(logIpblurhat .* Kphat,'symmetric');
                    logIpblur = logIpblur./maskblur;
                    logIpblur(~maskp) = bc;
                    sliceView(logIpblur)
                    drawnow;
                end
                % highpass
                logIph = logIp - logIpblur;
                Iph = exp(logIph);
                Iph(~maskp) = min(Iph(maskp>0));
                
                
                % remove the padding and reassign
                im_ = Iph(topad(2)+1:end-topad(2),topad(1)+1:end-topad(1),topad(3)+1:end-topad(3));
                sliceView(x,y,z,im_)
                drawnow;
                figure;sliceView(x,y,z,im)
                
                if a == 0 && (strcmpi(opts.which,'atlas') || strcmpi(opts.which,'both'))
                    I = im_;
                    saveas(gcf,[prefix_outdir 'corrected_I.png'])
                elseif a == 1 && (strcmpi(opts.which,'target') || strcmpi(opts.which,'both'))
                    J = im_;
                    saveas(gcf,[prefix_outdir 'corrected_J.png'])
                end
                
                
                

            end
            
        elseif strcmpi(name,'pad')
            disp('Pad')
            prefix_outdir = [out_dir 'preprocessing_' num2str(p,'%03d') '_pad/'];
            mkdir([prefix_outdir])

            if (strcmpi(opts.which,'atlas') || strcmpi(opts.which,'both'))
                I = padarray(I,opts.n*[1,1,1],0,'both');
                for i = 1 : opts.n
                    xI = [xI(1)-dxI(1), xI, xI(end)+dxI(1)];
                    yI = [yI(1)-dxI(2), yI, yI(end)+dxI(2)];
                    zI = [zI(1)-dxI(3), zI, zI(end)+dxI(3)];
                end
                danfigure(1)
                sliceView(xI,yI,zI,I)
                saveas(gcf,[prefix_outdir 'padded_I.png'])
            end
            if (strcmpi(opts.which,'target') || strcmpi(opts.which,'both'))
                J = padarray(J,opts.n*[1,1,1],0,'both');
                for i = 1 : opts.n
                    xJ = [xJ(1)-dxJ(1), xJ, xJ(end)+dxJ(1)];
                    yJ = [yJ(1)-dxJ(2), yJ, yJ(end)+dxJ(2)];
                    zJ = [zJ(1)-dxJ(3), zJ, zJ(end)+dxJ(3)];
                end
                danfigure(1)
                sliceView(xJ,yJ,zJ,J)
                saveas(gcf,[prefix_outdir 'padded_J.png'])
                
            end
        elseif strcmpi(name,'mask')
            disp('Mask')
            prefix_outdir = [out_dir 'preprocessing_' num2str(p,'%03d') '_mask/'];
            mkdir([prefix_outdir])
            
            if (strcmpi(opts.which,'atlas') || strcmpi(opts.which,'both'))
                opt_ = struct;
                if ~isfield(opts,'n')
                    opts.n = 2;
                end
                opt_.M = opts.n;
                mask = kmeans(I(:),opt_); % note this is my kmeans, not builtin
                mask = reshape(sum(mask(:,1,2:end),3),size(im)); % fg
                I = I .*mask;
                danfigure(1)
                sliceView(xI,yI,zI,I)
                saveas(gcf,[prefix_outdir 'masked_I.png'])
                
            end
            if (strcmpi(opts.which,'target') || strcmpi(opts.which,'both'))
                opt_ = struct;
                opt_.M = opts.n;
                mask = kmeans(J(:),opt_); % note this is my kmeans, not builtin
                mask = reshape(sum(mask(:,1,2:end),3),size(im)); % fg
                J = J .*mask;
                
                danfigure(1)
                sliceView(xJ,yJ,zJ,J)
                saveas(gcf,[prefix_outdir 'masked_J.png'])

            end
        end
        
        
    end % end of preprocessing
    
    
    %%    
    % mapping
    % mapping parameters
    opt = mapping_3D_options{m};
    downs = opt.downs;
    niters = opt.niters;
    
    eV = opt.eV;
    eA = opt.eA;



    % other params
    nt = opt.nt; % number of timesteps
    a = opt.a; % one voxel is 1-2 hundred micron
    apfactor = opt.apfactor; % 2 voxels
    p = opt.p; % power of linaer operators

    sigmaR = opt.sigmaR;

    % this are for z scored data
    sigmaM = opt.sigmaM;
    sigmaA = opt.sigmaA;
    sigmaB = opt.sigmaB;

    % contrast
    order = opt.order;

    A0_ = eye(4);
    prior = opt.prior;

    % let's look carefully at my images and summarize what I see
    danfigure(1)
    sliceView(xI,yI,zI,I)
    danfigure(2)
    sliceView(xJ,yJ,zJ,J)
    % first a lot of stuff has been masked out in I
    % this happened at the bias correction stage
    % 
    
    % run mapping
    OPT = struct;
    OPT.fAphiIgif = 1;
    OPT.errgif = 1;
    [phiiAix,phiiAiy,phiiAiz,Aphix,Aphiy,Aphiz,A,vtx,vty,vtz] = ThreeD_to_3D_registration(...
        xI,yI,zI,I,...
        xJ,yJ,zJ,J,...
        A0_,...
        nt,a,p,...
        sigmaR,sigmaM,sigmaA,sigmaB,prior,...
        order,...
        eA,eV*1,apfactor,...
        downs,niters,OPT);
    saveas(6,'energyout.png');
    
    % test
    % this worked pretty well
%     OPT = struct;
%     OPT.MI = 1;
%     [phiiAix,phiiAiy,phiiAiz,Aphix,Aphiy,Aphiz,A,vtx,vty,vtz] = ThreeD_to_3D_registration(...
%         xI,yI,zI,I,...
%         xJ,yJ,zJ,J,...
%         A0_,...
%         nt,a*0.5,p,...
%         sigmaR*5e6,sigmaM,sigmaA,sigmaB,prior,...
%         order,...
%         eA*1e5,eV*2e5,apfactor,...
%         downs,niters,OPT);

    
    % trying below with smaller a
    % this ended up really bad
%     OPT = struct;
%     OPT.MI = 1;
%     OPT.preconditioner = 0;
%     [phiiAix,phiiAiy,phiiAiz,Aphix,Aphiy,Aphiz,A,vtx,vty,vtz] = ThreeD_to_3D_registration(...
%         xI,yI,zI,I,...
%         xJ,yJ,zJ,J,...
%         A0_,...
%         nt,a*0.25,p,...
%         sigmaR*5e6,sigmaM,sigmaA,sigmaB,prior,...
%         order,...
%         eA*5e4,eV*5e6,apfactor,...
%         downs,niters,OPT);
    
    xV = xI;
    yV = yI;
    zV = zI;
    
    F = griddedInterpolant({yI,xI,zI},I,'linear','nearest');
    AphiI = F(phiiAiy,phiiAix,phiiAiz);
    AphiI_ = (AphiI - min(AphiI(:)))/(max(AphiI(:)) - min(AphiI(:)));
    J_ = (J - min(J(:)))/(max(J(:)) - min(J(:)));
    danfigure(999)
    sliceView(xJ,yJ,zJ,cat(4,AphiI_,J_,AphiI_))
    % we need to validate this carefully

    % the result came up wiggedy, e.g. the midline is wiggly
    % but it didn't look wiggly when I was registering
    % is there some problem? need to double check
    % is there a problem with downsampling? is xV really xI? even when I
    % downsampled? I think it is since it is isotropic downsampling and the
    % last number is 1
    % ok it was getting wigly on the midline
    
    
    
    %% 
    % swap transforms if necessary and save (no swapping supported yet because of order of A and phi)
    prefix_outdir = [out_dir 'mapping_3D_' num2str(m,'%03d') '_' atlas_name '_to_' target_name '/'];
    mkdir([prefix_outdir])
    save([prefix_outdir 'output.mat'],'phiiAix','phiiAiy','phiiAiz','Aphix','Aphiy','Aphiz','A','vtx','vty','vtz','xI','yI','zI','I','xJ','yJ','zJ','J','xV','yV','zV')
    
    %%
    
    % transform all data to common space (at native resolution)
    % make sure vtk data does not have any nans in it as this won't load in
    % other viewers
    % load it all
%     close all;
    prefix_outdir = [out_dir 'common_space_3D_data/'];
    mkdir([prefix_outdir])
    % note my transform is saved at coordinates xI
    % so loop through all the data in atlas space and apply transforms
    % how should I sample it?
    % atlas space but native resolution
    
    % note this should not be data_3D_files{m}, 
    % it should come from the pairs variable
    for i = 1 : length(data_3D_files{pairs(1)})
        name = data_3D_files{pairs(1)}{i}{1};
        fname = data_3D_files{pairs(1)}{i}{2};
        outname = [num2str(pairs(1),'%03d') '_' num2str(i,'%03d') '_' name '_common.vtk'];
        
        [x,y,z,im,title_,names] = read_vtk_image(fname);
        % sampling domain
        xS = min(xI) : x(2)-x(1) : max(xI);
        yS = min(yI) : y(2)-y(1) : max(yI);
        zS = min(zI) : z(2)-z(1) : max(zI);
        [XS,YS,ZS] = meshgrid(xS,yS,zS);
        % transform target back to atlas using Aphix
        F = griddedInterpolant({yV,xV,zV},Aphix,'linear','nearest');
        mapx = F(YS,XS,ZS);
        F = griddedInterpolant({yV,xV,zV},Aphiy,'linear','nearest');
        mapy = F(YS,XS,ZS);
        F = griddedInterpolant({yV,xV,zV},Aphiz,'linear','nearest');
        mapz = F(YS,XS,ZS);
        
        if contains(name,'seg') 
            F = griddedInterpolant({y,x,z},im,'nearest','none');
            mapim = F(mapy,mapx,mapz);
            % note this will not work with nans
            mapim(isnan(mapim)) = 0;
            if max(im(:)) > 2^16-1
                write_vtk_image(xS,yS,zS,uint32(mapim),[prefix_outdir outname],[title_ '_to_common'],names)
            else
                write_vtk_image(xS,yS,zS,uint16(mapim),[prefix_outdir outname],[title_ '_to_common'],names)
            end
            
            danfigure(111);
            sliceView(xS,yS,zS,mod(mapim,32));
            saveas(gcf,[prefix_outdir outname(1:end-4) '.png'])
        else
            F = griddedInterpolant({y,x,z},im,'linear','none');
            mapim = F(mapy,mapx,mapz);
            % note this will not work with nans
            mapim(isnan(mapim)) = 0;
            write_vtk_image(xS,yS,zS,mapim,[prefix_outdir outname],[title_ '_to_common'],names)
            danfigure(1111);
            sliceView(xS,yS,zS,mapim);
            saveas(gcf,[prefix_outdir outname(1:end-4) '.png'])
        end
    end
    
    
    
    
end

% copy over the atlas
prefix_outdir = [out_dir 'common_space_3D_data/'];
mkdir([prefix_outdir])
for i = 1 : length(data_3D_files{data_3D_common})
    name = data_3D_files{data_3D_common}{i}{1};
    fname = data_3D_files{data_3D_common}{i}{2};
    outname = [num2str(data_3D_common,'%03d') '_' num2str(i,'%03d') '_' name '_common.vtk'];
    copyfile(fname,[prefix_outdir outname])
end



%%
% now we begin 2D mapping

% create geometry file, % this should already be done
% % pixel_size = (0.46*2)*32
% pixel_size = (0.46*2)*64; % I think I need bigger
% note slice thickness is 20, slice spacing is 4 times as much
% create_location_csv_marmoset_coronal_v02(data_2D_directory, pixel_size, pixel_size, 20);


%%
% different stains
for i = 1 : length(data_2D)
%     close all;
    if i == data_2D_register
        continue
    end
%     align_2D_intermodality(data_2D_directory, data_2D{data_2D_register}{2}, data_2D{i}{2}, output_prefix);
    OPT = struct;
    if i == 2
        OPT.method = 'test2'; % this one wasn't working great for myelin, particularly bad in the cerebellum
    elseif i == 3
        OPT.method = 'mi';
    end
    align_2D_intermodality_test(data_2D_directory, data_2D{data_2D_register}{2}, data_2D{i}{2}, output_prefix, OPT);
end

%%
% find centers
close all;
find_centers_for_initialization_nissl(data_2D_directory, data_2D{data_2D_register}{2}, output_prefix);



%%
% slice to neighbor
% there seems to be an issue with downsampling here
% where if I downsample too much I get errors
close all;
r = 25;
r = 100;
r = 10;
r = 3;
% this downsampling and iterations is enough for a good initial guess
% not enough for a full accurate reconstruction
downs = [16,8,4];
niter = [40,20,10];
% niter = [40,20,10]*2;
downs = [16,8,4,2];
niter = [40,20,10,5]*2;

downs = [4,2];


skip_thick = -1; % no thick slices to be skipped
load_initializer = 1;
e = 0.05;
atlas_free_rigid_alignment_v02(data_2D_directory, data_2D{data_2D_register}{2}, output_prefix, r, downs, niter, e, skip_thick, load_initializer)


%%
qcreport = {};
if stop_for_qc
    % set up qc pipeline here

    % check all intermodality
    for i = 1 : length(data_2D)
        if i == data_2D_register
            continue
        end
        
        pattern0 = data_2D{data_2D_register}{2};        
        pattern1 = data_2D{i}{2};
        input_dir = data_2D_directory;

                
        
        % load the file info
        geometry_file = dir([data_2D_directory '*.csv']);
        fid = fopen([input_dir geometry_file(1).name],'rt');
        line = fgetl(fid); % ignore the first line
        csv_data = {};
        is_0 = [];
        is_1 = [];
        count = 0;
        while 1
            line = fgetl(fid);
            if line == -1
                break
            end

            count = count + 1;

            % check if it matches the pattern
            if (regexp(line,regexptranslate('wildcard',pattern0)))
                is_0(count) = 1;
                is_1(count) = 0;
                is_other(count) = 0;
            elseif (regexp(line,regexptranslate('wildcard',pattern1)))
                is_1(count) = 1;
                is_0(count) = 0;
                is_other(count) = 0;
            else
                is_other(count) = 1;
                is_0(count) = 0;
                is_1(count) = 0;
            end


            % process this line, splitting at commas
            csv_data(count,:) = strsplit(line,',');
            %     
        end
        fclose(fid);
        files = csv_data(:,1);
        nxJ0 = cellfun(@(x)str2num(x), csv_data(:,2:3));
        dxJ0 = cellfun(@(x)str2num(x), csv_data(:,5:6));
        zJ0 = cellfun(@(x)str2num(x), csv_data(:,10));

        
        
        
        

        for j = 1 : length(zJ0)
            % load the saved transforms, every time in case I update
            vars = load([output_prefix 'zero_to_one_' strrep(pattern0,'*','x') '_to_' strrep(pattern1,'*','x') '.mat']);
            if ~vars.is_1(j)
                continue
            end
            
            % find the nearest pattern 0
            cost = (zJ0(j) - zJ0).^2 - is_0(:)*1e10;
            ind = find(  cost == min(cost) ,1,'first');
            inds(j) = ind;
            % load the images
            I = imread([input_dir files{ind}]);
            if isa(I,'uint8')
                I = double(I)/255.0;    
            elseif isa(I,'uint16')
                I = double(I);
                for c = 1 : size(I,3)
                    I(:,:,c) = tiedrank(I(:,:,c));
                end
                I = I / max(I(:));
            end
            dxI = dxJ0(ind,:);
            xI = (0:size(I,2)-1)*dxI(1);xI = xI - mean(xI);
            yI = (0:size(I,1)-1)*dxI(2);yI = yI - mean(yI);

            J = imread([input_dir files{j}]);
            if isa(J,'uint8')
                J = double(J)/255.0;    
            elseif isa(J,'uint16')

        %         J = double(J);
        %         for c = 1 : size(I,3)
        %             J(:,:,c) = tiedrank(J(:,:,c));
        %         end
        %         J = J / max(J(:));

                J = double(J);
                J_ = J(:);
                J_ = tiedrank(J_);
                J = reshape(J_,size(J));
                J = J / max(J(:));

            end

            dxJ = dxJ0(j,:);
            xJ = (0:size(J,2)-1)*dxJ(1);xJ = xJ - mean(xJ);
            yJ = (0:size(J,1)-1)*dxJ(2);yJ = yJ - mean(yJ);
            
             
            
            Aqc = vars.zero_to_one(:,:,j);
            % transform I
            Ai = inv(Aqc);
            [XJ,YJ] = meshgrid(xJ,yJ);
            Xs = Ai(1,1)*XJ + Ai(1,2)*YJ + Ai(1,3);
            Ys = Ai(2,1)*XJ + Ai(2,2)*YJ + Ai(2,3);
            AI = zeros(size(J));
            for c = 1 : 3
                F = griddedInterpolant({yI,xI},I(:,:,c));
                AI(:,:,c) = F(Ys,Xs);
            end
            danfigure(88);
            imagesc(xI,yI,I)
            axis image
            title I;
            danfigure(89);
            imagesc(xJ,yJ,AI)
            axis image
            title AI;
            danfigure(90);
            imagesc(xJ,yJ,J)
            axis image
            title J;            
            % now we need some kind of error display
            tmp = tiedrank(AI(:));
            tmp = tmp / max(tmp);
            tmp = reshape(tmp,size(AI));            
            X = reshape(tmp,[],3);
            Y = reshape(J,[],3);
            X = [X,ones(size(X(:,1)))];
            c = (X'*X)\(X'*Y);
            danfigure(91);
            imagesc(xJ,yJ,reshape(X*c,size(J)))   
            title fAI
            axis image
            
            danfigure(92);
            imagesc(xJ,yJ, (reshape(X*c,size(J)) - J)*0.5 + 0.5 )
            axis image
            title error
            
            drawnow;
            while 1 % this is the interactive qc loop
                
                disp(['**********'])
                disp([input_dir files{j}])
                disp(['file ' num2str(j) ' of ' num2str(length(zJ0))])
                
                
                yesno = input(['Is this registration poor quality? (0=keep as is, 1=redo) >>']);
                qcreport{i}{j} = struct('initial_screen_good',yesno);
                if isempty(yesno) || yesno==0 
                    save([output_prefix 'intermodality_qcreport.mat'],'qcreport');
                    break
                end
                % otherwise
                R = manual_landmark_align_rigid_for_qc(xI,yI,I,xJ,yJ,J);
                qcreport{i}{j}.initial_alignment = Aqc;
                qcreport{i}{j}.new_alignment = R;
                vars.zero_to_one(:,:,j) = R;
                Ai = inv(R);
                Xs = Ai(1,1)*XJ + Ai(1,2)*YJ + Ai(1,3);
                Ys = Ai(2,1)*XJ + Ai(2,2)*YJ + Ai(2,3);
                AI = zeros(size(J));
                for c = 1 : 3
                    F = griddedInterpolant({yI,xI},I(:,:,c),'linear','nearest');
                    AI(:,:,c) = F(Ys,Xs);
                end
                tmp = tiedrank(AI(:));
                tmp = tmp / max(tmp);
                tmp = reshape(tmp,size(AI));       
                X = reshape(tmp,[],3);
                Y = reshape(J,[],3);
                X = [X,ones(size(X(:,1)))];
                c = (X'*X)\(X'*Y);
                
                
                danfigure(89)
                imagesc(xJ,yJ,AI)
                axis image
                title AI;
            
            
                danfigure(91);
                imagesc(xJ,yJ,reshape(X*c,size(J)))                
                title('fAI')
                axis image

                danfigure(92);
                imagesc(xJ,yJ, (reshape(X*c,size(J)) - J)*0.5 + 0.5 )
                title('error')
                axis image
                drawnow;
                
                yesno = input(['Do you want to initialize a rigid transform or are you happy with this (0=no keep as is, 1=yes run rigid mapping)']);
                if isempty(yesno) || yesno==0 
                    % save
                    save([output_prefix 'zero_to_one_' strrep(pattern0,'*','x') '_to_' strrep(pattern1,'*','x') '.mat'],'-struct','vars');
                    save([output_prefix 'intermodality_qcreport.mat'],'qcreport');
                    break
                end
                % run the registration
                downs = [1];
                niter = 1000;
    %             bins = linspace(0,1,10);
    %             sigma = 0.1;
    %             e = 5e2;
    %             R1 = slice_to_slice_rigid_alignment_MI(xI,yI,I,xJ,yJ,J,R,downs,niter,e,bins,sigma);
                e = 0.5;
                OPT = struct;
                OPT.sigmaM = 0.15;
                OPT.sigmaA = 0.5;
                OPT.draw = 0;
                OPT.draw = 1;
                R1 = slice_to_slice_rigid_alignment_GN_weight(xI,yI,I,xJ,yJ,J,R,downs,niter,e,OPT);
                qcreport{i}{j}.new_alignment = R1;
                vars.zero_to_one(:,:,j) = R1;
                
                Ai = inv(R1);
                Xs = Ai(1,1)*XJ + Ai(1,2)*YJ + Ai(1,3);
                Ys = Ai(2,1)*XJ + Ai(2,2)*YJ + Ai(2,3);
                AI = zeros(size(J));
                for c = 1 : 3
                    F = griddedInterpolant({yI,xI},I(:,:,c),'linear','nearest');
                    AI(:,:,c) = F(Ys,Xs);
                end
                tmp = tiedrank(AI(:));
                tmp = tmp / max(tmp);
                tmp = reshape(tmp,size(AI));       
                X = reshape(tmp,[],3);
                Y = reshape(J,[],3);
                X = [X,ones(size(X(:,1)))];
                c = (X'*X)\(X'*Y);
                danfigure(89)
                imagesc(xJ,yJ,AI)
                axis image
                title AI;
            
            
                danfigure(91);
                imagesc(xJ,yJ,reshape(X*c,size(J)))                
                title('fAI')
                axis image

                danfigure(92);
                imagesc(xJ,yJ, (reshape(X*c,size(J)) - J)*0.5 + 0.5 )
                title('error')
                axis image
                drawnow;
                
                
                
                yesno = input(['Is the quality still bad (0=no Im happy, 1=yes redo it)']);
                
                
                if isempty(yesno) || yesno==0 
                    save([output_prefix 'zero_to_one_' strrep(pattern0,'*','x') '_to_' strrep(pattern1,'*','x') '.mat'],'-struct','vars');
                    save([output_prefix 'intermodality_qcreport.mat'],'qcreport');
                    break
                end
                
            end

        end % of loop over slices

         
    end % of loop over data types

    

    
end
%%
close all;
%% qc check within modality
% the plan here is 
% 1. check if it is an outlier
% 2. decide which slice to align to
% 3. do the alignment
% 4. repeat if necessary
% for outlier checking, it is much better to sample on a uniform and consitent grid
% than to have the grid changing
% especially for the first few slices

qcreport = {};
if stop_for_qc

    % check all within modality, this is generally nissl to nissl
    for i = 1 : length(data_2D)
        if i ~= data_2D_register % we negate the condition from previous section
            continue
        end
        
        pattern0 = data_2D{data_2D_register}{2};        
        input_dir = data_2D_directory;

        % load the file info for nissl only (i.e. pattern 0)
        geometry_file = dir([data_2D_directory '*.csv']);
        fid = fopen([input_dir geometry_file(1).name],'rt');
        line = fgetl(fid); % ignore the first line
        csv_data = {};
        count = 0;
        while 1
            line = fgetl(fid);
            if line == -1
                break
            end
            % check if it matches the pattern
            if isempty(regexp(line,regexptranslate('wildcard',pattern0)))
                continue
            end
            

            % process this line, splitting at commas
            count = count + 1;
            csv_data(count,:) = strsplit(line,',');
            %     
        end
        fclose(fid);
        files = csv_data(:,1);
        nxJ0 = cellfun(@(x)str2num(x), csv_data(:,2:3));
        dxJ0 = cellfun(@(x)str2num(x), csv_data(:,5:6));
        % how about x0?
        zJ0 = cellfun(@(x)str2num(x), csv_data(:,10));
        
        % find a sampling grid
        % I'm having an issue where there are some outliers leading to a
        % big grid which crashes memory
        vars = load([output_prefix 'initializer_A.mat']);
        xmin = inf;
        xmax = -inf;
        ymin = inf;
        ymax = -inf;
        corners = [];
        for j = 1 : length(zJ0)
            xJ = (0:nxJ0(j,1)-1)*dxJ0(j,1);xJ = xJ - mean(xJ);
            yJ = (0:nxJ0(j,2)-1)*dxJ0(j,2);yJ = yJ - mean(yJ);
            lastwarn('');
            Ai = inv(vars.AJ(:,:,j));
            % we need to check the lastwarn to see if it is singular
            if contains(lastwarn,'singular')
                continue
            end
            % on the other hand, what if it is an outlier, but not singular
            % possibly 354 which has 1e5
            % also 300 is 1e6
            % 293 also, 273
%             if j == 300 || j == 18 || j == 3 % these are big ones with amplitude over 1e6
%                 continue
%             end
            
            corners_ = Ai(1:2,1:2)*[
                    xJ(1),yJ(1);
                    xJ(1),yJ(end);
                    xJ(end),yJ(1);
                    xJ(end),yJ(end);
                    ]' + Ai(1:2,end);
%             thisxmin = min(corners_(1,:));
%             thisxmax = max(corners_(1,:));
%             thisymin = min(corners_(2,:));
%             thisymax = max(corners_(2,:));
%             corners = [corners;[thisxmin,thisxmax,thisymin,thisymax]];
            xmin = min(min(corners_(1,:)),xmin);
            xmax = max(max(corners_(1,:)),xmax);
            ymin = min(min(corners_(2,:)),ymin);
            ymax = max(max(corners_(2,:)),ymax);
        end
%         xmin = quantile(corners(:,1),0.2);
%         xmax = quantile(corners(:,2),0.8);
%         ymin = quantile(corners(:,3),0.2);
%         ymax = quantile(corners(:,4),0.8);
%         xmin = -1e4;xmax=1e4;
%         ymin = -1e4;ymax=1e4;
        xg = xmin : dxJ0(1,1) : xmax;
        yg = ymin : dxJ0(1,2) : ymax;
        [XG,YG] = meshgrid(xg,yg);

        

        for j = 1 : length(zJ0)
            while 1 % this is the interactive qc loop
                
                disp(['**********'])
                disp([input_dir files{j}])
                disp(['file ' num2str(j) ' of ' num2str(length(zJ0))])
                
                
                
                % reload each time in case there are updates
                vars = load([output_prefix 'initializer_A.mat']);

                % now let's load this slice I
                I = imread([input_dir files{j}]);
                if isa(I,'uint8')
                    I = double(I)/255.0;    
                elseif isa(I,'uint16')
                    I = double(I);
                    for c = 1 : size(I,3)
                        I(:,:,c) = tiedrank(I(:,:,c));
                    end
                    I = I / max(I(:));
                end
                dxI = dxJ0(j,:);
                xI = (0:size(I,2)-1)*dxI(1);xI = xI - mean(xI);
                yI = (0:size(I,1)-1)*dxI(2);yI = yI - mean(yI);
                [XI,YI] = meshgrid(xI,yI);
                h = danfigure(1001);
                clf;
                pos = get(h,'position');
                % set width to three times the height
                pos(3) = pos(4)*4;
                set(h,'position',pos);
%                 subplot(1,5,3)
%                 imagesc(xI,yI,I)
%                 axis image
%                 title('This slice (3)')

                thisA = vars.AJ(:,:,j);


                % load 4 other slices, and transform to match this slice
                rshow = 3;
                for k = -rshow : rshow

                    try
                        J = imread([input_dir files{j+k}]);
                    catch
                        continue
                    end                    
                    if isa(J,'uint8')
                        J = double(J)/255.0;    
                    elseif isa(J,'uint16')
                        J = double(J);
                        for c = 1 : size(J,3)
                            J(:,:,c) = tiedrank(J(:,:,c));
                        end
                        J = J / max(J(:));
                    end
                    dxJ = dxJ0(j+k,:);
                    xJ = (0:size(J,2)-1)*dxJ(1);xJ = xJ - mean(xJ);
                    yJ = (0:size(J,1)-1)*dxJ(2);yJ = yJ - mean(yJ);
                    [XJ,YJ] = meshgrid(xJ,yJ);

                    

                    Ai= vars.AJ(:,:,j+k);
                    
                    Xs = Ai(1,1)*XG + Ai(1,2)*YG + Ai(1,3);
                    Ys = Ai(2,1)*XG + Ai(2,2)*YG + Ai(2,3);
                    J_ = zeros(size(XG,1),size(XG,2),size(J,3));
                    for c = 1 : size(J,3)
                        F = griddedInterpolant({yJ,xJ},J(:,:,c),'linear','nearest');
                        J_(:,:,c) = F(Ys,Xs);
                    end

                    danfigure(h);
                    subplot(1,rshow*2+1,k+rshow+1)
                    imagesc(xg,yg,J_)
                    axis image;
                    title(num2str(k+rshow+1))
                    if k == 0
                        title('This slice')
                    end

                end


                yesno = input(['Is this slice an outlier? (0=keep as is, 1=redo) >>']);
                qcreport{i}{j} = struct('initial_screen_good',yesno);
                if isempty(yesno) || yesno==0 
                    save([output_prefix 'intramodality_qcreport.mat'],'qcreport');
                    break
                end
                
                number = input(['Which slice do you want to align to? (number) >>']);                    
                ind = j + (number-rshow-1);
                

                % load the image
                J = imread([input_dir files{ind}]);
                if isa(J,'uint8')
                    J = double(J)/255.0;    
                elseif isa(J,'uint16')
                    J = double(J);
                    for c = 1 : size(J,3)
                        J(:,:,c) = tiedrank(J(:,:,c));
                    end
                    J = J / max(J(:));
                end
                dxJ = dxJ0(ind,:);
                xJ = (0:size(J,2)-1)*dxJ(1);xJ = xJ - mean(xJ);
                yJ = (0:size(J,1)-1)*dxJ(2);yJ = yJ - mean(yJ);
                
                
                % let's find a common space for resampling                                
                A_ = vars.AJ(:,:,ind);
                Ai = inv(A_);
                % computed corners on fixed grid above
%                 % find the action on the corners of J
%                 corners = [
%                     xJ(1),yJ(1);
%                     xJ(1),yJ(end);
%                     xJ(end),yJ(1);
%                     xJ(end),yJ(end);
%                     ]';
%                 newcorners = Ai(1:2,1:2)*corners + Ai(1:2,end);
%                 
%                 
%                 xg = min(newcorners(1,:)) : dxJ(1) : max(newcorners(1,:));
%                 yg = min(newcorners(2,:)) : dxJ(2) : max(newcorners(2,:));
%                 [XG,YG] = meshgrid(xg,yg);
                                
                Xs = A_(1,1)*XG + A_(1,2)*YG + A_(1,3);
                Ys = A_(2,1)*XG + A_(2,2)*YG + A_(2,3);
                J_ = zeros(size(XG,1),size(XG,2),size(J,3));
                for c = 1 : size(J,3)
                    F = griddedInterpolant({yJ,xJ},J(:,:,c),'linear','nearest');
                    J_(:,:,c) = F(Ys,Xs);
                end

                



                % now I have to get the transform matrices
                % this is more complicated because both images have a transform
                % that I need to worry about
                % I should combine them so I is in its native space or
                % something
                Aqc = vars.AJ(:,:,j);
                % transform I
                [XJ,YJ] = meshgrid(xJ,yJ);
                Xs = Aqc(1,1)*XG + Aqc(1,2)*YG + Aqc(1,3);
                Ys = Aqc(2,1)*XG + Aqc(2,2)*YG + Aqc(2,3);
                AI = zeros(size(XG,1),size(XG,2),size(J,3));
                for c = 1 : 3
                    F = griddedInterpolant({yI,xI},I(:,:,c),'linear','nearest');
                    AI(:,:,c) = F(Ys,Xs);
                end
                danfigure(88);
                imagesc(xI,yI,I)
                axis image
                title I;
                danfigure(89);
                imagesc(xg,yg,AI)
                axis image
                title AI;
                danfigure(90);
                imagesc(xg,yg,J_)
                axis image
                title J;
                % now we need some kind of error display
                danfigure(92);
                imagesc(xg,yg, (AI - J_)*0.5 + 0.5 );
                axis image
                title error
                drawnow;
                
                % when I'm doing qc
                % keep in mind that A is from common to the images
                % I will calculate image I to common though
                % so I will need to invert it
                yesno = input(['Is this registration poor quality? (0=keep as is, 1=redo) >>']);
                qcreport{i}{j} = struct('initial_screen_good',yesno);
                if isempty(yesno) || yesno==0 
                    save([output_prefix 'intramodality_qcreport.mat'],'qcreport');
                    break
                end
                % otherwise
                R = manual_landmark_align_rigid_for_qc(xI,yI,I,xg,yg,J_);
                % R is from image to common, but I want to save the inverse
                
                qcreport{i}{j}.initial_alignment = Aqc;
                qcreport{i}{j}.new_alignment = inv(R);
                vars.AJ(:,:,j) = inv(R);
                      
                Ai = inv(R);
                Xs = Ai(1,1)*XG + Ai(1,2)*YG + Ai(1,3);
                Ys = Ai(2,1)*XG + Ai(2,2)*YG + Ai(2,3);
                AI = zeros(size(J_));
                for c = 1 : 3
                    F = griddedInterpolant({yI,xI},I(:,:,c),'linear','nearest');
                    AI(:,:,c) = F(Ys,Xs);
                end
                
                
                danfigure(89)
                imagesc(xg,yg,AI)
                axis image
                title AI;
            
 
                danfigure(92);
                imagesc(xg,yg, (AI- J_)*0.5 + 0.5 )
                title('error')
                axis image
                drawnow;
                
                yesno = input(['Do you want to initialize a rigid transform or are you happy with this (0=no keep as is, 1=yes run rigid mapping)']);
                if isempty(yesno) || yesno==0 
                    % save
                    save([output_prefix 'initializer_A.mat'],'-struct','vars');
                    save([output_prefix 'intramodality_qcreport.mat'],'qcreport');
                    break
                end
                % run the registration
                downs = [1];
                niter = 1000;
                niter = 500; % less because within modality is easier
    %             bins = linspace(0,1,10);
    %             sigma = 0.1;
    %             e = 5e2;
    %             R1 = slice_to_slice_rigid_alignment_MI(xI,yI,I,xJ,yJ,J,R,downs,niter,e,bins,sigma);
                e = 0.5;
                OPT = struct;
                OPT.sigmaM = 0.15;
                OPT.sigmaA = 0.5;
                OPT.draw = 0;
                OPT.draw = 1;
                R1 = slice_to_slice_rigid_alignment_GN_weight(xI,yI,I,xg,yg,J_,R,downs,niter,e,OPT);
                qcreport{i}{j}.new_alignment = inv(R1);
                vars.AJ(:,:,j) = inv(R1);
                
                Ai = inv(R1);
                Xs = Ai(1,1)*XG + Ai(1,2)*YG + Ai(1,3);
                Ys = Ai(2,1)*XG + Ai(2,2)*YG + Ai(2,3);
                AI = zeros(size(J_));
                for c = 1 : 3
                    F = griddedInterpolant({yI,xI},I(:,:,c),'linear','nearest');
                    AI(:,:,c) = F(Ys,Xs);
                end
                danfigure(89)
                imagesc(xg,yg,AI)
                axis image
                title AI;
            

                danfigure(92);
                imagesc(xg,yg, (AI - J_)*0.5 + 0.5 )
                title('error')
                axis image
                drawnow;
                
                
                
                yesno = input(['Is the quality still bad (0=no Im happy, 1=yes redo it)']);
                
                
                if isempty(yesno) || yesno==0 
                    save([output_prefix 'initializer_A.mat'],'-struct','vars');
                    save([output_prefix 'intermodality_qcreport.mat'],'qcreport');
                    break
                end
                
            end

        end % of loop over slices

         
    end % of loop over data types

    

    
end
%%
% initial affine
close all;
downs = [8,4];
niter = 30;
% get atlas file

files = dir([out_dir 'common_space_3D_data/' num2str(mapping_3D_2D(1),'%03d') '_' num2str(mapping_3D_2D(2),'%03d') '*.vtk']);
atlas_file = [out_dir 'common_space_3D_data/' files(1).name];
affine_for_initial_alignment(atlas_file, data_2D_directory, data_2D{data_2D_register}{2}, output_prefix, downs, niter,A0)


%%
% 2D to 3D
% now 3D to 2D transformations for nissl
close all;
warning('off','MATLAB:griddedInterpolant:MeshgridEval2DWarnId')
warning('off','MATLAB:griddedInterpolant:MeshgridEval3DWarnId')

% here we will run registration at either only low res
% or only low and medium res
ThreeD_to_2D_registration(atlas_file, data_2D_directory, data_2D{data_2D_register}{2}, config_file, output_prefix)



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we are finished with registration
% the remainder is applying transforms and creating standard outputs
% we will do 2D data, then 3D data
% nothing below should reference varaibles from above (other than inputs 
% on the first cell)

%%
% now we begin creating the standard outputs
% all data in all spaces
% let's start with the 2D data
close all;

% load the transform
disp('Starting to load saved data')
Avars = load([output_prefix 'down1_A.mat']);
vvars = load([output_prefix 'down1_v.mat']);
disp('Finished loading saved data')
A = Avars.A;
AJ = Avars.AJ;
dxVJ = zeros(length(vvars.xJ),2);
for f = 1 : length(vvars.xJ)
    try
        dxVJ(f,:) = [vvars.xJ{f}(2)-vvars.xJ{f}(1),vvars.yJ{f}(2)-vvars.yJ{f}(1)];
    catch
        keyboard
        error('mismatch in number of files')
    end
end

vtx = vvars.vtx;
vty = vvars.vty;
vtz = vvars.vtz;

xV = vvars.xI;
yV = vvars.yI;
zV = vvars.zI;
[XV,YV,ZV] = meshgrid(xV,yV,zV);

nt = size(vtx,4);
dt = 1.0/nt;
vJtx = vvars.vJtx;
vJty = vvars.vJty;

% create 3D transform
disp('Starting to integrate velocity field')
phiinvx = XV;
phiinvy = YV;
phiinvz = ZV;
for t = 1 : nt
    disp(['Integrating v field step '  num2str(t) ' of ' num2str(nt)])
    % update phi
    Xs = XV - vtx(:,:,:,t)*dt;
    Ys = YV - vty(:,:,:,t)*dt;
    Zs = ZV - vtz(:,:,:,t)*dt;
    % subtract and add identity
    F = griddedInterpolant({yV,xV,zV},phiinvx-XV,'linear','nearest');
    phiinvx = F(Ys,Xs,Zs) + Xs;
    F = griddedInterpolant({yV,xV,zV},phiinvy-YV,'linear','nearest');
    phiinvy = F(Ys,Xs,Zs) + Ys;
    F = griddedInterpolant({yV,xV,zV},phiinvz-ZV,'linear','nearest');
    phiinvz = F(Ys,Xs,Zs) + Zs;
end

phix = XV;
phiy = YV;
phiz = ZV;
disp('Starting to integrate velocity field backward')
for t = nt : -1 : 1
    disp(['Integrating v field step '  num2str(t) ' of ' num2str(nt)])
    % update phi
    Xs = XV + vtx(:,:,:,t)*dt;
    Ys = YV + vty(:,:,:,t)*dt;
    Zs = ZV + vtz(:,:,:,t)*dt;
    % subtract and add identity
    F = griddedInterpolant({yV,xV,zV},phix-XV,'linear','nearest');
    phix = F(Ys,Xs,Zs) + Xs;
    F = griddedInterpolant({yV,xV,zV},phiy-YV,'linear','nearest');
    phiy = F(Ys,Xs,Zs) + Ys;
    F = griddedInterpolant({yV,xV,zV},phiz-ZV,'linear','nearest');
    phiz = F(Ys,Xs,Zs) + Zs;
end
Aphix = A(1,1)*phix + A(1,2)*phiy + A(1,3)*phiz + A(1,4);
Aphiy = A(2,1)*phix + A(2,2)*phiy + A(2,3)*phiz + A(2,4);
Aphiz = A(3,1)*phix + A(3,2)*phiy + A(3,3)*phiz + A(3,4);

%%
% now we will loop over all the slices
disp('Starting to load slice geometry info')
geometry_file = dir([data_2D_directory '*.csv']);
copyfile([data_2D_directory geometry_file(1).name],output_prefix)

fid = fopen([data_2D_directory geometry_file(1).name],'rt');
line = fgetl(fid); % ignore the first line
% it should say
% filename, nx, ny, nz, dx, dy, dz, x0, y0, z0

csv_data = {};
count = 0;
while 1
    line = fgetl(fid);
    if line == -1
        break
    end
    count = count + 1;
    % process this line, splitting at commas
    csv_data(count,:) = strsplit(line,',');
    %     
end
fclose(fid);
slice_filenames = csv_data(:,1);
nxJ = cellfun(@(x)str2num(x), csv_data(:,2:3));
dxJ = cellfun(@(x)str2num(x), csv_data(:,5:6));

x0J = cellfun(@(x)str2num(x), csv_data(:,8:9));
z0J = cellfun(@(x)str2num(x), csv_data(:,10));

zJ = z0J;
dzJ = cellfun(@(x) str2num(x), csv_data(:,7));

xJ = cell(0);
yJ = cell(0);
for f = 1 : length(slice_filenames)
    xJ{f} = x0J(f,1) + (0:nxJ(f,1)-1)*dxJ(f,1);
    yJ{f} = x0J(f,2) + (0:nxJ(f,2)-1)*dxJ(f,2);    
end
disp('finished loading slice geometry info')
%%
% combine the 2D rigid alignments so that they match the geometry file
disp('combining 2D transforms from different modalities')
close all;

AJAll = zeros(3,3,length(slice_filenames));
vvarsAll = vvars;
count = 0;
for f = 1 : length(slice_filenames)
    % check the pattern
    for i = 1 : length(data_2D)
        pattern = data_2D{i}{2};
        if ~isempty(regexp(slice_filenames{f},regexptranslate('wildcard',pattern)))
            break;            
        end
    end
    % now we have an index stored in i
    thisname = data_2D{i}{1};
    thispattern = data_2D{i}{2};
    if i == data_2D_register
        count = count + 1;
        AJAll(:,:,f) = AJ(:,:,count);
        vvarsAll.vJtx{f} = vvars.vJtx{count};
        vvarsAll.vJty{f} = vvars.vJty{count};
        vvarsAll.xJ{f} = vvars.xJ{count};
        vvarsAll.yJ{f} = vvars.yJ{count};
        continue
    end
    
    
    % otherwise load data
    vars = load([output_prefix 'zero_to_one_' strrep(data_2D{data_2D_register}{2},'*','x')  '_to_' strrep(thispattern,'*','x') '.mat']);
    ind = vars.inds(f); 
    % this stores the nearest register slice in file order
    % we still need to know its count
    count_ = 0;
    for f_ = 1 : ind
        if ~isempty(regexp(slice_filenames{f_},regexptranslate('wildcard',data_2D{data_2D_register}{2})))
            count_ = count_ + 1;
        end
    end

    
    AJAll(:,:,f) = AJ(:,:,count_) * vars.zero_to_one(:,:,f);
    vvarsAll.vJtx{f} = vvars.vJtx{count_};
    vvarsAll.vJty{f} = vvars.vJty{count_};
    vvarsAll.xJ{f} = vvars.xJ{count_};
    vvarsAll.yJ{f} = vvars.yJ{count_};


    
    
end

% combine and save as previous
AvarsAll = Avars;
AvarsAll.AJ = AJAll;
save([output_prefix 'combined_A.mat'],'-struct','AvarsAll')
save([output_prefix 'combined_v.mat'],'-struct','vvarsAll')
% and we need the v with the right index as well
disp('finished 2D transforms from different modalities')

%%
% get sampling grid
disp('Constructing uniform 2D sampling grid')
themin = inf;
themax = -inf;
for d = 1 : length(data_3D_files)    
    for di = 1 : length(data_3D_files{d})
        thisname = data_3D_files{d}{di}{1};
        % load the image
        common_dir = [out_dir 'common_space_3D_data/'];
        files = dir( [common_dir num2str(d,'%03d') '_' num2str(di,'%03') '*.vtk']);
        fname = [common_dir files(1).name];
        [xI,yI,zI,I,titleI,namesI] = read_vtk_image(fname);
        themin = min([xI,yI,zI,themin]);
        themax = max([xI,yI,zI,themax]);
        
    end
end
xg = themin:dxJ(1,1):themax;
yg = themin:dxJ(1,2):themax;
[Xg,Yg] = meshgrid(xg,yg);
disp('Finished constructing uniform 2D sampling grid')

%%
% in this example we bring each 3D dataset 
% from the COMMON SPACE to the RECONSTRUCTED SPACE
% note that from the common space to the reconstructed space is only one
% transform
% for other spaces there could be other transforms
% NOTE, currently this is getting data from the resampled space which is
% not ideal
% better would be to get it from the native space, composing transforms,
% that is a to do
% march 6
% I need to get atlas overlay and recolor it
% I will need to chose d=1,di=3,f=458 (for slice 160)

close all;


for native = [0] % if native = 1, don't apply rigid
    if native == 0
        output_dir = [output_prefix '2D_' data_2D{data_2D_register}{1} '_space/'];
    else
        output_dir = [output_prefix '2D_' data_2D{data_2D_register}{1} '_native_space/'];
    end
    mkdir(output_dir)

%     for d = 1%1 : length(data_3D_files)
    for d = 1 : length(data_3D_files)
        disp(['Space ' num2str(d)])

%         for di = 3%1 : length(data_3D_files{d})
        for di = 1 : length(data_3D_files{d})
            thisname = data_3D_files{d}{di}{1};
            % load the image
            common_dir = [out_dir 'common_space_3D_data/'];
            files = dir( [common_dir num2str(d,'%03d') '_' num2str(di,'%03d') '*.vtk']);
            fname = [common_dir files(1).name];

            [xI,yI,zI,I,titleI,namesI] = read_vtk_image(fname);

            danfigure(999);
            if contains(thisname,'seg')
                sliceView(xI,yI,zI,mod(I,32))
            else
                sliceView(xI,yI,zI,I)
            end
            title(thisname)


            % loop through the slices that were used for registration
            count = 0;
%             for f = 458%1 : length(slice_filenames)
            for f = 1 : length(slice_filenames)
                % I only need to load the 2D images the first time to build
                % them in reconstructed space
                
                first_time = (d == 1 && di == 1) ;
                if first_time
                    
                    
                    % load it
                    J = imread([data_2D_directory slice_filenames{f}]);
                    if strcmp(class(J) ,'uint8')
                        J = double(J)/255.0;

                    elseif strcmp(class(J),'uint16')  
                        J = double(J);
        %                 for c = 1 : size(J,3)
        %                     J(:,:,c) = tiedrank(J(:,:,c));
        %                 end
        %                 J = J / max(J(:));     
                        for c = 1 : size(J,3)
                            J(:,:,c) = J(:,:,c) - min(min(J(:,:,c)));
                            J(:,:,c) = J(:,:,c) / max(max(J(:,:,c)));
                        end
                    end


                    danfigure(2);
                    imagesc(xJ{f},yJ{f},J);
                    axis image

                    % apply its rigid motion and determine a sampling grid
                    % the sampling grid should be large to allow padding
                    A_ = AJAll(:,:,f);
                    if native == 1
                        A_ = eye(3);
                    end
                    Xs = A_(1,1)*Xg + A_(1,2)*Yg + A_(1,3);
                    Ys = A_(2,1)*Xg + A_(2,2)*Yg + A_(2,3);
                    AiJ = zeros(size(Xs,1),size(Xs,2),size(J,3));
                    for c = 1 : size(J,3)
                        F = griddedInterpolant({yJ{f},xJ{f}},J(:,:,c),'linear','none');
                        AiJ(:,:,c) = F(Ys,Xs);
                    end

                    danfigure(3);
                    imagesc(xg,yg,AiJ)
                    axis image
                    
                    % write out AiJ
                    [a,b,c] = fileparts(slice_filenames{f});
                    write_vtk_image(xg,yg,[zJ(f),zJ(f)+dzJ(f)],AiJ,[ output_dir b '.vtk'],thisname);
                    
                    % and write out the transform
                    Zg = zeros(size(Xg)) + zJ(f);
                    Zs = Zg;
                    write_vtk_image(xg,yg,[zJ(f),zJ(f)+dzJ(f)],cat(4,Xs-Xg,Ys-Yg,Zs-Zg),[output_dir b  '_to_' 'input' '_displacement.vtk'],thisname);
                    
                    
                    % write out AJ
                    fid = fopen([output_dir b  '_to_' 'input' '_matrix.txt'],'wt');
                    fprintf(fid,'%f %f %f\n',AJAll(:,:,f)');                    
                    fclose(fid);

                end
                % now transform I to match
                % I need Ai and phii
                Zg = zeros(size(Xg)) + zJ(f);
                Ai = inv(A);
                if native == 1
                    tmp = inv(AJAll(:,:,f));
                    Ai = eye(4);
                    Ai(1:2,1:2) = tmp(1:2,1:2);
                    Ai(1:2,end) = tmp(1:2,end);
                    Ai = inv(A) * Ai;
                end
                Xs = Ai(1,1)*Xg + Ai(1,2)*Yg + Ai(1,3)*Zg + Ai(1,4);
                Ys = Ai(2,1)*Xg + Ai(2,2)*Yg + Ai(2,3)*Zg + Ai(2,4);
                Zs = Ai(3,1)*Xg + Ai(3,2)*Yg + Ai(3,3)*Zg + Ai(3,4);
                % now use phii
                F = griddedInterpolant({yV,xV,zV},phiinvx-XV,'linear','nearest');
                Xs_ = F(Ys,Xs,Zs) + Xs;
                F = griddedInterpolant({yV,xV,zV},phiinvy-YV,'linear','nearest');
                Ys_ = F(Ys,Xs,Zs) + Ys;
                F = griddedInterpolant({yV,xV,zV},phiinvz-ZV,'linear','nearest');
                Zs_ = F(Ys,Xs,Zs) + Zs;
                % now get the transformed image
                if contains(thisname,'seg') 
                    F = griddedInterpolant({yI,xI,zI},I,'nearest','none');
                else
                    F = griddedInterpolant({yI,xI,zI},I,'linear','none');
                end
                I_ = F(Ys_,Xs_,Zs_);
                danfigure(4);
                if contains(thisname,'seg') 
                    imagesc(xg,yg,mod(I_,32))
                    if max(I_) > 2^16 - 1
                        I_ = uint32(I_);
                    else
                        I_ = uint16(I_);
                    end
                else
                    imagesc(xg,yg,I_)
                end
                axis image

%                 % for bingxing
%                 imagesc(xg,yg,double(mod(I_,32))/32.0 + AiJ)
%                 axis image

                % write out I and J as vtk files
                [a,b,c] = fileparts(slice_filenames{f});
                write_vtk_image(xg,yg,[zJ(f),zJ(f)+dzJ(f)],I_,[output_dir thisname '_to_' b '.vtk'],thisname);
                
                
%                 [x,y,z,I,title_,names] = read_vtk_image([output_dir thisname '_to_' b '.vtk']); % test
                
                % if this is the first time we can write out the
                % displacement fields to the COMMON space     
                if first_time
                    write_vtk_image(xg,yg,[zJ(f),zJ(f)+dzJ(f)],cat(4,Xs_-Xg,Ys_-Yg,Zs_-Zg),[output_dir b  '_to_' data_3D_spaces{data_3D_common} '_displacement.vtk'],thisname);
                    % to do, write out coordinate grids as geogson
                    % write out mean coordinates
                    fid = fopen([output_dir b  '_to_' data_3D_spaces{data_3D_common} '_mean_xyz.txt'],'wt');
                    fprintf(fid,'%f\n',[mean(Xs_(:)),mean(Ys_(:)),mean(Zs_(:))]);
                    fclose(fid);
                    

                    % write out grids (common space coordinates)
                    % calculating curves in this section is quite slow
                    
                    tic
                    levels = min([xI,yI,zI]):dx_contour:max([xI,yI,zI]); 
                    
                    contoursx = {};
                    contoursy = {};
                    contoursz = {};
                    for l = 1 : length(levels)
                    contoursx{l} = contourc(xg,yg,Xs_,levels(l)*[1,1]);
                    contoursy{l} = contourc(xg,yg,Ys_,levels(l)*[1,1]);
                    contoursz{l} = contourc(xg,yg,Zs_,levels(l)*[1,1]);
                    end
                    
                    
                   
%                     danfigure(321654);
%                     imagesc(xg,yg,AiJ)
%                     hold on;
%                     contour(xg,yg,Xs_,levels,'r');
%                     contour(xg,yg,Ys_,levels,'g');
%                     contour(xg,yg,Zs_,levels,'b');
%                     hold off
%                     axis image;
%                     title(['z = ' num2str(zJ(f))])
%                     axis image;
%                     xlabel x;
%                     ylabel y;
%                     legend('common x','common y', 'common z')
                    

                    
                    % do x contours
                    contoursx_ = {};            
                    namesx_ = {};
                    contoursy_ = {};            
                    namesy_ = {};
                    contoursz_ = {};            
                    namesz_ = {};
                    for l = 1 : length(levels)
                        namesx_{l} = num2str(levels(l));
                        sc = contoursx{l};   
                        count = 0;
                        while 1
                            if isempty(sc)
                                break
                            end
                            count = count + 1;
                            n = sc(2,1);

                            xdata = sc(1,2:n+1);
                            ydata = sc(2,2:n+1);

                            contoursx_{l}{count} = [xdata;ydata];                    
                            sc = sc(:,n+2:end);

                        end
                    end
                    empty = cellfun(@(x) isempty(x),contoursx_);
                    contoursx_ = contoursx_(~empty);
                    namesx_ = namesx_(~empty);

                    for l = 1 : length(levels)
                        namesy_{l} = num2str(levels(l));
                        sc = contoursy{l};                
                        count = 0;
                        while 1
                            if isempty(sc)
                                break
                            end
                            count = count + 1;
                            n = sc(2,1);

                            xdata = sc(1,2:n+1);
                            ydata = sc(2,2:n+1);

                            contoursy_{l}{count} = [xdata;ydata];                    
                            sc = sc(:,n+2:end);

                        end
                    end
                    empty = cellfun(@(x) isempty(x),contoursy_);
                    contoursy_ = contoursy_(~empty);
                    namesy_ = namesy_(~empty);
                    for l = 1 : length(levels)
                        namesz_{l} = num2str(levels(l));
                        sc = contoursz{l};  
                        count = 0;
                        while 1
                            if isempty(sc)
                                break
                            end
                            count = count + 1;
                            n = sc(2,1);

                            xdata = sc(1,2:n+1);
                            ydata = sc(2,2:n+1);

                            contoursz_{l}{count} = [xdata;ydata];                    
                            sc = sc(:,n+2:end);

                        end
                    end
                    empty = cellfun(@(x) isempty(x),contoursz_);
                    contoursz_ = contoursz_(~empty);
                    namesz_ = namesz_(~empty);

                    % write out the json
                    contournamex = [output_dir b '_to_' data_3D_spaces{data_3D_common}  '_grid_x.geojson'];
                    contournamey = [output_dir b '_to_' data_3D_spaces{data_3D_common}  '_grid_y.geojson'];
                    contournamez = [output_dir b '_to_' data_3D_spaces{data_3D_common}  '_grid_z.geojson'];

                    write_geojson_grid(contournamex,namesx_,contoursx_);
                    write_geojson_grid(contournamey,namesy_,contoursy_);
                    write_geojson_grid(contournamez,namesz_,contoursz_);
                    
                    toc
                end

                if contains(thisname,'seg')
                    % note this is a double interpolation which is especially
                    % bad for labels
                    % I really should compose the transformation
                    % that will be a TODO
                    % ALSO need a QC figure TODO
                    
                    % load target J for qc figures
                    J = imread([data_2D_directory slice_filenames{f}]);
                    if strcmp(class(J) ,'uint8')
                        J = double(J)/255.0;
                    elseif strcmp(class(J),'uint16')  
                        J = double(J);
                        for c = 1 : size(J,3)
                            J(:,:,c) = J(:,:,c) - min(min(J(:,:,c)));
%                             J(:,:,c) = J(:,:,c) / max(max(J(:,:,c)));
                            tmp = J(:,:,c);
                            J(:,:,c) = J(:,:,c)/ quantile(tmp(:),0.99);
                            
                        end
                    end


                   % apply its rigid motion and determine a sampling grid
                    % the sampling grid should be large to allow padding
                    A_ = AJAll(:,:,f);
                    if native == 1
                        A_ = eye(3);
                    end
                    Xs = A_(1,1)*Xg + A_(1,2)*Yg + A_(1,3);
                    Ys = A_(2,1)*Xg + A_(2,2)*Yg + A_(2,3);
                    AiJ = zeros(size(Xs,1),size(Xs,2),size(J,3));
                    for c = 1 : size(J,3)
                        F = griddedInterpolant({yJ{f},xJ{f}},J(:,:,c),'linear','none');
                        AiJ(:,:,c) = F(Ys,Xs);
                    end

                    danfigure(3);
                    imagesc(xg,yg,AiJ)
                    axis image
                    
                    labels = unique(I_(:));
                    structure_contours = cell(1,length(labels));
                    for l = 1 : length(labels)
                        % this part is fairly slow
                        structure_contours{l} = contourc(xg,yg,double(I_==labels(l)),[0.5,0.5]); % draw single level contours
                    end


                    hold on;
                    % plot them
                    ids = [];
                    contours = {};
                    for l = 1 : length(labels)
                        ids(l) = labels(l);
                        sc = structure_contours{l};
                        count = 0;
                        while 1
                            if isempty(sc)
                                break
                            end
                            count = count + 1;
                            n = sc(2,1);
                            xdata = sc(1,2:n+1);
                            ydata = sc(2,2:n+1);


                            s = 1 : (length(xdata)-1);
                            s = s - mean(s);
                            bb = exp(-s.^2/2.0/10^2); % with scale1 = 5
                            bb = exp(-s.^2/2.0/5^2); % with scale1 = 10
                            bb = exp(-s.^2/2.0/1^2);
                            % I just want to blur it
                            bb = bb / sum(bb);
                            bhat = fft(bb);
                            xdatanew_ = ifft( (fft(xdata(1:end-1)) .* bhat) )';
                            ydatanew_ = ifft( (fft(ydata(1:end-1)) .* bhat)  )';
    %                         xdatanew = xdatanew_;
    %                         ydatanew = ydatanew_;

                            % or forget smoothing
                            xdatanew = xdata(1:end-1)';
                            ydatanew = ydata(1:end-1)';

                            sc = sc(:,n+2:end);




                            % if it is too small, do not add it
                            if length(xdatanew) < 20
                                count = count - 1;
                                continue
                            end
                            danfigure(3);
                            plot([xdatanew; xdatanew(1)],[ydatanew; ydatanew(1)],'k','linewidth',2)
                            plot([xdatanew; xdatanew(1)],[ydatanew; ydatanew(1)],'w','linewidth',1)
%                             plot([xdatanew; xdatanew(1)],[ydatanew; ydatanew(1)],'k','linewidth',1.5) % grant figure
                            children = get(gca,'children');
                            % first time through the children are line line
                            % image                            
                            children = [children(1);children(3:end);children(2);children(end)];
                            set(gca,'children',children)
    %                         drawnow


                            contours{l}{count} = [[xdatanew; xdatanew(1)],[ydatanew; ydatanew(1)]]';


                        end

                    end
                    hold off;
                    
                    
                    
                    % write out figure
                    contourname = [output_dir thisname '_to_' b '_QC.png'];
                    saveas(3,contourname)
%                     return

                    % write out
                    contourname = [output_dir thisname '_to_' b '.geojson'];

                    % note it is possible that ids will be longer than contours
                    % if the contours are empty
                    write_geojson_ontology(contourname,ids(1:length(contours)),contours,0);
                    
                    


                    

                end
                
                
                    
                    
                drawnow
                
                


            end
            
            


        end
    end
end


%%
% last is the 3D data
% remember we want all data in all spaces
% this includes the 2D summary data
% the common space we should already have, but it was native resolution
% we want same resolution for voxel to voxel comparisons
%
% I have saved maps from common to each space
% I need to use the forward and inverse composed
% 

close all;

for d0 = 1:length(data_3D_files)
    disp(['**********'])    
    disp(['Putting data in space ' num2str(d0) ': ' data_3D_spaces{d0}])    
    % load the first image from this space and display it
    [xI,yI,zI,I,titleI,namesI] = read_vtk_image(data_3D_files{d0}{1}{2});
    danfigure(1);
    sliceView(xI,yI,zI,I)
    title(['Image 1 from ' data_3D_spaces{d0} ' space'])
    [XI,YI,ZI] = meshgrid(xI,yI,zI);
    
    % find the map from common to this space
    % is this the common space?
    if d0 == data_3D_common
        % if so load identity
        vars0 = struct;
        vars0.Aphix = XI;
        vars0.Aphiy = YI;
        vars0.Aphiz = ZI;
        vars0.phiiAix = XI;
        vars0.phiiAiy = YI;
        vars0.phiiAiz = ZI;
        vars0.xI = xI;
        vars0.yI = yI;
        vars0.zI = zI;

        vars0.xJ = xI; % there is no xJ
        vars0.yJ = yI;
        vars0.zJ = zI;
        
    else
        % otherwise load the data
        m = find(mapping_3D_pairs(:,1) == d0);
        pairs = mapping_3D_pairs(m,:);

        atlas_name = data_3D_files{data_3D_common}{pairs(3)}{1};
        target_name = data_3D_files{pairs(1)}{pairs(2)}{1};

        prefix_outdir = [out_dir 'mapping_3D_' num2str(m,'%03d') '_' atlas_name '_to_' target_name '/'];
        vars0 = load([prefix_outdir 'output.mat']);
    end
    

    % create an output dir
    prefix_outdir_ = [out_dir  data_3D_spaces{d0} '_space/'];
    if ~exist(prefix_outdir_,'dir')
        mkdir(prefix_outdir_);
    end
    
    
    % now we loop over all the other spaces
    for d1 = 1 : length(data_3D_files)
        if d0 == d1
            continue
        end
        disp(['Collecting data from space ' num2str(d1) ': ' data_3D_spaces{d1}])
        
        % load the first image from this space and display it
        [xJ,yJ,zJ,J,titleJ,namesJ] = read_vtk_image(data_3D_files{d1}{1}{2});
        danfigure(2);
        sliceView(xJ,yJ,zJ,J)
        [XJ,YJ,ZJ] = meshgrid(xJ,yJ,zJ);

        % to chose what data to load, we exampine mapping_3D_pairs and data_3D_common
        % if it is the common space, we just load identity
        if d1 == data_3D_common
             % if so load identity
            vars1 = struct;
            vars1.Aphix = XJ;
            vars1.Aphiy = YJ;
            vars1.Aphiz = ZJ;
            vars1.phiiAix = XJ;
            vars1.phiiAiy = YJ;
            vars1.phiiAiz = ZJ;
            vars1.xJ = xJ;
            vars1.yJ = yJ;
            vars1.zJ = zJ;
            vars1.xI = xJ;
            vars1.yI = yJ;
            vars1.zI = zJ;
        else
            % otherwise load the data
            m = find(mapping_3D_pairs(:,1) == d1);
            pairs = mapping_3D_pairs(m,:);

            atlas_name = data_3D_files{data_3D_common}{pairs(3)}{1};
            target_name = data_3D_files{pairs(1)}{pairs(2)}{1};

            prefix_outdir = [out_dir 'mapping_3D_' num2str(m,'%03d') '_' atlas_name '_to_' target_name '/'];
            vars1 = load([prefix_outdir 'output.mat']);
        end

        % we want to take the coordinates in I
        [XI,YI,ZI] = meshgrid(xI,yI,zI);
        % then apply phiiAi0 to them (map the coordinates in I to the
        % coordinates in common)
        F = griddedInterpolant({vars0.yJ,vars0.xJ,vars0.zJ},vars0.phiiAix,'linear','nearest');
        Tx = F(YI,XI,ZI);
        F = griddedInterpolant({vars0.yJ,vars0.xJ,vars0.zJ},vars0.phiiAiy,'linear','nearest');
        Ty = F(YI,XI,ZI);
        F = griddedInterpolant({vars0.yJ,vars0.xJ,vars0.zJ},vars0.phiiAiz,'linear','nearest');
        Tz = F(YI,XI,ZI);

        % then apply Aphi1 to them (map the coordinates in common to the
        % coordinates in J)
        F = griddedInterpolant({vars1.yI,vars1.xI,vars1.zI},vars1.Aphix,'linear','nearest');
        TTx = F(Ty,Tx,Tz);
        F = griddedInterpolant({vars1.yI,vars1.xI,vars1.zI},vars1.Aphiy,'linear','nearest');
        TTy = F(Ty,Tx,Tz);
        F = griddedInterpolant({vars1.yI,vars1.xI,vars1.zI},vars1.Aphiz,'linear','nearest');
        TTz = F(Ty,Tx,Tz);
        
        % we write out the transformation
        outname = [prefix_outdir_  data_3D_spaces{d0} '_to_' data_3D_spaces{d1}  '_displacement.vtk'];
        write_vtk_image(xI,yI,zI,cat(4,TTx-XI,TTy-YI,TTz-ZI),outname,[data_3D_spaces{d0} '_to_' data_3D_spaces{d1}  '_displacement'])
        

        
            
        for di = 1 : length(data_3D_files{d1})
            thisname = data_3D_files{d1}{di}{1};
            disp(thisname)
            
             % load the first image from this space and display it
            [xJ,yJ,zJ,J,titleJ,namesJ] = read_vtk_image(data_3D_files{d1}{di}{2});

            
            
            % now transform the image            
            if contains(thisname,'seg') 
                F = griddedInterpolant({yJ,xJ,zJ},J,'nearest','none');
                TJ = F(TTy,TTx,TTz);
                if max(J(:)) > 2^16-1
                    TJ = uint32(TJ);
                else
                    TJ = uint16(TJ);
                end
                TJs = double(mod(TJ,7))/7;
            else
                F = griddedInterpolant({yJ,xJ,zJ},J,'linear','none');
                TJ = F(TTy,TTx,TTz);
                TJs = TJ;
            end
            danfigure(33);
            sliceView(xI,yI,zI,TJs)
            danfigure(34);
            sliceView(xI,yI,zI,cat(4,I,TJs,I))
            
            
            % write out the image
            outname = [prefix_outdir_  thisname '_to_' data_3D_spaces{d0}  '.vtk'];
            write_vtk_image(xI,yI,zI,TJ,outname,[thisname '_to_' data_3D_spaces{d0}])

            
            saveas(33,[outname(1:end-4) '.png'])
            saveas(34,[outname(1:end-4) '_error.png'])
    
    
            
        end % end of loop over images
    end % end of loop over "from" spaces
    
    
    
    % last, we will reconstruct the slices and map its summary into the 3D
    % spaces    

    
    % we have to start with the 3D coordinates
    % XI
    % then we go from 3D to common, so we apply vars0.phiiAix
    % then we go from common to recone, so we apply Aphix
    F = griddedInterpolant({vars0.yJ,vars0.xJ,vars0.zJ},vars0.phiiAix,'linear','nearest');
    Tx = F(YI,XI,ZI);
    F = griddedInterpolant({vars0.yJ,vars0.xJ,vars0.zJ},vars0.phiiAiy,'linear','nearest');
    Ty = F(YI,XI,ZI);
    F = griddedInterpolant({vars0.yJ,vars0.xJ,vars0.zJ},vars0.phiiAiz,'linear','nearest');
    Tz = F(YI,XI,ZI);
    F = griddedInterpolant({yV,xV,zV},Aphix,'linear','nearest');
    TTx = F(Ty,Tx,Tz);
    F = griddedInterpolant({yV,xV,zV},Aphiy,'linear','nearest');
    TTy = F(Ty,Tx,Tz);
    F = griddedInterpolant({yV,xV,zV},Aphiz,'linear','nearest');
    TTz = F(Ty,Tx,Tz);

    
    
    for d1 = 1 : length(data_2D)
        thisname = data_2D{d1}{1};
        disp(['Collecting 2D data ' thisname])
        % we will do this as follows
        % 1. construct a 3D volume the size of I
        % 2. load slices and downsample to about the resolution
        % 3. Add them and keep track of increments
        J3D = zeros(size(I,1),size(I,2),size(I,3),3); % rgb
        O3D = zeros(size(I,1),size(I,2),size(I,3)); % accumulator
        % for transformations we need
        %
        
        % loop through the slices in the geometry file
        count = 0;
        for f = 1 : length(slice_filenames)
            thisfile = slice_filenames{f};
            % check pattern
            pattern = data_2D{d1}{2};
            if isempty(regexp(thisfile,regexptranslate('wildcard',pattern)))
                continue 
            end
            
            J = imread([data_2D_directory thisfile]);
            J = double(J);
            xJ = (0:nxJ(f,1)-1)*dxJ(f,1) + x0J(f,1);
            yJ = (0:nxJ(f,2)-1)*dxJ(f,2) + x0J(f,2);
            danfigure(44);
            imagesc(xJ,yJ,J/quantile(J(:),0.9))
            axis image
            
            xJ_ = [xJ(1)-dxJ(f,1),xJ,xJ(end)+dxJ(f,1)];
            yJ_ = [yJ(1)-dxJ(f,2),yJ,yJ(end)+dxJ(f,2)];
            zJ_ = z0J(f) + [-1,0,1]*dzJ(f)*4;
            [XJ,YJ,ZJ] = meshgrid(xJ_,yJ_,zJ_);

            A_ = AJAll(:,:,f);
            A_ = [A_(1:2,1:2),[0;0],A_(1:2,end);0,0,1,0;0,0,0,1];
            TTTx = A_(1,1)*TTx + A_(1,2)*TTy + A_(1,3)*TTz + A_(1,4);
            TTTy = A_(2,1)*TTx + A_(2,2)*TTy + A_(2,3)*TTz + A_(2,4);
            TTTz = A_(3,1)*TTx + A_(3,2)*TTy + A_(3,3)*TTz + A_(3,4);
                
            % since this is so slow, can I cut out a bounding box?
            ind = (TTTx >= min(xJ_)) .* (TTTx <= max(xJ_)) .* (TTTy >= min(yJ_)) .* (TTTy <= max(yJ_)) .* (TTTz >= min(zJ_)) .* (TTTz <= max(zJ_));
            indz = sum(sum(ind,1),2);
            indz0 = find(indz,1,'first');
            indz1 = find(indz,1,'last');
            indx = sum(sum(ind,1),3);
            indx0 = find(indx,1,'first');
            indx1 = find(indx,1,'last');            
            indy = sum(sum(ind,2),3);
            indy0 = find(indy,1,'first');
            indy1 = find(indy,1,'last');
            
            for c = 1 : size(J,3)
                J_ = J(:,:,c);
                O_ = J_*0+1;
                % make this a 3D image
                J_ = padarray(J_,[1,1,1],0,'both');
                O_ = padarray(O_,[1,1,1],0,'both');
                
                
                
                F = griddedInterpolant({yJ_,xJ_,zJ_},J_,'linear','nearest');
                J3D(indy0:indy1,indx0:indx1,indz0:indz1,c) = J3D(indy0:indy1,indx0:indx1,indz0:indz1,c) + F(TTTy(indy0:indy1,indx0:indx1,indz0:indz1),TTTx(indy0:indy1,indx0:indx1,indz0:indz1),TTTz(indy0:indy1,indx0:indx1,indz0:indz1));
                if c == 1
                    F = griddedInterpolant({yJ_,xJ_,zJ_},O_,'linear','nearest');
                    O3D(indy0:indy1,indx0:indx1,indz0:indz1) = O3D(indy0:indy1,indx0:indx1,indz0:indz1) + F(TTTy(indy0:indy1,indx0:indx1,indz0:indz1),TTTx(indy0:indy1,indx0:indx1,indz0:indz1),TTTz(indy0:indy1,indx0:indx1,indz0:indz1));                
                end
                
                
            end
            
            if ~mod(count,10)
                danfigure(556);
    %             sliceView(xI,yI,zI,J3D./O3D,5,[min(J(:)),max(J(:))])
                sliceView(xI,yI,zI,J3D./O3D,5,quantile(J(:),[0.1,0.99]))
                drawnow
            end
            
            count = count + 1;             
            
        end % end of slices for this stain
        J3D_ = J3D./O3D;
        danfigure(556);
        sliceView(xI,yI,zI,J3D_,5,quantile(J3D_(:),[0.1,0.99]))
        drawnow
        

        % now I should write out the image
        outname = [prefix_outdir_  thisname '_to_' data_3D_spaces{d0}  '.vtk'];
        write_vtk_image(xI,yI,zI,J3D_,outname,[thisname '_to_' data_3D_spaces{d0}])
        saveas(556,[outname(1:end-4) '.png'])
        
        
    end % end of this stain
    
    
    
end % end of loop over "to" spaces