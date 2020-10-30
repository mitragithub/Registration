function ThreeD_to_2D_registration(template_name, target_dir, pattern, config_file, output_dir, nonrigid_thick_only)
% map 3D atlas to 2D slices
% only use files that match pattern (with wildcards)
% In addition to normal registration, we want to have a second pass of
% registration
% here we want to incorporate some manually painted labels
% to do this we will need to have labelled slices as an input
% and also we'll need to know which atlas label they correspond to
% let's do this in a simple manner with filenames
% a label image should be a binary valued tiff
% it should have the same filename
% but it should say L for label instead of N for nissl
% after, we add what labels it should correspond to
% separated by underscores, so hypens can stil be used to separate parts
% e.g.
% PMD1028&1027-N50-2012.04.19-18.40.56_PMD1027_1_0148.tif
% becomes
% PMD1028&1027-L50_0_1_2-2012.04.19-18.40.56_PMD1027_1_0148.tif
% for a binary label that corresponds to CCF 0,1,2
% we can have another one
% PMD1028&1027-L50_5-2012.04.19-18.40.56_PMD1027_1_0148.tif
% for a binary label that corresponds to CCF 5 only
% we will call this "edit mode"
% and enable it by inputing template name as a cell array (image, labels)


tic


if nargin < 6
    % for rnaseq,data, we want thick slices, >20um, to be nonrigidly
    % registered.  Otherwise slices are transformed rigidly only
    nonrigid_thick_only = 0;
end
thick_cutoff = 20;
if contains(template_name,'marmoset')
    thick_cutoff = 80;
end

addpath Functions/plotting
addpath Functions/downsample
addpath Functions/vtk
addpath Functions/frame2Gif
addpath Functions/gradient


% we will generally use an initializer for affine transformations
initializer = [output_dir 'initializer_A.mat'];
if ~exist(initializer,'file')
    initializer = '';
end

% edit mode
% if template is a cell
edit_mode = 0;
if iscell(template_name) && length(template_name) ==2
    edit_mode = 1;
    label_name = template_name{2};
    template_name = template_name{1};
    disp('Edit mode enabled')
end




%%
% make directory for outputs
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

%%
% get config options
config = config_reader(config_file);

%%
% first thing is now to get slice thicknesses and location
geometry_file = dir([target_dir '*.csv']);
fid = fopen([target_dir geometry_file(1).name],'rt');
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
    % check if it matches the pattern
    if isempty(regexp(line,regexptranslate('wildcard',pattern)))
        continue
    end
    
    % check pattern
    count = count + 1;
    % process this line, splitting at commas
    csv_data(count,:) = strsplit(line,',');
    %
end
fclose(fid);
files = csv_data(:,1);
nxJ0 = cellfun(@(x)str2num(x), csv_data(:,2:3));
dxJ0 = cellfun(@(x)str2num(x), csv_data(:,5:6));
zJ0 = cellfun(@(x)str2num(x), csv_data(:,10));
isthicks = cellfun(@(x) str2num(x), csv_data(:,7)) > thick_cutoff;
% for now don't worry about offset, I'll just calculate it to be the center


center = zeros(length(files),2);% this will be for centering my priors, if there is an initializer, it will get updated
%%
% for each slice
% we need a cell array

edit_info = cell(length(files),1);
if edit_mode
    
    % if you are in edit mode, get corresponding label images
    % note there may be more than one label per image
    edit_files = dir([target_dir '*-L*_*.tif']);
    label_strs = {};
    for f = 1 : length(edit_files)
        % find the ids and the labels
        fname = edit_files(f).name;
        parts = strsplit(fname,'-');
        % use last four digits to get number
        num = fname(end-7:end-4);
        % get the slice id
        tmp = strfind(files,['_' num '.tif']);
        slice_id = find(cellfun(@(x) ~isempty(x),tmp),1,'first');
        % what labels
        ls = strsplit(parts{2},'_');
        ls = sort(cellfun(@(x)str2num(x),ls(2:end)));
        
        
        edit_info{slice_id}(end+1) = struct('labels',ls,'filename',fname);
        
        label_strs{end+1} = num2str(ls);
        
    end
    
    
    
    % find what labels I'll need to extract in terms of binary images
    tmp = unique(label_strs);
    tmp = tmp(cellfun(@(x) ~isempty(x),tmp));
    label_sets = {};
    nL = length(tmp);
    for l = 1 : nL
        label_sets{l} = sscanf(tmp{l}, '%d');
    end
    % loop back through and assign a channel to each label
    for f = 1 : length(edit_info)
        for ff = 1 : length(edit_info{f})
            for fff = 1 : length(label_sets)
                if length(edit_info{f}(ff).labels) == length(label_sets{fff}) && all(edit_info{f}(ff).labels == label_sets{fff}')
                    edit_info{f}(ff).channel = fff;
                end
            end
%             
        end
    end
end

%%
for downloop = 1 : 3
%%    
    
    if downloop == 1
        % now we're not doing cropping and not doing downsampling beforehand
        downI = str2double(config.DOWNLOOP1.downI);
        downJ = str2double(config.DOWNLOOP1.downJ);
        blurI = str2double(config.DOWNLOOP1.blurI);
        blurJ = str2double(config.DOWNLOOP1.blurJ);
    elseif downloop == 2
        % now down by 2
        downI = str2double(config.DOWNLOOP2.downI);
        downJ = str2double(config.DOWNLOOP2.downJ);
        blurI = str2double(config.DOWNLOOP2.blurI);
        blurJ = str2double(config.DOWNLOOP2.blurJ);

    elseif downloop == 3
        % now down by 1
        downI = str2double(config.DOWNLOOP3.downI);
        downJ = str2double(config.DOWNLOOP3.downJ);
        blurI = str2double(config.DOWNLOOP3.blurI);
        blurJ = str2double(config.DOWNLOOP3.blurJ);
    end
    
    
    
    
    
    
    [xI,yI,zI,I,title_,names] = read_vtk_image(template_name);
    I = double(I);
    % scale it for numerical stability, since its scale doesn't matter
    Imean = mean(I(:));
    Istd = std(I(:));
    
    I = I - mean(I(:));
    I = I/std(I(:));
    
    if edit_mode
        % load labels and construct the one hot image for each set
        [~,~,~,L_,~,~] = read_vtk_image(label_name);
        L = zeros(size(L_,1),size(L_,2),size(L_,3),nL);
        for l = 1 : nL
            for ll = 1 : length(label_sets{l})
                L(:,:,:,l) = L(:,:,:,l) | (L_==label_sets{l}(ll));
            end
        end
    end
    
    if blurI>0
        x_ = [-round(3*blurI):round(3*blurI)];
        [X_,Y_,Z_] = meshgrid(x_);
        K_ = exp(-(X_.^2 + Y_.^2 + Z_.^2)/2/blurI^2);
        K_ = K_/sum(K_(:));
        n_ = size(K_);
        
        Ipad = padarray(I,[n_(1)-1,n_(2)-1,n_(3)-1]/2,'both');
        O = padarray(ones(size(I)),[n_(1)-1,n_(2)-1,n_(3)-1]/2,'both');
        
        K_pad = padarray(K_,size(Ipad)-size(K_),0,'post');
        K_pad = circshift(K_pad,-[1,1,1]*round(3*blurI));
        
        % really need to blur with better boundary conditions!
        % this would probably be to pad with reflections
        Ob = ifftn(fftn(O).*fftn(K_pad),'symmetric');
        tmp = ifftn(fftn(Ipad).*fftn(K_pad),'symmetric')./Ob;
        %     I = ifftn(fftn(I).*fftn(K_pad),'symmetric');
        I = tmp((n_(1)-1)/2+1:end-(n_(1)-1)/2, (n_(2)-1)/2+1:end-(n_(2)-1)/2, (n_(3)-1)/2+1:end-(n_(3)-1)/2);
        
        if edit_mode
            Lpad = padarray(L,[n_(1)-1,n_(2)-1,n_(3)-1,0]/2,'both');
            tmp = bsxfun(@times, fft(fft(fft(Lpad,[],1),[],2),[],3), K_pad);
            tmp = bsxfun(@rdivide,ifft(ifft(ifft(tmp,[],1),[],2),[],3,'symmetric'),Ob);
            L = tmp((n_(1)-1)/2+1:end-(n_(1)-1)/2, (n_(2)-1)/2+1:end-(n_(2)-1)/2, (n_(3)-1)/2+1:end-(n_(3)-1)/2,:);
        end
    end
    I = padarray(I,[1,1,1]*downI,-Imean/Istd,'both'); % zero pad (actually with - mean)
    [~,~,~,I] = downsample(1:size(I,2),1:size(I,1),1:size(I,3),I,downI*[1,1,1]);
    if edit_mode
        L_ = padarray(L,[1,1,1,0]*downI,0,'both');
        for l = 1 : nL
            [~,~,~,tmp] = downsample(1:size(L_,2),1:size(L_,1),1:size(L_,3),L_(:,:,:,l),downI*[1,1,1]);
            if l == 1
                L = zeros(size(tmp,1),size(tmp,2),size(tmp,3),nL);
            end
            L(:,:,:,l) = tmp;
        end
    end
    
    
    
    % note this currently assumes zero centered
    dxI = [xI(2)-xI(1), yI(2)-yI(1), zI(2)-zI(1)]*downI;
    nxI = [size(I,2),size(I,1),size(I,3)];
    xI = (0:nxI(1)-1)*dxI(1);
    yI = (0:nxI(2)-1)*dxI(2);
    zI = (0:nxI(3)-1)*dxI(3);
    xI = xI - mean(xI);
    yI = yI - mean(yI);
    zI = zI - mean(zI);
    
    danfigure(1);
    sliceView(xI,yI,zI,I);
    [XI,YI,ZI] = meshgrid(xI,yI,zI);
    fxI = (0:nxI(1)-1)/nxI(1)/dxI(1);
    fyI = (0:nxI(2)-1)/nxI(2)/dxI(2);
    fzI = (0:nxI(3)-1)/nxI(3)/dxI(3);
    [FXI,FYI,FZI] = meshgrid(fxI,fyI,fzI);
    
    
    
    
    
    %%
    % load all the files
    danfigure(2);
    dxJ = dxJ0*downJ;
    zJ = zJ0;
    
    for f = 1 : length(files)
        
        J{f} = double(imread([target_dir files{f}]))/255.0;
        % note this normalization may not be appropriate for different
        % datasets
        if contains(pattern,'DK39')
            % maybe I should detect padding as black pixels?
            J{f} = normcdf(bsxfun(@minus,J{f}, mean(mean(J{f},1),2))/std(J{f}(:)));
        end
        
        % first thing is to get a mask on pixels that are pure white
%         danfigure(2729);
%         WMask{f} = 1.0 - double(all(J{f}==1,3));
%         subplot(1,2,1)
%         imagesc(WMask{f})
%         axis image

        % this step may not be appropriate if we're not dealing with nissl
        % it will detect white backgrounds that span an entire row
        WMask{f} = 1.0 - detect_padding_in_nissl(J{f});
%         subplot(1,2,2);
%         imagesc(WMask{f})
%         axis image


        if edit_mode
            % load some label files
            LJ{f} = [];
            nLJ(f) = 0;
            if ~isempty(edit_info{f})
                for ff = 1 : length(edit_info{f})
                    tmp = double(imread([target_dir edit_info{f}(ff).filename]));
                    LJ{f}(:,:,ff) = tmp;
                end
                nLJ(f) = size(LJ{f},3);
            end
        end
        
       
        
        % blur it
        if blurJ>0
            x_ = [-round(3*blurJ):round(3*blurJ)];
            [X_,Y_] = meshgrid(x_);
            K_ = exp(-(X_.^2 + Y_.^2 )/2/blurJ^2);
            K_ = K_/sum(K_(:));
            n_ = size(K_);
            Jpad = padarray(J{f},[n_(1)-1,n_(2)-1,0]/2,'both');
            Maskpad = padarray(WMask{f},[n_(1)-1,n_(2)-1,0]/2,'both');
            O = padarray(ones(size(J{f},1),size(J{f},2)),[n_(1)-1,n_(2)-1,0]/2,'both');
            
            % need to pad to size N + M - 1
            K_pad = padarray(K_,[size(Jpad,1),size(Jpad,2)]-size(K_),0,'post');
            K_pad = circshift(K_pad,-[1,1,1]*round(3*blurJ));
            K_pad_hat = fftn(K_pad);
            Ob = ifftn(fftn(O).*K_pad_hat,'symmetric');
            for c = 1 : size(J{f},3)
                % we really should pad the data before bluring it
                tmp = ifftn(fftn(Jpad(:,:,c)).*K_pad_hat,'symmetric')./Ob;
                J{f}(:,:,c) = tmp((n_(1)-1)/2 + 1 : end - (n_(1)-1)/2,(n_(2)-1)/2 + 1 : end - (n_(2)-1)/2);
            end
            tmp = ifftn(fftn(Maskpad).*K_pad_hat,'symmetric')./Ob;
            WMask{f} = tmp((n_(1)-1)/2 + 1 : end - (n_(1)-1)/2,(n_(2)-1)/2 + 1 : end - (n_(2)-1)/2);
            if edit_mode && ~isempty(edit_info{f})
                LJ_pad = padarray(LJ{f},[n_(1)-1,n_(2)-1,0]/2,'both');
                for l = 1 : size(LJ{f},3)
                    tmp = ifftn(fftn(LJ_pad(:,:,l)).*K_pad_hat,'symmetric')./Ob;
                    LJ{f}(:,:,l) = tmp((n_(1)-1)/2 + 1 : end - (n_(1)-1)/2,(n_(2)-1)/2 + 1 : end - (n_(2)-1)/2);
                end
            end
        end
        
        tmp = [];
        for c = 1 : size(J{f},3)
            [~,~,tmp(:,:,c)] = downsample2D(1:size(J{f},2),1:size(J{f},1),J{f}(:,:,c),downJ*[1,1]);
        end
        J{f} = tmp;
        [~,~,WMask{f}] = downsample2D(1:size(WMask{f},2),1:size(WMask{f},1),WMask{f},downJ*[1,1]);
        if edit_mode && ~isempty(edit_info{f})
            tmp = [];
            for l = 1 : size(LJ{f},3)
                [~,~,tmp(:,:,l)] = downsample2D(1:size(LJ{f},2),1:size(LJ{f},1),LJ{f}(:,:,l),downJ*[1,1]);
            end
            LJ{f} = tmp;
        end
        
        
        nxJ{f} = [size(J{f},2),size(J{f},1)];
        xJ{f} = (0 : nxJ{f}(1)-1)*dxJ(1);
        yJ{f} = (0 : nxJ{f}(2)-1)*dxJ(2);
        xJ{f} = xJ{f} - mean(xJ{f});
        yJ{f} = yJ{f} - mean(yJ{f});
        [XJ{f},YJ{f}] = meshgrid(xJ{f},yJ{f});
        AJ(:,:,f) = eye(3); % affine matrix
        
        if ~mod(f-1,10)
            danfigure(2);
            imagesc(xJ{f},yJ{f},J{f})
            axis image
            drawnow
        end
    end
    
    
    
    
    %%
    % settings
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if downloop >= 1
        start_3d_diffeo = str2double(config.DOWNLOOP1.start_3d_diffeo);
        start_2d_diffeo = str2double(config.DOWNLOOP1.start_2d_diffeo);
        start_2d_affine = str2double(config.DOWNLOOP1.start_2d_affine);
        niter = str2double(config.DOWNLOOP1.niter);
        
        if ~isempty(initializer)
            Afilename = initializer;
        else
            Afilename = '';
        end
        vfilename = '';
        prefix = [output_dir 'down4_'];
        
        sigmaR = str2double(config.DEFAULT.sigmaR);
        sigmaRJ = str2double(config.DEFAULT.sigmaRJ);
        sigmaM = str2double(config.DEFAULT.sigmaM);
        sigmaB = str2double(config.DEFAULT.sigmaB);
        sigmaA = str2double(config.DEFAULT.sigmaA);
        
        
        piM_width = str2double(config.DEFAULT.piM_width);
        piM = str2double(config.DEFAULT.piM);
        piB = str2double(config.DEFAULT.piB);
        piA = str2double(config.DEFAULT.piA);
        
        
        
        % step sizes
        eVJ = str2double(config.DOWNLOOP1.eVJ);
        eVI = str2double(config.DOWNLOOP1.eVI);
        eTI = str2double(config.DOWNLOOP1.eTI);
        eLI = str2double(config.DOWNLOOP1.eLI);
        eLJ = str2double(config.DOWNLOOP1.eLJ);
        eTJ = str2double(config.DOWNLOOP1.eTJ);
        post_affine_reduce = str2double(config.DOWNLOOP1.post_affine_reduce);
        
        
        % bounds
        VImax_fac = str2double(config.DEFAULT.VImax_fac);
        VImax = VImax_fac*dxI(1);
        VJmax_fac = str2double(config.DEFAULT.VJmax_fac);
        VJmax = VJmax_fac*dxJ(1);
        LJmax = str2double(config.DEFAULT.LJmax);
        LImax = str2double(config.DEFAULT.LImax);
        TJmax_fac = str2double(config.DEFAULT.TJmax_fac);
        TJmax = dxJ(1)*TJmax_fac;
        TImax_fac = str2double(config.DEFAULT.TImax_fac);
        TImax = dxI(1)*TImax_fac;
        
        
        
        % smoothness
        p = str2double(config.DEFAULT.p);
        a = str2double(config.DEFAULT.a);
        pJ = str2double(config.DEFAULT.pJ);
        aJ = str2double(config.DEFAULT.aJ);
        ap_fac = str2double(config.DEFAULT.ap_fac);
        ap = dxI(1)*ap_fac;
        aJp_fac = str2double(config.DEFAULT.aJp_fac);
        aJp = dxJ(1)*aJp_fac;

        
        % flow
        nt = str2double(config.DEFAULT.nt);
        ntJ = str2double(config.DEFAULT.ntJ);

        % polynomial (order 4 is cubic)
        order = str2double(config.DEFAULT.order);
        
        
    end
    
    
    if downloop >= 2
        % now down 2
        Afilename = [output_dir 'down4_A.mat'];
        vfilename = [output_dir 'down4_v.mat'];
        prefix = [output_dir 'down2_'];
        
        eVJ = str2double(config.DOWNLOOP2.eVJ);
        eVI = str2double(config.DOWNLOOP2.eVI);
        post_affine_reduce = str2double(config.DOWNLOOP2.post_affine_reduce);
        start_3d_diffeo = str2double(config.DOWNLOOP2.start_3d_diffeo);
        start_2d_diffeo = str2double(config.DOWNLOOP2.start_2d_diffeo);
        start_2d_affine = str2double(config.DOWNLOOP2.start_2d_affine);
        niter = str2double(config.DOWNLOOP2.niter);

    end
    
    if downloop >= 3
        % now down 1
        eVI = str2double(config.DOWNLOOP2.eVI);
        eVJ = str2double(config.DOWNLOOP2.eVJ);
        post_affine_reduce = str2double(config.DOWNLOOP3.post_affine_reduce);
        
        prefix = [output_dir 'down1_'];
        Afilename = [output_dir 'down2_A.mat'];
        vfilename = [output_dir 'down2_v.mat'];
        
        start_3d_diffeo = str2double(config.DOWNLOOP3.start_3d_diffeo);
        start_2d_diffeo = str2double(config.DOWNLOOP3.start_2d_diffeo);
        start_2d_affine = str2double(config.DOWNLOOP3.start_2d_affine);
        niter = str2double(config.DOWNLOOP3.niter);

    end
    if edit_mode
        if downloop < 3
            disp(['running in edit mode only on high res'])
            continue
        else
            sigmaL = 0.1;
            prefix = [output_dir 'down1edit_'];
            Afilename = [output_dir 'down1_A.mat'];
            vfilename = [output_dir 'down1_v.mat'];
        end
    end
    
    % hack if stepsizes are too big
    % % I want to reduce them all by 2 since twice as many slices in 787
    % eVJ = eVJ/2; % this is reduced from before
    % eVI = eVI/2;
    % eTI = eTI/2;
    % eLI = eLI/2;
    % eLJ = eLJ/2;
    % eTJ = eTJ/2;
    
    
    
    % energy
    EAsave = [];
    EBsave = [];
    frameJS = [];
    framefIS = [];
    frameWS = [];
    vmaxsave = [];
    Asave = [];
    framePhiI = [];
    frameErrAll = [];
    frameWeightAll = [];
    frameIAll = [];
    if edit_mode
        frameLabelAll = [];
    end
    ERJsave = zeros(length(files),niter);
    Esave = zeros(1,niter);
    EMsave = zeros(1,niter);
    ERsave = zeros(1,niter);
    EMJsave = zeros(length(files),niter);
    EMall = zeros(length(files),1);
    
    
    % note about operators
    % I want to define the length scales a in a way that does not depend on
    % downsampling
    % but
    % I want conditioning in a way that does
    % flow for atlas
    dt = 1.0/nt;
    vtx = zeros([nxI([2,1,3]),nt]);
    vty = zeros([nxI([2,1,3]),nt]);
    vtz = zeros([nxI([2,1,3]),nt]);
    LL = ( 1.0 - a^2*2*( (cos(2.0*pi*FXI*dxI(1))-1.0)/dxI(1)^2 + (cos(2.0*pi*FYI*dxI(2))-1.0)/dxI(2)^2 + (cos(2.0*pi*FZI*dxI(3))-1.0)/dxI(3)^2  ) ).^(2*p);
    K = 1.0./LL;
    % Kp for preconditioner
    LLp = ( 1.0 - (ap)^2*2*( (cos(2.0*pi*FXI*dxI(1))-1.0)/dxI(1)^2 + (cos(2.0*pi*FYI*dxI(2))-1.0)/dxI(2)^2 + (cos(2.0*pi*FZI*dxI(3))-1.0)/dxI(3)^2  ) ).^(2*p);
    Kp = LL./LLp;
    
    % 2d flow for target
    dtJ = 1.0/ntJ;
    for f = 1 : length(files)
        vJtx{f} = zeros([nxJ{f}([2,1]),ntJ]);
        vJty{f} = zeros([nxJ{f}([2,1]),ntJ]);
        if nonrigid_thick_only && ~isthicks(f) % for non-thick slices no flow
            vJtx{f} = zeros([nxJ{f}([2,1]),0]);
            vJty{f} = zeros([nxJ{f}([2,1]),0]);
        end
        fxJ = (0:nxJ{f}(1)-1)/nxJ{f}(1)/dxJ(1);
        fyJ = (0:nxJ{f}(2)-1)/nxJ{f}(2)/dxJ(2);
        [FXJ,FYJ] = meshgrid(fxJ,fyJ);
        
        LLJ{f} = ( 1.0 - aJ^2*2*( (cos(2.0*pi*FXJ*dxJ(1))-1.0)/dxJ(1)^2 + (cos(2.0*pi*FYJ*dxJ(2))-1.0)/dxJ(2)^2 )).^(2*pJ);
        KJ{f} = 1.0./LLJ{f};
        LLJp{f} = ( 1.0 - (aJp)^2*2*( (cos(2.0*pi*FXJ*dxJ(1))-1.0)/dxJ(1)^2 + (cos(2.0*pi*FYJ*dxJ(2))-1.0)/dxJ(2)^2 )).^(2*pJ);
        KJp{f} = LLJ{f}./LLJp{f};
    end
    
    % affine for atlas
    A = eye(4);
    % permute xy, and flip z
    % this is only valid for Allen vtk atlas
    A = [0,1,0,0;
        1,0,0,0;
        0,0,-1,0;
        0,0,0,1]*A;
    
    
    
    if ~isempty(Afilename)
        load(Afilename);
        if strcmp(Afilename, initializer)
            % first time we have to do something about the mean because I'm
            % only working with zero mean AJ
            % find the mean
            logmAJ = zeros(size(AJ));
            AJ(isnan(AJ)) = 0; % get rid of any nans
            for i = 1 : size(logmAJ,3)                
                center(i,:) = AJ(1:2,3,i);
                tmp = logm(AJ(:,:,i));
                if any(imag(tmp(:))) % if any are imaginary, we just set to identity and look at translation only
                    AJ(1:2,1:2,i) = eye(2);
                    tmp = logm(AJ(:,:,i));
                end
                logmAJ(:,:,i) = tmp;
            end
            themean = mean(logmAJ,3);
            themeanAJ = expm(themean);
            % we want to "subtract" the mean from AJ and "add" it to A
            for i = 1 : size(logmAJ,3)
                AJ(:,:,i) = AJ(:,:,i)/themeanAJ; % inverse on the right
            end
            A = [themeanAJ(1,1:2),0,themeanAJ(1,3);
                themeanAJ(2,1:2),0,themeanAJ(2,3);
                0,0,1,0;
                0,0,0,1]*A; % forward on the left
        end
        
    end
    if ~isempty(vfilename)
        % do not load the xI's because they have been redefined above
        load(vfilename,'vtx', 'vty', 'vtz', 'vJtx', 'vJty')
    end
    
    % upsample the vs
    if any (size(vtx(:,:,:,1)) < size(I))
        vtx_ = vtx;
        vty_ = vty;
        vtz_ = vtz;
        vtx = zeros([size(I),nt]);
        vty = zeros([size(I),nt]);
        vtz = zeros([size(I),nt]);
        for t = 1 : nt
            tmp = upsample(vtx_(:,:,:,t),size(I));
            vtx(:,:,:,t) = tmp;
            tmp = upsample(vty_(:,:,:,t),size(I));
            vty(:,:,:,t) = tmp;
            tmp = upsample(vtz_(:,:,:,t),size(I));
            vtz(:,:,:,t) = tmp;
        end
    end
    % now go through every slice
    for i = 1 : length(files)
        if size(vJtx{i},3) > 0 && any(size(vJtx{i}(:,:,1)) < size(J{i}(:,:,1)))
            vtx_ = vJtx{i};
            vty_ = vJty{i};
            vJtx{i} = zeros([size(J{i},1),size(J{i},2),ntJ]);
            vJty{i} = zeros([size(J{i},1),size(J{i},2),ntJ]);
            for t = 1 : ntJ
                tmp = upsample(vtx_(:,:,t),[size(J{i},1),size(J{i},2),1]);
                vJtx{i}(:,:,t) = tmp;
                tmp = upsample(vty_(:,:,t),[size(J{i},1),size(J{i},2),1]);
                vJty{i}(:,:,t) = tmp;
                % multiply by 0 here if I don't want to do any 2d diffeo
            end
        end
        
    end
    
    % initialize norm grad variables for display
    normLJsave = zeros(length(files),niter);
    normLJstepsave = zeros(length(files),niter);
    normTJsave = zeros(length(files),niter);
    normTJstepsave = zeros(length(files),niter);
    normLIsave = zeros(1,niter);
    normLIstepsave = zeros(1,niter);
    normTIsave = zeros(1,niter);
    normTIstepsave = zeros(1,niter);
    
    [I_x,I_y,I_z] = gradient3d(I,dxI(1),dxI(2),dxI(3),1);
    if edit_mode
        L_x = zeros(size(L,1),size(L,2),size(L,3),nL);
        L_y = zeros(size(L,1),size(L,2),size(L,3),nL);
        L_z = zeros(size(L,1),size(L,2),size(L,3),nL);
        for l = 1 : nL
            [L_x(:,:,:,l),L_y(:,:,:,l),L_z(:,:,:,l)] = gradient3d(L(:,:,:,l),dxI(1),dxI(2),dxI(3),1);
        end
    end
    
    % if running in edit mode
    % we should load a saved velocity here and skip the low res versions
    
    

    for it = 1 : niter
        
        save_frames = it < 10 || ~mod(it-1,11) || it == niter;
        
        %% first the deformation
        phiinvx = XI;
        phiinvy = YI;
        phiinvz = ZI;
        It = repmat(I,[1,1,1,nt]);
        if edit_mode
            Lt = repmat(L,[1,1,1,1,nt]);
        end
        for t = 1 : nt*(it>start_3d_diffeo)
            % deform image
            if t > 1
                F = griddedInterpolant({yI,xI,zI},I,'linear','nearest');
                It(:,:,:,t) = F(phiinvy,phiinvx,phiinvz);
                if edit_mode
                    for l = 1 : nL
                        F = griddedInterpolant({yI,xI,zI},L(:,:,:,l),'linear','nearest');
                        Lt(:,:,:,l,t) = F(phiinvy,phiinvx,phiinvz);
                    end
                end
            end
            % update phi
            Xs = XI - vtx(:,:,:,t)*dt;
            Ys = YI - vty(:,:,:,t)*dt;
            Zs = ZI - vtz(:,:,:,t)*dt;
            % subtract and add identity
            F = griddedInterpolant({yI,xI,zI},phiinvx-XI,'linear','nearest');
            phiinvx = F(Ys,Xs,Zs) + Xs;
            F = griddedInterpolant({yI,xI,zI},phiinvy-YI,'linear','nearest');
            phiinvy = F(Ys,Xs,Zs) + Ys;
            F = griddedInterpolant({yI,xI,zI},phiinvz-ZI,'linear','nearest');
            phiinvz = F(Ys,Xs,Zs) + Zs;
            
        end
        % cost
        vtxhat = fft(fft(fft(vtx,[],1),[],2),[],3);
        vtyhat = fft(fft(fft(vty,[],1),[],2),[],3);
        vtzhat = fft(fft(fft(vtz,[],1),[],2),[],3);
        ER = sum(sum(sum(sum( bsxfun(@times, (abs(vtxhat).^2 + abs(vtyhat).^2 + abs(vtzhat).^2) ,  LL ) ))))*prod(dxI)*dt/numel(vtx(:,:,:,1))/sigmaR^2/2;
        ERsave(it) = ER;
        
        
        % I'd like to sample this in the same space so I can take its gradient
        % for later
        F = griddedInterpolant({yI,xI,zI},I,'linear','nearest');
        phiI = F(phiinvy,phiinvx,phiinvz);
%         [phiI_x,phiI_y,phiI_z] = gradient3d(phiI,dxI(1),dxI(2),dxI(3));
        % note this may not be accurate if there are boundary issues
        F = griddedInterpolant({yI,xI,zI},I_x,'linear','nearest');
        phi_I_x = F(phiinvy,phiinvx,phiinvz);
        F = griddedInterpolant({yI,xI,zI},I_y,'linear','nearest');
        phi_I_y = F(phiinvy,phiinvx,phiinvz);
        F = griddedInterpolant({yI,xI,zI},I_z,'linear','nearest');
        phi_I_z = F(phiinvy,phiinvx,phiinvz);
        [phiinvx_x, phiinvx_y, phiinvx_z] = gradient3d(phiinvx,dxI(1),dxI(2),dxI(3));
        [phiinvy_x, phiinvy_y, phiinvy_z] = gradient3d(phiinvy,dxI(1),dxI(2),dxI(3));
        [phiinvz_x, phiinvz_y, phiinvz_z] = gradient3d(phiinvz,dxI(1),dxI(2),dxI(3));
        phiI_x = phi_I_x.*phiinvx_x + phi_I_y.*phiinvy_x + phi_I_z.*phiinvz_x;
        phiI_y = phi_I_x.*phiinvx_y + phi_I_y.*phiinvy_y + phi_I_z.*phiinvz_y;
        phiI_z = phi_I_x.*phiinvx_z + phi_I_y.*phiinvy_z + phi_I_z.*phiinvz_z;
        
        
        
        danfigure(333);
        sliceView(xI,yI,zI,phiI)
        
        if edit_mode
            phiL = zeros([size(phiI),nL]);
            phiL_x = zeros([size(phiI),nL]);
            phiL_y = zeros([size(phiI),nL]);
            phiL_z = zeros([size(phiI),nL]);
            for l = 1 : nL
                F = griddedInterpolant({yI,xI,zI},L(:,:,:,l),'linear','nearest');
                phiL(:,:,:,l) = F(phiinvy,phiinvx,phiinvz);
                %         [phiI_x,phiI_y,phiI_z] = gradient3d(phiI,dxI(1),dxI(2),dxI(3));
                % note this may not be accurate if there are boundary issues
                F = griddedInterpolant({yI,xI,zI},L_x(:,:,:,l),'linear','nearest');
                phi_L_x = F(phiinvy,phiinvx,phiinvz);
                F = griddedInterpolant({yI,xI,zI},L_y(:,:,:,l),'linear','nearest');
                phi_L_y = F(phiinvy,phiinvx,phiinvz);
                F = griddedInterpolant({yI,xI,zI},L_z(:,:,:,l),'linear','nearest');
                phi_L_z = F(phiinvy,phiinvx,phiinvz);
                
                phiL_x(:,:,:,l) = phi_L_x.*phiinvx_x + phi_L_y.*phiinvy_x + phi_L_z.*phiinvz_x;
                phiL_y(:,:,:,l) = phi_L_x.*phiinvx_y + phi_L_y.*phiinvy_y + phi_L_z.*phiinvz_y;
                phiL_z(:,:,:,l) = phi_L_x.*phiinvx_z + phi_L_y.*phiinvy_z + phi_L_z.*phiinvz_z;
            end
        end
        
        %% next sample onto each slice
        % the sequence of transformations is
        % 1. phi
        % 2. A (3D pose)
        % 3. phiJ (2d diffeo if applicable)
        % 4. AJ (slice specific)
        % I want to save the norm of the gradients of each slice
        B = inv(A);
        
        % initialize the gradient. each iteration of the loop wlil add to it
        gradA = zeros(4,4);
        gradI = zeros(size(I));
        if edit_mode
            gradL = zeros(size(L));
        end
        
        for f = 1 : length(files)
            % parfor f = 1 : length(files) % haven't tried yet
            edit_this_slice = edit_mode && ~isempty(edit_info{f});
            if edit_this_slice
                channels = [edit_info{f}.channel];
            end
            % nonrigid deformation on this slice
            % if it is thick (always)
            % or if you are not using the nonrigid_thick_only option
            this_slice_nonrigid = (~nonrigid_thick_only || isthicks(f));            
            AJf = AJ(:,:,f);
            %%
            % now a diffeo in 2D
            phiJinvx = XJ{f};
            phiJinvy = YJ{f};
            phiJAphiIt = zeros([size(XJ{f}),ntJ]);
            if edit_this_slice
                phiJAphiLt = zeros([size(XJ{f}),nLJ(f),ntJ]);
            end
            
            for t = 1 : ntJ*(it>start_2d_diffeo)*this_slice_nonrigid
                % sample the image at this point
                % phiJ . A . phi . I
                % or I(phi^{-1}(A^{-1}(phiJ^{-1}(x))))
                % consider phiJ to be identity and constant in z
                % first we apply A
                AiPhiJiX = B(1,1)*phiJinvx + B(1,2)*phiJinvy + B(1,3)*zJ(f) + B(1,4);
                AiPhiJiY = B(2,1)*phiJinvx + B(2,2)*phiJinvy + B(2,3)*zJ(f) + B(2,4);
                AiPhiJiZ = B(3,1)*phiJinvx + B(3,2)*phiJinvy + B(3,3)*zJ(f) + B(3,4);
                % then apply phiinv
                F = griddedInterpolant({yI,xI,zI},phiinvx-XI,'linear','nearest');
                phiiAiPhiJiX = F(AiPhiJiY,AiPhiJiX,AiPhiJiZ) + AiPhiJiX;
                F = griddedInterpolant({yI,xI,zI},phiinvy-YI,'linear','nearest');
                phiiAiPhiJiY = F(AiPhiJiY,AiPhiJiX,AiPhiJiZ) + AiPhiJiY;
                F = griddedInterpolant({yI,xI,zI},phiinvz-ZI,'linear','nearest');
                phiiAiPhiJiZ = F(AiPhiJiY,AiPhiJiX,AiPhiJiZ) + AiPhiJiZ;
                F = griddedInterpolant({yI,xI,zI},I,'linear','nearest');
                phiJAphiIt(:,:,t) = F(phiiAiPhiJiY,phiiAiPhiJiX,phiiAiPhiJiZ);
                if edit_this_slice
                    for l = 1 : length(channels)
                        F = griddedInterpolant({yI,xI,zI},L(:,:,:,channels(l)),'linear','nearest');
                        phiJAphiLt(:,:,l,t) = F(phiiAiPhiJiY,phiiAiPhiJiX,phiiAiPhiJiZ);
                    end
                end
                
                % update phi
                Xs = XJ{f} - vJtx{f}(:,:,t)*dtJ;
                Ys = YJ{f} - vJty{f}(:,:,t)*dtJ;
                % subtract and add identity
                F = griddedInterpolant({yJ{f},xJ{f}},phiJinvx-XJ{f},'linear','nearest');
                phiJinvx = F(Ys,Xs) + Xs;
                F = griddedInterpolant({yJ{f},xJ{f}},phiJinvy-YJ{f},'linear','nearest');
                phiJinvy = F(Ys,Xs) + Ys;
            end
            % and the last timepoint (I can use it for the AJ gradient)
            % first we apply A
            AiPhiJiX = B(1,1)*phiJinvx + B(1,2)*phiJinvy + B(1,3)*zJ(f) + B(1,4);
            AiPhiJiY = B(2,1)*phiJinvx + B(2,2)*phiJinvy + B(2,3)*zJ(f) + B(2,4);
            AiPhiJiZ = B(3,1)*phiJinvx + B(3,2)*phiJinvy + B(3,3)*zJ(f) + B(3,4);
            % then apply phiinv
            F = griddedInterpolant({yI,xI,zI},phiinvx-XI,'linear','nearest');
            phiiAiPhiJiX = F(AiPhiJiY,AiPhiJiX,AiPhiJiZ) + AiPhiJiX;
            F = griddedInterpolant({yI,xI,zI},phiinvy-YI,'linear','nearest');
            phiiAiPhiJiY = F(AiPhiJiY,AiPhiJiX,AiPhiJiZ) + AiPhiJiY;
            F = griddedInterpolant({yI,xI,zI},phiinvz-ZI,'linear','nearest');
            phiiAiPhiJiZ = F(AiPhiJiY,AiPhiJiX,AiPhiJiZ) + AiPhiJiZ;
            F = griddedInterpolant({yI,xI,zI},I,'linear','nearest');
            phiJAphiI = F(phiiAiPhiJiY,phiiAiPhiJiX,phiiAiPhiJiZ);
            if edit_this_slice
                phiJAphiL = zeros([size(phiJAphiI),nL]);
                for l = 1 : nL
                    F = griddedInterpolant({yI,xI,zI},L(:,:,:,l),'linear','nearest');
                    phiJAphiL(:,:,l) = F(phiiAiPhiJiY,phiiAiPhiJiX,phiiAiPhiJiZ);
                end
            end
            
            % cost
            if size(vJtx{f},3) > 1
                vJtxhat = fft(fft(vJtx{f},[],1),[],2);
                vJtyhat = fft(fft(vJty{f},[],1),[],2);
                
                ERJ(f) = sum(sum(sum( bsxfun(@times, (abs(vJtxhat).^2 + abs(vJtyhat).^2 ) ,  LLJ{f} ) )))*prod(dxJ([1,2]))*dtJ/numel(vJtx{f}(:,:,1))/sigmaRJ^2/2;
            else
                ERJ(f) = 0;
            end
            %             ERJsave = [ERJsave,ERJ];
            ERJsave(f,it) = ERJ(f);
            
            
            
            % now the 2D affine transform
            % we apply
            % x -> AJ phiJ A phi x
            % so the inverse is
            % phi^{-1}(A^{-1}(phiJ^{-1}(AJ^{-1}(x))))
            % first AJ
            BJ = inv(AJf);
            BJ_ = [BJ(1:2,1:2),[0;0],BJ(1:2,3);0,0,1,0;0,0,0,1];
            AJiX = BJ_(1,1)*XJ{f} + BJ_(1,2)*YJ{f} + BJ_(1,3)*zJ(f) + BJ_(1,4);
            AJiY = BJ_(2,1)*XJ{f} + BJ_(2,2)*YJ{f} + BJ_(2,3)*zJ(f) + BJ_(2,4);
            AJiZ = BJ_(3,1)*XJ{f} + BJ_(3,2)*YJ{f} + BJ_(3,3)*zJ(f) + BJ_(3,4); % just ZJ constant
            % now phiJ
            F = griddedInterpolant({yJ{f},xJ{f}},phiJinvx-XJ{f},'linear','nearest');
            phiJiAJiX = F(AJiY,AJiX) + AJiX;
            F = griddedInterpolant({yJ{f},xJ{f}},phiJinvy-YJ{f},'linear','nearest');
            phiJiAJiY = F(AJiY,AJiX) + AJiY;
            phiJiAJiZ = AJiZ; % identity
            % now A
            AiPhiJiAJiX = B(1,1)*phiJiAJiX + B(1,2)*phiJiAJiY + B(1,3)*phiJiAJiZ + B(1,4);
            AiPhiJiAJiY = B(2,1)*phiJiAJiX + B(2,2)*phiJiAJiY + B(2,3)*phiJiAJiZ + B(2,4);
            AiPhiJiAJiZ = B(3,1)*phiJiAJiX + B(3,2)*phiJiAJiY + B(3,3)*phiJiAJiZ + B(3,4);
            % now phi
            F = griddedInterpolant({yI,xI,zI},phiinvx-XI,'linear','nearest');
            phiiAiPhiJiAJiX = F(AiPhiJiAJiY,AiPhiJiAJiX,AiPhiJiAJiZ) + AiPhiJiAJiX;
            F = griddedInterpolant({yI,xI,zI},phiinvy-YI,'linear','nearest');
            phiiAiPhiJiAJiY = F(AiPhiJiAJiY,AiPhiJiAJiX,AiPhiJiAJiZ) + AiPhiJiAJiY;
            F = griddedInterpolant({yI,xI,zI},phiinvz-ZI,'linear','nearest');
            phiiAiPhiJiAJiZ = F(AiPhiJiAJiY,AiPhiJiAJiX,AiPhiJiAJiZ) + AiPhiJiAJiZ;
            
            
            % now get the image
            F = griddedInterpolant({yI,xI,zI},I,'linear','nearest');
            SAphiI{f} = F(phiiAiPhiJiAJiY,phiiAiPhiJiAJiX,phiiAiPhiJiAJiZ);
            if edit_this_slice
                SAphiL{f} = zeros([size(SAphiI{f}),nLJ(f)]);
                for l = 1 : nLJ(f)
                    F = griddedInterpolant({yI,xI,zI},L(:,:,:,channels(l)),'linear','nearest');
                    SAphiL{f}(:,:,l) = F(phiiAiPhiJiAJiY,phiiAiPhiJiAJiX,phiiAiPhiJiAJiZ);
                end
            end
            
            %% now contrast mapping
            % basis functions and their derivative
            D = zeros(numel(SAphiI{f}),order); % data to map
            DD = zeros(numel(SAphiI{f}),order); % derivative for adjoint calculations
            for o = 1 : order
                D(:,o) = SAphiI{f}(:).^(o-1);
                if o > 1
                    DD(:,o) = SAphiI{f}(:).^(o-2)*(o-1);
                end
            end
            
            
            if (it == 1)% this is just for a first initialization!
                WM{f} = piM*ones(size(J{f},1),size(J{f},2));
                WA{f} = piA*ones(size(J{f},1),size(J{f},2));
                WB{f} = piB*ones(size(J{f},1),size(J{f},2));
                
                % get coeffs first time, otherwise I'll get them at the end
                % probably I should have a way of loading them
                coeffs{f} = (D'*bsxfun(@times,D,WM{f}(:).*WMask{f}(:)) + eye(order)*1e-3)\(D'*bsxfun(@times,reshape(J{f},numel(SAphiI{f}),size(J{f},3)),WM{f}(:).*WMask{f}(:)));
                if any(isnan(coeffs{f}(:)))
                    coeffs{f} = zeros(size(coeffs{f}));
                end
                CA{f} = sum(sum(bsxfun(@times, J{f},WA{f}.*WMask{f}),1),2)/sum(WA{f}(:).*WMask{f}(:));
                if any(isnan(CA{f}))
                CA{f} = zeros(size(CA{f}));
                end
                CB{f} = sum(sum(bsxfun(@times, J{f},WB{f}.*WMask{f}),1),2)/sum(WB{f}(:).*WMask{f}(:));
                if any(isnan(CB{f}))
                CB{f} = zeros(size(CB{f}));
                end
            end
            
            
            
            fSAphiI{f} = reshape(D*coeffs{f},size(J{f}));

            %%
            % update weights (I need to calculate fSAphiI first)
            % I like updating the weights all the time for the rigid stage
            % note
            % what happens if a slice is entirely out of bounds?
            % then we have to predict with three constants
            % it seems like I run into trouble and end up with NaNs
            % why?
            % well, if the image is constant
            % then I can predict it exactly using the background channel
            % then the matching channel will have zero weight
            % then when I estimate coeffs I get nans
            % I think that might be the issue
            
            if 1 % for now I'm doing this every time, later I should do it every "mstep" steps
                WM{f} = 1/sqrt(2*pi*sigmaM^2)^(size(J{f},3))*exp(-sum((fSAphiI{f} - J{f}).^2,3)/2/sigmaM^2);
                WA{f} = 1/sqrt(2*pi*sigmaA^2)^(size(J{f},3))*exp(-sum(bsxfun(@minus,CA{f},J{f}).^2,3)/2/sigmaA^2);
                WB{f} = 1/sqrt(2*pi*sigmaB^2)^(size(J{f},3))*exp(-sum(bsxfun(@minus,CB{f},J{f}).^2,3)/2/sigmaB^2);
                %                 WSum = WM{f} + WA{f} + WB{f};
                %                 WM{f} = WM{f}./WSum;
                %                 WA{f} = WA{f}./WSum;
                %                 WB{f} = WB{f}./WSum;
                % let's use priors, see if this helps
                % I don't think this is working very well
                % somethings that are obviously missing tissue are not seen
                % that way
                % basically waht I'm saying is really try to match!
                % we'll see how this goes
                % maybe a spatial prior is better
                spatial_prior = exp(-((XJ{f}-center(f,1)).^2 + (YJ{f}-center(f,2)).^2)/2/piM_width^2);
                piM_ = piM*spatial_prior;
                piB_ = piB*(1-spatial_prior);
                piA_ = piA;
                
                %                 WSum = piM.*WM{f} + piA.*WA{f} + piB.*WB{f};
                %                 WM{f} = piM.*WM{f}./WSum;
                %                 WA{f} = piA.*WA{f}./WSum;
                %                 WB{f} = piB.*WB{f}./WSum;
                
                WSum = piM_.*WM{f} + piA_.*WA{f} + piB_.*WB{f};
                WSum = WSum + max(WSum(:))*1e-6; % try to avoid dividing by zero
                WM{f} = piM_.*WM{f}./WSum;
                WA{f} = piA_.*WA{f}./WSum;
                WB{f} = piB_.*WB{f}./WSum;
                
            end
                        
            
            %% now evaluate errors
            err{f} = fSAphiI{f} - J{f};
            EMall(f) = sum(sum( sum(err{f}.^2,3).*WM{f}.*WMask{f}))/2/sigmaM^2*prod(dxJ(f,1:2));
            
            if edit_this_slice
                errL{f} = SAphiL{f} - LJ{f};
                % note LJ{f} is a weight as well as a target
                EMall(f) = EMall(f) + sum(sum(sum(sum( errL{f}.^2.*LJ{f},3 ) .*WMask{f})))/2/sigmaL^2*prod(dxJ(f,1:2));
            end
            
            %% now the adjoint with respect to the deformation
            % I is transformed by f as
            % fI_i = f(I) = sum_j c_{ij} I^(j)
            % so perturbations are transformed by
            % consider perturbations I \mapsto I + e dI
            % d_de f(I + e dI)_i |_{e = 0} = d_de sum_j c_{ij}(I + e dI)^j |_{e=0}
            % = (sum_j c_{ij} j I^{j-1})dI % linear in dI (of course)
            % to find the adjoint
            % we consider
            % sum(err * f_*(I)[dI]) = sum( f^*(I)[err] * dI )
            tmp = sum((DD*coeffs{f}).*reshape(err{f},size(J{f},1)*size(J{f},2),size(J{f},3)),2);
            ferr{f} = reshape(tmp,size(SAphiI{f}));
            ferrW{f} = ferr{f}.*WM{f}.*WMask{f};
            if edit_this_slice
                % no contrast transform for the labels
                errLW{f} = bsxfun(@times,errL{f}.*LJ{f},WMask{f});
            end            
            %% now evaluate 2D linear transform gradient
            % we have the following cost
            % E = int |f(I(phi^{-1}(A^{-1}(phiJ^{-1}(AJ^{-1}(x)))))) - J(x)|^2 dx
            % let I' = I(phi^{-1}(A^{-1}(phiJ^{-1}(x))))
            % E = int |f(I'(AJ^{-1}(x)) - J(x)|^2 dx
            % take perturbation
            % d_de int |f(I'(exp(-e dA)(AJ^{-1}(x))) - J(x)|^2 dx e=0
            % = int 2(f(I'(AJ^{-1}(x))) - J(x))^T Df(I'(AJ^{-1}(x)))DI'(AJ^{-1}(x))(-1)dA AJ^{-1}x  dx
            % = int 2(f(I'(AJ^{-1}(x))) - J(x))^T Df(I'(AJ^{-1}(x)))D[I'(AJ^{-1}(x))]AJ(-1)dA AJ^{-1}x  dx
            %             [I_x_ , I_y_ ] = gradient(SAphiI{f},dxJ(1),dxJ(2));
            
            [I_x_ , I_y_ ] = gradient2d(SAphiI{f},dxJ(1),dxJ(2));
            if edit_this_slice
                L_x_ = zeros([size(I_x_),nLJ(f)]);
                L_y_ = zeros([size(I_x_),nLJ(f)]);
                for l = 1 : nLJ(f)
                    [L_x_(:,:,l) , L_y_(:,:,l) ] = gradient2d(SAphiL{f}(:,:,l),dxJ(1),dxJ(2));
                end
            end
            gradAJ = zeros(4,4);
            AJ_ = [AJf(1:2,1:2),[0;0],AJf(1:2,3);0,0,1,0;0,0,0,1];
            for r = 1 : 2
                for c = [1,2,4]
                    dA = double((1:4)==r)' * double((1:4)==c);
                    Aall = AJ_*dA*BJ_;
                    Xs_ = Aall(1,1)*XJ{f} + Aall(1,2)*YJ{f} + Aall(1,3)*zJ(f)+ Aall(1,4);
                    Ys_ = Aall(2,1)*XJ{f} + Aall(2,2)*YJ{f} + Aall(2,3)*zJ(f)+ Aall(2,4);
                    tmp = ferrW{f}.*(I_x_.*Xs_ + I_y_.*Ys_);
                    
                    gradAJ(r,c) = sum(sum(tmp))*prod(dxJ(1:2))*(-1)/sigmaM^2;
                    if edit_this_slice
                        tmpL = sum(errLW{f} .* L_x_ ,3).*Xs_ + sum(errLW{f} .* L_y_ ,3).*Ys_;
                        gradAJ(r,c) = gradAJ(r,c) + sum(sum(tmpL))*prod(dxJ(1:2))*(-1)/sigmaL^2;
                    end
                end
            end
            grad_AJ_f = gradAJ(1:3,[1,2,4]);
            % if rigid
            grad_AJ_f(1:2,1:2) = grad_AJ_f(1:2,1:2) - grad_AJ_f(1:2,1:2)';
            gradAJall(:,:,f) = grad_AJ_f;
            
            
            
            
            %%
            % now find the gradient of the deformation
            % we've got some error
            % and we need to pull it back
            % using AJ and adjoint of phiJ
            % start my flow
            ferrWpad = zeros(nxJ{f}(2)+2, nxJ{f}(1)+2);
            ferrWpad(2:end-1,2:end-1) = ferrW{f};            
            xJpad = [xJ{f}(1)-dxJ(1), xJ{f}, xJ{f}(end)+dxJ(1)];
            yJpad = [yJ{f}(1)-dxJ(2), yJ{f}, yJ{f}(end)+dxJ(2)];
            if edit_this_slice
                errLWpad = zeros(nxJ{f}(2)+2, nxJ{f}(1)+2,nLJ(f));
                errLWpad(2:end-1,2:end-1,:) = errLW{f};  
            end
            phiJ1tix = XJ{f};
            phiJ1tiy = YJ{f};
            detjacJ = ones(size(XJ{f}));
            for t = ntJ*(it>start_2d_diffeo)*this_slice_nonrigid : -1 : 1
                % update diffeo (note plus)
                Xs = XJ{f} + vJtx{f}(:,:,t)*dtJ;
                Ys = YJ{f} + vJty{f}(:,:,t)*dtJ;
                
                F = griddedInterpolant({yJ{f},xJ{f}},phiJ1tix-XJ{f},'linear','nearest');
                phiJ1tix = F(Ys,Xs) + Xs;
                F = griddedInterpolant({yJ{f},xJ{f}},phiJ1tiy-YJ{f},'linear','nearest');
                phiJ1tiy = F(Ys,Xs) + Ys;
                
                % determinant of jacobian
                [phiJ1tix_x,phiJ1tix_y] = gradient2d(phiJ1tix,dxJ(1),dxJ(2));
                [phiJ1tiy_x,phiJ1tiy_y] = gradient2d(phiJ1tiy,dxJ(1),dxJ(2));
                detjacJ = phiJ1tix_x.*phiJ1tiy_y - phiJ1tix_y.*phiJ1tiy_x;
                
                % affine
                AJphiJ1tinvx = AJf(1,1)*phiJ1tix + AJf(1,2)*phiJ1tiy + AJf(1,3);
                AJphiJ1tinvy = AJf(2,1)*phiJ1tix + AJf(2,2)*phiJ1tiy + AJf(2,3);
                
                % resample error
                F = griddedInterpolant({yJpad,xJpad},ferrWpad,'linear','nearest');
                lambdat = F(AJphiJ1tinvy,AJphiJ1tinvx).*detjacJ*abs(det(AJf))*(-1)/sigmaM^2;
                
                % gradient of It
                [It_x,It_y] = gradient2d(phiJAphiIt(:,:,t),dxJ(1),dxJ(2));
                
                % matching function gradient
                gradx = It_x.*lambdat;
                grady = It_y.*lambdat;
                
                if edit_this_slice
                    warning('updating 2D deformations in edit mode not yet implemented')
                    % note not tested because MBA does not use slice
                    % deformation
                    
                    lambdat = zeros([size(lambdat),nLJ(f)]);
                    for l = 1 : nLJ(f)
                        F = griddedInterpolant({yJpad,xJpad},errLWpad(:,:,l),'linear','nearest');
                        lambdat(:,:,l) = F(AJphiJ1tinvy,AJphiJ1tinvx).*detjacJ*abs(det(AJf))*(-1)/sigmaL^2;
                    end
                    
                    % gradient of It
                    Lt_x = zeros([size(lambdat)]);
                    Lt_y = zeros([size(lambdat)]);
                    for l = 1 : nLJ(f)
                        [Lt_x(:,:,l),Lt_y(:,:,l)] = gradient2d(phiJAphiLt(:,:,l,t),dxJ(1),dxJ(2));
                    end
                    % matching function gradient
                    gradx = gradx + sum(Lt_x.*lambdat,3);
                    grady = grady + sum(Lt_y.*lambdat,3);
                end
                
                % regularization ( I could add more smoothness here if I wanted)
                gradx = ifftn(   KJp{f}.*(fftn(gradx).*KJ{f} + vJtxhat(:,:,t)/sigmaRJ^2)   ,'symmetric');
                grady = ifftn(   KJp{f}.*(fftn(grady).*KJ{f} + vJtyhat(:,:,t)/sigmaRJ^2)   ,'symmetric');
                
                stepx = gradx*eVJ;
                stepy = grady*eVJ;
                
                % squash it
                normstep = sqrt(stepx.^2 + stepy.^2);
                if ~(all(normstep(:))==0) % make sure its not all zero, if it is we're just gonna get nans here
                    squashnorm = VJmax*normstep./(VJmax + normstep);
                    stepx = stepx./normstep.*squashnorm;
                    stepy = stepy./normstep.*squashnorm;
                    
                    vJtx{f}(:,:,t) = vJtx{f}(:,:,t) - stepx;
                    vJty{f}(:,:,t) = vJty{f}(:,:,t) - stepy;
                end
            end
            
            
            
            %%
            % now update the affine for this slice (maybe should be at end)
            % update affine
            e = ([1;1;0]*[1,1,0]*eLJ + [1;1;0]*[0,0,1]*eTJ);
            if it > min(start_2d_diffeo,start_3d_diffeo)
                e = e*post_affine_reduce;
            end
            %             AJ(:,:,f) = AJ(:,:,f)*expm(-gradAJall(:,:,f).*e);
            % normalize, use a simple squashing function with reasonable
            % max, which is identity near 0
            % how about this, really simple!
            % f(x) = ax/(a + x)
            % for large x it goes to a
            % for small x it is x
            stepAJ = grad_AJ_f.*e;
            normLJ = sqrt(sum(sum(stepAJ(1:2,1:2).^2)));
            normLJsave(f,it) = sqrt(sum(sum(grad_AJ_f(1:2,1:2).^2)));
            normLJstepsave(f,it) = normLJ;
            if normLJ > 0
                normLJsquash = LJmax*normLJ/(LJmax + normLJ);
                stepAJ(1:2,1:2) = stepAJ(1:2,1:2)/normLJ*normLJsquash;
            end
            
            normTJ = sqrt(sum(stepAJ(1:2,3).^2));
            normTJsave(f,it) = sqrt(sum(sum(grad_AJ_f(1:2,3).^2)));
            normTJstepsave(f,it) = normTJ;
            if normTJ > 0
                normTJsquash = TJmax*normTJ/(TJmax + normTJ);
                stepAJ(1:2,3) = stepAJ(1:2,3)/normTJ*normTJsquash;
            end
            if it > start_2d_affine
                AJ(:,:,f) = AJ(:,:,f)*expm(-stepAJ);
            end
            
            % update coeffs
            coeffs{f} = (D'*bsxfun(@times,D,WM{f}(:).*WMask{f}(:)) + eye(order)*1e-3)\(D'*bsxfun(@times,reshape(J{f},numel(SAphiI{f}),size(J{f},3)),WM{f}(:).*WMask{f}(:)));
            if any(isnan(coeffs{f}(:)))
                coeffs{f} = zeros(size(coeffs{f}));
            end
            CA{f} = sum(sum(bsxfun(@times, J{f},WA{f}.*WMask{f}),1),2)/sum(WA{f}(:).*WMask{f}(:));
            if any(isnan(CA{f}))
                CA{f} = zeros(size(CA{f}));
            end
            CB{f} = sum(sum(bsxfun(@times, J{f},WB{f}.*WMask{f}),1),2)/sum(WB{f}(:).*WMask{f}(:));
            if any(isnan(CB{f}))
                CB{f} = zeros(size(CB{f}));
            end

            
            
            %             % display weights
            %             danfigure(7);
            %             imagesc(xJ{f},yJ{f},cat(3,WM{f},WA{f},WB{f}))
            %             axis image
            
            
            
            % we need to do a nice big figure with lots of slices so I can
            % see that it is actually working
            if it == 1
                ntoshow = 10;
                %                 ntoshow = 15;
                %                 ntoshow = 21;
                showlist = round(linspace(1 , length(files), ntoshow));
                tmp = 1 : length(files);
                thickones = tmp(isthicks>0);
                showlist =  unique([showlist, thickones]);
                ntoshow = length(showlist);
            end
            if any(f==showlist)
                if f == showlist(1);
                    showcount = 1;
                end
                danfigure(4554);
%                 pos = get(4554,'position');
%                 pos(3) = 960*2;
%                 set(4554,'position',pos);
                % show the target
                subplotdan(4,ntoshow,showcount);
                imagesc(xJ{f},yJ{f},J{f})
                axis image
                axis off
                % show the atlas
                subplotdan(4,ntoshow,showcount+ntoshow);
                imagesc(xJ{f},yJ{f},fSAphiI{f});
                axis image
                axis off;
                % now the error
                subplotdan(4,ntoshow,showcount+ntoshow*2);
                imagesc(xJ{f},yJ{f},err{f}/2*3 + 0.5);
                axis image
                axis off;
                % now the weight
                subplotdan(4,ntoshow,showcount+ntoshow*3);
                imagesc(xJ{f},yJ{f},cat(3,WM{f}.*WMask{f},WA{f}.*WMask{f},WB{f}.*WMask{f}));
                axis image
                axis off;
                
                showcount = showcount + 1;
                
                %                 if f == showlist(end)
                %                     drawnow;
                %                 end
            end
            
            % show the error for all
            % this is actually really slow since its so many plots
            % even making the subplots is slow
            if save_frames % only draw the figures if I'm gonna save them
                danfigure(4559);
                if it == 1
                    subplotdan(ceil(sqrt(length(files))),ceil(sqrt(length(files))),f);
                    herrAll(f) = imagesc(xJ{f},yJ{f},err{f}/2*3 + 0.5);
                    % combine these axis image and off
                    axis image off;
                else
                    set(herrAll(f),'cdata',err{f}/2*3 + 0.5);
                end
                
                danfigure(6666);
                if it == 1
                    subplotdan(ceil(sqrt(length(files))),ceil(sqrt(length(files))),f);
                    hJAll(f) = imagesc(xJ{f},yJ{f},J{f});
                    % combine these axis image and off
                    axis image off;
                end
                danfigure(6667);
                if it == 1
                    subplotdan(ceil(sqrt(length(files))),ceil(sqrt(length(files))),f);
                    hIAll(f) = imagesc(xJ{f},yJ{f},fSAphiI{f});
                    % combine these axis image and off
                    axis image off;
                else
                    set(hIAll(f),'cdata',fSAphiI{f});
                end
                
                
                % weight for all
                danfigure(4560);
                if it == 1
                    subplotdan(ceil(sqrt(length(files))),ceil(sqrt(length(files))),f);
                    hWeightAll(f) = imagesc(xJ{f},yJ{f},cat(3,WM{f}.*WMask{f},WA{f}.*WMask{f},WB{f}.*WMask{f}));
                    axis image off;
                else
                    set(hWeightAll(f),'cdata',cat(3,WM{f}.*WMask{f},WA{f}.*WMask{f},WB{f}.*WMask{f}))
                end
                
                % seg matching
                danfigure(4561);
                if edit_this_slice
                    show = (SAphiL{f} - LJ{f})/2 + 0.5;
                    if size(show,3) > 3
                        show = show(:,:,1:3);
                    elseif size(show,3) == 2
                        show = cat(3,show,show(:,:,1));
                    end
                    if it == 1
                        subplotdan(ceil(sqrt(length(files))),ceil(sqrt(length(files))),f);                        
                        hLabelAll(f) = imagesc(xJ{f},yJ{f},show,[0,1]);
                        axis image off;
                    else
                        set(hLabelAll(f),'cdata',show)
                    end
                end
            end
            
            % one example figure that I can see
            if f == 50
                danfigure(234598+1);
                imagesc(xJ{f},yJ{f},err{f}/2*3 + 0.5);
                axis image
                danfigure(234598+2);
                imagesc(xJ{f},yJ{f},fSAphiI{f});
                axis image
                danfigure(234598+3);
                imagesc(xJ{f},yJ{f},J{f});
                axis image
            end
            
            
            %%
            % now this slice's contribution of error to the 3D volume
            % d_de sum_i int |f(I(phi^{-1}(exp(-edA)A^{-1}(phiJ^{-1}(AJ^{-1}(x)))))) - J(x)|^2 dx e=0
            % let I' = f(I(phi^{-1}(x)))
            % so we get
            % d_de sum_i int |I'(exp(-edA)A^{-1}(phiJ^{-1}(AJ^{-1}(x)))) - J(x)|^2 dx e=0
            % = sum_i int (I'(A^{-1}(phiJ^{-1}(AJ^{-1}(x)))) - J(x))^T DI'(A^{-1}(phiJ^{-1}(AJ^{-1}(x)))) (-1)dA A^{-1}(phiJ^{-1}(AJ^{-1}(x)))    dx
            % this term has been calculated A^{-1}(phiJ^{-1}(AJ^{-1}(x))):
            % AiPhiJiAJiX, AiPhiJiAJiY, AiPhiJiAJiZ
            F = griddedInterpolant({yI,xI,zI},phiI_x,'linear','nearest');
            phiI_x_2d = F(AiPhiJiAJiY,AiPhiJiAJiX,AiPhiJiAJiZ);
            F = griddedInterpolant({yI,xI,zI},phiI_y,'linear','nearest');
            phiI_y_2d = F(AiPhiJiAJiY,AiPhiJiAJiX,AiPhiJiAJiZ);
            F = griddedInterpolant({yI,xI,zI},phiI_z,'linear','nearest');
            phiI_z_2d = F(AiPhiJiAJiY,AiPhiJiAJiX,AiPhiJiAJiZ);
            % note that here there is a possible issue if the image crosses
            % the boundary
            % we end up with a lot of edge artifacts on these derivatives
            % these edge artifacts may contribute to problems with the
            % affine
            % would there be a better way to do this using the derivative
            % of the original image?
            % DI' = D[I(phii(x))] = DI(phii(x))Dphii(x)
            % I could do that above
            
            if edit_this_slice                
                phiL_x_2d = zeros([size(phiI_x_2d),nLJ(f)]);
                phiL_y_2d = zeros([size(phiI_x_2d),nLJ(f)]);
                phiL_z_2d = zeros([size(phiI_x_2d),nLJ(f)]);
                for l = 1 : nLJ(f)
                    F = griddedInterpolant({yI,xI,zI},phiL_x(:,:,:,channels(l)),'linear','nearest');
                    phiL_x_2d(:,:,l) = F(AiPhiJiAJiY,AiPhiJiAJiX,AiPhiJiAJiZ);
                    F = griddedInterpolant({yI,xI,zI},phiL_y(:,:,:,channels(l)),'linear','nearest');
                    phiL_y_2d(:,:,l) = F(AiPhiJiAJiY,AiPhiJiAJiX,AiPhiJiAJiZ);
                    F = griddedInterpolant({yI,xI,zI},phiL_z(:,:,:,channels(l)),'linear','nearest');
                    phiL_z_2d(:,:,l) = F(AiPhiJiAJiY,AiPhiJiAJiX,AiPhiJiAJiZ);
                end
            end
            
            gradA_f = zeros(4);
            for r = 1 : 3
                for c = 1 : 4
                    dA = double((1:4)==r)' * double((1:4)==c);
                    
                    Xs_ = dA(1,1)*AiPhiJiAJiX + dA(1,2)*AiPhiJiAJiY + dA(1,3)*AiPhiJiAJiZ + dA(1,4);
                    Ys_ = dA(2,1)*AiPhiJiAJiX + dA(2,2)*AiPhiJiAJiY + dA(2,3)*AiPhiJiAJiZ + dA(2,4);
                    Zs_ = dA(3,1)*AiPhiJiAJiX + dA(3,2)*AiPhiJiAJiY + dA(3,3)*AiPhiJiAJiZ + dA(3,4);

                    tmp = ferrW{f}.*(phiI_x_2d.*Xs_ + phiI_y_2d.*Ys_ + phiI_z_2d.*Zs_);
                    gradA_f(r,c) = sum(sum(tmp))*prod(dxJ(f,1:2))*(-1)/sigmaM^2;
                    
                    if edit_this_slice
                        tmp = sum(errLW{f}.*phiL_x_2d,3).*Xs_ + sum(errLW{f}.*phiL_y_2d,3).*Ys_ + sum(errLW{f}.*phiL_z_2d,3).*Zs_;
                        gradA_f(r,c) = gradA_f(r,c) + sum(sum(tmp))*prod(dxJ(f,1:2))*(-1)/sigmaL^2;
                    end
                end
            end
            gradA = gradA + gradA_f;
            
            
            %%
            % now we need to find the contribution to the velocity gradient
            % from this slice
            % calculate the variation with respect to I1 = phiI
            % then I know how to map this back to I
            %
            % d_de sum_i int |f([I1 + edI1](A^{-1}(AJi^{-1}(xi)))) - J(xi)|^2 dxi e=0
            %
            % let's write the integral in a smarter way, and the 2d-3d in a
            % smarter way
            % I will do 2D 3D with delta functions, \Delta
            % they will actually be triangles with a width of 1 voxel (linear
            % interpolation)
            % the transformed image will be
            % I(x,y,z,1)
            % to
            % I(phi^{-1}(x,y,z))
            % I(phi^{-1}A^{-1}(x,y,z))
            % \int \Delta(z - zi) I(phi^{-1}A^{-1}(x,y,z)) dz
            % now we've got a 2D image that's a function of x,y only
            % then we have one more map
            % \int \Delta(z - zi) I(phi^{-1}A^{-1}(AJi^{-1}(x,y),z)) dz
            % this notation is not that nice
            % we can let AJi be from r3 to r3, and just set it to identiy for z
            % \int \Delta(z - zi) I(phi^{-1}(A^{-1}(AJi^{-1}(x,y,z)))) dz
            % remember this is NOT a funciton of z
            % so now the error is
            % \int \Delta(z - zi) I(phi^{-1}(A^{-1}(AJi^{-1}(x,y,z)))) dz - J(x,y)
            % we can of course consider J as a constant function of z
            % \int \Delta(z - zi) [I(phi^{-1}(A^{-1}(AJi^{-1}(x,y,z)))) - J(x,y,z)]dz
            % so the square error becomes
            % |\int \Delta(z - zi) [I(phi^{-1}(A^{-1}(AJi^{-1}(x,y,z))))  - J(x,y,z)]dz|^2
            % we can expand the square
            % int int  \Delta(z - zi)\Delta(z' - zi) [I(phi^{-1}(A^{-1}(AJi^{-1}(x,y,z))))  - J(x,y,z)] [I(phi^{-1}(A^{-1}(AJi^{-1}(x,y,z'))))  - J(x,y,z')] dz dz'
            % we assumew that the delta function is narrow
            % see latex figure
            % I will map the data for each slice using A and AJi
            
            % now we're going to do this by stacking my error to give proper
            % boundary conditions
            % what is the correct normalization here?
            % my image err(x,y) is multiplied by a delta function
            % the delta function should integrate to 1
            % in my case the function is just [0,1,0]
            % but it should be [0,1,0]/dz
            % I added this above
            
            % is that right? this is a triangle
            % if its width is 2dz, and its height is 1/dz,
            % then its area is 1/2 2dz 1/dz = 1
            % okay
            
            %         divide by dz to get correct normalization, see note above
            % note padarray is ending up quite slow
            errpad = zeros(nxJ{f}(2)+2,nxJ{f}(1)+2,3); % 3 is for 3 slices with appropriate zero padding, NOT RGB
            errpad(2:end-1,2:end-1,2) = ferrW{f}/dxI(3);
            xpad = [xJ{f}(1)-dxJ(1), xJ{f}, xJ{f}(end)+dxJ(1)];
            ypad = [yJ{f}(1)-dxJ(2), yJ{f}, yJ{f}(end)+dxJ(2)];
            zpad = zJ(f) + [-1,0,1]*dxI(3); % this will define a triangle function of the appropriate size for slicing
            
            if edit_this_slice
                errLpad = zeros(nxJ{f}(2)+2,nxJ{f}(1)+2,3,nLJ(f));
                errLpad(2:end-1,2:end-1,2,:) = errLW{f};
            end
            
            % now how are we gonna sample this?
            % we've got to pull back with AJ, then phiJ, then A
            % I"ve got this: AJphiJ1tinvx
            % so we start with XI
            % then map to A XI
            % then apply phiJ
            AX = A(1,1)*XI + A(1,2)*YI + A(1,3)*ZI + A(1,4);
            AY = A(2,1)*XI + A(2,2)*YI + A(2,3)*ZI + A(2,4);
            AZ = A(3,1)*XI + A(3,2)*YI + A(3,3)*ZI + A(3,4);
            % then apply phiJ, this is composition
            % we are assuming identity in Z!
            F = griddedInterpolant({yJ{f},xJ{f}},phiJ1tix-XJ{f}, 'linear','nearest');
            phiJAX = F(AY,AX) + AX;
            F = griddedInterpolant({yJ{f},xJ{f}},phiJ1tiy-YJ{f}, 'linear','nearest');
            phiJAY = F(AY,AX) + AY;
            phiJAZ = AZ;
            % then apply AJ, matrix multiplication
            % we are assuming identity in Z
            AJphiJAX = AJf(1,1)*phiJAX + AJf(1,2)*phiJAY + AJf(1,3);
            AJphiJAY = AJf(2,1)*phiJAX + AJf(2,2)*phiJAY + AJf(2,3);
            AJphiJAZ = phiJAZ;
            % and we need the detjac
            % the largest amount of time is spent calculating 3D gradients
            % it must be this!
            % there's got to be a better way, given that I already have the det
            % jacs (detjacJ and detJac) and that this is really just one slice!
            %
            %         [AJphiJAX_x,AJphiJAX_y,AJphiJAX_z] = gradient(AJphiJAX,dxI(1),dxI(2),dxI(3));
            %         [AJphiJAY_x,AJphiJAY_y,AJphiJAY_z] = gradient(AJphiJAY,dxI(1),dxI(2),dxI(3));
            %         [AJphiJAZ_x,AJphiJAZ_y,AJphiJAZ_z] = gradient(AJphiJAZ,dxI(1),dxI(2),dxI(3));
            %         detjacslice = AJphiJAX_x.*(AJphiJAY_y.*AJphiJAZ_z - AJphiJAY_z.*AJphiJAZ_y) ...
            %             - AJphiJAX_y.*(AJphiJAY_x.*AJphiJAZ_z - AJphiJAY_z.*AJphiJAZ_x) ...
            %             + AJphiJAX_z.*(AJphiJAY_x.*AJphiJAZ_y - AJphiJAY_y.*AJphiJAZ_x);
            %         F = griddedInterpolant({ypad,xpad,zpad},errpad,'linear','nearest');
            %         gradI = gradI + F(AJphiJAY,AJphiJAX,AJphiJAZ).*detjacslice;
            % the transformation I'm interested in is
            % x \mapsto AJ(phiJ(A(x)))
            % so the derivative is
            % DAJ(phiJ(A(x))) DphiJ(A(x)) DA(x)
            % AJ DphiJ(A(x)) A
            % the determinant is a product of matrices
            % |AJ| |DphiJ(A(x))| |A|
            % terms 1,3 I already have
            % the second term I need to calculate, but it should just be one
            % interpolation
            % interpolation requires at least two sample points in each
            % dimension
            if it <= start_2d_diffeo
                detjacslice = abs(det(AJ(:,:,f)))*abs(det(A));
            else
                F = griddedInterpolant({yJ{f},xJ{f},zpad},repmat(detjacJ,[1,1,3]),'linear','nearest');
                detjacslice = F(AY,AX,AZ)*abs(det(AJ(:,:,f)))*abs(det(A));
            end
            F = griddedInterpolant({ypad,xpad,zpad},errpad,'linear','nearest');
            gradI = gradI + F(AJphiJAY,AJphiJAX,AJphiJAZ).*detjacslice;
            if edit_this_slice
                for l = 1 : nLJ(f)
                    % map down to used channels          
                    F = griddedInterpolant({ypad,xpad,zpad},errLpad(:,:,:,l),'linear','nearest');
                    gradL(:,:,:,channels(l)) = gradL(:,:,:,channels(l)) + F(AJphiJAY,AJphiJAX,AJphiJAZ).*detjacslice;
                end
            end
            
        end % of file loop
        
        
        % now we start updates for the whole image
        % update A
        e = ([1;1;1;0]*[1,1,1,0]*eLI + [1;1;1;0]*[0,0,0,1]*eTI);
        if it > min(start_2d_diffeo,start_3d_diffeo)
            e = e*post_affine_reduce;
        end
        Asave = [Asave,A(:)];
        
        paramfig = 46;
        danfigure(paramfig);
        subplot(2,2,1)
        plot(Asave([1,2,3,5,6,7,9,10,11],:)')
        title('Linear')
        subplot(2,2,2)
        plot(Asave([13:15],:)')
        title('Translation')
        
        % let's normalize
        %     A = A*expm(-gradA.*e);
        stepA = gradA.*e;
        normLI = sqrt(sum(sum(stepA(1:3,1:3).^2)));
        normLIsave(it) = sqrt(sum(sum(gradA(1:3,1:3).^2)));
        normLIstepsave(it) = normLI;
        
        if normLI > 0
            normLIsquash = LImax*normLI/(LImax + normLI);
            stepA(1:3,1:3) = stepA(1:3,1:3)/normLI*normLIsquash;
        end
        normTI = sqrt(sum(stepA(1:3,4).^2));
        normTIsave(it) = sqrt(sum(sum(gradA(1:3,4).^2)));
        normTIstepsave(it) = normTI;
        if normTI > 0
            normTIsquash = TImax*normTI/(TImax + normTI);
            stepA(1:3,4) = stepA(1:3,4)/normTI*normTIsquash;
        end
        A = A*expm(-stepA);
        
        
        
        % now update deformation
        gradIpad = padarray(gradI,[1,1,1],0,'both');
        xpad = [xI(1)-dxI(1),xI,xI(end)+dxI(1)];
        ypad = [yI(1)-dxI(2),yI,yI(end)+dxI(2)];
        zpad = [zI(1)-dxI(3),zI,zI(end)+dxI(3)];
        if edit_mode
            gradLpad = padarray(gradL,[1,1,1,0],0,'both');
        end
        % start my flow
        phi1tinvx = XI;
        phi1tinvy = YI;
        phi1tinvz = ZI;
        for t = nt*(it>start_3d_diffeo) : -1 : 1
            % update diffeo (note plus)
            Xs = XI + vtx(:,:,:,t)*dt;
            Ys = YI + vty(:,:,:,t)*dt;
            Zs = ZI + vtz(:,:,:,t)*dt;
            F = griddedInterpolant({yI,xI,zI},phi1tinvx-XI,'linear','nearest');
            phi1tinvx = F(Ys,Xs,Zs) + Xs;
            F = griddedInterpolant({yI,xI,zI},phi1tinvy-YI,'linear','nearest');
            phi1tinvy = F(Ys,Xs,Zs) + Ys;
            F = griddedInterpolant({yI,xI,zI},phi1tinvz-ZI,'linear','nearest');
            phi1tinvz = F(Ys,Xs,Zs) + Zs;
            % determinant of jacobian
            [phi1tinvx_x,phi1tinvx_y,phi1tinvx_z] = gradient3d(phi1tinvx,dxI(1),dxI(2),dxI(3));
            [phi1tinvy_x,phi1tinvy_y,phi1tinvy_z] = gradient3d(phi1tinvy,dxI(1),dxI(2),dxI(3));
            [phi1tinvz_x,phi1tinvz_y,phi1tinvz_z] = gradient3d(phi1tinvz,dxI(1),dxI(2),dxI(3));
            detjac = phi1tinvx_x.*(phi1tinvy_y.*phi1tinvz_z - phi1tinvy_z.*phi1tinvz_y) ...
                - phi1tinvx_y.*(phi1tinvy_x.*phi1tinvz_z - phi1tinvy_z.*phi1tinvz_x) ...
                + phi1tinvx_z.*(phi1tinvy_x.*phi1tinvz_y - phi1tinvy_y.*phi1tinvz_x);
            % resample
            F = griddedInterpolant({ypad,xpad,zpad},gradIpad,'linear','nearest');
            lambdat = F(phi1tinvy,phi1tinvx,phi1tinvz).*detjac*(-1)/sigmaM^2;
            
            % gradient of It
            [It_x,It_y,It_z] = gradient3d(It(:,:,:,t),dxI(1),dxI(2),dxI(3));
            
            % matching function gradient
            gradx = It_x.*lambdat;
            grady = It_y.*lambdat;
            gradz = It_z.*lambdat;
            
            if edit_mode
                lambdat = zeros([size(lambdat),nL]);
                for l = 1 : nL
                    F = griddedInterpolant({ypad,xpad,zpad},gradLpad(:,:,:,l),'linear','nearest');
                    lambdat(:,:,:,l) = F(phi1tinvy,phi1tinvx,phi1tinvz).*detjac*(-1)/sigmaL^2;
                end
                
                % gradient of It
                Lt_x = zeros([size(lambdat)]);
                Lt_y = zeros([size(lambdat)]);
                Lt_z = zeros([size(lambdat)]);
                for l = 1 : nLJ(f)
                    [Lt_x(:,:,:,l),Lt_y(:,:,:,l),Lt_z(:,:,l)] = gradient3d(lambdat(:,:,:,l),dxI(1),dxI(2),dxI(3));
                end
                % matching function gradient
                gradx = gradx + sum(Lt_x.*lambdat,4);
                grady = grady + sum(Lt_y.*lambdat,4);
                gradz = gradz + sum(Lt_z.*lambdat,4);
            end
            
            % regularization ( I could add more smoothness here if I wanted)
            gradx = ifftn(   Kp.*(fftn(gradx).*K + vtxhat(:,:,:,t)/sigmaR^2)   ,'symmetric');
            grady = ifftn(   Kp.*(fftn(grady).*K + vtyhat(:,:,:,t)/sigmaR^2)   ,'symmetric');
            gradz = ifftn(   Kp.*(fftn(gradz).*K + vtzhat(:,:,:,t)/sigmaR^2)   ,'symmetric');
            
            stepx = eVI*gradx;
            stepy = eVI*grady;
            stepz = eVI*gradz;
            
            % squash it
            % maybe the squashing is really upping the reg energy, not sure
            normstep = sqrt(stepx.^2 + stepy.^2 + stepz.^2);
            squashnorm = VImax*normstep./(VImax + normstep);
            stepx = stepx./normstep.*squashnorm;
            stepy = stepy./normstep.*squashnorm;
            stepz = stepz./normstep.*squashnorm;
            
            vtx(:,:,:,t) = vtx(:,:,:,t) - stepx;
            vty(:,:,:,t) = vty(:,:,:,t) - stepy;
            vtz(:,:,:,t) = vtz(:,:,:,t) - stepz;
        end
        
        
        
        EM = sum(EMall);
        EMsave(it) = EM;
        EMJsave(:,it) = EMall;
        E = ER + EM + sum(ERJ);
        Esave(it) = E;
        subplot(2,2,3)
        plot([EMsave(1:it);ERsave(1:it);sum(ERJsave(:,1:it),1);Esave(1:it)]');

        saveas(paramfig,[prefix 'energy.png'])
        disp(['Iter: ' num2str(it) ', energy: ' num2str(E) ', matching energy: ' num2str(EM) ', reg energy: ' num2str(ER) ', regJ energy: ' num2str(sum(ERJ))] )
        
        %%
        danfigure(323);
        sliceView(xI,yI,zI,gradI);
        
        % a figure for norm of gradient
        danfigure(113);
        subplot(2,2,1)
        imagesc(normLJsave);colorbar
        title('norm grad LJ')
        subplot(2,2,2)
        imagesc(normTJsave);colorbar
        title('norm grad of TJ')
        subplot(2,2,3)
        plot(normLIsave)
        title('norm grad of LI')
        subplot(2,2,4)
        plot(normTIsave)
        title('norm grad of TI')
        
        
        danfigure(114);
        subplot(2,2,1)
        imagesc(normLJstepsave);colorbar
        title('norm grad step LJ')
        subplot(2,2,2)
        imagesc(normTJstepsave);colorbar
        title('norm grad of step TJ')
        subplot(2,2,3)
        plot(normLIstepsave)
        title('norm grad of step LI')
        subplot(2,2,4)
        plot(normTIstepsave)
        title('norm grad of step TI')
        
        
        % save some frames and write gifs
        if it == 1
            frameSlices = [];
            frameErr = [];
        end
        if save_frames
            frameSlices = [frameSlices,getframe(4554)];
            frame2Gif(frameSlices,[prefix 'slices.gif'])
            frameErr = [frameErr,getframe(323)];
            frame2Gif(frameErr,[prefix 'error3d.gif']);
            framePhiI = [framePhiI,getframe(333)];
            frame2Gif(framePhiI,[prefix 'phiI.gif']);
            frameErrAll = [frameErrAll,getframe(4559)];
            frame2Gif(frameErrAll,[prefix 'errAll.gif']);
            frameWeightAll = [frameWeightAll,getframe(4560)];
            frame2Gif(frameWeightAll,[prefix 'weightAll.gif']);
            frameIAll = [frameIAll, getframe(6667)];
            frame2Gif(frameIAll,[prefix 'IAll.gif'])
            if edit_mode
                frameLabelAll = [frameLabelAll,getframe(4561)];
                frame2Gif(frameLabelAll,[prefix 'LAll.gif'])
            end
        end
        
        
        
        
        
        % apply this constraint, average should be 0 (expressed in the gobal affine)
        logAJ = zeros(size(AJ));
        for f = 1 : length(files)
            logAJ(:,:,f) = real(logm(AJ(:,:,f)));
        end
        meanLogAJ = mean(logAJ,3);
%         logAJ = bsxfun(@minus,logAJ,meanLogAJ);

        % maybe what would be better is if they are not zero mean subtract the
        % mean from them and add it to A
%         for f = 1 : length(files)
%             AJ(:,:,f) = expm(logAJ(:,:,f));
%         end
        
        % it looks like all the "low contrast" slices are moving up for no
        % apparent reason
        % could it be that that are sliding up because the others are being
        % pushed down?  I believe the below approach can fix this
        
        
        % this is the approach used in initialization
        themeanAJ = expm(meanLogAJ);
        % we want to "subtract" the mean from AJ and "add" it to A
        for i = 1 : length(files)
            AJ(:,:,i) = AJ(:,:,i)/themeanAJ; % inverse on the right
        end
        A = [themeanAJ(1,1:2),0,themeanAJ(1,3);
            themeanAJ(2,1:2),0,themeanAJ(2,3);
            0,0,1,0;
            0,0,0,1]*A; % forward on the left

            
        % note that there's another issue here
        % which is that if the first moment of these translations in Z is not zero, it could
        % be modelled as an affine in A
        % so I still am overparameterized
        % so maybe I should subtract the first moment also
        
        save([prefix 'A.mat'],'AJ','A')
        
        drawnow;
        
        if ~mod(it-1,100) 
            save([prefix 'v.mat'],'vtx','vty','vtz', 'vJtx', 'vJty','xJ','yJ','zJ','xI','yI','zI')
            save([prefix 'coeffs.mat'],'coeffs','CA','CB')
        end
        
        
    end
    save([prefix 'v.mat'],'vtx','vty','vtz', 'vJtx', 'vJty','xJ','yJ','zJ','xI','yI','zI')
    save([prefix 'coeffs.mat'],'coeffs','CA','CB')
    
    
    toc

    
    % in example, quit before high res.
    if ~isempty(strfind(target_dir,'example')) && downloop == 1
        break        
    end
    
end % of downloop