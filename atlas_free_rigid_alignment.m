function atlas_free_rigid_alignment(target_dir, pattern, output_dir, r, downs, niter, et_factor, etheta_factor, skip_thick, load_initializer)

% atlas free slice alignment 
% rigid alignment of a stack
% each slice is iteratively aligned to a weighted average of its neighbors
% this can be used as a standalone program
% or as an initial guess atlas based slice alignment
%
% TODO, deal with slice spacing in physical units
% TODO, deal with slices left out
% skip_thick is a number that we don't update if it is bigger than this

addpath Functions/plotting
addpath Functions/downsample

if nargin < 4
    r = 25;
end
if nargin < 5
    downs = [32,16,8,4];
end
if nargin < 6
    niter = 200;
end
if nargin < 8
    % step sizes
    etheta_factor = 5e-11;
    et_factor = 5e-4;
    
    % on july 17, half the step size and double the iters
    etheta_factor = 2e-11;
    et_factor = 2e-4;
    etheta_factor = 1e-11;
    et_factor = 1e-4;


end

if nargin < 9
    skip_thick = -1;
end
if nargin < 10
    load_initializer=0;
end







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
dz0 = cellfun(@(x)str2num(x), csv_data(:,7));









n = length(files);


% reconstruction params
theta = zeros(1,n) ;
tx = zeros(1,n);
ty = zeros(1,n);

% AJ is the transform from RECONSTRUCTED to OBSERVED
% images deform with inverse, so we will use this to transform OBSERVED to
% RECONSTRUCTED
AJ = zeros(3,3,length(files));
% instead of zero let's initialize it
if exist([output_dir 'initializer_A.mat'],'file') && load_initializer
    vars = load([output_dir 'initializer_A.mat']);
    AJ = vars.AJ;

    % get angle from x,y with correct sign
    theta = squeeze(atan2(AJ(1,2,:),AJ(1,1,:)))';
    tx = squeeze(AJ(1,3,:))';
    ty = squeeze(AJ(2,3,:))';
end


Cost = zeros(1,10000);

itercount = 0;


for downcount = 1 : length(downs)
    down = downs(downcount);
    
    % step sizes
    % note about scaling
    % if I downsample by 2, we expect the changes to be dominated by edges
    % whose gradient will get two times smaller
    % but
    % we also must consider the value of the error at the edge
    % due to blurring, this may also get two times smaller
    % so maybe there should be a factor of down squared?
    etheta = etheta_factor*down;
    et = et_factor*down;
    
    
    % first pass, loop through files
    % get max size
    nx = [0,0];
    for i = 1 : length(files)
        I_ = imread([target_dir files{i}]);
        I_ = double(I_)/255.0;
        I_ = mean(I_,3);
        
        
        
        % for padding
        val(i) = mode(I_(I_(:)~=1));
        if isnan(val(i))
            val(i) = val(i-1); % entirely missing slice?
        end
        I_(I_==1) = val(i);
        
        
        % for thick images, just ignore them
        if dz0(i) > skip_thick
            if i > 1
                I_ = ones(size(I_))*val(i-1);
            else
                I_ = zeros(size(I_));
            end
        end
        
        
        [~,~,I_] = downsample2D(1:size(I_,2),1:size(I_,1),I_,down*[1 1]);
        
        Iload{i} = I_;
        nx = max(nx,[size(I_,2),size(I_,1)]);
    end
    nx = [nx,length(files)];
    dx = [down*dxJ0(1,1),down*dxJ0(1,2),1/n*nx(1)*dxJ0(1,1)*down]; 
    % the dz component doesn't matter, just use something to display as a square
    
    x = (0 : nx(1)-1)*dx(1);x = x - mean(x);
    y = (0 : nx(2)-1)*dx(2);y = y - mean(y);
    z = (0 : nx(3)-1)*dx(3);z = z - mean(z);
    
    [XS,YS] = meshgrid(x,y);
    
    % load data again and pad to a fixed size
    % we will need to use zero centered coordinates
    I = zeros(nx(2),nx(1),nx(3));
    for i = 1 : length(files)
        tmp = I(1:size(Iload{i},1),1:size(Iload{i},2),i);
        tmp(tmp==1) = val(i);
        I(:,:,i) = val(i);
        offr = round(size(I,1)/2-size(tmp,1)/2);
        offc = round(size(I,2)/2-size(tmp,2)/2);
        I((1:size(Iload{i},1))+offr, (1:size(Iload{i},2))+offc, i) = Iload{i};
    end
    I = (I - mean(I(:)))/std(I(:));
    qlim = [0.001,0.999];
    clim = quantile(I(:),qlim);
    
    %%
    % weights for averaging
    tmp = exp(-[-r:r].^2/2/(r/3)^2);
    tmp = tmp / sum(tmp);
    CC = ([-r:r]==0) - tmp;
    CC = CC*0.2;
    nmax = length(CC);
    %%
    
    [I_x,I_y,I_z] = gradient(I,dx(1),dx(2),1);
    for iter = 1 : niter
        I_R = zeros(size(I));
        I_x_R = zeros(size(I));
        I_y_R = zeros(size(I));
        % deform image and its derivatives
        for i = 1 : nx(3)
            Xs = cos(theta(i))*XS - sin(theta(i))*YS + tx(i);
            Ys = sin(theta(i))*XS + cos(theta(i))*YS + ty(i);
            F = griddedInterpolant({y,x},I(:,:,i),'linear','nearest');
            I_R(:,:,i) = F(Ys,Xs);
            F = griddedInterpolant({y,x},I_x(:,:,i),'linear','nearest');
            I_x_R(:,:,i) = F(Ys,Xs);
            F = griddedInterpolant({y,x},I_y(:,:,i),'linear','nearest');
            I_y_R(:,:,i) = F(Ys,Xs);            
        end
        
        danfigure(1);
        sliceView(x,y,z,I_R,7,clim)
        
        
        
        % apply A to R
        % 0 pad
        off = (nmax-1)/2;
        I_R_pad = padarray(I_R,[0,0,off],'symmetric','both');
        AI_R_pad = zeros(size(I_R_pad));
        for i = 1 : length(CC)
            %         i-off-1
            AI_R_pad = AI_R_pad + CC(i)*circshift(I_R_pad,[0,0,i-off-1]);
        end
        
        AI_R = AI_R_pad(:,:,1+off:end-off);
        danfigure(2);
        sliceView(x,y,z,AI_R);
        
        
        % get the cost
        E = sum(AI_R(:).*I_R(:)/2*prod(dx(1:3)));
        disp(['Down ' num2str(down) ', it ' num2str(iter) ', E ' num2str(E) ])
        itercount = itercount + 1;
        Cost(itercount) = E;
        danfigure(33);
        plot(Cost(1:itercount))
        
        % now calculate gradient on each slice
        thetagrad = zeros(size(theta));
        txgrad = zeros(size(tx));
        tygrad = zeros(size(ty));
        for i = 1 : nx(3)
            if dz0(i) > skip_thick
                continue
            end
            thetagrad(i) = sum(sum(AI_R(:,:,i) .* (I_x_R(:,:,i).*(-sin(theta(i))*XS - cos(theta(i))*YS) + I_y_R(:,:,i).*(cos(theta(i))*XS + -sin(theta(i))*YS)   ) )) * prod(dx(1:2));
            txgrad(i) = sum(sum(AI_R(:,:,i).*I_x_R(:,:,i)))*prod(dx(1:2));
            tygrad(i) = sum(sum(AI_R(:,:,i).*I_y_R(:,:,i)))*prod(dx(1:2));
        end
        theta = theta - etheta * thetagrad;
        tx = tx - et*txgrad;
        ty = ty - et*tygrad;
        
        % make them zero maen
        tx = tx - mean(tx);
        ty = ty - mean(ty);
        theta = theta - mean(theta);
        
%         disp(tx(14));disp(ty(14));disp(theta(14));
        
        
        drawnow;
    end % of loop over iterations
    
    
end % end of down loop


% save data in required format
for i = 1 : n
    % note I apply transforms to each slice like this
%     Xs = cos(theta(i))*XS - sin(theta(i))*YS + tx(i);
%     Ys = sin(theta(i))*XS + cos(theta(i))*YS + ty(i);
%     F = griddedInterpolant({y,x},I(:,:,i),'linear','nearest');

    A = [cos(theta(i)),-sin(theta(i)),tx(i);
        sin(theta(i)),cos(theta(i)),ty(i);
        0,0,1];

    AJ(:,:,i) = A;
end
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end
save([output_dir 'initializer_A.mat'],'AJ');


return
keyboard


%%

for i = 1 : n
    danfigure(44);
    imagesc(I_R(:,:,i))
    axis image
    title(['slice ' num2str(i) ' of ' num2str(n)])
    drawnow
%     pause
end