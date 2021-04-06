function A = ThreeD_to_3D_affine_registration_GN(xI,yI,zI,I,xJ,yJ,zJ,J,A0,downs,niter,e)
% this function will register two volumes based on affine
% the idea is to construct a volume from the aligned slices and do an
% initial affine map onto it
%
% affine registration between I and J
% run at several different stages of downsampling
%
% Now using Gauss newton optimization which is much more robust than
% gradient descent
%
% we will consider vector perturbations, not group
% and we will optimize over Ai (inverse)
% d_de 1/2 int (I(Aix + edAix) - J(x))^2 dx |_{e=0}
% = int (I(Ai x) - J(x)) DI(Ai x) dAi x dx
% = int (I(Ai x) - J(x)) D[I(Ai x)] A dAi x dx
% this factors into 
%
% keyboard
addpath Functions/downsample
addpath Functions/plotting
addpath Functions/gradient
if nargin == 0
    addpath /cis/home/dtward/Functions/avwQuiet
    [avw,nx,dx,xI,yI,zI] = avw_img_read_domain(['/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/STPT_data/STP_feedback_Brains/180722_50']);
    I = double(avw.img);
    xI = xI - mean(xI);
    yI = yI - mean(yI);
    zI = zI - mean(zI);

    [avw,nx,dx,xJ,yJ,zJ] = avw_img_read_domain(['/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/STPT_data/STP_feedback_Brains/190123_50']);
    J = double(avw.img);
    xJ = xJ - mean(xJ);
    yJ = yJ - mean(yJ);
    zJ = zJ - mean(zJ);
    
    downs = [4,2];
    niter = 10;
    e = 0.1;
end
if nargin < 9
    A0 = eye(4);
end
if nargin < 12
    e = 0.5;
end
if length(niter) == 1
    niter = ones(size(downs))*niter;
end


I = (I - mean(I(:))) / std(I(:));
J = (J - mean(J(:))) / std(J(:));

% get ready for downsampling
I0 = I;
xI0 = xI;
yI0 = yI;
zI0 = zI;

J0 = J;
xJ0 = xJ;
yJ0 = yJ;
zJ0 = zJ;


% initialize
A = A0;
Ai = inv(A);

E = zeros(1,sum(niter));
% loop over downsamples
itercount = 0;
for downloop = 1 : length(downs)
    down = downs(downloop);
    
    % downsample the images
    [xI,yI,zI,I] = downsample(xI0,yI0,zI0,I0,[1,1,1]*down);
    dxI = [xI(2)-xI(1), yI(2)-yI(1), zI(2)-zI(1)];
    
    danfigure(1);
    sliceView(xI,yI,zI,I);
    title('atlas')
    
    [xJ,yJ,zJ,J] = downsample(xJ0,yJ0,zJ0,J0,[1,1,1]*down);
    dxJ = [xJ(2)-xJ(1), yJ(2)-yJ(1), zJ(2)-zJ(1)];    
    [XJ,YJ,ZJ] = meshgrid(xJ,yJ,zJ);
    danfigure(2);
    sliceView(xJ,yJ,zJ,J);
    title('target')
    
    % start gradient descent loop
    
    for it = 1 : niter
        itercount = itercount + 1;
        % transform I
        Xs = Ai(1,1)*XJ + Ai(1,2)*YJ + Ai(1,3)*ZJ + Ai(1,4);
        Ys = Ai(2,1)*XJ + Ai(2,2)*YJ + Ai(2,3)*ZJ + Ai(2,4);
        Zs = Ai(3,1)*XJ + Ai(3,2)*YJ + Ai(3,3)*ZJ + Ai(3,4);
        F = griddedInterpolant({yI,xI,zI},I,'linear','nearest');
        AI = F(Ys,Xs,Zs);
        
        
        danfigure(3);
        sliceView(xJ,yJ,zJ,AI)
        
        % linear prediction
        COV = cov(AI(:),J(:));
        fAI = (AI - mean(AI(:)))/COV(1,1)*COV(1,2) + mean(J(:));
        
        % error
        err = (fAI - J);
        danfigure(4);
        sliceView(xJ,yJ,zJ,err)
        
        % cost
        E(itercount) = sum(abs(err(:)).^2)/2*prod(dxJ);
        danfigure(5);
        plot(E(1:itercount))
        disp(['Down loop ' num2str(downloop) ', iteration '  num2str(it) ', energy ' num2str(E(itercount))])
        
        % now we take gradient
        % to do this use right perturbation
        % int (I(Ai x) - J(x))^2 dx
        % d_de int (I(exp(-e dA )Ai x) - J(x))^2 dx e = 0
        % = int (I(Ai x) - J(x)) DI (Ai x) (-1) dA Ai x dx
        % = int (I(Ai x) - J(x)) D[I (Ai x)] A (-1) dA Ai x dx
        
        [fAI_x,fAI_y,fAI_z] = gradient3d(fAI,dxJ(1),dxJ(2),dxJ(3));
        
        Jerr = zeros(size(J,1),size(J,2),size(J,3),12);
        count = 0;
        for r = 1 : 3
            for c = 1 : 4
                dA = double((1:4==r))' * double((1:4==c));
                AdAAi = A*dA;
                Xs = AdAAi(1,1)*XJ + AdAAi(1,2)*YJ + AdAAi(1,3)*ZJ + AdAAi(1,4);
                Ys = AdAAi(2,1)*XJ + AdAAi(2,2)*YJ + AdAAi(2,3)*ZJ + AdAAi(2,4);
                Zs = AdAAi(3,1)*XJ + AdAAi(3,2)*YJ + AdAAi(3,3)*ZJ + AdAAi(3,4);
                count = count + 1;
                Jerr(:,:,:,count) = (bsxfun(@times, fAI_x,Xs) + bsxfun(@times, fAI_y,Ys) + bsxfun(@times, fAI_z,Zs));
            end
        end
        JerrJerr = squeeze(sum(sum(sum(bsxfun(@times, permute(Jerr,[1,2,3,4,5]) , permute(Jerr,[1,2,3,5,4])),3),2),1));
        
        % step
        step = JerrJerr \ squeeze(sum(sum(sum(bsxfun(@times, Jerr, err),3),2),1));
        step = reshape(step,4,3)';
        Ai(1:3,1:4) = Ai(1:3,1:4) - e * step;
        A = inv(Ai);
        
        

        drawnow;
        
    end % end of gradient descent iteration loop it
    
    
    
    
end % of downsampling loop (downloop)

A = inv(Ai);