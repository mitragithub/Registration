function A = slice_to_slice_rigid_alignment_GN_weight(xI,yI,I,xJ,yJ,J,A0,downs,niter,e)
% rigidly align two slices, possibly different modalities
% we wish to align image I to image J rigidly
% we use exponential parameterization
% minimize
% int 1/2 |f(I(Ai x)) - J(x))|^2 dx
% let fI = I'
% int 1/2 |I'(Ai x) - J(x))|^2 dx
% right perturbation
% A \to A exp(e dA)
% or Ai \to  exp(-edA)Ai
% d_ee int 1/2 |I'(exp(-edA) Ai x) - J(x))|^2 dx |_e=0
% = \int (I'(Ai x) - J(x))^T DI'(Ai x)  (-1) dA Ai x dx
%= (-1) \int (I'(Ai x) - J(x))^T D[I'(Ai x)] A  dA Ai x dx
%
% now for the gauss newton version
% we will optimize over the matrix B=Ai for simplicity
% we will not use group perturbation, instead, just project onto rigid
% d_ee int 1/2 |I'(B x + edB x) - J(x))|^2 dx |_e=0
% = int (I'(Bx) - J(x)) DI'(Bx) dB x dx
% = int (I'(Bx) - J(x)) D[I'(Bx)] Bi dB x dx
% 
% need to incorporate weight (done)
% in this code we set artifact channel value to 1 and do not update
%
if nargin == 0
    % run an example
    
    I = imread('/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/Xu2Daniel/PMD964/PMD965&964-N33-2012.06.21-00.30.07_PMD964_3_0099.tif');
    I = double(I)/255.0;
    dxI = [14.72 14.72];
    xI = (0:size(I,2)-1)*dxI(1);xI = xI - mean(xI);
    yI = (0:size(I,1)-1)*dxI(2);yI = yI - mean(yI);
    
    J = imread('/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/Xu2Daniel/PMD964/PMD965&964-F33-2012.06.17-14.22.42_PMD964_3_0099.tif');
    J = double(J)/255.0;
    dxJ = [14.72 14.72];
    xJ = (0:size(J,2)-1)*dxJ(1);xJ = xJ - mean(xJ);
    yJ = (0:size(J,1)-1)*dxJ(2);yJ = yJ - mean(yJ);
    
    order = 2;
    keyboard
end




addpath Functions/plotting
addpath Functions/gradient
addpath Functions/downsample

% I = mean(I,3);
% J = mean(J,3);
qlim = [0.01,0.99];
climI = quantile(I(:),qlim);
climJ = quantile(J(:),qlim);
danfigure(1);
imagesc(xI,yI,I)
axis image

danfigure(2);
imagesc(xJ,yJ,(J-climJ(1))/diff(climJ))
axis image

CI = size(I,3);
CJ = size(J,3);

[XJ,YJ] = meshgrid(xJ,yJ);
dxI = [xI(2) - xI(1), yI(2) - yI(1)];
dxJ = [xJ(2) - xJ(1), yJ(2) - yJ(1)];

if nargin < 7
    A0 = eye(3);
end
if nargin < 8
    downs = [8,4,2];
end
if nargin < 9
    niter = 100;
end
if length(niter) == 1
    niter = repmat(niter,length(downs));
end
if nargin < 10
    e = 0.5;
end


sigmaM = 0.05;
sigmaA = sigmaM*5;
muA = cat(3,1,1,1);


%%
A = A0;
Ai = inv(A);
ndraw = 21;
for downloop = 1 : length(downs)
    d = downs(downloop);

    
    
    Id = [];
    for c = 1 : CI
        [~,~,Id(:,:,c)] = downsample2D(1:size(I,2),1:size(I,1),I(:,:,c),[1,1]*d);
    end
    dxId = dxI*d;
    xId = (1 : size(Id,2))*dxId(1); xId = xId - mean(xId);
    yId = (1 : size(Id,1))*dxId(2); yId = yId - mean(yId);
    Jd = [];
    for c = 1 : CJ
        [~,~,Jd(:,:,c)] = downsample2D(1:size(J,2),1:size(J,1),J(:,:,c),[1,1]*d);
    end
    dxJd = dxJ*d;
    xJd = (1 : size(Jd,2))*dxJd(1); xJd = xJd - mean(xJd);
    yJd = (1 : size(Jd,1))*dxJd(2); yJd = yJd - mean(yJd);
    
    WM = ones(size(Jd(:,:,1)))*0.9;
    WA = ones(size(Jd(:,:,1)))*0.1;
    
    %
    [XJd,YJd] = meshgrid(xJd,yJd);
    for it = 1 : niter(downloop)
        
        % transform I
        AId = zeros(size(Jd,1),size(Jd,2),CI);
        Xs = Ai(1,1)*XJd + Ai(1,2)*YJd + Ai(1,3);
        Ys = Ai(2,1)*XJd + Ai(2,2)*YJd + Ai(2,3);
        for c = 1 : CI
            F = griddedInterpolant({yId,xId},Id(:,:,c),'linear','nearest');
            AId(:,:,c) = F(Ys,Xs);
        end
        
        
        
        
        % now intensity
        % start with just linear
        D = cat(3,ones(size(AId,1),size(AId,2)),AId);
        Dvec = reshape(D,[],size(D,3));
        coeffs = (bsxfun(@times,Dvec,WM(:))'*Dvec) \ (bsxfun(@times,Dvec,WM(:))' * reshape(Jd,[],CJ));

        if any(isnan(coeffs))
            D = cat(3,ones(size(AId,1),size(AId,2)));
            Dvec = reshape(D,[],size(D,3));
            coeffs = (Dvec'*Dvec) \ (Dvec' * reshape(Jd,[],CJ));
        end
        
        fAId = reshape(Dvec*coeffs,size(Jd));
        
        
        
        % now get the cost
        err = fAId - Jd;                                                
        E = sum(sum(sum(err.^2,3).*WM))/2/sigmaM^2*prod(dxJd);
        
        % update weights
        WM = exp(-sum(err.^2,3)/2/sigmaM^2)/sqrt(2*pi*sigmaM^2)^size(Jd,3);
        WA = exp(-sum((bsxfun(@minus,Jd, muA)).^2,3)/2/sigmaA^2)/sqrt(2*pi*sigmaA^2)^size(Jd,3);
        WSum = WM+WA;
        WM = WM./WSum;
        WA = WA./WSum;
        WM(WSum==0) = 0;
        WA(WSum==0) = 0;
        
        
        
        if ~mod(it-1,ndraw) || it == niter(downloop)
            danfigure(1);
        imagesc(xId,yId,Id)
        axis image
        title('I')
        
        danfigure(2);
        imagesc(xJd,yJd,(Jd-climJ(1))/diff(climJ))
        axis image
        title('J')
        
        danfigure(3);
        imagesc(xJd,yJd,AId);
        axis image
        title('AI')
        
        danfigure(4);
        imagesc(xJd,yJd,(fAId-climJ(1))/diff(climJ));
        axis image
        title('fAI')
        danfigure(5);
        imagesc(xJd,yJd,err/2.0/diff(climJ) + 0.5)
        axis image
        title('err')
        disp(['Iteration ' num2str(it) '/' num2str(niter(downloop)) ', energy ' num2str(E)]);
        
        end
        

        % and the gradient
        fAId_x = zeros(size(fAId));
        fAId_y = zeros(size(fAId));
        for c = 1 : CJ
            [fAId_x(:,:,c),fAId_y(:,:,c)] = gradient2d(fAId(:,:,c),dxJd(1),dxJd(2));
        end

        
        
        Jerr = zeros(size(Jd,1),size(Jd,2),size(Jd,3),6);
        count = 0;
        for r = 1 : 2
            for c = 1 : 3
                dA = double((1:3==r))' * double((1:3==c));
                AdAAi = A*dA;
                Xs = AdAAi(1,1)*XJd + AdAAi(1,2)*YJd + AdAAi(1,3);
                Ys = AdAAi(2,1)*XJd + AdAAi(2,2)*YJd + AdAAi(2,3);
                count = count + 1;
                Jerr(:,:,:,count) = (bsxfun(@times, fAId_x,Xs) + bsxfun(@times, fAId_y,Ys));
            end
        end
        JerrJerr = squeeze(sum(sum(sum(bsxfun(@times, permute(bsxfun(@times,Jerr,WM),[1,2,3,4,5]) , permute(Jerr,[1,2,3,5,4])),3),2),1));
        
        
        % update
        step = JerrJerr \ squeeze(sum(sum(sum(bsxfun(@times, bsxfun(@times,Jerr,WM), err),3),2),1));
        step = reshape(step,3,2)';
        if any(isnan(step))
            step = zeros(size(step));
        end
        Ai(1:2,1:3) = Ai(1:2,1:3) - e * step;
        
        % make it rigid
        [U,S,V] = svd(Ai(1:2,1:2));
        Ai(1:2,1:2) = U*V';
        A = inv(Ai);
        
        drawnow
        
        
    end % of iter
    
end % of downsample
% output
A = inv(Ai);