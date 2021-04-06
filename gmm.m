% gaussian mixture model

function [L,mu] = gmm(I,OPT)
if nargin == 1
    OPT = struct;
end
%default options
image = 0; % show an image, dimensions are N1xN2xC channels, versus NxC
niter = 10;
M = 2; % components
draw = 1;

% update options
if isfield(OPT,'image')
    image = OPT.image;
end
if isfield(OPT,'niter')
    niter = OPT.niter;
end
if isfield(OPT,'M')
    M = OPT.M;
end
if isfield(OPT,'draw')
    draw = OPT.draw;
end


if image
    imagesize = size(I);
    N = imagesize(1)*imagesize(2);
    C = size(I,3);
    I = reshape(I,N,C);
else
    N = size(I,1);
    C = size(I,2);
end






% initialize parameters
% I will rearrange my data as follows
% first axis realizations
% second axis vector components
% third axis mixture component
% or for covariance, the first and second are vector components

p = ones(1,1,M); p = p / sum(p);
mu = zeros(1,C,M);
Sigma = zeros(1,C,C,M);
for c = 1 : C
    qs = linspace(min(I(:,c)),max(I(:,c)),M+2);
    mu(1,c,:) = qs(2:end-1);
end
if isfield(OPT,'mu')
    mu = OPT.mu;
%     draw = 1;
%     keyboard
end

for m = 1 : M
    Sigma(1,:,:,m) = eye(C)*0.1;
end
if isfield(OPT,'Sigma')
    Sigma = OPT.Sigma;
end






if image
    if draw
%         h = danfigure();
        h = danfigure(37);
        subplot(1,2,1)
        imagesc(reshape(I,imagesize))
        axis image
    end
end

for it = 1 : niter
    % likelihood
    % first invert sigma
    Sigmai = zeros(size(Sigma));
    detSigma = zeros(1,1,M);
    for m = 1 : M
        Sigmai(1,:,:,m) = inv(squeeze(Sigma(1,:,:,m)));
        detSigma(m) = det(squeeze(Sigma(1,:,:,m)));
    end
    
    I0 = bsxfun(@minus,I,mu);
    SigmaiI0 = squeeze(sum(  bsxfun(@times, Sigmai, reshape(I0,N,1,C,M) ), 3 ));
    I0SigmaiI0 = sum(bsxfun(@times, I0, SigmaiI0),2);
    
    % compute conditional likelihood
    l = bsxfun(@rdivide,bsxfun(@times,exp(-I0SigmaiI0/2),p*0+1),sqrt(2*pi*detSigma).^C); % do I need a p in here?
    
    % update probs
    psample = bsxfun(@rdivide, l, sum(l,3));
    n = sum(psample,1); % effective sample size
    p = bsxfun(@rdivide,n,sum(n(:)));
    
    % update params, weight by psample
    mu = bsxfun(@rdivide, sum(bsxfun(@times,I,psample),1) , n);
    I0 = bsxfun(@minus,I,mu);
    Sigma = bsxfun(@rdivide, sum(bsxfun(@times, bsxfun(@times, reshape(I0,N,C,1,M), reshape(I0,N,1,C,M)  ), reshape(psample,N,1,1,M)), 1), reshape(n,1,1,1,M));
    
    
    if image
        if draw
            danfigure(h);
            subplot(1,2,2);
            if M == 2
                imagesc(reshape(cat(3,psample(:,:,1),psample(:,:,2),psample(:,:,1)),imagesize(1),imagesize(2),3))
            else
                imagesc(reshape(psample(:,:,1:3),imagesize(1),imagesize(2),3))
            end
            axis image
            drawnow
        end
    end
    
    
end

L = psample;
if image
    L = reshape(psample,imagesize(1),imagesize(2),M);
end