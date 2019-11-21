function mask = getmask(I,flag)
% segment image into 3 classes
% then

if nargin == 1
    flag = 0; % dark background for fluoro, else bright background for nissl
end

% flag = 1 for bright background
OPT = struct;
OPT.image = 1;
OPT.M = 3;
if flag == 0 % fluoro
    OPT.mu = cat(3,[0,0,0],[0.15,0.15,0],[1,1,1]);
    OPT.Sigma = reshape(cat(3,eye(3)*0.5,eye(3)*0.5,eye(3)*2),1,3,3,3);
    OPT.draw = 36;
    W = 0;
end
if flag == 1 % nissl
    OPT.mu = cat(3,[0.7,0.7,0.8],[0.9,0.8,0.9],[1,1,1]);
    OPT.draw = 37;
    
    W = detect_padding_in_nissl(I);
end
% [L,mu] = gmm(I,OPT);
[L,mu] = kmeans(I,OPT);
L = bsxfun(@times, L, (1-W));
% we need to find the two largest labels, and assign one background and
% one forground
n = squeeze(sum(sum(L,1),2));
[n,perm] = sort(n);
mu = mu(:,:,perm);
L = L(:,:,perm);
mu = squeeze(mean(mu,2));
% consider only biggest 2
mu = mu(end-1:end);
L = L(:,:,end-1:end);
if flag == 1
    % out of the two biggest, which is bg?  It is the BRIGHTER
    ind = find(mu == max(mu));
    if ind == length(mu)
        fgind = length(mu)-1;
    elseif ind == length(mu)-1
        fgind = length(mu);
    else
        warning('Could not identify components')
        fgind = length(mu);
    end
else
    ind = find(mu == min(mu));
    if ind == length(mu)
        fgind = length(mu)-1;
    elseif ind == length(mu)-1
        fgind = length(mu);
    else
        warning('Could not identify components')
        fgind = length(mu);
    end
end



% mask = (max(L,[],3) == L(:,:,fgind)) .* (1-W); % for em
mask = L(:,:,fgind); % for kmeans
try
mask = bwareafilt(mask>0,2);
catch
    keyboard
end
