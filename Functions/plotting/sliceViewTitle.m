function sliceView_title(x,y,z,I,n)
n0 = 5;
if nargin == 1
    I = x;
    x = 1 : size(I,2);
    y = 1 : size(I,1);
    z = 1 : size(I,3);
    n = n0;
end
if nargin == 2
    I = x;
    n = y;
    x = 1 : size(I,2);
    y = 1 : size(I,1);
    z = 1 : size(I,3);
end
if nargin == 4
    n = n0;
end


if length(size(I)) == 3 || size(I,4) == 1
    nx = [size(I,2),size(I,1),size(I,3)];
    
    qlim = [0.01,0.99];
    clim = quantile(I(:),qlim);
    
    % last index fixed
    slices = linspace(1,nx(3),n+2);
    slices = slices(2:end-1);
    slices = round(slices);
    for i = 1 : length(slices)
        s = slices(i);
        subplotdan(3,n,i);
        imagesc(x,y,squeeze(I(:,:,s)),clim);
        axis image;
        axis off;
        title(['z = ' num2str(z(s))])
    end
    
    
    % second last index fixed
    slices = linspace(1,nx(1),n+2);
    slices = slices(2:end-1);
    slices = round(slices);
    for i = 1 : length(slices)
        s = slices(i);
        subplotdan(3,n,i+n);
        imagesc(z,y,squeeze(I(:,s,:)),clim);
        axis image;
        axis off;
        title(['x = ' num2str(x(s))])
    end
    
    % first index fixed
    slices = linspace(1,nx(2),n+2);
    slices = slices(2:end-1);
    slices = round(slices);
    for i = 1 : length(slices)
        s = slices(i);
        subplotdan(3,n,i+n+n);
        imagesc(z,x,squeeze(I(s,:,:)),clim);
        axis image;
        axis off;
        title(['y = ' num2str(y(s))])
    end
end




% color
if length(size(I))==4
    nx = [size(I,2),size(I,1),size(I,3)];
    
    qlim = [0.01,0.99];
    clim = quantile(I(:),qlim);
    
    if size(I,4) == 2
        I = cat(4,I,zeros(nx(2),nx(1),nx(3)));
    elseif size(I,4) > 3
        I = I(:,:,:,1:3);        
    end
    
    I = (I - clim(1))/diff(clim);
    
    
    % last index fixed
    slices = linspace(1,nx(3),n+2);
    slices = slices(2:end-1);
    slices = round(slices);
    for i = 1 : length(slices)
        s = slices(i);
        subplotdan(3,n,i);
        imagesc(x,y,squeeze(I(:,:,s,:)));
        axis image;
        axis off;
    end
    
    % second last index fixed
    slices = linspace(1,nx(1),n+2);
    slices = slices(2:end-1);
    slices = round(slices);
    for i = 1 : length(slices)
        s = slices(i);
        subplotdan(3,n,i+n);
        imagesc(z,y,squeeze(I(:,s,:,:)));
        axis image;
        axis off;
    end
    
    % first index fixed
    slices = linspace(1,nx(2),n+2);
    slices = slices(2:end-1);
    slices = round(slices);
    for i = 1 : length(slices)
        s = slices(i);
        subplotdan(3,n,i+n+n);
        imagesc(z,x,squeeze(I(s,:,:,:)));
        axis image;
        axis off;
    end
    
end


colormap gray;