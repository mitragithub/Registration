% this script will make a figure of one slice (or a set of specified
% slices) from a given subject
% it will show
% 2D image with deformed grid
% atlas mapped onto the image with distortion
%
% how to do this
% first focus on the slice
% we need to load it
% then we need to apply the input to registered transform

clear all;
close all;
fclose all;
addpath Functions/vtk/
addpath Functions/plotting/



slices = [50 100 150 200 250 300];
output_dir = '3317_manuscript_figs/';
input_dir = 'PMD3317_test_00/';
target_dir = '/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/Xu2Daniel/PMD3317/';


seg_file = '/cis/home/dtward/Documents/ARA/Mouse_CCF/vtk/annotation_50.vtk';
atlas_file = '/cis/home/dtward/Documents/ARA/Mouse_CCF/vtk/ara_nissl_50.vtk';
contour_spacing = 1000; % contour distance in um (don't change)
xyrange = [-6000,6000]; % fix bounds of image (don't change)
maxscalerange = 2; % for detjac (don't change)

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end



[xI,yI,zI,I,title_,names] = read_vtk_image(atlas_file);
[xS,yS,zS,S,title_,names] = read_vtk_image(seg_file);



% get geometry file
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
    count = count + 1;
    % process this line, splitting at commas
    csv_data(count,:) = strsplit(line,',');
    %
end
fclose(fid);
files = csv_data(:,1);
nxJ = cellfun(@(x)str2num(x), csv_data(:,2:3));
dxJ = cellfun(@(x)str2num(x), csv_data(:,5:6));

x0J = cellfun(@(x)str2num(x), csv_data(:,8:9));
z0J = cellfun(@(x)str2num(x), csv_data(:,10));

zJ = z0J;
dzJ = cellfun(@(x) str2num(x), csv_data(:,7));

for f = 1 : length(files)
    xJ{f} = x0J(f,1) + (0:nxJ(f,1)-1)*dxJ(f,1);
    yJ{f} = x0J(f,2) + (0:nxJ(f,2)-1)*dxJ(f,2);
end



for f = slices
    
    [dir_,fname_,ext_] = fileparts(files{f});
    qlim = [0.001,0.999];
    for atloop = 1 : 2
        if atloop == 1
            J = imread([target_dir fname_ ext_]);
            J = double(J);
            [XJ,YJ] = meshgrid(xJ{f},yJ{f});
            % scale it
            clim = quantile(J(:),qlim);
            J = J - clim(1);
            J = J / (clim(2) - clim(1));
            
            
            % load the input to registered transform
            [x,y,z,tform,title_,names] = read_vtk_image([[input_dir 'registered_to_input_displacement_' fname_ '.vtk']]);
            [X,Y] = meshgrid(x,y);
            tform(:,:,:,1) = tform(:,:,:,1) + X;
            tform(:,:,:,2) = tform(:,:,:,2) + Y;
            
            % transform the image to straighen it
            Jreg = zeros(size(tform,1),size(tform,2),size(J,3));
            for c = 1 : size(J,3)
                F = griddedInterpolant({yJ{f},xJ{f}},J(:,:,c),'linear','none');
                Jreg(:,:,c) = F(tform(:,:,:,2),tform(:,:,:,1));
            end
            atstring = 'target';
        elseif atloop == 2
            % do the same thing but with atlas
            [x,y,z,tform,title_,names] = read_vtk_image([[input_dir 'registered_to_atlas_displacement_' fname_ '.vtk']]);
            [X,Y,Z] = meshgrid(x,y,z);
            tform(:,:,:,1) = tform(:,:,:,1) + X;
            tform(:,:,:,2) = tform(:,:,:,2) + Y;
            tform(:,:,:,3) = tform(:,:,:,3) + Z;
            F = griddedInterpolant({yI,xI,zI},I,'nearest','none');
            Idef = F(tform(:,:,:,2),tform(:,:,:,1),tform(:,:,:,3));
            clim = quantile(Idef(:),qlim);
            Idef = Idef - clim(1);
            Idef = Idef / (clim(2) - clim(1));
            
            Jreg = repmat(Idef,[1,1,3]);
            
            atstring = 'atlas';
        end
        
        
        
        
        danfigure(1);
        clf
        set(1,'defaultaxesfontsize',18)
        set(1,'defaultaxeslinewidth',2)
        set(1,'paperpositionmode','auto')
        imagesc(x/1000,y/1000,Jreg)
        axis image
        xlabel ('mm')
        ylabel ('mm')
        set(gca,'xlim',xyrange/1000,'ylim',xyrange/1000)
        h = colorbar();
        set(h,'visible','off')
        saveas(gcf,[output_dir fname_ '_' atstring '.png'])
        
        
        
        
        
        % okay now the next thing is to load the transform
        [x,y,z,tform,title_,names] = read_vtk_image([[input_dir 'registered_to_atlas_displacement_' fname_ '.vtk']]);
        [X,Y,Z] = meshgrid(x,y,z);
        tform(:,:,:,1) = tform(:,:,:,1) + X;
        tform(:,:,:,2) = tform(:,:,:,2) + Y;
        tform(:,:,:,3) = tform(:,:,:,3) + Z;
        
        
        
        % plot the grid
        danfigure(1);
        clf
        imagesc(x/1000,y/1000,Jreg)
        axis image
        xlabel ('mm')
        ylabel ('mm')
        set(gca,'xlim',xyrange/1000,'ylim',xyrange/1000)
        
        hold on;
        levels = 0 : contour_spacing : 10000;
        levels = [-flip(levels) levels(2:end)];
        contour(x/1000,y/1000,tform(:,:,:,1),levels,'k','linewidth',3)
        contour(x/1000,y/1000,tform(:,:,:,2),levels,'k','linewidth',3)
        contour(x/1000,y/1000,tform(:,:,:,3),levels,'k','linewidth',3)
        
        contour(x/1000,y/1000,tform(:,:,:,1),levels,'w','linewidth',1)
        contour(x/1000,y/1000,tform(:,:,:,2),levels,'w','linewidth',1)
        contour(x/1000,y/1000,tform(:,:,:,3),levels,'w','linewidth',1)
        
        h = colorbar();
        set(h,'visible','off')
        saveas(gcf,[output_dir fname_ '_' atstring '_grid' '.png'])
        
        
        
        
        
        % plot the detjac
        [x,y,z,detjac,title_,names] = read_vtk_image([[input_dir 'atlas_to_registered_detjac_' fname_ '.vtk']]);
        % here's what I want to do
        % first take log
        scale = log(detjac.^(1/3));
        % map it to 0,1
        scalerange = [-log(maxscalerange),log(maxscalerange)];
        scale = (scale - scalerange(1))/diff(scalerange);
        danfigure(1);
        clf
        Jreg(Jreg<0) = 0;
        Jreg(Jreg>1) = 1;
        Jreg(isnan(Jreg)) = 0;
        imagesc(x/1000,y/1000,Jreg)
        
        alpha = (((scale-median(scale(:))).^2).^0.25)*2;
        alpha(alpha > 1) = 1;
        alpha = abs(scale-0.5)*2;
        
        colors = jet(256);
        % map scale to colors
        F = griddedInterpolant(linspace(0,1,256),colors(:,1),'linear','nearest');
        scaleR = reshape(F(scale(:)),size(scale));
        F = griddedInterpolant(linspace(0,1,256),colors(:,2),'linear','nearest');
        scaleG = reshape(F(scale(:)),size(scale));
        F = griddedInterpolant(linspace(0,1,256),colors(:,3),'linear','nearest');
        scaleB = reshape(F(scale(:)),size(scale));
        scaleRGB = cat(3,scaleR,scaleG,scaleB);
        
        %     alpha = 0.5;
        Ishow = bsxfun(@times,Jreg,(1-alpha)) + bsxfun(@times, scaleRGB, (alpha));
        imagesc(x/1000,y/1000,Ishow)
        colormap jet
        set(gca,'clim',scalerange)
        h = colorbar;
        nticks = 5;
        ticks = linspace(scalerange(1),scalerange(2),nticks);
        set(h,'ticks',ticks)
        ticklabels = {};
        for i = 1 : nticks
            ticklabels{i} = num2str(exp(ticks(i)),'%2.1f');
        end
        set(h,'ticklabels',ticklabels)
        axis image
        xlabel ('mm')
        ylabel ('mm')
        
        set(gca,'xlim',xyrange/1000,'ylim',xyrange/1000)
        saveas(gcf,[output_dir fname_ '_' atstring '_detjac' '.png'])
        
        
        
        
        % plot the segs
        [x,y,z,tform,title_,names] = read_vtk_image([[input_dir 'registered_to_atlas_displacement_' fname_ '.vtk']]);
        [X,Y,Z] = meshgrid(x,y,z);
        tform(:,:,:,1) = tform(:,:,:,1) + X;
        tform(:,:,:,2) = tform(:,:,:,2) + Y;
        tform(:,:,:,3) = tform(:,:,:,3) + Z;
        F = griddedInterpolant({yS,xS,zS},S,'nearest','none');
        Sdef = F(tform(:,:,:,2),tform(:,:,:,1),tform(:,:,:,3));
        Sdef(isnan(Sdef)) = 0;
        rng(1);
        colors = rand(256,3);
        colors(1,:) = 0;
        Sdef = mod(Sdef,256);
        SdefR = reshape(colors(Sdef(:)+1,1),size(Sdef));
        SdefG = reshape(colors(Sdef(:)+1,2),size(Sdef));
        SdefB = reshape(colors(Sdef(:)+1,3),size(Sdef));
        SdefRGB = cat(3,SdefR,SdefG,SdefB);
        
        danfigure(1);
        clf
        alpha = 0.5;
        imagesc(x/1000,y/1000,Jreg*(1-alpha) + SdefRGB*alpha)
        axis image
        xlabel ('mm')
        ylabel ('mm')
        set(gca,'xlim',xyrange/1000,'ylim',xyrange/1000)
        h = colorbar;
        set(h,'visible','off')
        saveas(gcf,[output_dir fname_ '_' atstring '_segs' '.png'])
        
        
        
    end
end