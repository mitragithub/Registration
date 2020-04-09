function combine_nissl_and_fluoro_transforms(detailed_output_dir)

% we load the following
% from the detailed output dir
% we load down1_v.mat
% and down1_A.mat
% and also
% NtoF.mat

% if doing example, just load low res
if ~isempty(strfind(detailed_output_dir,'example'))
varsA = load([detailed_output_dir 'down4_A.mat']);
varsV = load([detailed_output_dir 'down4_v.mat']);
else
varsA = load([detailed_output_dir 'down1_A.mat']);
varsV = load([detailed_output_dir 'down1_v.mat']);
if exist([detailed_output_dir 'down1edit_A.mat'],'file')
    disp(['Using edited files to combine nissl fluoro'])
    varsA = load([detailed_output_dir 'down1edit_A.mat']);
    varsV = load([detailed_output_dir 'down1edit_v.mat']);
end
end
varsNtoF = load([detailed_output_dir 'NtoF.mat']);

%%
% check that we have the right sizes
sum(varsNtoF.is_nissl)
size(varsA.AJ,3)
vJtx = {};
vJty = {};
xJ = {};
yJ = {};

% loop through the slices
count = 0;
for i = 1 : length(varsNtoF.is_fluoro)
    if varsNtoF.is_nissl(i)
        count = count + 1;
        % for a nissl stain, use existing values
        AJ(:,:,i) = varsA.AJ(:,:,count);
        vJtx{i} = varsV.vJtx{count};
        vJty{i} = varsV.vJty{count};
        xJ{i} = varsV.xJ{count};
        yJ{i} = varsV.yJ{count};

        
    else
        % for a fluoro stain, assign value of nearest slice
        ind = varsNtoF.inds(i);
        all_inds = cumsum(varsNtoF.is_nissl);
        AJ(:,:,i) = varsNtoF.NtoF(:,:,i) * varsA.AJ(:,:,all_inds(ind));
        
        vJtx{i} = varsV.vJtx{all_inds(ind)};
        vJty{i} = varsV.vJty{all_inds(ind)};

        % need to save x and y, they are used for interpretting v
        % but their values do not matter because this is rigid only
        xJ{i} = [0,1];
        yJ{i} = [0,1];

        
    end
end

% other vars
vtx = varsV.vtx;
vty = varsV.vty;
vtz = varsV.vtz;
xI = varsV.xI;
yI = varsV.yI;
zI = varsV.zI;
A = varsA.A;

%%
% save data
save([detailed_output_dir 'combined_A.mat'],'A','AJ')
save([detailed_output_dir 'combined_v.mat'],'xI','yI','zI','vtx','vty','vtz','vJtx','vJty','xJ','yJ')