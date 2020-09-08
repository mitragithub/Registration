% search_mba_transformjp2.m
% This script searches for the correct directory that contains the
% final transformation output images
% Input:
%   - imgdir0: a string of parent directory on M drive that contains all the
%   registration results (e.g. 'M32/RegistrationData/Data/')
%   - brainID: a string of brain ID (e.g. 'PMD1234')
% Output:
%   - imgdir: a string that contains the correct directory with transformed
%   images.
function imgdir=search_mba_transformjp2(imgdir0,brainID,transformdir)
if nargin<3
    % identify transformation output dir
    transformdir=[imgdir0,'/',brainID,'/Transformation_OUTPUT/'];
end
% identify the image directories
imgdir=[transformdir,brainID,'_img/reg_high_tif_pad_jp2/'];
jp2files=dir([imgdir,'*.jp2']);
if ~exist(imgdir,'dir') || isempty(jp2files)
    imgdir=[transformdir,brainID,'_img/'];
    jp2files=dir([imgdir,'*.jp2']);
    if ~exist(imgdir,'dir') || isempty(jp2files)
        imgdir=[transformdir,'/reg_high_tif_pad_jp2/'];
        jp2files=dir([imgdir,'*.jp2']);
        if ~exist(imgdir,'dir') || isempty(jp2files)
            imgdir=[];
            warning('Cannot find registered images!')
        end
    end
end
end

function tf=dircheck(imgdir)
if exist(imgdir,'dir')
    A=dir([imgdir,'/*.jp2']);
    if ~isempty(A)
        tf=true;
    else
        tf=false;
    end
else
    tf=false;
end
end