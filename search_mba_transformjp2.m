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
function imgdir=search_mba_transformjp2(imgdir0,brainID)
% identify the image directories
imgdir=[imgdir0,'/',brainID,'/Transformation_OUTPUT/',brainID,'_img/reg_high_tif_pad_jp2/'];
if ~dircheck(imgdir)
    imgdir=[imgdir0,'/',brainID,'/Transformation_OUTPUT/',brainID,'_img/'];
    if ~dircheck(imgdir)
        imgdir=[imgdir0,'/',brainID,'/Transformation_OUTPUT/reg_high_tif_pad_jp2/'];
        if ~dircheck(imgdir)
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