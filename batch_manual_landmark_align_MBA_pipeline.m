function batch_manual_landmark_align_MBA_pipeline(datadir,animalID,faillist)
% datadir='~/CSHLservers/mitragpu3/disk132/main/RegistrationData/';
% animalID='PTM849';
% faillist files are in mitramba1:/usr/local/MBA_QC/static/mnt/12
% faillist=['/Users/bhuo/CSHLservers/qcfiles/',animalID,'.txt'];
inputdir=[datadir,'/Data/',animalID,'/INPUT_DATA/'];
outputdir=[datadir,'/Data/',animalID,'/Registration_OUTPUT/'];
fid=fopen(faillist);
qcfilelist=textscan(fid,'%q');
qcfilelist=qcfilelist{1};
%% correct individual sections
% qcfilelist
for i=1:length(qcfilelist)
    [~,qcfilename,ext]=fileparts(qcfilelist{i});
    qcskip=0;
    if contains(qcfilename,'-F')
        disp(['Processing ',qcfilename,'...'])
        % search for Nissl failures
        k=strfind(qcfilename,'_preview');
        secnum=qcfilename(k-4:k-1);
        secind=find(contains(qcfilelist,secnum));
        if length(secind)>1
            if contains([qcfilelist{secind};],'-N')
                qcskip=1;
                disp(['Skipped ',qcfilename,' due to Nissl failure.'])
            end
        end
    else
        qcskip=1; % correct only fluorescent sections
    end
    if qcskip==0
        qcfile=[outputdir,qcfilename,ext];
        disp(qcfile)
        manual_landmark_align_MBA_pipeline(inputdir,qcfile);
        disp('Done.')
    end
end
%% apply to the stack
seg_file=[datadir,'/ATLAS/annotation_50.vtk'];
atlas_file = [datadir,'/ATLAS/ara_nissl_50.vtk'];
tic;
apply_deformation({seg_file,atlas_file},inputdir,outputdir,outputdir);
toc;