function batch_manual_landmark_align_MBA_pipeline(datadir,animalID,faillist)
% datadir='~/CSHLservers/mitragpu3/disk132/main/RegistrationData/';
% animalID='PTM849';
% faillist='~/CSHLservers/mitragpu3/disk132/main/temp/PMD963_failed_qc.txt';
inputdir=[datadir,'/Data/',animalID,'/INPUT_DATA/'];
outputdir=[datadir,'/Data/',animalID,'/Registration_OUTPUT/'];
fid=fopen(faillist);
qcfilelist=textscan(fid,'%q');
qcfilelist=qcfilelist{1};
%% correct individual sections
% qcfilelist
for i=1:length(qcfilelist)
    [~,qcfilename,ext]=fileparts(qcfilelist{i});
    if contains(qcfilename,'-F')
        qcskip=0;
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
    end
    if qcskip==0
        qcfile=[outputdir,qcfilename,ext];
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