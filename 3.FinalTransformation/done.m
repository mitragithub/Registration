function done(pmdno)
    filename=['/sonas-hs/mitra/hpc/home/xli/RegistrationPipelineV3/Data/' pmdno '/transform_done.txt'];
    fid=fopen(filename, 'wt');
    fprintf(fid, 'Astra');
    fclose(fid);
end
