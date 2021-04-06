function done(pmdno)
    filename=['/grid/mitra/home/xli/RegistrationPipelineV3/Data/' pmdno '/registration_done.txt'];
    fid=fopen(filename, 'wt');
    fprintf(fid, 'Astra');
    fclose(fid);
end
