function savevtk(I,filename);

s=size(I);
%this is vector data
if (length(s)==4);
    tline{1}=sprintf('# vtk DataFile Version 3.0\n');
    tline{2}=sprintf('lddmm %d %d %d %d %d %d %d %d %d %d %d %d %d\n',8,0,0,s(4)-1,s(4)-1,0,0,s(3)-1,s(3)-1,0,0,s(2)-1,s(2)-1);
    tline{3}=sprintf('BINARY\n');
    tline{4}=sprintf('DATASET STRUCTURED_POINTS\n');
    tline{5}=sprintf('DIMENSIONS %d %d %d\n',s(2),s(3),s(4));
    tline{6}=sprintf('SPACING %d %d %d\n',1,1,1);
    tline{7}=sprintf('ORIGIN %d %d %d\n',0,0,0);
    tline{8}=sprintf('POINT_DATA %d\n',s(2)*s(3)*s(4));
    tline{9}=sprintf('VECTORS Vectors_ double\n');

%this is scalar data
elseif (length(s)==3);
    tline{1}=sprintf('# vtk DataFile Version 3.0\n');
    tline{2}=sprintf('lddmm %d %d %d %d %d %d %d %d %d %d %d %d %d\n',8,0,0,s(3)-1,s(3)-1,0,0,s(2)-1,s(2)-1,0,0,s(1)-1,s(1)-1);
    tline{3}=sprintf('BINARY\n');
    tline{4}=sprintf('DATASET STRUCTURED_POINTS\n');
    tline{5}=sprintf('DIMENSIONS %d %d %d\n',s(1),s(2),s(3));
    tline{6}=sprintf('SPACING %d %d %d\n',1,1,1);
    tline{7}=sprintf('ORIGIN %d %d %d\n',0,0,0);
    tline{8}=sprintf('POINT_DATA %d\n',s(1)*s(2)*s(3));
    tline{9}=sprintf('SCALARS Scalars_ double\n');
end;

    fid=fopen(filename,'w','ieee-be');
    for i=1:9;
        fprintf(fid,tline{i});
    end;

    I=I(:);
    count = fwrite(fid,I,'double');
    fclose(fid);

