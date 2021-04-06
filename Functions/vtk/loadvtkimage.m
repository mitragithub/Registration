function [I,tline]=loadvtk(filename);

%opening file
fid=fopen(filename,'r','ieee-be');

%reading header information
for i=1:9;
	tline{i} = fgetl(fid);
end;


%reading data dimension
A=sscanf(tline{5},'%s%d%d%d');
Nx=A(11);   Ny=A(12);   Nz=A(13);

A=sscanf(tline{9},'%s');

%reading vector data
if (strcmp(A(1:7),'VECTORS')==1);
	[I, count] = fread(fid,Nx*Ny*Nz*3,'double');
	I=reshape(I,3,Nx,Ny,Nz);

%reading scalar data    
elseif (strcmp(A(1:7),'SCALARS')==1);
	[I, count] = fread(fid,Nx*Ny*Nz,'double');
	I=reshape(I,Nx,Ny,Nz);       

else;
	disp('ERROR: UNKNOWN DATA FORMAT')

end;

fclose(fid);
