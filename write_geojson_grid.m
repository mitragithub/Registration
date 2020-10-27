function write_geojson_grid(fname,names,curves)
% ids is an array of ints
% curves is nested cells
% curves{i} is the multi polygon for name {i}
% curves{i}{j} is the jth polygon in multi polygon i
% note for polygons the first and last point is the same
%
% note about names
% for allen atlas contours Cold Spring Harbor is doing this
% {"type":"Feature","id":"84","properties":{"name":"Prelimbic area, layer 6a","acronym":"PL6a"}
% note about style
% all braces and brackets start on the next line
% with same indentation as their key
%
% note download allen ontology instructions
% Download the "Mouse Brain Atlas" ontology as a hierarchically structured json file
% http://api.brain-map.org/api/v2/structure_graph_download/1.json



fid = fopen(fname,'wt');
fprintf(fid,'{\n'); % start of the file, no leading space
fprintf(fid,'  "type": "FeatureCollection",\n'); % start of features, 2 leading spaces
fprintf(fid,'  "features":\n'); % start of features, 2 leading spaces
fprintf(fid,'  [\n'); % open bracket, same indent as above



n = length(names);
for i = 1 : n

fprintf(fid,'    {\n'); % start this feature, 4 spaces
fprintf(fid,'      "type": "Feature",\n');
%fprintf(fid,['      "id": "' ids(i) '",\n']);


% now properties
fprintf(fid,'      "properties":\n');
fprintf(fid,'      {\n');
fprintf(fid,['        "name": "' names{i} '",\n']);
%fprintf(fid,['        "acronym": "' acronym '"\n']);
fprintf(fid,'      },\n'); % end of properties




fprintf(fid,'      "geometry":\n');
fprintf(fid,'      {\n');
fprintf(fid,'        "type": "MultiLineString",\n');
fprintf(fid,'        "coordinates":\n');
fprintf(fid,'        [\n');
% now loop through the coordinates
ni = length(curves{i});
for j = 1 : ni
fprintf(fid,'          [\n');
% fprintf(fid,'            [\n'); 
% loop over the points
fprintf(fid,'               ');     
for k = 1 : size(curves{i}{j},2)
fprintf(fid,['[' num2str(curves{i}{j}(1,k)) ', '  num2str(curves{i}{j}(2,k))  ']']);
if k < size(curves{i}{j},2)
    fprintf(fid,[', ']);
end
end
fprintf(fid,'\n');
% fprintf(fid,'            ]\n'); % end of inner brackets for this polygon
fprintf(fid,'          ]'); % end of this polygon
if j < ni
    fprintf(fid,',\n');
else
    fprintf(fid,'\n');
end
end


fprintf(fid,'        ]\n'); % end of coordinates
fprintf(fid,'      }\n'); % end of geometry






fprintf(fid,'    }'); % end of this feature, 4 leading spaces
if i < n
    fprintf(fid,',\n');
else
    fprintf(fid,'\n');
end

end



fprintf(fid,'  ]\n'); % end of features, 2 leading spaces
fprintf(fid,'}\n'); % end of file, no leading spaces


fclose(fid);