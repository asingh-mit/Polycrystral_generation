function [Filename] = GMSHFileGeneration(co_or,CN,mattype)


Filename='GMSH_file.msh';
fileID = fopen(Filename, 'w');


fprintf(fileID,'$MeshFormat\n');
fprintf(fileID,'2.2 0 8\n');
fprintf(fileID,'$EndMeshFormat\n');


fprintf(fileID,'$MeshVersion\n');
fprintf(fileID,'$2.2.2\n');
fprintf(fileID,'$EndMeshVersion\n');


fprintf(fileID,'$Nodes\n');
fprintf(fileID,'%d\n',length(co_or));
 for i=1:1:length(co_or)
     fprintf(fileID,'%d %22.15E %22.15E %22.15E \n',i,co_or(i,:),0); 
 end
fprintf(fileID,'$EndNodes\n');


 fprintf(fileID,'$Elements\n');
 fprintf(fileID,'%d\n',length(CN));

 for i = 1:length(CN)
    fprintf(fileID,'%d %d %d %d %d %d %d %d %d %d \n',i,3,3,mattype(i),mattype(i),0,CN(i,1:4));
 end

fprintf(fileID,'$EndElements\n');

fclose(fileID);