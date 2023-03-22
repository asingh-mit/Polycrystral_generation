function [Filename,Scale] = Moose_Mesh_Section(mattype)

Filename='MESH.i';
fileID = fopen(Filename, 'w');

Scale = 11e-6;
Scaled_domain = [Scale Scale 0];

fprintf(fileID,'[Mesh]\n');
fprintf(fileID,'[initial_mesh]\n');
fprintf(fileID,'type = FileMeshGenerator\n');
fprintf(fileID,'file = GMSH_file.msh\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./scale]\n');
fprintf(fileID,'type = TransformGenerator\n');
fprintf(fileID,'input = initial_mesh\n');
fprintf(fileID,'transform = SCALE\n');
fprintf(fileID,[strcat('vector_value=','''',num2str(Scaled_domain)),'''' '\n']);
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./nodesets]\n');
fprintf(fileID,'type = SideSetsFromNormalsGenerator\n');
fprintf(fileID,'input = scale\n');
fprintf(fileID,'normals = ''-1 0 0\n');
fprintf(fileID,'1 0 0\n');
fprintf(fileID,'0 -1 0\n');
fprintf(fileID,'0 1 0''\n');
fprintf(fileID,'fixed_normal = true\n');
fprintf(fileID,'new_boundary = ''left right bottom top'' \n'); 
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./rename_blocks]\n');
fprintf(fileID,'type = RenameBlockGenerator\n');
fprintf(fileID,'input = nodesets\n');
fprintf(fileID,'old_block = ''');
fprintf(fileID,[num2str(1:1:mattype)]); 
fprintf(fileID,'''\n');
fprintf(fileID,'new_block = ''');
fprintf(fileID,[num2str(0:1:mattype-1)]); 
fprintf(fileID,'''\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./pmi_boundary]\n');
fprintf(fileID,'type = SideSetsAroundSubdomainGenerator\n');
fprintf(fileID,'input = rename_blocks\n');
fprintf(fileID,'block=''');
fprintf(fileID,[num2str(0:1:mattype-2)]);
fprintf(fileID,'''\n');
fprintf(fileID,'new_boundary = ''PMI''\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');
 

for i = 1:mattype-1
     if i == 1
     fprintf(fileID,[strcat('[./','gb_boundary',num2str(i),']') '\n']); 
     fprintf(fileID,'type = SideSetsBetweenSubdomainsGenerator\n');
     fprintf(fileID,'input = pmi_boundary\n');
     fprintf(fileID,[strcat('primary_block=',num2str(i-1)) '\n']);
     ind = 0:1:mattype-2;
     ind(ind == i-1) = []; 
     fprintf(fileID,[strcat('paired_block=','''',num2str(ind),'''') '\n']);
     fprintf(fileID,'new_boundary = ''GB1''\n');
     fprintf(fileID,'[]\n');
     fprintf(fileID,'\n');
     else
     fprintf(fileID,[strcat('[./','gb_boundary',num2str(i),']') '\n']); 
     fprintf(fileID,'type = SideSetsBetweenSubdomainsGenerator\n');
     fprintf(fileID,[strcat('input=gb_boundary',num2str(i-1)) '\n']);
     fprintf(fileID,[strcat('primary_block=',num2str(i-1)) '\n']);
     ind = 0:1:mattype-2;
     ind(ind == i-1) = []; 
     fprintf(fileID,[strcat('paired_block=','''',num2str(ind),'''') '\n']);
     fprintf(fileID,[strcat('new_boundary=','''','GB',num2str(i)),'''' '\n']); 
     fprintf(fileID,'[]\n');
     fprintf(fileID,'\n');
     end
 end

fprintf(fileID,'[./left_bottom]\n');
fprintf(fileID,'type = ExtraNodesetGenerator\n');
fprintf(fileID,'new_boundary = ''left_bottom''\n');
fprintf(fileID,[strcat('coord=','''',num2str([0 0 0])),'''' '\n']);
fprintf(fileID,[strcat('input=''gb_boundary',num2str(mattype-1)),'''' '\n']); 
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./left_top]\n');
fprintf(fileID,'type = ExtraNodesetGenerator\n');
fprintf(fileID,'new_boundary = ''left_top''\n');
fprintf(fileID,[strcat('coord=','''',num2str([0 Scale 0])),'''' '\n']);
fprintf(fileID,'input = left_bottom\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./right_bottom]\n');
fprintf(fileID,'type = ExtraNodesetGenerator\n');
fprintf(fileID,'new_boundary = ''right_bottom''\n');
fprintf(fileID,[strcat('coord=','''',num2str([Scale 0 0])),'''' '\n']);
fprintf(fileID,'input = left_top\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./right_top]\n');
fprintf(fileID,'type = ExtraNodesetGenerator\n');
fprintf(fileID,'new_boundary = ''right_top''\n');
fprintf(fileID,[strcat('coord=','''',num2str([Scale Scale 0])),'''' '\n']);
fprintf(fileID,'input = right_bottom\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'construct_side_list_from_node_list=true\n');
fprintf(fileID,'\n');

fprintf(fileID,'[]\n');