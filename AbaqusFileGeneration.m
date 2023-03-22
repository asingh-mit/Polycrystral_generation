function [Filename] = AbaqusFileGeneration(co_or,CN,mattype,NSETS,PMI,GB)


Filename='Abaqus_file.inp';
fileID = fopen(Filename, 'w');

for i = 1:length(CN)
    Elements_Sets{i}.Name=strcat('EB',num2str(i));           
    Elements_Sets{i}.Elements_Type='S4R';
end


for i = 1:length(NSETS)
    NSET{i}.Name=strcat('NS',num2str(i));           
end

for i = 1:length(unique(mattype))
    Elemennt_CN{i} = CN(find(i == mattype),:);
end



fprintf(fileID,'********************************** P A R T S **********************************\n');
fprintf(fileID,'*PART, NAME=Part-Default\n');
fprintf(fileID,'**\n');



fprintf(fileID,'********************************** N O D E S **********************************\n');
fprintf(fileID,'*NODE, NSET=ALLNODES\n');

 for i=1:1:length(co_or)
     fprintf(fileID,[num2str(i) ',' num2str(co_or(i,1)) ',' num2str(co_or(i,2)) ',' num2str(0.0) '\n']); 
 end



 fprintf(fileID,'**\n');
 fprintf(fileID,'********************************** E L E M E N T S ****************************\n');
 

% fprintf(fileID,strcat('*ELEMENT, TYPE=',Elements_Sets{1}.Elements_Type,', ELSET=',Elements_Sets{1}.Name,'\n'));
 count = 1;
 for i = 1:length(unique(mattype))
     fprintf(fileID,strcat('*ELEMENT, TYPE=',Elements_Sets{i}.Elements_Type,', ELSET=',Elements_Sets{i}.Name,'\n'));
 

    for j = 1:length(Elemennt_CN{i}(:,1))
      fprintf(fileID,[num2str(count) ',' num2str(Elemennt_CN{i}(j,1)) ',' num2str(Elemennt_CN{i}(j,2)) ',' num2str(Elemennt_CN{i}(j,3)) ',' num2str(Elemennt_CN{i}(j,4)) '\n']); 
    count = count+1;
    end

 end
 fprintf(fileID,'**\n');


 fprintf(fileID,'********************************** N O D E S E T S **********************************\n');

  fprintf(fileID,'*NSET, NSET=left_1\n');
    for j = 1:length(NSETS{1})
      fprintf(fileID,[num2str(NSETS{1}(j)) ',' ]); 
    end

  fprintf(fileID,'\n*NSET, NSET=right_1\n');
    for j = 1:length(NSETS{2})
      fprintf(fileID,[num2str(NSETS{2}(j)) ',' ]); 
    end

  fprintf(fileID,'\n*NSET, NSET=bottom_1\n');
    for j = 1:length(NSETS{3})
      fprintf(fileID,[num2str(NSETS{3}(j)) ',' ]); 
    end

  fprintf(fileID,'\n*NSET, NSET=top_1\n');
    for j = 1:length(NSETS{4})
      fprintf(fileID,[num2str(NSETS{4}(j)) ',' ]); 
    end

   fprintf(fileID,'\n*NSET, NSET=GB\n');
    for j = 1:length(GB)
      fprintf(fileID,[num2str(GB(j)) ',' ]); 
    end

   fprintf(fileID,'\n*NSET, NSET=PMI\n');
    for j = 1:length(PMI)
      fprintf(fileID,[num2str(PMI(j)) ',' ]); 
    end

        fprintf(fileID,'\n*NSET, NSET=corner1\n');
     for j = 1:length(NSETS{5})
       fprintf(fileID,[num2str(NSETS{5}(j)) ',' ]); 
     end

 fprintf(fileID,'\n**\n');

fclose(fileID);