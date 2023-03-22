function [co_or,CN,mattype,NSETS] = MeshGeneration(nnx,nny,img)

for count=1
    for j=1:nny-1
        for i=1:nnx-1
            kN=[i+nnx*(j-1) i+nnx*(j-1)+1 nnx+i+1+nnx*(j-1) nnx+i+nnx*(j-1)];
            CN(count,:)=kN; % storing connectivity values for each element
            count=count+1;
        end
    end
end

for count=1
for j=0:nny-1
    for i=0:nnx-1
        coor=[i/(nnx-1) j/(nny-1)];
        co_or(count,:)=coor;
        count=count+1;
    end
end
end


%==========================================================================
%==================Generating node sets ===================================
%==========================================================================
for i=1:nny
       NS_left(i) = 1 + (i-1)*nny;
end

for i=1:nny
       NS_right(i) = i*nny;
end

for i=1:nnx
       NS_bottom(i) = i;
end

for i=1:nnx
       NS_top(i) = nnx*nny - nnx + i;
end

NSETS{1} = NS_left;
NSETS{2} = NS_right;
NSETS{3} = NS_bottom;
NSETS{4} = NS_top;

NSETS{5} = [1; nnx; 1+(nny-1)*nny; nny*nny];
%==========================================================================
%==========================================================================


numnp=nnx*nny;
numelv=((nnx-1)*(nny-1));
numatv=1;

filename='Output';

vtk = fopen(strcat(filename,'.vtk'),'w');

fprintf('Output.vtk \n');

fprintf(vtk,'# vtk DataFile Version 2.0 \n');
fprintf(vtk,'Fracture Distribution \n');
fprintf(vtk,'ASCII \n \n');
fprintf(vtk,'DATASET UNSTRUCTURED_GRID\n');

fprintf(vtk,'POINTS  %d  double \n',numnp);
for i=1:numnp
    fprintf(vtk,'%22.15E %22.15E %22.15E \n',co_or(i,:),0);
end

fprintf(vtk,'\nCELLS   %5d %5d\n',numelv,5*numelv);
for i=1:length(CN(:,1))
    % paraview is 0-indexed
            fprintf(vtk,'%9d %8d %8d %8d %8d\n',4,CN(i,1:4)-1);
end

fprintf(vtk,'\nCELL_TYPES %7d\n',length(CN));
for i=1:length(CN)
    % use hexahedron elements
    fprintf(vtk,'%4d\n',9);
end
%==============================================================

fprintf(vtk,'\nCELL_DATA %4d\n',numelv);
fprintf(vtk,'SCALARS cell_scalars int 1\n');
fprintf(vtk,'LOOKUP_TABLE default\n');

img = imrotate(img,-90);

count =1;
for i = 1:length(img)-1
    for j = 1:length(img)-1
        mattype(count) = img(j,i);
        count = count + 1;
    end
end

for i=1:numelv
    fprintf(vtk,'%4d\n',mattype(i));
end

end