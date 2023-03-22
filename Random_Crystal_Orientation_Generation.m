function [Filename] = Random_Crystal_Orientation_Generation(mattype,Crystal_orien)


Filename='Random_Crystal_Orientation.txt';
fileID = fopen(Filename, 'w');

% r = 0 + (0+180)*rand(mattype,1);

for i = 1:1:mattype
    if i<mattype
     fprintf(fileID,'%4d %4d %4d\n',Crystal_orien(i),0,0);
    elseif i==mattype
     fprintf(fileID,'%4d %4d %4d\n',0,0,0);
    end
end

