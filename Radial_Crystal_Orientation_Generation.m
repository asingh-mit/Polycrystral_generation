function [Filename] = Radial_Crystal_Orientation_Generation(mattype,Theta_radial)


Filename='Radial_Crystal_Orientation.txt';
fileID = fopen(Filename, 'w');

for i = 1:1:mattype
    if i<mattype
     fprintf(fileID,'%4d %4d %4d\n',Theta_radial(i),0,0);
    elseif i==mattype
     fprintf(fileID,'%4d %4d %4d\n',0,0,0);
    end
end