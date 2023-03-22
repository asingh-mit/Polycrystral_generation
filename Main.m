% MATLAB code to generate the 2D microstructure and structured mesh with 
% 1. Different aspect ratio of grains with radial morphology
% 2. Different crystal orientation of grains with mean and std values
% 3. Generate .gmsh and .inp mesh file to directly load into the MOOSE
% 4. Outputs .vtk file for visualization
% 5. Generate MOOSE mesh input files

% Developed by Avtar singh (MIT and NREL) and Hongyi Xu (Uconn)

clear; clc; close all;

img_generation_switch = 1; % 1 - if want to generate new microstructure
                           % 0 - if want to load previously generated image


 Radial_crystal_orien_mean =  0;
 Radial_crystal_orien_std =  5;

if img_generation_switch == 1

    num = 83;           % The number of grains
    L = 201;            % Length of the image space. It must be an odd number!
    R = (L-(L-1)/6)/2;  % The radius of the particle. The coordinate of center point is [R+1, R+1]

    void_d_min = 0; % I see a void in the middle of the grain from the EBSD. This is the min radius of the void. To be obtained by image characterization
    glr = 5;        % Mean of grain "radius".   It should be obtained from image characterization
                    % Increasing glr increases the length of the grains
                    % Minimum value can be 0.01
    glrMIN = 1;     % minimum grain "radius". It should be obtained from image characterization
                    % Increasing glrMIN increases the length of the grains
                    % Minimum value = 0.0


% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% R
% In this script, the so-called "radius" and "diameter" are all along
% the radial direction (from image center to the edge).
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% -------------------------------------------------------------------------
% cdist = normrnd( R/2+10, (R/2)/5, [num,1] ); % For the grain center
% distance to img center, sampling in the linear space will not work. It is
% because we will get much denser grain centers near the img center. The
% resultant granular crystalline image looks weird. To resolve this issue,
% we propose to sample the distance d in the d^nd space, where nd > 1 
% nd can be used as an microtructure descriptor of the spatial
% dispersion of grain centers. This idea is inspired by:
% http://kanga.usask.ca/Geant/Exercises/sampling2.html 
                                                  
ang = lhsdesign(num,1) * 360; % In degree
nd = 1.5;      % Controls the length of the grain, for no effect put nd = 1
cdistmp = rand(num*2, 1);   % Generate more grain centers than we need, then pick num centers, which meet other requirements, randomly.
cdistmp = cdistmp.^(1/nd);
cdistmp = cdistmp*R;

cdist = cdistmp( find( cdistmp>void_d_min+glrMIN   &   cdistmp<R-glrMIN   ) );

if length(cdist) < num
    error('The number of sample points is fewer than the required number. Increase the factor value in cdistmp = rand(num*2, 1)');
else
    cdist = cdist(1:num);
end

figure(); 
plot(cdist.*cos(ang), cdist.*sin(ang),'.'); hold on;
axis equal; 
title('Grain center distribution')
% -------------------------------------------------------------------------



% Assume the range of the grain length is +-glr pixels from the grain center, uniform distribution. 
% The number & distribution func should be obtained by image characterization.
cdist_inner = cdist - rand(num,1)*glr;  % Distance to the image center, inner end
cdist_outer = cdist + rand(num,1)*glr;  % Distance to the image center, outer end
figure(); 
plot(cdist_inner); hold on; 
plot(cdist_outer); 
title('Inner-outer distances to image center');

for ii = 1:1:num
    % Check the inner end of each line segment is within a reasonable range
    if cdist_inner(ii) < void_d_min  
        cdist_inner(ii) = void_d_min;
    end

    if cdist_inner(ii) > R-glr
        cdist_inner(ii) = R-glr;
    end

    % Check the outer end of each line segment is within a reasonable range
    if cdist_outer(ii) < cdist_inner(ii) + glrMIN*2
        cdist_outer(ii) = cdist_inner(ii) + glrMIN*2;
    end

    if cdist_outer(ii) > R
        cdist_outer(ii) = R;
    end
end

figure(); plot(cdist_outer-cdist_inner, 'o-');
 
% Convert to the XY coordinate and round the number
x_inner = [];   y_inner = [];
x_outer = [];   y_outer = [];
x_cen = [];     y_cen = [];  % Record the center coordinate for debug use
for ii = 1:1:num
    theta = ang(ii,1) /360 * 2*pi;
    x_inner = [ x_inner; cdist_inner(ii)*cos(theta) + (L-1)/2];
    y_inner = [ y_inner; cdist_inner(ii)*sin(theta) + (L-1)/2];
    
    x_outer = [ x_outer; cdist_outer(ii)*cos(theta) + (L-1)/2];
    y_outer = [ y_outer; cdist_outer(ii)*sin(theta) + (L-1)/2];
    
    x_cen = [ x_cen; cdist(ii)*cos(theta) + (L-1)/2];
    y_cen = [ y_cen; cdist(ii)*sin(theta) + (L-1)/2];
end


%======================================================================
%===============Control the random or radial orientation===============
%======================================================================
for i=1:length(x_inner)
X = [x_inner(i) x_outer(i)];
Y = [y_inner(i) y_outer(i)];
angle_limit = 0;
ang = randi([0 angle_limit],1,1);  % degrees
Xc = x_cen(i);
Yc = y_cen(i);
Xrot =  (X-Xc)*cosd(ang) + (Y-Yc)*sind(ang) + Xc;
Yrot = -(X-Xc)*sind(ang) + (Y-Yc)*cosd(ang) + Yc;

x_inner(i) = Xrot(1,1);
y_inner(i) = Yrot(1,1);

x_outer(i) = Xrot(1,2);
y_outer(i) = Yrot(1,2);
end


figure(); hold on;
for ii = 1:1:num
    plot( [x_inner(ii), x_outer(ii)], [y_inner(ii), y_outer(ii)] ); hold on;
end
plot(x_cen, y_cen, 'o'); axis equal; title('Grain center distribution')
title('Center lines of all grains');

%==========================================================================
%==========================================================================




% What I did here is just to traverse all pixel in the image and calculate
% each pixel's distances to all line sections, and then pick the shortest
% one. As two traversals are used here, it is not efficient. There should
% be a way to improve it.

 img = ones( L, L )*(num+1); % Center point of the grain center line
for xx = 1:1:L
    disp(L-xx+1);
    for yy = 1:1:L

        pnt = [xx,yy];
        
        dist_rec = [];
        for ii = 1:1:num
            v1 = [ x_inner(ii), y_inner(ii) ];
            v2 = [ x_outer(ii), y_outer(ii) ];
            dist = point_to_line(pnt, v1, v2);
            dist_rec = [ dist_rec; dist];
        end

        idx = find(dist_rec == min(dist_rec)); % Location of the min distance value, also the number of the grain/center line
        if length(idx)>1
            idx = idx(1);
        end
        
        if (xx-(L-1)/2)^2 + (yy-(L-1)/2)^2 >= void_d_min^2   &&   (xx-(L-1)/2)^2 + (yy-(L-1)/2)^2 <= R^2
            img(xx,yy) = idx;  % -2 means void. Positive number i means the i-th particle
        end

    end
end

 figure(); 
 imagesc(img); 
 title('Grid points of the grains');
 axis equal;

%==========================================================================
%==========================================================================

    save('img.mat','img');

end


if img_generation_switch == 0
    example = matfile('img.mat');
    img = example.img;
    L = length(img);
    R = (L-(L-1)/6)/2;  % The radius of the particle. 
                         % The coordinate of center point is [R+1, R+1]
end

 %=========================================================================
 %================Plot boundary of the grains==============================
 %=========================================================================
 figure()
for i=1:length(unique(img))-1

     [temp_x,temp_y] = find(i==img(:,:));

angc = 90;  % degrees
Xcc = (L-1)/2;  % X cener of rotation
Ycc = (L-1)/2;  % Y center of rotation
% Shift X/Y to the rotation center
Xshiftc = temp_x - Xcc;        
Yshiftc = temp_y - Ycc;
% Rotate the coordinates
Xsrotc =  Xshiftc*cosd(angc) + Yshiftc*sind(angc);
Ysrotc = -Xshiftc*sind(angc) + Yshiftc*cosd(angc);
% Shift the rotated coordinates back to the original reference center
coor_x{i} = Xsrotc + Xcc;
coor_y{i} = Ysrotc + Ycc;


 k{i} = boundary(coor_x{i},coor_y{i});

 plot(coor_x{i}(k{i}),coor_y{i}(k{i})); hold on
  end
  title('Boundary of the grains');
  axis equal
%==========================================================================
%==========================================================================





%==========================================================================
%============Fitting grains into an ellipse================================
%=======for calculating the aspect ratio of the grains=====================
%==========================================================================
figure()
Color = jet(length(unique(img))-1);

for i = 1:length(unique(img))-1
N = length(coor_x{i});
X = [coor_x{i}, coor_y{i}];

G = i*ones(N,1);
 
gscatter(X(:,1), X(:,2), G, Color(i,:),'o',1.0); hold on;

     Mu = mean(X);
     X0 = bsxfun(@minus, X, Mu);
 
STD = 1.25;                     %# 2 standard deviations
conf = 2*normcdf(STD)-1;        %# covers around 95% of population
scale = chi2inv(conf,2);        %# inverse chi-squared with dof=#dimensions

Cov = cov(X0) * scale;
[V D] = eig(Cov);

    t = linspace(0,2*pi,100);
    e = [cos(t) ; sin(t)];        %# unit circle
    VV = V*sqrt(D);               %# scale eigenvectors
    e = bsxfun(@plus, VV*e, Mu'); %#' project circle back to orig space

%---- plot cov and major/minor axes-----
   plot(e(1,:), e(2,:), 'Color','k','linewidth',1.5); hold on;
    quiver(Mu(1),Mu(2), VV(1,1),VV(2,1), 'Color','k','linewidth',1); hold on;
    quiver(Mu(1),Mu(2), VV(1,2),VV(2,2), 'Color','k','linewidth',1); hold on;
axis equal

AR(i) = max(max(sqrt(D)))/min(max(sqrt(D)));

center(i,1) = Mu(:,1);
center(i,2) = Mu(:,2);

plot(Mu(:,1),Mu(:,2),'o'); hold on;

if sqrt(D(1,1))>sqrt(D(2,2))
    Phi_pt_x(i) = Mu(1) + VV(1,1);
    Phi_pt_y(i) = Mu(2) + VV(2,1);
else
    Phi_pt_x(i) =  Mu(1) + VV(1,2);
    Phi_pt_y(i) =  Mu(2) + VV(2,2);
end

end
title('Fitting grain data points into an ellipse');

figure()
histogram(AR,floor(length(AR)/4),'facecolor','r','Edgecolor','k');
xlabel('Aspect ratio');
ylabel('Number of grains');
set(gca,'linewidth',2,'fontsize',15,'fontname','arial','fontweight','normal','tickdir','out','fontangle','italic');


Mean_AR = sum(AR)/length(AR);
std_AR = std(AR);
%==========================================================================
%==========================================================================






%==========================================================================
%============Calculating the crystal morphological direction===============
%==========================================================================
acolor = jet(length(center));
figure()
x2 = L;
y2 = (L-1)/2;

x1 = (L-1)/2;
y1 = (L-1)/2;

for i = 1:length(center)

x3 = center(i,1);
y3 = center(i,2);

vector1 = [x2,y2,0] - [x1,y1,0];
vector2 = [x3,y3,0] - [x1,y1,0];

vc = cross(vector1, vector2);

if vc(3)<0
    factor = -1;
elseif vc(3)>0
    factor = 1;
end

Theta(i) = factor*atan2d(norm(cross(vector1, vector2)), dot(vector1, vector2));

x4 = (L-1)/2 + ((L-1)/2)*cosd(Theta(i));
y4 = (L-1)/2 + ((L-1)/2)*sind(Theta(i));

plot(center(i,1),center(i,2),'o'); hold on;
plot([x1 x2],[y1 y2]); hold on;
plot([x1 x3],[y1 y3]); hold on;
plot([x3 Phi_pt_x(i)],[y3 Phi_pt_y(i)],'Color',acolor(i,:)); hold on;
%text(x4,y4,num2str(Theta(i)));
text(x4,y4,num2str(i));


P0 = [x3, y3];
P1 = [x1, y1];
P2 = [Phi_pt_x(i), Phi_pt_y(i)];
n1 = (P2 - P0) / norm(P2 - P0);  % Normalized vectors
n2 = (P1 - P0) / norm(P1 - P0);
angle3(i) = atan2(norm(det([n2; n1])), dot(n1, n2))*180/pi;

if angle3(i) < 90
    phi(i) = angle3(i);
text(x3,y3,num2str(angle3(i)));
else
    phi(i) = 180 - angle3(i);
text(x3,y3,num2str(180-angle3(i)));
end

end


figure()
histogram(phi,floor(length(phi)/16),'facecolor','r','Edgecolor','k');
xlabel('Morphological crystal orientation');
ylabel('Number of grains');
set(gca,'linewidth',2,'fontsize',15,'fontname','arial','fontweight','normal','tickdir','out','fontangle','italic');
xticks([0 30 60 90])
xlim([0 90])

Mean_phi = sum(phi)/length(phi);
std_phi = std(phi);
acolor = jet(length(center));
%=========================================================================
%=========================================================================



%=========================================================================
%=====================Calculating the crystal orientation=================
%=========================================================================
figure()
x2 = L;
y2 = (L-1)/2;

x1 = (L-1)/2;
y1 = (L-1)/2;

for i = 1:length(center)

x3 = center(i,1);
y3 = center(i,2);

vector1 = [x2,y2,0] - [x1,y1,0];
vector2 = [x3,y3,0] - [x1,y1,0];

vc = cross(vector1, vector2);

if vc(3)<0
    factor = -1;
    Theta_radial(i) = factor*atan2d(norm(cross(vector1, vector2)), dot(vector1, vector2));
elseif vc(3)>0
    factor = 1;
    Theta_radial(i) = factor*atan2d(norm(cross(vector1, vector2)), dot(vector1, vector2));
end


%--------------------------------------------------------------------------
Crystal_orien_min_random = -180;
Crystal_orien_max_random = 180;
Crystal_orien_random(i) = Theta_radial(i) + randi([Crystal_orien_min_random Crystal_orien_max_random],1,1);
%--------------------------------------------------------------------------


plot([x1 x2],[y1 y2]); hold on;
plot([x1 x3],[y1 y3]); hold on;
text(x3,y3,num2str(Theta_radial(i)));
end
axis off

Crystal_orien_radial = Theta_radial + normrnd(Radial_crystal_orien_mean,Radial_crystal_orien_std,[1 length(Theta_radial)]);

%==========================================================================
%=========================Random crystal orientation=======================
%==========================================================================

for i = 1:length(Theta_radial)
   crystal_direction_temp = Crystal_orien_random(i) - Theta_radial(i);
   crystal_direction_random(i) = crystal_direction_temp;
end


Mean_theta_random = sum((crystal_direction_random))/length(crystal_direction_random);
std_theta_random = std(crystal_direction_random);

img_crys = ones( L, L )*(0);
 for i = 1:length(crystal_direction_random)
     [temp_x,temp_y] = find(i==img(:,:));
      
     for j = 1:length(temp_y)
         img_crys(temp_x(j),temp_y(j)) = crystal_direction_random(i);
     end
 end

figure(); 
imagesc(img_crys); 
colormap(jet);
clim([-180 180]);
axis equal
axis off
h = colorbar;
h.Ticks = linspace(-180,180,7);

set(h,'FontSize',12);
set(get(h,'label'),'string','Crystal Orientation (\theta)','FontSize',13);

figure()
histogram(crystal_direction_random,floor(length(crystal_direction_random)/4),'facecolor','r','Edgecolor','k');
xlabel('Crystal orientation');
ylabel('Number of grains');
set(gca,'linewidth',2,'fontsize',15,'fontname','arial','fontweight','normal','tickdir','out','fontangle','italic');
xticks([-180 -120 -60 0 60 120 180])
xlim([-180 180])
%==========================================================================
%==========================================================================



%==========================================================================
%=======================Radial Orientation=================================
%==========================================================================
for i = 1:length(Theta_radial)
   crystal_direction_temp = Crystal_orien_radial(i)-Theta_radial(i);
   crystal_direction_radial(i) = crystal_direction_temp;

end

 Mean_crystal_direction_radial = sum((crystal_direction_radial))/length(crystal_direction_radial);
 std_crystal_direction_radial = std(crystal_direction_radial);
 
 img_crys = ones( L, L )*(0);
  for i = 1:length(crystal_direction_radial)
      [temp_x,temp_y] = find(i==img(:,:));
       
      for j = 1:length(temp_y)
          img_crys(temp_x(j),temp_y(j)) = crystal_direction_radial(i);
      end
  end
 
 figure(); 
 imagesc(img_crys); 
 colormap(jet);
 clim([-180 180]);
 axis equal
 axis off
 h = colorbar;
 h.Ticks = linspace(-180,180,7);
 set(h,'FontSize',12);
 set(get(h,'label'),'string','Radial Crystal Orientation (\theta)','FontSize',13);

figure()
histogram(crystal_direction_radial,floor(max(crystal_direction_radial)/4),'facecolor','r','Edgecolor','k');
xlabel('Radial crystal orientation');
ylabel('Number of grains');
set(gca,'linewidth',2,'fontsize',15,'fontname','arial','fontweight','normal','tickdir','out','fontangle','italic');
xticks([-180 -120 -60 0 60 120 180])
xlim([-180 180])
%==========================================================================
%==========================================================================



Data_stored = [phi;Theta_radial;Crystal_orien_random;Crystal_orien_radial;crystal_direction_random;crystal_direction_radial];
save('Data_stored.mat','Data_stored');



%==========================================================================
%==========================Meshing the domain and==========================
%==========================generating the vtk file=========================
%==========================================================================
nnx = L; % Number of nodes in x direction
nny = L; % Number of nodes in y direction

[co_or,CN,mattype,NSETS] = MeshGeneration(nnx,nny,img);
%==========================================================================
%==========================================================================



%==========================================================================
%============Finding grain boundaries and particle matrix interfaces=======
%==========================================================================
for i = 1:length(unique(mattype))
    Elemennt_CN{i} = CN(find(i == mattype),:);
end

PMI = [];
for i = 1:length(unique(mattype))-1
    PMI = [PMI; intersect(Elemennt_CN{i},Elemennt_CN{length(unique(mattype))})];
end

GB = [];
for i = 1:length(unique(mattype))-1
    for j = 1:length(unique(mattype))-1
        if i ~= j
        GB = [GB; intersect(Elemennt_CN{i},Elemennt_CN{j})];
        end
    end
end
%==========================================================================
%==========================================================================


%==========================================================================
%===================Writing the abaqus file (.inp)=========================
%==========================================================================
[Filename1] = AbaqusFileGeneration(co_or,CN,mattype,NSETS,PMI,GB);
%==========================================================================
%==========================================================================

[Filename2] = GMSHFileGeneration(co_or,CN,mattype);

numat = length(unique(mattype));

[Filename3,Domain_Scale] = Moose_Mesh_Section(numat);

lc_frac = 0.1e-6; % Length scale for interface 
[Filename5] = Sub_app_Section(numat,lc_frac);

[Filename6] = Random_Crystal_Orientation_Generation(numat,Crystal_orien_random);

[Filename7] = Radial_Crystal_Orientation_Generation(numat,Crystal_orien_radial);