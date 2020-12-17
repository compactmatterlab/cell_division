% Code written by Nadia Caroline Beydoun
% Cellular Rounding of a Single Cell During Mitosis
% Date last edited: 12/10/2020

clc; clear all;
tic
global dist0 Fmagpr storedistc storeFcph storeperipheryforce storecentforce cph gammax gammay gammaz resize Fmagk rad N_p Fmagd Fmagcr no_o no_c no_d no_k no_cr o spe spa ac acr  N_points pcell Fmagat Fmagatt  DT volume1 c2cd Fmagce Tfb Fmagpa Fmagpe storeFc2p storeFrepc storeFrepcr storeFcr2c storeFrpa storeFrpe storeFatp storeFcont storeF_pe storeSA storeSAt
rng('shuffle')

%% Boundary conditions 
 
%periphery points
N_points = 50; % test

%number of dynein   
no_d = 240;

%number of kinesin
no_k = 160;

%number of Chromosomes
no_cr = 23;    

%number of centrosomes
no_c = 2;    

%periphery sphere radius(microns)
rad=5;

%Centrosome Sphere Radius (microns) 
ac=1.8;                  

%Chromosome Sphere Radius (microns)
acr=1.2;                   

%smallest displacement between chromosomes (perpandicular to centrosomes)(microns)  
spe=0.6; 

%smallest displacement between chromosomes (parallel to centrosomes)(microns)                 
spa=0.5;

%smallest distance between centrosomes
c2cd=4;% microns

%coverplate height
cph=5.0; 

%Damping Constant
gammax=20000;                 % pN*s/um (Caragine et al., 2018)
gammay=20000;
gammaz=20000;

%% Force Magnitudes 

%Intracellular Pressure of cell (pN/um^2)
pcell=13;     

%Force Magnitude Kinesin (pN)            
Fmagk=10;                  % Force of Kinesin =6pN (Belyy et al. 2016)

%Force Magnitude Dynein (pN)    
Fmagd=8;                   % Force of dynein =4.3pN (Belyy et al. 2016)

%Force Magnitude Periphery attraction to plate  (pN) 
Fmagat=0.00000;

%Contractile force magnitude periphery points(pN)
Fmagatt=300;

%Force magnitude Chromosomes onto centrosomes  (pN)
Fmagcr=11;                  
           
%Force Magnitude Repulsion centrosomes to centrosomes (pN)
Fmagce=100;

%Force Magnitude Repulsion chromosome to chromosome in parallel (pN)
Fmagpa=10;

%Force Magnitude Repulsion chromosome to chromosome in perpandicular (pN)
Fmagpe=100;

%Force Magnitude Repulsion chromosome to periphery (pN)
Fmagpr=100;

%% Equidistant lateral periphery points Solver (periphery.m)
%Code from Brian Z Bentz (2020). mySphere(N) (https://www.mathworks.com/matlabcentral/fileexchange/57877-mysphere-n), MATLAB Central File Exchange. Retrieved December 10, 2020.

[X,Y,Z,N_n] = periphery(N_points);   %(calling periphery function)

switch 1      %scaling the radius and diplacement
    case 1
        % size of radius
        X = X.*rad;
        Y = Y.*rad;
        Z = Z.*rad;
        
    case 2
        % positions the sphere
        X = X + 0;
        Y = Y + 0;
        Z = Z + 0;        
end

nn=[ X' Y' Z'];
N_pp=size(nn,1);

%% Initial Triangulation of sphere code from mathworks.com /help/matlab/ref/triangulation.vertexnormal.html?s_tid=doc_ta

DT=delaunayTriangulation(nn(:,1),nn(:,2),nn(:,3));

[Tfb,Xfb] = freeBoundary(DT);           % triangulate the boundary (surface) 

%% Distance between points on the sphere
for n=1:N_pp
    for m=1:N_pp
        
        Np2Np2(n,m,:)=(nn(n,:)-nn(m,:)); %Np2Np2 are the position vectors between each periphery point
        %Normalize (Np2Np)
        test0(:)=Np2Np2(n,m,:);
        dist0(n,m)=norm(test0);          % dist0= distance between periphery points on the sphere
        
    end
end


%% Hemisphere Creation
for iii=1:N_n
    if Z(iii)<0
    Z(iii)=0;
    end
end


%% Periphery Coordinates

zz=[ X' Y' Z'];
z0=zz(:);
N_p=size(zz,1);            % number of periphery points after periphery function and hemisphere creation. (can vary from N_points)

%% Convex Hull (volume of initial hemisphere)

[K1,volume1] = convhulln(zz);

%% Centroid of Hemisphere

Mcent=[mean(zz(:,1)),mean(zz(:,2)),mean(zz(:,3))];

%% Point of Origin, same as centroid

no_o=1;                  
o = zeros(no_o,3);
o=o+Mcent;

%% Positions Centrosomes and Chromosomes and Origin 
%random numbers within sphere code from https://www.mathworks.com/help/matlab/math/numbers-placed-randomly-within-volume-of-sphere.html

%Centrosome positions (microns)

rvals = 2*rand(no_c,1)-1; 
elevation = asin(rvals); 
azimuth = 2*pi*rand(no_c,1); 
[xr,yr,zr] = sph2cart(azimuth,elevation,ac);
c=[ xr+Mcent(:,1), yr+Mcent(:,2) zr+Mcent(:,3)];

%Chromosome positions (microns)

ravals = 2*rand(no_cr,1)-1; 
elevations = asin(ravals); 
azi = 2*pi*rand(no_cr,1); 
radiii = acr*(rand(no_cr,1).^(1/3)); 
unit=1;
[xa,ya,za] = sph2cart(azi,elevations,radiii);
cr= [xa+Mcent(:,1), ya+Mcent(:,2), za+Mcent(:,3)];




%% Original Starting positions of each point chromosomes, centrosomes, centroid and periphery

startingposcr=cr;     
startingposc=c;
startingposp=zz;          
Mcentroid=Mcent;

%% Size of matrix

v=vertcat(c, cr, zz);      % Concatenates matrices vertically
v0=v(:);                   % Colon operator
resize = size(v);          


%% Time Interval

t0=0;                      %Initial time (seconds)                      
%tf=1800;                   %Final time (seconds)
tf=30;

% time step interval(seconds)
deltat=10; 

tspan=t0:1:tf;            % interval for results in ode45 https://www.mathworks.com/matlabcentral/answers/92961-how-do-i-use-a-fixed-step-size-with-ode23-and-ode45-in-matlab
 


%% Storage for all forces
storeperipheryforce=NaN(4,tf+1,N_p,3); %(maximum number of forces on each periphery point,final time, N_p, xyz)
storecentforce=NaN(4,tf+1,no_c,3);
storeFc2p=NaN(no_c,N_p,3);
storeFcr2c=NaN(no_c,no_cr,3);
storeFrepc=NaN(no_c,no_c,3);
storeFrepcr=NaN(no_cr,3);
storeFrpa=NaN(no_cr,no_cr,3);
storeFrpe=NaN(no_cr,no_cr,3);
storeFatp=NaN(N_p,3);
storeFcont=NaN(N_p,3);
storeF_pe=NaN(N_p,3);
storeFcph=NaN(N_p,3);


%% Solver with Movie (using dynamics function)
for total_t = deltat:deltat:tf
    
[t,v1]=ode45('dynamics',(t0:total_t),v0);  % calling dynamics function

t0=total_t;                                % resetting start point. Going from previous step to the next step. 

v0(:)=v1(length(v1(:,1)),:)';              % v0 is every 10(deltat) seconds
Vstore(:,total_t/deltat)=v0(:);            % Vstore(:,total_t/deltat) increases by each iteration total_t/deltat 

%Vstored(:,t0/deltat:total_t/deltat)=v1(:); %Stores every iteration
end

%% reshape matrix
v = reshape(Vstore(:,length(Vstore(1,:))),resize); % reshape(first size, final size)   

%% Solver without Movie (using dynamics function)

%[t,v1]=ode45('dynamics2_V2_02_11_2020',tspan,v0);  % calling dynamics function 
%v=reshape(v1(size(v1,1),:)); 

                         
%% Function    

c(1:no_c,:)=v(1:no_c,:);                                
cr(1:no_cr,:)=v(no_c+1:no_c+no_cr,:);                   
zz(1:N_p,:)=v(no_c+no_cr+1:no_c+no_cr+N_p,:);

N_p=size(zz,1);

%Centroid
Mcent=[mean(zz(:,1)),mean(zz(:,2)),mean(zz(:,3))];  

%Convex Hull/ Final volume of cell
[K2,volume2] = convhulln(zz);  
                               
                               
%P cell ratio
P_cell=pcell*volume1/volume2;   % multiply pcell constant by ratio (initial volume)/ (final volume) 


% Final Triangulation
TR = triangulation(Tfb,zz);  

%Unit normal vectors to the vertices
VN = vertexNormal(TR);          % mathworks.com/help/matlab/ref/triangulation.vertexnormal.html?s_tid=doc_ta           


%% Plot of starting positions

figure
hold on

plot3(startingposp(:,1), startingposp(:,2), startingposp(:,3), 'rx');

plot3(startingposcr(:,1), startingposcr(:,2), startingposcr(:,3), 'c*');

plot3(startingposc(:,1), startingposc(:,2), startingposc(:,3), 'g+');

plot3(o(:,1), o(:,2), o(:,3), 'ms');

axis([-10 10 -10 10 -10 10])

title('Initial Positions')

savefig('InitialPositions')

hold off

%% Plot of Forces on centrosomes
figure
hold on
for tt=1:tf
    for m=2%centrosome 1
        
        xt(:)=storecentforce(1,tt,m,:); %Fc2p
        xy=sqrt(xt(1)^2+xt(2)^2+xt(3)^2);
        scatter(tt/60,xy,'mo');
        
        xt(:)=storecentforce(2,tt,m,:); %Fcr2c
        xy=sqrt(xt(1)^2+xt(2)^2+xt(3)^2);
        scatter(tt/60,xy,'ro');
        
        xt(:)=storecentforce(3,tt,m,:); %Frepc
        xy=sqrt(xt(1)^2+xt(2)^2+xt(3)^2);
        scatter(tt/60,xy,'ko');
    end
end

title('Forces on Centrosomes')

savefig('ForcesCentrosomes')

hold off

%% Plot of Forces on periphery 
figure
hold on
for tt=1:tf
    for nn=1  %point 13
        xx(:)=storeperipheryforce(1,tt,nn,:); %Fatp
        xxy=sqrt(xx(1)^2+xx(2)^2+xx(3)^2);
        scatter(tt/60,xxy,'mo');
            
        xx(:)=storeperipheryforce(2,tt,nn,:); %Fcont
        xxy=sqrt(xx(1)^2+xx(2)^2+xx(3)^2);
        scatter(tt/60,xxy,'ro');
        
        xx(:)=storeperipheryforce(3,tt,nn,:); %F_pe
        xxy=sqrt(xx(1)^2+xx(2)^2+xx(3)^2);
        scatter(tt/60,xxy,'ko');
        
        xx(:)=storeperipheryforce(4,tt,nn,:); %Fcph
        xxy=sqrt(xx(1)^2+xx(2)^2+xx(3)^2);
        scatter(tt/60,xxy,'go');
        
    end
end
title('Forces on Periphery')
savefig('ForcesPeriphery')
hold off

%% Force Vectors on Periphery
figure
hold on

plot3(zz(:,1), zz(:,2), zz(:,3), 'rx'); 

plot3(cr(:,1), cr(:,2), cr(:,3), 'c*');

plot3(c(:,1), c(:,2), c(:,3), 'g+'); 

plot3(o(:,1), o(:,2), o(:,3), 'ms');

quiver3(zz(:,1), zz(:,2), zz(:,3), -storeFc2p(1,:,1)', -storeFc2p(1,:,2)', -storeFc2p(1,:,3)','b') %blue= Force vectors centrosome 1

quiver3(zz(:,1), zz(:,2), zz(:,3), -storeFc2p(2,:,1)', -storeFc2p(2,:,2)', -storeFc2p(2,:,3)','k') %black= Force vectors centrosome 2

title('Force Vectors Fc2p')
savefig('ForceVectorsFc2p')
hold off

%% centrosomes and chromosomes starting and ending positions
hold on

figure                                                    %shows before and after of chromosome position 
plot3(cr(:,1), cr(:,2), cr(:,3), 'c*');
hold on
plot3(startingposcr(:,1),startingposcr(:,2),startingposcr(:,3),'mo');

plot3(c(:,1), c(:,2), c(:,3), 'k*');

plot3 (startingposc(:,1), startingposc(:,2), startingposc(:,3), 'rx');

title('starting and ending position of chromosomes and centrosomes')
savefig('chromosomecentrosome') 
hold off


%% Force Vectors from Pressure
figure
hold on 
quiver3(zz(:,1),zz(:,2),zz(:,3),storeF_pe(:,1),storeF_pe(:,2),storeF_pe(:,3), 'Color','m');   %vector of pressure in pink
plot3(zz(:,1),zz(:,2),zz(:,3), 'r*');

trisurf(TR,'FaceColor',[0.8 0.8 1.0]);

axis equal
% quiver3(Xfb(:,1),Xfb(:,2),Xfb(:,3), ...
%      VN(:,1),VN(:,2),VN(:,3),0.5,'Color','b');  %vector normal to triangulation 
plot3(TR.Points(:,1), TR.Points(:,2), TR.Points(:,3),'ks') %why TR.Points?
plot3(cr(:,1), cr(:,2), cr(:,3), 'c*');
plot3(c(:,1), c(:,2), c(:,3), 'm+');
plot3(o(:,1), o(:,2), o(:,3), 'gs');

alpha 0.2   %adds transparency 

title ('Force vectors from pressure')
savefig('rounding')
hold off

%% Initial and Final triangulation
figure
%Plot of traingulation and unit normal vectors and periphery points  : mathworks.com/help/matlab/ref/triangulation.vertexnormal.html?s_tid=doc_ta
trisurf(TR,'FaceColor',[0.8 0.8 1.0]);
axis equal
hold on
quiver3(Xfb(:,1),Xfb(:,2),Xfb(:,3), ...
    VN(:,1),VN(:,2),VN(:,3),0.5,'Color','b');
plot3(TR.Points(:,1), TR.Points(:,2), TR.Points(:,3),'r*') %why TR.Points?

%Volume of traingulation Convex hull of triangulation. Code from https://www.mathworks.com/help/matlab/ref/delaunaytriangulation.convexhull.html
[CVH,vol]=convexHull(DT);  % vol is the volume in microns of the shape

trisurf(CVH,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3), ...
    'FaceColor','cyan')
alpha 0.2   %adds transparency

title ('initial and final triangulation')
savefig('triangulation')
hold off
   
%% Plot final cell 
figure
hold on 

plot3(zz(:,1),zz(:,2),zz(:,3), 'r*');

trisurf(TR,'FaceColor',[0.8 0.8 1.0]);

axis equal
hold on

plot3(TR.Points(:,1), TR.Points(:,2), TR.Points(:,3),'ks') %why TR.Points?
plot3(cr(:,1), cr(:,2), cr(:,3), 'c*');
plot3(c(:,1), c(:,2), c(:,3), 'm+');
plot3(Mcent(:,1), Mcent(:,2), Mcent(:,3), 'gs');
   
alpha 0.2   %adds transparency 

title ('Final Cell Shape')
savefig('Finalcell')
hold off

%% Periphery Force Diagram
figure
hold on
axis equal
plot3(zz(:,1), zz(:,2), zz(:,3), 'rx'); 
plot3(cr(:,1), cr(:,2), cr(:,3), 'c*');
plot3(c(:,1), c(:,2), c(:,3), 'k+'); 
plot3(Mcent(:,1), Mcent(:,2), Mcent(:,3), 'ks');

%Force Vectors on periphery

quiver3(zz(:,1), zz(:,2), zz(:,3), storeF_pe(:,1),storeF_pe(:,2),storeF_pe(:,3), 'Color','b');
quiver3(zz(:,1), zz(:,2), zz(:,3), storeFcont(:,1),storeFcont(:,2),storeFcont(:,3), 'Color',[0.8500, 0.3250, 0.0980]); %orange
quiver3(zz(:,1), zz(:,2), zz(:,3), storeFatp(:,1),storeFatp(:,2),storeFatp(:,3), 'Color','k'); 
quiver3(zz(:,1), zz(:,2), zz(:,3), storeFcph(:,1),storeFcph(:,2),storeFcph(:,3), 'Color', [0, 0.5, 0]); %green
title('Force Diagram')

savefig('ForceDiagram')

hold off

%% Intracellular Force Diagram
figure
hold on

axis equal
plot3(zz(:,1), zz(:,2), zz(:,3), 'kx'); 
plot3(cr(:,1), cr(:,2), cr(:,3), 'c*');
plot3(c(:,1), c(:,2), c(:,3), 'm+'); 
plot3(Mcent(:,1), Mcent(:,2), Mcent(:,3), 'gs');

%Force Vectors on periphery

myquiv=quiver3(zz(:,1), zz(:,2), zz(:,3), -storeFc2p(1,:,1)', -storeFc2p(1,:,2)', -storeFc2p(1,:,3)','Color',[0, 0.5, 0]); %green= Force vectors centrosome 1
myquiv2=quiver3(zz(:,1), zz(:,2), zz(:,3), -storeFc2p(2,:,1)', -storeFc2p(2,:,2)', -storeFc2p(2,:,3)','Color',[0.8500, 0.3250, 0.0980]); %orange= Force vectors centrosome 2
quiver3(cr(:,1), cr(:,2), cr(:,3), storeFcr2c(1,:,1)', storeFcr2c(1,:,2)', storeFcr2c(1,:,3)','Color',[0.4940, 0.1840, 0.5560]); %purple
quiver3(cr(:,1), cr(:,2), cr(:,3), storeFcr2c(2,:,1)', storeFcr2c(2,:,2)', storeFcr2c(2,:,3)','Color',[0.4940, 0.1840, 0.5560]); %purple
myquiv.ShowArrowHead = 'off';
myquiv2.ShowArrowHead = 'off';

title('Intracellular Force Diagram')

savefig('Intracellularforcediagram')

hold off


%% MOVIE of Final Cell
Fmov=VideoWriter('video');
Fmov.Quality=95;
open(Fmov)
for moviet=1:length(Vstore(1,:))
   
v = reshape(Vstore(:,moviet),resize);


c(1:no_c,:)=v(1:no_c,:);                                
cr(1:no_cr,:)=v(no_c+1:no_c+no_cr,:);                   
zz(1:N_p,:)=v(no_c+no_cr+1:no_c+no_cr+N_p,:);

N_p=size(zz,1);

%Centroid
Mcent=[mean(zz(:,1)),mean(zz(:,2)),mean(zz(:,3))];  

% Initial Triangulation
TR = triangulation(Tfb,zz);  

%Unit normal vectors to the vertices
VN = vertexNormal(TR);     

                                
if rem(moviet,1)==0 %looks at every time step

%figure             %uncomment if you want to capture each frame 

hold off 
%quiver3(zz(:,1),zz(:,2),zz(:,3),storeF_pe(:,1),storeF_pe(:,2),storeF_pe(:,3), 'Color','m');   %uncomment if you want to plot force due to pressure  
plot3(zz(:,1),zz(:,2),zz(:,3), 'r*');

hold on

plot3(cr(:,1), cr(:,2), cr(:,3), 'c*');
plot3(c(:,1), c(:,2), c(:,3), 'm+');
plot3(Mcent(:,1), Mcent(:,2), Mcent(:,3),'ks');

trisurf(TR,'FaceColor',[0.8 0.8 1.0]);

axis equal
 
plot3(TR.Points(:,1), TR.Points(:,2), TR.Points(:,3),'ks') 

alpha 0.2   %adds transparency 

axis([-10 10 -10 10 -10 10])

view([180,10])

text(9,9,-9,strjoin([string(moviet*deltat/60),'mins'])); % (12,12,-12) position text is found. moviet=iterations, deltat=time step, /60 makes it in minutes. moviet*deltat=seconds passed. string=line of text instead of a number. 

title ('force from pressure')

frame=getframe;

writeVideo(Fmov,frame);

end


end

% %% MOVIE of Force Vectors
% Fmov=VideoWriter('video');
% Fmov.Quality=95;
% open(Fmov)
% for moviet=1:length(Vstore(1,:))
%    
% v = reshape(Vstore(:,moviet),resize);
% 
% 
% c(1:no_c,:)=v(1:no_c,:);                                %c(2 centrosomes,xyz)=v(2 centrosomes,xyz) creating a 2x3
% cr(1:no_cr,:)=v(no_c+1:no_c+no_cr,:);                   %creating a 3x5 matrix of cr in v
% zz(1:N_p,:)=v(no_c+no_cr+1:no_c+no_cr+N_p,:);
% 
% N_p=size(zz,1);
% 
%                                      
% if rem(moviet,1)==0 %looking at every time step
% 
%  figure
% %Plot of final traingulation and unit normal vectors and periphery points  : mathworks.com/help/matlab/ref/triangulation.vertexnormal.html?s_tid=doc_ta
% hold off
% %quiver3(zz(:,1),zz(:,2),zz(:,3),storeF_pe(:,1),storeF_pe(:,2),storeF_pe(:,3), 'Color','m');   %vector of pressure in pink
% hold on
% plot3(zz(:,1), zz(:,2), zz(:,3), 'kx'); 
% plot3(cr(:,1), cr(:,2), cr(:,3), 'c*');
% plot3(c(:,1), c(:,2), c(:,3), 'm+'); 
% plot3(Mcent(:,1), Mcent(:,2), Mcent(:,3), 'gs');
% 
% %Force Vectors on periphery
% 
% myquiv=quiver3(zz(:,1), zz(:,2), zz(:,3), -storeFc2p(1,:,1)', -storeFc2p(1,:,2)', -storeFc2p(1,:,3)','Color',[0, 0.5, 0]); %green= Force vectors centrosome 1
% myquiv2=quiver3(zz(:,1), zz(:,2), zz(:,3), -storeFc2p(2,:,1)', -storeFc2p(2,:,2)', -storeFc2p(2,:,3)','Color',[0.8500, 0.3250, 0.0980]); %orange= Force vectors centrosome 2
% myquiv3=quiver3(cr(:,1), cr(:,2), cr(:,3), storeFcr2c(1,:,1)', storeFcr2c(1,:,2)', storeFcr2c(1,:,3)','Color',[0.4940, 0.1840, 0.5560]); %purple
% myquiv4=quiver3(cr(:,1), cr(:,2), cr(:,3), storeFcr2c(2,:,1)', storeFcr2c(2,:,2)', storeFcr2c(2,:,3)','Color',[0.4940, 0.1840, 0.5560]); %purple
% myquiv.ShowArrowHead = 'off';
% myquiv2.ShowArrowHead = 'off';
% myquiv3.ShowArrowHead = 'off';
% myquiv4.ShowArrowHead = 'off';
% 
% % quiver3(zz(:,1), zz(:,2), zz(:,3), storeF_pe(:,1),storeF_pe(:,2),storeF_pe(:,3), 'Color','b');
% % quiver3(zz(:,1), zz(:,2), zz(:,3), storeFcont(:,1),storeFcont(:,2),storeFcont(:,3), 'Color',[0.8500, 0.3250, 0.0980]); %orange
% % quiver3(zz(:,1), zz(:,2), zz(:,3), storeFatp(:,1),storeFatp(:,2),storeFatp(:,3), 'Color','k'); 
% % quiver3(zz(:,1), zz(:,2), zz(:,3), storeFcph(:,1),storeFcph(:,2),storeFcph(:,3), 'Color', [0, 0.5, 0]); %green
% 
% axis equal
% hold on
% 
% 
%    
% alpha 0.2   %adds transparency 
% 
% axis([-10 10 -10 10 -10 10])
% 
% view([180,10])
% 
% text(9,9,-9,strjoin([string(moviet*deltat/60),'mins'])); % (12,12,-12) position text is found. moviet=iterations, deltat=time step, /60 makes it in minutes. moviet*deltat=seconds passed. string=line of text instead of a number. 
% 
% title ('force from pressure')
% 
% %close(figure)
% 
% frame=getframe;
% 
% writeVideo(Fmov,frame);
% 
% end
% 
% 
% end

close(Fmov)
save('zz.mat','zz','Vstore','volume1','volume2','N_p','storeF_pe','storeFatp','storeFc2p','storeFcont','storeFcr2c','storeFrepc','storeFrepcr','storeFrpa','storeFc2p','storeFrpe','storeperipheryforce','storedistc','storeFcph','storeSAt','storeSA','cr','c')

toc/60/60/60
