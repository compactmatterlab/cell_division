function [X,Y,Z,N_n] = periphery(N_points)
%function[output variables]= function name (input varables)

%Code from Brian Z Bentz (2020). mySphere(N) (https://www.mathworks.com/matlabcentral/fileexchange/57877-mysphere-n), MATLAB Central File Exchange. Retrieved December 10, 2020.
%Code is based on Deserno et al. Max-Planck-Institut. 2004. 
%"How to generate equidistributed points on the surface of a sphere". https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf



%% 3D Sphere with laterally equidistant points
% Assume unit radius
% Sphere Centered at (0,0,0)
% N_points: initial number of nodes
% N_n: final number of nodes (may differ from initial)
% X,Y,Z: column vectors of length N_new 
 
r_unit = 1;
Area = 4*pi*r_unit^2/N_points;     %Surface Area of each point domain
Distance = sqrt(Area);        
M_theta = round(pi/Distance); % round to nearest decimal or integer
d_theta = pi/M_theta;         
d_phi = Area/d_theta;         
N_n = 0;
for m = 0:M_theta-1
    
    Theta = pi*(m+0.5)/M_theta;
    M_phi = round(2*pi*sin(Theta)/d_phi); % not exact
    
    for n = 0:M_phi-1        
        Phi = 2*pi*n/M_phi;    
        
        N_n = N_n + 1;
        
        X(N_n) = sin(Theta)*cos(Phi);
        Y(N_n) = sin(Theta)*sin(Phi);
        Z(N_n) = cos(Theta);
        
    end
end



    