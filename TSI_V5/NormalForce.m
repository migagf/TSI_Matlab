function [Nk,a,b] = NormalForce(ind,deltadot,RailProps,deltadotn,cr,RyWheel,RyRail)
% NORMALFORCE Function calculates the normal force N (MN) between two objects using
% the Hertz theory of contact.

%   Input Parameters:
%   - ind (m): Indentation between the two surfaces
%   - E1 (MPa):  Elastic modulus of the first surface  
%   - E2 (MPa):  Elastic modulus of the second surface
%   - v1 (-):    Poisson's ratio of the first surface
%   - v2 (-):    Poisson's ratio of the second surface
%   - R1x (m):   Radius of curvature for the 1st surface in x direction
%   - R1y (m):   Radius of curvature for the 1st surface in y direction
%   - R1x (m):   Radius of curvature for the 1st surface in x direction
%   - R1x (m):   Radius of curvature for the 1st surface in y direction

E1 = RailProps.E1; E2 = RailProps.E2;
v1 = RailProps.v1; v2 = RailProps.v2; 
R1x = RailProps.R1x;
R1y = RyRail;

% R1y = min(RyRail,RailProps.R1y);
% R1y = RailProps.R1y;

R2x = RailProps.R2x; 
R2y = RyWheel;

% R2y = min(RyWheel,RailProps.R2y);
% R2y = RailProps.R2y;

% Intermediate calculations:
Es = ((1-v1^2)/E1+(1-v2^2)/E2)^-1;
A = 0.5*(1/R1y+1/R2y);
B = 0.5*(1/R1x+1/R2x);

theta = rad2deg(acos(abs(B-A)/(B+A)));

r_table = [0 5 10 30 60 90 120 150 170 175 180;
           0 0.2969 0.4280 0.7263 0.9376 1 0.9376 0.7263 0.4280 0.2969 0];

m_table = [0 5 10 30 60 90 120 150 170 175 180;
           10000 11.238 6.612 2.731 1.486 1 0.7171 0.4930 0.311 0.2381 0]; 

nm_table = [0 5 10 30 60 90 120 150 170 175 180;
           0 0.0212 0.0470 0.1806 0.4826 1 2.0720 5.5380 21.26 47.20 10000];

r  = interp1(r_table(1,:),r_table(2,:),theta);
m  = interp1(m_table(1,:),m_table(2,:),theta);
nm = interp1(nm_table(1,:),nm_table(2,:),theta);
n  = m*nm;

% Normal force
K = 4*Es/3*sqrt((1/(A+B))*(1/r)^3);
%Nk = 4*Es/3*sqrt((1/(A+B))*(ind/r)^3); 
%Nk = K*(ind)^(3/2); % MN

% Damping Coefficient
% Fixed
% D = 3000; % MN*s/m

% Damping Coefficient - Hysteretic
Xi = Hysteretic_Parameter(K,cr,deltadotn);
Nk = (K+Xi*deltadot)*(ind)^(3/2);

% Dimensions of the contact patch
a = m*((3/4)*(Nk/Es)*(1/(A+B)))^(1/3);
b = n*((3/4)*(Nk/Es)*(1/(A+B)))^(1/3);

end