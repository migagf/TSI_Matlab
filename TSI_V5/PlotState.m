function [F,NF_L,NF_R,vec,Ft,delta,Momt] = PlotState(Uw,Vw,Utr,Vtr,WheelGeom_pol,RailGeom_pol,RailProps,plotting,Vel,deltai,Momt,deltadotn_vec,U_car)
% CONTACTFORCE calculates the magnitude and components of the wheel/rail
% contact forces, using the Hertz theory and the geometry of the contact.
cr = 0.6;
% Input Parameters:
%   - Y_w (m): Lateral coordinate of the wheelset (in track coordinates)
%   - Z_w (m): Vertical coordinate of the wheelset (in track coordinates)
%   - phi_w (-): Rotation of the wheelset (in track coordinates)
%   - Wheel_geom: is a matrix containing the geometry of the wheelsets.

%% System DOF's
% Wheelset coordinates
Y_w = Uw(1);         % m
Z_w = Uw(2);         % m
phi_w = Uw(3);       % rad
R_wg = [Y_w Z_w];    % R vector, in global coordinate system

% Track coordinates
Y_tr = Utr(1);     % m
Z_tr = Utr(2);     % m
phi_tr = Utr(3);     % rad  
R_track = [Y_tr Z_tr];

R_w = R_wg-R_track; % Wheelset location in track coordinate system

% Left rail (All these values are provided in track coordinate system)
phi_Lr = phi_tr;                    % Rotation of the left rail
Z_Lr  = -3.124*25.4/1000; %-Z_tr;    % Vertical position of the left rail (m)   
Y_Lr  = -34.3355*25.4/1000; %-Y_tr;  % Horizontal position of the left rail (m)  
R_Lr = [Y_Lr Z_Lr];                 % R vector, in track coordinate system

% Right rail (All these values are provided in track coordinate system)
phi_Rr = phi_tr;                    % Rotation of the left rail
Z_Rr  = -3.124*25.4/1000; %-Z_tr;    % Vertical position of the left rail (m)   
Y_Rr  = 34.3355*25.4/1000; %-Y_tr;  % Horizontal position of the left rail (m)  
R_Rr = [Y_Rr Z_Rr];                 % R vector, in track coordinate system

%% Rail Geometry
% Rail geometry, heights measured from the center of the rail (m);
Rail_geom = 25.4/1000*[-1.3355 -1.3167 -1.1874 -0.8725 -0.6244 0.0000 0.6244 0.8725 1.1874 1.3167 1.3355;
                        2.3292  3.0933  3.4382  3.6292  3.6612 3.6751 3.6612 3.6292 3.4382 3.0933 2.3292]';

% Spline approximation
fgrid_n = 500; % Number of nodes for the fine grid
Rail_geom_fine = zeros(fgrid_n,2);

%Rail_geom_fine(:,1) = linspace(min(Rail_geom(:,1)),max(Rail_geom(:,1)),fgrid_n);

nodes = linspace(0,fgrid_n-1,fgrid_n)';
a = min(Rail_geom(:,1)); b = max(Rail_geom(:,1));
Rail_geom_fine(:,1) = 0.5*(a+b)+0.5*(b-a)*cos(nodes/fgrid_n*pi);

Rail_geom_fine(:,2) = ppval(RailGeom_pol,Rail_geom_fine(:,1)); % Evaluate the polynomial in the defined nodes
Rail_geom_fine(:,2) = -Rail_geom_fine(:,2); % Z coordinates are inverted.

%% Wheel Geometry
% Wheel_Geom is a matrix that contains points with the discretization of the
% Wheel geometry in local coordinate system. Loaded as input.

fgrid_n = 2*fgrid_n;
% Left Wheel Geometry and discretization
LWheel_geom_fine = zeros(fgrid_n,2);
LWheel_geom_fine(:,1) = linspace(min(WheelGeom_pol.breaks),max(WheelGeom_pol.breaks),fgrid_n);
%LWheel_geom_fine(:,2) = -pchip(WheelGeom(:,1),WheelGeom(:,2),LWheel_geom_fine(:,1));
LWheel_geom_fine(:,2) = -ppval(WheelGeom_pol,LWheel_geom_fine(:,1));

% Right Wheel Geometry
RWheel_geom_fine = zeros(fgrid_n,2);
RWheel_geom_fine(:,1) = -LWheel_geom_fine(:,1);
RWheel_geom_fine(:,2) = LWheel_geom_fine(:,2);

%% Initial contact velocities
deltadotn_L1 = deltadotn_vec(1); 
deltadotn_L2 = deltadotn_vec(2);
deltadotn_R1 = deltadotn_vec(3);
deltadotn_R2 = deltadotn_vec(4);

%% Left Wheel Contact
% -- Rail to track coordinate system --
phi_Lrail2track = phi_Lr - phi_tr;        % Angle between reference systems
T_Lrail2track = Tmatrix(phi_Lrail2track);    % Rotation Matrix
LRail_geom_Tr = (T_Lrail2track*Rail_geom_fine' + R_Lr')'; %Left rail in track coordinates

% -- Wheel to track coordinate system
Rlw_ws = 25.4/1000*[-34.3671 0];
LWheel_geom = LWheel_geom_fine + Rlw_ws;
R_ws2track = R_w;% - R_track;
phi_ws2track = phi_w - phi_tr;
LWheel_geom_Tr = (Tmatrix(-phi_ws2track)*LWheel_geom')' + R_ws2track;

% Calculate maximum indentation
[delta_L1,maxDz_L1,phi_cont_L1,delta_L2,maxDz_L2,phi_cont_L2] = FindContact(LWheel_geom_Tr,LRail_geom_Tr);
deltadot_L1 = delta_L1 - deltai(1);
deltadot_L2 = delta_L2 - deltai(2);

% Initial rate of indentation
deltadotn_L1 = deltadot_n(delta_L1,deltai(1),deltadot_L1,deltadotn_L1);
deltadotn_L2 = deltadot_n(delta_L2,deltai(2),deltadot_L2,deltadotn_L2);

[Nforce_L1,a_L1,b_L1] = NormalForce(delta_L1,deltadot_L1,RailProps,deltadotn_L1,cr,RailProps.R1y,RailProps.R1y);
[Nforce_L2,a_L2,b_L2] = NormalForce(delta_L2,deltadot_L2,RailProps,deltadotn_L2,cr,RailProps.R1y,RailProps.R1y);

%Forces in the global coordinate system
Ny_L1 = Nforce_L1*sin(phi_cont_L1); Ny_L2 = Nforce_L2*sin(phi_cont_L2);
Nz_L1 = -Nforce_L1*cos(phi_cont_L1); Nz_L2 = -Nforce_L2*cos(phi_cont_L2);
Mx_L1 = [Nz_L1 Ny_L1]*(maxDz_L1-R_w'); Mx_L2 = [Nz_L2 Ny_L2]*(maxDz_L2-R_w');

% Creepages and friction forces
R_L1 = maxDz_L1-R_w';
r_L1 = norm(R_L1);
psi_L1 = atan(R_L1(2)/R_L1(1)); % Angle point of contact and wheelset ??

Vy_L1 = -Vw(3)*r_L1*sin(psi_L1)+Vw(1)-Vtr(1); % Rel velocity in the Y direction in the point of contact
Vz_L1 = Vw(3)*r_L1*cos(psi_L1)+Vw(2)-Vtr(2);  % Rel velocity in the Z direction in the point of contact
Vn_L1 = Tmatrix(phi_cont_L1)*[Vy_L1 Vz_L1]'; 

% Normal plane relative velocity

Ft_L1 = creep_forces(a_L1,b_L1,Nforce_L1,Vn_L1(1),Vw(3)-Vtr(3),Vel,RailProps); % Tangential force due to creep
% Ft_L1/Nforce_L1

R_L2 = maxDz_L2-R_w';
r_L2 = norm(R_L2);
psi_L2 = atan(R_L2(2)/R_L2(1)); % Angle point of contact and wheelset ??

Vy_L2 = -Vw(3)*r_L2*sin(psi_L2)+Vw(1)-Vtr(1); % Rel velocity in the Y direction in the point of contact
Vz_L2 = Vw(3)*r_L2*cos(psi_L2)+Vw(2)-Vtr(2);  % Rel velocity in the Z direction in the point of contact

Vn_L2 = Tmatrix(phi_cont_L2)*[Vy_L2 Vz_L2]';
Ft_L2 = creep_forces(a_L2,b_L2,Nforce_L2,Vn_L2(1),Vw(3)-Vtr(3),Vel,RailProps);

Ft_Ly = Ft_L1*cos(phi_cont_L1) + Ft_L2*cos(phi_cont_L2);
Ft_Lz = Ft_L1*sin(phi_cont_L1) + Ft_L2*sin(phi_cont_L2);

Mt_L =  -Ft_L1*cos(phi_cont_L1)*R_L1(2)-Ft_L2*cos(phi_cont_L2)*R_L2(2)+...
    + Ft_L1*sin(phi_cont_L1)*R_L1(1) + Ft_L2*sin(phi_cont_L2)*R_L2(1);

% Mt_L =  -Ft_L1*cos(phi_cont_L1)*R_L1(2)+Ft_L1*sin(phi_cont_L1)*R_L1(1);

%% Right Wheel Contact
% -- Rail to track coordinate system --
phi_Rrail2track = phi_Rr - phi_tr;        % Angle between reference systems
T_Rrail2track = Tmatrix(phi_Rrail2track);    % Rotation Matrix
RRail_geom_Tr = (T_Rrail2track*Rail_geom_fine' + R_Rr')'; %Left rail in track coordinates

% -- Wheel to track coordinate system
Rrw_ws = 25.4/1000*[34.3671 0];
RWheel_geom = RWheel_geom_fine + Rrw_ws;
RWheel_geom_Tr = (Tmatrix(-phi_ws2track)*RWheel_geom')' + R_ws2track;

% Calculate maximum indentation
[delta_R1,maxDz_R1,phi_cont_R1,delta_R2,maxDz_R2,phi_cont_R2] = FindContact(RWheel_geom_Tr,RRail_geom_Tr);

deltadot_R1 = delta_R1 - deltai(3);
deltadot_R2 = delta_R2 - deltai(4);

% Initial rate of indentation
deltadotn_R1 = deltadot_n(delta_R1,deltai(3),deltadot_R1,deltadotn_R1);
deltadotn_R2 = deltadot_n(delta_R2,deltai(4),deltadot_R2,deltadotn_R2);

[Nforce_R1,a_R1,b_R1] = NormalForce(delta_R1,deltadot_R1,RailProps,deltadotn_R1,cr,RailProps.R1y,RailProps.R1y);
[Nforce_R2,a_R2,b_R2] = NormalForce(delta_R2,deltadot_R2,RailProps,deltadotn_R2,cr,RailProps.R1y,RailProps.R1y);

%Forces in the global coordinate system
Ny_R1 = Nforce_R1*sin(phi_cont_R1);
Ny_R2 = Nforce_R2*sin(phi_cont_R2);
Nz_R1 = -Nforce_R1*cos(phi_cont_R1); 
Nz_R2 = -Nforce_R2*cos(phi_cont_R2);
Mx_R1 = [Nz_R1 Ny_R1]*(maxDz_R1-R_w');
Mx_R2 = [Nz_R2 Ny_R2]*(maxDz_R2-R_w');

% Creepages and friction forces
R_R1 = maxDz_R1-R_w';
r_R1 = norm(R_R1);
psi_R1 = atan(R_R1(2)/R_R1(1));

Vy_R1 = -Vw(3)*r_R1*sin(psi_R1)+Vw(1)-Vtr(1); % Rel velocity in the Y direction in the point of contact
Vz_R1 = -Vw(3)*r_R1*cos(psi_R1)+Vw(2)-Vtr(2);  % Rel velocity in the Z direction in the point of contact
V_R1 = [Vy_R1 Vz_R1]';
Vn_R1 = Tmatrix(phi_cont_R1)*V_R1;

Ft_R1 = creep_forces(a_R1,b_R1,Nforce_R1,Vn_R1(1),Vw(3)-Vtr(3),Vel,RailProps);

R_R2 = maxDz_R2-R_w';
r_R2 = norm(R_R2);
psi_R2 = atan(R_R2(2)/R_R2(1));

Vy_R2 = -Vw(3)*r_R2*sin(psi_R2)+Vw(1)-Vtr(1); % Rel velocity in the Y direction in the point of contact
Vz_R2 = -Vw(3)*r_R2*cos(psi_R2)+Vw(2)-Vtr(2);  % Rel velocity in the Z direction in the point of contact
V_R2 = [Vy_R2 Vz_R2]';

Vn_R2 = Tmatrix(phi_cont_R1)*V_R2;
Ft_R2 = creep_forces(a_R2,b_R2,Nforce_R2,Vn_R2(1),Vw(3)-Vtr(3),Vel,RailProps);
Ft_Ry = Ft_R1*cos(phi_cont_R1) + Ft_R2*cos(phi_cont_R2);
Ft_Rz = Ft_R1*sin(phi_cont_R1) + Ft_R2*sin(phi_cont_R2);

% Data extraction during analysis
vec = [Ft_Ly Ft_Lz Ft_Ry Ft_Rz];% Tangential forces
%vec = [Nforce_L1 Nforce_L2 Nforce_R1 Nforce_R2];

Mt_R =  -Ft_R1*cos(phi_cont_R1)*R_R1(2)-Ft_R2*cos(phi_cont_R2)*R_R2(2)+...
   +Ft_R1*sin(phi_cont_R1)*R_R1(1) + Ft_R2*sin(phi_cont_R2)*R_R2(1);
% Mt_R =  -Ft_R1*cos(phi_cont_R1)*R_R1(2)+Ft_R1*sin(phi_cont_R1)*R_R1(1);

if plotting == 1
    figure(1)
    plot(LWheel_geom_Tr(:,1)+R_track(1),LWheel_geom_Tr(:,2)+R_track(2),'b',LRail_geom_Tr(:,1)+R_track(1),LRail_geom_Tr(:,2)+R_track(2),'k',RWheel_geom_Tr(:,1)+R_track(1),RWheel_geom_Tr(:,2)+R_track(2),'b',RRail_geom_Tr(:,1)+R_track(1),RRail_geom_Tr(:,2)+R_track(2),'k'), grid on, set(gca, 'YDir','reverse') %, ylim([-0.4 -0.1])%axis([-1.5 1.5 -0.5 0.0])
    draw_rectangle([0 + U_car(1), -((3.125+30+50+15)*2.54/100 + U_car(2))], 2.54/100*126, 105*2.54/100, U_car(3),[1 0 0],1)
    draw_rectangle([0 + U_car(4), -1.0 - U_car(5)], 3, 0.5, U_car(6),[1 0 0],1)
    drawnow
end

Fn = [Ny_L1+Ny_R1+Ny_L2+Ny_R2; Nz_L1+Nz_R1+Nz_L2+Nz_R2 ; Mx_L1+Mx_R1+Mx_L2+Mx_R2];
Ft = [Ft_Ry+Ft_Ly; Ft_Rz+Ft_Lz ; Mt_R+Mt_L];

F = Ft + Fn;
NF_L = Nforce_L1+Nforce_L2;
NF_R = Nforce_R1+Nforce_R2;

delta = [delta_L1;delta_L2;delta_R1;delta_R2];
Momt = [Momt, [Mx_L1 Mx_R1 Mx_L2 Mx_R2]'];

end





