function [F,NF_L,NF_R,vec,Ft,delta,Momt,Uplift] = ContactForce(...
    Uw,Vw, ...     % Wheelset displacements
    Utr, ...       % Track displacements
    Us,Vs,...      % Track displacements and velocities 
    WheelGeom_pol, RailGeom_pol,...
    RailGeom_cur, WheelGeom_cur, RailProps,...
    plotting, Vel, deltai, Momt, deltadotn_vec,...
    Rail_geom_fine, LWheel_geom_fine, RWheel_geom_fine)

% CONTACTFORCE calculates the magnitude and components of the wheel/rail
% contact forces, using the Hertz theory and the geometry of the contact.
cr = 1.0;

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
Y_tr = Utr(1);       % m
Z_tr = Utr(2);       % m
phi_tr = Utr(3);     % rad  
R_track = [Y_tr Z_tr];


R_w = R_wg - R_track;  % Wheelset location in track coordinate system


% Left rail (All these values are provided in track coordinate system)
phi_Lr = -phi_tr - Us(3);                       % Rotation of the left rail
Z_Lr  = - 3.124 * 25.4 /1000 + Us(2) - Z_tr; %-Z_tr;    % Vertical position of the left rail (m)   
Y_Lr  = -34.3355 * 25.4/1000 + Us(1) - Y_tr; %-Y_tr;  % Horizontal position of the left rail (m)

R_Lr = [Y_Lr Z_Lr];                 % R vector, in global coordinate system

% Right rail (All these values are provided in track coordinate system)
phi_Rr = -phi_tr - Us(6);                    % Rotation of the left rail
Z_Rr  = -3.124*25.4/1000 + Us(5) - Z_tr; %-Z_tr;    % Vertical position of the left rail (m)   
Y_Rr  = 34.3355*25.4/1000 + Us(4) - Y_tr; %-Y_tr;  % Horizontal position of the left rail (m)  

R_Rr = [Y_Rr Z_Rr];                 % R vector, in track coordinate system

%% Initial contact velocities
deltadotn_L1 = deltadotn_vec(1); 
deltadotn_L2 = deltadotn_vec(2);
deltadotn_R1 = deltadotn_vec(3);
deltadotn_R2 = deltadotn_vec(4);


%% Left Wheel Contact
% -- Rail to track coordinate system --
phi_Lrail2track = phi_Lr - phi_tr;        % Angle between reference systems
T_Lrail2track = Tmatrix(phi_Lrail2track);    % Rotation Matrix
LRail_geom_Tr = (T_Lrail2track * Rail_geom_fine' + R_Lr')'; %Left rail in track coordinates


% -- Wheel to track coordinate system --
Rlw_ws = 25.4/1000*[-34.3355 0];
LWheel_geom = LWheel_geom_fine + Rlw_ws;
R_ws2track = R_w;% - R_track;
phi_ws2track = phi_w - phi_tr;
LWheel_geom_Tr = (Tmatrix(-phi_ws2track)*LWheel_geom')' + R_ws2track;


% Calculate maximum indentation
[delta_L1,maxDz_L1,phi_cont_L1,delta_L2,maxDz_L2,phi_cont_L2,D_LW] = FindContact(LWheel_geom_Tr,LRail_geom_Tr);
deltadot_L1 = delta_L1 - deltai(1);
deltadot_L2 = delta_L2 - deltai(2);


% Initial rate of indentation
deltadotn_L1 = deltadot_n(delta_L1,deltai(1),deltadot_L1,deltadotn_L1);
deltadotn_L2 = deltadot_n(delta_L2,deltai(2),deltadot_L2,deltadotn_L2);


% Calculate curvatures at the points of contact
[Rx_lw_1,Rx_lr_1] = ContactCurvature(delta_L1,maxDz_L1,R_w,Rlw_ws,R_Lr,phi_w,phi_tr,WheelGeom_cur,RailGeom_cur,'l',RailProps);
[Rx_lw_2,Rx_lr_2] = ContactCurvature(delta_L2,maxDz_L2,R_w,Rlw_ws,R_Lr,phi_w,phi_tr,WheelGeom_cur,RailGeom_cur,'l',RailProps);


[Nforce_L1,a_L1,b_L1] = NormalForce(delta_L1,deltadot_L1,RailProps,deltadotn_L1,cr,Rx_lw_1,Rx_lr_1);
[Nforce_L2,a_L2,b_L2] = NormalForce(delta_L2,deltadot_L2,RailProps,deltadotn_L2,cr,Rx_lw_2,Rx_lr_2);


%Forces in the global coordinate system
Ny_L1 =  Nforce_L1*sin(phi_cont_L1+phi_tr); Ny_L2 =  Nforce_L2*sin(phi_cont_L2+phi_tr);
Nz_L1 = -Nforce_L1*cos(phi_cont_L1+phi_tr); Nz_L2 = -Nforce_L2*cos(phi_cont_L2+phi_tr);
Mx_L1 = [Nz_L1 Ny_L1]*(maxDz_L1-R_w'); Mx_L2 = [Nz_L2 Ny_L2]*(maxDz_L2-R_w');


% Creepages and friction forces
R_L1 = maxDz_L1-R_w';
r_L1 = norm(R_L1);
psi_L1 = atan(R_L1(2)/R_L1(1)); % Angle point of contact and wheelset ??


Vy_L1 = -Vw(3)*r_L1*sin(psi_L1)+Vw(1)-Vs(1); % Rel velocity in the Y direction in the point of contact
Vz_L1 = Vw(3)*r_L1*cos(psi_L1)+Vw(2)-Vs(2);  % Rel velocity in the Z direction in the point of contact
Vn_L1 = Tmatrix(phi_cont_L1)*[Vy_L1 Vz_L1]'; 


% Normal plane relative velocity
Ft_L1 = creep_forces(a_L1,b_L1,Nforce_L1,Vn_L1(1),Vw(3)+Vs(3),Vel,RailProps); % Tangential force due to creep


% Ft_L1/Nforce_L1

R_L2 = maxDz_L2-R_w';
r_L2 = norm(R_L2);
psi_L2 = atan(R_L2(2)/R_L2(1)); % Angle point of contact and wheelset ??


Vy_L2 = -Vw(3)*r_L2*sin(psi_L2)+Vw(1)-Vs(1); % Rel velocity in the Y direction in the point of contact
Vz_L2 = Vw(3)*r_L2*cos(psi_L2)+Vw(2)-Vs(2);  % Rel velocity in the Z direction in the point of contact


Vn_L2 = Tmatrix(phi_cont_L2)*[Vy_L2 Vz_L2]';
Ft_L2 = creep_forces(a_L2,b_L2,Nforce_L2,Vn_L2(1),Vw(3)+Vs(3),Vel,RailProps);


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
Rrw_ws = 25.4/1000*[34.3355 0];
RWheel_geom = RWheel_geom_fine + Rrw_ws;
RWheel_geom_Tr = (Tmatrix(-phi_ws2track)*RWheel_geom')' + R_ws2track;

% Calculate maximum indentation
[delta_R1,maxDz_R1,phi_cont_R1,delta_R2,maxDz_R2,phi_cont_R2,D_RW] = FindContact(RWheel_geom_Tr,RRail_geom_Tr);

deltadot_R1 = delta_R1 - deltai(3);
deltadot_R2 = delta_R2 - deltai(4);

% Initial rate of indentation
deltadotn_R1 = deltadot_n(delta_R1,deltai(3),deltadot_R1,deltadotn_R1);
deltadotn_R2 = deltadot_n(delta_R2,deltai(4),deltadot_R2,deltadotn_R2);

% Calculate curvatures at the points of contact
[Rx_rw_1,Rx_rr_1] = ContactCurvature(delta_R1,maxDz_R1,R_w,Rrw_ws,R_Rr,phi_w,phi_tr,WheelGeom_cur,RailGeom_cur,'r',RailProps);
[Rx_rw_2,Rx_rr_2] = ContactCurvature(delta_R2,maxDz_R2,R_w,Rrw_ws,R_Rr,phi_w,phi_tr,WheelGeom_cur,RailGeom_cur,'r',RailProps);

[Nforce_R1,a_R1,b_R1] = NormalForce(delta_R1,deltadot_R1,RailProps,deltadotn_R1,cr,Rx_rw_1,Rx_rr_1);
[Nforce_R2,a_R2,b_R2] = NormalForce(delta_R2,deltadot_R2,RailProps,deltadotn_R2,cr,Rx_rw_2,Rx_rr_2);

%Forces in the global coordinate system
Ny_R1 = Nforce_R1*sin(phi_cont_R1+phi_tr);
Ny_R2 = Nforce_R2*sin(phi_cont_R2+phi_tr);
Nz_R1 = -Nforce_R1*cos(phi_cont_R1+phi_tr); 
Nz_R2 = -Nforce_R2*cos(phi_cont_R2+phi_tr);
Mx_R1 = [Nz_R1 Ny_R1]*(maxDz_R1-R_w');
Mx_R2 = [Nz_R2 Ny_R2]*(maxDz_R2-R_w');

% Creepages and friction forces
R_R1 = maxDz_R1-R_w';
r_R1 = norm(R_R1);
psi_R1 = atan(R_R1(2)/R_R1(1));

Vy_R1 = -Vw(3)*r_R1*sin(psi_R1)+Vw(1)-Vs(4); % Rel velocity in the Y direction in the point of contact
Vz_R1 = -Vw(3)*r_R1*cos(psi_R1)+Vw(2)-Vs(5);  % Rel velocity in the Z direction in the point of contact
V_R1 = [Vy_R1 Vz_R1]';
Vn_R1 = Tmatrix(phi_cont_R1)*V_R1;

Ft_R1 = creep_forces(a_R1,b_R1,Nforce_R1,Vn_R1(1),Vw(3)+Vs(6),Vel,RailProps);

R_R2 = maxDz_R2-R_w';
r_R2 = norm(R_R2);
psi_R2 = atan(R_R2(2)/R_R2(1));

Vy_R2 = -Vw(3)*r_R2*sin(psi_R2)+Vw(1)-Vs(4); % Rel velocity in the Y direction in the point of contact
Vz_R2 = -Vw(3)*r_R2*cos(psi_R2)+Vw(2)-Vs(5);  % Rel velocity in the Z direction in the point of contact
V_R2 = [Vy_R2 Vz_R2]';

Vn_R2 = Tmatrix(phi_cont_R1)*V_R2;
Ft_R2 = creep_forces(a_R2,b_R2,Nforce_R2,Vn_R2(1),Vw(3)+Vs(6),Vel,RailProps);
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
    plot(LWheel_geom_Tr(:,1) + R_track(1), LWheel_geom_Tr(:,2) + R_track(2),'b',...
        LRail_geom_Tr(:,1) + R_track(1), LRail_geom_Tr(:,2) + R_track(2),'k',...
        RWheel_geom_Tr(:,1) + R_track(1), RWheel_geom_Tr(:,2) + R_track(2),'b',...
        RRail_geom_Tr(:,1) + R_track(1), RRail_geom_Tr(:,2) + R_track(2),'k'), ...
        grid on, set(gca, 'YDir','reverse')%, axis([-2 2 -0.5 0.2])
    drawnow
end

Fn = [Ny_L1 + Ny_L2, Ny_R2 + Ny_R1; 
      Nz_L1 + Nz_L2, Nz_R2 + Nz_R1; 
     (Mx_L1 + Mx_L2), (Mx_R2 + Mx_R1)];

Ft = [Ft_Ry, Ft_Ly; 
      Ft_Rz, Ft_Lz;
      Mt_R, Mt_L];

%Fn(1:2,1) = Tmatrix(-phi_tr) * Fn(1:2,1);
%Ft(1:2,1) = Tmatrix(-phi_tr) * Ft(1:2,1);

F = Ft + Fn;

NF_L = Nforce_L1 + Nforce_L2;
NF_R = Nforce_R1 + Nforce_R2;

delta = [delta_L1;delta_L2;delta_R1;delta_R2];
Momt = [phi_cont_L1 phi_cont_L1 phi_cont_L1 phi_cont_L1]';

Uplift = [D_LW , D_RW];

end


function [T] = Tmatrix(phi)
%Tmatrix creates the transformation matrix for 2 coordinate systems rotated
%in phi radians
%   phi is the angle between the 2 coordinate systems
%   T is the transformation matrix
T = [cos(phi) sin(phi); -sin(phi) cos(phi)];
end





