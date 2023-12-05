function [delta_z1,maxDz1,phi_cont1,delta_z2,maxDz2,phi_cont2,d] = FindContact(WheelGeom,RailGeom)
n_points = 50;

% Establish search limits
L_lim = max(min(WheelGeom(:,1)),min(RailGeom(:,1)));
R_lim = min(max(WheelGeom(:,1)),max(RailGeom(:,1)));
X_con = linspace(L_lim,R_lim,n_points); % Discretization of the contact area
Z_wheel = interp1(WheelGeom(:,1),WheelGeom(:,2),X_con);
Z_rail = interp1(RailGeom(:,1),RailGeom(:,2),X_con);

dZ = Z_rail - Z_wheel; % Vertical indentation

% Calculate points of minimal vertical distance
minima = dZ(1:end-1).*dZ(2:end);
X_contact = X_con(find(minima<0)+1);
Z_contact = Z_wheel(find(minima<0)+1);

R_wheel_con = [X_con',Z_wheel'];
R_rail_con = [X_con',Z_rail'];

if length(X_contact)==2 % If there is 1 contact patch
    
    X_cont1 = X_contact(1:2); Z_cont1 = Z_contact(1:2);
    
    % Inclination of the contact plane
    phi_cont1 = atan((Z_cont1(2)-Z_cont1(1))/(X_cont1(2)-X_cont1(1)));
    
    % semiaxis of contact point

    a1 = 0.5*abs(X_cont1(2)-X_cont1(1))/cos(phi_cont1);

    R_wheel_con_perp1 = (Tmatrix(phi_cont1)*(R_wheel_con' - [X_cont1(1);Z_cont1(1)]))';
    R_rail_con_perp1 = (Tmatrix(phi_cont1)*(R_rail_con' - [X_cont1(1);Z_cont1(1)]))';
    
    dZ_perp1 = R_wheel_con_perp1(:,2) - R_rail_con_perp1(:,2);
    [delta_z1 , X_coord1] = max(dZ_perp1);

    % Point of application of the force (At the point of max indentation)
    maxDz1 = Tmatrix(-phi_cont1)*[R_rail_con_perp1(X_coord1,1),0]'+[X_cont1(1);Z_cont1(1)];
    
    delta_z2 = 0;
    maxDz2 = [0;0];
    phi_cont2 = 0;
    %  disp('1 contact patches')

    % Distance between wheel and rail
    d = 0;

elseif length(X_contact)==4
    % Calculations for the first point of contact
    
    X_cont1 = X_contact(1:2); Z_cont1 = Z_contact(1:2);
    
    % Inclination of the contact plane
    phi_cont1 = atan((Z_cont1(2)-Z_cont1(1))/(X_cont1(2)-X_cont1(1)));
    
    R_wheel_con_perp1 = (Tmatrix(phi_cont1)*(R_wheel_con' - [X_cont1(1);Z_cont1(1)]))';
    R_rail_con_perp1 = (Tmatrix(phi_cont1)*(R_rail_con' - [X_cont1(1);Z_cont1(1)]))';
    
    dZ_perp1 = R_wheel_con_perp1(:,2) - R_rail_con_perp1(:,2);
    [delta_z1 , X_coord1] = max(dZ_perp1);
    
    % Point of application of the force (At the point of max indentation)
    maxDz1 = Tmatrix(phi_cont1)*[R_rail_con_perp1(X_coord1,1),R_wheel_con_perp1(X_coord1,2)]'+[X_cont1(1);Z_cont1(1)];
    
    % Calculations for the second point of contact

    X_cont2 = X_contact(3:4); 
    Z_cont2 = Z_contact(3:4);
    
    % Inclination of the contact plane
    phi_cont2 = atan((Z_cont2(2)-Z_cont2(1))/(X_cont2(2)-X_cont2(1)));
    
    % a2 = 0.5*abs(X_cont2(2)-X_cont2(1))/cos(phi_cont2);

    R_wheel_con_perp2 = (Tmatrix(phi_cont2)'*(R_wheel_con' - [X_cont2(1);Z_cont2(1)]))';
    R_rail_con_perp2 = (Tmatrix(phi_cont2)'*(R_rail_con' - [X_cont2(1);Z_cont2(1)]))';
    
    dZ_perp2 = R_wheel_con_perp2(:,2) - R_rail_con_perp2(:,2);
    [delta_z2 , X_coord2] = max(dZ_perp2);
    
    % Point of application of the force (At the point of max indentation)
    maxDz2 = Tmatrix(-phi_cont2)*[R_rail_con_perp2(X_coord2,2),0]'+[X_cont2(1);Z_cont2(1)];
%     disp('2 contact patches')
    
    % Distance between wheel and rail
    d = 0;

else % If there is no contact
    % Find Smallest Distance between Wheel and Rail
    d = sqrt((X_con - X_con.').^2 + (Z_wheel - Z_rail.').^2);
    d = d + diag(inf(length(X_con),1));
    d = min(min(d));

    delta_z1 = 0;
    maxDz1 = [0;0];
    phi_cont1 = 0;

    delta_z2 = 0;
    maxDz2 = [0;0];
    phi_cont2 = 0;
end

end