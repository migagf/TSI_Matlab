function [d] = UpliftValue(WheelGeom,RailGeom)
n_points = 20;
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

if isempty(X_contact) % If there is no contact
    
    % Calculate distance
    d = sqrt((X_con - X_con.').^2 + (Z_wheel - Z_rail.').^2);
    d = d + diag(inf(length(X_con),1));
    d = min(min(d));

else % If there is no contact
    d = 0;

end

end