function [ft] = creep_forces(a, b, fn, vtan, omega, speed, railprops, method)
% This function computes the tangential forces in the wheel-rail
% interaction using Kalker or Carter simplified models


% (1) Check if we have enough input parameters...
if nargin < 8
    % Method not specified in the input parameters, using Kalker as
    % default
    method = 'kalker';
end

% (2) Shear modulus of the materials
G = railprops.E1/(2*(1 + railprops.v1)); % (MPa)

% (3) Friction coefficient (update: make an input parameter)
mu = 0.35;

% (4) Method selector (Kalker or Carter equations)
switch method
    case 'kalker'
        c = sqrt(a*b); % Average size of the contact patch
        if c ~= 0
            % Linear and rotational creepages
            vy = vtan/speed;
            spin = omega/speed;
            
            c22 = 2.4014 + 2.3179/(b/a) - 0.02/((b/a)^2);
            c23 = 0.4147 + 1.0184/(b/a) + 0.065/((b/a)^2) - 0.0013/((b/a)^3);
            
            % Tangential force
            ft = -G*a*b*c22*vy + G*a*b*c23*c*spin;
            
            % If we reach saturation
            if abs(ft) > mu*fn 
                ft = -sign(vy)*mu*fn;
            end

        else
            % No contact (a or b or both are 0).
            ft = 0;
        end

    case 'carter'
        vy = vtan/speed;
        c = sqrt(a*b);

        R = (1/railprops.R2x + 1/railprops.R1y);

        if c ~= 0 && vy ~= 0

            % If we reach saturation, then use coulomb friction
            if abs(-(2 + R/(mu*b)*abs(vy))*R*vy/(mu*b)) >= 1
                ft = -sign(vy)*mu*fn;

            % Carter's equation
            else 
                ft = -(2+R/(mu*b)*abs(vy))*R*vy/(mu*b)*mu*fn;
            end
        else
            ft = 0;
        end
end
end