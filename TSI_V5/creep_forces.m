function [ft] = creep_forces(a, b, fn, vtan, omega, speed, railprops, method)
    % This function computes the tangential forces in the wheel-rail
    % interaction
    if nargin < 8
        % Method not specified in the input parameters
        method = 'kalker';
    end
    % Shear modulus of the materials
    G = railprops.E1/(2*(1 + railprops.v1)); % (MPa)

    % Friction coefficient
    mu = 0.35; 
    
    switch method
        % If using kalker equations
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