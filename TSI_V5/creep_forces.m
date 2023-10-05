function [Ft] = creep_forces(a,b,Fn,Vn,omega,Vel,RailProps)

G = RailProps.E1/(2*(1+RailProps.v1)); %MPa
mu = 0.3; % Friction Coefficient

method = 'carter';

switch method
    case 'kalker'
        c = sqrt(a*b);

        if c ~= 0
            vy = Vn/Vel; % Creepage
            spin = omega/Vel;
            
            c22 = 2.4014 + 2.3179/(b/a) - 0.02/((b/a)^2);
            c23 = 0.4147 + 1.0184/(b/a) + 0.065/((b/a)^2) - 0.0013/((b/a)^3);
            
            % Tangential force
            Ft = -G*a*b*c22*vy + G*a*b*c23*c*spin;

            if abs(Ft) > mu*Fn % Saturation
                Ft = -sign(vy)*mu*Fn;
            end

        else
            Ft = 0;
        end

    case 'carter'
        vy = Vn/Vel;
        c = sqrt(a*b);

        R = (1/RailProps.R2x+1/RailProps.R1y);
        if c ~= 0 && vy ~= 0
            if abs(R*vy/(mu*b)) >= 1
                Ft = -sign(vy)*mu*Fn;
            else
                Ft = -(2+R/(mu*b)*abs(vy))*R*vy/(mu*b)*mu*Fn;
            end
        else
            Ft = 0;
        end
        
    
end
end