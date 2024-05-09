clear, clc, close all

% Define the Ground Motion
tt = 0:0.01:5;
g = 9.81;
ugddot = 0.5 * g * sin(10 * tt);

% plot(t, ugddot)

yprime = @(t, y) dydt(y, t, tt, ugddot);
y0 = [0, 0, 0, 0]';

tspan = [0, 5];

[t, y] = ode45(yprime, tspan, y0);


%% Plot
plot(t, y(:, 1), '-', t, y(:,3), '--')
legend('dof1', 'dof2')


%% Functions
function f = dydt(y, t, tvec, ugddot)
    % This function generates y' = f(y,t)
    % y is the response vector y = [y1, y2, y3, y4]'
    % t is the time for the solution
    % tvec is the time vector for the ugddot
    % ugddot is the ground motion accelerations
    g = 9.81;

    % Model parameters
    gamma = 0.3;
    omega = 1.0;
    zeta = 0.2;
    mu = 0.2;

    % Find interpolated value of ugddot for time t
    ugddot_int = interp1(tvec, ugddot, t);
    
    a1 = (gamma + 1) / gamma * omega ^ 2;
    a2 = (gamma + 1) / gamma * (2 * zeta * omega);
    a3 = a1 * gamma;
    a4 = a2 * gamma;

    A = [
        0, 1, 0, 0;
        -a1, -a2, a1, a2;
        0, 0, 0, 1;
        a3, a4, -a3, -a4;
        ];
    
    b = [
        0;
        - ugddot_int - (gamma+1)/gamma*mu*g*sign(y(2));
        0;
        - ugddot_int;
        ];

    f = A*y + b;
    
end