function [u,v,a] = CA_script(m,k,z,dt,ag)
% Constant Acceleration Method
% This script calculates the response of a 1DOF with properties m (kips/g),
% k (kips/in) and damping z, under the action of accelerations from an
% earthquake ag.

% Input Parameters
%k = 12; % stifness (kips/in)
%Tn = 1; % natural period (sec)

% Initial conditions
u0 = 0; % initial displacement
v0 = 0; % initial velocity

% Prellocate the u and v vectors
u = zeros(length(ag),1);
v = zeros(length(ag),1);
a = zeros(length(ag),1);
p_hat = zeros(length(ag),1);
u(1) = u0; v(1) = v0;

%% Previous calculations
wn = sqrt(k/m); % natural freq.
T = 2*pi/wn;

c = 2*z*m*wn;

% Force function
p = -m*ag; % dynamic force

gamma = 1/2;
beta = 1/4;

a0 = (p(1) - k*u0 - c*v0)/m; a(1) = a0;
a1 = 1/(beta*dt^2)*m + gamma/(beta*dt)*c;
a2 = 1/(beta*dt)*m + (gamma/beta - 1)*c;
a3 = (1/(2*beta) - 1)*m + dt*(gamma/(2*beta) - 1)*c;
k_hat = k + a1;

%% Iterations

for i = 1:length(ag)-1
    p_hat(i+1) = p(i+1) + a1*u(i) + a2*v(i) + a3*a(i);
    u(i+1) = p_hat(i+1)/k_hat;
    v(i+1) = gamma/(beta*dt)*(u(i+1)-u(i)) + (1 - gamma/beta)*v(i) + dt*(1 - gamma/(2*beta))*a(i);
    a(i+1) = 1/(beta*dt^2)*(u(i+1) - u(i)) - 1/(beta*dt)*v(i) - (1/(2*beta)-1)*a(i);
end

end