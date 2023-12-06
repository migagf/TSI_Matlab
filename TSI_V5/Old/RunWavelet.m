%% Solution of the equation of motion
% This code solves the equation of motion of the system, using an explicit
% scheme for the integration
% clear, clc, close
showplot = 0;

%% Load data
load WheelGeom
load RailGeom
load RailProps

RailProps.R1y = 0.2;

%% Solution Parameters
dt = 1.0E-4; % time step (sec)
psi = 0.5;
phi = 0.5;
g = 9.81; % m/s2
Vel = 1000/3.6; % m/sec

% Initial conditions
load X_initial

X0 = X_initial; %[0 Z1c 0 0 Z1b 0 0 Z1w 0]';  % Initial positions
V0 = [0 0 0 0 0 0 0 0 0]';  % Initial velocities
A0 = [0 0 0 0 0 0 0 0 0]';  % Initial accelerations

X1 = X_initial; %[0 Z1c 0 0 Z1b 0 0 Z1w 0]';  % 2nd initial position
V1 = [0 0 0 0 0 0 0 0 0]';  % 2nd initial velocities
A1 = [0 0 0 0 0 0 0 0 0]';  % Initial accelerations

%% Mass, Stiffness and Damping Matrices
% Call the functions that create the mass, stiffness and damping matrices

[K] = StiffnessMatrix(); 
[C] = DampingMatrix();

[M,Mc,Mt,Mw] = MassMatrix();
% % 
% w1 = 2*pi*100; %rad/sec
% a1 = 2/w1;
% C = C + a1*K;

%% Iterations

F_ine = [0*Mc*g Mc*g 0 0 Mt*g 0 0 Mw*g 0]'; % N
%load Cont_Force %Load initial contact force
F_ext = F_ine; %+ Cont_Force;

WheelGeom = 25.4/1000*WheelGeom;
% Interpolation polynomials of the wheel and rail geometries
RailGeom_pol = pchip(RailGeom(:,1),RailGeom(:,2));
WheelGeom_pol = makima(WheelGeom(:,1),WheelGeom(:,2));

% Forcing parameters

fm = 4; w = 2*pi*fm; T = 1/fm;
Ag = 0;
a = Ag/(w^2); % Amplitude

%%
total = 3*T; %10*dt;
t = 0:dt:total;

x_track = a*(1-2*pi^2*fm^2*(t).^2).*exp(-pi^2*fm^2.*(t).^2);
x_track = [x_track(end:-1:1),x_track];
v_track = [0 (1/dt)*(x_track(2:end)-x_track(1:end-1))];

total = (length(x_track)-1)*dt;

%plot(t,x_track(1:length(t)))
%%
tic
vari = 1;
steps = round(total/dt);

X = zeros(9,steps); V = X; A = X;

X(:,1) = X0; X(:,2) = X1;
V(:,1) = V0; V(:,2) = V1;
A(:,1) = A0; A(:,2) = A1;
Momt = [0 0 0 0]';
Creepforces = zeros(4,steps);
delta = zeros(4,length(x_track));
deltadotn_vec = [0 0 0 0];

for i = 2:steps

    X2 = X1+V1*dt+(0.5+psi)*A1*dt^2-psi*A0*dt^2;
    V2 = V1+(1+phi)*A1*dt-phi*A0*dt;
    
    [F,NF_L,NF_R,vec,Ft,delta(:,i),Momt] = ...
        ContactForce(X2(7:9)',V2(7:9)',[x_track(i),0,0],...
        [v_track(i),0,0],WheelGeom_pol,RailGeom_pol,...
        RailProps,showplot,Vel,delta(:,i-1),Momt,deltadotn_vec);
    
    Cont_Force = 10^6*[0 0 0 0 0 0 F']';
    F_ext =  F_ine + Cont_Force; % N
    
    vrel(i) = V2(7)-v_track(i);
    vect(i,:) = Ft;
    A2 = M\(F_ext-K*X2-C*V2);
    Creepforces(:,i) = vec;

    Left_Cont(i) = 1000*NF_L; Right_Cont(i) = 1000*NF_R;
    
    X(:,i) = X2; V(:,i) = V2; A(:,i) = A2;
    
    X0 = X1; V0 = V1; A0 = A1;
    X1 = X2; V1 = V2; A1 = A2;

end

toc
%%  Post Processing
% plotting_animation
t = linspace(0,total,steps);

Ek = sum(M*V.^2)/2;
Ep = sum(K*X.^2)/2;

figure(1)
plot(Left_Cont)
hold on
plot(Right_Cont)
% plot(t(1:end-1),Ek+Ep,t(1:end-1),Ek,t(1:end-1),Ep)
grid on
xlabel('Time (sec)'), ylabel('Contact Force (kN)')

post_processing
% hold off
% figure(3)
% plot(t,vrel) %Relative velocity between wheelset and rail
% hold on
% plot(t,vect(:,1)) % Tangential force in the Y direction
% plot(t(1:end),Left_Cont/1000) % Left contact force
% plot(t(1:end),Right_Cont/1000) % Right contact force
% plot(t(1:end),A(7,:)/981)
% grid on
% legend('Relative velocity','Tangential force','Left normal force','Right normal force')
