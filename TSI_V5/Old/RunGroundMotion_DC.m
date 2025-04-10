%% Run Multi Cases of Bridge Only Simulation
% clear, clc, close
nfig = 0;
showplot = 0;

%% Load Train Data
load WheelGeom
load RailGeom
load RailProps

%% Solution Parameters
dtt  = 5.0E-4;   % time step (sec) of train

% Train Finite Difference Sol. Parameters
psi = 0.5;
phi = 0.5;

g   = 9.81;       % m/s2
% Vel = 0.01/3.6;   % m/sec MODIFY!!!!!
CF  = 0.0;        % 1 for coupled analysis, 0 for uncoupled

%% Load Train Data
CreateWheelRailGeom

%% Mass, Stiffness and Damping Matrices
% Call the functions that create the stiffness and damping matrices

[K] = StiffnessMatrix(); 
[C] = DampingMatrix();
[M,Mc,Mt,Mw] = MassMatrix();

AddRayleighDamp = false;
if AddRayleighDamp
    w1 = 2*pi*100; %rad/sec
    a1 = 2/w1;
    C = C + a1*K;
end

%% Forcing parameters
% SF = 1.0; % PUT IN THE CALLING FUNCTION

% Resampling of Acceleration Time History
tt = 0:dtt:tb(end);  % Time vector for train analysis

X_Track   = BridgeResponse.X_Track;
V_Track   = BridgeResponse.V_Track;
Phi_Track = -BridgeResponse.X(3,:);
X_Track(:,1) = X_Track(:,3);
X_Track(:,2) = X_Track(:,3);

tsteps = length(tt);

%% Initial conditions and prellocation of variables

X = zeros(9,tsteps); % Train Global Coordinates 
V = zeros(9,tsteps); % Train Velocities
A = zeros(9,tsteps); % Train Accelerations

load X_initial_DC

X0 = X_initial;             % Initial Global Coordinates
V0 = [0 0 0 0 0 0 0 0 0]';  % Initial Velocities
A0 = [0 0 0 0 0 0 0 0 0]';  % Initial Accelerations

X1 = X_initial;             % Initial Global Coordinates
V1 = [0 0 0 0 0 0 0 0 0]';  % Initial Velocities
A1 = [0 0 0 0 0 0 0 0 0]';  % Initial Accelerations

X(:,1) = X0; X(:,2) = X1;
V(:,1) = V0; V(:,2) = V1;
A(:,1) = A0; A(:,2) = A1;

% Extra output variables
Momt = [0 0 0 0]';
Creepforces = zeros(4,tsteps);
delta       = zeros(4,tsteps);
deltadotn_vec = [0 0 0 0];

%% Gravity Loading
F_ine = [0 Mc*g 0 0 Mt*g 0 0 Mw*g 0]'; % Inertial Forces (N)

%% Solution of the EOM
tic

disp(['ScaleFactor = ',num2str(SF),' \\ BridgeModel = ',num2str(BM)])

for it = 2:tsteps
    % Integration of Train EOM
    X2 = X1+V1*dtt+(0.5+psi)*A1*dtt^2-psi*A0*dtt^2;
    V2 = V1+(1+phi)*A1*dtt-phi*A0*dtt;
    
    % Update the force vector with Contact Algorithm
    [F,NF_L,NF_R,vec,Ft,delta(:,it-1),Momt] = ... 
        ContactForce(X2(7:9)', ...
                     V2(7:9)', ...
                     [X_Track(1,it-1), X_Track(2,it-1),Phi_Track(it-1)], ...
                     [V_Track(1,it-1), V_Track(2,it-1),Phi_Track(it-1)], ...
                     WheelGeom_pol, ...
                     RailGeom_pol,...
                     RailGeom_cur, ...
                     WheelGeom_cur, ...
                     RailProps, ...
                     showplot, ...
                     Vel, ...
                     delta(:,it-1), ...
                     Momt, ...
                     deltadotn_vec, ...
                     Rail_geom_fine, ...
                     LWheel_geom_fine, ...
                     RWheel_geom_fine);
    
    % Updated force vector
    Cont_Force = 10^6*[0 0 0 0 0 0 F']';
    F_ext =  F_ine + Cont_Force; % N
    
    % Variable storage
    % vrel(i) = V2(7)-vx_track(i); vect(i,:) = Ft;
    Left_Cont(it) = 1000*NF_L; Right_Cont(it) = 1000*NF_R;

    %Numerical integration of Train EOM
    A2 = M\(F_ext-K*X2-C*V2);
    X(:,it) = X2; V(:,it) = V2; A(:,it) = A2;
    
    % Update the state vectors
    X0 = X1; V0 = V1; A0 = A1;
    X1 = X2; V1 = V2; A1 = A2;
end

toc

%%  Post Processing

% PostProcessingGM

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
