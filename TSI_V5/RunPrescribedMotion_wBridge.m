%% RunPrescribedMotion_wBridge
% Runs wavelets on bridge-train system
clear, clc, close

%% Load data
load WheelGeom
load RailGeom
load RailProps

RailProps.R1y = 0.2;

RailProps.E1 = RailProps.E1/10;
RailProps.E2 = RailProps.E2/10;

%% Load Interface Data
CreateWheelRailGeom

%% Solution Parameters
dt = 1.0e-3;
dtt  = dt;   % time step (sec) of train
dtb  = dt;   % time step (sec) of bridge 
g = 9.81;
Vel = 80/3.6;   % m/sec
CF  = 1.0;        % 1 for coupled analysis, 0 for uncoupled

% Scale Factors
SFvals = 0.05:0.05:5;

% Frequency Values
freqvals = 0.02:0.02:1;

initID = 1;

ID = initID;

for freq = freqvals   % Factor that modifies the frequency of the pulse
for SF = SFvals

disp(['Running Analysis #',num2str(ID),' of ',num2str(length(SFvals)*length(freqvals))])

%% Loading parameters (wavelet definition)

% Create a wavelet for analysis of derailment
wp = 2*pi*freq;
[time, accel] = RickerWave(wp);
t = 0:dt:max(time);
ugddot = interp1(time,accel,t) * SF * g;

tt = t;
tb = t;
tsteps = length(tt);
bsteps = length(tb);

%% Initial conditions and prellocation of variables

X = zeros(9,tsteps); % Train Global Coordinates 
V = zeros(9,tsteps); % Train Velocities
A = zeros(9,tsteps); % Train Accelerations

load X_initial

X0 = X_initial;             % Initial Global Coordinates
V0 = [0 0 0 0 0 0 0 0 0]';  % Initial Velocities
A0 = [0 0 0 0 0 0 0 0 0]';  % Initial Accelerations

X1 = X_initial;             % Initial Global Coordinates
V1 = [0 0 0 0 0 0 0 0 0]';  % Initial Velocities
A1 = [0 0 0 0 0 0 0 0 0]';  % Initial Accelerations

X(:,1) = X0; X(:,2) = X1;
V(:,1) = V0; V(:,2) = V1;
A(:,1) = A0; A(:,2) = A1;

BridgeResponse.X_Track = zeros(2,bsteps);
BridgeResponse.V_Track = zeros(2,bsteps);

% Extra output variables
Momt = [0 0 0 0]';
Creepforces = zeros(4,tsteps);
delta       = zeros(4,tsteps);
deltadotn_vec = [0 0 0 0];

%% Mass, Stiffness and Damping Matrices
% Call the functions that create the stiffness and damping matrices

[K] = StiffnessMatrix(); 
[C] = DampingMatrix();
[M,Mc,Mt,Mw] = MassMatrix();

% AddRayleighDamp = false;
% if AddRayleighDamp
%     w1 = 2*pi*100; %rad/sec
%     a1 = 2/w1;
%     C = C + a1*K;
% end

% EigenValue Analysis
% Options for EigenValue Analysis
% showplot = false;
% EigenValueAnalysis;

%% Initialize Bridge Data

InitializeBridgeModel

%% Gravity Loading
F_ine = [0 Mc*g 0 0 Mt*g 0 0 Mw*g 0]'; % Inertial Forces (N)

%% Solution of the EOM
tic

ib = 2; tbridge = 0;
psi = 0.5;
phi = 0.5;

for it = 2:tsteps
    % Integration of Train EOM
    X2 = X1+V1*dtt+(0.5+psi)*A1*dtt^2-psi*A0*dtt^2;
    V2 = V1+(1+phi)*A1*dtt-phi*A0*dtt;
   
    % Check if derailment has been reached
    if (abs(X2(7)-BridgeResponse.X_Track(1,it-1)) > 0.1 && abs(X2(9)) < pi/30) || ...
            abs(X2(9)) > pi/2 || ...
            (abs(X2(7)-BridgeResponse.X_Track(1,it-1)) > 1.0 && abs(X2(9)) < pi/2)
        break
    else
        % Update the force vector with Contact Algorithm
        [F,NF_L,NF_R,vec,Ft,delta(:,it-1),Momt] = ContactForce(X2(7:9)', ...
                                                  V2(7:9)', ...
                                                  [BridgeResponse.X_Track(1,ib-1), BridgeResponse.X_Track(2,ib-1),0], ...
                                                  [BridgeResponse.V_Track(1,ib-1), BridgeResponse.V_Track(2,ib-1),0], ...
                                                  WheelGeom_pol, ...
                                                  RailGeom_pol,...
                                                  RailGeom_cur, ...
                                                  WheelGeom_cur, ...
                                                  RailProps, ...
                                                  0, ...
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
        
        tbridge = tbridge + dtt;
        
        % Integrate Bridge EOM
        if tbridge >= dtb
            BridgeTimeIncrement
            UpdateX_Track;
            ib = ib + 1;
            tbridge = 0;
        end
        
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
end

toc

% Save Results
save(['WaveletWBridge\Wave_ID_',num2str(ID,'%04.f'),'_wF_',num2str(freq*100),'_sF_',num2str(SF*100),'.mat'])
ID = ID + 1;

end
end

%% Post Processing Block
% Nothing here (yet)

%% Functions
% Ricker wavelet
function [tau, r] = RickerWave(wp)
    % Function creates a ricker wavelet with the specified frequency
    tau = linspace(-10/wp,10/wp,10000);
    r = (1-0.5*wp^2*tau.^2).*exp(-1/4*wp^2*tau.^2);
    tau = tau - min(tau);
    % plot(tau,r)
end