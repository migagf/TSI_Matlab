%% Solution of the equation of motion
% This code solves the equation of motion of the system, using an explicit
% scheme for the integration
clear, clc, close all
nfig = 0;
showplot = 0;
on_bridge = false;

% Factor for LE Bridge or NL Bridge
NL = 1.0;

SFvals   = 0.1:0.05:5.0;
freqvals = 0.1:0.05:5.0;

for SF = SFvals
for freq = freqvals

%% Load Data
load WheelGeom
load RailGeom
load RailProps
RailProps.E1 = RailProps.E1;
RailProps.E1 = RailProps.E1;


%% Solution Parameters
dtt = 5.0E-4;   % time step (sec) of train
dtb = 5.0E-4;   % time step (sec) of bridge 

% Train Finite Difference Sol. Parameters
psi = 0.5;
phi = 0.5;

g   = 9.81;       % m/s2
Vel = 0.01/3.6;   % m/sec
CF  = 1.0;        % 1 for coupled analysis, 0 for uncoupled

%% Forcing parameters
%SF = 1.0;
%freq = 1.0;

run_filename = strcat('Pulse_AG_SF', num2str(100 * SF, "%03.0f"), '_FQ', num2str(100 * freq, "%03.0f"));

disp(strcat("Running Case:", run_filename))

% Create a wavelet for analysis of derailment
wp = 2*pi*freq;

[trec, ugddot] = RickerWave(wp);    % Create Ricker Wavelet
ugddot = ugddot' * SF * g;           % Apply scaling

dtrec = trec(2);                             % Time step of record.
ugddot = [0*(0:dtrec:4)'; ugddot];           % Add dead time before applying ground motion
trec  = 0:dtrec:dtrec*(length(ugddot)-1);    % Time vector 

%ugdot  = [0 diff(urec')]/(dtrec);           % Velocity Time-History
%ugddot = [0 diff(ugdot)]/(dtrec);           % Acceleration Time-History

ugdot = cumtrapz(ugddot) * dtrec;
urec  = cumtrapz(ugdot) * dtrec;

% Resampling of Acceleration Time History
tt = 0:dtt:trec(end);  % Time vector for train analysis
tb = 0:dtb:trec(end);  % Time vector for bridge analysis

ugddot  = interp1(trec,ugddot,tb);  % Excitation is required only for bridge
ugdot   = interp1(trec,ugdot,tb);   % Excitation is required only for bridge
urec    = interp1(trec,urec,tb);    % Excitation is required only for bridge

ug = urec;

tsteps = length(tt);
bsteps = length(tb);

% % Plot EQ record
% nfig = nfig + 1;
% 
% figure() 
% subplot(3,1,1), plot(tt,ugddot/9.81), xlabel('Time (sec)'), ylabel('Accel. (g)'), title('Ground Motion (LPE)') 
% subplot(3,1,2), plot(tt,ugdot), xlabel('Time (sec)'), ylabel('Vel. (m/s)'), title('Ground Motion (LPE)')
% subplot(3,1,3), plot(tt,urec), xlabel('Time (sec)'), ylabel('Dis. (m)'), title('Ground Motion (LPE)')

%% Load Train Data
CreateWheelRailGeom

%% Initialize Bridge Data
InitializeBridgeModel

%% Initial conditions and prellocation of variables

X = zeros(9,tsteps); % Train Global Coordinates 
V = zeros(9,tsteps); % Train Velocities
A = zeros(9,tsteps); % Train Accelerations

%load X_initial

%X0 = X_initial;      % Initial Global Coordinates
% X0(1) = 0; X0(4) = 0; X0(7) = 0.0;

%X0 = 0.5655 * [0 -1 0 0 -1 0 0 -1 0]';

if ~on_bridge
    X0 = [-0.0119, -0.3943, -0.0022, -0.0082, -0.4216, -0.0022, -0.0079, -0.5573, -0.0022]' + ...
        0.01 * [0 -1 0 0 -1 0 0 -1 0]';
else
    X0 = [-0.0119, -0.3943, -0.0022, -0.0082, -0.4216, -0.0022, -0.0079, -0.5573, -0.0022]' + ...
        0.01 * [0 -1 0 0 -1 0 0 -1 0]';
end

V0 = [0 0 0 0 0 0 0 0 0]';  % Initial Velocities
A0 = [0 0 0 0 0 0 0 0 0]';  % Initial Accelerations

X1 = X0;             % Initial Global Coordinates
V1 = [0 0 0 0 0 0 0 0 0]';  % Initial Velocities
A1 = [0 0 0 0 0 0 0 0 0]';  % Initial Accelerations

X(:,1) = X0; X(:,2) = X1;
V(:,1) = V0; V(:,2) = V1;
A(:,1) = A0; A(:,2) = A1;

BridgeResponse.X_Track = zeros(2,bsteps);
BridgeResponse.V_Track = zeros(2,bsteps);

% Extra output variables
Momt = [0 0 0 0]';
Creepforces = zeros(4, tsteps);
delta       = zeros(4, tsteps);
deltadotn_vec = [0 0 0 0];
Left_Cont = zeros(1, tsteps);
Right_Cont =zeros(1, tsteps);


% Contact forces (storage)
QL = zeros(1, tsteps);
QR = zeros(1, tsteps);
YL = zeros(1, tsteps);
YR = zeros(1, tsteps);


%% Mass, Stiffness and Damping Matrices
% Call the functions that create the stiffness and damping matrices

[K] = StiffnessMatrix(); 
[C] = DampingMatrix();
[M,Mc,Mt,Mw] = MassMatrix();

AddRayleighDamp = false;
if AddRayleighDamp
    w1 = 2 * pi * 100; %rad/sec
    a1 = 2/w1;
    C = C + a1*K;
end

% % EigenValue Analysis
% % Options for EigenValue Analysis
% showplot = false;
% EigenValueAnalysis;

%% Gravity Loading
F_ine = [0 Mc*g 0 0 Mt*g 0 0 Mw*g 0]'; % Inertial Forces (N)

%% Solution of the EOM
tic

ib = 2; tbridge = 0;

for it = 2:tsteps
    % Integration of Train EOM
    X2 = X1 + V1 * dtt + (0.5 + psi) * A1 * dtt ^ 2 - psi * A0 * dtt ^ 2;
    V2 = V1 + (1 + phi) * A1 * dtt - phi * A0 * dtt;

    % Check for derailment
    if (abs(X2(7) - BridgeResponse.X_Track(1, it - 1)) > 0.1 && abs(X2(9)) < pi/30)
        dr_type = 1;   % Derailment in slide-off mode
        disp("Slide-off Derailment")
        break
    elseif abs(X2(9)) > pi/2
        dr_type = 2;   % Derailment in overturning mode
        disp("Overturning Derailment")
        break
    elseif (abs(X2(7) - BridgeResponse.X_Track(1, it - 1)) > 1.0 && abs(X2(9)) < pi/2)
        dr_type = 3;
        disp("Combined Derailment")
        break           % Derailment in combined mode
    else 
        dr_type = 0;
    end

    
    % Update the force vector with Contact Algorithm
    [F,NF_L,NF_R,vec,Ft,delta(:,it-1),Momt] = ContactForce(...
        X2(7:9)', ...
        V2(7:9)', ...
        [BridgeResponse.X_Track(1, ib-1), BridgeResponse.X_Track(2, ib-1), -BridgeResponse.X(3, ib-1)], ...
        BridgeResponse.Xt(end-5:end,ib-1), ...
        BridgeResponse.Xtdot(end-5:end,ib-1), ...
        WheelGeom_pol, ...
        RailGeom_pol, ...
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
        RWheel_geom_fine ...
        );
    
    % Force Vector to Structure
    Cont_Force = 10 ^ 6 * [0, 0, 0, -F(1, 1), -F(2, 1), -F(1, 1) * 0.1, -F(1, 2), -F(2, 2), -F(1, 2) * 0.1]';

    % Force vector to train
    Cforce = 10 ^ 6 * [0, 0, 0, 0, 0, 0, F(1, 1) + F(1, 2), F(2, 1) + F(2, 2), F(3, 1) + F(3, 2)]';
    
    F_ext =  F_ine + Cforce; % N
    
    tbridge = tbridge + dtt;

    if tbridge >= dtb
        BridgeTimeIncrement
        UpdateX_Track;
        ib = ib + 1;
        tbridge = 0;
    end
    
    % Variable storage
    % vrel(i) = V2(7)-vx_track(i); vect(i,:) = Ft;
    Left_Cont(it) = 1000 * NF_L; 
    Right_Cont(it) = 1000 * NF_R;

    % Contact forces (storage)
    QL(it) = F(2, 1);
    QR(it) = F(2, 2);
    YL(it) = F(1, 1);
    YR(it) = F(1, 2);

    %Numerical integration of Train EOM
    A2 = M \ (F_ext - K * X2 - C * V2);
    
    X(:,it) = X2; 
    V(:,it) = V2; 
    A(:,it) = A2;
    
    % Update the state vectors
    X0 = X1; V0 = V1; A0 = A1;
    X1 = X2; V1 = V2; A1 = A2;

end

toc

save(strcat("Pulses_NoBridge\", run_filename))
disp("Done...")

% AnimateEarthquake
clearvars -except CF dtb dtt freqvals g nfig NL on_bridge phi psi RailGeom RailProps SFvals showplot Vel WheelGeom SF freq

end
end
%%  Post Processing

PostProcessingGM

%% Functions
% Ricker wavelet
function [tau, r] = RickerWave(wp)
    % Function creates a ricker wavelet with the specified frequency
    tau = linspace(-10/wp,10/wp,10000);
    r = (1-0.5*wp^2*tau.^2).*exp(-1/4*wp^2*tau.^2);
    tau = tau - min(tau);
    %plot(tau,r)
end
