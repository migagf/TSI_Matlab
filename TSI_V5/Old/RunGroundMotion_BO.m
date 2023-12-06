%% RunGroundMotion No Train
% This code runs a time history analysis over the bridge only

%% Solution of the equation of motion
% This code solves the equation of motion of the system, using an explicit
% scheme for the integration

% nfig = 0;

%% Solution Parameters
% dtb = 5.0e-4; % time step (sec)
g = 9.81; % m/s2

%% Forcing parameters
% SF = 3.15; % Define Scale Factor
% EQname = 'RSN185_IMPVALLH2'; % Define Ground Motion
% AnalysisID = 1;

load([cd,'\GMs\Used Records\',EQname])
ugddot = SF*TimeAccelData(:,2)*9.81;        % Acceleration Time-History

dtrec = round(TimeAccelData(2,1),3);        % Time step of record
trec  = 0:dtrec:dtrec*(length(ugddot)-1);   % Time vector 

ugdot  = dtrec*cumtrapz(ugddot);          % Velocity Time-History
urec   = dtrec*cumtrapz(ugddot);          % Displacement Time-History

% Resample ground motion
tb = 0:dtb:trec(end);  % Time vector for bridge analysis

ugddot  = interp1(trec,ugddot,tb);  % Excitation is required only for bridge
ugdot   = interp1(trec,ugdot,tb);   % Excitation is required only for bridge
urec    = interp1(trec,urec,tb);    % Excitation is required only for bridge

bsteps = length(tb);

% nfig = nfig + 1;
% figure(nfig)
% subplot(3,1,1), plot(tb,ugddot/9.81), xlabel('Time (sec)'), ylabel('Accel. (g)'), title('Ground Motion (LPE)')
% subplot(3,1,2), plot(tb,ugdot), xlabel('Time (sec)'), ylabel('Vel. (m/sec)'), title('Ground Motion (LPE)')
% subplot(3,1,3), plot(tb,urec), xlabel('Time (sec)'), ylabel('Disp. (m)'), title('Ground Motion (LPE)')

%% Initialize Bridge Data
% BM = 0;
if BM == 1
    InitializeBridgeModel_DC
else
    InitializeBridgeModel_DC_LE
end

%% Initial conditions and prellocation of variables

BridgeResponse.X_Track = zeros(2,bsteps);
BridgeResponse.V_Track = zeros(2,bsteps);
CF = 0; Cont_Force = zeros(1,9);

%% Solution of the EOM

% Integration parameters
psi = 0.5;
phi = 0.5;

tic

for ib = 2:bsteps
    BridgeTimeIncrement;
    UpdateX_Track;
end

toc

%% Save Simulation Results
% FileName = ['SimResults_BO\SimRes_',num2str(AnalysisID,'%03.f'),EQname];
% save(FileName)

%%  Post Processing

% PostProcessingPM
% PostProcessingGM
% PostProcessingGMNoTrain
