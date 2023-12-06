%% Run Multiple Cases
EQname = 'RSN185_IMPVALLH2';
ScaleFactor = [3.49];
% ScaleFactor = 3.15*[0.21 0.74 1.0 1.38];

TrainSpeed  = 0.01; %[0.01 80 240]; % 80 240

Coupled     = [1]; % 1 = Coupled   ; 0 = Uncoupled
BridgeModel = 1; %[0 1]; %[0 1]; % 1 = Nonlinear ; 0 = Linear

% Select Time Steps
dtb = 5.0E-4;
dtt = 5.0E-4;

AnalysisID = 457;
for sf = 1:length(ScaleFactor)
% Define Scale Factor of Analysis
SF = ScaleFactor(sf);

for ts = 1:length(TrainSpeed)
% Define Train Velocity
Vel = TrainSpeed(ts);

for c = 1:length(Coupled)
% Define Coupled/Uncoupled Analysis
CF = Coupled(c);

for bm = 1:length(BridgeModel)
% Define Linear/Nonlinear Analysis
BM = BridgeModel(bm);

% Run Analysis
AnalysisID = AnalysisID + 1;
disp('\\---')
disp(['Running Analysis ',num2str(AnalysisID)])
disp(['ScaleFactor = ',num2str(SF),' \\ Coupled = ',num2str(CF),' \\ Speed = ',num2str(Vel),' \\ BridgeModel = ',num2str(BM)])
disp('---')
RunGroundMotion_Multi

% Store data
FileName = ['SimResults\SimRes_',num2str(AnalysisID,'%04.f'),EQname];

% Create Log File
% LogFile{AnalysisID,1} = FileName;   % Filename
% LogFile{AnalysisID,2} = dtb;        % Time step for bridge
% LogFile{AnalysisID,3} = dtt;        % Time step for train
% LogFile{AnalysisID,4} = SF;         % Scale Factor
% LogFile{AnalysisID,5} = EQname;     % Ground Motion ID
% LogFile{AnalysisID,6} = BM;         % Model Type (Linear Elastic / Nonlinear)
% LogFile{AnalysisID,7} = CF;         % Coupled / Uncoupled
% LogFile{AnalysisID,8} = Vel;        % Train Speed
% 
% % Main response parameters
% LogFile{AnalysisID,9}  = max(abs(BridgeResponse.X(1,:)));   % Max bridge lateral displacement
% LogFile{AnalysisID,10} = max(abs(BridgeResponse.X(2,:)));   % Max bridge hinge rotation
% LogFile{AnalysisID,11} = max(abs(ugddot));  % PGA
% LogFile{AnalysisID,12} = max(abs(ugdot));   % PGV
% LogFile{AnalysisID,13} = max(abs(urec));   % PGD

% Store simulation results
save(FileName)
disp('\\---')
end
end
end
end
