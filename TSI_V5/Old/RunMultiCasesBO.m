%% Run Multiple Bridge Only Analyses
TheFiles = dir('GMs\Used Records\*.mat');
GMScaleFactors = ...
    [3.49, 3.12, 3.11, 4.25, 3.28, 2.60, 2.95, 2.26, 2.06, 3.05,...
     3.08, 3.88, 2.71, 2.32, 1.81, 2.19, 1.24, 2.32, 2.65, 3.15];

AnalysisID = 20;

for file = 6:length(TheFiles)

EQname = TheFiles(file).name(1:end-4);
% EQname = 'RSN185_IMPVALLH2';
% ScaleFactor = [3.49];

ScaleFactor = GMScaleFactors(file)*[0.21 0.74 1.0 1.38];

% TrainSpeed  = [0.01 80 240]; % 80 240

BridgeModel = 1; % 1 = Nonlinear ; 0 = Linear

% Select Time Steps
dtb = 5.0E-4;

for sf = 1:length(ScaleFactor)
% Define Scale Factor of Analysis
SF = ScaleFactor(sf);

for bm = 1:length(BridgeModel)
% Define Linear/Nonlinear Analysis
BM = BridgeModel(bm);

% Run Analysis
AnalysisID = AnalysisID + 1;

% Run Bridge Only
disp('\\---')
disp(['Running Bridge Only Analysis #',num2str(AnalysisID)])
disp(EQname)

RunGroundMotion_BO

% Store data
FileName = ['SimResults_BO\SimRes_BO_',num2str(AnalysisID,'%04.f'),EQname];

% Create Log File
% LogFile{AnalysisID,1} = FileName;   % Filename
% LogFile{AnalysisID,2} = dtb;        % Time step for bridge
% LogFile{AnalysisID,3} = 'N/A';        % Time step for train
% LogFile{AnalysisID,4} = SF;         % Scale Factor
% LogFile{AnalysisID,5} = EQname;     % Ground Motion ID
% LogFile{AnalysisID,6} = BM;         % Model Type (Linear Elastic / Nonlinear)

% Main response parameters
% LogFile{AnalysisID,7}  = max(abs(BridgeResponse.X(1,:)));   % Max bridge lateral displacement
% LogFile{AnalysisID,8} = max(abs(BridgeResponse.X(2,:)));   % Max bridge hinge rotation
% LogFile{AnalysisID,9} = max(abs(ugddot));  % PGA
% LogFile{AnalysisID,10} = max(abs(ugdot));   % PGV
% LogFile{AnalysisID,11} = max(abs(urec));    % PGD

% Store simulation results
save(FileName)
disp('\\---')

end
end
end
disp('Complete')
