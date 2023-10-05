%% Run Multiple Bridge Only Analyses
TheFiles = dir('GMs\Used Records\*.mat');
GMScaleFactors = ...
    2.0*[3.49, 3.12, 3.11, 4.25, 3.28, 2.60, 2.95, 2.26, 2.06, 3.05,...
     3.08, 3.88, 2.71, 2.32, 1.81, 2.19, 1.24, 2.32, 2.65, 3.15];

AnalysisID = 0;
TrainSpeed  = [0.01 80 240];

for file = 1:length(TheFiles)

EQname = TheFiles(file).name(1:end-4);
% EQname = 'RSN185_IMPVALLH2';
% ScaleFactor = [3.49];

ScaleFactor = GMScaleFactors(file)*[0.21 0.74 1.0 1.38];

% ; % 80 240

BridgeModel = 1; % 1 = Nonlinear ; 0 = Linear

    for sf = 1:length(ScaleFactor)
        % Define Scale Factor of Analysis
        SF = ScaleFactor(sf);
        for Vel = TrainSpeed
        % Run Analysis
        AnalysisID = AnalysisID + 1;
        tic
        % Run Bridge Only
        disp('\\---')
        disp(['Running At Grade Analysis #',num2str(AnalysisID),' // Speed = ',num2str(Vel)])
        disp(EQname)
        
        RunGroundMotion_AG
        
        % Store data
        FileName = ['SimResults_AG\SimRes_AG_',num2str(AnalysisID,'%04.f'),EQname];
        
        % Store simulation results
        save(FileName)
        disp('\\---')
        toc
        end
    end
end

disp('Complete')
