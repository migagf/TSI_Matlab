%% Run Multiple Decoupled Train Analyses
% First, run Bridge-Only analyses. Then, run the decoupled analysis. 

TheFiles1 = dir('SimResults_BO/*.mat');

TrainSpeed  = [0.01 80 240]; % 80 240

% BridgeModel = [0 1]; %[0 1]; % 1 = Nonlinear ; 0 = Linear

% Select Time Step for train
dtt = 5.0E-4;

AnalysisIDc = 0;

for file = 1:length(TheFiles1)
% Define Scale Factor of Analysis
load(['SimResults_BO\',TheFiles1(file).name]);

    for ts = 1:length(TrainSpeed)
        % Define Linear/Nonlinear Analysis
        Vel = TrainSpeed(ts);
        
        % Run Analysis
        AnalysisIDc = AnalysisIDc + 1;
        
        % Run Bridge Only
        disp('\\---')
        disp(['Running Decoupled Analysis #',num2str(AnalysisIDc)])
        disp(['ScaleFactor = ',num2str(SF),' \\ BridgeModel = ',num2str(BM),...
              ' \\ Train Speed = ',num2str(Vel)])
        
        
        
        % Store data
        FileName = ['SimResults_DC\SimRes_DC_',num2str(AnalysisIDc,'%04.f'),EQname];
        disp(TheFiles1(file).name)
        try
            RunGroundMotion_DC
        catch
            disp(['Skipped Analysis',num2str(AnalysisIDc)])
        end
        % Store simulation results
        save(FileName)
        disp('\\---')
    
    end
end

disp('Complete')