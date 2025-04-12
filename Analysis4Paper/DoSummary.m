%% Create data matrix for analysis
clear, clc, close all

analysis_type = 'OB';  % Here select [OB: over bridge | AG: at grade]
saveName = strcat('SummTable', analysis_type);

sourceFolder = 'C:\Users\Miguel.MIGUEL-DESK\Documents\GitHub\TSI_Matlab\TSI_V5\Results\';
destinationFolder = 'C:\Users\Miguel.MIGUEL-DESK\Documents\GitHub\TSI_Matlab\TSI_V5\Results\Summary';

% Names and information of all runs
TheFiles = dir(strcat(sourceFolder, 'Runs_', analysis_type, '\*.mat'));

% Show how many files were loaded
disp(['Loaded ' num2str(length(TheFiles)) ' files.'])

SummMatrix = zeros(length(TheFiles), 22);
doPlots = false;
mode = 'piola';

% Load each run file and extract relevant information

for File = 1:length(SummMatrix)
    % Load DataSet
    [DataSet] = LoadData(strcat(TheFiles(File).folder,'\' ,TheFiles(File).name));
    
    % Get ground motion ID
    groundmotion = str2double(TheFiles(File).name(6:7));
    
    % Relative wheelset/track displacement
    RelDispl = 100 / 2.54 * (DataSet.TrainResponse.X(7, :) - DataSet.BridgeResponse.X_Track(1, :)); % inches
    
    % Rotation of the track
    TrackRotation = (DataSet.BridgeResponse.X(3,:) + (-DataSet.BridgeResponse.X(8,:) + DataSet.BridgeResponse.X(5,:)) / 1.7742);
    
    % Relative rotation track/wheelset
    RelRotation = rad2deg(abs(-DataSet.TrainResponse.X(9, :) - TrackRotation));

    % Check for derailment (if relative displacement at the end of the analysis is large, then that's derailment).
    if or(max(abs(RelDispl)) >= 5.0, max(abs(RelRotation)) > 45)
        DRCase = 1;
    else
        DRCase = 0;
    end
    
    % Now, check for derailment type and find derailment instant
    climb_strt = find(abs(RelDispl) > 1.0, 1);     % Start climbing instant
    climb_full = find(abs(RelDispl) > 3.0, 1);     % Full climbing instant
    overt_full = find(abs(RelRotation) > 30, 1);   % Overturning instant
    
    if DRCase == 1 % If derailment was identified
        if ~isempty(overt_full) && abs(DataSet.NormalForce(overt_full)) > 0
            if ~strcmp(mode, 'piola')
                disp('Overturning Derail')
            end
            firstDerail = overt_full;
            drType = 3;

        elseif ~isempty(climb_full) && abs(RelRotation(climb_full)) < 3.0
            if ~strcmp(mode, 'piola')
                disp('Sliding Derail')
            end
            firstDerail = climb_full;
            drType = 1;
            
        else
            if ~strcmp(mode, 'piola')
                disp('Combined Derail')
            end
            % Find whether there are 2000 datapoints with derailment (1
            % second of full uplift of both wheels means derailment)
            firstDerail = max(find(DataSet.NormalForce(50:end)==0, 4000)) - 4000;
            drType = 2;

        end
    else
        if ~strcmp(mode, 'piola')
            disp('No derailment')
        end
        firstDerail = length(RelDispl);
        drType = 0;
    end
    
    % Get max response of bridge (for the whole analysis)
    PeakDBridge = max(abs(DataSet.BridgeResponse.X(1,1:end)));
    PeakVBridge = max(abs(DataSet.BridgeResponse.Xtdot(1,1:end)));
    PeakABridge = max(abs(DataSet.BridgeResponse.Xddot(1,1:end) + DataSet.ugddot));
    PeakRBridge = max(abs(DataSet.BridgeResponse.X(2,1:end)));
    
    % Get max response of train (until derailment)
    PeakDCarBody = max(abs(DataSet.TrainResponse.X(1,1:firstDerail)));
    PeakDWheelSet = max(abs(DataSet.TrainResponse.X(7,1:firstDerail)));
    PeakRCarBody = max(abs(DataSet.TrainResponse.X(3,1:firstDerail)));
    PeakRWheelSet = max(abs(DataSet.TrainResponse.X(9,1:firstDerail)));
    
    % Get peak relative rotation between wheelset and track
    PeakRotUplift = max(abs(RelRotation(1:firstDerail)));    % Degrees

    % Below are just for checking
    % disp(PeakRotUplift/0.836)
    % disp(firstDerail)

    % Get peak response of train
    PeakAhCarBody = max(abs(DataSet.TrainResponse.A(1,1:firstDerail)));
    PeakAvCarBody = max(abs(DataSet.TrainResponse.A(2,1:firstDerail)));
    
    % Ground Motion parameters
    PGA = max(abs(DataSet.ugddot));
    PGV = max(abs(DataSet.ugdot));
    
    % Get max uplift value
    peakUplift = max([max(abs(DataSet.Uplift(1:firstDerail))), max(abs(DataSet.Uplift(2:firstDerail)))]);

    % Store Variables
    SummMatrix(File,:) = [...
        DataSet.HL, DataSet.CF, DataSet.Vel, DataSet.BM, PGA, PGV,...
        PeakDBridge, PeakVBridge, PeakABridge, PeakRBridge,...
        PeakDCarBody, PeakDWheelSet, PeakRCarBody, PeakRWheelSet,...
        PeakRotUplift, PeakAhCarBody, PeakAvCarBody,...
        abs(PeakDBridge - PeakDWheelSet), DRCase, groundmotion, drType, peakUplift];
    
    % If wanted, create plots
    if doPlots
        figure()
        plot(TrackRotation(1: firstDerail))
        hold on
        plot(-DataSet.TrainResponse.X(9, 1:firstDerail), ':')
        title('Wheelset rotation')
    
        figure()
        plot(RelRotation(1: firstDerail))
        title('Relative wheel/rail rotation')
    
        figure()
        plot(DataSet.TrainResponse.X(7, 1:firstDerail))
        hold on
        plot(DataSet.BridgeResponse.X_Track(1, 1:firstDerail))
        title('Lateral Displacements of train and track')
    
        figure()
        plot(RelDispl(1:firstDerail))
        title('Relative Displacements')
    end

    clear DataSet
end

SummTable = array2table(SummMatrix,...
    'VariableNames', ...
    {'hazlvl', 'cf', 'spd', 'bridgemodel', 'pga', 'pgv', ...
    'pbd', 'pbv', 'pba', 'pbrot', ... 
    'pcd', 'pwd', 'pcrot', 'pwrot', ...
    'purot', 'pcha', 'pcva', 'drdis', 'drcase', 'gm', 'drtype', 'peakUplift'});

save(saveName, "SummTable")
disp('Done! Summary matrix saved!')

%% Functions

function [Data] = LoadData(FileName)
    % This function loads analysis and extracts required data
    load(FileName, 'BridgeResponse', 'X', 'V', 'A', 'HL', 'Vel', ...
        'ugddot', 'ugdot', 'NL', 'CF', 'Left_Cont', 'Right_Cont', 'uplift');
    
    Data.HL = HL;
    Data.ugddot = ugddot;
    Data.ugdot = ugdot;
    Data.BM = NL;
    Data.BridgeResponse = BridgeResponse;
    Data.Vel = Vel;
    Data.TrainResponse.X = X;
    Data.TrainResponse.V = V;
    Data.TrainResponse.A = A;
    Data.CF = CF;
    Data.NormalForce = Left_Cont + Right_Cont;
    Data.Uplift = uplift;
end




