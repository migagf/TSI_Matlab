% Load each run and store global response parameters
TheFiles = dir("C:\Users\Miguel Gomez\Documents\PhD Files\TSI_Runs\Pulses_AG\*.mat");

for File = 1:length(TheFiles)
    % Load DataSet
    [DataSet] = LoadData(strcat(TheFiles(File).folder,'\' ,TheFiles(File).name));
    groundmotion = str2double(TheFiles(File).name(6:7));

    DRCase = DataSet.dr_type;

    % Obtain max response of bridge
    PeakDBridge = max(abs(DataSet.BridgeResponse.X(1,1:end)));
    PeakVBridge = max(abs(DataSet.BridgeResponse.Xdot(1,1:end)));
    PeakABridge = max(abs(DataSet.BridgeResponse.Xddot(1,1:end)));
    PeakRBridge = max(abs(DataSet.BridgeResponse.X(2,1:end)));
    

    % Obtain max response of train
    PeakDCarBody   = max(abs(DataSet.TrainResponse.X(1,1:firstDerail)));
    PeakDWheelSet  = max(abs(DataSet.TrainResponse.X(7,1:firstDerail)));
    PeakRCarBody   = max(abs(DataSet.TrainResponse.X(3,1:firstDerail)));
    PeakRWheelSet  = max(abs(DataSet.TrainResponse.X(9,1:firstDerail)));
    

    PeakRotUplift  = max(abs(DataSet.TrainResponse.X(9, 1:firstDerail) + DataSet.BridgeResponse.X(3, 1:firstDerail) + ...
        - (DataSet.BridgeResponse.X(8, 1:firstDerail) - DataSet.BridgeResponse.X(5, 1:firstDerail))/1.2));
    

    % plot(DataSet.TrainResponse.X(9,:)+DataSet.BridgeResponse.X(3,:))
    %PeakLftWUplift = max(abs(DataSet.TrainResponse.Uplift(:,1)));
    %PeakRgtWUplift = max(abs(DataSet.TrainResponse.Uplift(:,2)));
    PeakAhCarBody  = max(abs(DataSet.TrainResponse.A(1,1:firstDerail)));
    PeakAvCarBody  = max(abs(DataSet.TrainResponse.A(2,1:firstDerail)));
    
    % Ground Motion parameters
    PGA = max(abs(DataSet.ugddot));
    PGV = max(abs(DataSet.ugdot));
    
    sf = DataSet.sf;
    freq = DataSet.freq;

    % Store Variables
    SummMatrix(File,:) = [...
        DataSet.HL, DataSet.CF, DataSet.Vel, DataSet.BM, PGA, PGV,...
        PeakDBridge, PeakVBridge, PeakABridge, PeakRBridge,...
        PeakDCarBody, PeakDWheelSet, PeakRCarBody, PeakRWheelSet,...
        PeakRotUplift, PeakAhCarBody, PeakAvCarBody,...
        abs(PeakDBridge - PeakDWheelSet), DRCase, groundmotion, sf, freq];
    
   
    clear DataSet

end

SummTable = array2table(SummMatrix,...
    'VariableNames', ...
    {'hazlvl', 'cf', 'spd', 'bridgemodel', 'pga', 'pgv', ...
    'pbd', 'pbv', 'pba', 'pbrot', ... 
    'pcd', 'pwd', 'pcrot', 'pwrot', ...
    'purot', 'pcha', 'pcva', 'drdis', 'drcase', 'gm', 'sf', 'freq'});

save("SummTable", "SummTable")
% [DataSet] = LoadData(TheFiles(File).name);


%% In code function to extract selected response parameters

function [Data] = LoadData(FileName)%% Extract Data
    load(FileName, 'BridgeResponse', 'X', 'V', 'A', 'HL', 'Vel', 'ugddot', 'ugdot', 'NL', 'CF', 'dr_type', 'SF', 'freq');
    
    Data.HL = HL;
    Data.ugddot = ugddot;
    Data.ugdot = ugdot;
    Data.SF = SF;
    Data.freq = freq;
    Data.BM = NL;
    Data.BridgeResponse = BridgeResponse;
    Data.Vel = Vel;
    Data.TrainResponse.X = X;
    Data.TrainResponse.V = V;
    Data.TrainResponse.A = A;
    Data.CF = CF;
    Data.dr_type = dr_type;

end