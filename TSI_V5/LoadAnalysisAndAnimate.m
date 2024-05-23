% Load and animate analysis
analysisID = 127;

TheFiles = dir("C:\Users\Miguel.MIGUEL-DESK\Documents\PhD Files\TSI_Runs\Runs_AG\*.mat");

Folder = TheFiles(analysisID).folder;
% FileName = TheFiles(analysisID).name

FileName = "AG_1_09_6_1_1.mat";
%FileName = "AG_3_09_6_1_1.mat"
myFile = strcat(Folder,'\' ,FileName);

load(myFile)

AnimateEarthquake