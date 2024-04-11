clear, clc, close all

% freq = 0.3;
% sf = 1.0;
% 
% folder = "C:\Users\Miguel Gomez\Documents\PhD Files\TSI_Runs\Pulses_OB\";
% filename = strcat("Pulse_OB_SF", num2str(100*sf, "%03.0f"), "_FQ", num2str(100*freq, "%03.0f"), ".mat");
% address = strcat(folder, filename);
% 
% load(address)
% 
% PostProcessingGM

% Exploratory analysis...
filesdir = "C:\Users\Miguel.MIGUEL-DESK\Documents\PhD Files\TSI_Runs\Runs_OB";

% Select the case to be plotted
speed_case = 1;
ground_motion = 10;
scale_case = 4;
coupling = 2;
nonlinear = 1;

filename = strcat('OB_', num2str(speed_case), '_', num2str(ground_motion, '%02.0f'), '_', num2str(scale_case), '_', num2str(coupling), '_', num2str(nonlinear),'.mat');

load(strcat(filesdir, '\', filename))
PostProcessingGM

% Ground motion 2, scale 4
% Ground motion 6, scale 5







