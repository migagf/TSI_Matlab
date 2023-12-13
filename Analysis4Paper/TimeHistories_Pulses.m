clear, clc, close all

freq = 0.3;
sf = 1.0;

folder = "C:\Users\Miguel Gomez\Documents\PhD Files\TSI_Runs\Pulses_OB\";
filename = strcat("Pulse_OB_SF", num2str(100*sf, "%03.0f"), "_FQ", num2str(100*freq, "%03.0f"), ".mat");
address = strcat(folder, filename);

load(address)

PostProcessingGM