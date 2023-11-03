% Run Ground Motion
clear, clc, close

%% Load Data
load WheelGeom
load RailGeom
load RailProps

%% Solution Parameters
dt = 1.0e-4;
psi = 1/2;
phi = 1/2;

g = 9.81;
vel = 0.01/3.6;

coupled = true;

%% Ground motion

load("GMs/Used Records/RSN181_IMPVALLH1.mat")










