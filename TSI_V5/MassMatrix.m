function [M,Mc,Mt,Mw] = MassMatrix()% ----------
% Mass and inertia
% ----------
Mc = 29600;         % Car body mass (kg)
Mt = 2*3621;        % Bogie mass (kg)
Mw = 4*648;         % Wheelset mass (kg)

Icx = 3.2e4;        % Mass moment of inertia of car body about X axis (kg-m2)
Itx = 2*944;       % Mass moment of inertia of bogie about X axis (kg-m2)
Iwx = 4*362.41;       % Mass moment of inertia of wheelset about X axis (kg-m2)

%% Mass Matrix

M = diag([Mc Mc Icx Mt Mt Itx Mw Mw Iwx]);

end