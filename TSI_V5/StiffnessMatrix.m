function [K] = StiffnessMatrix()
% StiffnessMatrix creates the stiffness matrix of the train

% ----------
% Geometry
% ----------

Hcb = 52*2.54/100;  % 50 (in) but in (m)
Hbt = 15*2.54/100;  % 14 (in) but in (m)
Htw = 5*2.54/100;   %  1 (in) but in (m)

dtk = 21.4685*2.54/100;  % Half distance from springs in secondary suspension system
dpk = 21.4685*2.54/100;  % Half distance from springs in primary susp. system

% ----------
% Spring Constants
% ----------

Kty = 2*0.1019*2;  % Stiffness coefficient of primary suspension along Y axis (MN/m)
Ktz = 4*0.1664*2*4;  % Stiffness coefficient of primary suspension along Z axis (MN/m)

Kpy = 2*0.1019*2;      % Stiffness coefficient of secondary suspension along Y axis (MN/m)
Kpz = 4*0.1664*2;      % Stiffness coefficient of secondary suspension along Z axis (MN/m)

%% Stiffness Matrix Coefficients

k11 = 2*Kty;
k21 = 0;
k31 = -2*Kty*Hcb;
k41 = -2*Kty;
k51 = 0;
k61 = -2*Kty*Hbt;
k71 = 0;
k81 = 0;
k91 = 0;

k12 = k21;
k22 = 2*Ktz;
k32 = 0; 
k42 = 0;
k52 = -2*Ktz;
k62 = 0;
k72 = 0;
k82 = 0;
k92 = 0;

k13 = k31;
k23 = k32;
k33 = 2*Ktz*dtk^2 + 2*Kty*Hcb^2;
k43 = 2*Kty*Hcb;
k53 = 0;
k63 = 2*Kty*Hcb*Hbt - 2*Ktz*dtk^2;
k73 = 0;
k83 = 0;
k93 = 0;

k14 = k41;
k24 = k42;
k34 = k43;
k44 = 2*Kty + 2*Kpy;
k54 = 0;
k64 = 2*Kty*Hbt - 2*Kpy*Htw;
k74 = -2*Kpy;
k84 = 0;
k94 = 0;

k15 = k51;
k25 = k52;
k35 = k53;
k45 = k54;
k55 = 2*Ktz + 2*Kpz;
k65 = 0;
k75 = 0;
k85 = -2*Kpz;
k95 = 0;

k16 = k61;
k26 = k62;
k36 = k63;
k46 = k64;
k56 = k65;
k66 = 2*Kty*Hbt^2 + 2*Kpy*Htw^2 + 2*Ktz*dtk^2 + 2*Kpz*dpk^2;
k76 = 2*Kpy*Htw;
k86 = 0;
k96 = -2*Kpz*dpk^2;

k17 = k71;
k27 = k72;
k37 = k73;
k47 = k74;
k57 = k75;
k67 = k76;
k77 = 2*Kpy;
k87 = 0;
k97 = 0;

k18 = k81;
k28 = k82;
k38 = k83;
k48 = k84;
k58 = k85;
k68 = k86;
k78 = k87;
k88 = 2*Kpz;
k98 = 0;

k19 = k91;
k29 = k92;
k39 = k93;
k49 = k94;
k59 = k95;
k69 = k96;
k79 = k97;
k89 = k98;
k99 = 2*Kpz*dpk^2;

K = [k11 k12 k13 k14 k15 k16 k17 k18 k19;
     k21 k22 k23 k24 k25 k26 k27 k28 k29;
     k31 k32 k33 k34 k35 k36 k37 k38 k39;
     k41 k42 k43 k44 k45 k46 k47 k48 k49;
     k51 k52 k53 k54 k55 k56 k57 k58 k59;
     k61 k62 k63 k64 k65 k66 k67 k68 k69;
     k71 k72 k73 k74 k75 k76 k77 k78 k79;
     k81 k82 k83 k84 k85 k86 k87 k88 k89;
     k91 k92 k93 k94 k95 k96 k97 k98 k99];

K = K*1.0E6;

end
