function [C] = DampingMatrix()
% Damping Matrix

Hcb = 50*2.54/100;  % 50 (in) but in (m)
Hbt = 15*2.54/100;  % 14 (in) but in (m)
Htw =  5*2.54/100;   %  1 (in) but in (m)

dtc = 38.9685*2.54/100;  % Half distance from springs in secondary suspension system
dpc = 38.9685*2.54/100;  % Half distance from springs in primary susp. system

Cty = (17068*2)*2;     % Damping coefficient of secondary suspension along Y axis (N-s/m)
Ctz = (17068*2)*2;     % Damping coefficient of secondary suspension along Z axis (N-s/m)

Cpy = (17068*2)*2;     % Damping coefficient of primary suspension along Z axis (N-s/m)
Cpz = (17068*2)*2;     % Damping coefficient of primary suspension along Z axis (N-s/m)

%% Damping Matrix Assemble

c11 = 2*Cty;
c21 = 0;
c31 = -2*Cty*Hcb;
c41 = -2*Cty;
c51 = 0;
c61 = -2*Cty*Hbt;
c71 = 0;
c81 = 0;
c91 = 0;

c12 = c21;
c22 = 2*Ctz;
c32 = 0; 
c42 = 0;
c52 = -2*Ctz;
c62 = 0;
c72 = 0;
c82 = 0;
c92 = 0;

c13 = c31;
c23 = c32;
c33 = 2*Ctz*dtc^2 + 2*Cty*Hcb^2;
c43 = 2*Cty*Hcb;
c53 = 0;
c63 = 2*Cty*Hcb*Hbt - 2*Ctz*dtc^2;
c73 = 0;
c83 = 0;
c93 = 0;

c14 = c41;
c24 = c42;
c34 = c43;
c44 = 2*Cty + 2*Cpy;
c54 = 0;
c64 = 2*Cty*Hbt - 2*Cpy*Htw;
c74 = -2*Cpy;
c84 = 0;
c94 = 0;

c15 = c51;
c25 = c52;
c35 = c53;
c45 = c54;
c55 = 2*Ctz + 2*Cpz;
c65 = 0;
c75 = 0;
c85 = -2*Cpz;
c95 = 0;

c16 = c61;
c26 = c62;
c36 = c63;
c46 = c64;
c56 = c65;
c66 = 2*Cty*Hbt^2 + 2*Cpy*Htw^2 + 2*Ctz*dtc^2 + 2*Cpz*dpc^2;
c76 = 2*Cpy*Htw;
c86 = 0;
c96 = -2*Cpz*dpc^2;

c17 = c71;
c27 = c72;
c37 = c73;
c47 = c74;
c57 = c75;
c67 = c76;
c77 = 2*Cpy;
c87 = 0;
c97 = 0;

c18 = c81;
c28 = c82;
c38 = c83;
c48 = c84;
c58 = c85;
c68 = c86;
c78 = c87;
c88 = 2*Cpz;
c98 = 0;

c19 = c91;
c29 = c92;
c39 = c93;
c49 = c94;
c59 = c95;
c69 = c96;
c79 = c97;
c89 = c98;
c99 = 2*Cpz*dpc^2;

C = [c11 c12 c13 c14 c15 c16 c17 c18 c19;
     c21 c22 c23 c24 c25 c26 c27 c28 c29;
     c31 c32 c33 c34 c35 c36 c37 c38 c39;
     c41 c42 c43 c44 c45 c46 c47 c48 c49;
     c51 c52 c53 c54 c55 c56 c57 c58 c59;
     c61 c62 c63 c64 c65 c66 c67 c68 c69;
     c71 c72 c73 c74 c75 c76 c77 c78 c79;
     c81 c82 c83 c84 c85 c86 c87 c88 c89;
     c91 c92 c93 c94 c95 c96 c97 c98 c99];

end
