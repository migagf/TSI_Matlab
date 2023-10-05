% Initialize Structural Model 
%
% This code creates the objects for running an analysis with a structure

% ------------------------------------------------------------------------
% Structure's Parameters
% ------------------------------------------------------------------------

BridgePar.E = 44.8e9;         % Young's Modulus of Concrete (Pa) (GPa)
BridgePar.I = 0.9 * 0.1825;   % Moment of Inertia (m4)                                       
BridgePar.L = 6.0;            % Length of Column (m)

BridgePar.m = 262914;         % Mass of structure (kg)
BridgePar.Itheta = 2075226;   % Rotational Moment of Inertia of Structure (kg-m2)

% ------------------------------------------------------------------------
% Plastic Hinge Parameters (BWBN)
% ------------------------------------------------------------------------

PlasticHingePar.Fy       = 38432000;
PlasticHingePar.thetay   = 0.001;
PlasticHingePar.xu       = 10*PlasticHingePar.thetay;

PlasticHingePar.alpha    = 0.007;
PlasticHingePar.ko       = PlasticHingePar.Fy/PlasticHingePar.thetay;
PlasticHingePar.n        = 1.28;
PlasticHingePar.eta      = 0.0;
PlasticHingePar.beta     = 17000;

PlasticHingePar.rhoeps   = 0.0;
PlasticHingePar.rhox     = 0.25; % 0.11
PlasticHingePar.phi      = 1.6;
PlasticHingePar.deltak   = 15.0; % 47.0 Real = 40
PlasticHingePar.deltaf   = 0.50;  % 0.5

PlasticHingePar.sigma = 0.01;
PlasticHingePar.u     = 4.0;
PlasticHingePar.epsp  = 0.001;
PlasticHingePar.rhop  = 0.1;

PlasticHingePar.xmax  = 0;
PlasticHingePar.xmaxp = 0;

% ------------------------------------------------------------------------
% Interface Model Parameters
% ------------------------------------------------------------------------

kh = 10;    % Lateral Stiffness
kv = 10;    % Vertical Stiffness
kr = 100;   % Rotational Stiffness

% ------------------------------------------------------------------------
% Convergence Parameters
% ------------------------------------------------------------------------

PlasticHingePar.tolerance  = 0.001;
PlasticHingePar.maxNumIter = 1000;

% ------------------------------------------------------------------------
% Input time history of displacements
% ------------------------------------------------------------------------

% Prellocation of parameters
BridgeResponse.eps       = zeros(1,length(tb));
BridgeResponse.z         = zeros(1,length(tb));
BridgeResponse.Mtheta    = zeros(1,length(tb));
BridgeResponse.ktheta    = zeros(1,length(tb));
BridgeResponse.theta1    = zeros(1,length(tb));
BridgeResponse.theta1dot = zeros(1,length(tb));

% Initial Parameters
BridgeResponse.ktheta(1) = PlasticHingePar.ko; 

% ------------------------------------------------------------------------
% Stiffness/Mass/Damping
% ------------------------------------------------------------------------

% Form Stiffness Matrix
BridgePar.Kt = InitialKt(BridgePar.E,BridgePar.I,BridgePar.L,PlasticHingePar.ko);

% Form Mass Matrix
BridgePar.M = diag([BridgePar.m 1.0 BridgePar.Itheta]);

% Form Damping Matrix
zeta = 0.01;

% Calculate Mass and Stiffness Proportional Damping
[w2] = eig(BridgePar.M\BridgePar.Kt); w2 = w2(2:3);
wn = sqrt(w2); 
w1 = wn(1); T1 = 2*pi/w1;
w2 = wn(2); T2 = 2*pi/w2;
a1 = zeta*2.0*w1*w2/(w1 + w2);	% mass damping coefficient based on first and second modes
a0 = zeta*2.0/(w1 + w2);        % stiffness damping coefficient based on first and second modes

BridgePar.C = a0*BridgePar.Kt + a1*BridgePar.M; % Add Some Stiffness-Proportional Damping

%% Newmark's Integration Method

% Integration Parameters

IntegrationPar.Nbeta  = 1/4;
IntegrationPar.Ngamma = 1/2;

% Calculate Integration Constants

IntegrationPar.c0 = 1/(IntegrationPar.Nbeta*dtb^2);
IntegrationPar.c1 = 1/(IntegrationPar.Nbeta*dtb);
IntegrationPar.c2 = IntegrationPar.Ngamma/(IntegrationPar.Nbeta*dtb);
IntegrationPar.c3 = 1/(2*IntegrationPar.Nbeta)-1;
IntegrationPar.c4 = IntegrationPar.Ngamma/IntegrationPar.Nbeta-1;
IntegrationPar.c5 = IntegrationPar.Ngamma/(2*IntegrationPar.Nbeta)-1;
BridgePar.r  = [1 0 0]';

% Effective Stiffness Matrix
BridgePar.Keff = IntegrationPar.c0*BridgePar.M+IntegrationPar.c2*BridgePar.C+BridgePar.Kt;

% Initial Conditions
BridgeResponse.X     = zeros(3,length(tb));
BridgeResponse.Xdot  = zeros(3,length(tb));
BridgeResponse.Xddot = zeros(3,length(tb));
BridgeResponse.Peff  = zeros(3,length(tb));
BridgeResponse.Pr0   = zeros(3,length(tb));
BridgeResponse.X_Track = zeros(2,bsteps);
BridgeResponse.V_Track = zeros(2,bsteps);

load BridgeIC % Load Initial Conditions of Bridge
ApplyBridgeIC % Inserts the initial condition at t = 0 for each response par.

function [K] = InitialKt(E,I,L,ktheta)
    k11 = 12*E*I/L^3;  k12 = 6*E*I/L^2;       k13 = 6*E*I/L^2;
    k21 = k12;         k22 = 4*E*I/L+ktheta;  k23 = 2*E*I/L;
    k31 = k13;         k32 = k23;             k33 = 4*E*I/L;

    K = [k11 k12 k13; k21 k22 k23; k31 k32 k33];
end

function [K_full] = create_stiffness_all(ks, kh, kv, kr, l, h)
% This functions extends the structural stiffness matrix to convert it into
% the full stiffness matrix with the interface rail

K_full = [ks(1,1), ks(1,2), ks(1,3), -kh, 0, 0, -kh, 0, 0;
          ks(2,1), ks(2,2), ks(2,3), 0, kv * (l - h/2), -kr, 0, kv * (l + h/2), -kr;
          ks(3,1), ks(3,2), ks(3,3), 0, ]
end






