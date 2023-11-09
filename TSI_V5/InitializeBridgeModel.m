% Initialize Structural Model 

% ------------------------------------------------------------------------
% Structure's Parameters
% ------------------------------------------------------------------------

BridgePar.E   = 44.8e9;      % Young's Modulus of Concrete (Pa) (GPa)
BridgePar.I   = 0.9*0.1825;  % Moment of Inertia (m4)                                       
BridgePar.L   = 6.0;         % Length of Column (m)

BridgePar.mass = 262914;      % Mass of structure (kg)
BridgePar.itheta = 2075226;   % Rotational Moment of Inertia of Structure (kg-m2)

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

% Convergence Parameters

PlasticHingePar.tolerance  = 0.001;
PlasticHingePar.maxNumIter = 1000;

% ------------------------------------------------------------------------
% Input time history of displacements
% ------------------------------------------------------------------------

% Prellocation of parameters
BridgeResponse.eps = zeros(1,length(tb));
BridgeResponse.z = zeros(1,length(tb));
BridgeResponse.Mtheta = zeros(1,length(tb));
BridgeResponse.ktheta = zeros(1,length(tb));
BridgeResponse.theta1 = zeros(1,length(tb));
BridgeResponse.theta1dot = zeros(1,length(tb));

% Initial Parameters
BridgeResponse.ktheta(1) = PlasticHingePar.ko; 

% ------------------------------------------------------------------------
% Stiffness/Mass/Damping
% ------------------------------------------------------------------------

% Form Stiffness Matrix
[k_bridge, ks_bridge] = kbridge(BridgePar.E, BridgePar.I, BridgePar.L, PlasticHingePar.ko);

% Form Mass Matrix
m_bridge = diag([BridgePar.mass, BridgePar.itheta/10, BridgePar.itheta]);

% Form Damping Matrix
zeta = 0.01;

% Calculate Mass and Stiffness Proportional Damping
w2 = eig(m_bridge \ k_bridge); w2 = w2(2:3);

wn = sqrt(w2); 
w1 = wn(1); T1 = 2 * pi / w1;
w2 = wn(2); T2 = 2 * pi / w2;

a1 = zeta * 2.0 * w1 * w2 / (w1 + w2);	% mass damping coefficient based on first and second modes
a0 = zeta * 2.0 / (w1 + w2);            % stiffness damping coefficient based on first and second modes

c_bridge = a0 * k_bridge + a1 * m_bridge;                   % Some Stiffness-Proportional Damping

% ----------------------------
% Newmark's Integration Method
% ----------------------------

% Integration Parameters

IntegrationPar.Nbeta  = 1/4;
IntegrationPar.Ngamma = 1/2;

% Calculate Integration Constants

IntegrationPar.c0 = 1 / (IntegrationPar.Nbeta * dtb ^ 2);
IntegrationPar.c1 = 1 / (IntegrationPar.Nbeta * dtb);
IntegrationPar.c2 = IntegrationPar.Ngamma / (IntegrationPar.Nbeta * dtb);
IntegrationPar.c3 = 1 / (2 * IntegrationPar.Nbeta) - 1;
IntegrationPar.c4 = IntegrationPar.Ngamma / IntegrationPar.Nbeta - 1;
IntegrationPar.c5 = IntegrationPar.Ngamma /(2 * IntegrationPar.Nbeta) - 1;


% Force collocation vector
BridgePar.r  = [1, 0, 0, 1, 0, 0, 1, 0, 0]';    % 9 dof model


% Initial Conditions
BridgeResponse.X = zeros(9,length(tb));
BridgeResponse.Xdot = zeros(9,length(tb));
BridgeResponse.Xddot = zeros(9,length(tb));
BridgeResponse.Peff = zeros(9,length(tb));
BridgeResponse.Pr0 = zeros(9,length(tb));
BridgeResponse.X_Track = zeros(2,bsteps);
BridgeResponse.V_Track = zeros(2,bsteps);


%load BridgeIC % Load Initial Conditions of Bridge
%ApplyBridgeIC % Inserts the initial condition at t = 0 for each response par.


% k_rail(kh, kv, ko, l2, l3)
% l2: distance from the center of pier to the top of pier
% l3: distance from the center of pier to the base of the rail

kh = 1.0e6;        % N/m
kv = 1.0e7;        % N/m
ktheta = 1.0e7;    % N-m/rad

ch = 1.0e7;        % N/m
cv = 1.0e7;        % N/m
ctheta = 1.0e7;    % N-m/rad

krail1 = k_rail(kh, kv, ktheta, 0.5,  0.2);
krail2 = k_rail(kh, kv, ktheta, 0.5,  0.8);

crail1 = k_rail(ch, cv, ctheta, 0.5,  0.2);
crail2 = k_rail(ch, cv, ctheta, 0.5,  0.8);

kfull = fullstiffness(k_bridge, krail1, krail2);
ksfull = fullstiffness(ks_bridge, krail1, krail2);
cfull = fullstiffness(c_bridge, crail1, crail2);

mfull = diag([ ...
    BridgePar.mass, BridgePar.itheta/10, BridgePar.itheta, ...
    BridgePar.mass/100, BridgePar.mass/100, BridgePar.itheta/100, ...
    BridgePar.mass/100, BridgePar.mass/100, BridgePar.itheta/100]);

% Here, pay attention to the mass and inertia of the track assemblies
%[phi, tn, xin] = eigenvalue(mfull, kfull, cfull);


BridgePar.M = mfull;
BridgePar.Kt = kfull;
BridgePar.Kts = ksfull;
BridgePar.C = cfull;

% ------------------------------------------------------------------------
% Functions
% ------------------------------------------------------------------------

function [K, Ks] = kbridge(E,I,L,ktheta)

    k11 = 12*E*I/L^3;  k12 = 6*E*I/L^2;       k13 = 6*E*I/L^2;
    k21 = k12;         k22 = 4*E*I/L+ktheta;  k23 = 2*E*I/L;
    k31 = k13;         k32 = k23;             k33 = 4*E*I/L;
    
    k22s = 4*E*I/L;

    K = [k11 k12 k13; k21 k22 k23; k31 k32 k33];
    Ks = [k11 k12 k13; k21 k22s k23; k31 k32 k33];
end


function [k_r] = k_rail(kh, kv, ko, l2, l3)
    % K_RAIL
    % This function creates the stiffness matrix of one of the rails
    % and takes it into the global coordinate system
    
    k_e = [
        kh, 0, 0, -kh, 0, 0;
        0, kv, 0, 0, -kv, 0;
        0, 0, ko, 0, 0, -ko;
        -kh, 0, 0, kh, 0, 0;
        0, -kv, 0, 0, kv, 0;
        0, 0, -ko, 0, 0, ko
        ];

    Tmatrix = [
        1, -l2, 0, 0, 0;
        0, -l3, 0, 0, 0;
        0, 1, 0, 0, 0;
        0, 0, 1, 0, 0;
        0, 0, 0, 1, 0;
        0, 0, 0, 0, 1
        ];

    k_r = Tmatrix' * k_e * Tmatrix;
    
end


function [k_full] = fullstiffness(k_b, k_r1, k_r2)
    % FULLSTIFFNESS
    % This function assembles the stiffness matrix of the whole 
    % structure

    ndof_b = size(k_b, 1);           % Get the size of the structure's stiffness matrix
    k_f1 = zeros(ndof_b + 6);     % Added 6 dofs for the 
    k_f2 = zeros(ndof_b + 6);
    k_f3 = zeros(ndof_b + 6);

    i_b = [1, 2, 3];
    i_r1 = [1, 3, 4, 5, 6];
    i_r2 = [1, 3, 7, 8, 9];
    
    k_f1(i_b, i_b) = k_b;
    k_f2(i_r1, i_r1) = k_r1;
    k_f3(i_r2, i_r2) = k_r2;

    k_full = k_f1 + k_f2 + k_f3;
    
end


function [m_full] = fullmass(m_b, m_r1, m_r2)
    % FULLSTIFFNESS
    % This function assembles the stiffness matrix of the whole 
    % structure

    ndof_b = size(k_b, 1);           % Get the size of the structure's stiffness matrix
    m_full = zeros(ndof_b + 6);     % Added 6 dofs for the 
    
    i_b = [1, 2, 3];
    i_r1 = [1, 3, 4, 5, 6];
    i_r2 = [1, 3, 7, 8, 9];
    
    m_full(i_b, i_b) = m_b;
    m_full(i_r1, i_r1) = m_r1;
    m_full(i_r2, i_r2) = m_r2;

end


function [phi, tn, xin] = eigenvalue(M, K, C)
    [phi, w2] = eig(M \ K);

    w2 = diag(w2);
    wn = sqrt(w2);
    tn = 2 * pi ./ wn;
    
    % Modal masses
    Mn = diag(phi' * M * phi);
    Cn = diag(phi' * C * phi);

    xin = Cn./(2 * wn .* Mn);

end





