% RUN MULTI CASES
% Runs the time history analysis for multiple ground motions

clear, clc, close all
nfig = 0;

% Select showplot = 1 if want to plot the animation with the response
showplot = 0;

% true = over bridge | false = at grade -- AG still has structure rubber pads
on_bridge = false;

% Load all Ground Motions and the scale factors
AllGMs = dir("GMs\UsedRecords\*.mat");
load("GMs\ScaleFactors.mat");

NL_vals = [1.0 0.0];     % Linear/Nonlinear Analysis Cases [1 = nonlinear | 0 = linear]
CF_vals = [1.0 0.0];     % Coupling Factor Cases [1 = coupled | 0 = decoupled]
SP_vals = [0.01 40 80];  % Speed of train (kph)
SF_vals = [0.21 0.74 1.0 1.38 2.0 3.0 4.0 5.0];  % "Hazard" Cases

%% Load Data
load WheelGeom
load RailGeom
load RailProps
RailProps.E1 = RailProps.E1;
RailProps.E1 = RailProps.E1;

%% Load Train Data
CreateWheelRailGeom

%% Solution Parameters
dtt = 5.0E-4;   % time step (sec) of train
dtb = 5.0E-4;   % time step (sec) of bridge 

% Train Finite Difference Sol. Parameters
psi = 0.5;
phi = 0.5;
g = 9.81;       % m/s2

for SP_index = 1:3   % Speed
    for GM_index = 1:20   % GM index
        for SF_index = 1:8   % Scale Factor
            for CF_index = 1:1   % Coupling Factor (Here is just one, because no coupled analysis)
                for NL_index = 1:1   % Linear/Nonlinear
                    
                    % (1) Generate analysis code
                    AnalysisCode = strcat('AG_',num2str(SP_index),'_',num2str(GM_index,"%02.f"),'_',num2str(SF_index),'_',num2str(CF_index),'_',num2str(NL_index));
                    
                    % (2) Factor for LE Bridge or NL Bridge
                    NL = NL_vals(NL_index);           % NL = 1 --> Nonlinear | NL = 0 --> Linear
                    Vel = SP_vals(SP_index)/3.6;      % m/sec (turn speed from kph to m/s)
                    CF = CF_vals(CF_index);           % Coupling factor
                    HL = SF_vals(SF_index);           % Hazard level factor
                    SF = HL * ScaleFactors(GM_index); % Scale factor for the analysis
                    
                    % (3) Load ground motion from suite and print name of running GM
                    load("GMs\UsedRecords\"+AllGMs(GM_index).name);
                    disp(AllGMs(GM_index).name)
                    disp('SP_GM_SF_CF_NL')
                    disp(AnalysisCode)
                    disp('of_3_20_8_2_2')
                    
                    % (4) Proceessing of the ground motion
                    ugddot = SF*TimeAccelData(:,2)*9.81;         % Displacement Time-History (m)
                    dtrec = round(0.005,3);                      % Time step of record
                    ugddot = [0*(0:dtrec:2)'; ugddot];
                    trec  = 0:dtrec:dtrec*(length(ugddot)-1);    % Time vector 
                    
                    %ugdot  = [0 diff(urec')]/(dtrec);           % Velocity Time-History
                    %ugddot = [0 diff(ugdot)]/(dtrec);           % Acceleration Time-History
                    
                    ugdot = cumtrapz(ugddot) * dtrec;
                    urec  = cumtrapz(ugdot) * dtrec;
                    
                    % Plot EQ record
                    nfig = nfig + 1;
                    
                    % figure(nfig) 
                    % subplot(3,1,1), plot(trec,ugddot/9.81), xlabel('Time (sec)'), ylabel('Accel. (g)'), title('Ground Motion (LPE)') 
                    % subplot(3,1,2), plot(trec,ugdot), xlabel('Time (sec)'), ylabel('Vel. (m/s)'), title('Ground Motion (LPE)')
                    % subplot(3,1,3), plot(trec,urec), xlabel('Time (sec)'), ylabel('Dis. (m)'), title('Ground Motion (LPE)')
                    
                    % Resampling of Acceleration Time History
                    tt = 0:dtt:trec(end);  % Time vector for train analysis
                    tb = 0:dtb:trec(end);  % Time vector for bridge analysis
                    
                    ugddot = interp1(trec,ugddot,tb);  % Excitation is required only for bridge
                    ugdot = interp1(trec,ugdot,tb);   % Excitation is required only for bridge
                    urec = interp1(trec,urec,tb);    % Excitation is required only for bridge
                    
                    ug = urec;
                    
                    tsteps = length(tt);
                    bsteps = length(tb);
                    
                    % (5) Initialization of the model

                    InitializeBridgeModel
                    
                    % Initial conditions and prellocation of variables
                    
                    X = zeros(9,tsteps); % Train Global Coordinates 
                    V = zeros(9,tsteps); % Train Velocities
                    A = zeros(9,tsteps); % Train Accelerations
                    
                    % Deprecated
                    % load X_initial
                    % X0 = X_initial;      % Initial Global Coordinates
                    % X0(1) = 0; X0(4) = 0; X0(7) = 0.0;
                    
                    X0 = 0.5655 * [0 -1 0 0 -1 0 0 -1 0]';
                    V0 = [0 0 0 0 0 0 0 0 0]';  % Initial Velocities
                    A0 = [0 0 0 0 0 0 0 0 0]';  % Initial Accelerations
                    
                    X1 = X0;                    % Initial Global Coordinates
                    V1 = [0 0 0 0 0 0 0 0 0]';  % Initial Velocities
                    A1 = [0 0 0 0 0 0 0 0 0]';  % Initial Accelerations
                    
                    X(:,1) = X0; X(:,2) = X1;
                    V(:,1) = V0; V(:,2) = V1;
                    A(:,1) = A0; A(:,2) = A1;
                    
                    BridgeResponse.X_Track = zeros(2,bsteps);
                    BridgeResponse.V_Track = zeros(2,bsteps);
                    
                    % Extra output variables
                    Momt = [0 0 0 0]';
                    Creepforces = zeros(4, tsteps);
                    delta       = zeros(4, tsteps);
                    deltadotn_vec = [0 0 0 0];
                    Left_Cont = zeros(1, tsteps);
                    Right_Cont =zeros(1, tsteps);
                    
                    % Contact forces (storage)
                    QL = zeros(1, tsteps);
                    QR = zeros(1, tsteps);
                    YL = zeros(1, tsteps);
                    YR = zeros(1, tsteps);
                    
                    % (6) Mass, Stiffness and Damping Matrices                   
                    [K] = StiffnessMatrix(); 
                    [C] = DampingMatrix();
                    [M, Mc, Mt, Mw] = MassMatrix();
                    
                    AddRayleighDamp = true;
                    if AddRayleighDamp
                        w1 = 2 * pi * 100; %rad/sec
                        a1 = 2/w1;
                        C = C + a1*K;
                    end
                    
                    % EigenValue Analysis
                    % Options for EigenValue Analysis
                    % showplot = false;
                    % EigenValueAnalysis;
                    
                    % (7) Apply Gravity Loading
                    F_ine = [0 Mc*g 0 0 Mt*g 0 0 Mw*g 0]'; % Inertial Forces (N)
                    
                    % (8) Solution of the EOM
                    tic
                    
                    ib = 2; tbridge = 0;
                    
                    for it = 2:tsteps
                        
                        % Integration of Train EOM
                        X2 = X1 + V1 * dtt + (0.5 + psi) * A1 * dtt ^ 2 - psi * A0 * dtt ^ 2;
                        V2 = V1 + (1 + phi) * A1 * dtt - phi * A0 * dtt;
                        
                        try
                            % Update the force vector with Contact Algorithm
                            [F,NF_L,NF_R,vec,Ft,delta(:,it-1),Momt] = ...
                                ContactForce( ...
                                X2(7:9)', ...
                                V2(7:9)', ...
                                [BridgeResponse.X_Track(1, ib-1), BridgeResponse.X_Track(2, ib-1), -BridgeResponse.X(3, ib-1)], ...
                                BridgeResponse.Xt(end-5:end,ib-1), ...
                                BridgeResponse.Xtdot(end-5:end,ib-1), ...
                                WheelGeom_pol, ...
                                RailGeom_pol, ...
                                RailGeom_cur, ...
                                WheelGeom_cur, ...
                                RailProps, ...
                                showplot, ...
                                Vel, ...
                                delta(:,it-1), ...
                                Momt, ...
                                deltadotn_vec, ...
                                Rail_geom_fine, ...
                                LWheel_geom_fine, ...
                                RWheel_geom_fine ...
                                );
                            % If everything works ok, then no derailment
                            pderail = 0;
                        catch
                            pderail = 1.0;
                            disp('Problem in the integration scheme. Derailment is likely')
                            % break
                        end
                        
                        % Force Vector to Structure
                        Cont_Force = 10 ^ 6 * [0, 0, 0, -F(1, 1), -F(2, 1), -F(1, 1) * 0.1, -F(1, 2), -F(2, 2), -F(1, 2) * 0.1]';
                    
                        % Force vector to train
                        Cforce = 10 ^ 6 * [0, 0, 0, 0, 0, 0, F(1, 1) + F(1, 2), F(2, 1) + F(2, 2), F(3, 1) + F(3, 2)]';
                        
                        F_ext =  F_ine + Cforce; % N
                        
                        tbridge = tbridge + dtt;
                    
                        if tbridge >= dtb
                            BridgeTimeIncrement
                            UpdateX_Track;
                            ib = ib + 1;
                            tbridge = 0;
                        end
                        
                        % Variable storage
                        % vrel(i) = V2(7)-vx_track(i); vect(i,:) = Ft;
                        Left_Cont(it) = 1000 * NF_L; 
                        Right_Cont(it) = 1000 * NF_R;
                    
                        % Contact forces (storage)
                        QL(it) = F(2, 1);
                        QR(it) = F(2, 2);
                        YL(it) = F(1, 1);
                        YR(it) = F(1, 2);
                    
                        %Numerical integration of Train EOM
                        A2 = M \ (F_ext - K * X2 - C * V2);
                        
                        X(:,it) = X2; V(:,it) = V2; A(:,it) = A2;
                        
                        % Update the state vectors
                        X0 = X1; V0 = V1; A0 = A1;
                        X1 = X2; V1 = V2; A1 = A2;
                    
                    end
                    
                    close all

                    toc
                    disp("-------")
                    
                    save(strcat("Results\Runs_AG\", AnalysisCode))
                    % PostProcessingGM
                    % This is here to stop the code whenever needed
                    % myvar = input('code is paused');
                end
            end
        end
    end
end

disp('All Done')

%% Post Processing

% PostProcessingGM

% hold off
% figure(3)
% plot(t,vrel) %Relative velocity between wheelset and rail
% hold on
% plot(t,vect(:,1)) % Tangential force in the Y direction
% plot(t(1:end),Left_Cont/1000) % Left contact force
% plot(t(1:end),Right_Cont/1000) % Right contact force
% plot(t(1:end),A(7,:)/981)
% grid on
% legend('Relative velocity','Tangential force','Left normal force','Right normal force')
