%% RunPrescribedMotion
% Runs wavelets on train-only system
clear, clc, close

%% Load data
load WheelGeom
load RailGeom
load RailProps

RailProps.R1y = 0.2;

RailProps.E1 = RailProps.E1/10;
RailProps.E2 = RailProps.E2/10;

%% Load Train Data
CreateWheelRailGeom

%% Solution Parameters
dt  = 1.0E-3;                            % time step (sec)
g   = 9.81;                              % m/s2
Vel = 80 / 3.6;                          % m/sec

% Scale Factors
SFvals = 0.05:0.05:5.0;

% Frequency Values
freqvals = 0.02:0.02:1.0;

initID = 1;

ID = initID;
%%
for freq = freqvals   % Factor that modifies the frequency of the pulse
for SF = SFvals

disp(['Running Analysis #',num2str(ID),' of ',num2str(length(SFvals)*length(freqvals))])

% Initial conditions

load X_initial
X0 = X_initial;                         % [0 Z1c 0 0 Z1b 0 0 Z1w 0]';  % Initial positions
                                        % X0 = [0 -0.4 0 0 -0.43 0 0 -0.5 0]' - 0.065*[0 1 0 0 1 0 0 1 0]'
V0 = [0 0 0 0 0 0 0 0 0]';              % Initial velocities
A0 = [0 0 0 0 0 0 0 0 0]';              % Initial accelerations

X1 = X0;                                %[0 Z1c 0 0 Z1b 0 0 Z1w 0]';  % 2nd initial position
V1 = [0 0 0 0 0 0 0 0 0]';              % 2nd initial velocities
A1 = [0 0 0 0 0 0 0 0 0]';              % Initial accelerations

% Stiffness and Damping Matrices

[K] = StiffnessMatrix(); 
[C] = DampingMatrix();
[M,Mc,Mt,Mw] = MassMatrix();

% AddRayleighDamp = false;
% if AddRayleighDamp
%     w1 = 2*pi*100; %rad/sec
%     a1 = 2/w1;
%     C = C + a1*K;
% end
% 
% % EigenValue Analysis
% showplot = false;
% EigenValueAnalysis;

% Forcing parameters

% Create a wavelet for analysis of derailment
wp = 2*pi*freq;
[time, accel] = RickerWave(wp);
t = 0:dt:max(time);
ag = interp1(time,accel,t) * SF * g;

ax_track = ag;
az_track = 0*ag;

vx_track = dt*cumtrapz(ax_track);
vz_track = 0*vx_track;

dx_track = dt*cumtrapz(vx_track);
dz_track = 0*dx_track;

% nfig = nfig + 2;
% figure(nfig)
% subplot(1,3,1), plot(t,ax_track)
% subplot(1,3,2), plot(t,vx_track)
% subplot(1,3,3), plot(t,dx_track)
        
%% Iterations
% Initial definitions
psi = 0.5;
phi = 0.5;

F_ine = [0 Mc*g 0 0 Mt*g 0 0 Mw*g 0]'; % N

vari = 1;
steps = length(t);
X = zeros(9,steps); V = X; A = X;
X(:,1) = X0; X(:,2) = X1;
V(:,1) = V0; V(:,2) = V1;
A(:,1) = A0; A(:,2) = A1;

Momt = [0 0 0 0]';
Creepforces = zeros(4,steps);
delta = zeros(4,length(dx_track));
deltadotn_vec = [0 0 0 0];

tic
for i = 3:steps
    
    % Integration with explicit method (Part 1)
    X2 = X1 + V1 * dt + (0.5 + psi) * A1 * dt ^ 2 - psi * A0 * dt ^ 2;
    V2 = V1 + (1 + phi) * A1 * dt - phi * A0 * dt;
    
    % Check if derailment has been reached
    if (abs(X2(7)-dx_track(i)) > 0.1 && abs(X2(9)) < pi/30) || abs(X2(9)) > pi/2 || (abs(X2(7)-dx_track(i)) > 1.0 && abs(X2(9)) < pi/2)
        break
    else
        % Update the force vector with Contact Algorithm
        [F,NF_L,NF_R,vec,Ft,delta(:,i),Momt] = ... 
            ContactForce(X2(7:9)',V2(7:9)',[dx_track(i),dz_track(i),0],...
            [vx_track(i),vz_track(i),0],WheelGeom_pol,RailGeom_pol,...
            RailGeom_cur,WheelGeom_cur,RailProps, ...
            0,Vel,delta(:,i-1),Momt,deltadotn_vec, ...
            Rail_geom_fine,LWheel_geom_fine,RWheel_geom_fine);
        
        % Updated force vector
        Cont_Force = 10^6*[0 0 0 0 0 0 F']';
        F_ext =  F_ine + Cont_Force; %  N
        
        % Variable storage
        vrel(i) = V2(7)-vx_track(i); vect(i,:) = Ft;
        Left_Cont(i) = 1000*NF_L; Right_Cont(i) = 1000*NF_R;
        
        % Numerical integration of the EOM
        A2 = M\(F_ext-K*X2-C*V2);
        X(:,i) = X2; V(:,i) = V2; A(:,i) = A2;
        
        % Update the state vectors
        X0 = X1; V0 = V1; A0 = A1;
        X1 = X2; V1 = V2; A1 = A2;
    end
end
toc
save(['WaveletCases\Wave_ID_',num2str(ID,'%04.f'),'_wF_',num2str(freq*100),'_sF_',num2str(SF*100),'.mat'])
clearvars -except WheelGeom_pol WheelGeom_cur WheelGeom Vel SFvals RWheel_geom_fine RailProps RailGeom_pol RailGeom_cur RailGeom Rail_geom_fine Rail_geom nodes nfig LWheel_geom_fine g fgrid_n dt b a ID freqvals freq
ID = ID + 1;
end
end

%%  Post Processing
% nfig = 1
% PostProcessingPM

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


%% Functions
% Ricker wavelet
function [tau, r] = RickerWave(wp)
    % Function creates a ricker wavelet with the specified frequency
    tau = linspace(-10/wp,10/wp,10000);
    r = (1-0.5*wp^2*tau.^2).*exp(-1/4*wp^2*tau.^2);
    tau = tau - min(tau);
    % plot(tau,r)
end








