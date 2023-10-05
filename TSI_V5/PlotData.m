%% Plot Data
% clear, clc, close

%% Plot Response Spectrum
load GMs\Response_Spectrum.txt

T  = Response_Spectrum(:,1);
Sa = Response_Spectrum(:,2);
Sd = Sa.*(T/(2*pi)).^2*981;

figure(1)
subplot(1,2,1), plot(T,Sa), title('Sa Spectrum'), xlabel('Period T (sec)'), ylabel('Pseudo Acceleration (g)'), grid on
subplot(1,2,2), plot(T,Sd), title('Sd Spectrum'), xlabel('Period T (sec)'), ylabel('Displacement (cm)'), grid on

%% Plot OpenSees Model Response

load OpenSeesRes\BridgeResponseOS

figure(2)
subplot(1,2,1), plot(BridgeResponseOS(:,1),BridgeResponseOS(:,2)), title('Displacement Time History'), xlabel('Time (sec)'), ylabel('Top Mass Displacement (m)'), grid on
subplot(1,2,2), plot(BridgeResponseOS(:,2),BridgeResponseOS(:,3)), title('Moment-Deformation'), xlabel('Top Mass Displacement (m)'), ylabel('Base Moment (kN-m)'), grid on

%% Plot OpenSees Model Response vs Matlab Model Response
% load SimResults\SimRes_010_NLBridgeTrainTogether
% load SimResults\SimRes_010_NLBridgeOnly % Original Stiffness Degradation
% load SimResults\SimRes_012_NLBridgeOnly % Modified Reduced Stiffness Degradation
% load SimResults\SimRes_010_LEBridgeOnly
% load SimResults\SimRes_011_NLBridgeTrain
% load SimResults\SimRes_001_NLB_TC % Nonlinear Bridge / Coupled / SF = 1.0
% load SimResults\SimRes_002_NLB_TC % Nonlinear Bridge / Coupled / SF = 1.8
  load SimResults\SimRes_022_NLB_TO

figure(3)
subplot(1,2,1), plot(BridgeResponseOS(:,1),BridgeResponseOS(:,2),t,BridgeResponse.X(1,:)), legend('OpenSees', 'Matlab'), title('Displacement Time History'), xlabel('Time (sec)'), ylabel('Top Mass Displacement (m)'), grid on
subplot(1,2,2), plot(BridgeResponseOS(:,2),BridgeResponseOS(:,3),-BridgeResponse.X(1,:),BridgeResponse.Mtheta/1000), legend('OpenSees', 'Matlab'), title('Moment-Deformation'), xlabel('Top Mass Displacement (m)'), ylabel('Base Moment (kN-m)'), grid on

%% Plot Comparison Train and No Train
% 
% load SimResults\SimRes_011_NLBridgeTrain
% 
% figure(3)
% 
% subplot(1,2,1), hold on, plot(t,BridgeResponse.X(1,:)), title('Displacement Time History'), xlabel('Time (sec)'), ylabel('Top Mass Displacement (m)'), grid on
% subplot(1,2,2), hold on, plot(-BridgeResponse.X(1,:),BridgeResponse.Mtheta/1000), title('Moment-Deformation'), xlabel('Top Mass Displacement (m)'), ylabel('Base Moment (kN-m)'), grid on
% 
% load SimResults\SimRes_011_NLBridgeTrainSeparate
% 
% figure(3)
% subplot(1,2,1), hold on, plot(t,BridgeResponse.X(1,:)), title('Displacement Time History'), xlabel('Time (sec)'), ylabel('Top Mass Displacement (m)'), grid on
% subplot(1,2,2), hold on, plot(-BridgeResponse.X(1,:),BridgeResponse.Mtheta/1000), title('Moment-Deformation'), xlabel('Top Mass Displacement (m)'), ylabel('Base Moment (kN-m)'), grid on
% 







