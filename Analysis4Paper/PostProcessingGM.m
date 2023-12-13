%% Post Processing
tmin = 0;
tmax = max(trec);

%% Normal Forces
% nfig = nfig+1; 
% 
% figure(nfig)

% %% Energy Calculations
% Ek = sum(M*V.^2)/2;
% Ep = sum(K*X.^2)/2;
% plot(t(1:end-1),Ek+Ep,t(1:end-1),Ek,t(1:end-1),Ep)

%% Plot translations and rotations
figure()
subplot(3,1,1), plot(tt,X(7,:),tt,X(4,:),tt,X(1,:),tb,BridgeResponse.X_Track(1,:))
ylabel('Horiz. trans. (m)'), xlabel('Time (s)')
xlim([tmin tmax]), grid on
legend('Wheelset','Bogie','Car','location','eastoutside')

subplot(3,1,2), plot(tt,X(9,:),tt,X(6,:),tt,X(3,:))
ylabel('Rotation (-)'), xlabel('Time (s)')
xlim([tmin tmax]), grid on
legend('Wheelset','Bogie','Car','location','eastoutside')

subplot(3,1,3)
plot(tt,Left_Cont, ...
     tt, Right_Cont)
xlim([tmin tmax])
grid on
xlabel('Time (sec)'), ylabel('Contact Force (kN)')
legend('L. Wheel','R. Wheel','location','eastoutside')

% %% Plot accelerations
% nfig = nfig+1; 
% figure(nfig)
% 
% load GroundMotions\ax_track.mat
% ax_track_interp = interp1(ax_track(:,1),ax_track(:,2),tt);
% % ax_track_interp = ugddot/9.81;
% % Acceleration of the top mass (car)
% Uddot_car = A(1,:);
% plot(tt,Uddot_car/9.81+ax_track_interp), xlabel('Time (s)'), ylabel('Car Acceleration (g)'), grid on

%% Plot Moment-Rotation of Plastic Hinge
if on_bridge
    figure()
    subplot(2,2,1), plot(tb,BridgeResponse.X(1,:)), xlabel('Time (s)'), ylabel('Top Lateral Disp (m)'), grid on
    subplot(2,2,3), plot(tb,BridgeResponse.X(2,:),tb,BridgeResponse.X(3,:)), xlabel('Time (s)'), ylabel('Rotation (rad)'), legend('Plastic Hinge','Top of the Column'), grid on
    subplot(2,2,2), plot(BridgeResponse.X(2,:),BridgeResponse.Mtheta/1000), xlabel('Plastic Hinge Rotation (rad)'), ylabel('Base Moment (kN-m)'), grid on
    subplot(2,2,4), plot(-BridgeResponse.X(1,:),BridgeResponse.Mtheta/1000), xlabel('Top Lateral Disp (m)'), ylabel('Base Moment (kN-m)'), grid on
end

%% Nadal Index
close all
NadalL = (movmedian(YL, 150)./movmedian(QL, 150));
NadalR = (movmedian(YR, 150)./movmedian(QR, 150));

tmin = 4;

figure()
subplot(3,1,1), plot(tt - tmin, NadalL, 'r'), xlim([0, tmax - tmin]), ylim([-5.0, 5.0]), ylabel('Nadal Left Wheel'), hold on %, ylim([-100 100])
subplot(3,1,1), plot(tt - tmin, NadalR, 'k'), xlim([0, tmax - tmin]), ylim([-5.0, 5.0]), ylabel('Nadal Right Wheel'), grid on %, ylim([-100 100])
subplot(3,1,2), plot(tt - tmin, -BridgeResponse.Xt(1,:), tt - tmin, -ug, tt - tmin, -X(7,:)), xlim([0, tmax - tmin]), ylabel('Bridge Vel. (m/s)'), grid on
subplot(3,1,3), plot(tt - tmin, BridgeResponse.X(3,:), tt - tmin, -X(9,:)), xlim([0, tmax - tmin]), ylabel('Bridge Vel. (m/s)'), grid on


% Build scalograms
NadalL(isnan(NadalL)) = 0;
NadalR(isnan(NadalR)) = 0;

% figure()
% cwt(BridgeResponse.Xtdot(1, :), 2000)
% figure()
% cwt(NadalR, 2000)




