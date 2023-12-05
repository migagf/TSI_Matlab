%% Post Processing
tmin = 2;
tmax = 14;

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
subplot(4,1,1)
plot(tt,-X(8,:)+X(8,1),tt,-X(5,:)+X(5,1),tt,-X(2,:)+X(2,1),tb,-BridgeResponse.X_Track(2,:))
ylabel('Vertical trans. (m)'), xlabel('Time (s)')
axis([tmin tmax -0.05 0.05]), grid on
legend('Wheelset','Bogie','Car','location','eastoutside')

subplot(4,1,2), plot(tt,X(7,:),tt,X(4,:),tt,X(1,:),tb,BridgeResponse.X_Track(1,:))
ylabel('Horiz. trans. (m)'), xlabel('Time (s)')
axis([tmin tmax -0.20 0.20]), grid on
legend('Wheelset','Bogie','Car','location','eastoutside')

subplot(4,1,3), plot(tt,X(9,:),tt,X(6,:),tt,X(3,:))
ylabel('Rotation (-)'), xlabel('Time (s)')
axis([tmin tmax -0.05 0.05]), grid on
legend('Wheelset','Bogie','Car','location','eastoutside')

subplot(4,1,4)
plot(tt,Left_Cont)
hold on
plot(tt,Right_Cont), axis([tmin tmax 0 500])
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
NadalL = (movmedian(YL, 200)./movmedian(QL, 200));
NadalR = (movmedian(YR, 200)./movmedian(QR, 200));

tmin = 0;
tmax = 30;

figure()
subplot(3,1,1), plot(tt, NadalL, 'r'), xlim([tmin, tmax]), ylim([-2.0, 2.0]), ylabel('Nadal Left Wheel') %, ylim([-100 100])
subplot(3,1,2), plot(tt, NadalR, 'k'), xlim([tmin, tmax]), ylim([-2.0, 2.0]), ylabel('Nadal Right Wheel') %, ylim([-100 100])
subplot(3,1,3), plot(tt, BridgeResponse.Xtdot(1,:), tt, ugdot), xlim([tmin, tmax]), ylabel('Bridge Vel. (m/s)')


% Build scalograms
NadalL(isnan(NadalL)) = 0;
NadalR(isnan(NadalR)) = 0;

figure()
cwt(BridgeResponse.Xtdot(1, :), 2000)
figure()
cwt(NadalR, 2000)




