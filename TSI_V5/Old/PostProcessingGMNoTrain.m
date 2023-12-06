%% Post Processing No Train

%% Normal Forces

%% Plot translations and rotations

%% Plot accelerations

%% Plot Moment-Rotation of Plastic Hinge
nfig = nfig+1; 
figure(nfig)
subplot(2,2,1), plot(tb,BridgeResponse.X(1,:)), xlabel('Time (s)'), ylabel('Top Lateral Disp (m)'), grid on
subplot(2,2,3), plot(tb,BridgeResponse.X(2,:),tb,BridgeResponse.X(3,:)), xlabel('Time (s)'), ylabel('Rotation (rad)'), legend('Plastic Hinge','Top of the Column'), grid on
subplot(2,2,2), plot(BridgeResponse.X(2,:),BridgeResponse.Mtheta/1000), xlabel('Plastic Hinge Rotation (rad)'), ylabel('Base Moment (kN-m)'), grid on
subplot(2,2,4), plot(BridgeResponse.X(1,:),-BridgeResponse.Mtheta/1000), xlabel('Top Lateral Disp (m)'), ylabel('Base Moment (kN-m)'), grid on


