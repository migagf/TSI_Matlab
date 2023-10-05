%% Post Processing
t = 0:dt:(total/dt-1)*dt;
%% Normal Forces
nfig = nfig+1; 

figure(nfig)
plot(t,Left_Cont)
hold on
plot(t,Right_Cont)
grid on
xlabel('Time (sec)'), ylabel('Contact Force (kN)')

% %% Energy Calculations
% Ek = sum(M*V.^2)/2;
% Ep = sum(K*X.^2)/2;
% plot(t(1:end-1),Ek+Ep,t(1:end-1),Ek,t(1:end-1),Ep)

%% Plot translations and rotations
nfig = nfig+1; 
figure(nfig)

subplot(3,1,1)
plot(t,-X(8,:)+X(8,1),t,-X(5,:)+X(5,1),t,-X(2,:)+X(2,1),t,-z_track(1:end-1)), ylabel('Vertical trans. (m)'), xlabel('Time (s)');
legend('Wheelset','Bogie','Car')
subplot(3,1,2), plot(t,X(7,:),t,X(4,:),t,X(1,:),t,x_track(1:end-1)), ylabel('Horiz. trans. (m)'), xlabel('Time (s)');
legend('Wheelset','Bogie','Car')
subplot(3,1,3), plot(t,X(9,:),t,X(6,:),t,X(3,:)), ylabel('Rotation (-)'), xlabel('Time (s)');
legend('Wheelset','Bogie','Car')

%% Plot accelerations
nfig = nfig+1; 
figure(nfig)

load GroundMotions\ax_track.mat
ax_track_interp = interp1(ax_track(:,1),ax_track(:,2),t);
% Acceleration of the top mass (car)
Uddot_car = A(1,:);
plot(t,Uddot_car/9.81), xlabel('Time (s)'), ylabel('Car Acceleration (g)'), grid on


