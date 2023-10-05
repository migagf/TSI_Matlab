%% Post Processing
%

t = 0:dt:(total/dt-1)*dt;
% Plot vertical translation of the wheelset
figure(2)
subplot(3,1,1)
plot(t,-X(8,:)+X(8,1),t,-X(5,:)+X(5,1),t,-X(2,:)+X(2,1)), ylabel('Vertical trans. (m)'), xlabel('Time (s)');
legend('Wheelset','Bogie','Car')
% Plot lateral translation of the wheelset
subplot(3,1,2), plot(t,X(7,:),t,X(4,:),t,X(1,:),t,x_track(1:end-2)), ylabel('Horiz. trans. (m)'), xlabel('Time (s)');

legend('Wheelset','Bogie','Car')
% Plot rotation of the wheelset
subplot(3,1,3), plot(t,X(9,:),t,X(6,:),t,X(3,:)), ylabel('Rotation (-)'), xlabel('Time (s)');
legend('Wheelset','Bogie','Car')


