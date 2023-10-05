%% Animation
close all
% Use 500 for EQ
% Use 25 for pulse
cont = 1;

for i = 1:5:length(X(1,:))
    X_ref = X(:,i);
    V_ref = V(:,i);
    % Activate this one for 5-pulse TH
    F = PlotState(X_ref(7:9)',V_ref(7:9)',[dx_track(i),0,0],[0,0,0],WheelGeom_pol,RailGeom_pol,RailProps,1,Vel,delta,[],[0 0 0 0],X_ref');
    % frictioncoeff(cont) = F(1)/F(2);
    
    % Activate for prescribed motion
    %     F = ContactForce(X_ref(7:9)',V_ref(7:9)',[x_track(i),z_track(i),0],[0,0,0],WheelGeom_pol,RailGeom_pol,RailProps,0,Vel,delta,[]);
    %     frictioncoeff(cont) = F(1)/F(2);
    if (abs(X_ref(7)-dx_track(i)) > 0.1 && abs(X_ref(9)) < pi/30) || abs(X_ref(9)) > pi/2 || (abs(X_ref(7)-dx_track(i)) > 1.0 && abs(X_ref(9)) < pi/2)
        input('Analysis Paused, press any key ')
    end
    cont = cont + 1;
    % Activate this one for EQ analysis (coupled structure + train)
    % ContactForce(X_ref(7:9)',V_ref(7:9)',[Y_tr(i),Z_tr(i),phi_tr(i)],[Ydot_tr(i),Zdot_tr(i),phidot_tr(i)],WheelGeom_pol,RailGeom_pol,RailProps,1,Vel);  
end