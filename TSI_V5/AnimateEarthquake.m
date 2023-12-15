%% Animation
% Use 500 for EQ
% Use 25 for pulse

cont = 1;
imax = length(BridgeResponse.X_Track(1,:));

if isempty(find(abs(X(7,:)-BridgeResponse.X_Track(1,:))>2*2.54/100))
    imax = length(BridgeResponse.X_Track(1,:));
else
    imax = min(find(abs(X(7,:)-BridgeResponse.X_Track(1,:))>2*2.54/100));
end

for i = 1:10:imax
    x_track(i) = BridgeResponse.X_Track(1,i);
    z_track(i) = BridgeResponse.X_Track(2,i)+0.01;
    phi_track(i) = -BridgeResponse.X(3,i);

    X_ref = X(:,i);
    V_ref = V(:,i);
    % Activate this one for 5-pulse TH
    % F = ContactForce(X_ref(7:9)',V_ref(7:9)',[x_track(i),0,0],[0,0,0],WheelGeom_pol,RailGeom_pol,RailProps,1,Vel,delta,[])
    % frictioncoeff(cont) = F(1)/F(2);
    
    % Activate for prescribed motion

    F = PlotState(X_ref(7:9)',V_ref(7:9)',[x_track(i),z_track(i),phi_track(i)],[0,0,0],WheelGeom_pol,RailGeom_pol,RailProps,1,Vel,delta,[],[0 0 0 0],X_ref(1:6)');
    frictioncoeff(cont) = F(1)/F(2);
    cont = cont + 1;
    % Activate this one for EQ analysis (coupled structure + train)
    % ContactForce(X_ref(7:9)',V_ref(7:9)',[Y_tr(i),Z_tr(i),phi_tr(i)],[Ydot_tr(i),Zdot_tr(i),phidot_tr(i)],WheelGeom_pol,RailGeom_pol,RailProps,1,Vel);  
end

