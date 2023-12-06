%% Assume Bridge Response

BridgeResponse.eps(ip) = BridgeResponse.eps(ip-1);
BridgeResponse.z(ip)   = BridgeResponse.z(ip-1);

BridgeResponse.Mtheta(ip)    = BridgeResponse.Mtheta(ip-1);
BridgeResponse.ktheta(ip)    = BridgeResponse.ktheta(ip-1);
BridgeResponse.theta1(ip)    = BridgeResponse.theta1(ip-1);
BridgeResponse.theta1dot(ip) = BridgeResponse.theta1dot(ip-1);

BridgeResponse.X(:,ip)     = BridgeResponse.X(:,ip-1);
BridgeResponse.Xdot(:,ip)  = BridgeResponse.Xdot(:,ip-1); 
BridgeResponse.Xddot(:,ip) = BridgeResponse.Xddot(:,ip-1);

BridgeResponse.Peff(:,ip) = BridgeResponse.Peff(:,ip-1);
BridgeResponse.Pr0(:,ip)  = BridgeResponse.Pr0(:,ip);

BridgeResponse.X_Track(:,ip) = BridgeResponse.X_Track(:,ip-1);
BridgeResponse.V_Track(:,ip) = BridgeResponse.V_Track(:,ip-1);
