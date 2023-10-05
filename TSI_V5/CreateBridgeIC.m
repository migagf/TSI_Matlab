
BridgeIC.eps        = BridgeResponse.eps(end)          ;   
BridgeIC.z          = BridgeResponse.z(end)            ;
BridgeIC.Mtheta     = BridgeResponse.Mtheta(end)       ;  
BridgeIC.ktheta     = BridgeResponse.ktheta(end)       ;  
BridgeIC.theta1     = BridgeResponse.theta1(end)       ;  
BridgeIC.theta1dot  = BridgeResponse.theta1dot(end)    ;     
BridgeIC.X          = BridgeResponse.X(:,end)          ;
BridgeIC.Xdot       = BridgeResponse.Xdot(:,end)       ;  
BridgeIC.Xddot      = BridgeResponse.Xddot(:,end)      ;   
BridgeIC.Peff       = BridgeResponse.Peff(:,end)       ;  
BridgeIC.Pr0        = BridgeResponse.Pr0(:,end)        ; 
BridgeIC.X_Track    = BridgeResponse.X_Track(:,end)    ;     
BridgeIC.V_Track    = BridgeResponse.V_Track(:,end)    ; 

save('BridgeIC','BridgeIC')