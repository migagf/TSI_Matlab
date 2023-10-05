% This subroutine calculates XTrack and VTrack

BridgeResponse.X_Track(1,ib) = BridgeResponse.X(1,ib-1) - 1.5*BridgeResponse.X(3,ib-1)*sin(pi/10);
BridgeResponse.X_Track(2,ib) = -1.5*BridgeResponse.X(3,ib-1)*cos(pi/10);
    
BridgeResponse.V_Track(1,ib) = BridgeResponse.Xdot(1,ib-1) - 1.5*BridgeResponse.Xdot(3,ib-1)*sin(pi/10);
BridgeResponse.V_Track(2,ib) = -1.5*BridgeResponse.Xdot(3,ib-1)*cos(pi/10);