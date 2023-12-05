% This subroutine calculates XTrack and VTrack

if on_bridge
    BridgeResponse.X_Track(1,ib) = BridgeResponse.Xt(1,ib-1) - 0.5 * BridgeResponse.Xt(3,ib-1);
    BridgeResponse.X_Track(2,ib) = -0.6 * BridgeResponse.Xt(3,ib-1);
    
    BridgeResponse.V_Track(1,ib) = BridgeResponse.Xtdot(1,ib-1) - 0.5 * BridgeResponse.Xtdot(3,ib-1);
    BridgeResponse.V_Track(2,ib) = -0.6 * BridgeResponse.Xdot(3,ib-1);

else
    %disp(BridgeResponse.X(1, ib - 1))
    BridgeResponse.X_Track(1,ib) = BridgeResponse.Xt(1, ib - 1);
    BridgeResponse.X_Track(2,ib) = BridgeResponse.Xt(2, ib - 1);
    
    BridgeResponse.V_Track(1,ib) = BridgeResponse.Xtdot(1, ib - 1);
    BridgeResponse.V_Track(2,ib) = BridgeResponse.Xtdot(2, ib - 1);

end


