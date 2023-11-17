function [RyWheel,RyRail] = ContactCurvature(...
    delta,Rc_tr,Rw_tr,Rw_ws,Rr_tr,phi_w,phi_tr,WheelGeom_cur,RailGeom_cur,side,RailProps)
% Definition of input parameters:
% - Rc_tr  - Location of point of contact in 
% - Rw_tr  - Wheelset location w/r to track coordinates
% - Rw_ws  - Location of wheel (left or right) in wheelset coordinates
% - Rr_tr  - Location of rail (left or right) in track coordinates
% - phi_w  - Rotation of wheelset in global coordinates
% - phi_tr - Rotation of track in global coordinates
% - WheelGeom_cur - Values of curvature in local coordinates
% - side   - 'l' for left side, 'r' for right side.

% Point of contact in local rail coordinates
%disp(side)
Rc_r = Rc_tr - Rr_tr';

% Point of contact in local wheel coordinates
Rc_w = Tmatrix(phi_w-phi_tr)*(Rc_tr - Rw_tr') - Rw_ws';

% Evaluate the curvatures
switch side
    case 'l'
        CurvWheel = abs(ppval(WheelGeom_cur,Rc_w(1)));
        CurvRail  = abs(ppval(RailGeom_cur,Rc_r(1)));
    case 'r'
        CurvWheel = abs(ppval(WheelGeom_cur,-Rc_w(1)));
        CurvRail  = abs(ppval(RailGeom_cur,Rc_r(1)));
end

if CurvWheel ~= 0
    RyWheel = 1/CurvWheel;
else
    RyWheel = 1000000; %Arbitrarily large number (plane)
end

if CurvRail ~= 0
    RyRail = 1/CurvRail;
else
    RyRail = 1000000; %Arbitrarily large number (plane)
end

if delta == 0
    RyWheel = RailProps.R2y;
    RyRail  = RailProps.R1y;
end

%RyWheel = RyWheel;
%RyRail  = RyRail ;

end