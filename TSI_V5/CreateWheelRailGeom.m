%% Rail and Wheel Geometries
WheelGeom = 25.4/1000*WheelGeom;

% Interpolation polynomials of the wheel and rail geometries

RailGeom_pol = pchip(RailGeom(:,1),RailGeom(:,2));
RailGeom_cur = pchip(RailGeom(:,1),RailGeom(:,3)); % Curvatures of the rail geom

WheelGeom_pol = pchip(WheelGeom(:,1),WheelGeom(:,2));
WheelGeom_cur = pchip(WheelGeom(:,1),WheelGeom(:,3)); % Curvatures of the wheel geom

% Rail geometry, heights measured from the center of the rail (m);
Rail_geom = 25.4/1000*[-1.3355 -1.3167 -1.1874 -0.8725 -0.6244 0.0000 0.6244 0.8725 1.1874 1.3167 1.3355;
                        2.3292  3.0933  3.4382  3.6292  3.6612 3.6751 3.6612 3.6292 3.4382 3.0933 2.3292]';

% Spline approximation
fgrid_n = 1000; % Number of nodes for the fine grid
Rail_geom_fine = zeros(fgrid_n,2);

%Rail_geom_fine(:,1) = linspace(min(Rail_geom(:,1)),max(Rail_geom(:,1)),fgrid_n);

nodes = linspace(0,fgrid_n-1,fgrid_n)';
a = min(Rail_geom(:,1)); b = max(Rail_geom(:,1));
Rail_geom_fine(:,1) = 0.5*(a+b)+0.5*(b-a)*cos(nodes/fgrid_n*pi);

Rail_geom_fine(:,2) = ppval(RailGeom_pol,Rail_geom_fine(:,1)); % Evaluate the polynomial in the defined nodes
Rail_geom_fine(:,2) = -Rail_geom_fine(:,2); % Z coordinates are inverted.

%% Wheel Geometry
% Wheel_Geom is a matrix that contains points with the discretization of the
% Wheel geometry in local coordinate system. Loaded as input.

fgrid_n = 2*fgrid_n;
% Left Wheel Geometry and discretization
LWheel_geom_fine = zeros(fgrid_n,2);
LWheel_geom_fine(:,1) = linspace(min(WheelGeom_pol.breaks),max(WheelGeom_pol.breaks),fgrid_n);
%LWheel_geom_fine(:,2) = -pchip(WheelGeom(:,1),WheelGeom(:,2),LWheel_geom_fine(:,1));
LWheel_geom_fine(:,2) = -ppval(WheelGeom_pol,LWheel_geom_fine(:,1));

% Right Wheel Geometry
RWheel_geom_fine = zeros(fgrid_n,2);
RWheel_geom_fine(:,1) = -LWheel_geom_fine(:,1);
RWheel_geom_fine(:,2) = LWheel_geom_fine(:,2);
