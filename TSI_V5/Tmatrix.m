function [T] = Tmatrix(phi)
%Tmatrix creates the transformation matrix for 2 coordinate systems rotated
%in phi radians
%   phi is the angle between the 2 coordinate systems
%   T is the transformation matrix
T = [cos(phi) sin(phi); -sin(phi) cos(phi)];
end