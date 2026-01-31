function [X,Y,Z] = llh2xyz(lat,long, h, a, e) 
% Convert lat, long, height data to ECEF X,Y,Z.
%
% [X,Y,Z] = llh2xyz(lat,long, h, a, e)
%
%   Inputs:
%       lat         - latitude, degrees
%       lon         - longitude, degrees
%       h           - altitude, meters
%       a           - semi-major axis for ellipsoid (default is 1980 Geodetic
%                   Reference System ellipsoid) 
%       e           - eccentricity of ellipsoid (default is 1980 Geodetic Reference
%                   System ellipsoid)
%
%   Outputs:
%       X           - X in ECEF system, meters
%       Y           - Y in ECEF system, meters
%       Z           - Z in ECEF system, meters

%--------------------------------------------------------------------------
%                             MJTL_SDR v1.0
% Copyright (C) 2026 J. Sun, Z. Tang, J. Wei and J. Zhou.
% Developed for GPS L1/L5 multi-frequency joint tracking by J. Sun, 
% Z. Tang, J. Wei and J. Zhou.
% -------------------------------------------------------------------------


if nargin <= 3
    ellipsoid = almanac('earth', 'ellipsoid', 'm');
    a = ellipsoid(1);
    e = ellipsoid(2);
end

lat = lat/180*pi;   % converting to radians
long = long/180*pi; % converting to radians

e2 = e^2;

chi = sqrt(1-e2*(sin(lat)).^2); 
X = (a./chi +h).*cos(lat).*cos(long); 
Y = (a./chi +h).*cos(lat).*sin(long); 
Z = (a*(1-e2)./chi + h).*sin(lat);