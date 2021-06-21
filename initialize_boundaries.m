function c2_boundary =initialize_boundaries()

%Boundaries are defined as a set of contours e.g. [x1 NaN x2 ; y1 NaN y2]

%Adding walls around the domain
% Lx=65e4; Ly=65e4;
% x=[-1 -1 1 1 -1]*Lx; 
% y=[-1 1 1 -1 -1]*Ly;

%%Inertial
% R = 65e4;
% t = 0:pi/50:2*pi;
% x = R*cos(t); y = R*sin(t);

%%Nares
Lx=2e5; Ly=1e6;
x=[-1 -1 1 1 -1]*Lx/2; 
y=[-1 1 1 -1 -1]*Ly/2;

c2_boundary = [x; y];

% include something intresting like WA state boundaries
include_WashingtonState=0; % set to 0 if you don't want this extra boundary

if include_WashingtonState
    load('WA_LatLon.mat');
    y=(Lat-nanmean(Lat))/max(abs(Lat-nanmean(Lat)))*63e3;
    x=(Lon-nanmean(Lon))/max(abs(Lon-nanmean(Lon)))*63e3;
    x=x(1:end-140); y=y(1:end-140);
    c2_boundary = [ c2_boundary(1,:) NaN x; c2_boundary(2,:) NaN y];    
end


end