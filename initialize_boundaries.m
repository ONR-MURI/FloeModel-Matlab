function c2_boundary =initialize_boundaries()

%Boundaries are defined as a set of contours e.g. [x1 NaN x2 ; y1 NaN y2]

%Adding walls around the domain
L_boundary=65e3;
x=[-1 -1 1 1 -1]*L_boundary; 
y=[-1 1 1 -1 -1]*L_boundary;
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