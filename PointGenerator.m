function [XX,YY] = PointGenerator(N,c2_boundary,SD,N0)
%This function creates the initial distribution of points throughout the
%floe domain

%Randomly distribute the points clusters will be generated around
Lx = max(c2_boundary(1,:))-min(c2_boundary(1,:));
Ly = max(c2_boundary(2,:))-min(c2_boundary(2,:));
X = min(c2_boundary(1,:)) + Lx*rand(1,N0);
Y = min(c2_boundary(2,:)) +Ly*rand(1,N0);
XY = randn(fix(N/N0),fix(N/N0));
XX = []; YY = [];

%Create points around those clusters based upon a specified standard
%deviation
for ii = 1:N0
    k = randi(length(XY));
    XX = [XX SD*XY(k,:)+X(ii)];
    YY = [YY; SD*XY(:,k)+Y(ii)];
end

end