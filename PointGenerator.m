function [XX,YY] = PointGenerator(N,c2_boundary,SD,N0)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%N0 = fix(N^(2/3));
Lx = max(c2_boundary(1,:));
Ly = max(c2_boundary(2,:));
X = -Lx + 2*Lx*rand(1,N0);
Y = -Ly + 2*Ly*rand(1,N0);
XY = randn(fix(N/N0),fix(N/N0));
XX = []; YY = [];
for ii = 1:N0
    k = randi(length(XY));
    XX = [XX SD*XY(k,:)+X(ii)];
    YY = [YY; SD*XY(:,k)+Y(ii)];
end
% [in] = inpolygon(XX,YY,c2_boundary(1,:),c2_boundary(2,:));
% XX(~in) = [];
% YY(~in) = [];
% close all
% plot(XX,YY','ok','linewidth',2)
% hold on
% plot(X,Y','^r','linewidth',2)
% xlim([-Lx Lx])
% ylim([-Ly Ly])
% figure
% voronoi(XX,YY')
% xlim([-Lx Lx])
% ylim([-Ly Ly])
end

