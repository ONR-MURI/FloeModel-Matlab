function [Sigma] = Calc_Stress(eularian_data,dt, c2_boundary)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% find current terms
U = squeeze(eularian_data.u);
V = squeeze(eularian_data.v);
dU = squeeze(eularian_data.du);
dV = squeeze(eularian_data.dv);
Fx = squeeze(eularian_data.force_x);
Fy = squeeze(eularian_data.force_y);
SigXX = squeeze(eularian_data.stressxx);
SigYX = squeeze(eularian_data.stressyx);
SigXY = squeeze(eularian_data.stressxy);
SigYY = squeeze(eularian_data.stressyy);

%Add ghost cells for current time step
[Ny,Nx] = size(Fx);
Nx = Nx+2;
Ny = Ny+2;
u(2:Ny-1,2:Nx-1) = U;
u(1,:) = u(2,:);
u(Ny,:) = u(Ny-1,:);
u(:,1) = u(:,2);
u(:,Nx) = u(:,Nx-1);
v(2:Ny-1,2:Nx-1) = V;
v(1,:) = v(2,:);
v(Ny,:) = v(Ny-1,:);
v(:,1) = v(:,2);
v(:,Nx) = v(:,Nx-1);
fx(2:Ny-1,2:Nx-1) = Fx;
fx(1,:) = fx(2,:);
fx(Ny,:) = fx(Ny-1,:);
fx(:,1) = fx(:,2);
fx(:,Nx) = fx(:,Nx-1);
fy(2:Ny-1,2:Nx-1) = Fy;
fy(1,:) = fy(2,:);
fy(Ny,:) = fy(Ny-1,:);
fy(:,1) = fy(:,2);
fy(:,Nx) = fy(:,Nx-1);
Sigxx(2:Ny-1,2:Nx-1) = SigXX;
Sigxx(1,:) = Sigxx(2,:);
Sigxx(Ny,:) = Sigxx(Ny-1,:);
Sigxx(:,1) = Sigxx(:,2);
Sigxx(:,Nx) = Sigxx(:,Nx-1);
Sigyx(2:Ny-1,2:Nx-1) = SigYX;
Sigyx(1,:) = Sigyx(2,:);
Sigyx(Ny,:) = Sigyx(Ny-1,:);
Sigyx(:,1) = Sigyx(:,2);
Sigyx(:,Nx) = Sigyx(:,Nx-1);
Sigxy(2:Ny-1,2:Nx-1) = SigXY;
Sigxy(1,:) = Sigxy(2,:);
Sigxy(Ny,:) = Sigxy(Ny-1,:);
Sigxy(:,1) = Sigxy(:,2);
Sigxy(:,Nx) = Sigxy(:,Nx-1);
Sigyy(2:Ny-1,2:Nx-1) = SigYY;
Sigyy(1,:) = Sigyy(2,:);
Sigyy(Ny,:) = Sigyy(Ny-1,:);
Sigyy(:,1) = Sigyy(:,2);
Sigyy(:,Nx) = Sigyy(:,Nx-1);

%% Define grid
x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/(Nx-2):max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/(Ny-2):max(c2_boundary(2,:));
delx = abs(x(2)-x(1));
dely = -abs(y(2)-y(1));

%% Shift points for calculating Stress Terms
Uxshift = u(2:Ny-1,2:Nx);
Uxshift(:,2:end-1) = 0.5*(u(2:Ny-1,2:Nx-2)+u(2:Ny-1,3:Nx-1));
Uyshift = u(2:Ny,2:Nx-1);
Uyshift(2:Ny-2,:) = 0.5*(u(2:Ny-2,2:Nx-1)+u(3:Ny-1,2:Nx-1));

Vxshift = v(2:Ny-1,2:Nx);
Vxshift(:,2:end-1) = 0.5*(v(2:Ny-1,2:Nx-2)+v(2:Ny-1,3:Nx-1));
Vyshift = v(2:Ny,2:Nx-1);
Vyshift(2:Ny-2,:) = 0.5*(v(2:Ny-2,2:Nx-1)+v(3:Ny-1,2:Nx-1));


SigXXxshift = Sigxx(2:Ny-1,2:Nx);
SigXXxshift(:,2:end-1) = 0.5*(Sigxx(2:Ny-1,2:Nx-2)+Sigxx(2:Ny-1,3:Nx-1));
SigYXyshift = Sigyx(2:Ny,2:Nx-1);
SigYXyshift(2:Ny-2,:) = 0.5*(Sigyx(2:Ny-2,2:Nx-1)+Sigyx(3:Ny-1,2:Nx-1));

SigXYxshift = Sigxy(2:Ny-1,2:Nx);
SigXYxshift(:,2:end-1) = 0.5*(Sigxy(2:Ny-1,2:Nx-2)+Sigxy(2:Ny-1,3:Nx-1));
SigYYyshift = Sigyy(2:Ny,2:Nx-1);
SigYYyshift(2:Ny-2,:) = 0.5*(Sigyy(2:Ny-2,2:Nx-1)+Sigyy(3:Ny-1,2:Nx-1));

FXxshift = fx(2:Ny-1,2:Nx);
FXxshift(:,2:end-1) = 0.5*(fx(2:Ny-1,2:Nx-2)+fx(2:Ny-1,3:Nx-1));
FXyshift = fx(2:Ny,2:Nx-1);
FXyshift(2:Ny-2,:) = 0.5*(fx(2:Ny-2,2:Nx-1)+fx(3:Ny-1,2:Nx-1));

FYxshift = fy(2:Ny-1,2:Nx);
FYxshift(:,2:end-1) = 0.5*(fy(2:Ny-1,2:Nx-2)+fy(2:Ny-1,3:Nx-1));
FYyshift = fy(2:Ny,2:Nx-1);
FYyshift(2:Ny-2,:) = 0.5*(fy(2:Ny-2,2:Nx-1)+fy(3:Ny-1,2:Nx-1));

%% Find Stress Tensor
for jj = 1:Nx-2
    for ii = 1:Ny-2
        Sigma.XX(ii,jj) = (FXxshift(ii,jj)-FXxshift(ii,jj+1))/dely;
        Sigma.XY(ii,jj) = (FYxshift(ii,jj)-FYxshift(ii,jj+1))/dely;
        Sigma.YY(ii,jj) = (FYyshift(ii,jj)-FYyshift(ii+1,jj))/delx;
        Sigma.YX(ii,jj) = (FXyshift(ii,jj)-FXyshift(ii+1,jj))/delx;
    end
end
rho_ice = 920;
DivSig1 = rho_ice*dU-Fx*rho_ice;
DivSig2 = rho_ice*dV-Fy*rho_ice;

% DivSig1 = rho_ice*dU+rho_ice*(u(2:Ny-1,2:Nx-1).*diff(Uxshift,1,2)/delx + v(2:Ny-1,2:Nx-1).*diff(Uyshift,1,1)/dely)-Fx*rho_ice;
% DivSig2 = rho_ice*dV+rho_ice*(u(2:Ny-1,2:Nx-1).*diff(Vxshift,1,2)/delx + v(2:Ny-1,2:Nx-1).*diff(Vyshift,1,1)/dely)-Fy*rho_ice;

DivSigX = diff(SigXXxshift,1,2)/delx + diff(SigYXyshift,1,1)/dely;
DivSigY = diff(SigXYxshift,1,2)/delx + diff(SigYYyshift,1,1)/dely;

xx = 1;
xx(1) = [1 2];
end

