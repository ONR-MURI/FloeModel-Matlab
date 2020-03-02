Nx = 40;

X = 0:2*pi/Nx:2*pi;
X = (X(1:end-1)+X(2:end))/2;
Y = fliplr(X)';

dt = 1e-6;

%% find current terms
U = cos(Y)*sin(X);
V = sin(Y)*cos(X);
vdcurrent = cos(Y)*cos(X);

%Add ghost cells for current time step
[Ny,Nx] = size(U);
Nx = Nx+4;
Ny = Ny+4;
u(3:Ny-2,3:Nx-2) = U;
u(1:2,:) = [u(3,:); u(3,:)];
u(Ny-1:Ny,:) = [u(Ny-2,:);u(Ny-2,:)];
u(:,1:2) = [u(:,3) u(:,3)];
u(:,Nx-1:Nx) = [u(:,Nx-2) u(:,Nx-2)];
v(3:Ny-2,3:Nx-2) = V;
v(1:2,:) = [v(3,:); v(3,:)];
v(Ny-1:Ny,:) = [v(Ny-2,:);v(Ny-2,:)];
v(:,1:2) = [v(:,3) v(:,3)];
v(:,Nx-1:Nx) = [v(:,Nx-2) v(:,Nx-2)];
Vdcurrent(3:Ny-2,3:Nx-2) = vdcurrent;
Vdcurrent(1:2,:) = [Vdcurrent(3,:); Vdcurrent(3,:)];
Vdcurrent(Ny-1:Ny,:) = [Vdcurrent(Ny-2,:);Vdcurrent(Ny-2,:)];
Vdcurrent(:,1:2) = [Vdcurrent(:,3) Vdcurrent(:,3)];
Vdcurrent(:,Nx-1:Nx) = [Vdcurrent(:,Nx-2) Vdcurrent(:,Nx-2)];


%define strength of diffusion
diffusion = 1e4;

Dissolved_new = zeros(Ny,Nx);
Dissolved_new(3:Ny-2,3:Nx-2) = dt*(2*diffusion*cos(Y)*cos(X)+cos(Y).^2*cos(2*X)+cos(2*Y)*(cos(X).^2));

%% Define grid
x = X;
y = Y;
delx = abs(x(2)-x(1));
dely = abs(y(2)-y(1));
x = [min(x)-2*delx min(x)-delx x max(x)+delx max(x)+2*delx];
y = [max(y)+2*dely; max(y)+dely; y; min(y)-dely; min(y)-2*dely];

%% Find differentiation matrix for x
d2x = zeros(Nx,Nx);
b = 1/(delx^2);
d2x(2:Nx+1:Nx*Nx) = b;
d2x(2+Nx:Nx+1:Nx*Nx) = -2*b;
d2x(2+2*Nx:Nx+1:Nx*Nx) = b;
d2x(1,1) = b/(3/2);
d2x(1,2) = -2*b/(3/2);
d2x(1,3) = b/(3/2);
d2x(Nx,Nx) = b/(3/2);
d2x(Nx,Nx-1) = -2*b/(3/2);
d2x(Nx,Nx-2) = b/(3/2);

%Find differentiation matrix for y
d2y = zeros(Ny,Ny);
b = 1/(dely^2);
d2y(2:Ny+1:Ny*Ny) = b;
d2y(2+Ny:Ny+1:Ny*Ny) = -2*b;
d2y(2+2*Ny:Ny+1:Ny*Ny) = b;
d2y(1,1) = b/(3/2);
d2y(1,2) = -2*b/(3/2);
d2y(1,3) = b/(3/2);
d2y(Ny,Ny) = b/(3/2);
d2y(Ny,Ny-1) = -2*b/(3/2);
d2y(Ny,Ny-2) = b/(3/2);

Ix = eye(Nx);
Iy = eye(Ny);

%% Shift points for calculating Advective Terms
Vxshift = zeros(Ny-4,Nx-3);
Vxshift(:,2:end-1) = 0.5*(Vdcurrent(3:Ny-2,3:Nx-3)+Vdcurrent(3:Ny-2,4:Nx-2));
Ushift = zeros(Ny-4,Nx-3);
Ushift(:,2:end-1) = 0.5*(u(3:Ny-2,3:Nx-3)+u(3:Ny-2,4:Nx-2));
Vyshift = zeros(Ny-3,Nx-4);
Vyshift(2:end-1,:) = 0.5*(Vdcurrent(3:Ny-3,3:Nx-2)+Vdcurrent(4:Ny-2,3:Nx-2));
Vshift = zeros(Ny-3,Nx-4);
Vshift(2:end-1,:) = 0.5*(v(3:Ny-3,3:Nx-2)+v(4:Ny-2,3:Nx-2));
%% Time step new equation for dissolved ice. Adams Basheforth for NL terms and Crank-Nicolsen for linear terms
%Calculate RHS first without advection terms
for ii = 1:1e5
    RHS = Vdcurrent(:)+Dissolved_new(:)+dt*diffusion*(kron(Ix,d2y)+kron(d2x,Iy))*Vdcurrent(:);
    Vd_new = reshape(RHS,Ny,Nx);
    
    %calculate advection terms
    Advec = diff(Ushift.*Vxshift,1,2)/delx+ diff(Vshift.*Vyshift,1,1)/(-dely);
    Vd_new(3:Ny-2,3:Nx-2) = Vd_new(3:Ny-2,3:Nx-2)-dt*Advec;
    
    err = (Vdcurrent(3:Ny-2,3:Nx-2)-Vd_new(3:Ny-2,3:Nx-2))./Vdcurrent(3:Ny-2,3:Nx-2);
    max(max(err))
    pause
    
    Vdcurrent(3:Ny-2,3:Nx-2) = Vd_new(3:Ny-2,3:Nx-2);
    Vdcurrent(1:2,:) = [Vdcurrent(3,:); Vdcurrent(3,:)];
    Vdcurrent(Ny-1:Ny,:) = [Vdcurrent(Ny-2,:);Vdcurrent(Ny-2,:)];
    Vdcurrent(:,1:2) = [Vdcurrent(:,3) Vdcurrent(:,3)];
    Vdcurrent(:,Nx-1:Nx) = [Vdcurrent(:,Nx-2) Vdcurrent(:,Nx-2)];
end
%end
Vd_new = Vd_new(3:Ny-2,3:Nx-2);