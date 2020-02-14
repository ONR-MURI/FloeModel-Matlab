function Vd_new = Dissolved_Ice(Vd,coarseMean,im_num,dissolved_new,c2_boundary,dt)
%% find current terms
U = squeeze(coarseMean(2,:,:,im_num));
V = squeeze(coarseMean(3,:,:,im_num));
vdcurrent = Vd(:,:,1);

%Add ghost cells for current time step
[Ny,Nx] = size(U);
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
Vdcurrent(2:Ny-1,2:Nx-1) = vdcurrent;
Vdcurrent(1,:) = Vdcurrent(2,:);
Vdcurrent(Ny,:) = Vdcurrent(Ny-1,:);
Vdcurrent(:,1) = Vdcurrent(:,2);
Vdcurrent(:,Nx) = Vdcurrent(:,Nx-1);

%Set ghost cells for old terms
if im_num==1
    uold = u;
    vold = v;
    Vdold = Vdcurrent;
else
    Uold = squeeze(coarseMean(2,:,:,im_num-1));
    Vold = squeeze(coarseMean(3,:,:,im_num-1));
    vdold = Vd(:,:,2);
end
uold(2:Ny-1,2:Nx-1) = Uold;
uold(1,:) = uold(2,:);
uold(Ny,:) = uold(Ny-1,:);
uold(:,1) = uold(:,2);
uold(:,Nx) = uold(:,Nx-1);
vold(2:Ny-1,2:Nx-1) = Vold;
vold(1,:) = vold(2,:);
vold(Ny,:) = vold(Ny-1,:);
vold(:,1) = vold(:,2);
vold(:,Nx) = vold(:,Nx-1);
Vdold(2:Ny-1,2:Nx-1) = vdold;
Vdold(1,:) = Vdold(2,:);
Vdold(Ny,:) = Vdold(Ny-1,:);
Vdold(:,1) = Vdold(:,2);
Vdold(:,Nx) = Vdold(:,Nx-1);
Dissolved_new = zeros(Ny,Nx);
Dissolved_new(2:Ny-1,2:Nx-1) = dissolved_new;

%define strength of diffusion
diffusion = 1e4;

%% Define grid
x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/(Nx-2):max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/(Ny-2):max(c2_boundary(2,:));
delx = abs(x(2)-x(1));
dely = abs(y(2)-y(1));
x = [min(x)-delx x max(x)+delx];
y = [min(y)-dely y max(y)+dely];

%% Find differentiation matrix for x
dx = zeros(Nx,Nx);
a = 1/(2*delx);
dx(Nx+1:1+Nx:Nx*Nx) = a;
dx(2:1+Nx:Nx*Nx-Nx) = -a;
dx(1,2) = 2*a;
dx(1,1) = -2*a;
dx(Nx,Nx-1) = -2*a;
dx(Nx,Nx) = 2*a;
d2x = dx^2;

dxfor = zeros(Nx,Nx);
a = 1/(delx);
dxfor(1:1+Nx:Nx*Nx) = -a;
dxfor(Nx+1:1+Nx:Nx*Nx-1) = a;
dxfor(Nx,Nx) = a;
dxfor(Nx,Nx-1) = -a;

dxback = zeros(Nx,Nx);
a = 1/(delx);
dxback(Nx+2:1+Nx:Nx*Nx) = a;
dxback(2:1+Nx:Nx*Nx-Nx) = -a;
dxback(1,2) = a;
dxback(1,1) = -a;

%Find differentiation matrix for y
dy = zeros(Ny,Ny);
b = 1/(2*dely);
dy(Ny+1:1+Ny:Ny*Ny) = b;
dy(2:1+Ny:Ny*Ny-Ny) = -b;
dy(1,2) = 2*b;
dy(1,1) = -2*b;
dy(Ny,Ny-1) = -2*b;
dy(Ny,Ny) = 2*b;
d2y = dy^2;

dyfor = zeros(Ny,Ny);
a = 1/(dely);
dyfor(1:1+Ny:Ny*Ny) = -a;
dyfor(Ny+1:1+Ny:Ny*Ny-1) = a;
dyfor(Ny,Ny) = a;
dyfor(Ny,Ny-1) = -a;

dyback = zeros(Ny,Ny);
a = 1/(dely);
dyback(Ny+2:1+Ny:Ny*Ny) = a;
dyback(2:1+Ny:Ny*Ny-Ny) = -a;
dyback(1,2) = a;
dyback(1,1) = -a;

%% Use Kroenicker product to allow for differentiation in both x and y direction
Ix = eye(Nx);
Iy = eye(Ny);
Dy = kron(Ix,dy);
Dyback = kron(Ix,dyback);
Dyfor = kron(Ix,dyfor);
Dx = kron(dx,Iy);
Dxback = kron(dxback,Iy);
Dxfor = kron(dxfor,Iy);
L = kron(Iy,Ix)-diffusion*dt/2*(kron(Ix,d2y)+kron(d2x,Iy));

%Add Boundary conditions
%Dy = 0 at top and bottom 
for ii = Ny:Ny:Nx*Ny
    L(ii,ii-1) = -2*a;
    L(ii,ii) = 2*a;
end
for ii = 0:Ny:Nx*Ny-Ny
    L(ii+1,ii+1) = -2*a;
    L(ii+1,ii+2) = 2*a;
end
%Dx = 0 at left and right
for ii = 1:Ny
    L(ii,ii) = -2*a;
    L(ii,ii+Ny) = 2*a;
    L(Nx*Ny-Ny+ii,Nx*Ny-2*Ny + ii) = -2*a;
    L(Nx*Ny-Ny+ii,Nx*Ny-Ny+ii) = 2*a;
end
%t = 0;
%% Time step new equation for dissolved ice. Adams Basheforth for NL terms and Crank-Nicolsen for linear terms
%while t < tnew
    RHS = Vdcurrent(:)+Dissolved_new(:)+dt*(diffusion/2*(kron(Ix,d2y)+kron(d2x,Iy))*Vdcurrent(:)...
        -1/2*(3*(heaviside(u(:)).*Dxback*(u(:).*Vdcurrent(:))+heaviside(-u(:)).*Dxfor*(u(:).*Vdcurrent(:)))...
        -(heaviside(uold(:)).*Dxback*(uold(:).*Vdold(:))+heaviside(-uold(:)).*Dxfor*(uold(:).*Vdold(:))))...
        -1/2*(3*(heaviside(v(:)).*Dyback*(v(:).*Vdcurrent(:))+heaviside(-v(:)).*Dyfor*(v(:).*Vdcurrent(:)))...
        -(heaviside(vold(:)).*Dyback*(vold(:).*Vdold(:))+heaviside(-vold(:)).*Dyfor*(vold(:).*Vdold(:)))));
        %-1/2*(3*Dx*(u(:).*Vdcurrent(:))-Dx*(uold(:).*Vdold(:))) - 1/2*(3*Dy*(v(:).*Vdcurrent(:))-Dy*(vold(:).*Vdold(:))));     
    %Setting RHS to 0 for BCs
    RHS(1:Ny:Nx*Ny-1) = 0;
    RHS(Ny:Ny:Nx*Ny) = 0;
    RHS(1:Ny) = 0;
    RHS(Nx*Ny-Ny:Nx*Ny) = 0;
    
    %Invert, reshape and take care of corners
    VdNEW = L\RHS;
    Vd_new = reshape(VdNEW,Ny,Nx);
    Vd_new(1,1) = mean([Vd_new(1,2),Vd_new(2,1)]);
    Vd_new(1,Nx) = mean([Vd_new(1,Nx-1),Vd_new(2,Nx)]);
    Vd_new(Ny,1) = mean([Vd_new(Ny,2),Vd_new(Ny-1,1)]);
    Vd_new(Ny,Nx) = mean([Vd_new(Ny,Nx-1),Vd_new(Ny-1,Nx)]);
%    Vdold = Vdcurrent;
%    Vdcurrent = Vd_new;
%    t = t+dt;
%end
Vd_new = Vd_new(2:Ny-1,2:Nx-1);
end