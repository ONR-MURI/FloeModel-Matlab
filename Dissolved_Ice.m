function Vd_new = Dissolved_Ice(Vd,coarseMean,im_num,dissolved_new,c2_boundary,dt)
%% find current and old terms
u = squeeze(coarseMean(2,:,:,im_num));
v = squeeze(coarseMean(3,:,:,im_num));
Vdcurrent = Vd(:,:,1);
[Ny,Nx] = size(u);
if im_num==1
    uold = u;
    vold = v;
    Vdold = Vdcurrent;
else
    uold = squeeze(coarseMean(2,:,:,im_num-1));
    vold = squeeze(coarseMean(3,:,:,im_num-1));
    Vdold = Vd(:,:,2);
end

%define strength of diffusion
diffusion = 0.5;

%% Define grid
x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));

%% Find differentiation matrix for x
delx = abs(x(2)-x(1));
dx = zeros(Nx,Nx);
a = 1/(2*delx);
dx(Nx+1:1+Nx:Nx*Nx) = a;
dx(2:1+Nx:Nx*Nx-Nx) = -a;
dx(1,2) = 2*a;
dx(1,1) = -2*a;
dx(Nx,Nx-1) = -2*a;
dx(Nx,Nx) = 2*a;
d2x = dx^2;

%Find differentiation matrix for y
dely = abs(y(2)-y(1));
dy = zeros(Ny,Ny);
b = 1/(2*dely);
dy(Ny+1:1+Ny:Ny*Ny) = b;
dy(2:1+Ny:Ny*Ny-Ny) = -b;
dy(1,2) = 2*b;
dy(1,1) = -2*b;
dy(Ny,Ny-1) = -2*b;
dy(Ny,Ny) = 2*b;
d2y = dy^2;

%% Use Kroenicker product to allow for differentiation in both x and y direction
Ix = eye(Nx);
Iy = eye(Ny);
Dy = kron(Ix,dy);
Dx = kron(dx,Iy);
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

%% Time step new equation for dissolved ice. Adams Basheforth for NL terms and Crank-Nicolsen for linear terms
RHS = Vdcurrent(:)+dt*(dissolved_new(:) + diffusion/2*(kron(Ix,d2y)+kron(d2x,Iy))*Vdcurrent(:)...
    -1/2*(3*Dx*(u(:).*Vdcurrent(:))-Dx*(uold(:).*Vdold(:))) - 1/2*(3*Dy*(v(:).*Vdcurrent(:))-Dy*(vold(:).*Vdold(:))));

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
end