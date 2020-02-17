
for ii = 1:1e5
Dissolved_new = 0;
Vxshift = zeros(Ny-4,Nx-3);
Vxshift(:,2:end-1) = 0.5*(Vdcurrent(3:Ny-2,3:Nx-3)+Vdcurrent(3:Ny-2,4:Nx-2));
Ushift = zeros(Ny-4,Nx-3);
Ushift(:,2:end-1) = 0.5*(u(3:Ny-2,3:Nx-3)+u(3:Ny-2,4:Nx-2));
Vyshift = zeros(Ny-3,Nx-4);
Vyshift(2:end-1,:) = 0.5*(Vdcurrent(3:Ny-3,3:Nx-2)+Vdcurrent(4:Ny-2,3:Nx-2));
Vshift = zeros(Ny-3,Nx-4);
Vshift(2:end-1,:) = 0.5*(v(3:Ny-3,3:Nx-2)+v(4:Ny-2,3:Nx-2));
Ax = diff(Ushift.*Vxshift,1,2)/delx;
Ay = diff(Vshift.*Vyshift,1,1)/dely;
    RHS = Vdcurrent(:)+dt*(diffusion*(kron(Ix,d2y)+kron(d2x,Iy))*Vdcurrent(:));%...
%         -1/2*(3*(heaviside(-u(:)).*Dxback*(u(:).*Vdcurrent(:))+heaviside(u(:)).*Dxfor*(u(:).*Vdcurrent(:)))...
%         -(heaviside(-uold(:)).*Dxback*(uold(:).*Vdold(:))+heaviside(uold(:)).*Dxfor*(uold(:).*Vdold(:))))...
%         -1/2*(3*(heaviside(v(:)).*Dyback*(v(:).*Vdcurrent(:))+heaviside(-v(:)).*Dyfor*(v(:).*Vdcurrent(:)))...
%         -(heaviside(vold(:)).*Dyback*(vold(:).*Vdold(:))+heaviside(-vold(:)).*Dyfor*(vold(:).*Vdold(:)))));
%         -(heaviside(-u(:)).*(Dxback*(u(:).*Vdcurrent(:)))+heaviside(u(:)).*(Dxfor*(u(:).*Vdcurrent(:))))...
%         -(heaviside(-v(:)).*(Dyback*(v(:).*Vdcurrent(:)))+heaviside(v(:)).*(Dyfor*(v(:).*Vdcurrent(:)))));
        %-1/2*(3*Dx*(u(:).*Vdcurrent(:))-Dx*(uold(:).*Vdold(:))) - 1/2*(3*Dy*(v(:).*Vdcurrent(:))-Dy*(vold(:).*Vdold(:))));     
    %Setting RHS to 0 for BCs:: 
%     RHS(1:Ny:Nx*Ny-1) = 0;
%     RHS(Ny:Ny:Nx*Ny) = 0;
%     RHS(1:Ny) = 0;
%     RHS(Nx*Ny-Ny:Nx*Ny) = 0;
    
    Vd_new = reshape(RHS,Ny,Nx);
    Vd_new(3:Ny-2,3:Nx-2) = Vd_new(3:Ny-2,3:Nx-2)+dt*(Ay+Ax);
%     Vd_new(1,1) = mean([Vd_new(1,2),Vd_new(2,1)]);
%     Vd_new(1,Nx) = mean([Vd_new(1,Nx-1),Vd_new(2,Nx)]);
%     Vd_new(Ny,1) = mean([Vd_new(Ny,2),Vd_new(Ny-1,1)]);
%     Vd_new(Ny,Nx) = mean([Vd_new(Ny,Nx-1),Vd_new(Ny-1,Nx)]);
    
    Viscous = diffusion*(kron(Ix,d2y)+kron(d2x,Iy))*Vdcurrent(:);
    VF= reshape(Viscous,Ny,Nx);
    
    Advec1 = -1/2*(3*(heaviside(v(:)).*Dyback*(v(:).*Vdcurrent(:))+heaviside(-v(:)).*Dyfor*(v(:).*Vdcurrent(:)))...
        -(heaviside(vold(:)).*Dyback*(vold(:).*Vdold(:))+heaviside(-vold(:)).*Dyfor*(vold(:).*Vdold(:))));
    A1= reshape(Advec1,Ny,Nx);
    Advec1 = A1(3:Ny-2,3:Nx-2);
    ATrack1(ii) = sum(Ax(:)+Ay(:));
    
    Advec2 = -(heaviside(-u(:)).*(Dxback*(u(:).*Vdcurrent(:)))+heaviside(u(:)).*(Dxfor*(u(:).*Vdcurrent(:))))...
        -(heaviside(v(:)).*(Dyback*(v(:).*Vdcurrent(:)))+heaviside(-v(:)).*(Dyfor*(v(:).*Vdcurrent(:))));
    A2= reshape(Advec2,Ny,Nx);
    Advec2 = A2(3:Ny-2,3:Nx-2);
    ATrack2(ii) = sum(Advec2(:));
    
    Advec3 = -1/2*(3*(heaviside(-u(:)).*Dxback*(u(:).*Vdcurrent(:))+heaviside(u(:)).*Dxfor*(u(:).*Vdcurrent(:)))...
        -(heaviside(-uold(:)).*Dxback*(uold(:).*Vdold(:))+heaviside(uold(:)).*Dxfor*(uold(:).*Vdold(:))))...
        -1/2*(3*(heaviside(v(:)).*Dyback*(v(:).*Vdcurrent(:))+heaviside(-v(:)).*Dyfor*(v(:).*Vdcurrent(:)))...
        -(heaviside(vold(:)).*Dyback*(vold(:).*Vdold(:))+heaviside(-vold(:)).*Dyfor*(vold(:).*Vdold(:))));
    A3= reshape(Advec3,Ny,Nx);
    Advec3 = A3(3:Ny-2,3:Nx-2);
    ATrack3(ii) = sum(Advec3(:));
    
    Viscous = VF(3:Ny-2,3:Nx-2);
    VisTrack(ii) = sum(Viscous(:));
    
    Vdold = Vdcurrent;
    Vdcurrent = Vd_new;
    Vdcurrent(1:2,:) = [Vdcurrent(3,:); Vdcurrent(3,:)];
    Vdcurrent(Ny-1:Ny,:) = [Vdcurrent(Ny-2,:);Vdcurrent(Ny-2,:)];
    Vdcurrent(:,1:2) = [Vdcurrent(:,3) Vdcurrent(:,3)];
    Vdcurrent(:,Nx-1:Nx) = [Vdcurrent(:,Nx-2) Vdcurrent(:,Nx-2)];
    
    VT = Vdcurrent(3:Ny-2,3:Nx-2);
    Vtrack(ii) = sum(VT(:));
    
end