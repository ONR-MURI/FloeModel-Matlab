function [floe,Fx,Fy] = calc_trajectory_int(dt,ocean,floe,HFo)

ext_force=floe.collision_force;
ext_torque=floe.collision_torque;
HFo = mean(HFo(:));
if ~isempty(floe.interactions)
    a=floe.interactions;
    r=[floe.Xi floe.Yi];
    floe.Stress =1/(2*floe.area*floe.h)*([sum((a(:,4)-r(1)).*a(:,2)) sum((a(:,5)-r(2)).*a(:,2)); sum((a(:,4)-r(1)).*a(:,3)) sum((a(:,5)-r(2)).*a(:,3))]...
        +[sum(a(:,2).*(a(:,4)-r(2))) sum(a(:,3).*(a(:,4)-r(1))); sum(a(:,2).*(a(:,5)-r(2))) sum(a(:,3).*(a(:,5)-r(2)))]);
end

if length(ext_force) == 1
    ext_force = [0 0];
end

if floe.h > 10
    xx = 1;
    xx(1) =[1 2];
end
while max((abs(ext_force))) > floe.mass/(5*dt)
    ext_force = ext_force/10;
    ext_torque = ext_torque/10;
    if ~isempty(floe.interactions); a = a/10; end
end
% end

Xo=ocean.Xo;
Yo=ocean.Yo;
%[Xocn, Yocn]=meshgrid(Xo,Yo);

% Spin up of an ice floe by an ocean eddy;
% ice floe is a disk with a radius R_floe;
% ocean eddy has a uniform vorticity and is
% centered at the center of the ice floe (hence no drift for now).

% ice-water drag is parameterized as tau=rho0*Cd*|Uocean-Uice|*(Uocean-Uice)
% no turning angle is used for now


%% ice floe params

floe_mass=floe.mass; % total mass
h = floe.h;
floe_inertia_moment=floe.inertia_moment; % moment of inertia
dh = HFo*dt./h;
floe_mass = (h+dh)./h.*floe_mass; floe.mass = floe_mass;
floe_inertia_moment = (h+dh)./h.*floe_inertia_moment;
floe.inertia_moment = floe_inertia_moment;

%%

if isnan(floe.Xi), disp('Ice floe sacked: out of ocean grid bounds!'); xx = 1; xx(1) = [1 2]; floe=[];
else
    
    if  ( max(max(floe.c_alpha(1,:)))+floe.Xi>max(Xo) || min(min(floe.c_alpha(1,:)))+floe.Xi<min(Xo) || max(max(floe.c_alpha(2,:)))+floe.Yi>max(Yo) || min(min(floe.c_alpha(2,:)))+floe.Yi<min(Yo)   )
        disp('Ice floe sacked: out of ocean grid bounds!'); floe=[];
    elseif floe.alive == 1
        
        
        %Using 2nd order time-stepping here, utilizing tendencies calculated at
        %the previos time steps d = 1.5*dt*(d/dt)-0.5*dt*(d/dt)_previous
        
        % updating the ice floe coordinates with velocities
        dx = 1.5*dt*floe.Ui -0.5*dt*floe.dXi_p; dy = 1.5*dt*floe.Vi -0.5*dt*floe.dYi_p;
        floe.Xi=floe.Xi+dx;  floe.dXi_p=floe.Ui;
        floe.Yi=floe.Yi+dy;  floe.dYi_p=floe.Vi;
        floe.alpha_i=floe.alpha_i+1.5*dt*floe.ksi_ice-0.5*dt*floe.dalpha_i_p; floe.dalpha_i_p=floe.ksi_ice;
        
        
        % updating the ice floe velocities with mean forces and torques
        dUi_dt=ext_force(1)/floe_mass;
        frac = [];
        if abs(dt*dUi_dt) > 0.5
            dUi_dt = sign(dUi_dt)*0.5/dt;
            frac = dUi_dt/(ext_force(1))/floe_mass;
        end
        floe.Ui=floe.Ui+1.5*dt*dUi_dt-0.5*dt*floe.dUi_p;  floe.dUi_p=dUi_dt;
        if abs(floe.Ui) > 5
            xx = 1;
            xx(1) = [1 2];
        end
        
        dVi_dt=ext_force(2)/floe_mass;
        if abs(dt*dVi_dt) > 0.5
            dVi_dt = sign(dVi_dt)*0.5/dt;
            frac = dVi_dt/(ext_force(2))/floe_mass;
        end
        floe.Vi=floe.Vi+1.5*dt*dVi_dt - 0.5*dt*floe.dVi_p;  floe.dVi_p=dVi_dt;
        if abs(floe.Vi) > 5
            xx = 1;
            xx(1) = [1 2];
        end
        
        dksi_ice_dt=ext_torque/floe_inertia_moment;
        if ~isempty(frac)
            dksi_ice_dt = frac*dksi_ice_dt;
        end
        ksi_ice=floe.ksi_ice+1.5*dt*dksi_ice_dt - 0.5*dt*floe.dksi_ice_p;
        if abs(ksi_ice) > 1e-4
            xx = 1;
            xx(1) = [1 2];
        end
        floe.ksi_ice=ksi_ice;
        floe.dksi_ice_p=dksi_ice_dt;
        
        A_rot=[cos(floe.alpha_i) -sin(floe.alpha_i); sin(floe.alpha_i) cos(floe.alpha_i)]; %rotation matrix
        floe.c_alpha=A_rot*floe.c0; %rotate floe contour
        
        floe.strain = 1/2*([dUi_dt*dt/dx dVi_dt*dt/dx; dUi_dt*dt/dy dVi_dt*dt/dy] + [dUi_dt*dt/dx dUi_dt*dt/dy; dVi_dt*dt/dx dVi_dt*dt/dy]);
    end
end

if isnan(floe.Ui) || isnan(floe.ksi_ice) || isnan(floe.Vi) || isinf(floe.ksi_ice)
    xx =1 ;
    xx(1) = [1 2];
end

Fx = ext_force(1)/floe_mass;
Fy = ext_force(2)/floe_mass;
end


