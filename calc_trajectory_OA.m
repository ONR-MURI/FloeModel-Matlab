function [floe, fracture,Fx,Fy,Torque] =calc_trajectory_OA(dt,ocean,winds,floe,doInt)

dt = dt*doInt.step;

Xo=ocean.Xo;
Yo=ocean.Yo;
Uocn=ocean.Uocn;
Vocn=ocean.Vocn;
%[Xocn, Yocn]=meshgrid(Xo,Yo);
dXo=Xo(2)-Xo(1);
fracture = 0;

Xi=floe.Xi;
Yi=floe.Yi;
Fx = 0; Fy = 0;

% Spin up of an ice floe by an ocean eddy;
% ice floe is a disk with a radius R_floe;
% ocean eddy has a uniform vorticity and is
% centered at the center of the ice floe (hence no drift for now).

% ice-water drag is parameterized as tau=rho0*Cd*|Uocean-Uice|*(Uocean-Uice)
% no turning angle is used for now

% ice-ocean parameters

rho_ice=920; % kg/m3
rho0=1027;   % ocean density
Cd=3e-3;

% ice-water drag coefficient
rho_air=1.2;
Cd_atm=1e-3;

fc=ocean.fCoriolis; %coriolis parameter

%% ice floe params

floe_area=floe.area;
floe_mass=floe.mass; % total mass
floe_inertia_moment=floe.inertia_moment; % moment of inertia
R_floe=sqrt(2)*floe.rmax;%sqrt(max(floe.Xg)^2+max(floe.Yg)^2);

%%
U10=winds(1); % atmospheric winds
V10=winds(2); % constant here

if isnan(floe.Xi), disp('Ice floe sacked: out of ocean grid bounds!'); xx = 1; xx(1) = [1 2]; floe=[];
else
    
    x = floe.X;
    y = floe.Y;
    A_rot=[cos(floe.alpha_i) -sin(floe.alpha_i); sin(floe.alpha_i) cos(floe.alpha_i)]; %rotation matrix
    xr = A_rot*[x';y'];
    floe_mask = floe.A;%inpolygon(xr(1,:),xr(2,:),floe.c_alpha(1,:),floe.c_alpha(2,:));
    
    if sum(floe_mask(:))==0
        floe.alive = 0;
    end
    
    if  ( max(max(floe.c_alpha(1,:)))+floe.Xi>max(Xo) || min(min(floe.c_alpha(1,:)))+floe.Xi<min(Xo) || max(max(floe.c_alpha(2,:)))+floe.Yi>max(Yo) || min(min(floe.c_alpha(2,:)))+floe.Yi<min(Yo)   )       
        disp('Ice floe sacked: out of ocean grid bounds!'); floe=[];        
    elseif floe.alive == 1
        
        Xg = xr(1,:)+Xi; Yg = xr(2,:)+Yi;
        
        [theta,rho] = cart2pol(xr(1,:),xr(2,:));
        
        
        Uice=floe.Ui-rho*floe.ksi_ice.*sin(theta); % X-dir floe velocity (variable within the ice floe)
        Vice=floe.Vi+rho*floe.ksi_ice.*cos(theta); % Y-dir velocity
        
        if max(abs(Vice(:)))<100 && max(abs(Uice(:)))<100
            
            % interpolating ocean currents onto ice floe grid.
            x_ind=logical((Xo <= Xi+R_floe+2*dXo).*(Xo >= Xi-R_floe-2*dXo));
            y_ind=logical((Yo <= Yi+R_floe+2*dXo).*(Yo >= Yi-R_floe-2*dXo));
            
            Uocn_interp=interp2(Xo(x_ind),Yo(y_ind), Uocn(y_ind,x_ind),Xg,Yg);
            Vocn_interp=interp2(Xo(x_ind),Yo(y_ind), Vocn(y_ind,x_ind),Xg,Yg);
            
            Fx_atm=rho_air*Cd_atm*sqrt(U10^2+V10^2)*U10;
            Fy_atm=rho_air*Cd_atm*sqrt(U10^2+V10^2)*V10;
            
            Fx_pressureGrad=-(floe_mass/floe_area)*fc*Vocn_interp; % SSH tilt term
            Fy_pressureGrad=+(floe_mass/floe_area)*fc*Uocn_interp;        
        
            du=Uocn_interp-Uice; dv=Vocn_interp-Vice;        
        
            tau_ocnX=rho0*Cd*sqrt(du.^2+dv.^2).*( cos(ocean.turn_angle)*du+sin(ocean.turn_angle)*dv); % ocean stress with the turning angle
            tau_ocnY=rho0*Cd*sqrt(du.^2+dv.^2).*(-sin(ocean.turn_angle)*du+cos(ocean.turn_angle)*dv);
            
            Fx=tau_ocnX+Fx_atm+Fx_pressureGrad; 
            Fy=tau_ocnY+Fy_atm+Fy_pressureGrad;
            
            if max(isinf(Fx(:))) || max(isnan(Fx(:)))
                xx = 1;
                xx(1) = [1 2];
            end
            
            % updating the ice floe vorticity with averaged torques over the ice floe area
            torque=(-Fx.*sin(theta)+Fy.*cos(theta)).*rho;  % torque
            
            
            %Using 2nd order time-stepping here, utilizing tendencies calculated at
            %the previos time steps d = 1.5*dt*(d/dt)-0.5*dt*(d/dt)_previous
            
            % updating the ice floe coordinates with velocities
            dx = 1.5*dt*floe.Ui -0.5*dt*floe.dXi_p; dy = 1.5*dt*floe.Vi -0.5*dt*floe.dYi_p;
            floe.Xi=floe.Xi+dx;  floe.dXi_p=floe.Ui;
            floe.Yi=floe.Yi+dy;  floe.dYi_p=floe.Vi;
            floe.alpha_i=floe.alpha_i+1.5*dt*floe.ksi_ice-0.5*dt*floe.dalpha_i_p; floe.dalpha_i_p=floe.ksi_ice;
            
            
            % updating the ice floe velocities with mean forces and torques
            dUi_dt=(mean(Fx(floe_mask))*floe_area)/floe_mass;
            floe.Ui=floe.Ui+1.5*dt*dUi_dt-0.5*dt*floe.dUi_p;  
            if abs(floe.Ui) > 5
                xx = 1;
                xx(1) = [1 2];
            end
            floe.dUi_p=dUi_dt;
            
            dVi_dt=(mean(Fy(floe_mask))*floe_area)/floe_mass;
            floe.Vi=floe.Vi+1.5*dt*dVi_dt - 0.5*dt*floe.dVi_p;  
            if abs(floe.Vi) > 5
                xx = 1;
                xx(1) = [1 2];
            end
            floe.dVi_p=dVi_dt;
            
            dksi_ice_dt=(mean(torque(floe_mask))*floe_area)/floe_inertia_moment;
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
        else
            fracture = 1; % if large stress then designate floe to fracture releasing energy
        end
    end
end

if isnan(floe.Ui) || isnan(floe.ksi_ice) || isnan(floe.Vi) || isinf(floe.ksi_ice)
    xx =1 ;
    xx(1) = [1 2];
end

Fx = mean(Fx(floe_mask))*floe_area/floe_mass;
Fy = mean(Fy(floe_mask))*floe_area/floe_mass;
Torque = mean(torque(floe_mask))*floe_area;
end


