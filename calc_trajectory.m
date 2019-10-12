function floe=calc_trajectory(dt,ocean,winds,floe)

ext_force=floe.collision_force;
ext_torque=floe.collision_torque;

Xo=ocean.Xo;
Yo=ocean.Yo;
Uocn=ocean.Uocn;
Vocn=ocean.Vocn;
%[Xocn, Yocn]=meshgrid(Xo,Yo);
dXo=Xo(2)-Xo(1);


Xi=floe.Xi;
Yi=floe.Yi;


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

%% ice floe params

floe_area=floe.area;
floe_mass=floe.mass; % total mass
floe_inertia_moment=floe.inertia_moment; % moment of inertia
R_floe=floe.rmax;

%%
U10=winds(1); % atmospheric winds
V10=winds(2); % constant here

if isnan(floe.Xi), disp('Ice floe sacked: out of ocean grid bounds!'); floe=[];
else
    
    A_alpha=imrotate(floe.A,-floe.alpha_i/pi*180,'bilinear','crop');
    floe_mask=(A_alpha==1);
    
    if  ( max(max(floe.X(floe_mask)))+floe.Xi>max(Xo) || min(min(floe.X(floe_mask)))+floe.Xi<min(Xo) || max(max(floe.Y(floe_mask)))+floe.Yi>max(Yo) || min(min(floe.Y(floe_mask)))+floe.Yi<min(Yo)   )       
        disp('Ice floe sacked: out of ocean grid bounds!'); floe=[];        
    else
        
        
        Xg=floe.X+floe.Xi; % grid centered around the ice floe
        Yg=floe.Y+floe.Yi;
        
        [theta,rho] = cart2pol(Xg-Xi,Yg-Yi);
        
        
        Uice=floe.Ui-rho*floe.ksi_ice.*sin(theta); % X-dir floe velocity (variable within the ice floe)
        Vice=floe.Vi+rho*floe.ksi_ice.*cos(theta); % Y-dir velocity
        
        % interpolating ocean currents onto ice floe grid.
        %   floe_mask_ocn_grid=logical( (Xocn <= Xi+R_floe+3*dXo).*(Xocn >= Xi-R_floe-3*dXo).*(Yocn <= Yi+R_floe+3*dXo).*(Yocn >= Yi-R_floe-3*dXo) );
        x_ind=logical((Xo <= Xi+R_floe+dXo).*(Xo >= Xi-R_floe-dXo));
        y_ind=logical((Yo <= Yi+R_floe+dXo).*(Yo >= Yi-R_floe-dXo));
        
        Uocn_interp=interp2(Xo(x_ind),Yo(y_ind), Uocn(y_ind,x_ind),Xg,Yg);
        Vocn_interp=interp2(Xo(x_ind),Yo(y_ind), Vocn(y_ind,x_ind),Xg,Yg);
        
        Fx_atm=rho_air*Cd_atm*sqrt(U10^2+V10^2)*U10;
        Fy_atm=rho_air*Cd_atm*sqrt(U10^2+V10^2)*V10;
        
        Fx=rho0*Cd*sqrt((Uocn_interp-Uice).^2+(Vocn_interp-Vice).^2).*(Uocn_interp-Uice)+Fx_atm; % ice-water drag forces in X and Y-dir
        Fy=rho0*Cd*sqrt((Uocn_interp-Uice).^2+(Vocn_interp-Vice).^2).*(Vocn_interp-Vice)+Fy_atm;
        
        % updating the ice floe vorticity with averaged torques over the ice floe area
        torque=(-Fx.*sin(theta)+Fy.*cos(theta)).*rho;  % torque
        
        %choosing a time step based on advection and vorticity criteria
        %   dt_max=min([ abs(0.05*max(1e-5,abs(floe.ksi_ice))/(mean(torque(floe_mask))*floe_area/floe_inertia_moment))   0.35*dXo/sqrt(floe.Ui^2+floe.Vi^2)]);
        %dt=max([100 dt]);
        
        %   dt=300; % fixed dt for all floes in order to synchronize!
        %   if (dt>dt_max); display('warning: timestep may be too large'); end
        
        
        %Using 2nd order time-stepping here, utilizing tendencies calculated at
        %the previos time steps d = 1.5*dt*(d/dt)-0.5*dt*(d/dt)_previous
        
        % updating the ice floe coordinates with velocities
        floe.Xi=floe.Xi+1.5*dt*floe.Ui -0.5*dt*floe.dXi_p;  floe.dXi_p=floe.Ui;
        floe.Yi=floe.Yi+1.5*dt*floe.Vi -0.5*dt*floe.dYi_p;  floe.dYi_p=floe.Vi;
        floe.alpha_i=floe.alpha_i+1.5*dt*floe.ksi_ice-0.5*dt*floe.dalpha_i_p; floe.dalpha_i_p=floe.ksi_ice;
        
        % updating the ice floe velocities with mean forces and torques
        dUi_dt=(mean(Fx(floe_mask))*floe_area+ext_force(1))/floe_mass;
        floe.Ui=floe.Ui+1.5*dt*dUi_dt-0.5*dt*floe.dUi_p;  floe.dUi_p=dUi_dt;
        
        dVi_dt=(mean(Fy(floe_mask))*floe_area+ext_force(2))/floe_mass;
        floe.Vi=floe.Vi+1.5*dt*dVi_dt - 0.5*dt*floe.dVi_p;  floe.dVi_p=dVi_dt;
        
        dksi_ice_dt=(mean(torque(floe_mask))*floe_area+ext_torque)/floe_inertia_moment;
        floe.ksi_ice=floe.ksi_ice+1.5*dt*dksi_ice_dt - 0.5*dt*floe.dksi_ice_p; floe.dksi_ice_p=dksi_ice_dt;
        
        A_rot=[cos(floe.alpha_i) -sin(floe.alpha_i); sin(floe.alpha_i) cos(floe.alpha_i)]; %rotation matrix
        floe.c_alpha=A_rot*floe.c0; %rotate floe contour
    end
end


end


