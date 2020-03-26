function floe=calc_trajectory(dt,ocean,winds,floe)

ext_force=floe.collision_force;
ext_torque=floe.collision_torque;

Xo=ocean.Xo;
Yo=ocean.Yo;
Uocn=ocean.Uocn;
Vocn=ocean.Vocn;
%[Xocn, Yocn]=meshgrid(Xo,Yo);
dXo=Xo(2)-Xo(1);

if isempty(floe.SubFloes)
    SUBFLOES = true;
    floe = initialize_floe_values(floe.poly,SUBFLOES);
end

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

if isnan(floe.Xi)||isnan(floe.alpha_i)||isnan(floe.ksi_ice), disp('Ice floe sacked: NaN state vars.'); floe=[];
else
     
%     if  ( floe.Xi+floe.rmax>max(Xo) || floe.Xi-floe.rmax<min(Xo) || floe.Yi+floe.rmax>max(Yo) || floe.Yi-floe.rmax<min(Yo)   )       
%         disp('Ice floe sacked: out of ocean grid bounds!'); floe=[];        
%     else
        
    if  ( max(floe.poly.Vertices(:,1))>max(Xo) || min(floe.poly.Vertices(:,1))<min(Xo) || max(floe.poly.Vertices(:,2))>max(Yo) || min(floe.poly.Vertices(:,2))<min(Yo)   )
%checking out of bounds only for Y-direction as X-direction is periodic.
%     if  ~PERIODIC && (max(floe.poly.Vertices(:,2))>max(Yo) || min(floe.poly.Vertices(:,2))<min(Yo)   )
%     if (max(floe.poly.Vertices(:,2))>max(Yo) || min(floe.poly.Vertices(:,2))<min(Yo)   )
        disp('Ice floe sacked: out of ocean grid bounds!'); floe=[];        
    else
        
%        A_alpha=imrotate(floe.A,-floe.alpha_i/pi*180,'bilinear','crop');

%        floe_mask=(A_alpha==1);
        dX=500; % resolution of the grid inside the flow
        n=(fix(floe.rmax/dX)+1); n=dX*(-n:n);
        [Xg, Yg]= meshgrid(n+floe.Xi, n+floe.Yi);
        
       % in=inpolygon(Xg(:),Yg(:),floe.poly.Vertices(:,1),floe.poly.Vertices(:,2)); 
        in=inpoly2([Xg(:) Yg(:)],floe.poly.Vertices); 
        floe_mask=reshape(in,length(Xg),length(Xg));

        if sum(floe_mask(:)) < 2
            floe.alive = 0;
        else
        
            [theta,rho] = cart2pol(Xg-Xi,Yg-Yi);
            
            
            Uice=floe.Ui-rho*floe.ksi_ice.*sin(theta); % X-dir floe velocity (variable within the ice floe)
            Vice=floe.Vi+rho*floe.ksi_ice.*cos(theta); % Y-dir velocity
            
            % interpolating ocean currents onto ice floe grid.
            %   floe_mask_ocn_grid=logical( (Xocn <= Xi+R_floe+3*dXo).*(Xocn >= Xi-R_floe-3*dXo).*(Yocn <= Yi+R_floe+3*dXo).*(Yocn >= Yi-R_floe-3*dXo) );
            x_ind=logical((Xo <= Xi+R_floe+2*dXo).*(Xo >= Xi-R_floe-2*dXo));
            y_ind=logical((Yo <= Yi+R_floe+2*dXo).*(Yo >= Yi-R_floe-2*dXo));
            
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
            dalpha=1.5*dt*floe.ksi_ice-0.5*dt*floe.dalpha_i_p;
            floe.alpha_i=floe.alpha_i+dalpha; floe.dalpha_i_p=floe.ksi_ice;
            
            A_rot=[cos(dalpha) -sin(dalpha); sin(dalpha) cos(dalpha)]; %rotation matrix
%             Vertices=(A_rot*(floe.poly.Vertices - [floe.Xi floe.Yi])')'; %rotate floe contour around its old center of mass
            
            floe.poly=rotate(floe.poly,dalpha*180/pi,[floe.Xi, floe.Yi]);
            for ii = 1:length(floe.SubFloes)
                floe.SubFloes(ii).poly=rotate(floe.SubFloes(ii).poly,dalpha*180/pi,[floe.Xi, floe.Yi]);
            end
            vorVert = (A_rot*([floe.vorX floe.vorY] - [floe.Xi floe.Yi])')';
            floe.vorX = vorVert(:,1)+floe.Xi; floe.vorY = vorVert(:,2)+floe.Yi;
            dx=1.5*dt*floe.Ui-0.5*dt*floe.dXi_p;
            floe.Xi=floe.Xi+dx; floe.Xm=floe.Xm+dx;  floe.dXi_p=floe.Ui;
            floe.vorX = floe.vorX + dx; floe.vorbox(:,1) = floe.vorbox(:,1) + dx;
            
            dy=1.5*dt*floe.Vi-0.5*dt*floe.dYi_p;
            floe.Yi=floe.Yi+dy;  floe.Ym=floe.Ym+dy;floe.dYi_p=floe.Vi;
            for ii = 1:length(floe.SubFloes)
                floe.SubFloes(ii).poly = translate(floe.SubFloes(ii).poly,[dx dy]);
            end
            floe.vorY = floe.vorY + dy; floe.vorbox(:,2) = floe.vorbox(:,2) + dy;
            
            floe.poly.Vertices= floe.poly.Vertices + [dx dy];  % shift by its updated center coordinate
            
            
            % updating the ice floe velocities with mean forces and torques
            dUi_dt=(mean(Fx(floe_mask))*floe_area+ext_force(1))/floe_mass;
            du=1.5*dt*dUi_dt-0.5*dt*floe.dUi_p;
            floe.Ui=floe.Ui+du;  floe.dUi_p=dUi_dt;
            
            dVi_dt=(mean(Fy(floe_mask))*floe_area+ext_force(2))/floe_mass;
            dv=1.5*dt*dVi_dt - 0.5*dt*floe.dVi_p;
            floe.Vi=floe.Vi+dv;  floe.dVi_p=dVi_dt;
            
            dksi_ice_dt=(mean(torque(floe_mask))*floe_area+ext_torque)/floe_inertia_moment;
            dksi=1.5*dt*dksi_ice_dt - 0.5*dt*floe.dksi_ice_p;
            floe.ksi_ice=floe.ksi_ice+dksi; floe.dksi_ice_p=dksi_ice_dt;
            

        end
    end
    
end


end


