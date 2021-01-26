function [floe, fracture,FxOA,FyOA] =calc_trajectory_Nares(dt,ocean,winds,floe,HFo, doInt, c2_boundary)


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
elseif floe.mass<0
    floe.mass = 1e3;
    floe.alive=0;
end
while max((abs(ext_force))) > floe.mass/(5*dt)
    ext_force = ext_force/10;
    ext_torque = ext_torque/10;
    if ~isempty(floe.interactions); a = a/10; end
end
% end

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
h = floe.h;
floe_inertia_moment=floe.inertia_moment; % moment of inertia
dh = HFo*dt./h;
floe_mass = (h+dh)./h.*floe_mass; floe.mass = floe_mass;
floe_inertia_moment = (h+dh)./h.*floe_inertia_moment;
floe.inertia_moment = floe_inertia_moment;
R_floe=sqrt(2)*floe.rmax;%sqrt(max(floe.Xg)^2+max(floe.Yg)^2);

%%
U10=winds(1); % atmospheric winds
V10=winds(2); % constant here

if isnan(floe.Xi), disp('Ice floe sacked: out of ocean grid bounds!'); xx = 1; xx(1) = [1 2]; floe=[];
else
    
%     A_alpha=imrotate(floe.A,-floe.alpha_i/pi*180,'bilinear','crop');
%     floe_mask=(A_alpha==1);
    if doInt.flag || isempty(floe.FxOA) || floe.h < 0.1
        x = floe.X;
        y = floe.Y;
        A_rot=[cos(floe.alpha_i) -sin(floe.alpha_i); sin(floe.alpha_i) cos(floe.alpha_i)]; %rotation matrix
        xr = A_rot*[x';y'];
        %     A = sum(in)/length(x)*4*R_floe^2;
        
        if sum(floe.A(:))==0
            floe.alive = 0;
        end
    end
%     dX=250; % resolution of the grid inside the floe
%     n=(fix(floe.rmax/dX)+1); n=dX*(-n:n);
%     [Xg, Yg]= meshgrid(n, n);
% 
%     [in] = inpolygon(Xg(:), Yg(:),floe.c0(1,:),floe.c0(2,:));
%     floe_mask=reshape(in,length(Xg),length(Xg));
    
    if  ( max(max(floe.c_alpha(1,:)))+floe.Xi>max(Xo) || min(min(floe.c_alpha(1,:)))+floe.Xi<min(Xo) || max(max(floe.c_alpha(2,:)))+floe.Yi>max(Yo) || min(min(floe.c_alpha(2,:)))+floe.Yi<min(Yo)   )       
        disp('Ice floe sacked: out of ocean grid bounds!'); floe=[];        
    elseif floe.alive == 1
        
        
%         Xg=floe.X+floe.Xi; % grid centered around the ice floe
%         Yg=floe.Y+floe.Yi;
%         
%         [theta,rho] = cart2pol(Xg-Xi,Yg-Yi);
        
        if doInt.flag || isempty(floe.FxOA) || floe.h < 0.1
            Xg = xr(1,:)+Xi; Yg = xr(2,:)+Yi;
            
            [theta,rho] = cart2pol(xr(1,:),xr(2,:));
            
            
            Uice=floe.Ui-rho*floe.ksi_ice.*sin(theta); % X-dir floe velocity (variable within the ice floe)
            Vice=floe.Vi+rho*floe.ksi_ice.*cos(theta); % Y-dir velocity
        
%         if max(abs(Vice(:)))<100 && max(abs(Uice(:)))<100
            
            % interpolating ocean currents onto ice floe grid.
            %   floe_mask_ocn_grid=logical( (Xocn <= Xi+R_floe+3*dXo).*(Xocn >= Xi-R_floe-3*dXo).*(Yocn <= Yi+R_floe+3*dXo).*(Yocn >= Yi-R_floe-3*dXo) );
            x_ind=logical((Xo <= Xi+R_floe+2*dXo).*(Xo >= Xi-R_floe-2*dXo));
            y_ind=logical((Yo <= Yi+R_floe+2*dXo).*(Yo >= Yi-R_floe-2*dXo));
            
            Uocn_interp=interp2(Xo(x_ind),Yo(y_ind), Uocn(y_ind,x_ind),Xg,Yg);
            Vocn_interp=interp2(Xo(x_ind),Yo(y_ind), Vocn(y_ind,x_ind),Xg,Yg);
%             for ii = 1:length(Xg)
%                 fU(ii) = real(sum(ocean.Uon(:).*exp(1i*(Xg(ii)*ocean.k(:)+Yg(ii)*ocean.l(:)))));
%                 fV(ii) = real(sum(ocean.Von(:).*exp(1i*(Xg(ii)*ocean.k(:)+Yg(ii)*ocean.l(:)))));
%             end
%             xx = 1; xx(1) = [1 2];
            
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
                save('floefail.mat','dt','ocean','winds','floe','HFo', 'doInt', 'c2_boundary')
                xx = 1;
                xx(1) = [1 2];
            end
            
            % updating the ice floe vorticity with averaged torques over the ice floe area
            torque=(-Fx.*sin(theta)+Fy.*cos(theta)).*rho;  % torque
            
            floe.FxOA = mean(Fx(floe.A));
            floe.FyOA = mean(Fy(floe.A));
            floe.torqueOA = mean(torque(floe.A));
        end
            
        %choosing a time step based on advection and vorticity criteria
        %   dt_max=min([ abs(0.05*max(1e-5,abs(floe.ksi_ice))/(mean(torque(floe_mask))*floe_area/floe_inertia_moment))   0.35*dXo/sqrt(floe.Ui^2+floe.Vi^2)]);
        %dt=max([100 dt]);
        
        %   dt=300; % fixed dt for all floes in order to synchronize!
        %   if (dt>dt_max); display('warning: timestep may be too large'); end
        
        
        %Using 2nd order time-stepping here, utilizing tendencies calculated at
        %the previos time steps d = 1.5*dt*(d/dt)-0.5*dt*(d/dt)_previous
        
        % updating the ice floe coordinates with velocities
        dx = 1.5*dt*floe.Ui -0.5*dt*floe.dXi_p; dy = 1.5*dt*floe.Vi -0.5*dt*floe.dYi_p;
        floe.Xi=floe.Xi+dx;  floe.dXi_p=floe.Ui;
        floe.Yi=floe.Yi+dy;  floe.dYi_p=floe.Vi;
        floe.alpha_i=floe.alpha_i+1.5*dt*floe.ksi_ice-0.5*dt*floe.dalpha_i_p; floe.dalpha_i_p=floe.ksi_ice;
            
        % check if Floe has gone out bottom boundary
        if min(floe.c_alpha(2,:))+floe.Yi<min(c2_boundary(2,:))
            floe.alive = 0;
        end
        
        % updating the ice floe velocities with mean forces and torques
        dUi_dt=(floe.FxOA*floe_area+ext_force(1))/floe_mass;
        frac = [];
        if abs(dt*dUi_dt) > 0.5
            dUi_dt = sign(dUi_dt)*0.5/dt;
            frac = dUi_dt/(floe.FxOA*floe_area+ext_force(1))/floe_mass;
        end
        floe.Ui=floe.Ui+1.5*dt*dUi_dt-0.5*dt*floe.dUi_p;  
        if abs(floe.Ui) > 5
            xx = 1;
            xx(1) = [1 2];
        end
        floe.dUi_p=dUi_dt;
            
        dVi_dt=(floe.FyOA*floe_area+ext_force(2))/floe_mass;
        if abs(dt*dVi_dt) > 0.5
            dVi_dt = sign(dVi_dt)*0.5/dt;
            frac = dVi_dt/(floe.FyOA*floe_area+ext_force(2))/floe_mass;
        end
        floe.Vi=floe.Vi+1.5*dt*dVi_dt - 0.5*dt*floe.dVi_p;  
        if abs(floe.Vi) > 5
            xx = 1;
            xx(1) = [1 2];
        end
        floe.dVi_p=dVi_dt;
        
        dksi_ice_dt=(floe.torqueOA*floe_area+ext_torque)/floe_inertia_moment;
        if ~isempty(frac)
            dksi_ice_dt = frac*dksi_ice_dt;
        end
        ksi_ice=floe.ksi_ice+1.5*dt*dksi_ice_dt - 0.5*dt*floe.dksi_ice_p;
        if abs(ksi_ice) > 1e-5
            ksi_ice = floe.ksi_ice;
        end
        floe.ksi_ice = ksi_ice;
        floe.dksi_ice_p=dksi_ice_dt;
        
        A_rot=[cos(floe.alpha_i) -sin(floe.alpha_i); sin(floe.alpha_i) cos(floe.alpha_i)]; %rotation matrix
        floe.c_alpha=A_rot*floe.c0; %rotate floe contour
            
        floe.strain = 1/4*([dUi_dt*dt*dy/floe_area dVi_dt*dt*dy/floe_area; dUi_dt*dt*dx/floe_area dVi_dt*dt*dx/floe_area] + [dUi_dt*dt*dy/floe_area dUi_dt*dt*dx/floe_area; dVi_dt*dt*dy/floe_area dVi_dt*dt*dx/floe_area]);
        %floe.strain = 1/2*([dUi_dt*dt/dx dVi_dt*dt/dx; dUi_dt*dt/dy dVi_dt*dt/dy] + [dUi_dt*dt/dx dUi_dt*dt/dy; dVi_dt*dt/dx dVi_dt*dt/dy]);
    end
end

if isnan(floe.Ui) || isnan(floe.ksi_ice) || isnan(floe.Vi) || isinf(floe.ksi_ice)
    xx =1 ;
    xx(1) = [1 2];
end

% if floe.alive == 0 && floe.area > 1e8
%     xx = 1;
%     xx(1) =[1 2];
% end
FxOA = floe.FxOA*floe_area/floe_mass;
FyOA = floe.FyOA*floe_area/floe_mass;
end


