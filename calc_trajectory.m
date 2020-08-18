function floe=calc_trajectory(dt,ocean,winds,floe,heat_flux,c2_boundary_poly, SUBFLOES)
%%This function advances the individual floes in time calculating the new
%%trajectory based upon the calculated forces

ext_force=floe.collision_force;
ext_torque=floe.collision_torque;

Xo=ocean.Xo;
Yo=ocean.Yo;
[Xocn,Yocn] = meshgrid(Xo,Yo);
Uocn=ocean.Uocn;
Vocn=ocean.Vocn;
dXo=Xo(2)-Xo(1);

Xa=winds.X;
Ya=winds.Y;
[Xatm,Yatm] = meshgrid(Xa,Ya);
Uatm=winds.U;
Vatm=winds.V;
dXa=Xa(2)-Xa(1);

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
Cd=5.5e-3;

% ice-water drag coefficient
rho_air=1.2;
Cd_atm=1.6e-3;

fc=ocean.fCoriolis; %coriolis parameter

%% ice floe params

[ka,~] = dsearchn([Xatm(:),Yatm(:)],[Xi,Yi]);
[ko,~] = dsearchn([Xocn(:),Yocn(:)],[Xi,Yi]);
HFa = heat_flux.atm(ka);
HFo = heat_flux.oc(ko);
                        
floe_area=floe.area;
R_floe=floe.rmax;
N = length(floe.SubFloes);
areaS = zeros(N,1);
inertia = zeros(N,1);
centers = zeros(N,2);
for ii = 1:N
    floe.SubFloes(ii).h = floe.SubFloes(ii).h-HFo*dt/(floe.SubFloes(ii).h);
    areaS(ii) = area(floe.SubFloes(ii).poly);
    if floe.SubFloes(ii).poly.NumHoles > 0
        breaks = isnan(floe.SubFloes(ii).poly.Vertices(:,1));
        I = find(breaks == 1);
        I = [0 I' length(breaks)+1];
        inertia(ii) = 0;
        for jj = length(I) -1
            inertia(ii) = inertia(ii) + PolygonMoments(floe.SubFloes(ii).poly.Vertices(I(jj)+1:I(jj+1)-1,:),floe.SubFloes(ii).h);
        end
    else
        inertia(ii) = PolygonMoments(floe.SubFloes(ii).poly.Vertices,floe.SubFloes(ii).h);
    end
    [Xi,Yi] = centroid(floe.SubFloes(ii).poly);
    centers(ii,:) = [Xi,Yi];
end
floe.mass = sum(rho_ice*areaS.*cat(1,floe.SubFloes.h));
floe_mass=floe.mass; % total mass
floe.Xm = sum(rho_ice*areaS.*cat(1,floe.SubFloes.h).*centers(:,1))./floe.mass;
floe.Ym = sum(rho_ice*areaS.*cat(1,floe.SubFloes.h).*centers(:,2))./floe.mass;
floe.h = sum(rho_ice*areaS.*cat(1,floe.SubFloes.h).*cat(1,floe.SubFloes.h))./floe.mass;
floe.inertia_moment = sum(inertia+cat(1,floe.SubFloes.h).*sqrt((centers(:,1)-floe.Xm).^2+(centers(:,2)-floe.Ym).^2));
if floe.h <= 0
    floe.alive = 0;
end


%%

if isnan(floe.Xi)||isnan(floe.alpha_i)||isnan(floe.ksi_ice), disp('Ice floe sacked: NaN state vars.'); floe=[];
%     x = 1;
%     x(1) = [1 2];
    
else
    
    
    if  ( max(floe.poly.Vertices(:,1))>max(Xo) || min(floe.poly.Vertices(:,1))<min(Xo) || max(floe.poly.Vertices(:,2))>max(Yo) || min(floe.poly.Vertices(:,2))<min(Yo)   )
        disp('Ice floe sacked: out of ocean grid bounds!');
%         x = 1;
%         x(1) = [1 2];
        floe=[];
    else
        %Advance the floe in time
        
        dX=500; % resolution of the grid inside the flow
        n=(fix(floe.rmax/dX)+1); n=dX*(-n:n);
        [Xg, Yg]= meshgrid(n+floe.Xi, n+floe.Yi);
        
        in=inpoly2([Xg(:) Yg(:)],floe.poly.Vertices);
        floe_mask=reshape(in,length(Xg),length(Xg));
        
        if sum(floe_mask(:)) < 2
            floe.alive = 0; %If floe is to small mark floe to stop being tracked
        else
        
            
            [theta,rho] = cart2pol(Xg-Xi,Yg-Yi); % rho is the radius from the center of the floe; theta is the angle
            
            
            Uice=floe.Ui-rho*floe.ksi_ice.*sin(theta); % X-dir floe velocity (variable within the ice floe)
            Vice=floe.Vi+rho*floe.ksi_ice.*cos(theta); % Y-dir velocity
            
            % interpolating ocean and atm currents onto ice floe grid.
            x_ind=logical((Xo <= Xi+R_floe+2*dXo).*(Xo >= Xi-R_floe-2*dXo));
            y_ind=logical((Yo <= Yi+R_floe+2*dXo).*(Yo >= Yi-R_floe-2*dXo));
            
            Uocn_interp=interp2(Xo(x_ind),Yo(y_ind), Uocn(y_ind,x_ind),Xg,Yg);
            Vocn_interp=interp2(Xo(x_ind),Yo(y_ind), Vocn(y_ind,x_ind),Xg,Yg);
            
            Uatm_interp=interp2(Xo(x_ind),Yo(y_ind), Uatm(y_ind,x_ind),Xg,Yg);
            Vatm_interp=interp2(Xo(x_ind),Yo(y_ind), Vatm(y_ind,x_ind),Xg,Yg);
            
            % all forces are per unit area, units of N/m^2;
            Fx_atm = rho_air*Cd_atm*sqrt(Uatm_interp.^2.+Vatm_interp.^2).*Uatm_interp;  % wind drag (no turning angle here!)
            Fy_atm = rho_air*Cd_atm*sqrt(Uatm_interp.^2.+Vatm_interp.^2).*Vatm_interp;
            
            Fx_pressureGrad=-(floe_mass/floe_area)*fc*Vocn_interp; % SSH tilt term
            Fy_pressureGrad=+(floe_mass/floe_area)*fc*Uocn_interp;        
        
            du=Uocn_interp-Uice; dv=Vocn_interp-Vice;        
        
            tau_ocnX=rho0*Cd*sqrt(du.^2+dv.^2).*( cos(ocean.turn_angle)*du+sin(ocean.turn_angle)*dv); % ocean stress with the turning angle
            tau_ocnY=rho0*Cd*sqrt(du.^2+dv.^2).*(-sin(ocean.turn_angle)*du+cos(ocean.turn_angle)*dv);
        
            Fx=tau_ocnX+Fx_atm+Fx_pressureGrad; % adding up all forces except the Coriolis force
            Fy=tau_ocnY+Fy_atm+Fy_pressureGrad;
        
            % updating the ice floe vorticity with averaged torques over the ice floe area
            torque=(-Fx.*sin(theta)+Fy.*cos(theta)).*rho;  % torque
        
            %adding the remaining Coriolis force; it has no torque.
            Fx=Fx+(floe_mass/floe_area)*fc*floe.Vi;
            Fy=Fy-(floe_mass/floe_area)*fc*floe.Ui;
            
            
            
            %Using 2nd order time-stepping here, utilizing tendencies calculated at
            %the previos time steps d = 1.5*dt*(d/dt)-0.5*dt*(d/dt)_previous
            
            % updating the ice floe coordinates with velocities
            dalpha=1.5*dt*floe.ksi_ice-0.5*dt*floe.dalpha_i_p;
            floe.alpha_i=floe.alpha_i+dalpha; floe.dalpha_i_p=floe.ksi_ice;
            if isinf(floe.alpha_i)
                floe.alive = 0;
%                 if floe.area > 1e7
%                     xx = 1;
%                     xx(1) = [1 2];
%                 end
            end
            
            A_rot=[cos(dalpha) -sin(dalpha); sin(dalpha) cos(dalpha)]; %rotation matrix
            
            
            floe.poly=rotate(floe.poly,dalpha*180/pi,[floe.Xi, floe.Yi]);

            
            vorVert = (A_rot*([floe.vorX floe.vorY] - [floe.Xi floe.Yi])')';
            floe.vorX = vorVert(:,1)+floe.Xi; floe.vorY = vorVert(:,2)+floe.Yi;
            dx=1.5*dt*floe.Ui-0.5*dt*floe.dXi_p;
            floe.Xi=floe.Xi+dx; floe.Xm=floe.Xm+dx;  floe.dXi_p=floe.Ui;
            floe.vorX = floe.vorX + dx; floe.vorbox(:,1) = floe.vorbox(:,1) + dx;
            
            dy=1.5*dt*floe.Vi-0.5*dt*floe.dYi_p;
            floe.Yi=floe.Yi+dy;  floe.Ym=floe.Ym+dy;floe.dYi_p=floe.Vi;
            floe.vorY = floe.vorY + dy; floe.vorbox(:,2) = floe.vorbox(:,2) + dy;
            
            floe.poly = translate(floe.poly,[dx,dy]);
            polyout = sortregions(floe.poly,'area','descend');
            R = regions(polyout);
            floe.poly = R(1);
            
            if SUBFLOES
                for ii = 1:length(floe.SubFloes)
                    floe.SubFloes(ii).poly=rotate(floe.SubFloes(ii).poly,dalpha*180/pi,[floe.Xi, floe.Yi]);
                    floe.SubFloes(ii).poly = translate(floe.SubFloes(ii).poly,[dx dy]);
                    floe.SubFloes(ii).h = floe.SubFloes(ii).h + dh;
                    if floe.SubFloes(ii).h > 30; floe.SubFloes(ii).h = 30; end;
                end
            else
                floe.SubFloes = [];
                floe.SubFloes.poly = floe.poly;
                floe.SubFloes.h = floe.h;
            end
            
            % updating the ice floe velocities with mean forces and torques
            dUi_dt=(mean(Fx(floe_mask))*floe_area+ext_force(1))/floe_mass;
            du=1.5*dt*dUi_dt-0.5*dt*floe.dUi_p;
            if abs(du) > 1
                du = sign(du);
                dUi_dt=(du+0.5*dt*floe.dUi_p)/(1.5*dt);
%             elseif isnan(du)
%                 xx = 1;
%                 xx(1) = [1 2];
            end
            floe.Ui=floe.Ui+du;  floe.dUi_p=dUi_dt;
            if abs(floe.Ui) > 5 
                floe.alive = 0;
%                 if floe.h > 0.5 && floe.area > 1e7
%                     xx = 1;
%                     xx(1) = [1 2];
%                 end
            end
            
            dVi_dt=(mean(Fy(floe_mask))*floe_area+ext_force(2))/floe_mass;
            dv=1.5*dt*dVi_dt - 0.5*dt*floe.dVi_p;
            if abs(dv) > 1
                dv = sign(dv);
                dVi_dt=(dv+0.5*dt*floe.dVi_p)/(1.5*dt);
%             elseif isnan(du)
%                 xx = 1;
%                 xx(1) = [1 2];
            end
            floe.Vi=floe.Vi+dv;  floe.dVi_p=dVi_dt;
            if  abs(floe.Vi) > 5 && floe.area > 1e7
                floe.alive = 0;
%                 if floe.h > 0.5
%                     xx = 1;
%                     xx(1) = [1 2];
%                 end
            end
                       
            dksi_ice_dt=(mean(torque(floe_mask))*floe_area+ext_torque)/floe.inertia_moment;
            dksi=1.5*dt*dksi_ice_dt - 0.5*dt*floe.dksi_ice_p;
            if abs(dksi) > 1e-4
                dksi = sign(dksi)*1e-4;
                dksi_ice_dt=(dksi+0.5*dt*floe.dksi_ice_p)/(1.5*dt);
%             elseif isnan(dksi)
%                 xx = 1;
%                 xx(1) = [1 2];
            end
            floe.ksi_ice=floe.ksi_ice+dksi; floe.dksi_ice_p=dksi_ice_dt;
            if  abs(floe.ksi_ice) > 5 && floe.area > 1e7
                floe.alive = 0;
%                 if floe.h > 0.5
%                     xx = 1;
%                     xx(1) = [1 2];
%                 end
            end
            
            
            if ~isempty(floe.potentialInteractions)
                if area(intersect(floe.poly,union([floe.potentialInteractions.c])))/floe.area > 0.5
                    floe.alive=0;
                end
            end
            if isinf(floe.ksi_ice)
                floe.alive = 0;
%                 if floe.area > 1e7
%                     xx = 1;
%                     xx(1) = [1 2];
%                 end
            end
%             if isnan(floe.Xi) || isnan(floe.Yi)
%                 xx = 1;
%                 xx(1) = [1 2];
%             elseif isnan(floe.ksi_ice)
%                 xx = 1;
%                 xx(1) = [1 2];
%             end
        end
    end
    
end


end


