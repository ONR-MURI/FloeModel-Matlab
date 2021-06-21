function [floenew] = FuseFloes(floe1,floe2)
%This function takes two input floes and fuses them together while
%conserving mass and momentum
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
rho_ice = 920;
mass = cat(1,floe2.mass);

%Create and initialize the new shape
polynew = union([floe1.poly floe2.poly]);
polyout = sortregions(polynew,'area','descend');
R = regions(polyout);
poly1new = R(1);
poly1new = rmholes(poly1new);
poly = simplify(poly1new);
floenew.mass = floe1.mass + sum(mass);
height.mean = floenew.mass/(area(poly)*rho_ice);
height.delta = 0;
floenew = initialize_floe_values(poly, height);


%% Calculate new properties to conserve the momemtum of the existing floes

floenew.Ui = (floe1.Ui*floe1.mass + sum(cat(1,floe2.Ui).*mass))/(floenew.mass);
floenew.Vi = (floe1.Vi*floe1.mass + sum(cat(1,floe2.Vi).*mass))/(floenew.mass);
floenew.dUi_p = (floe1.dUi_p*floe1.mass + sum(cat(1,floe2.dUi_p).*mass))/(floenew.mass);
floenew.dVi_p = (floe1.dVi_p*floe1.mass + sum(cat(1,floe2.dVi_p).*mass))/(floenew.mass);
floenew.ksi_ice = (floe1.ksi_ice*floe1.inertia_moment + sum(cat(1,floe2.ksi_ice).*cat(1,floe2.inertia_moment)))/(floenew.inertia_moment);%use inertia moment instead of mass
floenew.dXi_p = (floe1.dXi_p*floe1.mass + sum(cat(1,floe2.dXi_p).*mass))/(floenew.mass);
floenew.dYi_p = (floe1.dYi_p*floe1.mass + sum(cat(1,floe2.dYi_p).*mass))/(floenew.mass);
floenew.dksi_ice_p = (floe1.dksi_ice_p*floe1.inertia_moment + sum(cat(1,floe2.dksi_ice_p).*cat(1,floe2.inertia_moment)))/(floenew.inertia_moment);
floenew.Fx = (floe1.Fx*floe1.mass + sum(cat(1,floe2.Fx).*mass))/(floenew.mass);
floenew.Fy = (floe1.Fy*floe1.mass + sum(cat(1,floe2.Fy).*mass))/(floenew.mass);

floenew.Stress = floe1.Stress*floe1.mass;
for ii = 1:length(floe2)
    floenew.Stress = floenew.Stress+floe2(ii).Stress*floe2(ii).mass;
end
floenew.Stress = floenew.Stress/floenew.mass;

if isinf(floenew.ksi_ice) || isnan(floenew.ksi_ice)
    xx = 1;
    xx(1) = [1 2];
end

for ii = 1:length(floenew)
    if abs(floenew(ii).area/area(polyshape(floenew(ii).c_alpha'))-1)>1e-3
        xx = 1;
        xx(1) =[1 2];
    end
end

% h = [floe1.h floe2.h];
% if floenew.h/max(h)-1 > 0.1
%     xx = 1;
%     xx(1) =[1 2];
% end

end

