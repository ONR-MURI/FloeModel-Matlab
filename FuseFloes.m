function [floenew] = FuseFloes(floe1,floe2)
%This function takes two input floes and fuses them together while
%conserving mass and momentum
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
rho_ice = 920;

%Create and initialize the new shape
polynew = union(floe1.poly, floe2.poly);
polyout = sortregions(polynew,'area','descend');
R = regions(polyout);
poly1new = R(1);
poly1new = rmholes(poly1new);
poly = simplify(poly1new);
floenew.mass = floe1.mass + floe2.mass;
height.mean = floenew.mass/(area(poly)*rho_ice);
height.delta = 0;
floenew = initialize_floe_values(poly, height);


%% Calculate new properties to conserve the momemtum of the existing floes

floenew.Ui = (floe1.Ui*floe1.mass + floe2.Ui*floe2.mass)/(floenew.mass);
floenew.Vi = (floe1.Vi*floe1.mass + floe2.Vi*floe2.mass)/(floenew.mass);
floenew.dUi_p = (floe1.dUi_p*floe1.mass + floe2.dUi_p*floe2.mass)/(floenew.mass);
floenew.dVi_p = (floe1.dVi_p*floe1.mass + floe2.dVi_p*floe2.mass)/(floenew.mass);
floenew.ksi_ice = (floe1.ksi_ice*floe1.inertia_moment + floe2.ksi_ice*floe2.inertia_moment)/(floenew.inertia_moment);%use inertia moment instead of mass
floenew.dXi_p = (floe1.dXi_p*floe1.mass + floe2.dXi_p*floe2.mass)/(floenew.mass);
floenew.dYi_p = (floe1.dYi_p*floe1.mass + floe2.dYi_p*floe2.mass)/(floenew.mass);
floenew.dksi_ice_p = (floe1.dksi_ice_p*floe1.inertia_moment + floe2.dksi_ice_p*floe2.inertia_moment)/(floenew.inertia_moment);

if isinf(floenew.ksi_ice) || isnan(floenew.ksi_ice)
    xx = 1;
    xx(1) = [1 2];
end

end

