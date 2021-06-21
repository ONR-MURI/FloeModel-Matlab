function [Floe] = FracMohr(Floe,Nb,min_floe_size)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
        q = 5.2; SigC = 250e3;
        A = cat(1,Floe.area);
        Sig1 = (1/q+1)*SigC/(1/q-q);
        Sig2 = q*Sig1+SigC;
        Sig11 = 1e8;
        Sig22 = q*Sig11+SigC;
        MohrX = [Sig1; Sig11; Sig22];
        MohrY = [Sig2; Sig22; Sig11];
        Mohr = polyshape(MohrX,MohrY);
        
        for ii = 1:length(Floe)
            Stress = eig(Floe(ii).Stress);
            Princ1(ii) = max(Stress);
            Princ2(ii) = min(Stress);
        end
        %[d_min, ~, ~] = p_poly_dist(Princ1, Princ2, Mohr.Vertices(:,1), Mohr.Vertices(:,2),true);
        in = inpolygon(Princ1,Princ2,MohrX, MohrY);
        %keep=rand(length(Floe),1)>d_min;
        keep = zeros(length(Floe),1);
        keep(in) = 1;
        keep(A<min_floe_size)=1;
        keep(1:Nb) = ones(Nb,1);
        keep = logical(keep);
        fracturedFloes=fracture_floe(Floe(~keep),3);
        if ~isempty(fracturedFloes)
            Floe=[Floe(keep) fracturedFloes];
        end

end

