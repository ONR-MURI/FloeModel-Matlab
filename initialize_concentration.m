function [Floe] = initialize_concentration(c,c2_boundary,ocean,SUBFLOES,PERIODIC,height,NumFloes)
%% 
SHIFT = true;
PERIODIC = false;
[Ny, Nx] = size(c);
N = floor(NumFloes/sum(sum(c)));
c = flipud(c);
x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
xf = [-1 -1 1 1 -1];
yf = [-1 1 1 -1 -1];
Floe = [];
for jj = 1:Ny
    for ii = 1:Nx
        polynew=polyshape(mean([x(ii),(ii+1)]) +xf,mean([y(jj),y(jj+1)]) + yf);
        Floe0 = initialize_floe_values(polynew, height,SHIFT,SUBFLOES);
        Floe0.alive = 0;
        boundary = [x(ii) x(ii) x(ii+1) x(ii+1) x(ii); y(jj) y(jj+1) y(jj+1) y(jj) y(jj)];
        [dFloe]= FloeGeneratorConcentration(Floe0,boundary,c(jj,ii),N,SUBFLOES,height);
        figure
        plot([dFloe.poly])
        if sum(cat(1,dFloe.area))/area(polyshape(boundary(1,:),boundary(2,:)))<c(jj,ii)
            [dFloe,~] = pack_ice(dFloe,boundary,1.1,0,c(jj,ii),ocean,height,SUBFLOES, PERIODIC);
        end
        figure
        plot([dFloe.poly])
        Floe = [Floe dFloe];
    end
end

end

