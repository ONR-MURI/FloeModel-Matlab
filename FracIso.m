function [Floe] = FracIso(Floe,Ny,Nx,Nb,c2_boundary,Sig)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Create coarse grids
x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
Lx = (x(2)-x(1))/2; Ly = (y(2)-y(1))/2;

%Find floes and create bins 
Xi=cat(1,Floe.Xi);
Yi=cat(1,Floe.Yi);
Binx = fix((Xi-min(x))/(max(x)-min(x))*Nx+1);
Biny = fix((Yi-min(y))/(max(y)-min(y))*Ny+1);

% Idenfity floes that are alive
live = cat(1,Floe.alive);
for ii =1:length(Floe)
    Stress(ii) = max(eig(Floe(ii).Stress));
end
Stress = abs(Stress);
Stress(Stress<1000)= 0;
StressN = zeros(length(Stress),1);
Sig = flipud(abs(Sig));

%Place all floes in correct bins
for ii = 1:Nx
    for jj = 1:Ny
        if Sig(jj,ii)>0
            test = Stress(live == 1 & Binx == ii & Biny == jj)/Sig(jj,ii);
            StressN(live == 1 & Binx == ii & Biny == jj) = Stress(live == 1 & Binx == ii & Biny == jj)/(Sig(jj,ii)*max(test));
            xi = (x(ii)+x(ii+1))/2; yi = (y(jj)+y(jj+1))/2;
            X = Lx*rand(1);Y = Ly*rand(1);
            m = 1;
            y1 = Y+m*(Lx-X); y2 = Y+m*(-Lx-X);
            x1 = X+(Ly-Y)/m; x2 = X+(-Ly-Y)/m;
            xx = [x1 x2 Lx -Lx]; yy = [Ly -Ly y1 y2];
            bound = [x(ii) x(ii) x(ii+1) x(ii+1) x(ii);y(ii) y(ii+1) y(ii+1) y(ii) y(ii)];
            [in] = inpolygon(xx+xi,yy+yi,bound(1,:)',bound(2,:)');
            var = min([Lx-abs(X)+xi, Ly-abs(Y)+yi]);
            mid = var*rand(2,1);
            X = X+mid(1)+xi; Y = Y+mid(2)+yi;
            xc = xx(in)+xi; yc = yy(in)+yi;
            Pc = cutpolygon(bound', [xc(1) yc(1); xc(2) yc(2)],1);
            Pc(1,:) = [];
            while Lx+Ly-sum(abs(Pc(1,:)))>0
                Pc=[Pc(2:end,:);Pc(1,:)];
            end
            [k,dist] = dsearchn(Pc,[xc(1) yc(1); xc(2) yc(2)]);
            Pc1 = [Pc(1:min(k),:); X Y; Pc(min(k)+1:end,:)];
            xx = 1; xx(1) = [1 2];
        end
    end
end


if max(StressN>0)
    keep=8*rand(length(Floe),1)>StressN.^2;
    %keep=3*rand(length(Floe),1)>Stress'/max(Stress);
%     xx = 1; xx(1) =[1 2];
else
    overlapArea=cat(1,Floe.OverlapArea)./cat(1,Floe.area);
    keep=rand(length(Floe),1)>2*overlapArea;
end
keep(1:Nb) = ones(Nb,1);
fracturedFloes=fracture_floe(Floe(~keep),3);
if ~isempty(fracturedFloes)
    Floe=[Floe(keep) fracturedFloes];
end
end

