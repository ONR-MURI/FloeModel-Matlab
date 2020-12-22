function [Floe] = FracLeads(Floe,Ny,Nx,Nb,c2_boundary,eularian_data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Create coarse grids
x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
Lx = (x(2)-x(1))/2; Ly = (y(2)-y(1))/2;

%Find floes and create bins 
Xi=cat(1,Floe.Xi);
Yi=cat(1,Floe.Yi);
Ri=cat(1,Floe.rmax);
Binx = fix((Xi-min(x))/(max(x)-min(x))*Nx+1);
Biny = fix((Yi-min(y))/(max(y)-min(y))*Ny+1);

% Idenfity floes that are alive
live = cat(1,Floe.alive);
SigXX = flipud(squeeze(eularian_data.stressxx)); SigYX = flipud(squeeze(eularian_data.stressyx));
SigXY = flipud(squeeze(eularian_data.stressxy)); SigYY = flipud(squeeze(eularian_data.stressyy));

%Place all floes in correct bins
fracturedFloes = [];
% xx = 1; xx(1) =[1 2];
for ii = 1:Nx
    for jj = 1:Ny
        if abs(SigXX(jj,ii))>0
            alive = live(live == 1 & Binx == ii & Biny == jj);
            floenew = Floe(live == 1 & Binx == ii & Biny == jj);
            [~,V] = eig([SigXX(jj,ii) SigYX(jj,ii);SigXY(jj,ii) SigYY(jj,ii)]);
            r = randi([1 2],1,1);
            m = vecnorm(V(1,r)/V(2,r));
            if isinf(m)
                m = Ly;
            end
            xi = (x(ii)+x(ii+1))/2; yi = (y(jj)+y(jj+1))/2;
            x0 = Xi(live == 1 & Binx == ii & Biny == jj)-xi;
            y0 = Yi(live == 1 & Binx == ii & Biny == jj)-yi;
            Lx = max(abs(x0)+Ri(live == 1 & Binx == ii & Biny == jj));
            Ly = max(abs(y0)+Ri(live == 1 & Binx == ii & Biny == jj));
            X = Lx*rand(1);Y = Ly*rand(1);
            y1 = Y+m*(Lx-X); y2 = Y+m*(-Lx-X);
            x1 = X+(Ly-Y)/m; x2 = X+(-Ly-Y)/m;
            if m == 0
                x1 = 1e16; x2 = 1e16;
            end
            xx = [x1 x2 Lx -Lx]; yy = [Ly -Ly y1 y2];
            %bound = [x(ii) x(ii) x(ii+1) x(ii+1) x(ii);y(jj) y(jj+1) y(jj+1) y(jj) y(jj)];
            bound = [-Lx -Lx Lx Lx -Lx; -Ly Ly Ly -Ly -Ly]+[xi;yi];
            [in] = inpolygon(xx+xi,yy+yi,bound(1,:)',bound(2,:)');
            var = min([Lx-abs(X), Ly-abs(Y)]);
            mid = var*rand(2,1);
            Xn = X+mid(1)+xi; Yn = Y+mid(2)+yi;
            xc = xx(in)+xi; yc = yy(in)+yi;
            Pc = cutpolygon(bound', [xc(1) yc(1); xc(2) yc(2)],1);
            Pc(1,:) = [];
            while Lx+Ly-sum(abs(Pc(1,:)-[xi,yi]))>100
                Pc=[Pc(2:end,:);Pc(1,:)];
            end
            [k,~] = dsearchn(Pc,[xc(1) yc(1); xc(2) yc(2)]);
            Pc1 = polyshape([Pc(1:min(k),:); Xn Yn; Pc(min(k)+1:end,:)]);
            c2 = [xc(1) Xn xc(2); yc(1) Yn yc(2)]; 
            for kk = 1:length(alive)
                c1=[floenew(kk).c_alpha(1,:)+floenew(kk).Xi; floenew(kk).c_alpha(2,:)+floenew(kk).Yi];
                P=InterX(c1,c2);
                if ~isempty(P)
                    Poly = polyshape(floenew(kk).c_alpha'+[floenew(kk).Xi floenew(kk).Yi]);
                    polynew = subtract(Poly,Pc1);
%                 if area(polynew)/floenew(kk).area < 1
                    polynew2 = subtract(Poly,polynew);
                    poly = [polynew polynew2];
                    Floes=fractures(floenew(kk),poly);
                    fracturedFloes = [fracturedFloes Floes];
                    alive(kk) = 0;
                end
            end
%             for kk =1:length(floenew)
%                 if alive(kk)
%                     poly(kk) = polyshape(floenew(kk).c_alpha'+[floenew(kk).Xi floenew(kk).Yi]);
%                 end
%             end
%             for kk =1:length(fracturedFloes)
%                 polyF(kk) = polyshape(fracturedFloes(kk).c_alpha'+[fracturedFloes(kk).Xi fracturedFloes(kk).Yi]);
%             end
            live(live == 1 & Binx == ii & Biny == jj) = alive;
        end
    end
end
Floe = [Floe(logical(live)) fracturedFloes];
for kk =1:length(Floe)
    polyF(kk) = polyshape(Floe(kk).c_alpha'+[Floe(kk).Xi Floe(kk).Yi]);
end

end

