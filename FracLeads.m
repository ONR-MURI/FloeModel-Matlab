function [Floe] = FracLeads(Floe,Ny,Nx,Nb,c2_boundary,eularian_data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Create coarse grids
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)

for kk = 1:length(Floe)
    polyO(kk) = polyshape(Floe(kk).c_alpha' +[Floe(kk).Xi Floe(kk).Yi]);
end

x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
dx = (x(2)-x(1))/2; dy = (y(2)-y(1))/2;
box = [-dx -dx dx dx -dx; -dy dy dy -dy -dy];

%Find floes and create bins 
Xi=cat(1,Floe.Xi);
Yi=cat(1,Floe.Yi);
Ri=cat(1,Floe.rmax);
Binx = fix((Xi-min(x))/(max(x)-min(x))*Nx+1);
Biny = fix((Yi-min(y))/(max(y)-min(y))*Ny+1);

% Idenfity floes that are alive
live = cat(1,Floe.alive);
bound = ones(length(Floe),1); bound(1:Nb) = zeros(1:Nb);
SigXX = flipud(squeeze(eularian_data.stressxx)); SigYX = flipud(squeeze(eularian_data.stressyx));
SigXY = flipud(squeeze(eularian_data.stressxy)); SigYY = flipud(squeeze(eularian_data.stressyy));

%Place all floes in correct bins
fracturedFloes = [];
% xx = 1; xx(1) =[1 2];
for ii = 1:Nx
    for jj = 1:Ny
        if abs(SigXX(jj,ii))>0
            alive = live(live == 1 & Binx == ii & Biny == jj & bound == 1);
            floenew = Floe(live == 1 & Binx == ii & Biny == jj & bound == 1);
            shift =[(x(ii)+x(ii+1))/2,(y(jj)+y(jj+1))/2];
            dx = zeros(1,length(floenew)); dy = dx;
            for kk = 1:length(floenew)
                dx(kk) = max(floenew(kk).c_alpha(1,:)-(floenew(kk).Xi-shift(1)));
                dy(kk) = max(floenew(kk).c_alpha(2,:)-(floenew(kk).Yi-shift(2)));
            end
            dx = max(dx); dy = max(dy);
            box = [-dx -dx dx dx -dx; -dy dy dy -dy -dy];
            
            thet2 = atan(2*SigXY(jj,ii)/(SigXX(jj,ii)-SigYY(jj,ii)));
%             m = [tan(thet2/2) -tan(thet2/2)];
            r = randi([0 1],1,1);
            thet = (-1)^r*thet2;
%             [V,~] = eig([SigXX(jj,ii) SigYX(jj,ii);SigXY(jj,ii) SigYY(jj,ii)]);
%             r = randi([1 2],1,1);
%             thet = atan(V(2,r)/V(1,r));
            Rot = [cos(thet) -sin(thet); sin(thet) cos(thet)];
            P = randomwalk(21,1); xx = -10:10;
            walk = [xx;P'-P(11)]; P2 = ([dx; dx/5]*ones(1,length(xx))).*(Rot*walk)/5;
            x1 = dx*(rand(1)-0.5); y1 = dy*(rand(1)-0.5);
            P2 = [P2(1,:)+x1; P2(2,:)+y1];
            crossing=InterX(box,P2);
            [in] = inpolygon(P2(1,:)',P2(2,:)',box(1,:)',box(2,:)');
            [m,n] = size(crossing);
            if n==2
                if abs(min(crossing(1,:))+dx)<1 && abs(max(crossing(1,:))-dx)<1
                    [~,i1] = min(crossing(1,:)); [~,i2] = max(crossing(1,:));
                    lead = [crossing(:,i1) P2(:,in) crossing(:,i2) box(:,4:5)];
                    polylead = polyshape(lead');
                    if polylead.NumRegions >1
                        lead = [crossing(:,i2) P2(:,in) crossing(:,i1) box(:,4:5)];
                        polylead = polyshape(lead');
                    end
                elseif abs(min(crossing(1,:))+dx)<1 && abs(max(crossing(2,:))-dy)<1
                    [~,i1] = min(crossing(1,:)); [~,i2] = max(crossing(2,:));
                    lead = [crossing(:,i1) P2(:,in) crossing(:,i2) box(:,2)];
                    polylead = polyshape(lead');
                    if polylead.NumRegions >1
                        lead = [crossing(:,i2) P2(:,in) crossing(:,i1) box(:,2)];
                        polylead = polyshape(lead');
                    end
                elseif abs(min(crossing(1,:))+dx)<1 && abs(min(crossing(2,:))+dy)<1
                    [~,i1] = min(crossing(1,:)); [~,i2] = max(crossing(1,:));
                    lead = [crossing(:,i1) P2(:,in) crossing(:,i2) box(:,1)];
                    polylead = polyshape(lead');
                    if polylead.NumRegions >1
                        lead = [crossing(:,i2) P2(:,in) crossing(:,i1) box(:,1)];
                        polylead = polyshape(lead');
                    end
                elseif abs(min(crossing(2,:))+dy)<1 && abs(max(crossing(1,:))-dx)<1
                    [~,i1] = min(crossing(1,:)); [~,i2] = max(crossing(1,:));
                    lead = [crossing(:,i1) P2(:,in) crossing(:,i2) box(:,4)];
                    polylead = polyshape(lead');
                    if polylead.NumRegions >1
                        lead = [crossing(:,i2) P2(:,in) crossing(:,i1) box(:,4)];
                        polylead = polyshape(lead');
                    end
                elseif abs(max(crossing(2,:))-dy)<1 && abs(max(crossing(1,:))-dx)<1
                    [~,i1] = min(crossing(1,:)); [~,i2] = max(crossing(1,:));
                    lead = [crossing(:,i1) P2(:,in) crossing(:,i2) box(:,3)];
                    polylead = polyshape(lead');
                    if polylead.NumRegions >1
                        lead = [crossing(:,i2) P2(:,in) crossing(:,i1) box(:,3)];
                        polylead = polyshape(lead');
                    end
                elseif abs(min(crossing(2,:))+dy)<1 && abs(max(crossing(2,:))-dy)<1
                    [~,i1] = min(crossing(2,:)); [~,i2] = max(crossing(2,:));
                    lead = [crossing(:,i1) P2(:,in) crossing(:,i2) box(:,3:4)];
                    polylead = polyshape(lead');
                    if polylead.NumRegions >1
                        lead = [crossing(:,i2) P2(:,in) crossing(:,i1) box(:,4) box(:,3)];
                        polylead = polyshape(lead');
                    end
                end
                
                polylead = rmholes(polylead);
                if polylead.NumRegions >1
                    xx = 1; xx(1) =[1 2];
                end
                
                %            m = vecnorm(V(1,r)/V(2,r));
                %             if m>0 && m<Inf
                %                 xx = 1; xx(1) =[1 2];
                %             end
                if isinf(m)
                    m = Ly;
                end
    
                %             polylead = translate(polylead,[xi,yi]);
%                 x0 = Xi(live == 1 & Binx == ii & Biny == jj)-xi;
%                 y0 = Yi(live == 1 & Binx == ii & Biny == jj)-yi;
                for kk = 1:length(alive)
                    P=InterX(floenew(kk).c_alpha+[floenew(kk).Xi-shift(1);floenew(kk).Yi-shift(2)],P2);
                    [~,m] = size(P);
                    if m>1
                        Poly = polyshape(floenew(kk).c_alpha'+[floenew(kk).Xi-shift(1),floenew(kk).Yi-shift(2)]);
                        polynew = subtract(Poly,polylead);
                        polynew2 = subtract(Poly,polynew);
%                         polynew = translate(polynew,[floenew(kk).Xi,floenew(kk).Yi]);
%                         polynew2 = translate(polynew2,[floenew(kk).Xi,floenew(kk).Yi]);
                        poly = [regions(polynew); regions(polynew2)];
                        Floes=fractures(floenew(kk),poly,[floenew(kk).Xi-shift(1),floenew(kk).Yi-shift(2)]);
                        fracturedFloes = [fracturedFloes Floes];
                        alive(kk) = 0;
                    end
                end

                live(live == 1 & Binx == ii & Biny == jj) = alive;
            end
        end
    end
end
Floe = [Floe(logical(live)) fracturedFloes];
for kk =1:length(Floe)
    polyF(kk) = polyshape(Floe(kk).c_alpha'+[Floe(kk).Xi Floe(kk).Yi]);
end
warning('on',id)
warning('on',id3)

end

