function [Floe,Vd] = pack_ice(Floe,c2_boundary,dhdt,Vd,target,ocean, height, SUBFLOES, PERIODIC)
%UNTITLED2 Summary of this function goes here
%% 
id = 'MATLAB:polyshape:tinyBoundaryDropped';
SHIFT = false;
warning('off',id);

FloeOld = Floe;

N0=length(Floe);
Lx= max(c2_boundary(1,:));
Ly= max(c2_boundary(2,:));%c2 must be symmetric around x=0 for channel boundary conditions.

Xo=ocean.Xo;
Yo=ocean.Yo;
[Xocn,Yocn] = meshgrid(Xo,Yo);
Uocn=ocean.Uocn;
Vocn=ocean.Vocn;

if PERIODIC
    
    ghostFloeX=[];
    ghostFloeY=[];
    parent=[];
    
    x=cat(1,Floe.Xi);
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        
        %   if alive(i) && (x(i)>Lx-rmax(i)) || (x(i)<-Lx+rmax(i))
        if alive(i) && (max(abs(Floe(i).poly.Vertices(:,1)))>Lx)
            
            ghostFloeX=[ghostFloeX  Floe(i)];
            ghostFloeX(end).poly=translate(Floe(i).poly,[-2*Lx*sign(x(i)) 0]);
            for ii = 1:length(Floe(i).SubFloes)
                ghostFloeX(end).SubFloes(ii).poly=translate(Floe(i).SubFloes(ii).poly,[-2*Lx*sign(x(i)) 0]);
            end
            ghostFloeX(end).Xi=Floe(i).Xi-2*Lx*sign(x(i));
            ghostFloeX(end).vorX=Floe(i).vorX-2*Lx*sign(x(i));
            ghostFloeX(end).vorbox(:,1)=Floe(i).vorbox(:,1)-2*Lx*sign(x(i));
            parent=[parent  i];
            
        end
        
        
    end
    
    Floe=[Floe ghostFloeX];
    
    x=cat(1,Floe.Xi);
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        
        %   if alive(i) && (x(i)>Lx-rmax(i)) || (x(i)<-Lx+rmax(i))
        if alive(i) && (max(abs(Floe(i).poly.Vertices(:,2)))>Ly)
            
            ghostFloeY=[ghostFloeY  Floe(i)];
            ghostFloeY(end).poly=translate(Floe(i).poly,[0 -2*Ly*sign(y(i))]);
            for ii = 1:length(Floe(i).SubFloes)
                ghostFloeY(end).SubFloes(ii).poly=translate(Floe(i).SubFloes(ii).poly,[0 -2*Ly*sign(y(i))]);
            end
            ghostFloeY(end).Yi=Floe(i).Yi-2*Ly*sign(y(i));
            ghostFloeY(end).vorY=Floe(i).vorY-2*Ly*sign(y(i));
            ghostFloeY(end).vorbox(:,2)=Floe(i).vorbox(:,2)-2*Ly*sign(y(i));
            parent=[parent  i];
            
        end
        
    end
    
    Floe=[Floe ghostFloeY];
    
end

floe2 = [];
floenew = [];
rho_ice = 920;
[Ny,Nx,~] = size(Vd);
x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
y = fliplr(y);
dx = abs(x(2)-x(1));
dy = abs(y(2)-y(1));
[xx,yy] = meshgrid(0.5*(x(1:end-1)+x(2:end)),0.5*(y(1:end-1)+y(2:end)));
xf = cat(1,Floe.Xi);
yf = cat(1,Floe.Yi);
r_max = sqrt((dx/2)^2+(dy/2)^2);
rmax = cat(1,Floe.rmax);
potentialInteractions = zeros(Ny,Nx,length(Floe));
for kk = 1:length(Floe)
    pint = sqrt((xx-xf(kk)).^2+(yy-yf(kk)).^2)-(rmax(kk)+r_max);
    pint(pint>0) = 0;
    pint(pint<0) = 1;
    potentialInteractions(:,:,kk) = pint;
end
%[eularian_data] = calc_eulerian_data2(Floe,Nx,Ny,c2_boundary);

%make actual function where i can input dhdt to get valuef
ramp = @(dhdt) heaviside(dhdt)*dhdt;
%% 

for ii = 1:Nx
    for jj = 1:Ny
        
        p = rand(1);
        if p <ramp(dhdt)
            bound = [x(ii) x(ii) x(ii+1) x(ii+1) x(ii);y(jj) y(jj+1) y(jj+1) y(jj) y(jj)];
            box = polyshape(bound(1,:), bound(2,:));
            k = find(logical(potentialInteractions(jj,ii,:))==1);
            if isempty(k)
                poly = [];
                A = 0;
                in = [];
            else
                poly = intersect(box,[Floe(logical(potentialInteractions(jj,ii,:))).poly]);
                A = area(poly);
                in = k(A>0);
            end
            
            if isempty(in)
                polyu = [];
            else
                polyu = union([Floe(in).poly]);
            end
%             polyu = union([poly(A>0)]);
            c = sum(A)/area(box);
            if c<0.99*target
                atarget = (target*area(box)-sum(A));
                if isempty(polyu)
                    floe.poly = box;
                else
                    floe.poly = subtract(box,polyu);
                end
%                 polyout = sortregions(polyout,'area','descend');
%                 R = regions(polyout);
%                 poly1new = R(1);
%                 floe.poly = rmholes(poly1new);
                anew = 0;
                floe.rmax = r_max;
                clear areas
                [Xi,Yi] = centroid(box);
                floe.Xi = Xi; floe.Yi = Yi;
                N = 10; N2 = 1;
                while N2 > 0.5
                    [subfloes,~] = create_subfloes(floe,N,false);
                    if sum([subfloes.NumRegions]) == length(subfloes) %&& sum([subfloes.NumHoles])<0.5
                        N2 = 0;
%                     elseif sum([subfloes.NumRegions]) == length(subfloes) && sum([subfloes.NumHoles])>0.5
%                         floes = subfloes(cat(1,subfloes.NumHoles)>0);
%                         k = find(logical(potentialInteractions(jj,ii,:))==1);
%                         for kk = 1:length(floes)
%                             poly2 = rmholes(floes(kk));
%                             polyO = intersect(poly2,poly);
%                             A2 = area(polyO);
%                             polyin = A2./A;
%                             in = k(polyin>0.99);
%                             figure
%                             plot(floes(kk))
%                             hold on
%                             plot([Floe(in).poly])
%                             xx = 1;
%                             xx(1) = [1 2];
%                             Floe(i) = FuseFloes(Floe(i),Floe(j),SUBFLOES);
%                             Floe(j) = [];
%                         end
                    else 
                        N = N+5;
                    end
                end
%                 subfloes = subtract([subfloes],polyu);
                areas = area(subfloes);
                N = length(areas(areas>3500));
                count = 1;
                while anew < atarget
                    
                    [~,I] = max(areas);
                    if areas(I) > 3500
                        
                        floe2 = initialize_floe_values(subfloes(I),height, SHIFT, SUBFLOES);
                        [k,~] = dsearchn([Xocn(:),Yocn(:)],[floe2.Xi,floe2.Yi]);
                        floe2.Ui = Uocn(k); floe2.Vi = Vocn(k);
                        if floe2.poly.NumHoles>0
                            k = find(logical(potentialInteractions(jj,ii,:))==1);
%                         for kk = 1:length(floes)
                            poly2 = rmholes(floe2.poly);
                            polyO = intersect(poly2,poly);
                            A2 = area(polyO);
                            polyin = A2./A;
                            in = k(polyin>0.99);
%                             figure
%                             plot(floe2.poly)
%                             hold on
%                             plot([Floe(in).poly])
                            in = flipud(in);
                            for kk = 1:length(in)
                                floe2 = FuseFloes(floe2,Floe(in(kk)),SUBFLOES);
                                Floe(in(kk)).alive = 0;
                                polyu = subtract(polyu,Floe(in(kk)).poly);
                            end
%                             figure
%                             plot(floe2.poly)
%                             xx = 1;
%                             xx(1) = [1 2];
                        else
                            floe2 = rmfield(floe2, 'potentialInteractions');
                            if floe2.area/area(box)>0.75
                                xx = 1;
                                xx(1) = [1 2];
                            end
                        end
                        Vd(jj,ii,1) = Vd(jj,ii,1)-floe2.h*floe2.area*rho_ice;
                        anew = anew + areas(I);
                        areas(I) = 0;
                        
                        if floe2.poly.NumHoles > 0
                            Holes= floe2.poly.NumHoles;
                            poly2 = rmholes(floe2.poly);
                            polyO = intersect(poly2,poly);
                            A2 = area(polyO);
                            for kk = 1:Holes
                                polyin = A2./abs(area(floe2.poly,kk+1));
                                in = k(polyin>0.99);                                                          
                                floe2 = FuseFloes(floe2,Floe(in),SUBFLOES);
                                Floe(in(kk)).alive = 0;
                                polyu = subtract(polyu,Floe(in).poly);
                            end
                            if floe2.poly.NumHoles > 0
                                xx = 1;
                                xx(1) = [1 2];
                            end
                            %  floe2.poly = rmholes(floe2.poly);
                        elseif isempty(floe2.SubFloes)
                            xx = 1;
                            xx(1) = [1 2];
                        elseif floe2.poly.NumRegions>1
                            xx = 1;
                            xx(1) = [1 2];
                        elseif ~isempty(polyu)
                            if area(intersect(floe2.poly,polyu))/floe2.area > 0.25
                                xx = 1;
                                xx(1) = [1 2];
                            end
                        end
                        
                        floenew = [floenew floe2];
                        clear floe2;
                    end
                    if count >= N
                        anew = atarget;
                    end
                    count = count+1;
                end
               
            end
                        
        end
    end
end
Floe=Floe(1:N0);
alive = cat(1,Floe.alive);
Floe(~logical(alive)) = [];
SimpMin = @(A) 15*log10(A);%15+(A-1e4)*(1e9-1e4)/(200-15);
for ii = 1:length(floenew)
    floe = floenew(ii);
    if length(floenew(ii).poly.Vertices) > SimpMin(floenew(ii).area)
        Floe1 = FloeSimplify(floenew(ii), 250,SUBFLOES);
        if length(Floe1) >1
            floenew(ii) = Floe1(1);
            floenew = [floenew Floe1(2:end)];
        else
            floenew(ii) = Floe1;
        end
    end
    if floenew(ii).poly.NumRegions>1
        xx = 1;
        xx(1) = [1 2];
    end
end
if numel(fieldnames(Floe))>31
    Floe=rmfield(Floe,{'potentialInteractions'});
end
Floe = [Floe floenew];
warning('on',id)
end