function [Floe,Vd] = pack_ice(Floe,c2_boundary,dhdt,Vd,target,ocean, height, min_floe_size, PERIODIC)
%% This function takes in the existing floe state and creates new thin floes in the open space to match the input target concentration
id = 'MATLAB:polyshape:tinyBoundaryDropped';
warning('off',id);
id2 ='MATLAB:polyshape:repairedBySimplify';
warning('off',id2)

N0=length(Floe);
kill = zeros(N0,1);
Lx= max(c2_boundary(1,:));
Ly= max(c2_boundary(2,:));%c2 must be symmetric around x=0 for channel boundary conditions.

Xo=ocean.Xo;
Yo=ocean.Yo;
[Xocn,Yocn] = meshgrid(Xo,Yo);
Uocn=ocean.Uocn;
Vocn=ocean.Vocn;
hO = height;
hvsd = @(x) [0.5*(x == 0) + (x > 0)];

%If the floe is periodic populate the ghost floes
N0=length(Floe);
if PERIODIC
    
    ghostFloeX=[];
    ghostFloeY=[];
    parent=[];
    translation = [];
    
    x=cat(1,Floe.Xi);
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        poly = polyshape(Floe(i).c_alpha'+[x(i) y(i)]);
        %   if alive(i) && (x(i)>Lx-rmax(i)) || (x(i)<-Lx+rmax(i))
        if alive(i) && (max(abs(poly.Vertices(:,1)))>Lx)
            
            ghostFloeX=[ghostFloeX  Floe(i)];
            ghostFloeX(end).Xi=Floe(i).Xi-2*Lx*sign(x(i));
            parent=[parent  i];
            translation = [translation; -2*Lx*sign(x(i)) 0];
            
        end
        
        
    end
    
    Floe=[Floe ghostFloeX];
    
    x=cat(1,Floe.Xi);
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        
        %   if alive(i) && (x(i)>Lx-rmax(i)) || (x(i)<-Lx+rmax(i))
        if alive(i) && (max(abs(poly.Vertices(:,2)))>Ly)
            
            ghostFloeY=[ghostFloeY  Floe(i)];
            ghostFloeY(end).Yi=Floe(i).Yi-2*Ly*sign(y(i));
            parent=[parent  i];
            translation = [translation; 0 -2*Ly*sign(y(i))];
            
        end
        
    end
    
    Floe=[Floe ghostFloeY];
    
end

%Caclulate the coarse grid and any potential interactions for the different
%regions of the floe where ice needs to be created
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
    Floe(kk).poly = polyshape(Floe(kk).c_alpha'+[Floe(kk).Xi Floe(kk).Yi]);
    pint = sqrt((xx-xf(kk)).^2+(yy-yf(kk)).^2)-(rmax(kk)+r_max);
    pint(pint>0) = 0;
    pint(pint<0) = 1;
    potentialInteractions(:,:,kk) = pint;
end

%function to determine the probablity that ice gets created
ramp = @(dhdt) hvsd(dhdt)*dhdt;
SimpMin = @(A) 3*log10(A);

%% Loop through all regions creating sea ice in each if probability criteria is met

for j = 1:Nx*Ny
    ii = ceil(j/Ny);
    jj = mod(j,Ny);
    if jj == 0;  jj = Ny; end
    
    p = rand(1);
    if p <ramp(dhdt)
        %Find coverage of sea ice in region of interest and calculate
        %concentration
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
            polyu = intersect(polyu,box);
        end
        R = regions(polyu);
        FloeArea = sum(area(R));
        c = FloeArea/area(box);
        
        %if concentration is below target then create new sea ice
        if c<0.99*target
            atarget = (target*area(box)-FloeArea);
            floe.poly = box;
            
            if isempty(polyu)
                Op = box;
            else
                Op = subtract(box,polyu);
            end
            
            
            %Break up open water into chunks that will be the new floes
            anew = 0;
            floe.rmax = r_max;
            clear areas
            [Xi,Yi] = centroid(box);
            floe.Xi = Xi; floe.Yi = Yi;
            N = ceil(atarget/min_floe_size);
            if N < 3
                N = 3;
            end
            
            X = Xi+floe.rmax*(2*rand(N,1)-1);
            Y = Yi+floe.rmax*(2*rand(N,1)-1);
            boundingbox=[-1 ,-1; 1,-1; 1,1; -1 ,1]*floe.rmax+ [Xi Yi];
            [~, b,~,~,~] = polybnd_voronoi([X Y],boundingbox);
            
            for kk = 1:length(b)
                if ~isnan(b{kk})
                    subfloes(kk) = polyshape(b{kk});
                end
            end
            
            new = intersect(Op,subfloes);
            subfloes = [];
            for kk = 1:length(new)
                R = regions(new(kk));
                subfloes = [subfloes; R(area(R)>min_floe_size)];
            end
            
            areas = area(subfloes);
            Iold = [];
            N = length(areas(areas>min_floe_size));%only keep floes above a certain area
            count = 1;
            
%             xx = 1; xx(1) =[1 2];
            
            %Add in new floes until we reach target concentration
            while anew < atarget
                
                [~,I] = max(areas);
                if areas(I) > min_floe_size
                    %Cacluate properties of the new floes
                    if sum(size(height.mean)) > 2
                        [Xi,Yi] = centroid(subfloes(I));
                        [ko,~] = dsearchn([Xocn(:),Yocn(:)],[Xi,Yi]);
                        height.mean = hO.mean(ko);
                    end
                    floe2 = initialize_floe_values(subfloes(I),height);
%                     if length(floe2new.poly.Vertices) > SimpMin(floe2new.area)
%                         floeS = FloeSimplify(floe2new,true);
%                         AreaSimp = cat(1,floeS.area);
%                         [~,Imax] = max(AreaSimp);
%                         floe2 = floeS(Imax);
% %                         floe2 = [];
% %                         if length(floeS) >1
% %                             floe2 = floeS(1);
% %                             xx = 1; xx(1) =[1 2];
% %                         else
% %                             floe2 = floeS;
% %                         end
%                     else
%                         floe2 = floe2new;
%                     end
                    Vd(jj,ii,1) = Vd(jj,ii,1)-floe2.h*floe2.area*rho_ice;
%                     [k,~] = dsearchn([Xocn(:),Yocn(:)],[floe2.Xi,floe2.Yi]);
%                     floe2.Ui = Uocn(k); floe2.Vi = Vocn(k);
                    
                    if floe2.poly.NumHoles > 0
                        floe2holes = floe2;
                        if isnan(floe2.ksi_ice)
                            xx = 1;
                            xx(1) = [1 2];
                        end
                        Holes= floe2.poly.NumHoles;
                        poly2 = rmholes(floe2.poly);
                        polyO = intersect(poly2,poly);
                        A2 = area(polyO);
                        k2 = find(logical(potentialInteractions(jj,ii,:))==1);
                        for kk = 1:Holes
                            polyin = A2./abs(area(floe2holes.poly,kk+1));
                            polyin2 = abs(area(floe2holes.poly,kk+1))./A2;
                            polyin2(isinf(polyin2)) = 0;
                            in = unique([k2(polyin>0.99); k2(polyin2>0.99)]);
                            test1 = initialize_floe_values(rmholes(floe2.poly),height);
                            for k = 1:length(in)
                                if Floe(in(k)).alive == 1
                                    V2 = area(intersect(test1.poly,Floe(in(k)).poly))*Floe(in(k)).h;
                                    [test1, test2] = ridge_values_update(test1,Floe(in(k)), V2);
                                    polyu = subtract(polyu,union(test1.poly));
%                                     test = FuseFloes(floe2,FloeIn);
                                    if isempty(test2)
                                        Floe(in(k)).alive = 0;
                                        kill(in(k)) = 1;
                                    elseif length(test2)>1
                                        Floe(in(k)) = test2(1);
                                        floenew = [floenew test2(2:end)];
                                    else
                                        Floe(in(k)) = test2;
                                    end
                                    if length(test1) > 1
                                        floe2 = test1(1);
                                        floenew =[floenew test1(2:end)];
                                    else
                                        floe2 = test1;
                                    end
                                end
                            end
                        end
                        if isnan(floe2.ksi_ice)
                            xx = 1;
                            xx(1) = [1 2];
                        end
                    end
                    
                    if isnan(floe2.ksi_ice)
                        xx = 1;
                        xx(1) = [1 2];
                    end
                    
                    %Keep track of how much area is now covered by sea
                    %ice
                    anew = anew + areas(I);
                    areas(I) = 0;
                    Iold = [Iold I];
                    
                    %Now check for any overlap in holes
                    overlapS = intersect(floe2.poly,[subfloes]);
                    areas(area(overlapS)./areas>0.9) = 0;
                    
                    if floe2.poly.NumHoles > 0
                        xx = 1;
                        xx(1) = [1 2];
                        poly2 = rmholes(floe2.poly);
                        if area(intersect(poly2,polyu))>0
                            k2 = find(logical(potentialInteractions(jj,ii,:))==1);
                            polyO = intersect(poly2,poly);
                            A2 = area(polyO);
                            polyin = A2./A;
                            in = k2(polyin>0.99);
                            in = flipud(in);
                            for kk = 1:length(in)
                                floe2 = FuseFloes(floe2,Floe(in(kk)));
                                Floe(in(kk)).alive = 0;
                                kill = 1;
                                polyu = subtract(polyu,Floe(in(kk)).poly);
                            end
                        end
                    end
                    if isnan(floe2.ksi_ice)
                        xx = 1;
                        xx(1) = [1 2];
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
Floe=Floe(1:N0);
Floe(logical(kill))=[];
alive = cat(1,Floe.alive);
Floe(~logical(alive)) = [];

% %Simplify any new floes that are to complicated
% SimpMin = @(A) 3*log10(A);
% for ii = 1:length(floenew)
%     floe = floenew(ii);
%     if length(floenew(ii).poly.Vertices) > SimpMin(floenew(ii).area)
%         Floe1 = FloeSimplify(floenew(ii));
%         if length(Floe1) >1
%             floenew(ii) = Floe1(1);
%             floenew = [floenew Floe1(2:end)];
%         else
%             floenew(ii) = Floe1;
%         end
%     end
% end

% [eularian_data] = calc_eulerian_data(Floe,10,10,0,c2_boundary,PERIODIC);
% if max(max(eularian_data.c))-1 > 0.01
%     xx = 1;
%     xx(1) = [1 2];
% end

Floe = [Floe floenew];
warning('on',id)
warning('on',id2)

Floe=rmfield(Floe,{'poly'});
Vd(Vd<0)=0;

for ii = 1:length(Floe)
    if isnan(Floe(ii).ksi_ice)
        xx = 1;
        xx(1) = [1 2];
    end
end

end