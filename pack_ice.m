function [Floe,Vd] = pack_ice(Floe,c2_boundary,dhdt,Vd,target,ocean, height, SUBFLOES, PERIODIC)
%% This function takes in the existing floe state and creates new thin floes in the open space to match the input target concentration
id = 'MATLAB:polyshape:tinyBoundaryDropped';
SHIFT = false;
warning('off',id);

N0=length(Floe);
Lx= max(c2_boundary(1,:));
Ly= max(c2_boundary(2,:));%c2 must be symmetric around x=0 for channel boundary conditions.

Xo=ocean.Xo;
Yo=ocean.Yo;
[Xocn,Yocn] = meshgrid(Xo,Yo);
Uocn=ocean.Uocn;
Vocn=ocean.Vocn;

%If the floe is periodic populate the ghost floes
if PERIODIC
    
    ghostFloeX=[];
    ghostFloeY=[];
    parent=[];
    
    x=cat(1,Floe.Xi);
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        
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
    pint = sqrt((xx-xf(kk)).^2+(yy-yf(kk)).^2)-(rmax(kk)+r_max);
    pint(pint>0) = 0;
    pint(pint<0) = 1;
    potentialInteractions(:,:,kk) = pint;
end

%function to determine the probablity that ice gets created
ramp = @(dhdt) heaviside(dhdt)*dhdt;

%% Loop through all regions creating sea ice in each if probability criteria is met
for ii = 1:Nx
    for jj = 1:Ny
        
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
            end
            c = sum(A)/area(box);
            
            %if concentration is below target then create new sea ice
            if c<0.99*target
                atarget = (target*area(box)-sum(A));
                if isempty(polyu)
                    floe.poly = box;
                else
                    floe.poly = subtract(box,polyu);
                end

                %Break up open water into chunks that will be the new floes
                anew = 0;
                floe.rmax = r_max;
                clear areas
                [Xi,Yi] = centroid(box);
                floe.Xi = Xi; floe.Yi = Yi;
                N = 10; N2 = 1;
                while N2 > 0.5
                    [subfloes,~] = create_subfloes(floe,N,false);
                    if sum([subfloes.NumRegions]) == length(subfloes) 
                        N2 = 0;
                    else 
                        N = N+5;
                    end
                end
                areas = area(subfloes);
                Iold = [];
                N = length(areas(areas>3500));%only keep floes above a certain area
                count = 1;
                
                %Add in new floes until we reach target concentration
                while anew < atarget
                    
                    [~,I] = max(areas);
                    if areas(I) > 3500
                        %Cacluate properties of the new floes
                        floe2 = initialize_floe_values(subfloes(I),height, SHIFT, SUBFLOES);
                        Vd(jj,ii,1) = Vd(jj,ii,1)-floe2.h*floe2.area*rho_ice;
                        [k,~] = dsearchn([Xocn(:),Yocn(:)],[floe2.Xi,floe2.Yi]);
                        floe2.Ui = Uocn(k); floe2.Vi = Vocn(k);
                        if floe2.poly.NumHoles > 0
                            if area(intersect(floe2.poly,polyu))>0
                                k2 = find(logical(potentialInteractions(jj,ii,:))==1);
                                poly2 = rmholes(floe2.poly);
                                polyO = intersect(poly2,poly);
                                A2 = area(polyO);
                                polyin = A2./A;
                                in = k2(polyin>0.99);
                                in = flipud(in);
                                for kk = 1:length(in)
                                    floe2 = FuseFloes(floe2,Floe(in(kk)),SUBFLOES);
                                    Floe(in(kk)).alive = 0;
                                    polyu = subtract(polyu,Floe(in(kk)).poly);
                                end
                            end
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
                            floe2holes = floe2;
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
                                for k = 1:length(in)
                                    test = FuseFloes(floe2,Floe(in(k)),SUBFLOES);
                                    if isempty(test)
                                        xx = 1;
                                        xx(1) = [1 2];
                                    end
                                    floe2 = test;
                                    Floe(in(k)).alive = 0;
                                    polyu = subtract(polyu,Floe(in(k)).poly);
                                end
                            end
                        end
                        
                        if numel(fieldnames(floe2))>31
                            floe2=rmfield(floe2,{'potentialInteractions'});
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

%Simplify any new floes that are to complicated
SimpMin = @(A) 15*log10(A);
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
end
if numel(fieldnames(Floe))>31
    Floe=rmfield(Floe,{'potentialInteractions'});
end
Floe = [Floe floenew];
warning('on',id)
end