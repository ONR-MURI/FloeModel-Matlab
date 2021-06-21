function [Floe, Nb] = initial_concentration(c2_boundary,target_concentration,height, NumFloes, min_floe_size)
%% This function is used to generate the initial floe field

%Identify the grids to align with the concentrations specified by the input
[Ny, Nx] = size(target_concentration);
c = flipud(target_concentration);
x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
c2_boundary_poly = polyshape(c2_boundary');
dx = x(2)-x(1);
dy = y(2)-y(1);

%Create floes that act as boundaries and wont move
% x1 = [Lx/2 Lx/2 1.5e4 1.5e4]; y1 = [-Ly/2 -1e4 -Ly/2+3e5 -Ly/2];
% B1 = polyshape(x1, y1);
% B2 = polyshape(-x1, y1);
% x2 = [Lx/2 Lx/2 -Lx/2 -Lx/2]; y2 = [-1e5 -Ly/2 -Ly/2 -1e5];
% B3 = polyshape(x2,y2);
% Floe1 = initialize_floe_values(B1,height,0,SUBFLOES);
% Floe2 = initialize_floe_values(B2,height,0,SUBFLOES);
% bound = subtract(c2_boundary_poly, B1);
% bound = subtract(bound, B2);
%Floe = [Floe1 Floe2];

nx=40; ny=4;%fix(Nx*LyO/LxO);
xc = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/nx:max(c2_boundary(1,:));
yc = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/ny:max(c2_boundary(2,:));
Xc = (xc(1:end-1)+xc(2:end))/2; Yc = -(yc(1:end-1)+yc(2:end))/2;
[xx,yy] = meshgrid(Xc,Yc); X = xx(:); Y = yy(:);%[yc(1) yc(2) yc(2) yc(1)];

Nb = 0;
Floe = [];

%Loop through all the regions of the domain to create new floes
for jj = 1:Ny
    for ii = 1:Nx
        if c(jj,ii)>0
            boundary = polyshape([x(ii) x(ii) x(ii+1) x(ii+1)], [y(jj) y(jj+1) y(jj+1) y(jj)]);
            boundary = intersect(boundary,c2_boundary_poly);
            N = ceil(NumFloes*area(boundary)/area(c2_boundary_poly)/c(jj,ii));
%             poly = intersect(bound,boundary); %Use these when having
%             boundaries
%             N = 4*ceil(NumFloes*area(poly)/area(bound)/c(jj,ii)); %Use these when having
%             boundaries
            X = 0.95*dx/2*(2*rand(N,1)-1)+(x(ii)+x(ii+1))/2;
            Y = 0.95*dy/2*(2*rand(N,1)-1)+(y(jj)+y(jj+1))/2;
            in = inpolygon(X,Y,boundary.Vertices(:,1),boundary.Vertices(:,2));
            X = X(in); Y = Y(in);
%             for i = 1:2%nx
%                 X = [(-1)^i*dx/2 (-1)^i*dx/2 0 0];%xx(:,i:i+1); X = X(:);
%                 Y = [-dy/2 dy/2 dy/2 -dy/2];%Y = yy(:,i:i+1); Y = Y(:);
%                 b{i} = [X',Y'];
%             end
           [~, b,~,~,~] = polybnd_voronoi([X Y],boundary.Vertices);
%             [~, b,~,~,~] = polybnd_voronoi([X Y],boundary.Vertices); %%Use these when having
%             boundaries
            Nf = 1:length(b);%randperm(length(b));
            Atot = 0;
            count = 1;
            while Atot/area(boundary)<=c(jj,ii)
                if ~isnan(b{Nf(count)})
                    poly = polyshape(b{Nf(count)});
                    floenew = initialize_floe_values(poly,height);
                    Floe = [Floe floenew];
                    Atot = Atot+area(poly);
                end
                count = count+1;
                if count > length(Nf)
                    Atot = area(boundary)+1;
                end
            end
        end
    end
end

areas = cat(1,Floe.area);
Floe(areas<min_floe_size)=[];
end

