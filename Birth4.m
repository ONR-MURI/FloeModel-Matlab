function [FloeNEW]= BirthDeath4(Floe,c2_boundary,Target)

Area = cFine0;
Area(cFine0<1)=0;
rho_ice=920;
mean_thickness = 3;
ddx = 250;
[m,n] = size(Area);
cnow = mean(mean(cFine0));
N = 8e2;
Xv = 1.5*x;
Yv = 1.5*y;
X = Xv(1)+rand(1,N)*(max(Xv)-min(Xv));
Y = Yv(1)+rand(1,N)*(max(Yv)-min(Yv));
for ii = length(X):-1:1
    [Mx,jj] = min(abs(x(1:n-1)-X(ii)));
    [My,kk] = min(abs(y(1:m-1)-Y(ii)));
    if Area(kk,jj) > 0
        X(ii)=[];
        Y(ii) = [];
    end
end
dt = delaunayTriangulation(X',Y');
[Vold,C] = voronoiDiagram(dt);
V = Vold;
[vxold,vyold]=voronoi(X,Y);
[vx, vy] = vertexcrop(x,y,vxold,vyold);
for ii = 1:length(vx)
    [M(ii),jj] = min(abs(V(:,1)-vxold(1,ii)));
    if abs(M(ii))==0
        V(jj,1)=vx(1,ii);
        V(jj,2)=vy(1,ii);
    end
end

Atot = zeros(m,n);
%htot = zeros(m,n);
h = 0.5;
E = 1e9;
sigma_m = 0.4e6;
Atot2 = zeros(m,n);
R = randperm(length(C));
hbar = normrnd(mean_thickness,0.2,[1,length(C)]);
%% Merging Floes
% P = [X;Y];
% PQ = [X(1,R(1));Y(1,R(1))];
% P2 = P;
% P2(:,R(1))=[];
% k = dsearchn(P2',PQ');
% [M,Ind] = min(abs(P(1,:)-P2(1,k)));
% 
% ii = 1;
% ind = 1;
% ind2 = [];
% Merg = [R(1) Ind];
% Merg = sort(Merg,'descend');
% for jj = 1:length(Merg)
%     if C{Merg(jj)}(1) ==1
%         C{Merg(jj)}(1)=[];
%         C{Merg(jj)} = C{Merg(jj)}(2:end);
%     end
%     Vnow1 = V(C{Merg(jj)},1);
%     Vnow1(isnan(Vnow1))=[];
%     Vnow2 = V(C{Merg(jj)},2);
%     Vnow2(isnan(Vnow2))=[];
%     if jj == 1
%         poly1 = polyshape(Vnow1,Vnow2);
%     else
%         poly2 = polyshape(Vnow1,Vnow2);
%     end
% end
% polyout = union(poly1,poly2);
% Vnow1 = polyout.Vertices(:,1);
% Vnow2 = polyout.Vertices(:,2);
% Vnow1 = V(C{Merg(jj)},1);
% Vnow1(isnan(Vnow1))=[];
% Vnow2 = V(C{Merg(jj)},2);
% Vnow2(isnan(Vnow2))=[];
% BW = poly2mask((Vnow1+max(x))/ddx,(Vnow2+max(y))/ddx,m, n);
% Anew = double(BW)-double(Area)-double(Atot);
% Anew(Anew<0)=0;
% if sum(sum(Anew)) > 0
%     stats = regionprops(Anew);
%     [areaPoly,jj] = max(cat(1,stats.Area));
%     if isempty(jj) == 0
%         Xi(ind,1) = stats(jj).Centroid(1)*ddx-max(x)-ddx;
%         Yi(ind,1) = stats(jj).Centroid(2)*ddx-max(y)-ddx;
%         total_multiple_contour_floes=0;
%         
%         Anew = logical(Anew);
%         xx=1:size(Anew,2); yy=1:size(Anew,1);
%         [xx,yy]=meshgrid(ddx*xx,ddx*yy); %% ddx is the resolution size of model
%         xc=mean(xx(Anew)); yc=mean(yy(Anew));
%         r=sqrt((xx-xc).^2+(yy-yc).^2);
%         r_max(ind,1)=max(r(Anew))+2*ddx;
%         Ngr=2*fix(r_max(ind,1)/ddx); % number of grid boxes over floe
%         
%         Yg1=floor(stats(jj).BoundingBox(2)+stats(jj).BoundingBox(4)/2-Ngr/2):floor(stats(jj).BoundingBox(2)+stats(jj).BoundingBox(4)/2+Ngr/2);
%         Xg1=floor(stats(jj).BoundingBox(1)+stats(jj).BoundingBox(3)/2-Ngr/2):floor(stats(jj).BoundingBox(1)+stats(jj).BoundingBox(3)/2+Ngr/2);
%         Yg1(Yg1<1) = [];
%         Yg1(Yg1>m)=[];
%         Xg1(Xg1<1) = [];
%         Xg1(Xg1>n) = [];
%         A = Anew(Yg1,Xg1);
%         
%         AreaArray(ind,1) = sum(sum(A))*ddx^2;
%         Mass(ind,1) = AreaArray(ind)*thickness*rho_ice;
%         inertia_moment(ind,1)=sum(r(A).^2)*ddx^2*thickness*rho_ice;
%         
%         Atot = Atot+Anew;
%         cFine = Atot+cFine0;
%         cFine(cFine>1)=1;
%         cnow = mean(mean(cFine));
%         ind2(ind) = ii;
%         ii = ii+1;
%         ind = ind+1;
%         Merged = 1;
%         
%     end
% else
    ii = 1;
    ind = 1;
    Merged = 0;
% end

%% Loop if Concentration is 1
if Target == 1
for ii = 1:length(C)
    if ii>length(R)
        Anew = 0;
    elseif R(ii) > length(C)
        Anew = 0;
    else
        if C{R(ii)}(1) ==1
            C{R(ii)}(1)=[];
            C{R(ii)} = C{R(ii)}(2:end);
        end
        Vnow1 = V(C{R(ii)},1);
        Vnow1(isnan(Vnow1))=[];
        Vnow2 = V(C{R(ii)},2);
        Vnow2(isnan(Vnow2))=[];
        BW = poly2mask((Vnow1+max(x))/ddx,(Vnow2+max(y))/ddx,m, n);
        Anew = double(BW)-double(Area)-double(Atot);
        Anew(Anew<0)=0;
    end
    if sum(sum(Anew)) > 0
        stats = regionprops(Anew);
        [areaPoly,jj] = max(cat(1,stats.Area));
        if isempty(jj) == 0
            Xi(ind,1) = stats(jj).Centroid(1)*ddx-max(x)-ddx;
            Yi(ind,1) = stats(jj).Centroid(2)*ddx-max(y)-ddx;
            total_multiple_contour_floes=0;
            
            Anew = logical(Anew);
            xx=1:size(Anew,2); yy=1:size(Anew,1);
            [xx,yy]=meshgrid(ddx*xx,ddx*yy); %% ddx is the resolution size of model
            xc=mean(xx(Anew)); yc=mean(yy(Anew));
            r=sqrt((xx-xc).^2+(yy-yc).^2); 
            r_max(ind,1)=max(r(Anew))+2*ddx;
            Ngr=2*fix(r_max(ind,1)/ddx); % number of grid boxes over floe
            
            Yg1=floor(stats(jj).BoundingBox(2)+stats(jj).BoundingBox(4)/2-Ngr/2):floor(stats(jj).BoundingBox(2)+stats(jj).BoundingBox(4)/2+Ngr/2);
            Xg1=floor(stats(jj).BoundingBox(1)+stats(jj).BoundingBox(3)/2-Ngr/2):floor(stats(jj).BoundingBox(1)+stats(jj).BoundingBox(3)/2+Ngr/2);
            Yg1(Yg1<1) = [];
            Yg1(Yg1>m)=[];
            Xg1(Xg1<1) = [];
            Xg1(Xg1>n) = [];
            A = Anew(Yg1,Xg1);
            r = r(Yg1,Xg1);
            
            %thickness=2; % ice floe thickness in meters; save for all floes.
            d = logical(abs(double(A)-1));
            d = bwdist(d)-1;
%             h = sin(pi*d)./(pi*d)+hbar(ii);
%             %h = r_max(ind,1)*sin(pi*(r-r_max(ind,1)))./(pi*(r-r_max(ind,1)))+hbar(ii);
%             h(isinf(1./A)) = 0;
%             h(isnan(h)) = 1 + hbar(ii);
%             thickness = abs(h)+rand(length(Yg1),length(Xg1))*(min(abs(h(A))));
%             thickness(isinf(1./A))=0;
            
            AreaArray(ind,1) = sum(sum(A))*ddx^2;
            Mass(ind,1) = AreaArray(ind).*h*rho_ice;
            inertia_moment(ind,1)=sum(r(A).^2)*ddx^2.*h*rho_ice;
            hArray(ind,1) = h;
            Earray(ind,1) = E;
            sig_array(ind,1) = sigma_m;
            
            Atot = Atot+Anew;
%             htot(Yg1,Xg1) = htot(Yg1,Xg1)+Atot(Yg1,Xg1).*thickness;
            cFine = Atot+cFine0;
            cFine(cFine>1)=1;
            cnow = mean(mean(cFine));
            ind2(ind) = ii;
            ind = ind+1;
        end
    end
    if ii >= length(C)
        cnow = 1;
    end
    ii = ii+1;    
end
end
%% Loop if Concentration < 1
if Target < 1
while cnow < Target
    if ii>length(R)
        Anew = 0;
    elseif R(ii) > length(C)
        Anew = 0;
    else
        if C{R(ii)}(1) ==1
            C{R(ii)}(1)=[];
            C{R(ii)} = C{R(ii)}(2:end);
        end
        Vnow1 = V(C{R(ii)},1);
        Vnow1(isnan(Vnow1))=[];
        Vnow2 = V(C{R(ii)},2);
        Vnow2(isnan(Vnow2))=[];
        BW = poly2mask((Vnow1+max(x))/ddx,(Vnow2+max(y))/ddx,m, n);
        Anew = double(BW)-double(Area)-double(Atot);
        Anew(Anew<0)=0;
    end
    if sum(sum(Anew)) > 0
        stats = regionprops(Anew);
        [areaPoly,jj] = max(cat(1,stats.Area));
        if isempty(jj) == 0
            Xi(ind,1) = stats(jj).Centroid(1)*ddx-max(x)-ddx;
            Yi(ind,1) = stats(jj).Centroid(2)*ddx-max(y)-ddx;
            total_multiple_contour_floes=0;
            
            Anew = logical(Anew);
            xx=1:size(Anew,2); yy=1:size(Anew,1);
            [xx,yy]=meshgrid(ddx*xx,ddx*yy); %% ddx is the resolution size of model
            xc=mean(xx(Anew)); yc=mean(yy(Anew));
            r=sqrt((xx-xc).^2+(yy-yc).^2); 
            r_max(ind,1)=max(r(Anew))+2*ddx;
            Ngr=2*fix(r_max(ind,1)/ddx); % number of grid boxes over floe
            
            Yg1=floor(stats(jj).BoundingBox(2)+stats(jj).BoundingBox(4)/2-Ngr/2):floor(stats(jj).BoundingBox(2)+stats(jj).BoundingBox(4)/2+Ngr/2);
            Xg1=floor(stats(jj).BoundingBox(1)+stats(jj).BoundingBox(3)/2-Ngr/2):floor(stats(jj).BoundingBox(1)+stats(jj).BoundingBox(3)/2+Ngr/2);
            Yg1(Yg1<1) = [];
            Yg1(Yg1>m)=[];
            Xg1(Xg1<1) = [];
            Xg1(Xg1>n) = [];
            A = Anew(Yg1,Xg1);
            r = r(Yg1,Xg1);
            
            %thickness=2; % ice floe thickness in meters; save for all floes.
            d = logical(abs(double(A)-1));
            d = bwdist(d)-1;
%             h = sin(pi*d)./(pi*d)+hbar(ii);
%             %h = r_max(ind,1)*sin(pi*(r-r_max(ind,1)))./(pi*(r-r_max(ind,1)))+hbar(ii);
%             h(isinf(1./A)) = 0;
%             h(isnan(h)) = 1 + hbar(ii);
%             thickness = abs(h)+rand(length(Yg1),length(Xg1))*(min(abs(h(A))));
%             thickness(isinf(1./A))=0;
            
            AreaArray(ind,1) = sum(sum(A))*ddx^2;
            Mass(ind,1) = sum(sum(A.*h))*rho_ice;%AreaArray(ind).*thickness*rho_ice;
            inertia_moment(ind,1)=sum(h*r(A).^2)*ddx^2.*rho_ice;%sum(r(A).^2)*ddx^2.*thickness*rho_ice;
            hArray(ind,1) = h;
            Earray(ind,1) = E;
            sig_array(ind,1) = sigma_m;
            
            Atot = Atot+Anew;
%             htot(Yg1,Xg1) = htot(Yg1,Xg1)+Atot(Yg1,Xg1).*thickness;
            cFine = Atot+cFine0;
            cFine(cFine>1)=1;
            cnow = mean(mean(cFine));
            ind2(ind) = ii;
            ind = ind+1;
        end
    end
    if ii >= length(C)
        cnow = 1;
    end
    ii = ii+1;    
end
end

%% Create Structure
if isempty(ind2) == 0
    col0 = logical(zeros(length(ind2),1));
    col1 = logical(ones(length(ind2),1));
    T = table(AreaArray,Mass,inertia_moment,col0,col0,col0,r_max,col0,col0,col0,col0,Xi,Yi,col0,col0,col0,col0,col1,...
        col0,col0,col0,col0,col0,col0,col0,col0,col0,col0,hArray,Earray,sig_array,...
        'VariableNames',{'area','mass','inertia_moment','c0','c_alpha','A','rmax','Xg','Yg',...
        'X','Y','Xi','Yi','alpha_i','Ui','Vi','ksi_ice','alive','dXi_p','dYi_p','dUi_p','dVi_p',...
        'dalpha_i_p','dksi_ice_p','interactions','potentialInteractions','collision_force','collision_torque','h','E','sigma_m'});
    FloeNEW = table2struct(T);
    FloeNEW = FloeNEW';
    clear r_max
    clear T
end
%% Fill in structure
Atot = zeros(m,n);
length(ind2);
if isempty(ind2) == 0
    if isnan(length(ind2)) == 0
        for ii = 1:length(ind2)
            if Merged ==1
                if ii == 1
                    Vnow1 = polyout.Vertices(:,1);
                    Vnow2 = polyout.Vertices(:,2);
                else
                    Vnow1 = V(C{R(ind2(ii))},1);
                    Vnow1(isnan(Vnow1))=[];
                    Vnow2 = V(C{R(ind2(ii))},2);
                    Vnow2(isnan(Vnow2))=[];
                end
            else
                Vnow1 = V(C{R(ind2(ii))},1);
                Vnow1(isnan(Vnow1))=[];
                Vnow2 = V(C{R(ind2(ii))},2);
                Vnow2(isnan(Vnow2))=[];
            end
            BW = poly2mask((Vnow1+max(x))/ddx,(Vnow2+max(y))/ddx,m, n);
            Anew = double(BW)-double(Area)-double(Atot);
            Anew(Anew<0)=0;
            
            if sum(sum(Anew)) > 0
                stats = regionprops(Anew);
                if isempty(stats.Area) == 0
                %[areaPoly,jj] = max(cat(1,stats.Area));
                
                    total_multiple_contour_floes=0;
                    
                    Anew = logical(Anew);
                    xx=1:size(Anew,2); yy=1:size(Anew,1);
                    [xx,yy]=meshgrid(ddx*xx,ddx*yy); %% floePixelSize m is the pixel size from Earthdata NASA webpage
                    xc=mean(xx(Anew)); yc=mean(yy(Anew));
                    r=sqrt((xx-xc).^2+(yy-yc).^2); r_max=max(r(Anew))+2*ddx; % pad with extra row/column of zeros around the floe
                    
                    %h = r_max*sin(10*pi*(r-r_max)/r_max)./(pi*(r-r_max))+1;                  
                    
                    Ngr=2*fix(r_max/ddx); % number of grid boxes over floe
                    
                    [XX,YY]=meshgrid((-1:2/Ngr:1)*r_max, (-1:2/Ngr:1)*r_max); % floe grid
                    Xg=(-1:2/Ngr:1)*r_max; Yg=Xg;
                    
                    A=interp2(xx-xc,yy-yc,double(Anew),XX,YY); A(isnan(A))=0;
                    A = logical(A);
                    %thickness = interp2(xx-xc,yy-yc,h,XX,YY); thickness(isinf(1./A))=0;
                    
                    
                    c0=contourc(Xg,Yg,double(A),[0.5 0.5]);
                    
                    if isempty(c0) == 1
                        FloeNEW(ii).alive=0;
                    else
                        if c0(2,1)< size(c0,2)-1,
                            multiple_contours=1;
                            total_multiple_contour_floes=total_multiple_contour_floes+1;
                        else multiple_contours=0;
                        end
                        
                        c0(:,c0(1,:)==0.5)=NaN;
                        c0=c0(:,2:end);
                        
                        alpha_ini=0;
                        A_rot=[cos(alpha_ini) -sin(alpha_ini); sin(alpha_ini) cos(alpha_ini)]; %rotation matrix
                        c_alpha=A_rot*c0;
                        
                        FloeNEW(ii).c0=c0; %contour
                        FloeNEW(ii).c_alpha=c_alpha;
                        FloeNEW(ii).A=uint8(A);
                        FloeNEW(ii).Xg=Xg;
                        FloeNEW(ii).Yg=Yg;
                        FloeNEW(ii).X=XX;
                        FloeNEW(ii).Y=YY;
                    end
                end
            end
        end
    else
        FloeNEW = [];
    end
    else
    FloeNEW = [];
end