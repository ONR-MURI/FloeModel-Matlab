function [ force_1, pcontact, overlap,Force_factor] = floe_interactions_poly_con(floe1, floe2, c2_boundary,PERIODIC,Modulus)
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)

%Modulus = 8e7; 
h1 = floe1.h; h2 = floe2.h;
r1 = sqrt(floe1.area); r2 = sqrt(floe2.area);
% Modulus = Modulus*10;
Force_factor=Modulus*(h1*h2)/(h1*r2+h2*r1); overlap = 0;
if isfield(floe2,'alive')
    Force_factor=Modulus*h1/r1;
end
% xx = 1; xx(1) =[1 2];
% if abs(Force_factor-1.5e3)/1.5e3 > 0.1
%     xx = 1; xx(1) =[1 2];
% end
% load F_dir
% Force_factor = 1.5e3;
% load En;
% Force_factor = kmod;
G = 3.8e9; mu = 0.25;
c1=[floe1.c_alpha(1,:)+floe1.Xi; floe1.c_alpha(2,:)+floe1.Yi];
if isfield(floe2,'c')
    c2=floe2.c;
    poly2 = polyshape(c2');
    boundary = 0;
else
    c2 = [floe2.poly.Vertices; floe2.poly.Vertices(1,:)]';
    poly2 = floe2.poly;
    boundary = 1;
end
poly1 = polyshape(c1');
polyA = area(poly1);
polyI = intersect(poly1,poly2);

if  ( max(c1(1,:))<max(c2_boundary(1,:)) && min(c1(1,:))>min(c2_boundary(1,:)) && max(c1(2,:))<max(c2_boundary(2,:)) && min(c1(2,:))>min(c2_boundary(2,:))|| area(poly2)<0.95*area(polyshape(c2_boundary')) || PERIODIC) 
    if area(polyI)/area(poly1) > 0.6
        overlap = Inf;
    elseif area(polyI)/area(poly2) > 0.6
        overlap = -Inf;
    end
    
end
if norm(c1(:,1)-c1(:,end))> 1
    c1(:,length(c1)+1) = c1(:,1);
end
if norm(c2(:,1)-c2(:,end))> 1
    c2(:,length(c2)+1) = c2(:,1);
end

P=InterX(c1,c2);

if isempty(P) || size(P,2)<2 || isinf(overlap)
    force_1=[0 0];
    pcenter=[0 0];
    pcontact=[0 0];    
else
    
    r = regions(polyI); N1 = length(poly1.Vertices); N2 = length(poly2.Vertices);
    Amin =  min([N1,N2])*100/1.75;%min([area(poly1) area(poly2)]);
%     Amin = 1.8e-4*Amin;
%     if Amin < 1000; Amin = 1000; end
    r = r(area(r)>Amin);%1e-4*Amin
    N_contact=length(r);
    
    force_1=zeros(N_contact,2);
    
    pcenter=zeros(N_contact,2);
    pcontact=zeros(N_contact,2);
    
        
    for k=1:N_contact
        
        poly = r(k);
        [cx, cy] = centroid(poly);
%         if isfield(floe2,'alive')  
%             if cy<0 && abs(cx)<9e4
%                 xx = 1; xx(1) =[1 2];
%             end
%         end
        [verts,dist] = dsearchn(poly.Vertices,P');
        p = poly.Vertices(verts(dist<1),:);
        [m,~] = size(p);
        
        if m == 2
            contact = mean(p,1);% [~, x_d_min, y_d_min] = p_poly_dist(contact(1), contact(2), r(k).Vertices(:,1), r(k).Vertices(:,2),true);
            pcontact(k,:)= [cx, cy];% contact;%[x_d_min, y_d_min];
            xv = [p(1,1); p(2,1)];
            xgh = xv(2:end)-xv(1:end-1); xm = (xv(2:end)+xv(1:end-1))/2;
            yv = [p(1,2); p(2,2)];
            ygh = yv(2:end)-yv(1:end-1); ym = (yv(2:end)+yv(1:end-1))/2;
            b = sqrt(xgh.^2+ygh.^2); force_dir = [-ygh./b; xgh./b];
        else
            xv = [poly.Vertices(:,1); poly.Vertices(1,1)];
            xgh = xv(2:end)-xv(1:end-1); xm = (xv(2:end)+xv(1:end-1))/2;
            yv = [poly.Vertices(:,2); poly.Vertices(1,2)];
            ygh = yv(2:end)-yv(1:end-1); ym = (yv(2:end)+yv(1:end-1))/2;
            b = sqrt(xgh.^2+ygh.^2); n = [-ygh./b xgh./b];
            xt = xm+n(:,1)/100; yt = ym+n(:,2)/100;
            in = inpolygon(xt,yt,xv,yv);
            n(~in,:) = -n(~in,:);
            Fn = -Force_factor*(b*ones(1,2)).*n;
            [d_min1] = p_poly_dist(xm, ym,poly1.Vertices(:,1), poly1.Vertices(:,2));
            on = logical(abs(d_min1)<1e-8);
            if length(on)<length(d_min1)
                f_dir = sum(Fn(on,:),1);
                force_dir=f_dir'/sqrt(f_dir*f_dir');
            else
                force_dir = [0; 0];
            end
            pcontact(k,:) = [cx, cy];
        end

        poly1new = translate(poly1,force_dir');
        polyInew = intersect(poly1new,poly2);
        rnew = regions(polyInew); Anew = intersect(r(k),rnew);
        [~,Imax] = max(area(Anew));
        if area(rnew(Imax))/area(r(k))-1 > 0
            force_dir = -force_dir;
            poly1new = translate(poly1,force_dir');
            polyInew = intersect(poly1new,poly2);
            rnew = regions(polyInew); Anew = intersect(r(k),rnew);
            [~,Imax2] = max(area(Anew));
            if Imax2 > Imax
                force_dir = -force_dir;
            end
        end
        
        force=force_dir*area(r(k))*Force_factor; %proportional to the overlap area
                
%         dir_t = [-force_dir(2) force_dir(1)];
        v1 = ([floe1.Ui floe1.Vi]+ floe1.ksi_ice*(pcenter(k,:)-[floe1.Xi floe1.Yi]));
        v2 = ([floe2.Ui floe2.Vi]+ floe2.ksi_ice*(pcenter(k,:)-[floe2.Xi floe2.Yi]));
        dir_t = [-force_dir(2) force_dir(1)];
        v_t = dot(v1-v2,dir_t);
        v_t = (v1-v2);
        if max(abs(v_t)) == 0
            dir_t =[0 0];
        else
            dir_t = v_t/vecnorm(v_t);
        end
        dl = mean(b);
%         dir_t = sign(v_t)*dir_t;
%         v_r = (v1+v2)/2;
%       force_t = dl*G*v_t*sign(v_r-v1).*abs(dir_t);
      force_t = -dl*G*vecnorm(v_t)*dir_t;
%         force_t = [0 0];%-area(r(k))*(v1-v2)*Force_factor/4;%[0 0];%

        if vecnorm(force_t)>mu*vecnorm(force)
            force_t = -mu*vecnorm(force)*dir_t;
        end

        if isnan(force(1)) || isnan(force(2))
            xx = 1;
            xx(1) = [1 2];
        end
        
        force_1(k,:)= force'+force_t;
        overlap(k) = area(r(k));
        
    end

end


warning('on',id)
warning('on',id3)

end

