function [ force_1, pcontact, overlap,Force_factor] = floe_interactions_poly(floe1, floe2, c2_boundary,PERIODIC,Modulus)
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
G = 3.8e9; mu = 0.45;
c1=[floe1.c_alpha(1,:)+floe1.Xi; floe1.c_alpha(2,:)+floe1.Yi];
if isfield(floe2,'c')
    c2=floe2.c;
    poly2 = polyshape(c2');
    boundary = 0;
else
    c2 = c2_boundary;
    poly2 = floe2.poly;
    boundary = 1;
end
poly1 = polyshape(c1');
polyA = area(poly1);
polyI = intersect(poly1,poly2);
% poly2 = polyshape(c2');
% polyB = polyshape(c2_boundary');
% if area(intersect(polyB,poly2))/area(polyB)<0.95
%     polyI = intersect(poly1,poly2);
%     if area(polyI)/area(poly1) > 0.6
%         overlap = Inf;
%     elseif area(polyI)/area(poly2) > 0.6
%         overlap = Inf;
%     end
% end
if  ( max(c1(1,:))<max(c2_boundary(1,:)) && min(c1(1,:))>min(c2_boundary(1,:)) && max(c1(2,:))<max(c2_boundary(2,:)) && min(c1(2,:))>min(c2_boundary(2,:))|| area(poly2)<0.95*area(polyshape(c2_boundary')) || PERIODIC) 
    if area(polyI)/area(poly1) > 0.6
        overlap = Inf;
%         xx = 1;
%         xx(1) =[1 2];
    elseif area(polyI)/area(poly2) > 0.6
        overlap = -Inf;
%         xx = 1;
%         xx(1) =[1 2];
    end
    
%     if area(polyI)/area(poly1) > 0.25 && area(polyI)/area(poly2) > 0.25
%         xx = 1;
%         xx(1) =[1 2];
%     end
end
if norm(c1(:,1)-c1(:,end))> 1
    c1(:,length(c1)+1) = c1(:,1);
end
if norm(c2(:,1)-c2(:,end))> 1
    c2(:,length(c2)+1) = c2(:,1);
end

P=InterX(c1,c2);
% xx = 1; xx(1) =[1 2];
% [xint,yint] = polyxpoly(c1(1,:),c1(2,:),c2(1,:),c2(2,:),'unique');
% P = [xint'; yint'];

if isempty(P) || size(P,2)<2 || isinf(overlap)
    force_1=[0 0];
    pcenter=[0 0];
    pcontact=[0 0];    
else
    
%     [P, worked] = sort_contact_points(P, c1, c2 );
    r = regions(polyI);
    r = r(area(r)>100);
    N_contact=length(r);%size(P,2)-1;
%     if abs(N_contact-length(fdir))>0
%         clear fdir
%         fdir.dir =[0;0];fdir.loc =[0 0];
%     end
%     Force_factor = Force_factor*N_contact;
    
    force_1=zeros(N_contact,2);
    
    pcenter=zeros(N_contact,2);
    pcontact=zeros(N_contact,2);
    
    same_dir=zeros(1,N_contact);
        
    for k=1:N_contact
        [d_min1] = p_poly_dist(P(1,:)', P(2,:)', r(k).Vertices(:,1), r(k).Vertices(:,2));
        [~,I] = mink(d_min1,2);
        p = P(:,I);
%         p=P(1:2,k:k+1);

        dl=sqrt((p(1,1)-p(1,2))^2+(p(2,1)-p(2,2))^2);
        
        contact = mean(p,2); [~, x_d_min, y_d_min] = p_poly_dist(contact(1), contact(2), r(k).Vertices(:,1), r(k).Vertices(:,2),true);
        pcenter(k,:)=[x_d_min, y_d_min]; dp=p(:,1)-p(:,2);
        
        d=1e5; %checking intersections with countour at +-d distance from the center point.
        
        dirr=dp(1)/dp(2);
        if abs(dirr)>1e3;
            
            x1=pcenter(k,1);
            y1=pcenter(k,2) - d;
            
            x2=pcenter(k,1);
            y2=pcenter(k,2) + d;
            
        else
            x1=pcenter(k,1) + d/sqrt(1+dirr^2);
            y1=pcenter(k,2) - d/sqrt(1+dirr^2)*dirr;
            
            x2=pcenter(k,1) - d/sqrt(1+dirr^2);
            y2=pcenter(k,2) + d/sqrt(1+dirr^2)*dirr;
        end
        
        pp1=InterX(c1,[x1 x2; y1 y2]);   pp2=InterX(c2,[x1 x2; y1 y2]);

        %[min_dist_center1, i_min] = min( sqrt( ( (pp1(1,:)-pcenter(k,1)).^2 + (pp1(2,:)-pcenter(k,2)).^2)));
        if ~isempty(pp1) && ~isempty(pp2)
            [sorted_dist_center1, i_min] = sort( sqrt( ( (pp1(1,:)-pcenter(k,1)).^2 + (pp1(2,:)-pcenter(k,2)).^2)));
            min_dist_center1=sorted_dist_center1(1);
            
            f_dir=pp1(:,i_min(end))-pp1(:,i_min(1)); %from point with min distance to the furthest point.
            pp1=pp1(:,i_min(1));
            
            %[min_dist_center2, i_min] = min(sqrt( ( (pp2(1,:)-pcenter(k,1)).^2 + (pp2(2,:)-pcenter(k,2)).^2)));
            [sorted_dist_center2, i_min] = sort(sqrt( ( (pp2(1,:)-pcenter(k,1)).^2 + (pp2(2,:)-pcenter(k,2)).^2)));
            min_dist_center2=sorted_dist_center2(1);
            pp2=pp2(:,i_min(1));
            
            if max(abs(pp2-pp1))==0
                force_1(k,:) = [0 0];
                overlap(k) = 0;
%                 pp2=InterX(c2,[x1 x2; y1 y2]);
%                 if size(pp2,2)>1
%                     [sorted_dist_center2, i_min] = sort(sqrt( ( (pp2(1,:)-pcenter(k,1)).^2 + (pp2(2,:)-pcenter(k,2)).^2)));
%                     min_dist_center2=sorted_dist_center2(1);
%                     pp2=pp2(:,i_min(2));
%                 else
%                     pp1=InterX(c1,[x1 x2; y1 y2]);
%                     [sorted_dist_center1, i_min] = sort( sqrt( ( (pp1(1,:)-pcenter(k,1)).^2 + (pp1(2,:)-pcenter(k,2)).^2)));
%                     min_dist_center1=sorted_dist_center1(1);
%                     f_dir=pp1(:,i_min(end))-pp1(:,i_min(1)); %from point with min distance to the furthest point.
%                     pp1=pp1(:,i_min(2));
%                     [sorted_dist_center2, i_min] = sort(sqrt( ( (pp2(1,:)-pcenter(k,1)).^2 + (pp2(2,:)-pcenter(k,2)).^2)));
%                     min_dist_center2=sorted_dist_center2(1);
%                     pp2=pp2(:,i_min(1));
%                 end
%                 if max(abs(pp2-pp1))==0
%                     xx=1;
%                     xx(1) =[1 2];
%                 end
            else
                same_dir(k)=((pp2-pp1)'*f_dir)>0;
                [cx,cy] = centroid(r(k));
%                 cx = pcenter(k,1); cy = pcenter(k,2);
                pcontact(k,:) = [cx, cy];
                xp = cx; yp = cy;
                xv = [p(1,1) p(1,2)]; yv =[p(2,1) p(2,2)];
                [d_min, x_d_min, y_d_min] =   p_poly_dist(xp, yp, xv, yv);
                f_dir = ([cx cy]-[x_d_min, y_d_min])';
                check1 = [];
                
                if d_min<1000 && boundary == 0
                    worked1 = 1;
                    v1 = ([floe1.Ui floe1.Vi]+ floe1.ksi_ice*(pcenter(k,:)-[floe1.Xi floe1.Yi]));
                    v2 = ([floe2.Ui floe2.Vi]+ floe2.ksi_ice*(pcenter(k,:)-[floe2.Xi floe2.Yi]));
                    poly = r(k);
                    c0 = poly.Vertices;
                    [d_min1] = p_poly_dist(c0(:,1), c0(:,2),poly1.Vertices(:,1), poly1.Vertices(:,2));
                    check1 = zeros(length(c0),1);
                    check1(abs(d_min1)<1e-5) = 1;
                    [d_min2] = p_poly_dist(poly1.Vertices(:,1), poly1.Vertices(:,2),c0(:,1), c0(:,2));
                    check2 = zeros(length(poly1.Vertices),1);
                    check2(abs(d_min2)<1e-3) = 1;
                    if sum(check1)<length(check1)
                        while check1(1)+check1(end)>1.5 || check1(1)<0.5
                            check1 = [check1(2:end); check1(1)];
                            poly.Vertices = [poly.Vertices(2:end,:);poly.Vertices(1,:)];
                        end
                    elseif sum(check1) == length(check1) && boundary ==1
                        worked1 = 0; check1 = [];
                    else
                        [d_min1] = p_poly_dist(c0(:,1), c0(:,2),poly2.Vertices(:,1), poly2.Vertices(:,2));
                        check1 = zeros(length(c0),1);
                        check1(abs(d_min1)<1e-3) = 1;
                        if sum(check1)<length(check1) && sum(check1) > 1
                            while check1(1)+check1(end)>1.5 || check1(1)<0.5
                                check1 = [check1(2:end); check1(1)];
                                poly.Vertices = [poly.Vertices(2:end,:);poly.Vertices(1,:)];
                            end
                        else
%                             [d_min1] = p_poly_dist(P(1,:), P(2,:),c0(:,1)', c0(:,2)');
%                             check1 = zeros(length(c0),1);
%                             check1(abs(d_min1)<1e-3) = 1;
%                             while check1(1)+check1(2)<1.5
%                                 check1 = [check1(2:end); check1(1)];
%                                 P = [P(:,2:end) P(:,1)];
%                             end
                            worked1 = 0;
%                             xx =1; xx(1) =[1 2];
                        end
                    end
                    check1 = logical(check1);
                    check2 = logical(check2);
                    if worked1 == 1
                        verts = poly.Vertices(check1,:);
                        tangent = [verts(1,1)-verts(end,1); verts(1,2)-verts(end,1)];%[(p(1,2)-p(1,1)); (p(2,2)-p(2,1))];
                        f_dir = [-tangent(2); tangent(1)];
%                         xx = 1; xx(1) =[1 2];
%                         P = poly1.Vertices(check2,:);
                        
%                         [k,~] = dsearchn(P,[cx, cy]);
%                         x_d_min = P(k,1); y_d_min = P(k,2);
%                         xv = poly.Vertices(check1,1)'; yv =poly.Vertices(check1,2)';
%                         per = sqrt((xv(2:end)-xv(1:end-1)).^2+(yv(2:end)-yv(1:end-1)).^2)/sum(sqrt((xv(2:end)-xv(1:end-1)).^2+(yv(2:end)-yv(1:end-1)).^2));
%                         xc = (xv(2:end)+xv(1:end-1))/2; yc = (yv(2:end)+yv(1:end-1))/2; n = length(xc);
%                         if n<3
%                             [~,I2] = max((poly2.Vertices(:,1)-cx).^2+(poly2.Vertices(:,2)-cy).^2);
%                             xc(3) = poly2.Vertices(I2,1); yc(3) = poly2.Vertices(I2,2);
%                         end
%                         c = floe2.c;
%                         box = [min(c(1,:)) min(c(1,:)) max(c(1,:)) max(c(1,:)) min(c(1,:)); min(c(2,:)) max(c(2,:)) max(c(2,:)) min(c(2,:)) min(c(2,:))];
%                         [~, b,~,~,~] = polybnd_voronoi([xc; yc;]',[box(1,:); box(2,:)]');
%                         for ii = 1:length(b)
%                             polyb = polyshape(b{ii}); verts = b{ii};
%                             [in] = inpolygon(xc',yc',verts(:,1),verts(:,2));
%                             polynew = intersect(polyb,poly);
%                             [cxi, cyi] = centroid(polynew);
%                             per(ii) = area(polynew)/area(poly);
%                             if ~isnan(cxi)
%                                 F_dir(:,ii) = per(ii)*[xc(in)-cxi; yc(in)-cyi];
%                             else
%                                 F_dir(:,ii) =[0;0];
%                             end
%                         end
                        %f_dir = sum(F_dir,2);
%                         xx = 1; xx(1) =[1 2];
                        [~, x_d_min, y_d_min] =   p_poly_dist(xp, yp, xv, yv,0);                        
                    elseif worked1 == 0 && boundary == 0
                        xv = poly1.Vertices(:,1)'; yv =poly1.Vertices(:,2)';
                        [~, x_d_min, y_d_min] =   p_poly_dist(xp, yp, xv, yv);
                        f_dir = ([cx cy]-[x_d_min, y_d_min])';
%                         load F_dir
%                         xv = P(1,check1); yv = P(2,check1);
                    end
                    
%                     if abs (abs(f_dir2(1))-abs(f_dir(1)))/abs(f_dir(1))>0.001
%                         xx =1; xx(1) =[1 2];
%                     end
                end
%                 f_dir = [floe1.Ui-floe2.Ui; floe1.Vi-floe2.Vi];%(v1-v2)'; 

%                 alpha_i = 0;
%                 A_rot=[cos(alpha_i) -sin(alpha_i); sin(alpha_i) cos(alpha_i)]; %rotation matrix
%                 if abs(f_dir(1))<abs(f_dir(2))
%                     xx = 1; xx(1) =[1 2];
%                 end
                Px1=InterX([xv; yv],[cx floe1.Xi; cy floe1.Yi]);
                Px2=InterX([xv; yv],[cx floe2.Xi; cy floe2.Yi]);
%                 if ~isempty(Px1)
%                     f_dir = ([cx; cy]-Px1);
%                 elseif ~isempty(Px2)
%                     f_dir = ([cx; cy]-Px2);
%                 else
%                     f_dir = ([cx cy]-pcenter(k,:))';
%                 end
                floe_Rforce1=[cx cy]-[floe1.Xi floe1.Yi];
                floe_Rforce2=[cx cy]-[floe2.Xi floe2.Yi];
%                 f_dir = [floe_Rforce1(1)-floe_Rforce2(1); floe_Rforce1(2)-floe_Rforce2(2)];
%                 f_dir = ([cx cy]-[p(1,1) p(2,1)])';

                dist_eff=(min_dist_center2+min_dist_center1)/2;
                if dist_eff*Force_factor > 1e6 && polyA > 1e7 % ensure that forces are somewhat bounded
                    %display(dist_eff*Force_factor);
%                     xx = 1; xx(1) =[1 2];
                    dist_eff=1e6/Force_factor;
                elseif dist_eff*Force_factor > 1e5 && polyA <= 1e7
                    dist_eff=10;
                end
                
                if sum(check1)
                    f_dir=pp2-pp1; % if sorting worked use the local force criteria, if didn't work push the floe towards the furthest intersect point.
                end
%                 f_dir = sign(f_dir(1))*f_dir; %f_dir2 = sign(f_dir2(1))*f_dir2;
%                 CosTheta = max(min(dot(f_dir,f_dir2)/(norm(f_dir)*norm(f_dir2))),-1);
%                 ThetaInDegrees = real(acosd(CosTheta));
%                 if vecnorm(fdir(1).dir)==0
%                     force_dir=f_dir/sqrt(f_dir'*f_dir); %always push the flow away from contact point!
%                     fdir.dir = force_dir; fdir.loc =[cx cy];
%                 else
%                     worked = 0;
% %                     xx = 1; xx(1) =[1 2];
%                     for ii = 1:length(fdir)
%                         if sqrt((cx-fdir(ii).loc(1))^2+(cy-fdir(ii).loc(2))^2)<1
%                             force_dir = fdir(ii).dir; fdir(ii).loc = [cx cy];
%                             worked = 1;
%                         end
%                     end
%                     if worked == 0
%                         force_dir=f_dir/sqrt(f_dir'*f_dir); %always push the flow away from contact point!
%                         fdir(length(fdir)+1).dir = force_dir; fdir(length(fdir)).loc =[cx cy];
%                     end
%                 end
                force_dir=f_dir/sqrt(f_dir'*f_dir); %always push the flow away from contact point!
%                 save('F_dir.mat','fdir','f_dir')
                        
%                     save('F_dir.mat','force_dir');
%                 end
                
%                 poly1new = translate(poly1,force_dir');
%                 polyInew = intersect(poly1new,poly2);
%                 if area(polyInew)/area(polyI)-1 > 0
%                     force_dir = -force_dir;
%                 end
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
                save('F_dir.mat','force_dir')
%                 xx = 1; xx(1) =[1 2];
%                 load F_dir
                floe_Rforce1=[cx cy]-[floe1.Xi floe1.Yi];
                floe_torque1=cross([floe_Rforce1 0], [force_dir' 0]);
                floe_Rforce2=[cx cy]-[floe2.Xi floe2.Yi];
                floe_torque2=cross([floe_Rforce2 0], [-force_dir' 0]);
                if abs(floe_torque1(3)+floe_torque2(3)) > 1
%                     f_dir = ([cx cy]-[p(1,2) p(2,2)])';
%                     force_dir=f_dir/sqrt(f_dir'*f_dir); %always push the flow away from contact point!
%                     poly1new = translate(poly1,force_dir');
%                     polyInew = intersect(poly1new,poly2);
%                     rnew = regions(polyInew); Anew = intersect(r(k),rnew);
%                     [~,Imax] = max(area(Anew));
%                     if area(rnew(Imax))/area(r(k))-1 > 0
%                         force_dir = -force_dir;
%                     end
% %                     floe_Rforce1=[cx cy]-[floe1.Xi floe1.Yi];
% %                     floe_torque1=cross([floe_Rforce1 0], [force_dir' 0]);
% %                     floe_Rforce2=[cx cy]-[floe2.Xi floe2.Yi];
% %                     floe_torque2=cross([floe_Rforce2 0], [-force_dir' 0]);
%                     if sign(floe_torque1(3))*sign(floe_torque2(3)) > 0
%                         xx = 1; xx(1)=[1 2];
%                     end
                end
%                 [m,n] = size(f_hist_t); f_hist_t(:,n+1) = force_dir';
               % f_hist(:,n+1) = force_dir';
%                 save('fhist.mat','f_hist','f_hist_t')

%                 force=force_dir*(dl*dist_eff)*Force_factor; %proportional to the overlap area
%                 force=force_dir*area(r(k))*Force_factor/cos(CosTheta); %proportional to the overlap area
                force=force_dir*area(r(k))*Force_factor; %proportional to the overlap area
                
                v1 = ([floe1.Ui floe1.Vi]+ floe1.ksi_ice*(pcenter(k,:)-[floe1.Xi floe1.Yi]));
                v2 = ([floe2.Ui floe2.Vi]+ floe2.ksi_ice*(pcenter(k,:)-[floe2.Xi floe2.Yi]));
                dir_t = [-force_dir(2) force_dir(1)];
                v_t = dot(v1-v2,dir_t);
                v_r = (v1+v2)/2;
%                 force_t = dl*G*v_t*sign(v_r-v1).*abs(dir_t);
%                 force_t = dl*G*v_t*dir_t;
                force_t = [0 0];%-area(r(k))*(v1-v2)*Force_factor/4;%[0 0];%
%                 vnew = [floe1.Ui floe1.Vi]+force_t/norm(force_t)*(0.1*norm([floe1.Ui floe1.Vi]));
%                 if sum(vnew.*vnew)>floe1.Ui^2+floe1.Vi^2
%                     dir_t = -dir_t;
%                     force_t = -force_t;%dl*G*v_t*dir_t;
%                     vnew2= [floe1.Ui floe1.Vi]+force_t/norm(force_t)*(0.1*norm([floe1.Ui floe1.Vi]));
%                     if sum(vnew2.*vnew2)>floe1.Ui^2+floe1.Vi^2
%                         xx = 1; xx(1) =[1 2];
%                     end
%                 end
%                 else
%                     force_t = dl*G*v_t*sign(v2).*abs(dir_t);
%                 end
                if vecnorm(force_t)>mu*vecnorm(force)
%                     if abs(v1)>0
                        force_t = mu*vecnorm(force)*sign(v_r-v1).*abs(dir_t);
%                     else
%                         force_t = mu*vecnorm(force)*sign(v2).*abs(dir_t);
%                     end
                end
%                 xx = 1;
%                 xx(1) = [1 2];

                if isnan(force(1)) || isnan(force(2))
                    xx = 1;
                    xx(1) = [1 2];
                end
                    
                force_1(k,:)= force'+force_t;
                overlap(k) = area(r(k));
                
            end
            
        else
            force_1(k,:) = [0 0];
            overlap(k) = 0;
        end

        
    end
%     if N_contact > 1
%         xx = 1; xx(1) =[1 2];
%     end
end


warning('on',id)
warning('on',id3)

end

