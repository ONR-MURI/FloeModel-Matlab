function [ force_1, pcenter, overlap] = floe_interactions(floe1, floe2, c2_boundary,PERIODIC)
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)

Modulus = 1e7; h1 = floe1.h; h2 = floe2.h;
r1 = sqrt(floe1.area); r2 = sqrt(floe2.area);
Force_factor=Modulus*(h1*h2)/(h1*r2+h2*r1); overlap = 0;
G = 3.8e9; mu = 0.45;
c1=[floe1.c_alpha(1,:)+floe1.Xi; floe1.c_alpha(2,:)+floe1.Yi];
if isfield(floe2,'c')
    c2=floe2.c;
    poly2 = polyshape(c2');
    border = 0;
else
    c2 = c2_boundary;
    poly2 = floe2.poly;
    border = 1;
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
if  ( max(c1(1,:))<max(c2_boundary(1,:)) && min(c1(1,:))>min(c2_boundary(1,:)) && max(c1(2,:))<max(c2_boundary(2,:)) && min(c1(2,:))>min(c2_boundary(2,:)) || PERIODIC)
    
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
% if border && max(abs(P(1,:))) == max(abs(c2_boundary(1,:))) && max(abs(P(2,:))) == max(abs(c2_boundary(2,:)))
% %     xx = 1; xx(1) =[1 2];
%     P(:,end+1) = [sign(floe1.Xi)*max(abs(c2_boundary(1,:))); sign(floe1.Yi)*max(abs(c2_boundary(2,:)))];
% end
% [xint,yint] = polyxpoly(c1(1,:),c1(2,:),c2(1,:),c2(2,:),'unique');
% P = [xint'; yint'];

if isempty(P) || size(P,2)<2 || isinf(overlap)
    force_1=[0 0];
    pcenter=[0 0];
else
    
    [P, worked] = sort_contact_points(P, c1, c2 );
    
    N_contact=size(P,2)-1;
    
    force_1=zeros(N_contact,2);
    
    pcenter=zeros(N_contact,2);
    
    same_dir=zeros(1,N_contact);
    
    for k=1:N_contact
        
        p=P(1:2,k:k+1);
        
        dl=sqrt((p(1,1)-p(1,2))^2+(p(2,1)-p(2,2))^2);
        
        pcenter(k,:)=mean(p,2); dp=p(:,1)-p(:,2);
        
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
                
                dist_eff=(min_dist_center2+min_dist_center1)/2;
                
                if dist_eff*Force_factor > 1e6 && polyA > 1e7 % ensure that forces are somewhat bounded
                    %display(dist_eff*Force_factor);
                    %                     xx = 1; xx(1) =[1 2];
                    dist_eff=1e6/Force_factor;
                elseif dist_eff*Force_factor > 1e5 && polyA <= 1e7
                    dist_eff=10;
                end
                
                
                f_dir=pp2-pp1; % if sorting worked use the local force criteria, if didn't work push the floe towards the furthest intersect point.
                
                force_dir=f_dir/sqrt(f_dir'*f_dir); %always push the flow away from contact point!
                
                poly1new = translate(poly1,force_dir');
                polyInew = intersect(poly1new,poly2);
                if area(polyInew)/area(polyI)-1 > 0 %&& ~border
                    force_dir = -force_dir;
%                 elseif border
%                     force_dir = -sign(pcenter(k,:)).*abs(force_dir);
                end
                
                force=force_dir*(dl*dist_eff)*Force_factor; %proportional to the overlap area
                
                v1 = ([floe1.Ui floe1.Vi]+ floe1.ksi_ice*(pcenter(k,:)-[floe1.Xi floe1.Yi]));
                v2 = ([floe2.Ui floe2.Vi]+ floe2.ksi_ice*(pcenter(k,:)-[floe2.Xi floe2.Yi]));
                dir_t = [-force_dir(2) force_dir(1)];
                v_t = dot(v1-v2,dir_t);
                v_r = (v1+v2)/2;
                %                 if abs(v1)>0
                force_t = dl*G*v_t*sign(v_r-v1).*abs(dir_t);
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
                overlap(k) = area(polyI);
                
            end
            
        else
            force_1(k,:) = [0 0];
            overlap(k) = 0;
        end
        
        
    end
    
end

warning('on',id)
warning('on',id3)

end

