function [ force_1, pcenter, overlap] = floe_interactions(c1, c2)
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)

Force_factor=1e3; overlap = 0;
polyA = area(polyshape(c1'));
% c1=[Floe1.c_alpha(1,:)+Floe1.Xi; Floe1.c_alpha(2,:)+Floe1.Yi];
%
% if ~isa(Floe2,'struct')
%     c2=Floe2; % interaction with a boundary contour
% else
%     c2=[Floe2.c_alpha(1,:)+Floe2.Xi; Floe2.c_alpha(2,:)+Floe2.Yi];
% end
if norm(c1(:,1)-c1(:,end))> 1
    c1(:,length(c1)+1) = c1(:,1);
end
if norm(c2(:,1)-c2(:,end))> 1
    c2(:,length(c2)+1) = c2(:,1);
end

P=InterX(c1,c2);
% [xint,yint] = polyxpoly(c1(1,:),c1(2,:),c2(1,:),c2(2,:),'unique');
% P = [xint'; yint'];

if isempty(P) || size(P,2)<2
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
                
                if dist_eff*Force_factor > 1e5 && polyA > 1e7 % ensure that forces are somewhat bounded
                    %display(dist_eff*Force_factor);
                    dist_eff=100;
                elseif dist_eff*Force_factor > 1e5 && polyA <= 1e7
                    dist_eff=10;
                end
                
                
                f_dir=pp2-pp1; % if sorting worked use the local force criteria, if didn't work push the floe towards the furthest intersect point.
                
                force_dir=f_dir/sqrt(f_dir'*f_dir); %always push the flow away from contact point!
                
                
                force=force_dir*(dl*dist_eff)*Force_factor; %proportional to the overlap area
                
                
                if isnan(force(1)) || isnan(force(2))
                    xx = 1;
                    xx(1) = [1 2];
                end
                
                force_1(k,:)= force;
                overlap(k) = dl*dist_eff;
                
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

