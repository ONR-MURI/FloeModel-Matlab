function [ force_1, pcenter, worked] = floe_interactions_poly(c1, c2)

Force_factor=1e3;

% c1=[Floe1.c_alpha(1,:)+Floe1.Xi; Floe1.c_alpha(2,:)+Floe1.Yi];
% 
% if ~isa(Floe2,'struct')
%     c2=Floe2; % interaction with a boundary contour
% else
%     c2=[Floe2.c_alpha(1,:)+Floe2.Xi; Floe2.c_alpha(2,:)+Floe2.Yi];
% end

c1=c1.Vertices'; c1=[c1 c1(:,1)];
c2=c2.Vertices'; c2=[c2 c2(:,1)]; %adding the first point to the end to make sure lines are closed! polygons do not include the last point!

P=InterX(c1,c2);


if isempty(P) || size(P,2)<2
    force_1=[0 0];
    pcenter=[0 0];
    worked=1;
else
    
    
    [P, worked] = sort_contact_points(P, c1, c2 );
    
    N_contact=size(P,2)/2;
    
    force_1=zeros(N_contact,2);

    pcenter=zeros(N_contact,2);
    
    same_dir=zeros(1,N_contact);
    
    for k=1:N_contact
        
        p=P(1:2,2*k-1:2*k);
        
        dl=sqrt((p(1,1)-p(1,2))^2+(p(2,1)-p(2,2))^2);
        
        pcenter(k,:)=mean(p,2); dp=p(:,1)-p(:,2);
        
        d=1e5; %checking intersections with countour at +-d distance from the center point.
        
        dirr=dp(1)/dp(2);
        if abs(dirr)>1e3
            
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
        
        pp1=InterX(c1,[x1 x2; y1 y2]);
        %[min_dist_center1, i_min] = min( sqrt( ( (pp1(1,:)-pcenter(k,1)).^2 + (pp1(2,:)-pcenter(k,2)).^2)));
        [sorted_dist_center1, i_min] = sort( sqrt( ( (pp1(1,:)-pcenter(k,1)).^2 + (pp1(2,:)-pcenter(k,2)).^2)));
        min_dist_center1=sorted_dist_center1(1);
        
        f_dir=pp1(:,i_min(end))-pp1(:,i_min(1)); %from point with min distance to the furthest point.
        pp1=pp1(:,i_min(1));
        
        pp2=InterX(c2,[x1 x2; y1 y2]);
        %[min_dist_center2, i_min] = min(sqrt( ( (pp2(1,:)-pcenter(k,1)).^2 + (pp2(2,:)-pcenter(k,2)).^2)));
        [sorted_dist_center2, i_min] = sort(sqrt( ( (pp2(1,:)-pcenter(k,1)).^2 + (pp2(2,:)-pcenter(k,2)).^2)));
        min_dist_center2=sorted_dist_center2(1);
        pp2=pp2(:,i_min(1));
        
        same_dir(k)=((pp2-pp1)'*f_dir)>0;
        
        dist_eff=(min_dist_center2+min_dist_center1)/2;
        
%         dist_center2center1=sqrt((Floe1.Xi-pcenter(k,1))^2+(Floe1.Yi-pcenter(k,2))^2);
%         
%         if ~isa(Floe2,'struct'), 
%             dist_ratio=dist_eff/dist_center2center1;
%         else
%             dist_center2center2=sqrt((Floe2.Xi-pcenter(k,1))^2+(Floe2.Yi-pcenter(k,2))^2);
%             dist_ratio=max(dist_eff/dist_center2center1,dist_eff/dist_center2center2);
%         end
%         
%         
%         if dist_ratio > 1, %increase forces for floes that overlap too much.
%             %display(dist_ratio),
%             Force_factor=Force_factor*dist_ratio;
%         end
%         
        
         if dist_eff*Force_factor > 1e5, % ensure that forces are somewhat bounded
             %display(dist_eff*Force_factor);
             dist_eff=100;
         end
        
        if worked, f_dir=pp2-pp1; end % if sorting worked use the local force criteria, if didn't work push the floe towards the furthest intersect point.
        force_dir=f_dir/sqrt(f_dir'*f_dir); %always push the flow away from contact point!
        
        force=force_dir*(dl*dist_eff)*Force_factor; %proportional to the overlap area
        
        force_1(k,:)= force;
%       force_2(k,:)=-force;
%        
%         r1=[Floe1.Xi; Floe1.Yi];
%         if ~isa(Floe2,'struct'); r2=NaN; else r2=[Floe2.Xi; Floe2.Yi]; end
%         
%         torque=cross([pcenter(k,:)-r1' 0], [force_1(k,:) 0]);torque_1(k)=torque(3);
%         torque=cross([pcenter(k,:)-r2' 0], [force_2(k,:) 0]);torque_2(k)=torque(3);
%         
%                         figure;
%                         plot(c1(1,:),c1(2,:),'r'); hold on;
%                         plot(c2(1,:),c2(2,:),'b');
%                         plot(P(1,2*k-1:2*k),P(2,2*k-1:2*k),'x');
%                         quiver(pcenter(k,1),pcenter(k,2),force_1(k,1),force_1(k,2),1/Force_factor,'r')
%                         quiver(pcenter(k,1),pcenter(k,2),force_2(k,1),force_2(k,2),1/Force_factor,'b')
%         
%                         plot(Floe1.Xi,Floe1.Yi,'r.');
%                         if isa(Floe2,'struct'); plot(Floe2.Xi,Floe2.Yi,'b.'); end
%                         plot(pp1(1),pp1(2),'ro');
%                         plot(pp2(1),pp2(2),'bo');
%                         plot(pcenter(k,1),pcenter(k,2),'*')
        
    end
    
    
    
end


end

