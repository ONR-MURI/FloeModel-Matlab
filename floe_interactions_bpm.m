function [ force_1, pcenter, worked] = floe_interactions_bpm(c1, c2)
%% 
% Find Vertices of overlapping polygons
Force_factor=1e3;
polyout = intersect(c1,c2);
areaPoly = area(polyout);
poly1 = c1;
c1=c1.Vertices'; c1=[c1 c1(:,1)];
c0(:,1) = polyout.Vertices(:,1); c0(:,2) = polyout.Vertices(:,2);
c0 = c0';

%Set up arrays such that they can be read for different overlapping areas
breaks = isnan(polyout.Vertices(:,1));
I = find(breaks == 1);
[m,n] = size(I);
if m  == 1
    I = [0 I length(breaks)+1];
elseif n == 1
    I = [0 I' length(breaks)+1];
end
worked=1;
%% 

%Calculate forces on different overlapping areas
if areaPoly==0
    force_1=[0 0];
    pcenter=[0 0];   
else    
    N_contact = polyout.NumRegions;    
    force_1=zeros(N_contact,2);
    [pcenter(:,1),pcenter(:,2)] = centroid(polyout, [1:polyout.NumRegions]);
    
    for k=1:N_contact
        
        %find vertices of points inside poly1
        in1 = c0(:,I(k)+1:I(k+1)-1)';
        [~,d] = dsearchn(c1',in1);
        
        %calculate direction of force by calculating direction from center
        %of each interior edge pointing to center of exterior edges of
        %overlapping area
        %EdgeCenter = mean(in1(d==0,:));
        if max(max(in1)) == 1
            in1 = in1(d>0,:);
        end
        for ii = 1:length(in1)-1
            %dir(ii,:) = EdgeCenter-mean(in1(ii:ii+1,:));
            normal = [(in1(ii+1,2)-in1(ii,2)) (in1(ii,1)-in1(ii+1,1))];
            dir(ii,:) = normal/norm(normal);
            checkpoint = mean(in1(ii:ii+1,:))+dir(ii,:);
            [check] =inpolygon(checkpoint(1),checkpoint(2),in1(:,1),in1(:,2));
            %dir(ii,:) = EdgeCenter-mean(in1(ii:ii+1,:));
            if check == 0
                dir(ii,:) = -dir(ii,:);    
            end
            dl(ii) = sqrt((in1(ii,1)-in1(ii+1,1))^2+(in1(ii,2)-in1(ii+1,2))^2);
        end 
        f_dir = -sum(dir.*[dl' dl']);          
        force_dir=f_dir/norm(f_dir);
        
        force_1(k,:)=force_dir*area(polyout,k)*Force_factor;        
    end  
end
%% 

if max(max(force_1))>0
    displace = 10*sum(force_1,1)/norm(sum(force_1,1));
    poly1.Vertices(:,1) = poly1.Vertices(:,1) + displace(1);
    poly1.Vertices(:,2) = poly1.Vertices(:,2)+displace(2);
    polynew = intersect(poly1,c2);
    if area(polynew)>areaPoly
        force_1 = -force_1;
        worked = 0;
    end
end

end