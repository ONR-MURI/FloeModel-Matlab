function [ force_1, pcenter, worked] = floe_interactions_bpm2(c1, c2)
%% 
% Find Vertices of overlapping polygons
Force_factor=500;%1e3;
polyout = intersect(c1,c2);
% if polyout.NumRegions>0
%     polyout2 = sortregions(polyout,'area','descend');
%     R = regions(polyout2);
%     R = R(area(R)>100);
%     if ~isempty(area(R)); polyout = union(R); areaPoly = area(polyout); else; areaPoly = 0; end
% else
%     areaPoly = 0;
% end
areaPoly = area(polyout);
if areaPoly/area(c1) >0.75 && area(c1) > 5e7
    X = 1;
    X(1) = [1 2];
end
worked=1;
%% 

%Calculate forces on different overlapping areas
if areaPoly==0
    force_1=[0 0];
    pcenter=[0 0];   
else    
    N_contact = polyout.NumRegions;
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
    node1 = c1.Vertices;
    edge1 = [1:length(node1); 2:length(node1) 1]';
    [check1,~] = inpoly2(c0',node1,edge1) ;
    node2 = c1.Vertices;
    edge2 = [1:length(node2); 2:length(node2) 1]';
    [check2,~] = inpoly2(c0',node2,edge2) ;
%     check1 = inpolygon(c0(:,1),c0(:,2),c2.Vertices(:,1),c2.Vertices(:,2));
%     check2 = inpolygon(c0(:,1),c0(:,2),c1.Vertices(:,1),c1.Vertices(:,2));
    force_1=zeros(N_contact,2);
    [pcenter(:,1),pcenter(:,2)] = centroid(polyout, 1:polyout.NumRegions);
    
    
    for k=1:N_contact       
        %find vertices of points inside poly1
        in2 = c0(:,I(k)+1:I(k+1)-1)';
%         in2 = R(k).Vertices;
        if sum(check1(I(k)+1:I(k+1)-1))>sum(check2(I(k)+1:I(k+1)-1))
            [~,d] = dsearchn(c2.Vertices,in2);
            if length(d(d>0))>2
                in1 = in2(d>0,:);
                [~,d2] = dsearchn(c1.Vertices,in1);
                if max(d2) > 0 && min(abs(d2))==0
                    while abs(d2(end))==0 || abs(d2(end-1))> 0
                        in1 = [in1(2:end,:);in1(1,:)];
                        d2 = [d2(2:end);d2(1)];
                    end
                elseif max(d2) > 0 && min(abs(d2)) > 0
                    while abs(d2(end))==0
                        in1 = [in1(2:end,:);in1(1,:)];
                        d2 = [d2(2:end);d2(1)];
                    end
                end
            else
                in1 = in2;
            end
        else
            [~,d] = dsearchn(c1.Vertices,in2);
            if length(d(d>0))>2
                in1 = in2(d>0,:);
                [~,d2] = dsearchn(c2.Vertices,in1);
                if max(d2) > 0 && min(abs(d2))==0
                    while abs(d2(end))==0 || abs(d2(end-1))> 0
                        in1 = [in1(2:end,:);in1(1,:)];
                        d2 = [d2(2:end);d2(1)];
                    end
                elseif max(d2) > 0 && min(abs(d2)) > 0
                    while abs(d2(end))==0 
                        in1 = [in1(2:end,:);in1(1,:)];
                        d2 = [d2(2:end);d2(1)];
                    end
                end
            else
                in1 = in2;    
            end
        end
        
        
        %calculate direction of force by calculating direction from center
        %of each interior edge pointing to center of exterior edges of
        %overlapping area
        normal = [in1(2:end,2)-in1(1:end-1,2) in1(1:end-1,1)-in1(2:end,1)];
        dir = [normal(:,1)./vecnorm(normal,2,2) normal(:,2)./vecnorm(normal,2,2)];
        checkpoint = [0.5*(in1(1:end-1,1)+in1(2:end,1)) 0.5*(in1(1:end-1,2)+in1(2:end,2))] + 1e-4*dir;
        [check] = inpolygon(checkpoint(:,1),checkpoint(:,2),in2(:,1),in2(:,2));
        dir(check == 0,:) = -dir(check == 0,:); 
        dl = sqrt((in1(1:end-1,1)-in1(2:end,1)).^2 + (in1(1:end-1,2)-in1(2:end,2)).^2);
        f_dir = -sum(dir.*[dl dl]);          
        force_dir=f_dir/norm(f_dir);
        
        force_1(k,:)=force_dir*area(polyout,k)*Force_factor;        
    end  
end
%% 

if max(max(abs(force_1)))>0
    displace = 0.1*sum(force_1,1)/norm(sum(force_1,1));
    c1.Vertices(:,1) = c1.Vertices(:,1) + displace(1);
    c1.Vertices(:,2) = c1.Vertices(:,2)+displace(2);
    polynew = intersect(c1,c2);
    if area(polynew)>areaPoly
        force_1 = -force_1;
        worked = 0;
    end
end
end