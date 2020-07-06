function [ force_1, pcenter, worked,live] = floe_interactions_poly(c1, c2)
%% This function calculates the forces between two interaction floes

% Find overlapping areas between polygons
live = [1 1];
Force_factor=1.1e3;
polyout = intersect(c1,c2);
if polyout.NumRegions>0
    polyout2 = sortregions(polyout,'area','descend');
    R = regions(polyout2);
    R = R(area(R)>10);
    if ~isempty(area(R));  areaPoly = max(area(R)); else; areaPoly = 0; end
else
    areaPoly = 0;
end

worked=1;

%Make sure one polygon is not entirely within the other
if areaPoly/area(c1)>0.99 && areaPoly/area(c2)>0.99
    x = 1;
    x(1) = [1 2];
elseif areaPoly/area(c1)>0.5 && area(c1)/area(c2)<0.05
    force_1=[0 0];
    pcenter=[0 0];  
    areaPoly = 0;
    live(1) = 0;
elseif areaPoly/area(c2)>0.5 && area(c2)/area(c2)>0.05
    force_1=[0 0];
    pcenter=[0 0];
    areaPoly = 0;
    live(2) = 0;
end

%If overlapping area is to small then set forces to zero
if areaPoly<10
    force_1=[0 0];
    pcenter=[0 0];   
elseif areaPoly/area(c1)>0.75 && c2.NumRegions < 1
    force_1=[0 0];
    pcenter=[0 0]; 
elseif areaPoly/area(c2)>0.75 
    force_1=[0 0];
    pcenter=[0 0]; 
else
    %Idenfity the number of overlapping regions between the polygons
    N_contact = length(R);
    force_1=zeros(N_contact,2);
    pcenter = zeros(N_contact,2);
    
    %Calculate the force from each of these overlaps
    for k=1:N_contact      
        poly = rmholes(R(k));
        c0 = poly.Vertices;
        [pcenter(k,1),pcenter(k,2)] = centroid(R(k));

        %Idenfity the vertices of overlapping polygons
        [d_min1] = p_poly_dist(c0(:,1), c0(:,2), c1.Vertices(:,1), c1.Vertices(:,2));
        [d_min2] = p_poly_dist(c0(:,1), c0(:,2), c2.Vertices(:,1), c2.Vertices(:,2));        
        check1 = zeros(length(c0),1); check2 = check1;
        check1(abs(d_min1)<1e-3) = 1;
        check2(abs(d_min2)<1e-3) = 1;
        [checkin1] = inpolygon(c0(:,1),c0(:,2),c1.Vertices(:,1),c1.Vertices(:,2));
        [checkin2] = inpolygon(c0(:,1),c0(:,2),c2.Vertices(:,1),c2.Vertices(:,2));
        checkin1(checkin1-check1<1) = 0;
        checkin2(checkin2-check2<1) = 0;
        in2 = c0;
        if sum(checkin1)>sum(checkin2)
            in = logical(checkin1);
        else
            in = logical(checkin2);
        end
        
        %Check to see if there are any holes in the polyshape that need to
        %be dealt with
        if c2.NumHoles>0
            holes = isnan(c2.Vertices(:,1));
            II = find(holes==1);
            [d_min2] = p_poly_dist(c0(:,1), c0(:,2), c2.Vertices(II+1:end,1), c2.Vertices(II+1:end,2));
            check2(abs(d_min2)<1e-3) = 1;
            [checkin2] = inpolygon(c0(:,1),c0(:,2),c2.Vertices(:,1),c2.Vertices(:,2));
            checkin2(checkin2-check2<1) = 0;
            in = logical(checkin2);
            %Force_factor=1.5e3;
        end
        
        %Identify points that are on edges versus in the polygons
        on = logical(abs(in-1));
        ii = 1;
        if sum(in) < 0.5 || sum(on) < 1.5
            force_1(k,:)= [0 0];
            overlap = false;
        else
            while in(1) > 0.5 || in(1)+in(2)< 0.5
                in = [in(2:end); in(1)];
                c0 = [c0(2:end,:); c0(1,:)];
                if ii == length(in)+1
                    xx = 1;
                    xx(1) = [1 2];
                end
                ii = ii+1;
            end
            overlap = true;

        end
        
        if overlap
            numReg = 0;
            clear BreakStart; clear BreakEnd
            for ii=2:length(in)
                if in(ii) > in(ii-1)
                    numReg = numReg + 1;
                    BreakStart(numReg) = ii;
                elseif in(ii)<in(ii-1)
                    BreakEnd(numReg) = ii-1;
                end
            end
            if in(end)>in(1)
                BreakEnd(numReg) = length(in);
            end
        
            for jj = 1:numReg
                if BreakEnd(jj) == length(in)
                    in1 = [c0(BreakStart(jj)-1:BreakEnd(jj),:); c0(1,:)];
                else
                    in1 = c0(BreakStart(jj)-1:BreakEnd(jj)+1,:);
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
                f_dir(jj,:) = -sum(dir.*[dl dl]);
            end
            
            %Find direction from all regions that contribute to this one
            f_dir = sum(f_dir,1);
            if norm(f_dir) < 0.001
                force_dir = 0;
            else
                force_dir=f_dir/norm(f_dir);
            end
            
            %Check to make sure force will actually push the two floes
            %apart
            if max(max(abs(force_dir)))>0
                displace = 5*force_dir;
                C1(:,1) = c1.Vertices(:,1) + displace(1);
                C1(:,2) = c1.Vertices(:,2)+displace(2);
                C2(:,1) = c1.Vertices(:,1) - displace(1);
                C2(:,2) = c1.Vertices(:,2)-displace(2);
                polynew1 = intersect(polyshape(C1),R(k));
                polynew2 = intersect(polyshape(C2),R(k));
                if area(polynew1)>area(polynew2)
                    force_dir = -force_dir;
                    worked = 0;
                end
                
            end
            
            force_1(k,:)=force_dir*area(R(k))*Force_factor;
        end
    end  
end
%% 

if isnan(force_1)
    xx = 1;
    xx(1) = [1 2];
end

end