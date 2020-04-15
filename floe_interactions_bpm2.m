function [ force_1, pcenter, worked] = floe_interactions_bpm2(c1, c2)
%% 
% Find Vertices of overlapping polygons
Force_factor=5e2;
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
%% 

%Calculate forces on different overlapping areas
if areaPoly<10
    force_1=[0 0];
    pcenter=[0 0];   
elseif areaPoly/area(c1)>0.75 || areaPoly/area(c2)>0.75
    force_1=[0 0];
    pcenter=[0 0]; 
else
    %% 
    
    N_contact = length(R);
    force_1=zeros(N_contact,2);
    pcenter = zeros(N_contact,2);
    
    %% 
    
    for k=1:N_contact       
        c0 = R(k).Vertices;
        [pcenter(k,1),pcenter(k,2)] = centroid(R(k));

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
        
        ii = 1;
        if sum(in) < 0.5
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
        
        if overlap
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
            f_dir = sum(f_dir,1);
            force_dir=f_dir/norm(f_dir);
                        
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