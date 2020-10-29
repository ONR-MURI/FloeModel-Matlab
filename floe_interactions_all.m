function Floe = floe_interactions_all(Floe, ocean, winds,c2_boundary, dt, PERIODIC, RIDGING)

id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)

Lx= max(c2_boundary(1,:));
Ly= max(c2_boundary(2,:));
height.mean = 1; height.delta = 0;
c2_boundary_poly = polyshape(c2_boundary');
FloeB = initialize_floe_values(c2_boundary_poly, height);
%Floe(i).interactions=[floeNumber Fx Fy px py torque];
%Floe(i).potentialInteractions(j).floeNum
%Floe(i).potentialInteractions(j).c_alpha=Floe(floeNum).c_alpha.

N0=length(Floe);
if PERIODIC
    
    ghostFloeX=[];
    ghostFloeY=[];
    parent=[];
    translation = [];
    
    x=cat(1,Floe.Xi);
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        poly = polyshape(Floe(i).c_alpha'+[x(i) y(i)]);
        %   if alive(i) && (x(i)>Lx-rmax(i)) || (x(i)<-Lx+rmax(i))
        if alive(i) && (max(abs(poly.Vertices(:,1)))>Lx)
            
            ghostFloeX=[ghostFloeX  Floe(i)];
            ghostFloeX(end).Xi=Floe(i).Xi-2*Lx*sign(x(i));
            parent=[parent  i];
            translation = [translation; -2*Lx*sign(x(i)) 0];
            
        end
        
        
    end
    
    Floe=[Floe ghostFloeX];
    
    x=cat(1,Floe.Xi);
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        
        %   if alive(i) && (x(i)>Lx-rmax(i)) || (x(i)<-Lx+rmax(i))
        if alive(i) && (max(abs(poly.Vertices(:,2)))>Ly)
            
            ghostFloeY=[ghostFloeY  Floe(i)];
            ghostFloeY(end).Yi=Floe(i).Yi-2*Ly*sign(y(i));
            parent=[parent  i];
            translation = [translation; 0 -2*Ly*sign(y(i))];
            
        end
        
    end
    
    Floe=[Floe ghostFloeY];
    
end

%Find length of new Floe variable including the ghost floes
N=length(Floe);
x=cat(1,Floe.Xi);
y=cat(1,Floe.Yi);
rmax=cat(1,Floe.rmax);
alive=cat(1,Floe.alive);

for i=1:N  %do interactions with boundary in a separate parfor loop
    
    Floe(i).interactions=[];
    
    Floe(i).OverlapArea = 0;
    
    Floe(i).potentialInteractions=[];
    
    Floe(i).collision_force=[0 0];
    
    Floe(i).collision_torque=0;
    
    k=1;
    
    if ( alive(i) && ~isnan(x(i)) )
        for j=1:N
            %display(j);
            if j>i && alive(j) && sqrt((x(i)-x(j))^2 + (y(i)-y(j))^2)<(rmax(i)+rmax(j)) % if floes are potentially overlapping
                Floe(i).potentialInteractions(k).floeNum=j;
                Floe(i).potentialInteractions(k).c=[Floe(j).c_alpha(1,:)+x(j); Floe(j).c_alpha(2,:)+y(j)];
                k=k+1;
            end
            
        end
        
    end
    
end


parfor i=1:N  %now the interactions could be calculated in a parfor loop!
        
    c1=[Floe(i).c_alpha(1,:)+x(i); Floe(i).c_alpha(2,:)+y(i)];
    
    if ~isempty(Floe(i).potentialInteractions)
        
        for k=1:length(Floe(i).potentialInteractions)
            
            floeNum=Floe(i).potentialInteractions(k).floeNum;
            
            c2=Floe(i).potentialInteractions(k).c;
            
            [force_j,P_j, overlap] = floe_interactions(c1,c2);
            
            %if ~worked, disp(['contact points issue for (' num2str(i) ',' num2str(floeNum) ')' ]); end
            
            if sum(abs(force_j))~=0
                Floe(i).interactions=[Floe(i).interactions ; floeNum*ones(size(force_j,1),1) force_j P_j zeros(size(force_j,1),1) overlap'];
                Floe(i).OverlapArea = sum(overlap)+Floe(i).OverlapArea;
            end
            
        end
        
    end
    
    if ~PERIODIC
        [force_b, P_j, overlap] = floe_interactions(c1, c2_boundary);
        in = inpolygon(x(i),y(i),c2_boundary(1,:)',c2_boundary(2,:)');
        if ~in
            xx = 1;
            xx(1) = [1 2];
        end
        
%     if ~worked, display(['contact points issue for (' num2str(i) ', boundary)' ]); end
        if sum(abs(force_b))~=0,
            force_b = [-sign(x(i)) -sign(y(i))].*abs(force_b);
            % boundary will be recorded as floe number Inf;
            Floe(i).interactions=[Floe(i).interactions ; Inf*ones(size(force_b,1),1) force_b P_j zeros(size(force_b,1),2)];
            Floe(i).OverlapArea = sum(overlap)+Floe(i).OverlapArea;
            Floe(i).potentialInteractions(end+1).floeNum = Inf;
            Floe(i).potentialInteractions(end).c = c2_boundary;
        end
    end
    
end



%Fill the lower part of the interacton matrix (floe_i,floe_j) for floes with j<i
for i=1:N %this has to be done sequentially
      
    if ~isempty(Floe(i).interactions)
        
        a=Floe(i).interactions;
        
        indx=a(:,1);
        
        for j=1:length(indx)
            
            if indx(j)<=N && indx(j)>i
                Floe(indx(j)).interactions=[Floe(indx(j)).interactions; i -a(j,2:3) a(j,4:5) 0 a(j,7)];   % 0 is torque here that is to be calculated below
                Floe(indx(j)).OverlapArea = Floe(indx(j)).OverlapArea + a(j,7);
            end
            
        end
        
    end
    
end

% calculate all torques from forces
if PERIODIC
    
    parfor i=N0+1:N
        
        if ~isempty(Floe(i).interactions)
            
            a=Floe(i).interactions;
            r=[x(i) y(i)];
            for k=1:size(a,1)
                floe_Rforce=a(k,4:5);
                floe_force=a(k,2:3);
                floe_torque=cross([floe_Rforce-r 0], [floe_force 0]);
                Floe(i).interactions(k,6)=floe_torque(3);
            end
            
            Floe(i).collision_force=sum(Floe(i).interactions(:,2:3),1);
            Floe(i).collision_torque=sum(Floe(i).interactions(:,6),1);
            
        end
        
    end
     %add forces and torques from ghost floes to their parents; ghost floes
    %begin with the index N0+1
    for i=1:length(parent)
        Floe(parent(i)).collision_force =Floe(parent(i)).collision_force +Floe(N0+i).collision_force;
        Floe(parent(i)).collision_torque=Floe(parent(i)).collision_torque+Floe(N0+i).collision_torque;
    end
end

for i=1:N0
    
    if ~isempty(Floe(i).interactions)
        
       a=Floe(i).interactions;
       r=[x(i) y(i)];
        for k=1:size(a,1)
            floe_Rforce=a(k,4:5);
            floe_force=a(k,2:3);
            floe_torque=cross([floe_Rforce-r 0], [floe_force 0]);
            Floe(i).interactions(k,6)=floe_torque(3);
        end
        
       Floe(i).collision_force=sum(Floe(i).interactions(:,2:3),1)+Floe(i).collision_force;
       Floe(i).collision_torque=sum(Floe(i).interactions(:,6),1)+Floe(i).collision_torque;
        
    end
    
   %Do the timestepping now that forces and torques are known.
    if Floe(i).alive,
        tmp=calc_trajectory(dt,ocean, winds,Floe(i)); % calculate trajectory
        if (isempty(tmp) || isnan(x(i)) ), Floe(i).alive=0; else Floe(i)=tmp; end
    end
        
end

%timer1 = tic;
floenew = [];
if RIDGING
    %Create a function to control probability that ridging will occur
    overlapArea=cat(1,Floe.OverlapArea)./cat(1,Floe.area);
    keep=rand(length(Floe),1)>overlapArea;
%     keep=rand(length(Floe),1)<0.5;
    for i=1:N0

        if Floe(i).alive && ~keep(i)
            c2 = 0;
            if ~isempty(Floe(i).potentialInteractions)
                c = cat(1,Floe(i).potentialInteractions.floeNum);
                if isinf(max(c))
                    c2 = 1;
                end
                c = c(c<Inf);
                for ii = 1:length(c)
                    [Floe1, Floe2] = ridging(Floe(i),Floe(c(ii)),c2_boundary_poly,true);
                    if length(Floe1) > 1
                        Floe(i) = Floe1(1);
                        floenew = [floenew Floe1(2:end)];
                    else
                        Floe(i) = Floe1;
                    end
                    if length(Floe2) > 1
                        Floe(c(ii)) = Floe2(1);
                        floenew = [floenew Floe2(2:end)];
                    else
                        Floe(c(ii)) = Floe2;
                    end
                    
                end
                if c2 == 1
                    [Floe1, ~] = ridging(Floe(i),FloeB,c2_boundary_poly,false);
                    if length(Floe1) > 1
                        Floe(i) = Floe1(1);
                        floenew = [floenew Floe1(2:end)];
                    else
                        Floe(i) = Floe1;
                    end
                end
            end
        end
    end
end
%toc(timer1)
Floe=Floe(1:N0); % ditch the ghost floes.

for ii = 1:length(Floe)
    if isnan(Floe(ii).ksi_ice)
        xx = 1;
        xx(1) = [1 2];
    end
end

warning('on',id)
warning('on',id3)
end