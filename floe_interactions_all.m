function Floe = floe_interactions_all(Floe, ocean, winds,c2_boundary, dt)


N=length(Floe);
%Floe(i).interactions=[floeNumber Fx Fy px py torque];
%Floe(i).potentialInteractions(j).floeNum
%Floe(i).potentialInteractions(j).c_alpha=Floe(floeNum).c_alpha.
Vd = 0;
x = c2_boundary_poly.Vertices(:,1);
y = c2_boundary_poly.Vertices(:,2);
dx = max(x)-min(x);
dy = max(y)-min(y);
xb = [min(x)-dx min(x)-dx max(x)+dx max(x)+dx min(x)-dx];
yb = [min(y)-dy max(y)+dy max(y)+dy min(y)-dy min(y)-dy];
poly3 = polyshape(xb,yb);
c2_poly = subtract(poly3,c2_boundary_poly);
%%

x=cat(1,Floe.Xi);
y=cat(1,Floe.Yi);
rmax=cat(1,Floe.rmax);
alive=cat(1,Floe.alive);

for i=1:N  %do interactions with boundary in a separate parfor loop
    
    Floe(i).interactions=[];
    
    Floe(i).OverlapArea=0;
    
    Floe(i).potentialInteractions=[];
    
    Floe(i).collision_force=[0 0];
    
    Floe(i).collision_torque=0;
    
    k=1;
    
    if ( alive(i) && ~isnan(x(i)) )
        for j=1:N
            %display(j);
            if j>i && alive(j) && sqrt((x(i)-x(j))^2 + (y(i)-y(j))^2)<(rmax(i)+rmax(j)) % if floes are potentially overlapping
                Floe(i).potentialInteractions(k).floeNum=j;
%                Floe(i).potentialInteractions(k).c=[Floe(j).c_alpha(1,:)+Floe(j).Xi; Floe(j).c_alpha(2,:)+Floe(j).Yi];
                Floe(i).potentialInteractions(k).c=Floe(j).poly;
                k=k+1;
            end
            
        end
        
    end
    
end

%%

c2_boundary_poly=polyshape(c2_boundary(1,:),c2_boundary(2,:));

for i=1:N  %now the interactions could be calculated in a parfor loop!
        
%   c1=[Floe(i).c_alpha(1,:)+Floe(i).Xi; Floe(i).c_alpha(2,:)+Floe(i).Yi];
    c1=Floe(i).poly;
    
            
    if ~isempty(Floe(i).potentialInteractions)
        
        for k=1:length(Floe(i).potentialInteractions)
            
            floeNum=Floe(i).potentialInteractions(k).floeNum;
            
            Floe(i).potentialInteractions(k).overlapArea=0;
            
            
            %add ridging
            [Floe(i),Floe(Floe(i).potentialInteractions(k).floeNum),Vd] = ridging(Vd,Floe(i),Floe(Floe(i).potentialInteractions(k).floeNum));
            
            c1=Floe(i).poly;
            c2=Floe(i).potentialInteractions(k).c;
            
<<<<<<< Updated upstream
%             [force_j,P_j,worked] = floe_interactions_poly(c1,c2);
            [force_j,P_j,worked] = floe_interactions_bpm(c1,c2);
            [force_j,P_j,worked] = floe_interactions_bpm(c1,c2);
%             [force_j,P_j,worked] = floe_interactions_poly(c1,c2);
=======
%             [force_j,P_j,worked] = floe_interactions_bpm(c1,c2);
            [force_j,P_j,worked] = floe_interactions_poly(c1,c2);
>>>>>>> Stashed changes
            
            %if ~worked, disp(['potential contact issue for floes (' num2str(i) ',' num2str(floeNum) ')' ]); end
            
            if sum(abs(force_j))~=0
                Floe(i).interactions=[Floe(i).interactions ; floeNum*ones(size(force_j,1),1) force_j P_j zeros(size(force_j,1),1)];

                Floe(i).OverlapArea=Floe(i).OverlapArea+area(intersect(c1,c2));
           
                Floe(i).potentialInteractions(k).overlapArea=area(intersect(c1,c2));
            end
            
        end
        
    end
    
<<<<<<< Updated upstream
    [force_b, P_j, worked] = floe_interactions_bpm(c1, c2_poly);
%     [force_b, P_j, worked] = floe_interactions_poly(c1, c2_boundary_poly);
=======
    [force_b, P_j, worked] = floe_interactions_poly(c1, c2_boundary_poly);
%     [force_b, P_j, worked] = floe_interactions_bpm(c1, c2_poly);
>>>>>>> Stashed changes
    %if ~worked, disp(['potential contact issue for floes (' num2str(i) ', boundary)' ]); end
    if sum(abs(force_b))~=0
        % boundary will be recorded as floe number Inf;
        Floe(i).interactions=[Floe(i).interactions ; Inf*ones(size(force_b,1),1) force_b P_j zeros(size(force_b,1),1)];
        Floe(i).OverlapArea=Floe(i).OverlapArea+(area(c1)-area(intersect(c1,c2_boundary_poly)));
    end
    
end


%Floe=rmfield(Floe,{'potentialInteractions'});
% for ii = 1:N
%     if Floe(ii).alive == 1
%         if isempty(Floe(ii).potentialInteractions) ==0
%             for jj = 1:length(Floe(ii).potentialInteractions)
%                 if Floe(Floe(ii).potentialInteractions(jj).floeNum).alive == 1 && Floe(ii).alive == 1
%                     [Floe(ii),Floe(Floe(ii).potentialInteractions(jj).floeNum),Vd] = ridging(Vd,Floe(ii),Floe(Floe(ii).potentialInteractions(jj).floeNum));
%                 end
%             end
%         end
%     end
% end

%%
%Fill the lower part of the interacton matrix (floe_i,floe_j) for floes with j<i
for i=1:N %this has to be done sequentially
      
    if ~isempty(Floe(i).interactions)
        
        a=Floe(i).interactions;
        
        indx=a(:,1);
        
        for j=1:length(indx)
            
            if indx(j)<=N && indx(j)>i
                Floe(indx(j)).interactions=[Floe(indx(j)).interactions; i -a(j,2:3) a(j,4:5) 0];   % 0 is torque here that is to be calculated below
            end
            
        end
        
        
        for j=1:length(Floe(i).potentialInteractions)
            
            overlapArea=Floe(i).potentialInteractions(j).overlapArea;
            floeNum=Floe(i).potentialInteractions(j).floeNum;
            
            Floe(i).OverlapArea=Floe(i).OverlapArea+overlapArea;
            Floe(floeNum).OverlapArea=Floe(floeNum).OverlapArea+overlapArea;
        end
        
    end
    
end

%%

% calculate all torques from forces

for i=1:N
    
    if ~isempty(Floe(i).interactions)
        
       a=Floe(i).interactions;
       r=[Floe(i).Xi Floe(i).Yi];
        for k=1:size(a,1)
            floe_Rforce=a(k,4:5);
            floe_force=a(k,2:3);
            floe_torque=cross([floe_Rforce-r 0], [floe_force 0]);
            Floe(i).interactions(k,6)=floe_torque(3);
        end
        
       Floe(i).collision_force=sum(Floe(i).interactions(:,2:3),1);
       Floe(i).collision_torque=sum(Floe(i).interactions(:,6),1);
        
    end
    
   %Do the timestepping now that forces and torques are known.
    if Floe(i).alive
        tmp=calc_trajectory(dt,ocean, winds,Floe(i)); % calculate trajectory
        if (isempty(tmp) || isnan(Floe(i).Xi) ), Floe(i).alive=0; else Floe(i)=tmp; end
    end
    
end

for i=1:N
    
    if ~isempty(Floe(i).interactions)
        
        for k = 1:length(Floe(i).potentialInteractions)
        
            [Floe(i),Floe(Floe(i).potentialInteractions(k).floeNum),Vd] = ridging(Vd,Floe(i),Floe(Floe(i).potentialInteractions(k).floeNum));
    
        end
    end
end

Floe=rmfield(Floe,{'potentialInteractions'});

for i=1:N             
    if ~isempty(Floe(i).potentialInteractions)        
        NpotInter=length(Floe(i).potentialInteractions);                
        for k=1:NpotInter
            %[Floe(i),Floe(Floe(i).potentialInteractions(k).floeNum),Vd] = ridging(Vd,Floe(i),Floe(Floe(i).potentialInteractions(k).floeNum));
        end
    end
end
Floe=rmfield(Floe,{'potentialInteractions'});
end
