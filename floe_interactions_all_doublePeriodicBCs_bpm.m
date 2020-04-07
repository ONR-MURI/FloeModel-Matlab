function [Floe,dissolvedNEW, SackedOB] = floe_interactions_all_doublePeriodicBCs_bpm(Floe, ocean, winds,c2_boundary_poly, dt,dissolvedNEW,SackedOB,Nx,Ny, RIDGING, PERIODIC,SUBFLOES)

Lx= max(c2_boundary_poly.Vertices(:,1));
Ly= max(c2_boundary_poly.Vertices(:,2));%c2 must be symmetric around x=0 for channel boundary conditions.
x=[-1 -1 1 1 -1]*Lx*2; 
y=[-1 1 1 -1 -1]*Ly*2;
polybound = polyshape(x,y);
c2_poly = subtract(polybound,c2_boundary_poly);

%%
N0=length(Floe);

Lx= max(c2_boundary_poly.Vertices(:,1));
Ly= max(c2_boundary_poly.Vertices(:,2));%c2 must be symmetric around x=0 for channel boundary conditions.

if PERIODIC
    
    ghostFloeX=[];
    ghostFloeY=[];
    parent=[];
    
    x=cat(1,Floe.Xi);
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        
        %   if alive(i) && (x(i)>Lx-rmax(i)) || (x(i)<-Lx+rmax(i))
        if alive(i) && (max(abs(Floe(i).poly.Vertices(:,1)))>Lx)
            
            ghostFloeX=[ghostFloeX  Floe(i)];
            ghostFloeX(end).poly=translate(Floe(i).poly,[-2*Lx*sign(x(i)) 0]);
            for ii = 1:length(Floe(i).SubFloes)
                ghostFloeX(end).SubFloes(ii).poly=translate(Floe(i).SubFloes(ii).poly,[-2*Lx*sign(x(i)) 0]);
            end
            ghostFloeX(end).Xi=Floe(i).Xi-2*Lx*sign(x(i));
            ghostFloeX(end).vorX=Floe(i).vorX-2*Lx*sign(x(i));
            ghostFloeX(end).vorbox(:,1)=Floe(i).vorbox(:,1)-2*Lx*sign(x(i));
            parent=[parent  i];
            
        end
        
        
    end
    
    Floe=[Floe ghostFloeX];
    
    x=cat(1,Floe.Xi);
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        
        %   if alive(i) && (x(i)>Lx-rmax(i)) || (x(i)<-Lx+rmax(i))
        if alive(i) && (max(abs(Floe(i).poly.Vertices(:,2)))>Ly)
            
            ghostFloeY=[ghostFloeY  Floe(i)];
            ghostFloeY(end).poly=translate(Floe(i).poly,[0 -2*Ly*sign(y(i))]);
            for ii = 1:length(Floe(i).SubFloes)
                ghostFloeY(end).SubFloes(ii).poly=translate(Floe(i).SubFloes(ii).poly,[0 -2*Ly*sign(y(i))]);
            end
            ghostFloeY(end).Yi=Floe(i).Yi-2*Ly*sign(y(i));
            ghostFloeY(end).vorY=Floe(i).vorY-2*Ly*sign(y(i));
            ghostFloeY(end).vorbox(:,2)=Floe(i).vorbox(:,2)-2*Ly*sign(y(i));
            parent=[parent  i];
            
        end
        
    end
    
    Floe=[Floe ghostFloeY];
    
end

N=length(Floe);
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
            if j>i && alive(j)
                                
                if sqrt((x(i)-x(j))^2 + (y(i)-y(j))^2)<(rmax(i)+rmax(j)) % if floes are potentially overlapping
                Floe(i).potentialInteractions(k).floeNum=j;
%                Floe(i).potentialInteractions(k).c=[Floe(j).c_alpha(1,:)+Floe(j).Xi; Floe(j).c_alpha(2,:)+Floe(j).Yi];
           %     Floe(i).potentialInteractions(k).c=Floe(j).poly;
                k=k+1;
            
                end
            end
            
        end
        
    end
    
    if ~isempty(Floe(i).potentialInteractions)
        potFloes=[Floe([Floe(i).potentialInteractions.floeNum]).poly];
        a=area(intersect(Floe(i).poly,potFloes));
        Floe(i).potentialInteractions=Floe(i).potentialInteractions(a>0); %ditch floes that have zero overlap
        b=num2cell(a(a>0)); [Floe(i).potentialInteractions.overlapArea]=b{:}; %save the overlap area
        [Floe(i).potentialInteractions.c]= Floe([Floe(i).potentialInteractions.floeNum]).poly; %copy polyshapes for parfor calculation of forces.
    end
end

%%

for i=1:N  %now the interactions could be calculated in a parfor loop!
    
    c1=Floe(i).poly;
            
    if ~isempty(Floe(i).potentialInteractions)
        
        NpotInter=length(Floe(i).potentialInteractions);
        for k=1:NpotInter
            
            floeNum=Floe(i).potentialInteractions(k).floeNum;
            
            %Floe(i).potentialInteractions(k).overlapArea=0;
            
            c2=Floe(i).potentialInteractions(k).c;
            
            [force_j,P_j,worked] = floe_interactions_bpm2(c1,c2);
            
            %if ~worked, disp(['potential contact issue for floes (' num2str(i) ',' num2str(floeNum) ')' ]); end
            
            if sum(abs(force_j(:)))~=0 && abs(P_j(1))<Lx && abs(P_j(2))<Ly %make sure to count only interactions inside the domain, even for ghost floes
                Floe(i).interactions=[Floe(i).interactions ; floeNum*ones(size(force_j,1),1) force_j P_j zeros(size(force_j,1),1)];
              %  Floe(i).potentialInteractions(k).overlapArea=area(intersect(c1,c2));
            end
            
        end
%         if NpotInter > 1
%             c2 = union([Floe(i).potentialInteractions.c]);
%         else
%             c2 = Floe(i).potentialInteractions.c;
%         end
%         
%         [force_j,P_j,worked] = floe_interactions_bpm2(c1,c2);
%                 
%         if sum(abs(force_j))~=0
%             Floe(i).interactions=[Floe(i).interactions ; ones(size(force_j,1),1) force_j P_j zeros(size(force_j,1),1)];
%         end
        
    end
    
    if ~PERIODIC
        [force_b, P_j, worked] = floe_interactions_bpm2(c1, c2_poly);
        if sum(abs(force_b(:)))~=0 && abs(P_j(1))<(Lx-0.1) %; subtracting 0.1m to ensure that interactions with x=+-Lx boundary are excluded as this is a periodic boundary.
            % boundary will be recorded as floe number Inf;
            Floe(i).interactions=[Floe(i).interactions ; Inf*ones(size(force_b,1),1) force_b P_j zeros(size(force_b,1),1)];
            Floe(i).OverlapArea=area(c1)-area(intersect(c1,c2_boundary_poly)); % overlap area with the boundary only
        end
    end
end

%Floe=rmfield(Floe,{'potentialInteractions'});

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
% calculate all torques from forces including for ghost floes
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
        
%        if abs(Floe(i).collision_force(1)/Floe(i).mass) > 0.35
%            if Floe(i).area > 3500
%                xx = 1;
%                xx(1) = [1 2];
%            end
%        elseif abs(Floe(i).collision_force(2)/Floe(i).mass) > 0.35
%            if Floe(i).area > 3500
%                xx = 1;
%                xx(1) = [1 2];
%            end
%        end
    end
    
end

if PERIODIC
    %add forces and torques from ghost floes to their parents; ghost floes
    %begin with the index N0+1
    for i=1:length(parent)
        Floe(parent(i)).collision_force =Floe(parent(i)).collision_force +Floe(N0+i).collision_force;
        Floe(parent(i)).collision_torque=Floe(parent(i)).collision_torque+Floe(N0+i).collision_torque;
    end
end


%Do the timestepping for parent floes now that their forces and torques are known.
for i=1:N0
    
    if Floe(i).alive
        
        if PERIODIC
            
            floe = Floe(i);
            if abs(Floe(i).Xi)>Lx %if floe got out of periodic bounds, put it on the other end
                %
                Floe(i).poly=translate(Floe(i).poly,[-2*Lx*sign(Floe(i).Xi) 0]);
                Floe(i).Xm=Floe(i).Xm-2*Lx*sign(Floe(i).Xi);
                for ii = 1:length(Floe(i).SubFloes)
                    Floe(i).SubFloes(ii).poly=translate(Floe(i).SubFloes(ii).poly,[-2*Lx*sign(Floe(i).Xi) 0]);
                end
                Floe(i).vorX=Floe(i).vorX-2*Lx*sign(Floe(i).Xi);
                Floe(i).vorbox(:,1)=Floe(i).vorbox(:,1)-2*Lx*sign(Floe(i).Xi);
                Floe(i).Xi=Floe(i).Xi-2*Lx*sign(Floe(i).Xi);
            end
            
            if abs(Floe(i).Yi)>Ly %if floe got out of periodic bounds, put it on the other end
                %
                Floe(i).poly=translate(Floe(i).poly,[0 -2*Ly*sign(Floe(i).Yi)]);
                Floe(i).Ym=Floe(i).Ym-2*Ly*sign(Floe(i).Yi);
                for ii = 1:length(Floe(i).SubFloes)
                    Floe(i).SubFloes(ii).poly=translate(Floe(i).SubFloes(ii).poly,[0 -2*Ly*sign(Floe(i).Yi)]);
                end
                Floe(i).vorY=Floe(i).vorY-2*Ly*sign(Floe(i).Yi);
                Floe(i).vorbox(:,2)=Floe(i).vorbox(:,2)-2*Ly*sign(Floe(i).Yi);
                Floe(i).Yi=Floe(i).Yi-2*Ly*sign(Floe(i).Yi);
            end
            
        end
        

        
        tmp=calc_trajectory(dt,ocean, winds,Floe(i)); % calculate trajectory
        if (isempty(tmp) || isnan(Floe(i).Xi) )
            Floe(i).alive=0; 
            SackedOB = SackedOB +1;
        elseif Floe(i).alive == 0
            dissolvedNEW = calc_vol_dissolved(Floe(i),Nx,Ny,c2_boundary_poly);
            Floe(i) = [];
        else
            Floe(i)=tmp;
        end

    end
    
end

if RIDGING
    
    for i=1:N0
        floe = Floe(i);
        if ~isempty(Floe(i).interactions)
            
            if ~isempty(Floe(i).potentialInteractions)
                
                for k = 1:length(Floe(i).potentialInteractions)

                    [Floe(i),Floe(Floe(i).potentialInteractions(k).floeNum),dissolvedNEW] = ridging(dissolvedNEW,Floe(i),Floe(Floe(i).potentialInteractions(k).floeNum),Nx,Ny,c2_boundary_poly,PERIODIC,SUBFLOES);                  
                    
                    if Floe(i).poly.NumRegions > 1
                        polyout = sortregions(Floe(i).poly,'area','descend');
                        R = regions(polyout);
                        polynew = R(1);
                        Floe(i).poly = polynew;
                    elseif Floe(Floe(i).potentialInteractions(k).floeNum).poly.NumRegions > 1
                        polyout = sortregions(Floe(Floe(i).potentialInteractions(k).floeNum).poly,'area','descend');
                        R = regions(polyout);
                        polynew = R(1);
                        Floe(Floe(i).potentialInteractions(k).floeNum).poly=polynew;
                    end
                end
                
            end
        end
    end
    
    
end

Floe=Floe(1:N0); % ditch the ghost floes.

Floe=rmfield(Floe,{'potentialInteractions'});
end
