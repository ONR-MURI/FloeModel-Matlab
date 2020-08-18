function [Floe,dissolvedNEW, SackedOB] = floe_interactions_all(Floe, ocean, winds,heat_flux,c2_boundary_poly, dt,dissolvedNEW,SackedOB,Nx,Ny, RIDGING, PERIODIC,SUBFLOES)
%This function time steps the floe model forward by creating ghost floes if
%periodic. Then calculating the forces between all interacting floes and
%using the forces to time step each individual floes forward. Once these
%floes have been advanced if ridging is allowed this process is run


%Create a polyshape for the boundary 
Lx= max(c2_boundary_poly.Vertices(:,1));
Ly= max(c2_boundary_poly.Vertices(:,2));%c2 must be symmetric around x=0 for channel boundary conditions.
x=[-1 -1 1 1 -1]*Lx*2; 
y=[-1 1 1 -1 -1]*Ly*2;
polybound = polyshape(x,y);
c2_poly = subtract(polybound,c2_boundary_poly);

%% Create ghost floes in the case that periodicity is being used

%Identify the length of original floe variable before creating ghost floes
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
            translation = [translation; -2*Lx*sign(x(i)) 0];
            
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

%Idenfity all interacting floes
for i=1:N 
    
    Floe(i).interactions=[];
    
    Floe(i).OverlapArea=0;
    
    Floe(i).potentialInteractions=[];
    
    Floe(i).collision_force=[0 0];
    
    Floe(i).collision_torque=0;
    
    k=1;
    
    if ( alive(i) && ~isnan(x(i)) )
        for j=1:N
            if j>i && alive(j)
                                
                if sqrt((x(i)-x(j))^2 + (y(i)-y(j))^2)<(rmax(i)+rmax(j)) % if floes are potentially overlapping
                Floe(i).potentialInteractions(k).floeNum=j;
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

%% Calculate forces between interacting floes

parfor i=1:N  %now the interactions could be calculated in a parfor loop!
    
    c1=Floe(i).poly;
            
    if ~isempty(Floe(i).potentialInteractions)
        
        NpotInter=length(Floe(i).potentialInteractions);
        Floe(i).killed = zeros(1,NpotInter);
        for k=1:NpotInter
            
            %Create polyshape for interacting floes
            floeNum=Floe(i).potentialInteractions(k).floeNum;            
            c2=Floe(i).potentialInteractions(k).c;
            
            %Calcualte the forces between the two
            [force_j,P_j,worked,live] = floe_interactions_poly(c1,c2);
            
            %If any floes do not survive interactions put their mass in
            %dissolved variable
            if live(1) == 0
                Floe(i).alive = 0;
            elseif live(2) == 0
                Floe(i).killed(k) = floeNum;
            end
            
            
            if sum(abs(force_j(:)))~=0 && abs(P_j(1))<Lx && abs(P_j(2))<Ly %make sure to count only interactions inside the domain, even for ghost floes
                Floe(i).interactions=[Floe(i).interactions ; floeNum*ones(size(force_j,1),1) force_j P_j zeros(size(force_j,1),1)];
            end
            
        end
        
    end
    
    %If not periodic then calculate interactions with boundary
    if ~PERIODIC
        [force_b, P_j, worked,live] = floe_interactions_poly(c1, c2_poly);
        if live(1) == 0
            Floe(i).alive = 0;
        end
        if sum(abs(force_b(:)))~=0 %; subtracting 0.1m to ensure that interactions with x=+-Lx boundary are excluded as this is a periodic boundary.
            Floe(i).interactions=[Floe(i).interactions ; Inf*ones(size(force_b,1),1) force_b P_j zeros(size(force_b,1),1)];
            Floe(i).OverlapArea=area(c1)-area(intersect(c1,c2_boundary_poly)); % overlap area with the boundary only
        end
    end
end

if numel(fieldnames(Floe))>32
    for i = 1:N    
        if ~isempty(Floe(i).killed(Floe(i).killed>0))
            killed = Floe(i).killed(Floe(i).killed>0);
            killed = killed(killed>0);
            if ~isempty(killed)
                for j = 1:length(killed)
                    Floe(killed(j)).alive = 0;
                end
            end
        end
    end    
    Floe=rmfield(Floe,{'killed'});
end
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
parfor i=1:N
    
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
       
       if isinf(Floe(i).collision_torque) || isinf(max(abs(Floe(i).collision_force)))
           xx =1;
           xx(1) = [1 2];
       end
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

%% 

%Do the timestepping for parent floes now that their forces and torques are known.
parfor i=1:N0
    
    if Floe(i).alive
        
        if PERIODIC
            
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
       
        
        tmp=calc_trajectory(dt,ocean, winds,Floe(i),heat_flux,SUBFLOES); % calculate trajectory
        if (isempty(tmp) || isnan(Floe(i).Xi) )
            Floe(i).alive=0; 
            SackedOB = SackedOB +1;
        else
            Floe(i)=tmp;
        end
        

    end
    

    
end

%% Now performing ridging

floenew = [];
SimpMin = @(A) 3*log10(A);%create function for limit on number of vertices a floe can have

if RIDGING
    %Create a function to control probability that ridging will occur
    overlapArea=cat(1,Floe.OverlapArea)./cat(1,Floe.area);
    keep=rand(length(Floe),1)<5*overlapArea/max(overlapArea);
    for i=1:N0
        
        if Floe(i).alive && keep(i)
            
            if ~isempty(Floe(i).potentialInteractions)
                
                
                for k = 1:length(Floe(i).potentialInteractions)
               
                    
                    if Floe(i).potentialInteractions(k).floeNum<= N0
                        floeNum = Floe(i).potentialInteractions(k).floeNum;
                        [Floe1,Floe2] = ridging(Floe(i),Floe(Floe(i).potentialInteractions(k).floeNum),c2_boundary_poly,PERIODIC);
                    else
                        %Identify parent floes that can ridge and move them
                        %to ghost cell location
                        floeNum = parent(Floe(i).potentialInteractions(k).floeNum - N0);
                        trans = translation(Floe(i).potentialInteractions(k).floeNum-N0,:);
                        if floeNum > N0
                            trans = trans+translation(floeNum-N0,:);
                            floeNum = parent(floeNum-N0);
                        end
                        Floe2 = Floe(floeNum);
                        
                        Floe2.poly=translate(Floe2.poly,trans);
                        for ii = 1:length(Floe2.SubFloes)
                            Floe2.SubFloes(ii).poly=translate(Floe2.SubFloes(ii).poly,trans);
                        end
                        Floe2.vorX=Floe2.vorX+trans(1);
                        Floe2.vorbox(:,1)=Floe2.vorbox(:,1)+trans(1);
                        Floe2.Xm=Floe2.Xm+trans(1);
                        Floe2.Xi=Floe2.Xi+trans(1);
                        
                        Floe2.vorY=Floe2.vorY+trans(2);
                        Floe2.vorbox(:,2)=Floe2.vorbox(:,2)+trans(2);
                        Floe2.Ym=Floe2.Ym+trans(2);
                        Floe2.Yi=Floe2.Yi+trans(2);
                        
                        overlap = intersect(Floe2.poly,Floe(Floe(i).potentialInteractions(k).floeNum).poly);
                        if area(overlap)/Floe2.area>0.9
                            
                            [Floe1,Floe2] = ridging(Floe(i),Floe2,c2_boundary_poly,PERIODIC);
                        else
                            Floe1 = Floe(i);
                        end
                        
                        %Translate back to original floe location
                        if trans(1)
                            for jj = 1:length(Floe2)
                                Floe2(jj).poly=translate(Floe2(jj).poly,[-trans(1) 0]);
                                for ii = 1:length(Floe2(jj).SubFloes)
                                    Floe2(jj).SubFloes(ii).poly=translate(Floe2(jj).SubFloes(ii).poly,[-trans(1) 0]);
                                end
                                Floe2(jj).vorX=Floe2(jj).vorX-trans(1);
                                Floe2(jj).vorbox(:,1)=Floe2(jj).vorbox(:,1)-trans(1);
                                Floe2(jj).Xm=Floe2(jj).Xm-trans(1);
                                Floe2(jj).Xi=Floe2(jj).Xi-trans(1);
                            end
                        end
                        
                        if trans(2)
                            for jj = 1:length(Floe2)
                                Floe2(jj).poly=translate(Floe2(jj).poly,[0 -trans(2)]);
                                for ii = 1:length(Floe2(jj).SubFloes)
                                    Floe2(jj).SubFloes(ii).poly=translate(Floe2(jj).SubFloes(ii).poly,[0 -trans(2)]);
                                end
                                Floe2(jj).vorY=Floe2(jj).vorY-trans(2);
                                Floe2(jj).vorbox(:,2)=Floe2(jj).vorbox(:,2)-trans(2);
                                Floe2(jj).Ym=Floe2(jj).Ym-trans(2);
                                Floe2(jj).Yi=Floe2(jj).Yi-trans(2);
                            end
                        end
                        
                    end
                    
                    %Check if any floe simplification is necessary
                    for kk = 1:length(Floe1)
                        floe = Floe1(kk);
                        ddx = 100;
                        if length(floe.poly.Vertices) > SimpMin(floe.area)
                            floe = FloeSimplify(floe, ddx,SUBFLOES);
                        end
                        for jj = 1:length(floe)
                            
                            if kk == 1 && jj == 1
                                Floe(i) = floe(jj);
                            else
                                floenew = [floenew floe(jj)];
                            end
                        end
                    end
                    
                    for kk = 1:length(Floe2)
                        floe = Floe2(kk);
                        if length(floe.poly.Vertices) > SimpMin(floe.area)
                            floe = FloeSimplify(floe, ddx,SUBFLOES);
                        end
                        for jj = 1:length(floe)
                            
                            if kk == 1 && jj == 1
                                Floe(floeNum) = floe(jj);
                            else
                                floenew = [floenew floe(jj)];
                            end
                        end
                    end
                    
                    
                end

                
            end
            
        end
    end
    
    
end

Floe=Floe(1:N0); % ditch the ghost floes.
Floe = [Floe floenew]; %add on any new floes that were if ridging caused one floe to become multiple

Floe=rmfield(Floe,{'potentialInteractions'}); %remove potential interactions field
live = cat(1,Floe.alive);
for ii = N0:-1:1
    if live(ii)==0
        dissolvedNEW = dissolvedNEW+calc_dissolved_mass(Floe(ii),Nx,Ny,c2_boundary_poly);
    end
end
Floe(live==0)=[]; %remove any floes that got dissolved so they do not take up space


end