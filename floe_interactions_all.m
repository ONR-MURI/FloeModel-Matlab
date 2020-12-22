function [Floe,dissolvedNEW] = floe_interactions_all(Floe, ocean, winds,c2_boundary, dt, HFo, min_floe_size, Nx,Ny,Nb, dissolvedNEW,COLLISION, PERIODIC, RIDGING)

id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)

Lx= max(c2_boundary(1,:));
Ly= max(c2_boundary(2,:));
height.mean = 1; height.delta = 0;
c2_boundary_poly = polyshape(c2_boundary');
FloeB = initialize_floe_values(c2_boundary_poly, height);
live = cat(1,Floe.alive);
Floe(live==0)=[];
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
    
    Floe(i).Stress=zeros(2);
    
    Floe(i).collision_torque=0;
    
    k=1;
    
%     if abs(Floe(i).area/area(polyshape(Floe(i).c_alpha'))-1)>1e-3
%         xx = 1;
%         xx(1) =[1 2];
%     end
    
    if ( alive(i) && ~isnan(x(i)) ) && COLLISION
        for j=1:N
            %display(j);
            if j>i && alive(j) && sqrt((x(i)-x(j))^2 + (y(i)-y(j))^2)<(rmax(i)+rmax(j)) % if floes are potentially overlapping
                Floe(i).potentialInteractions(k).floeNum=j;
                Floe(i).potentialInteractions(k).c=[Floe(j).c_alpha(1,:)+x(j); Floe(j).c_alpha(2,:)+y(j)];
                Floe(i).potentialInteractions(k).Ui=Floe(j).Ui;
                Floe(i).potentialInteractions(k).Vi=Floe(j).Vi;
                Floe(i).potentialInteractions(k).Xi=x(j);
                Floe(i).potentialInteractions(k).Yi=y(j);
                Floe(i).potentialInteractions(k).ksi_ice = Floe(j).ksi_ice;
                k=k+1;
            end
            
        end
        
    end
    
end

weld = zeros(length(Floe),1);
kill = zeros(1,N0);

parfor i=1:N  %now the interactions could be calculated in a parfor loop!
        
    %c1=[Floe(i).c_alpha(1,:)+x(i); Floe(i).c_alpha(2,:)+y(i)];
    
    if ~isempty(Floe(i).potentialInteractions)
        
        for k=1:length(Floe(i).potentialInteractions)
            
            floeNum=Floe(i).potentialInteractions(k).floeNum;
            
            %c2=Floe(i).potentialInteractions(k).c;
            
            %[force_j,P_j, overlap] = floe_interactions(c1,c2,c2_boundary,PERIODIC);
            [force_j,P_j, overlap] = floe_interactions(Floe(i),Floe(i).potentialInteractions(k),c2_boundary,PERIODIC);
            
            %if ~worked, disp(['contact points issue for (' num2str(i) ',' num2str(floeNum) ')' ]); end
            
            if sum(abs(force_j))~=0
                Floe(i).interactions=[Floe(i).interactions ; floeNum*ones(size(force_j,1),1) force_j P_j zeros(size(force_j,1),1) overlap'];
                Floe(i).OverlapArea = sum(overlap)+Floe(i).OverlapArea;
            elseif isinf(overlap)
                if i <= N0 && sign(overlap)>0
                    kill(i) = i;%floeNum;
                elseif floeNum <= N0
                    kill(i) = floeNum;
                end
            end
            
        end
        
    end
    
    if ~PERIODIC
        c1=[Floe(i).c_alpha(1,:)+x(i); Floe(i).c_alpha(2,:)+y(i)];
        [force_b, P_j, overlap] = floe_interactions_bound(c1, c2_boundary,c2_boundary,PERIODIC);
        in = inpolygon(x(i),y(i),c2_boundary(1,:)',c2_boundary(2,:)');
        if ~in
            Floe(i).alive = 0;
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
% for i = 1+Nb:length(Floe)
%     if weld(i)>0 && Floe(i).alive > 0 && Floe(weld(i)).alive > 0
%         Floe(i).poly = polyshape(Floe(i).c_alpha'+[Floe(i).Xi Floe(i).Yi]);
%         Floe(weld(i)).poly = polyshape(Floe(weld(i)).c_alpha'+[Floe(weld(i)).Xi Floe(weld(i)).Yi]);
%         [floenew] = FuseFloes(Floe(i),Floe(weld(i)));
%         if floenew.h > 3
%             xx = 1;
%             xx(1) = [1 2];
%         elseif weld(i) > N0
%             xx = 1;
%             xx(1) = [1 2];
%         end
% %         if isfield(floenew,'poly')
% %             floenew=rmfield(floenew,{'poly'});
% %         end
%         Mass = [Floe(i).mass Floe(weld(i)).mass];
%         if Mass(1) > Mass(2)
%             floenew.interactions = Floe(i).interactions;
%             floenew.OverlapArea = Floe(i).OverlapArea;
%             floenew.potentialInteractions = Floe(i).potentialInteractions;
%             Floe(i) = floenew;
%             Floe(weld(i)).alive = 0; kill(weld(i)) = 0;
%             Floe(weld(i)).interactions = [];
%         else
%             floenew.interactions = Floe(weld(i)).interactions;
%             floenew.OverlapArea = Floe(weld(i)).OverlapArea;
%             floenew.potentialInteractions = Floe(weld(i)).potentialInteractions;
%             Floe(weld(i)) = floenew;
%             Floe(i).alive = 0; kill(i) = 0;
%             Floe(i).interactions = [];
%         end
%     end
% end
alive=cat(1,Floe.alive);
if isfield(Floe,'poly')
    Floe=rmfield(Floe,{'poly'});
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
    
    parfor i=N0+1:N %do this in parfor
        
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

keep = ones(1,N0);
parfor i=1+Nb:N0
    
    if ~isempty(Floe(i).interactions)
        
       a=Floe(i).interactions;
       r=[x(i) y(i)];
        for k=1:size(a,1)
            floe_Rforce=a(k,4:5);
            floe_force=a(k,2:3);
            floe_torque=cross([floe_Rforce-r 0], [floe_force 0]);
            Floe(i).interactions(k,6)=floe_torque(3);
        end
        
%        Floe(i).Stress =.5*([sum((a(:,4)-r(1)).*a(:,2)) sum((a(:,5)-r(2)).*a(:,2)); sum((a(:,4)-r(1)).*a(:,3)) sum((a(:,5)-r(2)).*a(:,3))]+[sum(a(:,2).*(a(:,4)-r(2))) sum(a(:,3).*(a(:,4)-r(2))); sum(a(:,2).*(a(:,5)-r(2))) sum(a(:,3).*(a(:,5)-r(2)))]);
       Floe(i).collision_force=sum(Floe(i).interactions(:,2:3),1)+Floe(i).collision_force;
       Floe(i).collision_torque=sum(Floe(i).interactions(:,6),1)+Floe(i).collision_torque;
        
    end
    
    if PERIODIC
            
        if abs(Floe(i).Xi)>Lx %if floe got out of periodic bounds, put it on the other end
            Floe(i).Xi=Floe(i).Xi-2*Lx*sign(Floe(i).Xi);
        end
        
        if abs(Floe(i).Yi)>Ly %if floe got out of periodic bounds, put it on the other end
            Floe(i).Yi=Floe(i).Yi-2*Ly*sign(Floe(i).Yi);
        end
        
    end
    
   %Do the timestepping now that forces and torques are known.
    if Floe(i).alive,
%         [tmp,frac,Fx,Fy]=calc_trajectory_Nares(dt,ocean, winds,Floe(i),HFo,c2_boundary); % calculate trajectory
        [tmp,frac,Fx,Fy]=calc_trajectory(dt,ocean, winds,Floe(i),HFo); % calculate trajectory
        if (isempty(tmp) || isnan(x(i)) ), kill(i)=i; elseif frac == 1, keep(i) = 0; else Floe(i)=tmp; Floe(i).Fx = Fx; Floe(i).Fy = Fy; end
    end
        
end

%timer1 = tic;
floenew = [];
if RIDGING
    %Create a function to control probability that ridging will occur
    overlapArea=cat(1,Floe.OverlapArea)./cat(1,Floe.area);
    keepR=rand(length(Floe),1)>10*overlapArea;
%     keep=rand(length(Floe),1)<0.5;
    for i=1+Nb:N0

        if Floe(i).alive && ~keepR(i)
            c2 = 0;
            if ~isempty(Floe(i).potentialInteractions)
                c = cat(1,Floe(i).potentialInteractions.floeNum);
                if isinf(max(c))
                    c2 = 1;
                end
                c = c(c<Inf);
                c = c(c>Nb);
                for ii = 1:length(c)
                    [Floe1, Floe2] = ridging(Floe(i),Floe(c(ii)),c2_boundary_poly,PERIODIC,min_floe_size);

                    if length(Floe1) > 1
                        Floe(i) = Floe1(1);
                        floenew = [floenew Floe1(2:end)];
                    else
                        Floe(i) = Floe1;
                        if Floe1.alive == 0
                            kill(i) = i;
                        end
                    end
                    if length(Floe2) > 1
                        Floe(c(ii)) = Floe2(1);
                        floenew = [floenew Floe2(2:end)];
                    else
                        Floe(c(ii)) = Floe2;
                        if Floe2.alive == 0 && c(ii) <= N0
                            kill(c(ii)) = c(ii);
                        end
                    end
                    
                end
                if c2 == 1
                    [Floe1, ~] = ridging(Floe(i),FloeB,c2_boundary_poly,PERIODIC,min_floe_size);
                    if length(Floe1) > 1
                        Floe(i) = Floe1(1);
                        floenew = [floenew Floe1(2:end)];
                    else
                        Floe(i) = Floe1;
                        if Floe1.alive == 0
                            kill(i) = i;
                        end
                    end
                end
            end
        end
    end
end
%toc(timer1)
Floe=Floe(1:N0); % ditch the ghost floes.

% for ii = 1+Nb:length(Floe)
%     if isnan(Floe(ii).ksi_ice)
%         xx = 1;
%         xx(1) = [1 2];
%     end
% end
if ~isempty(kill(kill>0)) 
    kill = unique(kill(kill>0));
    kill = sort(kill,'descend');
    for ii = length(kill)
        dissolvedNEW = dissolvedNEW+calc_dissolved_mass(Floe(kill(ii)),Nx,Ny,c2_boundary_poly);
        Floe(kill(ii)).alive = 0;
        %     end
    end
end
Floe(live==0)=[]; %remove any floes that got dissolved so they do not take up space
live = cat(1,Floe.alive);
notalive = logical(abs(live-1));
keep = logical(keep+notalive);
% fracturedFloes=fracture_floe(Floe(~keep),3);
% if ~isempty(fracturedFloes)
%     Floe=[Floe(keep) fracturedFloes];
% end
% for ii = 1:length(Floe)
%     if abs(Floe(ii).area/area(polyshape(Floe(ii).c_alpha'))-1)>1e-3
%         xx = 1;
%         xx(1) =[1 2];
%     end
% end


warning('on',id)
warning('on',id3)
end