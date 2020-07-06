function [fig, fig2] =plot_Floes(fig, fig2, Time,Floe,ocean,c2_boundary_poly, PERIODIC)
%This function creates plots of the floe state showing the stress and and
%thickness of the floes

showContactPoints=0;
showCenterOfMass=0;
stressColor=1;


%Populate the ghost floes
N0=length(Floe);

Lx= max(c2_boundary_poly.Vertices(:,1)); %c2 must be symmetric around x=0 for channel boundary conditions.
Ly= max(c2_boundary_poly.Vertices(:,2)); 

ghostFloeX=[];
ghostFloeY=[];
parent=[];

x=cat(1,Floe.Xi);
alive=cat(1,Floe.alive);

if PERIODIC
    
    for i=1:length(Floe)
        
        if alive(i) && (max(abs(Floe(i).poly.Vertices(:,1)))>Lx/2)
            
            ghostFloeX=[ghostFloeX  Floe(i)];
            ghostFloeX(end).poly=translate(Floe(i).poly,[-2*Lx*sign(x(i)) 0]);
            for ii = 1:length(Floe(i).SubFloes)
                Floe(i).SubFloes(ii).poly=translate(Floe(i).SubFloes(ii).poly,[-2*Lx*sign(Floe(i).Xi) 0]);
            end
            ghostFloeX(end).Xi=Floe(i).Xi-2*Lx*sign(x(i));
            parent=[parent  i];
            
        end
        
        
    end
    
    Floe=[Floe ghostFloeX];
    
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        
        if alive(i) && (max(abs(Floe(i).poly.Vertices(:,2)))>Ly/2)
            
            ghostFloeY=[ghostFloeY  Floe(i)];
            ghostFloeY(end).poly=translate(Floe(i).poly,[0 -2*Ly*sign(y(i))]);
            for ii = 1:length(Floe(i).SubFloes)
                Floe(i).SubFloes(ii).poly=translate(Floe(i).SubFloes(ii).poly,[0 -2*Ly*sign(Floe(i).Yi)]);
            end
            ghostFloeY(end).Yi=Floe(i).Yi-2*Ly*sign(y(i));
            parent=[parent  i];
            
        end
        
    end
    
    Floe=[Floe ghostFloeY];
    
end

%% Set Up The Plots

ratio=max(ocean.Yo)/max(ocean.Xo);
if (fig==0 || ~isvalid(fig))
    fig=figure('Position',[100 100 1000 1000*ratio],'visible','on');  
    set(fig,'PaperSize',12*[1 ratio],'PaperPosition',12*[0 0 1 ratio]);
end
figure(fig);
clf(fig);

dn=1; % plot every dn'th velocity vector
quiver(ocean.Xo(1:dn:end),ocean.Yo(1:dn:end),ocean.Uocn(1:dn:end,1:dn:end),ocean.Vocn(1:dn:end,1:dn:end));
hold on;

axis([ocean.Xo(1) ocean.Xo(end) ocean.Yo(1) ocean.Yo(end)]);

colormap('gray'); caxis([0 1]);

title(['Time = ' num2str(Time) ' s']);

if (fig2==0 || ~isvalid(fig2))
    fig2=figure('Position',[100 100 1000 1000*ratio],'visible','on');  
    set(fig2,'PaperSize',12*[1 ratio],'PaperPosition',12*[0 0 1 ratio]);
end
figure(fig2);
clf(fig2);

dn=1; % plot every dn'th velocity vector
quiver(ocean.Xo(1:dn:end),ocean.Yo(1:dn:end),ocean.Uocn(1:dn:end,1:dn:end),ocean.Vocn(1:dn:end,1:dn:end),'k');
hold on;

axis([ocean.Xo(1) ocean.Xo(end) ocean.Yo(1) ocean.Yo(end)]);

colormap('gray'); caxis([0 1]);

title(['Time = ' num2str(Time) ' s']);

for ii = 1:length(Floe)
    hmax(ii) = max(cat(1,Floe(ii).SubFloes.h));
    hmin(ii) = min(cat(1,Floe(ii).SubFloes.h));
end
Maxh = max(hmax);
Minh = min(hmin);
title(['Time = ' num2str(Time) ' s']);
hmax = log10(31);


%% Plot the Floes and Ghost Floes
figure(fig);
for j=1:length(Floe)
    if Floe(j).alive
        poly=intersect(Floe(j).poly,c2_boundary_poly);
        if PERIODIC, poly_ghost=subtract(Floe(j).poly,c2_boundary_poly); end

        
        if stressColor==1
            cFact=min(1,Floe(j).OverlapArea/Floe(j).area); 
            plot(poly,'FaceColor',[1 0 0]*cFact^0.25,'FaceAlpha',0.5);
            if PERIODIC, plot(poly_ghost,'FaceColor','k','FaceAlpha',0.2,'EdgeColor',[1 1 1]*0.4); end

        else
            plot(poly,'FaceColor','k'); 
            if PERIODIC,  plot(poly_ghost,'FaceColor','k','FaceAlpha',0.2,'EdgeColor',[1 1 1]*0.4); end

        end
    end
end

figure(fig2);
plot(c2_boundary_poly,'FaceColor',[0, 0.4470, 0.7410],'FaceAlpha',0.5);
for j=1:length(Floe)
    if Floe(j).alive
        for ii = 1:length(Floe(j).SubFloes)
            poly=intersect(Floe(j).SubFloes(ii).poly,c2_boundary_poly);
            
            if PERIODIC, poly_ghost=subtract(Floe(j).SubFloes(ii).poly,c2_boundary_poly); end
            
            if stressColor==1
                plot(poly,'FaceColor',[1 1 1]*log10(Floe(j).SubFloes(ii).h+1)/hmax,'FaceAlpha',0.75);
                if PERIODIC, plot(poly_ghost,'FaceColor','k','FaceAlpha',0.2,'EdgeColor',[1 1 1]*0.4); end
                
                
            else
                plot(poly,'FaceColor','k');
                if PERIODIC,  plot(poly_ghost,'FaceColor','k','FaceAlpha',0.2,'EdgeColor',[1 1 1]*0.4); end
            end
        end
    end
end

%% Add remaining details to plots

figure(fig);
if ~PERIODIC 
    xb=c2_boundary_poly.Vertices(:,1); xb(end+1)=xb(1);
    yb=c2_boundary_poly.Vertices(:,2); yb(end+1)=yb(1);
    plot(xb,yb, 'k-','linewidth',2);
end


if showCenterOfMass
    plot(cat(1,Floe.Xi),cat(1,Floe.Yi),'k.'); 
end

if showContactPoints
    a=cat(1,Floe(logical(cat(1,Floe.alive))).interactions);
    if ~isempty(a)
        plot(a(:,4),a(:,5),'ro','markersize',2,'MarkerFaceColor','r','lineWidth',1);
    end
end

colormap('gray'); caxis([0 1]);
axis([-Lx-Lx/10 Lx+Lx/10 -Ly-Ly/10 Ly+Ly/10])
xlabel('m');ylabel('m');
set(gca,'Ydir','normal');

drawnow;

figure(fig2);
if ~PERIODIC 
    xb=c2_boundary_poly.Vertices(:,1); xb(end+1)=xb(1);
    yb=c2_boundary_poly.Vertices(:,2); yb(end+1)=yb(1);
    plot(xb,yb, 'k-','linewidth',2);
end

if showCenterOfMass
    plot(cat(1,Floe.Xi),cat(1,Floe.Yi),'k.'); 
end

if showContactPoints
    a=cat(1,Floe(logical(cat(1,Floe.alive))).interactions);
    if ~isempty(a)
        plot(a(:,4),a(:,5),'ro','markersize',2,'MarkerFaceColor','r','lineWidth',1);
    end
end

colormap('gray'); caxis([0 1]);
axis([-Lx-Lx/10 Lx+Lx/10 -Ly-Ly/10 Ly+Ly/10])
xlabel('m');ylabel('m');
set(gca,'Ydir','normal');

drawnow;

mycolors = [0:0.01:1; 0:0.01:1; 0:0.01:1];
mycolors = mycolors';
colormap(mycolors)
caxis([0 1.5])
colorbar

%% Display number of sacked floes

sacked_floes=sum(~cat(1,Floe.alive));
if sacked_floes>0, display(['Sacked floes: ' num2str(sacked_floes)]); end


end