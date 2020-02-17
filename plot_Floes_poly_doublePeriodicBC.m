function fig=plot_Floes_poly_doublePeriodicBC(fig, Time,Floe,ocean,c2_boundary_poly)

showContactPoints=0;
showCenterOfMass=0;
stressColor=1;



N0=length(Floe);

Lx= max(c2_boundary_poly.Vertices(:,1)); %c2 must be symmetric around x=0 for channel boundary conditions.
Ly= max(c2_boundary_poly.Vertices(:,2)); 

ghostFloeX=[];
ghostFloeY=[];
parent=[];

x=cat(1,Floe.Xi);
alive=cat(1,Floe.alive);

for i=1:length(Floe)
    
%   if alive(i) && (x(i)>Lx-rmax(i)) || (x(i)<-Lx+rmax(i))
   if alive(i) && (max(abs(Floe(i).poly.Vertices(:,1)))>Lx/2)
   
    ghostFloeX=[ghostFloeX  Floe(i)];  
    ghostFloeX(end).poly=translate(Floe(i).poly,[-2*Lx*sign(x(i)) 0]);
    ghostFloeX(end).Xi=Floe(i).Xi-2*Lx*sign(x(i));
    parent=[parent  i];
   
   end
    
    
end

Floe=[Floe ghostFloeX];

y=cat(1,Floe.Yi);
alive=cat(1,Floe.alive);

for i=1:length(Floe)
    
%   if alive(i) && (x(i)>Lx-rmax(i)) || (x(i)<-Lx+rmax(i))
   if alive(i) && (max(abs(Floe(i).poly.Vertices(:,2)))>Ly/2)
   
    ghostFloeY=[ghostFloeY  Floe(i)];  
    ghostFloeY(end).poly=translate(Floe(i).poly,[0 -2*Ly*sign(y(i))]);
    ghostFloeY(end).Yi=Floe(i).Yi-2*Ly*sign(y(i));
    parent=[parent  i];
   
   end    
    
end

Floe=[Floe ghostFloeY];




%ratio=max(c2_boundary(:,1:end-1),[],2)-mean(c2_boundary(:,1:end-1),2); ratio=ratio(2)/ratio(1);

ratio=max(ocean.Yo)/max(ocean.Xo);
if (fig==0 || ~isvalid(fig))
    fig=figure('Position',[100 100 1000 1000*ratio],'visible','on');  
    set(fig,'PaperSize',12*[1 ratio],'PaperPosition',12*[0 0 1 ratio]);
end

clf(fig);

dn=1; % plot every dn'th velocity vector
quiver(ocean.Xo(1:dn:end),ocean.Yo(1:dn:end),ocean.Uocn(1:dn:end,1:dn:end),ocean.Vocn(1:dn:end,1:dn:end));
hold on;

%axis([-1 1 -1 1]*7e4);
%axis([nanmin(c2_boundary(1,:)) nanmax(c2_boundary(1,:)) nanmin(c2_boundary(2,:)) nanmax(c2_boundary(2,:))]);
axis([ocean.Xo(1) ocean.Xo(end) ocean.Yo(1) ocean.Yo(end)]);

colormap('gray'); caxis([0 1]);

title(['Time = ' num2str(Time) ' s']);

for j=1:length(Floe)
    if Floe(j).alive
        poly=intersect(Floe(j).poly,c2_boundary_poly);
        poly_ghost=subtract(Floe(j).poly,c2_boundary_poly);

        if stressColor==1
            cFact=min(1,Floe(j).OverlapArea/Floe(j).area); 
            plot(poly,'FaceColor',[1 0 0]*cFact^0.25,'FaceAlpha',0.5);
            plot(poly_ghost,'FaceColor','k','FaceAlpha',0.2,'EdgeColor',[1 1 1]*0.4);

        else
            plot(poly,'FaceColor','k'); 
            plot(poly_ghost,'FaceColor','k','FaceAlpha',0.2,'EdgeColor',[1 1 1]*0.4);

        end
    end
end

%plot([Floe(logical(cat(1,Floe.alive))).poly]);

xb=c2_boundary_poly.Vertices(:,1); xb(end+1)=xb(1); %elongate the boundaries to inf in the x direction such that the overlap areas and interactions with the x boundary do not count
yb=c2_boundary_poly.Vertices(:,2); yb(end+1)=yb(1);

%plot(xb,yb,'-','linewidth',2,'color', 'k');
%plot([min(ocean.Xo)  max(ocean.Xo)], [min(yb) min(yb)],'-','linewidth',2,'color', 'k')
%plot([min(ocean.Xo)  max(ocean.Xo)], [max(yb) max(yb)],'-','linewidth',2,'color', 'k')

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
axis([min(ocean.Xo) max(ocean.Xo) min(ocean.Yo) max(ocean.Yo)])
xlabel('m');ylabel('m');
set(gca,'Ydir','normal');

drawnow;

sacked_floes=sum(~cat(1,Floe.alive));
if sacked_floes>0, display(['Sacked floes: ' num2str(sacked_floes)]); end


end