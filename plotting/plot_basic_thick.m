function [fig] =plot_basic_thick(fig, Time,Floe,ocean,c2_boundary_poly,Nb)
%This function creates plots of the floe state showing the stress and and
%thickness of the floes
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)
Lx= max(c2_boundary_poly.Vertices(:,1)); %c2 must be symmetric around x=0 for channel boundary conditions.
Ly= max(c2_boundary_poly.Vertices(:,2)); 
live = cat(1,Floe.alive);
Floe(live == 0) = [];

%% Set Up The Plots

ratio=max(ocean.Yo)/max(ocean.Xo);
if (fig==0 || ~isvalid(fig))
    fig=figure('Position',[100 100 1000 1000*ratio],'visible','on');  
    set(fig,'PaperSize',12*[1 ratio],'PaperPosition',12*[0 0 1 ratio]);
end
figure(fig)
clf(fig);

dn=1; % plot every dn'th velocity vector
quiver(ocean.Xo(1:dn:end),ocean.Yo(1:dn:end),ocean.Uocn(1:dn:end,1:dn:end),ocean.Vocn(1:dn:end,1:dn:end),'w');
hold on;

axis([ocean.Xo(1) ocean.Xo(end) ocean.Yo(1) ocean.Yo(end)]);

colormap('gray'); caxis([0 1]);

set(gca,'Color','k')

title(['Time = ' num2str(Time/3600) ' hours'],'fontsize',24);
for ii =1:length(Floe)
    Floe(ii).poly = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
end
h = cat(1,Floe.h);
h(h>3)=3;

%% Plot the Floes and Ghost Floes
for i = 1:length(Floe)
    plot(Floe(i).poly,'FaceColor',[1 1 1]*h(i)/3,'EdgeColor',[1 1 1]*0.2,'FaceAlpha',0.3);
end

if Nb > 0
    plot([Floe(1:Nb).poly],'FaceColor','k','FaceAlpha',0.3,'EdgeColor',[1 1 1]*0.2);
end

set(0,'CurrentFigure',fig);
xb=c2_boundary_poly.Vertices(:,1); xb(end+1)=xb(1);
yb=c2_boundary_poly.Vertices(:,2); yb(end+1)=yb(1);
plot(xb,yb, 'k-','linewidth',2);

%colorbar
%colormap('gray'); caxis([0 1]);
axis([-Lx-Lx/10 Lx+Lx/10 -Lx-Lx/10 Lx+Lx/10])
%xlabel('m');ylabel('m');
set(gca,'Ydir','normal');
% set(gca,'xtick',[])
% set(gca,'ytick',[])

drawnow

warning('on',id)
warning('on',id3)
end