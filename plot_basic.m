function [fig] =plot_basic(fig, Time,Floe,ocean,c2_boundary_poly)
%This function creates plots of the floe state showing the stress and and
%thickness of the floes
Lx= max(c2_boundary_poly.Vertices(:,1)); %c2 must be symmetric around x=0 for channel boundary conditions.
Ly= max(c2_boundary_poly.Vertices(:,2)); 

%% Set Up The Plots

ratio=max(ocean.Yo)/max(ocean.Xo);
if (fig==0 || ~isvalid(fig))
    fig=figure('Position',[100 100 1000 1000*ratio],'visible','on');  
    set(fig,'PaperSize',12*[1 ratio],'PaperPosition',12*[0 0 1 ratio]);
end
figure(fig)
clf(fig);

dn=1; % plot every dn'th velocity vector
quiver(ocean.Xo(1:dn:end),ocean.Yo(1:dn:end),ocean.Uocn(1:dn:end,1:dn:end),ocean.Vocn(1:dn:end,1:dn:end));
hold on;

axis([ocean.Xo(1) ocean.Xo(end) ocean.Yo(1) ocean.Yo(end)]);

colormap('gray'); caxis([0 1]);

title(['Time = ' num2str(Time) ' s']);


%% Plot the Floes and Ghost Floes
plot([Floe.poly],'FaceColor','k','FaceAlpha',0.3,'EdgeColor',[1 1 1]*0.2);

set(0,'CurrentFigure',fig);
xb=c2_boundary_poly.Vertices(:,1); xb(end+1)=xb(1);
yb=c2_boundary_poly.Vertices(:,2); yb(end+1)=yb(1);
plot(xb,yb, 'k-','linewidth',2);


colormap('gray'); caxis([0 1]);
axis([-Lx-Lx/10 Lx+Lx/10 -Ly-Ly/10 Ly+Ly/10])
xlabel('m');ylabel('m');
set(gca,'Ydir','normal');

drawnow
end