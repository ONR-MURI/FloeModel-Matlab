function [fig] =plot_Arctic(fig, Time,Floe,ocean,c2_boundary,Nb)
%This function creates plots of the floe state showing the stress and and
%thickness of the floes

%% Set Up The Plots

% ratio=max(ocean.Yo)/max(ocean.Xo);
ratio=max(ocean.Yo)/max(ocean.Xo);
if (fig==0 || ~isvalid(fig))
    fig=figure('Position',[100 100 1000 1000*ratio],'visible','on');  
    set(fig,'PaperSize',12*[1 ratio],'PaperPosition',12*[0 0 1 ratio]);
end
figure(fig)
clf(fig);

dn=1; % plot every dn'th velocity vector
% quiver(ocean.Xo(1:dn:end),ocean.Yo(1:dn:end),ocean.Uocn(1:dn:end,1:dn:end),ocean.Vocn(1:dn:end,1:dn:end));
hold on;

axis([ocean.Xo(1) ocean.Xo(end) ocean.Yo(1) ocean.Yo(end)]);

colormap('gray'); caxis([0 1]);

% title(['Time = ' num2str(Time) ' s'],'fontsize',24);
for ii =1:length(Floe)
    Floe(ii).poly = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
    Stress(ii) = max(abs(eig(Floe(ii).Stress)));
end

A = cat(1,Floe.area);

plot([Floe(1+Nb:end).poly],'FaceColor','k','FaceAlpha',0.3,'EdgeColor',[1 1 1]*0.2);
if Nb > 0
    plot([Floe(1:Nb).poly],'FaceColor',[0 0.2 0],'FaceAlpha',0.75,'EdgeAlpha',0);
end

load('ArcticPolyshapes.mat')
Nx = 40; Ny = 40;
xc = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
yc = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
Xc = (xc(1:end-1)+xc(2:end))/2; Yc = (yc(1:end-1)+yc(2:end))/2;
[Xcg,Ycg] = meshgrid(Xc,Yc);
Xcg = Xcg(:); Ycg = Ycg(:);
[Uwinds, Vwinds] = ArcticWindGenerator;
in = inpolygon(Xcg,Ycg,polyAU.Vertices(:,1),polyAU.Vertices(:,2));
u = Uwinds(Xcg(:),Ycg(:)); v = Vwinds(Xcg(:),Ycg(:));
quiver(Xcg(~in),Ycg(~in),u(~in),v(~in));

% xx = 1; xx(1) =[1 2]
set(0,'CurrentFigure',fig);
xb=c2_boundary(1,:);
yb=c2_boundary(2,:);
% plot(xb,yb, 'r-','linewidth',2);
%title(['Time = ' num2str(Time) ' s']);

colormap('gray'); caxis([0 1]);
% axis([-Lx Lx -Ly Ly])
%xlabel('m');ylabel('m');
set(gca,'Ydir','normal');
% set(gca,'xtick',[])
% set(gca,'ytick',[])
ylim([min(yb) max(yb)])
xlim([min(xb) max(xb)])

drawnow
end