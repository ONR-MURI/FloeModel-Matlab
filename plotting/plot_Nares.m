function [fig] =plot_Nares(fig, Time,Floe,ocean,c2_boundary_poly,Nb,eularian_data)
%This function creates plots of the floe state showing the stress and and
%thickness of the floes
Lx= max(c2_boundary_poly.Vertices(:,1)); %c2 must be symmetric around x=0 for channel boundary conditions.
Ly= max(c2_boundary_poly.Vertices(:,2)); 

%% Set Up The Plots

% ratio=Ly/Lx;
ratio = abs(-4e5/2)/(Lx); %Ly/Lx;%
if (fig==0 || ~isvalid(fig))
    fig=figure('Position',[10 10 200 200*ratio],'visible','on');  
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

plot([Floe(1+Nb:end).poly],'FaceColor','k','FaceAlpha',0.3,'EdgeColor',[1 1 1]*0.2,'linewidth',0.2);
if Nb > 0
    plot([Floe(1:Nb).poly],'FaceColor',[0 0.2 0],'FaceAlpha',0.75,'EdgeColor',[1 1 1]*0.2);
end

[Ny, Nx]=size(eularian_data.u);%fix(Nx*LyO/LxO);
c2_boundary =c2_boundary_poly.Vertices';
xc = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
yc = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
Xc = (xc(1:end-1)+xc(2:end))/2; Yc = -(yc(1:end-1)+yc(2:end))/2;
inP= zeros(Ny,Nx); 
%for jj = 1:Nb
[Xcg,Ycg] = meshgrid(Xc,Yc);
for jj = 1:Nb
    in = inpolygon(Xcg(:),Ycg(:),Floe(jj).poly.Vertices(:,1),Floe(jj).poly.Vertices(:,2));
    inP(:) = inP(:) + in;
end
q1 = quiver(Xcg(inP==0),Ycg(inP==0),eularian_data.u(inP==0)*50000,eularian_data.v(inP==0)*50000,'color',[0.2549 0.4118 0.8824],'autoscale','off');
x1 = -9.45e4; y1 = -7.9e4; u = 0.75; v = 0;
% quiver(x1,y1,u*50000,v*50000,'color',[0.2549 0.4118 0.8824],'autoscale','off','linewidth',1.5,'maxheadsize',0.75);

set(0,'CurrentFigure',fig);
xb=c2_boundary_poly.Vertices(:,1); xb(end+1)=xb(1);
yb=c2_boundary_poly.Vertices(:,2); yb(end+1)=yb(1);
%plot(xb,yb, 'k-','linewidth',2);
%title(['Time = ' num2str(Time) ' s']);
% dim = [.125 .35 .425 .07];
% Time = round(10*Time/(24*3600))/10;
% %str = 'Time = ' num2str(Time) ' days';
% %annotation('textbox',dim,'String',['Time = ' num2str(Time) ' days'],'fontsize',16,'Color',[1, 1 ,1]);
% annotation('textbox',dim,'String',['Day ' num2str(Time)],'fontsize',14,'Color',[0, 0 ,0],'linestyle','none');%,'FitBoxToText','on');
% xu = [-9.65e4 -4.8e4]; yu = [-0.45e5 -0.45e5];
% plot(xu,yu,'k','linewidth',0.5)
% 
% %x = [0.13 0.32]; y = [0.1275 0.1275];
% %annotation('line',x,y,'linewidth',1.25);
% %annotation('doublearrow',x,y,'head1width',3,'head2width',3,'linewidth',1.25);
% 
% dim2 = [.15 .115 .425 .07];
% annotation('textbox',dim2,'String',['25 km'],'fontsize',8,'Color',[0, 0 ,0],'linestyle','none');%,'FitBoxToText','on');
% 
% xr = [-8.75e4 -6.25e4]; yr = [-1.2e5 -1.2e5];
% plot(xr,yr,'k','linewidth',2)
% 
% dim2 = [.127 .24 .425 .07];
% annotation('textbox',dim2,'String',['0.75 m/s'],'fontsize',8,'Color',[0, 0 ,0],'linestyle','none');%,'FitBoxToText','on');
% 
% dim = [.127 .85 .425 .07];
% % annotation('textbox',dim,'String',['A)'],'fontsize',18,'Color',[0, 0 ,0],'linestyle','none');

colormap('gray'); caxis([0 1]);
axis([-Lx Lx -Ly Ly])
%xlabel('m');ylabel('m');
set(gca,'Ydir','normal');
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
ylim([-2.5e5 1.5e5])
%legend(q1([1]),'Ice Velocity')

drawnow
end