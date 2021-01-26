function [fig] =plot_Nares(fig, Time,Floe,ocean,c2_boundary_poly,Nb)
%This function creates plots of the floe state showing the stress and and
%thickness of the floes
Lx= max(c2_boundary_poly.Vertices(:,1)); %c2 must be symmetric around x=0 for channel boundary conditions.
Ly= max(c2_boundary_poly.Vertices(:,2)); 

%% Set Up The Plots

% ratio=max(ocean.Yo)/max(ocean.Xo);
ratio = abs(-3e5)/(Lx); %Ly/Lx;%
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
    poly(ii) = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
    Stress(ii) = max(abs(eig(Floe(ii).Stress)));
end
% xx = 1; xx(1) =[1 2];
A = cat(1,Floe.area);
Stress7 = Stress(A<1e8); poly7 = poly(A<1e8);
Stress8 = Stress; Stress8(A<1e8) = 0; Stress8(A>1e9) = 0; 
poly8 = poly; 
poly8(Stress8==0) = []; Stress8(Stress8==0) = [];
Stress9 = Stress(A>1e9); poly9= poly(A>1e9);
%% Plot the Floes and Ghost Floes
% plot([Floe.poly],'FaceColor','k','FaceAlpha',0.3,'EdgeColor',[1 1 1]*0.2);
% if max(Stress)>0
%     for i = 1+Nb:length(Floe)
%         plot(Floe(i).poly,'FaceColor',[1 0 0]*Stress(i)/max(Stress),'EdgeColor',[1 1 1]*0.2);
%     end
% else
%     plot([Floe.poly],'FaceColor','r','FaceAlpha',0.3,'EdgeColor',[1 1 1]*0.2);
% end

if max(Stress7)>0
    for i = 1:length(poly7)
        plot(poly7(i),'FaceColor',[1 0 0]*Stress7(i)/max(Stress7),'EdgeColor',[1 1 1]*0.2);
    end
else
    plot(poly7,'FaceColor','r','FaceAlpha',0.3,'EdgeColor',[1 1 1]*0.2);
end
if max(Stress8)>0
    for i = 1:length(poly8)
        plot(poly8(i),'FaceColor',[1 0 0]*Stress8(i)/max(Stress8),'EdgeColor',[1 1 1]*0.2);
    end
else
    plot(poly8,'FaceColor','r','FaceAlpha',0.3,'EdgeColor',[1 1 1]*0.2);
end
if max(Stress9)>0
    for i = 1:length(poly9)
        plot(poly9(i),'FaceColor',[1 0 0]*Stress9(i)/max(Stress9),'EdgeColor',[1 1 1]*0.2);
    end
else
    plot(poly9,'FaceColor','r','FaceAlpha',0.3,'EdgeColor',[1 1 1]*0.2);
end
if Nb > 0
    plot(poly(1:Nb),'FaceColor','k','FaceAlpha',0.3,'EdgeColor',[1 1 1]*0.2);
end

set(0,'CurrentFigure',fig);
xb=c2_boundary_poly.Vertices(:,1); xb(end+1)=xb(1);
yb=c2_boundary_poly.Vertices(:,2); yb(end+1)=yb(1);
plot(xb,yb, 'k-','linewidth',2);


colormap('gray'); caxis([0 1]);
axis([-Lx Lx -Ly Ly])
%xlabel('m');ylabel('m');
set(gca,'Ydir','normal');
set(gca,'xtick',[])
set(gca,'ytick',[])
ylim([-3.5e5 -0.5e5])

drawnow
end