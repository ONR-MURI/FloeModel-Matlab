fig = 0;

close all
clear U
clear V
clear ocean
for ii = 2:91
    Time = dt*nDTOut*(ii-1);
%     ocean.Uocn=ocean.U*cos(ocean.fCoriolis*Time)*ones(size(ocean.Vocn));  
%     ocean.Vocn=-ocean.U*sin(ocean.fCoriolis*Time)*ones(size(ocean.Uocn));
    load(['./FloesCircle90/Floe' num2str(ii,'%07.f') '.mat']);
    [eularian_data] = calc_eulerian_data(Floe,Nx,Ny,c2_boundary,PERIODIC);    
    u =eularian_data.u;
    v=eularian_data.v;
    U(ii) = u(2,2);
    V(ii) = v(2,2);
    [fig]=plot_basic(fig, Time,Floe, ocean, c2_boundary_poly);
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    t = Time/3600;
    axes('Position',[.75 .75 .2 .2])
    box on
    hold on
    xlim([-0.6 0.6])
    ylim([-0.6 0.6])
    quiver(0,0,U(ii),V(ii),'color',[0 0 0],'autoscale',0,'linewidth',1.5)
    quiver(0,0,ocean.Uocn(1), ocean.Vocn(1),'color',[1 0 0],'autoscale',0,'linewidth',1.5)
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    plot(U(1:ii-1),V(1:ii-1),'ok','markerfacecolor','k')
    legend('Ice','Ocn')
    title(['Time = ' num2str(fix(t)) ' Hours'],'fontsize',24);
    set(0,'CurrentFigure',fig);
    saveas(fig,['./FigsVor/' num2str(ii,'%03.f') '.jpg'],'jpg');
end