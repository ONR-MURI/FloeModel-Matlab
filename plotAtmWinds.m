function fig = plotAtmWinds(fig,winds,Time)

clf(fig);

dn=4; % Plot every dn'th velocity vector
quiver(winds.X(1:dn:end),winds.Y(1:dn:end),winds.U(1:dn:end,1:dn:end),winds.V(1:dn:end,1:dn:end));

title(['Time = ' num2str(Time) ' s']);

%colormap summar
%colormap('gray'); caxis([0 1]);
axis([min(winds.X) max(winds.X) min(winds.Y) max(winds.Y)])
xlabel('m');ylabel('m');
set(gca,'Ydir','normal');

drawnow;

end
