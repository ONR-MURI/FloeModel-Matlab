%%Amp3
close all
Yequispaced = - (-0.25:0.01:0.25);
[Xeq, Yeq] = meshgrid(x4,Yequispaced);
[Xcheb,Ycheb] = meshgrid(x4,Y);
Ueq = interp2(Xcheb,Ycheb,U(:,:,3),Xeq,Yeq);
Veq = interp2(Xcheb,Ycheb,V(:,:,3),Xeq,Yeq);
quiver(x4(1:incx:Nx4),Yequispaced,Ueq(:,1:incx:Nx4),Veq(:,1:incx:Nx4),'k','linewidth',1);
ylim([-0.025 0.025])
xlim([0 0.5])
starty3 = [0 0 -0.005 -0.01 0.005 0.01];
startx3 = [0.1855 0.34 0.5 0.5 0.01 0.01];
streamline(x4,Yequispaced,Ueq,-Veq,startx3,starty3);

alabamahi = shaperead('usastatehi', 'UseGeoCoords', true,...
            'Selector',{@(name) strcmpi(name,'Alabama'), 'Name'});
bama = polyshape(alabamahi.Lon,alabamahi.Lat);
plot(bama)
%% Clean up
polynew = simplify(polynew);
polyout = sortregions(polynew,'area','descend');
R = regions(polyout);
polynew = R(1);
poly = rmholes(polynew);
States(i).poly = poly;
plot(States(i).poly)
i = i+1;