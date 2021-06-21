if FRACTURES
    overlapArea=cat(1,Floe.OverlapArea)./cat(1,Floe.area);
%     keep=rand(length(Floe),1)<5*overlapArea/max(overlapArea);
    fracturedFloes=fracture_floe_normal(Floe(~keep),[1 1]);
    if ~isempty(fracturedFloes), fracturedFloes=rmfield(fracturedFloes,'potentialInteractions');
        %             Floe=[Floe(keep) fracturedFloes];
    end
    
end
height.mean = 2;
height.delta = 0;
for ii =1:length(Floe)
%     poly = polyshape(b{ii});
%     poly(ii) = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
%     Floe(ii).angles = polyangles(poly.Vertices(:,1),poly.Vertices(:,2));
%     floe = initialize_floe_values(poly, height);
%     Floe0(ii) = floe;
    vec(ii) = length(Floe(ii).c0);
end
Floe = Floe0;
save('Floe0.mat','Floe')
max(vec)
[i1, i2] = max(vec)

Time = 0;
fig = 0;
for ii = 1:45
    load(['./FloesNares/Floe' num2str(ii,'%07.f') '.mat'],'Floe');
    [fig] =plot_Nares(fig, Time,Floe,ocean,c2_boundary_poly,Nb);
    ax = gca;
    exportgraphics(ax,['./Nares/' num2str(ii,'%03.f') '.jpg'],'Resolution',300)
end

close all;
dt = 30;  
Nb = 2;
fig = 0;
Lx=2e5; Ly=1e6;
x=[-1 -1 1 1 -1]*Lx/2; y=[-1 1 1 -1 -1]*Ly/2;
c2_boundary = [x; y]; c2_boundary_poly=polyshape(c2_boundary');
for jj = 1:168
    load(['./FloesNares/Floe' num2str(jj,'%07.f') '.mat'],'Floe');
    Time = dt*(jj-1)*100;
    clear poly
    for ii =1:length(Floe)
        poly(ii)= polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
        %Floe(ii).poly = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
    end
    plot(poly)
%     drawnow
    fig = plot_Nares(fig, Time,Floe,ocean,c2_boundary_poly,Nb);
    saveas(fig,['./Floes/' num2str(jj,'%03.f') '.jpg'],'jpg');
end

r = 10;
for m = 1:100
for n = 1:5
    x = r*(2*rand(10^n,1) - 1);
    y = r*(2*rand(10^n,1) - 1);
    in = inpolygon(x,y,poly.Vertices(:,1),poly.Vertices(:,2));
    A = sum(in)/length(x)*4*r^2;
    error(n) = 100*abs(area(poly)-A)/(area(poly));
end
E1(m,:) = error;
end
n = 1:5;
Error = mean(E1,1);
semilogx(10.^n,Error,'xk','linewidth',2)

figure
for n = 1:5
    dx = 20/sqrt(10^n);
    x = -10:dx:10;
    [xx,yy]=meshgrid(x,x);
    in = inpolygon(xx(:),yy(:),xc,yc);
    A = sum(in)/length(xx(:))*400;
    error2(n) = (100*pi-A)/(100*pi);
end
n = 1:5;
plot(10.^n,error2)
