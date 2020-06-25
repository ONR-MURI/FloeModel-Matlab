close all
for jj = 1:length(floe.SubFloes)
    clear x1; clear y1; clear z1;
    x1 = [floe.SubFloes(jj).poly.Vertices(:,1)' floe.SubFloes(jj).poly.Vertices(1,1); floe.SubFloes(jj).poly.Vertices(:,1)' floe.SubFloes(jj).poly.Vertices(1,1)];
    y1 = [floe.SubFloes(jj).poly.Vertices(:,2)' floe.SubFloes(jj).poly.Vertices(1,2); floe.SubFloes(jj).poly.Vertices(:,2)' floe.SubFloes(jj).poly.Vertices(1,2)];
    z1 = [-0.6*floe.SubFloes(jj).h*ones(1,length(x1)); 0.4*floe.SubFloes(jj).h*ones(1,length(x1))];
    h = cat(1,floe.SubFloes.h);
    Maxh = max(h);
    Minh = min(h);
    hmax = log10(1+max(Maxh));
    surface(x1,y1,z1,'FaceColor',[0 1 0]*(log10(1+floe.SubFloes(jj).h)/hmax),'linestyle','none')
    hold on
    h1 = patch(x1(1,:),y1(1,:),z1(1,:),'green','linewidth',0.5);
    set(h1,'FaceColor',[0 1 0]*(log10(1+floe.SubFloes(jj).h)/hmax));
    h2 = patch(x1(1,:),y1(1,:),z1(2,:),'green','linewidth',0.5);
    set(h2,'FaceColor',[0 1 0]*(log10(1+floe.SubFloes(jj).h)/hmax));
end