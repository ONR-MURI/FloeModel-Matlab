close all
floe = Floe(2);
x=cat(1,Floe.Xi);
y=cat(1,Floe.Yi);
N0 = length(Floe);
floenew = [];
for kk = 1:3000
    floe = FloeSimplify(floe);
    inter = floe.interactions(:,1);
    inter = inter(inter<=N0);
    Xi = x(inter);
    Yi = y(inter);
    d = sqrt((Xi-floe.Xi).^2-(Yi-floe.Yi).^2);
    
    angles=cat(1,floe.angles);
    Anorm = 180-360/length(angles);
    keep1=(rand(length(angles),1)>angles/Anorm);
    da = ones(length(angles),1);
    for jj = 1:length(angles)
        if ~keep1(jj)
            d2 = sqrt((Xi-floe.c_alpha(1,jj)).^2-(Yi-floe.c_alpha(2,jj)).^2);
            if max(d2-d)>0
                da(jj) = 1;
            end
        end
    end
    
    keep = logical(keep1 + da);
    
    if sum(keep)>0
        [~,I] = min(angles(keep));
        fracturedFloes = frac_corner(floe,I);
        floenew=[floenew fracturedFloes];
    end
    Areas = cat(1,fracturedFloes.area);
    [~,I] = max(Areas);
    fracturedFloes(I).interactions = floe.interactions;
    floe = fracturedFloes(I);
    poly(kk) = polyshape(floe.c_alpha'+[floe.Xi floe.Yi]);
%     plot(poly)
end
for ii = 0:50
    fig = plot(poly(ii*5+1));
    xlim([6e4 6.6e4])
    ylim([2.8e4 3.6e4])
    drawnow
    saveas(fig,['./figs/' num2str(ii+1,'%03.f') '.jpg'],'jpg');
end