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
    poly(ii) = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
%     Floe(ii).angles = polyangles(poly.Vertices(:,1),poly.Vertices(:,2));
%     floe = initialize_floe_values(poly, height);
%     Floe0(ii) = floe;
    %vec(ii) = length(Floe(ii).c0);
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

