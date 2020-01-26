function Floe=StochasticFloeRemoval(Floe)

p0=1; % probability of removing a fully overlapping floe during a timestep.

overlap=cat(1,Floe.OverlapArea)./cat(1,Floe.area);

keep=(rand(length(Floe),1) > p0*overlap); 

Floe=Floe(keep);

Nkilled=sum(~keep);

if Nkilled > 0, disp(['Fractured floes:' num2str(Nkilled)]); end


end