function plot_potentialInteractions(floe)


figure; 

plot(floe.poly,'FaceColor','k'); hold on; 

for k=1:length(floe.potentialInteractions),
    
    plot(floe.potentialInteractions(k).c); 

end


end

