time = 10*(1:im_num-1)*50;
edges = 10.^(3:fix(log10(Amax))+1);
figure
figure
plot(time,NumFloes(1:im_num-1),'k','linewidth',2)
set(gca,'fontsize',18)
ylabel('Number of Floes','fontsize',36)
xlabel('Time','fontsize',36)
for ii = 1:length(time)
    figure(1)
    A = histogram(FloeStats(ii).A/nDTOut,edges);
%     figure(2)
%     h = hisogram(FloeStats(ii).h/nDTOut,edges);
    p = polyfit(log10(edges(A.Values>0)),log10(A.Values(A.Values>0)), 1);
    m(ii) = p(1);
end
figure
plot(time,m,'xk','linewidth',2)
set(gca,'fontsize',18)
ylabel('Power Law Exponent','fontsize',36)
xlabel('Time','fontsize',36)
figure
A = cat(1,Floe.area);
histogram(A)
set(gca,'XScale','log','YScale','log')
set(gca,'fontsize',18)
ylabel('Number of Floes','fontsize',36)
xlabel('Floe Size','fontsize',36)
%loglog(A2,A1,'ok','linewidth',2)
