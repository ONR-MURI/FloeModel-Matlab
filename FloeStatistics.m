load FloeStats
% FloeStats = FloeStats
close all
time = 10*(0:length(FloeStats)-1)*nDTOut;
%edges = 6*10.^(1:fix(log10(sqrt(Amax)))+1);
figure
plot(time,cat(1,FloeStats.Num),'k','linewidth',2)
set(gca,'fontsize',18)
ylabel('Number of Floes','fontsize',36)
xlabel('Time','fontsize',36)
figure
figure
rho_ice = 920;
clear m
clear Vol
for ii = 2:length(time)
    figure(2)
    Vol(ii) = sum(FloeStats(ii).A(:,1).*FloeStats(ii).h(:,1))/nDTOut+FloeStats(ii).DissolvedMass/rho_ice;
    plot(sqrt(FloeStats(ii).A),FloeStats(ii).h,'ok','linewidth',2)
    figure(3)
    FloeStats(ii).A(sqrt(FloeStats(ii).A)<1e3) = [];
    A = histogram(sqrt(FloeStats(ii).A));
    bins = (A.BinEdges(1:end-1)+A.BinEdges(2:end))/2;
%     figure(3)
%     h = hisogram(FloeStats(ii).h/nDTOut,edges);
    p = polyfit(log10(bins(A.Values>0)),log10(A.Values(A.Values>0)), 1);
    m(ii) = p(1);
end
figure
plot(time,m,'xk','linewidth',2)
set(gca,'fontsize',18)
ylabel('Power Law Exponent','fontsize',36)
xlabel('Time','fontsize',36)
figure
A2 = histogram(sqrt(FloeStats(ii).A));
bins = (A2.BinEdges(1:end-1)+A2.BinEdges(2:end))/2;
figure
plot(bins(A2.Values>0),A2.Values(A2.Values>0),'ok','linewidth',2)
set(gca,'XScale','log','YScale','log')
set(gca,'fontsize',18)
ylabel('Number of Floes','fontsize',36)
xlabel('Floe Size','fontsize',36)
%loglog(A2,A1,'ok','linewidth',2)
figure
h2 = histogram(FloeStats(ii).h);
bins = (h2.BinEdges(1:end-1)+h2.BinEdges(2:end))/2;
figure
plot(bins(h2.Values>0),h2.Values(h2.Values>0),'ok','linewidth',2)
set(gca,'XScale','log','YScale','log')
set(gca,'fontsize',18)
ylabel('Number of Floes','fontsize',36)
xlabel('Floe Thickness','fontsize',36)
figure
plot(time(2:end),Vol(2:end),'k','linewidth',2)
set(gca,'fontsize',18)
ylabel('Volume of Ice','fontsize',36)
xlabel('Time','fontsize',36)
figure(2)
set(gca,'XScale','log','YScale','log')
set(gca,'fontsize',18)
xlabel('Floe size','fontsize',36)
ylabel('Floe Thickness','fontsize',36)
ylim([0.1 10])