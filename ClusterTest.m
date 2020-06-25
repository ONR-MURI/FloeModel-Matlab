clear all; close all;

SUBFLOES = false;

%Define ocean currents
[ocean, c2_boundary]=initialize_ocean_Gyre(1e4, 2e5, 1e5,4e3);

%Initialize Floe state
height.mean = 2;
height.delta = 0.5; %max difference between a flow thickness and the mean floe value

target_concentration=1; % could be a vector
N = 200:100:500;
N0 = fix(1200^(1/3));
for ii = 1:length(N)
    Floe0 = [];
    % Floe = initialize_concentration(target_concentration,c2_boundary,SUBFLOES, height, 50);
    [Floe]= FloeGeneratorConcentration(Floe0,c2_boundary,target_concentration,N(ii),SUBFLOES,height,N0);
    
    figure(1)
    A = sqrt(cat(1,Floe.area));
    Floe(A<sqrt(3500)) = [];
    A2 = histogram(sqrt(A));
    bins = (A2.BinEdges(1:end-1)+A2.BinEdges(2:end))/2;
    figure(2)
    subplot(2,2,ii)
    plot(bins(A2.Values>0),A2.Values(A2.Values>0),'ok','linewidth',2)
%     histogram(A)
    set(gca,'XScale','log','YScale','log')
    set(gca,'fontsize',18)
    title(['Total Floes = ' num2str(length(Floe))]);
    ylabel('Number of Floes','fontsize',36)
    xlabel('Floe Size','fontsize',36)
%     figure
%     plot([Floe.poly])
end