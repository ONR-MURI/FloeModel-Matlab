Nx=6; 
Lx=6e5; Ly=1.5e5;
x=[-1 -1 1 1 -1]*Lx;
y=[-1 1 1 -1 -1]*Ly;
c2_boundary = [x; y];
xc = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
Xc =  (xc(1:end-1)+xc(2:end))/2; Yc = 0;
dx = xc(2)-xc(1); y = [-65e4 65e4];
rho_ice=920; % kg/m3
dt = 1000;

min_floe_size = 4*Lx*Ly/20000;%4*Lx*Ly/25000;
xx=[-1 -1 1 1 -1]*dx; 
yy=[-1 1 1 -1 -1]*Ly;
c2_boundary = [xx; yy];
target_concentration = 1;
height.mean = 2;
height.delta = 0; %max difference between a flow thickness and the mean floe value
FloeNum = [75 100 125 150 250];
% mass = 0.5*6.6238e+13*ones(1,20);
jj = 1;
% for jj = 1:5
[Floe0, Nb] = initial_concentration(c2_boundary,target_concentration,height,FloeNum(jj),min_floe_size);
while length(Floe0)<FloeNum(jj)
    [Floe0, Nb] = initial_concentration(c2_boundary,target_concentration,height,FloeNum(jj),min_floe_size);
end
boundary0 = polyshape(c2_boundary');
boundary = boundary0;
N1(jj) = length(Floe0);
Floe = Floe0;
for ii = 2:20
%         ii = 10;
    close all
    A(ii-1) = area(boundary);
    load(['./FloesSubNew/Floe' num2str(ii,'%07.f') '.mat'],'U','SigXXa','mass');
    %     load(['./FloesNew/Floe' num2str(ii,'%07.f') '.mat'],'U','mass','SigXXa');
    Ua = (mass(1:2:length(mass)-1).*U(1:2:length(mass)-1)+mass(2:2:length(mass)).*U(2:2:length(mass)))./(mass(1:2:length(mass)-1)+mass(2:2:length(mass)));
%         load(['./FloesSpice/Macro3' num2str(ii-1,'%07.f') '.mat']);
    load(['./FloesSubNew/Floe' num2str(ii+1,'%07.f') '.mat'],'U','SigXXa','mass');
    SigXX = (mass(2:2:length(mass)-2).*SigXXa(2:2:length(mass)-2)+mass(3:2:length(mass)-1).*SigXXa(2:2:length(mass)-1))./(mass(2:2:length(mass)-2)+mass(3:2:length(mass)-1));
    [Sig,~,Floe,Hist, bound] = Subzero_Microscale(Ua(2:3),6,dt,Floe,boundary,SigXX);
    %     SigHist(jj,:) = Hist;
    %     N2(jj) = length(Floe);
    boundary = polyshape(bound');
    save(['./FloesSpice/Macro4' num2str(ii,'%07.f') '.mat'],'Floe', 'Sig','boundary');
end