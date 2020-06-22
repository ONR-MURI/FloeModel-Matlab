
close all
clear all

%%
A0=rgb2gray(imread('full_image.jpg'));
A=A0;

A(A0<10)=0; A(A0>=10)=1;
A = imfill(A,'holes');
%A = imbinarize(A,0.8);
A = imgaussfilt(double(A),0.5);
A = imbinarize(double(A),0.8);
A = bwareafilt(A,[50 1e10],8);

CC = bwconncomp(A,8);
FloeShapes = regionprops(CC,'Image','Centroid');

figure; 

subplot(1,2,1);
imagesc(A0);
axis xy;
colormap('gray')
subplot(1,2,2);
set(gca,'Color','k')
imagesc(A)
axis xy

%%
Floe=initialize_Floe(FloeShapes);

%%
figure; 

subplot(1,2,1);
imagesc(A0);
axis xy;
ylim([50 510]);
xlim([1 470])
colormap('gray')
subplot(1,2,2);
plot([Floe.poly],'FaceColor','k')
xlim([-0.9 1]*6e4);
ylim([-1 0.9]*6e4);

%%
RIDGING=false; 

FRACTURES=false;

PERIODIC=true;

SUBFLOES = false;

PACKING = false;


%Define ocean currents
Lx=14e4; Ly=20e4;
[ocean, c2_boundary]=initialize_ocean_Gyre(1e4, Lx, Ly,2e4);
c2_boundary_poly=polyshape(c2_boundary(1,:),c2_boundary(2,:));

ocean.Uocn=0.2*ones(size(ocean.Uocn));
ocean.Vocn=0.2*ones(size(ocean.Vocn));

%Define 10m winds
winds=[0 0];

%Initialize Floe state
height.mean = 2;
height.delta = 0; %max difference between a flow thickness and the mean floe value

c=0.9; % could be a vector
Floe1 = initialize_concentration(c,c2_boundary,SUBFLOES, height, 10000);

%%
fig=figure; 
subplot(1,2,1);
imagesc(A0);
axis xy;
ylim([50 510]);
xlim([1 470])
set(gca,'Xtick',[],'Ytick',[])
colormap('gray')

subplot(1,2,2);
plot([Floe1.poly Floe.poly],'FaceColor','w','FaceAlpha',1);
set(gca,'Color','k')
xlim([-0.9 1]*6e4);
ylim([-1 0.9]*6e4);
set(gca,'Xtick',[],'Ytick',[])
set(fig,'papersize',[12 5],'paperposition',[0 0 12 5]);
saveas(fig,'satelliteFloes.pdf','pdf');

