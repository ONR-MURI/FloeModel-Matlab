
close all
clear all


A=rgb2gray(imread('DSC01029_v2.JPG'));

A=imbinarize(A,0.85);
A = imfill(A,'holes');

A = bwareafilt(A,[2e1 5e4],8);

A = imgaussfilt(double(A),1);

A = imbinarize(A,0.5);

A = bwareafilt(A,[2e1 5e4],8);

figure; imagesc(A);


CC = bwconncomp(A,8);
FloeShapes = regionprops(CC,'Image','Centroid');

save('FloeShapes.mat','FloeShapes');