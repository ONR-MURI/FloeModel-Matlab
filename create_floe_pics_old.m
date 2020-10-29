
close all
clear all


A=rgb2gray(imread('full_image.jpg'));

A=imbinarize(A,0.75);
A = imfill(A,'holes');
A = bwareafilt(A,[5 1e5],8);

A = imgaussfilt(double(A),2);
A = imbinarize(A,0.65);
A = bwareafilt(A,[50 1e5],8);

CC = bwconncomp(A,8);
FloeShapes = regionprops(CC,'Image','Centroid');

figure; imagesc(A);

save('FloeShapes.mat','FloeShapes');