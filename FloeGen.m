close all
clear all
A=rgb2gray(imread('full_image.jpg'));
A = imfill(A,'holes');
A = imbinarize(A,0.8);
A = imgaussfilt(double(A),0.5);
A = imbinarize(double(A),0.8);
%A = bwareafilt(A,[50 1e5],8);
CC = bwconncomp(A,8);
FloeShapes = regionprops(CC,'Image','Centroid');
figure; imagesc(A);
save('FloeShapes.mat','FloeShapes');