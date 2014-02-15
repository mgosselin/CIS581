home;
clear;
tic;

img_input = cell(3,1);
img_input{1} = imread('quad08.jpg');
img_input{2} = imread('quad09.jpg');
% img_input{3} = imread('quad10.jpg');

img_mosaic = mymosaic(img_input);
imshow(img_mosaic);
toc;