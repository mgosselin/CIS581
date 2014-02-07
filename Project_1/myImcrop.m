function J = myImcrop(I);
figure(1) = image(I);
uiwait(msgbox('Please select two corners of the region you would like to crop-out of the image','User Prompt','modal'));
[X,Y] = ginput(2);
X = round(X);
Y = round(Y);
J = I(min(Y):max(Y),min(X):max(X),:);
close all;
figure(1) = image(J);

