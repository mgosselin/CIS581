function [p] = feat_desc(im,y,x);
%%

perim = 20;
[row,col] = size(im);

rowstart = 1+perim;
rowend = row+perim;
colstart = 1+perim;
colend = col+perim;

% create some black space around the picture
CMG = zeros(row+2*perim,col+2*perim);
CMG(rowstart:rowend,colstart:colend) = im;

% do all 4 edges first
for i = 1:perim
    % left edge
    CMG(rowstart:rowend,i) = CMG(rowstart:rowend,1+2*perim-i);
    % right edge
    CMG(rowstart:rowend,1+col+2*perim-i) = CMG(rowstart:rowend,i+col);
    % top edge
    CMG(i,colstart:colstart+col) = CMG(1+2*perim-i,colstart:colstart+col);
    % bottom edge
    CMG(1+row+2*perim-i,colstart:colstart+col) = CMG(i+row,colstart:colstart+col);
end

% do all the corners second
for i = 1:perim
    for j = 1:perim
        % top left
        CMG(i,j) = CMG(1+2*perim-i,1+2*perim-j);
        % bottom left
        CMG(1+row+2*perim-i,j) = CMG(row+i,1+2*perim-j);
        % top right
        CMG(rowstart-i,colend+j) = CMG(perim+i,colstart+col-j);
        % bottom right
        CMG(1+row+2*perim-i,1+col+perim*2-j) = CMG(i+row,j+col);
    end
end

%%

G = fspecial('gaussian',[3 3],0.60);
CMG = conv2(CMG,G,'same');

% imshow(CMG/255);

%%

N = length(y);

for j = 1:N
    pix = CMG((y(j)):(y(j)+40-1),(x(j)):(x(j)+40-1));
    pix = pix(:);
    pix = pix(1:25:end);
    pix = double(pix);
    pix = pix - mean(pix);
    pix = pix/std(pix);
    p(:,j) = pix;
end

%%