function [y,x,rmax] = anms(cimg,max_pts)
%%
% local maximum calculation, to get corner candidates

perim = 1;
[row,col] = size(cimg);

rowstart = 1+perim;
rowend = row+perim;
colstart = 1+perim;
colend = col+perim;

% create some black space around the picture
CMG = zeros(row+2*perim,col+2*perim);
CMG(rowstart:rowend,colstart:colend) = cimg;

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

[Y,X] = ndgrid(1:row,1:col);
Y = Y(:); 
X = X(:);
M = length(Y);
CORx = NaN(M,1);
CORy = NaN(M,1);

for k = 1:M
    if cimg(Y(k),X(k)) > CMG(Y(k),X(k)+1) &&...
       cimg(Y(k),X(k)) > CMG(Y(k)+1,X(k)+2) &&...
       cimg(Y(k),X(k)) > CMG(Y(k)+2,X(k)+1) &&...
       cimg(Y(k),X(k)) > CMG(Y(k)+1,X(k)) &&...
       cimg(Y(k),X(k)) > CMG(Y(k),X(k)) &&...
       cimg(Y(k),X(k)) > CMG(Y(k),X(k)+2) &&...
       cimg(Y(k),X(k)) > CMG(Y(k)+2,X(k)) &&...
       cimg(Y(k),X(k)) > CMG(Y(k)+2,X(k)+2)
    CORx(k) = X(k);
    CORy(k) = Y(k);
    end
end

%combine the x and y pixel coordinates for the corner candidates
CORx = CORx(~isnan(CORx));
CORy = CORy(~isnan(CORy));
COR = [CORx,CORy];

%%

% adaptive non-maximum suppression of corner candidates

% get the difference in the cornerliness scores for all the corner candidate pixels
CMC = cimg(sub2ind(size(cimg),COR(:,2),COR(:,1)));
difference = bsxfun(@minus,CMC(:)',CMC(:));

% figure out which candidates are bigger than each other, assign boolean
difference(difference <= 0) = 0;
difference(difference > 0) = 1;

% get all the distances between the corner candidates
D = pdist(COR,'euclidean');
D = squareform(D);
DIST = D.*difference;
DIST(DIST == 0) = NaN;
[RAD,IDX] = min(DIST,[],2);
RAD(isnan(RAD)) = -Inf;
[topRAD,topRADind] = sort(RAD,1,'descend');
rmax = topRAD(max_pts);
x = COR(topRADind(1:max_pts),1);
y = COR(topRADind(1:max_pts),2);

