function img_mosaic = mymosaic(img_input)

max_pts = 1500; % pick how many corners to compute
thresh = 0.5; % pick the threshold
images = length(img_input); % find out how many images there are
H = cell(images-1,1); % pre allocate memory for the homographies
% calculate the homographies
for i = 1:images - 1
    % pull the images out of the input array
    img1 = img_input{i};
    img2 = img_input{i+1};
    % convert to gray
    im1 = rgb2gray(img1);
    im2 = rgb2gray(img2);
    % get the maxima and the feature descriptors of image 1
    cimg1 = cornermetric(im1,'Harris');
    [y1,x1,rmax1] = anms(cimg1,max_pts);
    p1 = feat_desc(im1,y1,x1);
    % get the maxima and the feature descriptors of image 2
    cimg2 = cornermetric(im2,'Harris');
    [y2,x2,rmax2] = anms(cimg2,max_pts);
    p2 = feat_desc(im2,y2,x2);
    % get the feature match
    [m] = feat_match(p1,p2);
    %get the matched points
    pts1 = [x1(m~=-1),y1(m~=-1)];
    pts2 = [x2(m(m~=-1)),y2(m(m~=-1))];
    % calculate the homography
    [HH,inlier_ind] = ransac_est_homography(pts1(:,2),pts1(:,1),pts2(:,2),pts2(:,1),thresh);
    % save the homography in a new cell array for using later
    H{i} = HH;
end

%%
% assign the first image to the mosaic initially
img_mosaic = img_input{1};
shift = 0;
% make a loop to transform each image one by one and copy it to the output
for j = 1:images - 1
    % assign image '1' and image '2' ...not necessarily the 1st and 2nd
    img1 = img_mosaic;
    img2 = img_input{j+1};
    % convert to gray
    im1 = rgb2gray(img1);
    im2 = rgb2gray(img2);
    
    % compute the transformation matrix from image 'j' to the base frame
    HH = 1;
    for h = 1:j
        HH = H{h}*HH;
    end
    
    % get the corners of image 2 in the frame of image 1
    [row,col] = size(im2);
    corners22 = [1,row;col,row;col,1;1,1];
    [corners21(:,1),corners21(:,2)] = apply_homography(HH,corners22(:,1),corners22(:,2));

    % get the corners of image 1 in the frame of image 1
    [row,col] = size(im1);
    corners11 = [1,row;col,row;col,1;1,1];
    
    % compare corners of image 1 & 2 to get the dimensions to expand the bounding box
    % x-dim
    xmin = round(min([corners11(1,1),corners11(4,1),corners21(1,1),corners21(4,1)]));
    xmax = round(max([corners11(2,1),corners11(3,1),corners21(2,1),corners21(3,1)]));
    xdim = xmax-xmin;

    % y-dim
    ymin = round(min([corners11(3,2),corners11(4,2),corners21(3,2),corners21(4,2)]));
    ymax = round(max([corners11(1,2),corners11(2,2),corners21(1,2),corners21(2,2)]));
    ydim = ymax-ymin;
    
    % expand in x
    img_mosaic = [img_mosaic,uint8(zeros(row,(xdim - col),3))];
    
    % expand in y
    img_mosaic = [img_mosaic;uint8(zeros((ydim-row)/2,xdim,3))];
    img_mosaic = [uint8(zeros((ydim-row)/2,xdim,3));img_mosaic];
    
    % compute the transformation matrix to shift im2 down the right amount
    shift = shift + (ydim-row)/2;
    H1 = [1 0 0;0 1 shift;0 0 1];
    
    % get the indices of the box where the second image will go
    xmin = round(min([corners21(1,1),corners21(4,1)]));
    xmax = round(max([corners21(2,1),corners21(3,1)]));
    ymin = round(min([corners21(3,2),corners21(4,2)]));
    ymax = round(max([corners21(1,2),corners21(2,2)]));
%     xdim = xmax-xmin;
%     ydim = ymax-ymin;
    [Y,X] = ndgrid([ymin:ymax],[xmin:xmax]);
    Y = Y(:);
    X = X(:);
    N = length(Y);
    
    % convert these indices into the frame of image 2
    [X,Y] = apply_homography(inv(HH),X,Y);

    % throw away indices outside the pixel-space of image 2
    [row,col] = size(im2);
    X(X < 1 | X > col) = NaN;
    Y(Y < 1 | Y > row) = NaN;
    for k = 1:N
        if isnan(X(k)) | isnan(Y(k))
            X(k) = NaN;
            Y(k) = NaN;
        end
    end
    X = X(~isnan(X));
    Y = Y(~isnan(Y));
    
    % copy image 2 into the output mosaic using the transformed indices
    [X2,Y2] = apply_homography(H1*HH,X,Y);
    X = round(X);
    Y = round(Y);
    
    N = length(Y2);
    for l = 1:N
        x = round(X2(l));
        y = round(Y2(l));
        x = min(max(x,1),size(img_mosaic,2));
        y = min(max(y,1),size(img_mosaic,1));
        img_mosaic(y,x,:) = img2(Y(l),X(l),:);
    end
end
% get the corners of image 3 in the frame of image 2


