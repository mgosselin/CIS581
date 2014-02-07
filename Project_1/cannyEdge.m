function E = cannyEdge(I)
% This function should take in an image I (dim = 3) and output an edge map E (dim = 2).

I = im2double(I);

% use upper and lower bounds for hysteresis
upper = 0.10;
lower = 0.05;

%%

% first, dilate the image, surrounding it by zeros

perim = 3;
[row,col,nb] = size(I);

rowstart = 1+perim;
rowend = row+perim;
colstart = 1+perim;
colend = col+perim;

% create some black space around the picture
J = zeros(row+2*perim,col+2*perim,nb);
J(rowstart:rowend,colstart:colend,:) = I;
J = rgb2gray(J);

%%

% Perform mirroring to pad the region outside the image

% do all 4 edges first
for i = 1:perim
    % left edge
    J(rowstart:rowend,i) = J(rowstart:rowend,1+2*perim-i);
    % right edge
    J(rowstart:rowend,1+col+2*perim-i) = J(rowstart:rowend,i+col);
    % top edge
    J(i,colstart:colstart+col) = J(1+2*perim-i,colstart:colstart+col);
    % bottom edge
    J(1+row+2*perim-i,colstart:colstart+col) = J(i+row,colstart:colstart+col);
end

% do all the corners second
for i = 1:perim
    for j = 1:perim
        % top left
        J(i,j) = J(1+2*perim-i,1+2*perim-j);
        % bottom left
        J(1+row+2*perim-i,j) = J(row+i,1+2*perim-j);
        % top right
        J(rowstart-i,colend+j) = J(perim+i,colstart+col-j);
        % bottom right
        J(1+row+2*perim-i,1+col+perim*2-j) = J(i+row,j+col);
    end
end

%%

% make a smoothing filter and convolve it with the derivative...then use this to get the gradient vectors in the x and y coordinates

% make a 2-D guassian smoothing filter
G = fspecial('gaussian',[3 3],0.60);

% define the discrete-domain first-derivatives in x and y coordinates
dx = [1, -1];
dy = [1; -1];

% convolve the gaussian filter and the derivatives with the image
Ix = conv2(conv2(J(rowstart:rowend,colstart:colend),G,'same'),dx,'same');
Iy = conv2(conv2(J(rowstart:rowend,colstart:colend),G,'same'),dy,'same');

% make the gradient vectors into polar equivalents
[theta, rad] = cart2pol(Ix,Iy);

% imshow(rad);
% hold on
% quiver(1:col,1:row,Ix,Iy);

% E = rad;

%%

% use a variable to store the binary cadidates for local maxima for edge-linking in a later step 
candidate = logical(zeros(row,col));
rad = rad/max(rad(:));
% rad(rad < lower) = 0;

for i = 2:row-1
    for j = 2:col-1
        % use the angle of the gradient to get some pixels to interpolate
        run = cos(theta(i,j));
        rise = sin(theta(i,j));
        % figure out which quadrant the gradient vector is pointing
        if 0 < rise && 0 < run
            X = (j:j+1);
            Y = (i-1:i)';
        elseif 0 < rise && run <= 0
            X = (j-1:j);
            Y = (i-1:i)';
        elseif rise <= 0 && run <= 0
            X = (j-1:j);
            Y = (i:i+1)';
        elseif rise <= 0 && 0 < run
            X = (j:j+1);
            Y = (i:i+1)';
        end
        % get the values of the gradient magnitude at the four pixels in vectors X,Y
        V = rad(Y,X);
        % interpolate to get the pixel ahead (r) and the pixel behind (p) the current pixel            
        r = interp2(X,Y,V,j+cos(theta(i,j)),i-sin(theta(i,j)));
        % use the angle of the gradient to get some pixels to interpolate
        run = cos(theta(i,j)-pi);
        rise = sin(theta(i,j)-pi);
        % figure out which quadrant the gradient vector is pointing
        if 0 < rise && 0 < run
            X = (j:j+1);
            Y = (i-1:i)';
        elseif 0 < rise && run <= 0
            X = (j-1:j);
            Y = (i-1:i)';
        elseif rise <= 0 && run <= 0
            X = (j-1:j);
            Y = (i:i+1)';
        elseif rise <= 0 && 0 < run
            X = (j:j+1);
            Y = (i:i+1)';
        end
        % get the values of the gradient magnitude at the four pixels in vectors X,Y
        V = rad(Y,X);
        % interpolate to get the pixel ahead (r) and the pixel behind (p) the current pixel
        p = interp2(X,Y,V,j+cos(theta(i,j)-pi),i-sin(theta(i,j)-pi));      
        % check if the next pixel's gradient is bigger
        if (r < rad(i,j)) && (p < rad(i,j))
            % if current pixel is biggest, mark it as a local max
            candidate(i,j) = 1;
        end
    end
end

% E = candidate;

%%

% now do linking of the binary edge candidates

% use a variable to keep track of which pixels have been visited
visited = logical(zeros(row,col));

% now look for places where edges might start (where gradient magnitude is large)
starter = logical(zeros(row,col)) + (rad > upper);

% use a variable to store the binary edgemap for output at the end
output = logical(zeros(row,col));

% use a variable to switch between 'walking' mode and 'searching' mode
walk = 0;

for i = 2:row-1
    for j = 2:col-1
        % check for unvisited starter pixels which are also edge candidates
        if visited(i,j) == 0 && starter(i,j) == 1 && walk == 0 && candidate(i,j) == 1
            % mark starter pixel as 'visited'
            visited(i,j) = 1;
            % mark the starter pixel as a member of the edge
            output(i,j) = 1;
            % enter 'walking' mode
            walk = 1;
            walknum = 0;
            % look along the tangent to gradient vector, and get nearest indices
            next_i = round(i-sin(theta(i,j)+pi/2));
            next_i = min(row,next_i);
            next_i = max(1,next_i);
            next_j = round(j+cos(theta(i,j)+pi/2));
            next_j = min(col,next_j);
            next_j = max(1,next_j);
            % walk to the next edge candidate pixel and check if it's unvisited and if it is still above the lower threshold
            while walk == 1 && rad(next_i,next_j) > lower && visited(next_i,next_j) == 0 && candidate(next_i,next_j) == 1
                % mark the new pixel as visited and as an edge pixel
                output(i,j) = 1;
                visited(i,j) = 1;
                % turn the 'next' pixel into the 'current' pixel
                i = next_i; 
                j = next_j;
                % get the indices for the 'next' pixel
                next_i = round(i-sin(theta(i,j)+pi/2));
                next_i = min(row,next_i);
                next_i = max(1,next_i);
                next_j = round(j+cos(theta(i,j)+pi/2));
                next_j = min(col,next_j);
                next_j = max(1,next_j);
                % update the lower threshold
                lower = 0.05./walknum + 0.005;
                walknum = walknum + 1;
            end
            % exit 'walking' mode
            walk = 0;
            walknum = 0;
        end
    end
end

%%

E = output;

