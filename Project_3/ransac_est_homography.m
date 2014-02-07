function [H,inlier_ind] = ransac_est_homography(y1,x1,y2,x2,thresh)

pts1 = [x1,y1];
pts2 = [x2,y2];

inliers = logical(zeros(length(pts1),3000));

for i = 1:3000
    rset = randi(length(pts1),4,1);
    while length(rset) ~= length(unique(rset))
        rset = randi(length(pts1),4,1);
    end

    H = est_homography(x1(rset),y1(rset),x2(rset),y2(rset));
    [X2,Y2] = apply_homography(H,x2,y2);
    inliers(:,i) = (((Y2-y1).^2+(X2-x1).^2) < thresh);
end

[val,ind] = max(sum(inliers,1));
inlier_ind = find(inliers(:,ind));
x1 = x1(inliers(:,ind));
y1 = y1(inliers(:,ind));
x2 = x2(inliers(:,ind));
y2 = y2(inliers(:,ind));
H = est_homography(x1,y1,x2,y2);


