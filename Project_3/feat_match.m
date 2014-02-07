function [m] = feat_match(p1,p2)
%%
thresh = 0.3;
[M,N] = size(p1);

SSD = zeros(N,N);

for i = 1:N
    for j = 1:N        
        SSD(i,j) = sum((p1(:,i)-p2(:,j)).^2);        
    end
end

[DIST,IND] = sort(SSD,2,'ascend');
RATIO = DIST(:,1)./DIST(:,2);
m = IND(:,1);
m(RATIO >= thresh) = -1;
