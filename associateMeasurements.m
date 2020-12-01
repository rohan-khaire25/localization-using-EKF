function [v, H, R] = associateMeasurements(x, P, Z, R, M, g)
% [v, H, R] = associateMeasurements(x, P, Z, R, M, g) returns a set of
% innovation vectors and associated jacobians and measurement covariances
% by matching line features by Mahalanobis distance.
v = [];
H = zeros(2,3,5);
temp = [0.1000 0;0 0.1000];
md_val = [];
for index = 1:size(Z,2)
    for j = 1:size(M,2)
        [h, Hx] = measurementFunction(x,M(:,j));
        p = [Z(1,index) ; Z(2,index)] - [h(1);h(2)];
        O = Hx*P*Hx' + R(:,:,1);
        md = double(p'*inv(O)*p);
        md_val(index,j) = md;
        %mm = min(md_val(:,j));
        if md <= g^2
            if index >= 2 && md_val(index,j) < min(md_val(1:index-1,j))
                 v(:,j) = p;
            else
                 v(:,j) = p;
            end
            
            H(:,:,j) = Hx;
        end
    end
    
end
R = repmat(temp,1,1,size(H,3));
     
end

