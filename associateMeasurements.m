function [v, H, R] = associateMeasurements(x, P, Z, R, M, g)
% [v, H, R] = associateMeasurements(x, P, Z, R, M, g) returns a set of
% innovation vectors and associated jacobians and measurement covariances
% by matching line features by Mahalanobis distance.
v = zeros(2,5);
H = zeros(2,3,5);
for i = 1:size(Z,2)
  for j = 1:size(M,2)
    [h,Hx] = measurementFunction(x,M(:,j));
    p = [h(1); h(2)] - [Z(1,i); Z(2,i)];
    O = Hx*P*Hx' + R; %R(:,:,j)
    dist = p'*inv(O)*p;
    if dist<g^2
      v(:,j) = p;
      H(:,:,i) = Hx;
    end
  end
end
