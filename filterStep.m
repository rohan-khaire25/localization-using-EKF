function [x_posteriori, P_posteriori] = filterStep(x, P, u, Z, R, M, k, g, b)
% [x_posteriori, P_posteriori] = filterStep(x, P, u, z, R, M, k, g, b)
% returns an a posteriori estimate of the state and its covariance

% additional bells and whistles in case no line was detected, please
% incorporate this at a sensical position in your code
if size(Z,2) == 0
    x_posteriori = x_priori;
    P_posteriori = P_priori;
    return;
end
[f, F_x, F_u] = transitionFunction(x,u, b);
[v, H, R] = associateMeasurements(x, P, Z, R, M, g);
[h,Hx] = measurementFunction(x,M);
%p = [h(1); h(2)] - [Z(1,i); Z(2,i)];
O = Hx*P*Hx' + R %R(:,:,1);
K = P*Hx'*inv(O);
x_posteriori = f+K*v*0.01;
x_posteriori = x_posteriori(:,1);
P_posteriori = P-K*O*K';

