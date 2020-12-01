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
[x_priori, F_x, F_u] = transitionFunction(x, u, b);
Q = [k*abs(u(1,1)) 0; 0 k*abs(u(2,1))];
P_priori = F_x * P * F_x' + F_u * Q * F_u';

[v, H, R] = associateMeasurements(x_priori, P_priori, Z, R, M, g);

y = reshape(v, [], 1);
H = reshape(permute(H, [1,3,2]), [], length(x));
R = blockDiagonal(R);

% update state estimates (pp. 335)
S = H * P_priori * H' + R;
K = P_priori * H' * inv(S);

P_posteriori = (eye(3) - K*H)*P_priori;
x_posteriori = x_priori + K * y;
