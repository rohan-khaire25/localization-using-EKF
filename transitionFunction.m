function [f, F_x, F_u] = transitionFunction(x,u, b)
% [f, F_x, F_u] = transitionFunction(x,u,b) predicts the state x at time t given
% the state at time t-1 and the input u at time t. F_x denotes the Jacobian
% of the state transition function with respect to the state evaluated at
% the state and input provided. F_u denotes the Jacobian of the state
% transition function with respect to the input evaluated at the state and
% input provided.
% State and input are defined according to "Introduction to Autonomous Mobile Robots", pp. 337
%{
syms g y theta del_sl del_sr o

J = [((del_sl+del_sr)/2)*(cos(theta+((del_sr-del_sl)/(2*o)))); ...
                    ((del_sl+del_sr)/2)*(sin(theta+((del_sr-del_sl)/(2*o)))); ...
                        (del_sr-del_sl)/o];
temp1 = jacobian(J,[g y theta]);
temp2 = jacobian(J,[del_sl del_sr]);

g = x(1,1);
y = x(2,1);
theta = x(3,1);
del_sl = u(1,1);
del_sr = u(2,1);
o = b;
temp3 = subs(J);
F_x = subs(temp1);
F_u = subs(temp2);
f = [x(1,1); x(2,1); x(3,1)] + temp3;
%}
del_s = (u(2)+u(1))/2;
del_t = (u(2)-u(1))/b;
F_x = [0 0 -del_s*sin(x(3)+del_t/2); 0 0 del_s*cos(x(3)+del_t/2);0 0 0];
        
F_u = [0.5*cos(x(3)+del_t/2)+(del_s/(2*b))*sin(x(3)+del_t/2), 0.5*cos(x(3)+del_t/2)-(del_s/(2*b))*sin(x(3)+del_t/2); ...
       0.5*sin(x(3)+del_t/2)-(del_s/(2*b))*cos(x(3)+del_t/2), 0.5*sin(x(3)+del_t/2)+(del_s/(2*b))*cos(x(3)+del_t/2); ...
       -1/b, 1/b];
f = [x(1); x(2); x(3)] + [del_s*cos(x(3)+del_t/2);del_s*sin(x(3)+del_t/2);del_t];

