function V  = BellmanU( Par, be, bu, K, Z, Kp )
% V  = Bellman( Par, b, K, Z, Kp )
%   Evaluate the RHS of the Bellman equation
%
% Inputs
% Par   Parameter structure
% b     6 x 1 coefficients in polynomial for E[ V(K',Z') | Z ]
% K     n x 1 array of current capital
% Z     n x 1 array of current TFP
% Kp    n x 1 array of savings
%
% Output
% V     n x 1 array of value function
%

C = Par.b + K - Kp/(1+Par.r);
u = C.^(1-Par.gamma) / (1-Par.gamma);
V = u + Par.beta * (PolyBasis(Kp,Par.b)*(Par.q)*be + PolyBasis(Kp,Par.b)*(1-Par.q)*bu); %* b;
%V = u + Par.beta * PolyBasis(Kpu,Z) * b;

end