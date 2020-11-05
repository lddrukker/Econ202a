function y  = f(Par, K,Z )
% y  = f( K,Z )
%   Production function gross of undepreciated capital

y = (1/1+Par.r)*(K+Z - 1e-3); %exp(Z) .* K.^Par.alpha + (1-Par.delta)*K;
end