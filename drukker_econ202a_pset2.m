% Store the parameters in a structure
Par.beta = 0.99; %0.98;
Par.gamma = 2;
Par.alpha = 0.36;
Par.delta = 0.03;
Par.rho = 0.95;
Par.sigma = 0.007; %0.1
%Par.mu = 1;
%Par.b = 0.4;
%Par.p = 0.05;
%Par.q = 0.25;
%Par.r = 0.01;

%% Solve for the steady state
Zstar = 0;
Kstar = ((1/Par.beta - 1 + Par.delta)./(Par.alpha * exp(Zstar))).^(1/(Par.alpha-1));

%% Create a grid for Z
meanZ = 0;
Grid.nZ = 7;  % number of points in our grid for Z
numStdZ = 2;  % number of standard deviations to cover with the grid
[Grid.Z, Grid.PZ]  = tauchen(Grid.nZ, meanZ, Par.rho, Par.sigma, numStdZ);

Grid.PZ = Grid.PZ'; % this is a 7 x 7 transition matrix for which the columns sum to 1
% the (i,j) element is the probability of moving from j to i.

%% Create a grid for K
Grid.nK = 20;
Grid.K = linspace(0.75*Kstar, 1.25*Kstar,Grid.nK)';  % this is a 20 x 1 array of evenly spaced points

%% Create a product of the two grids
[ZZ,KK] =meshgrid(Grid.Z,Grid.K);
Grid.KK = KK(:);
Grid.ZZ = ZZ(:);



%% Updating guess of value function
%b = zeros(6,1);

%EV = reshape(V,Grid.nK,Grid.nZ) * Grid.PZ;

%b1 = PolyGetCoef(Grid.KK,Grid.ZZ,EV(:));

%% Putting the iteration in value function iteration
Kp0 = zeros(size(Grid.KK));
MAXIT = 2000;
for it = 1:MAXIT

    [V, Kp] = MaxBellman(Par,b,Grid);


    % take the expectation of the value function from the perspective of
    % the previous Z
    EV = reshape(V,Grid.nK,Grid.nZ) * Grid.PZ;

    % update our polynomial coefficients
    b = PolyGetCoef(Grid.KK,Grid.ZZ,EV(:));

    % see how much our policy rule has changed
    test = max(abs(Kp0 - Kp));
    Kp0 = Kp;

    disp(['iteration ' num2str(it) ', test = ' num2str(test)])
    if test < 1e-5
        break
    end
end

%% Results
DK = Grid.K/Kstar-1; % Capital grid as percent deviation from steady state

DKp = reshape(Kp,Grid.nK,Grid.nZ)./reshape(Grid.KK,Grid.nK,Grid.nZ) - 1;
  % savings policy rule as a 20 x 7 array expressed as a percent change from current K

plot(DK, DKp);  % plot the policy rule

hold on;        % next plots go on the same figure

plot(DK, zeros(Grid.nK,1), 'k--'); % add a zero line, k-- means black and dashsed

xlabel('K in % deviation from steady state')  % label the axes
ylabel('(K'' - K)/K')

%% Functions
function y  = f(Par, K,Z )
% y  = f( K,Z )
%   Production function gross of undepreciated capital

y =  exp(Z) .* K.^Par.alpha + (1-Par.delta)*K;
end

function V  = Bellman( Par, b, K, Z, Kp )
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

C = f(Par,K,Z) - Kp;
u = C.^(1-Par.gamma) / (1-Par.gamma);
V = u + Par.beta * PolyBasis(Kp,Z) * b;



end

function [V, Kp] = MaxBellman(Par,b,Grid)
% [V, Kp] = MaxBellman(Par,b,Grid)
%   Maximizes the RHS of the Bellman equation using golden section search
%
% Inputs
% Par       Parameter structure
% b     6 x 1 coefficients in polynomial for E[ V(K',Z') | Z ]
% Grid      Grid structure


p = (sqrt(5)-1)/2;

A = Grid.K(1) * ones(size(Grid.KK));
D = min(f(Par,Grid.KK,Grid.ZZ) - 1e-3, Grid.K(end)); % -1e-3 so we always have positve consumption.

B = p*A+(1-p)*D;
C = (1-p)*A + p * D;


fB = Bellman(Par,b,Grid.KK,Grid.ZZ,B);
fC = Bellman(Par,b,Grid.KK,Grid.ZZ,C);


MAXIT = 1000;
for it_inner = 1:MAXIT

    if all(D-A < 1e-6)
        break
    end

    I = fB > fC;

    D(I) = C(I);
    C(I) = B(I);
    fC(I) = fB(I);
    B(I) = p*C(I) + (1-p)*A(I);
    fB(I) = Bellman(Par,b,Grid.KK(I),Grid.ZZ(I),B(I));

    A(~I) = B(~I);
    B(~I) = C(~I);
    fB(~I) = fC(~I);
    C(~I) = p*B(~I) + (1-p)*D(~I);
    fC(~I) = Bellman(Par,b,Grid.KK(~I),Grid.ZZ(~I),C(~I));

end

% At this stage, A, B, C, and D are all within a small epsilon of one
% another.  We will use the average of B and C as the optimal level of
% savings.
Kp = (B+C)/2;

% evaluate the Bellman equation at the optimal policy to find the new
% value function.
V = Bellman(Par,b,Grid.KK,Grid.ZZ,Kp);

end



