% Solve the model using value function iteration.
% 
% Description of the model.
% V_u(A,Y) = max{A_u',Y'} {u(C_u) + beta((1-q)E[V_u(A_u',Y')] + qE[V_e(A_u',Y')])}
% V_e(A,Y) = max{A_e',Y'} {u(C_e) + beta((1-q)E[V_u(A_e',Y')] + qE[V_e(A_e',Y')])}
% C_u = A+b-(A_u')/(1+r)
% C_e = A+Y-(A_e')/(1+r)
% Y follows a stochastic process
% Y' = (1-rho)mu + rho*Y + epsilon' with epsilon ~  N(0,sigma^2)
% u(C) = C^{1-\gamma}/(1-gamma)
% A >= 0

% Store the parameters in a structure
Par.beta = 0.98;
Par.gamma = 2;
%Par.alpha = 0.36;
%Par.delta = 0.03;
Par.rho = 0.95;
Par.sigma = 0.1;
Par.mu = 1;
Par.b = 0.4;
Par.p = 0.05;
Par.q = 0.25;
Par.r = 0.01;

%% Solve for the steady state
%Zstar = 0;
%Kstar = ((1/Par.beta - 1 + Par.delta)./(Par.alpha * exp(Zstar))).^(1/(Par.alpha-1));

%% Create a grid for Z
% Note if unemployment didn't exist then in expectation,
% Y= (1-rho)mu + rho*Y + e => Y = mu + e/(1-rho) => E[Y] = mu + 0 = mu
meanZ = Par.mu; %(1-Par.rho)*Par.mu;
Grid.nZ = 7;  % number of points in our grid for Z
numStdZ = 3;  % number of standard deviations to cover with the grid
[Grid.Z, Grid.PZ]  = tauchen(Grid.nZ, meanZ, Par.rho, Par.sigma, numStdZ);
    %tauchen(N,mu,rho,sigma,m)
    
Grid.PZ = Grid.PZ'; % this is a 7 x 7 transition matrix for which the columns sum to 1
% the (i,j) element is the probability of moving from j to i.


%% Create a grid for K
Grid.nK = 20;
Grid.K = linspace(0,6,Grid.nK)';  % this is a 20 x 1 array of evenly spaced points

%% Create a product of the two grids
[ZZ,KK] =meshgrid(Grid.Z,Grid.K);
Grid.KK = KK(:);
Grid.ZZ = ZZ(:);

%% Updating guess of value function
be = zeros(6,1);
bu = zeros(6,1);

%% Putting the iteration in value function iteration
Kpe0 = zeros(size(Grid.KK));
Kpu0 = zeros(size(Grid.KK));
%Kpe = ones(size(Grid.KK));
%Kpu = ones(size(Grid.KK));

MAXIT = 2000;


for it = 1:MAXIT

    [Ve, Kpe] = MaxBellmanE(Par,be,bu,Grid);
    [Vu, Kpu] = MaxBellmanU(Par,be,bu,Grid);
        
    % take the expectation of the value function from the perspective of
    % the previous Z
    EVe = reshape(Ve,Grid.nK,Grid.nZ) * Grid.PZ;
    EVu = reshape(Vu,Grid.nK,Grid.nZ) * Grid.PZ; %Grid.PZ;

    % update our polynomial coefficients
    bu = PolyGetCoef(Grid.KK,Grid.ZZ,EVu(:));
    be = PolyGetCoef(Grid.KK,Grid.ZZ,EVe(:));
    
    
    if it > 1
        % see how much our policy rule has changed
        teste = max(abs(Kpe0 - Kpe));
        testu = max(abs(Kpu0 - Kpu));
        
        Kpe0 = Kpe;
        Kpu0 = Kpu;

        disp(['iteration ' num2str(it) ', teste = ' num2str(teste) ', testu = ' num2str(testu)])
        if teste < 1e-5 && testu < 1e-5
            break
        end
    end
end

%% Results
DK = Grid.K/Kp-1; % Capital grid as percent deviation from steady state

DKp = reshape(Kp,Grid.nK,Grid.nZ)./reshape(Grid.KK,Grid.nK,Grid.nZ) - 1;
  % savings policy rule as a 20 x 7 array expressed as a percent change from current K

plot(DK, DKp);  % plot the policy rule

hold on;        % next plots go on the same figure

plot(DK, zeros(Grid.nK,1), 'k--'); % add a zero line, k-- means black and dashsed

xlabel('K in % deviation from steady state')  % label the axes
ylabel('(K'' - K)/K')