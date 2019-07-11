%Replicator code

Niter = 2*10^5; % [-] number of iterations of the replicator equation
Iavg = 100000; % [-] How many of the last time steps do we save?
dtfact = 0.01; %Max percentage of change per time step

% Saved fitnesses for the last time steps
FitC = zeros(1,Iavg); % [day^-1] 
FitF = zeros(1,Iavg);
FitA = zeros(1,Iavg);
FitJ = zeros(1,Iavg);
FitM = zeros(1,Iavg);
FitB = zeros(1,Iavg);

%Coefficient to prevent extinction of strategies
coeff = 10^-7; % [-]
dC0 = coeff*P.C; % [gC m^-3] Minimum concentration of organisms in a strategy
dF0 = coeff*P.F;
dJ0 = coeff*P.J;
dA0 = coeff*P.A;
dM0 = coeff*P.M;
dB0 = coeff*P.B;


