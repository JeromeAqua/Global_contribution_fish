%Replicator code
global C F J M A B %the global variables are the proportions of each population using each strategy

Niter = 2*10^5; % [-] number of iterations of the replicator equation
Iavg = 100000; % [-] How many of the last time steps do we save?
dtfact = 0.01; %Max percentage of change per time step
P = Parameters();

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

%Initialization of the different strategy matrices
C = ones(P.n); C = C/sum(sum(C)); % [-] equal distribution at first - we can also use rand(n) to have random distribution
F = ones(P.n); F = F/sum(sum(F)); 
M = ones(P.n); M = M/sum(sum(M));
A = ones(P.n); A = A/sum(sum(A));
B = ones(P.n); B = B/sum(sum(B));
J = ones(P.n); J = J/sum(sum(J));






















function OUT = Clearance(predator,day,Q) %Clearance rate for each strategy (matrix) as a function of day - Q here is our Parameter file
OUT = zeros(n);

if strcmp(predator,'forage')
    if day==1
        Vis = P.RF*sqrt(P.LD./(P.KF+Q.LD))'; % [m] Depth-dependent visual range of forage fish during daytime
        OUT = P.gamma*pi*P.uF*P.MSDF.*repmat(Vis,1,P.n).^2; % [m^3 day^-1]
    elseif day==0
        Vis = P.RF*sqrt(P.LN./(P.KF+P.LN)); % [m] Depth-dependent visual range of forage fish during daytime
        OUT = P.gamma*pi*P.uF*P.MSDF.*repmat(Vis,P.n,1).^2; % [m^3 day^-1]
    end
    
    




end


