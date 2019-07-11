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


Cday = P.n*P.C*sum(C,2)'; % [gC m^-3] Average concentration in each layer during day for copepod
Cnight = P.n*P.C*sum(C,1); % [gC m^-3] Average concentration in each layer during night
Fday = P.n*P.F*sum(F,2)'; % [gC m^-3] Average concentration in each layer during day for forage fish
Fnight = P.n*P.F*sum(F,1); % [gC m^-3] Average concentration in each layer during night
Aday = P.n*P.A*sum(A,2)'; % [gC m^-3] Average concentration in each layer during day for top predator
Anight = P.n*P.A*sum(A,1); % [gC m^-3] Average concentration in each layer during night
Bday = P.n*P.B*sum(B,2)'; % [gC m^-3] Average concentration in each layer during day for bathypelagic fish
Bnight = P.n*P.B*sum(B,1); % [gC m^-3] Average concentration in each layer during night
Mday = P.n*P.M*sum(M,2)'; % [gC m^-3] Average concentration in each layer during day for mesopelagic fish
Mnight = P.n*P.M*sum(M,1); % [gC m^-3] Average concentration in each layer during night
Jday = P.n*P.J*sum(J,2)'; % [gC m^-3] Average concentration in each layer during day for tactile predator
Jnight = P.n*P.J*sum(J,1); % [gC m^-3] Average concentration in each layer during night


%Denominators for ingestion rates calculations
NF1 = P.IDF + P.EDF.*(pref('forage','detritus')*repmat(P.D',1,n)+pref('forage','copepod')*repmat(Cday',1,n)+...
                      pref('forage','benthos')*repmat(P.Benthos',1,n)+pref('forage','meso')*repmat(Mday',1,n)); % [gC day^-1] Denominator for ingestion function of forage fisg during day
NF0 = P.INF + P.ENF.*(pref('forage','detritus')*repmat(P.D,n,1)+pref('forage','copepod')*repmat(Cday,n,1)+...
                      pref('forage','benthos')*repmat(P.Benthos,n,1)+pref('forage','meso')*repmat(Mday,n,1)); % [gC day^-1] Denominator for the ingestion function
NA1 =
NA0 = 
                  
                  
                  
%Ingestion functions
IFC1 = P.IDF.*P.EDF*pref('forage','copepod')*repmat(Cday',1,n)./NF1; % [gC day^-1] Ingestion rate of copepods during daytime by forage fish
IFC0 = 
IFD1 = 
IFD0 = 
IFB1 = 
IFB0 =
IFM1 = 
IFM0 = ...;
