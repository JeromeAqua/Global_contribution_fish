%Replicator code
%global C F J M A B %the global variables are the proportions of each population using each strategy

Niter = 10^6; % [-] number of iterations of the replicator equation
Iavg = Niter; %100000; % [-] How many of the last time steps do we save?
dtfact = 0.01; %Max percentage of change per time step
P = Parameters();

% Saved fitnesses for the last time steps
FitC = zeros(1,Iavg); % [day^-1] 
FitF = zeros(1,Iavg);
FitA = zeros(1,Iavg);
FitJ = zeros(1,Iavg);
FitM = zeros(1,Iavg);
FitB = zeros(1,Iavg);

%Saved distributions for the last time steps
MAday = zeros(P.n,Iavg); % [-]
MAnight = MAday;
MBday = MAday; MBnight = MAday;
MCday = MAday; MCnight = MAday;
MJday = MAday; MJnight = MAday;
MMday = MAday; MMnight = MAday;
MFday = MAday; MFnight = MAday;

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

i = Niter;
notdone = 1;
tic
while notdone
    i = i - 1;

%Denominators for ingestion rates calculations
    NF1 = P.IDF + P.EDF.*(pref('forage','detritus')*repmat(P.D',1,P.n)+pref('forage','copepod')*repmat(Cday',1,P.n)+...
                          pref('forage','benthos')*repmat(P.Benthos',1,P.n)+pref('forage','meso')*repmat(Mday',1,P.n)); % [gC day^-1] Denominator for ingestion function of forage fish during day
    NF0 = P.INF + P.ENF.*(pref('forage','detritus')*repmat(P.D,P.n,1)+pref('forage','copepod')*repmat(Cnight,P.n,1)+...
                          pref('forage','benthos')*repmat(P.Benthos,P.n,1)+pref('forage','meso')*repmat(Mnight,P.n,1)); % [gC day^-1] Denominator for the ingestion function
    NA1 = P.IDA + P.EDA.*(pref('top','forage')*repmat(Fday',1,P.n)+pref('top','tactile')*repmat(Jday',1,P.n)+...
                          pref('top','bathy')*repmat(Bday',1,P.n)+pref('top','meso')*repmat(Mday',1,P.n)); % [gC day^-1] Denominator for ingestion function of top predator during day
    NA0 = P.INA + P.ENA.*(pref('top','forage')*repmat(Fnight,P.n,1)+pref('top','tactile')*repmat(Jnight,P.n,1)+...
                          pref('top','bathy')*repmat(Bnight,P.n,1)+pref('top','meso')*repmat(Mnight,P.n,1)); % [gC day^-1] Denominator for the ingestion function           
    NC1 = P.IDC + P.EDC.*(pref('copepod','phyto')*repmat(P.R',1,P.n)+pref('copepod','detritus')*repmat(P.D',1,P.n)); % [gC day^-1] Denominator for ingestion function of copepods during day
    NC0 = P.INC + P.ENC.*(pref('copepod','phyto')*repmat(P.R,P.n,1)+pref('copepod','detritus')*repmat(P.D,P.n,1)); % [gC day^-1] Denominator for the ingestion function                   
    %J: No denominator because functional response type I
    NM1 = P.IDM + P.EDM.*(pref('meso','detritus')*repmat(P.D',1,P.n)+pref('meso','copepod')*repmat(Cday',1,P.n)); % [gC day^-1] Denominator for ingestion function of mesopelagic during day
    NM0 = P.INM + P.ENM.*(pref('meso','detritus')*repmat(P.D,P.n,1)+pref('meso','copepod')*repmat(Cnight,P.n,1)); % [gC day^-1] Denominator for the ingestion function                   
    NB1 = P.IDB + P.EDB.*(pref('bathy','detritus')*repmat(P.D',1,P.n)+pref('bathy','benthos')*repmat(P.Benthos',1,P.n)+...
                          pref('bathy','copepod')*repmat(Cday',1,P.n)+pref('bathy','meso')*repmat(Mday',1,P.n)); % [gC day^-1] Denominator for ingestion function of top predator during day
    NB0 = P.INB + P.ENB.*(pref('bathy','detritus')*repmat(P.D,P.n,1)+pref('bathy','benthos')*repmat(P.Benthos,P.n,1)+...
                          pref('bathy','copepod')*repmat(Cnight,P.n,1)+pref('bathy','meso')*repmat(Mnight,P.n,1)); % [gC day^-1] Denominator for the ingestion function  


%Ingestion rates
    IFC1 = P.IDF.*P.EDF*pref('forage','copepod').* repmat(Cday' ,1,P.n)./NF1; % [gC day^-1] Ingestion rate of copepods during daytime by forage fish
    IFC0 = P.INF.*P.ENF*pref('forage','copepod').* repmat(Cnight,P.n,1)./NF0;
    IFD1 = P.IDF.*P.EDF*pref('forage','detritus').*repmat(P.D',1,P.n)  ./NF1;
    IFD0 = P.INF.*P.ENF*pref('forage','detritus').*repmat(P.D ,P.n,1)  ./NF0;
    IFb1 = P.IDF.*P.EDF*pref('forage','benthos') .*repmat(P.Benthos',1,P.n)./NF1; % [gC day^-1] small b to show that it is for benthos and not bathypelagic fish
    IFb0 = P.INF.*P.ENF*pref('forage','benthos') .*repmat(P.Benthos, P.n,1)./NF0;
    IFM1 = P.IDF.*P.EDF*pref('forage','meso')    .*repmat(Mday' ,1,P.n)./NF1;
    IFM0 = P.INF.*P.ENF*pref('forage','meso')    .*repmat(Mnight,P.n,1)./NF0;

    IAF1 = P.IDA.*P.EDA*pref('top','forage').* repmat(Fday' ,1,P.n)./NA1; % [gC day^-1] Ingestion rate of forage fish during daytime by top predators
    IAF0 = P.INA.*P.ENA*pref('top','forage').* repmat(Fnight,P.n,1)./NA0;
    IAJ1 = P.IDA.*P.EDA*pref('top','tactile').*repmat(Jday',1,P.n)  ./NA1;
    IAJ0 = P.INA.*P.ENA*pref('top','tactile').*repmat(Jnight ,P.n,1)  ./NA0;
    IAM1 = P.IDA.*P.EDA*pref('top','meso') .*repmat(Mday',1,P.n)./NA1; % [gC day^-1]
    IAM0 = P.INA.*P.ENA*pref('top','meso') .*repmat(Mnight, P.n,1)./NA0;
    IAB1 = P.IDA.*P.EDA*pref('top','bathy')    .*repmat(Bday' ,1,P.n)./NA1;
    IAB0 = P.INA.*P.ENA*pref('top','bathy')    .*repmat(Bnight,P.n,1)./NA0;

    ICR1 = P.IDC.*P.EDC*pref('copepod','phyto') .*repmat(P.R',1,P.n)./NC1; % [gC day^-1] Ingestion rate of phytoplankton during daytime by copepods
    ICR0 = P.INC.*P.ENC*pref('copepod','phyto') .*repmat(P.R, P.n,1)./NC0;
    ICD1 = P.IDC.*P.EDC*pref('copepod','detritus').*repmat(P.D',1,P.n)  ./NC1;
    ICD0 = P.INC.*P.ENC*pref('copepod','detritus').*repmat(P.D ,P.n,1)  ./NC0;

    IJC1 = P.EDJ*pref('tactile','copepod').*repmat(Cday',1,P.n); % [gC day^-1] Ingestion rate of copepods during daytime by tactile predators - note the Type I functional response
    IJC0 = P.ENJ*pref('tactile','copepod').*repmat(Cnight,P.n,1);
    IJM1 = P.EDJ*pref('tactile','meso').*repmat(Mday',1,P.n);
    IJM0 = P.ENJ*pref('tactile','meso').*repmat(Mnight,P.n,1);

    IMC1 = P.IDM.*P.EDM*pref('meso','detritus').*repmat(P.D',1,P.n)  ./NM1; % [gC day^-1] Ingestion rate of copepods during daytime by mesopelagic fish
    IMC0 = P.INM.*P.ENM*pref('meso','detritus').*repmat(P.D ,P.n,1)  ./NM0; 
    IMD1 = P.IDM.*P.EDM*pref('meso','copepod') .*repmat(Cday',1,P.n)./NM1; 
    IMD0 = P.INM.*P.ENM*pref('meso','copepod') .*repmat(Cnight, P.n,1)./NM0;

    IBD1 = P.IDB.*P.EDB*pref('bathy','detritus').* repmat(P.D' ,1,P.n)./NB1; % [gC day^-1] Ingestion rate of forage fish during daytime by top predators
    IBD0 = P.INB.*P.ENB*pref('bathy','detritus').* repmat(P.D,P.n,1)./NB0;
    IBb1 = P.IDB.*P.EDB*pref('bathy','benthos').*repmat(P.Benthos',1,P.n)  ./NB1;
    IBb0 = P.INB.*P.ENB*pref('bathy','benthos').*repmat(P.Benthos ,P.n,1)  ./NB0;
    IBC1 = P.IDB.*P.EDB*pref('bathy','copepod') .*repmat(Cday',1,P.n)./NB1; % [gC day^-1]
    IBC0 = P.INB.*P.ENB*pref('bathy','copepod') .*repmat(Cnight, P.n,1)./NB0;
    IBM1 = P.IDB.*P.EDB*pref('bathy','meso')    .*repmat(Mday' ,1,P.n)./NB1;
    IBM0 = P.INB.*P.ENB*pref('bathy','meso')    .*repmat(Mnight,P.n,1)./NB0;
    
    
%Remove all the NaN of ingestion rates when they do not feed at all
    IFC1(isnan(IFC1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by forage fish
    IFC0(isnan(IFC0)) = 0;
    IFD1(isnan(IFD1)) = 0;
    IFD0(isnan(IFD0)) = 0;
    IFb1(isnan(IFb1)) = 0; % [gC day^-1] small b to show that it is for benthos and not bathypelagic fish
    IFb0(isnan(IFb0)) = 0;
    IFM1(isnan(IFM1)) = 0;
    IFM0(isnan(IFM0)) = 0;
    IAF1(isnan(IAF1)) = 0; % [gC day^-1] Ingestion rate of forage fish during daytime by top predators
    IAF0(isnan(IAF0)) = 0;
    IAJ1(isnan(IAJ1)) = 0;
    IAJ0(isnan(IAJ0)) = 0;
    IAM1(isnan(IAM1)) = 0 ; % [gC day^-1]
    IAM0(isnan(IAM0)) = 0;
    IAB1(isnan(IAB1)) = 0;
    IAB0(isnan(IAB0)) = 0;
    ICR1(isnan(ICR1)) = 0; % [gC day^-1] Ingestion rate of phytoplankton during daytime by copepods
    ICR0(isnan(ICR0)) = 0;
    ICD1(isnan(ICD1)) = 0;
    ICD0(isnan(ICD0)) = 0;
    IJC1(isnan(IJC1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by tactile predators - note the Type I functional response
    IJC0(isnan(IJC0)) = 0;
    IJM1(isnan(IJM1)) = 0;
    IJM0(isnan(IJM0)) = 0;
    IMC1(isnan(IMC1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by mesopelagic fish
    IMC0(isnan(IMC0)) = 0; 
    IMD1(isnan(IMD1)) = 0; 
    IMD0(isnan(IMD0)) = 0;
    IBD1(isnan(IBD1)) = 0; % [gC day^-1] Ingestion rate of forage fish during daytime by top predators
    IBD0(isnan(IBD0)) = 0;
    IBb1(isnan(IBb1)) = 0;
    IBb0(isnan(IBb0)) = 0;
    IBC1(isnan(IBC1)) = 0; % [gC day^-1]
    IBC0(isnan(IBC0)) = 0;
    IBM1(isnan(IBM1)) = 0;
    IBM0(isnan(IBM0)) = 0 ;
    

%Assimilation rates
    IC = P.fC*(P.sigma*(ICR1+ICD1)+(1-P.sigma)*(ICR0+ICD0))/P.wC; % [day^-1] Total assimilation rate per individual per strategy for copepods
    IF = P.fF*(P.sigma*(IFD1+IFb1+IFC1+IFM1)+(1-P.sigma)*(IFD0+IFb0+IFC0+IFM0))/P.wF; % [day^-1] Total assimilation rate per individual per strategy for forage fish
    IA = P.fA*(P.sigma*(IAF1+IAJ1+IAM1+IAB1)+(1-P.sigma)*(IAF0+IAJ0+IAM0+IAB0))/P.wA; % [day^-1] Total assimilation rate per individual per strategy for top predator
    IJ = P.fJ*(P.sigma*(IJC1+IJM1)+(1-P.sigma)*(IJC0+IJM0))/P.wJ; % [day^-1] Total assimilation rate per individual per strategy for tactile predator
    IM = P.fM*(P.sigma*(IMD1+IMC1)+(1-P.sigma)*(IMD0+IMC0))/P.wM; % [day^-1] Total assimilation rate per individual per strategy for mesopelagic fish
    IB = P.fB*(P.sigma*(IBD1+IBb1+IBC1+IBM1)+(1-P.sigma)*(IBD0+IBb0+IBC0+IBM0))/P.wB; % [day^-1] Total assimilation rate per individual per strategy for bathypelagic fish

%Mortality rates due to predation - we calculate the mortality rate imposed by all strategies (predators) at each depth before redistributing it equally among prey
%Copepod
    mCday = IFC1*P.n^2*P.F.*F/P.wF + IJC1*P.n^2*P.J.*J/P.wJ + IMC1*P.n^2*P.M.*M/P.wM + IBC1*P.n^2*P.B.*B/P.wB; % [gC m^-3 day^-1] size n*n How much each strategy eats copepods during daytime
    mCD = (sum(mCday,2)./Cday')'; % [day^-1] size 1*n What is the mortality rate experienced at each depth during day
    MortDa = C.*repmat(mCD',1,P.n); % [day^-1] Mortality rate experienced by the different copepod strategies during daytime

    mCnight = IFC0*P.n^2*P.F.*F/P.wF + IJC0*P.n^2*P.J.*J/P.wJ + IMC0*P.n^2*P.M.*M/P.wM + IBC0*P.n^2*P.B.*B/P.wB; % [gC m^-3 day^-1]
    mCN = sum(mCnight,1)./Cnight; % [day^-1]
    MortNi = C.*repmat(mCN,P.n,1); % [day^-1] Mortality rate experienced by the different bathypelagic fish strategies during nighttime

    MortC = P.sigma*MortDa + (1-P.sigma)*MortNi; % [day^-1] Total mortality rate experienced by the different copepod strategies

%Forage fish
    mFday = IAF1*P.n^2*P.A.*A/P.wA; % [gC m^-3 day^-1] size n*n How much each strategy eats forage fish during daytime
    mFD = (sum(mFday,2)./Fday')'; % [day^-1] size 1*n What is the mortality rate experienced at each depth during day
    MortDa = F.*repmat(mFD',1,P.n); % [day^-1] Mortality rate experienced by the different forage fish strategies during daytime

    mFnight = IAF0*P.n^2*P.A*A/P.wA; % [gC m^-3 day^-1]
    mFN = sum(mFnight,1)./Fnight; % [day^-1]
    MortNi = F.*repmat(mFN,P.n,1); % [day^-1] Mortality rate experienced by the different forage fish strategies during nighttime

    MortF = P.sigma*MortDa + (1-P.sigma)*MortNi; % [day^-1] Total mortality rate experienced by the different forage fish strategies

%Tactile predator
    mJday = IAJ1*P.n^2*P.A.*A/P.wA; % [gC m^-3 day^-1] size n*n How much each strategy eats tactile pred during daytime
    mJD = (sum(mJday,2)./Jday')'; % [day^-1] size 1*n What is the mortality rate experienced at each depth during day
    MortDa = J.*repmat(mJD',1,P.n); % [day^-1] Mortality rate experienced by the different tactile pred strategies during daytime

    mJnight = IAJ0*P.n^2*P.A*A/P.wA; % [gC m^-3 day^-1]
    mJN = sum(mJnight,1)./Jnight; % [day^-1]
    MortNi = J.*repmat(mJN,P.n,1); % [day^-1] Mortality rate experienced by the different tactile pred strategies during nighttime

    MortJ = P.sigma*MortDa + (1-P.sigma)*MortNi; % [day^-1] Total mortality rate experienced by the different tactile predator strategies

%Mesopelagic fish
    mMday = IFM1*P.n^2*P.F.*F/P.wF + IJM1*P.n^2*P.J.*J/P.wJ + IAM1*P.n^2*P.A.*A/P.wA + IBM1*P.n^2*P.B.*B/P.wB; % [gC m^-3 day^-1] size n*n How much each strategy eats mesopelagic during daytime
    mMD = (sum(mMday,2)./Mday')'; % [day^-1] size 1*n What is the mortality rate experienced at each depth during day
    MortDa = M.*repmat(mMD',1,P.n); % [day^-1] Mortality rate experienced by the different mesopelagic fish strategies during daytime

    mMnight = IFM0*P.n^2*P.F.*F/P.wF + IJM0*P.n^2*P.J.*J/P.wJ + IAM0*P.n^2*P.A.*A/P.wA + IBM0*P.n^2*P.B.*B/P.wB; % [gC m^-3 day^-1]
    mMN = sum(mMnight,1)./Mnight; % [day^-1]
    MortNi = M.*repmat(mMN,P.n,1); % [day^-1] Mortality rate experienced by the different bathypelagic fish strategies during nighttime

    MortM = P.sigma*MortDa + (1-P.sigma)*MortNi; % [day^-1] Total mortality rate experienced by the different copepod strategies

%Bathypelagic fish
    mBday = IAB1*P.n^2*P.A.*A/P.wA; % [gC m^-3 day^-1] size n*n How much each strategy eats bathypelagic fish during daytime
    mBD = (sum(mBday,2)./Bday')'; % [day^-1] size 1*n What is the mortality rate experienced at each depth during day
    MortDa = B.*repmat(mBD',1,P.n); % [day^-1] Mortality rate experienced by the different bathypelagic fish strategies during daytime

    mBnight = IAB0*P.n^2*P.A*A/P.wA; % [gC m^-3 day^-1]
    mBN = sum(mBnight,1)./Bnight; % [day^-1]
    MortNi = B.*repmat(mBN,P.n,1); % [day^-1] Mortality rate experienced by the different bathypelagic fish strategies during nighttime

    MortB = P.sigma*MortDa + (1-P.sigma)*MortNi; % [day^-1] Total mortality rate experienced by the different bathypelagic fish strategies


%Fitnesses
    fitA = IA - 0 - P.CA/P.wA - P.metA - 0.01; % [day^-1] Fitness of top predator 0.05 is a background mortality rate
    fitC = IC - MortC - P.CC/P.wC - P.metC - 0.2 ; % [day^-1]
    fitJ = IJ - MortJ - P.CJ/P.wJ - P.metJ -0.1 ; % [day^-1]
    fitF = IF - MortF - P.CF/P.wF - P.metF -0.05; % [day^-1]
    fitM = IM - MortM - P.CM/P.wM - P.metM -0.05; % [day^-1]
    fitB = IB - MortB - P.CB/P.wB - P.metB -0.03; % [day^-1]

%%% NOW IS THE REPLICATOR PART
    FAmax = max(max(fitA)); FAmin = min(min(fitA));
    FBmax = max(max(fitB)); FBmin = min(min(fitB));
    FCmax = max(max(fitC)); FCmin = min(min(fitC));
    FJmax = max(max(fitJ)); FJmin = min(min(fitJ));
    FMmax = max(max(fitM)); FMmin = min(min(fitM));
    FFmax = max(max(fitF)); FFmin = min(min(fitF));


%the proportionality factor is dynamic so that maximum increase is at most 2% per time step
     factA = dtfact./max([FAmax, -FAmin]); % [day]
     factB = dtfact./max([FBmax, -FBmin]);
     factC = dtfact./max([FCmax, -FCmin]);
     factJ = dtfact./max([FJmax, -FJmin]);
     factM = dtfact./max([FMmax, -FMmin]);
     factF = dtfact./max([FFmax, -FFmin]);

%increment, the core of the replicator equation
    A = A.*(1 + factA*fitA); % [-] Proportion of all strategies, before renormalization
    B = B.*(1 + factB*fitB);
    C = C.*(1 + factC*fitC);
    J = J.*(1 + factJ*fitJ);
    M = M.*(1 + factM*fitM);
    F = F.*(1 + factF*fitF);
    
%no extinction, so that every strategy can still emerge
    A(A<dA0) = dA0; % [-]
    B(B<dB0) = dB0;
    C(C<dC0) = dC0;
    J(J<dJ0) = dJ0;
    M(M<dM0) = dM0;
    F(F<dF0) = dF0;
    
%Renormalization
    A = A/sum(sum(A)); % [-] Good matrices of strategies after each replicator time step
    B = B/sum(sum(B));
    C = C/sum(sum(C));
    M = M/sum(sum(M));
    J = J/sum(sum(J));
    F = F/sum(sum(F));

%Calculation of the day and night concentrations
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
     
     
%Save historic of convergence
    if i<Iavg
        MAday(:,Iavg-i) = Aday;
        MBday(:,Iavg-i) = Bday;
        MCday(:,Iavg-i) = Cday;
        MMday(:,Iavg-i) = Mday;
        MFday(:,Iavg-i) = Fday;
        MJday(:,Iavg-i) = Jday;
        MAnight(:,Iavg-i) = Anight;
        MBnight(:,Iavg-i) = Bnight;
        MCnight(:,Iavg-i) = Cnight;
        MMnight(:,Iavg-i) = Mnight;
        MFnight(:,Iavg-i) = Fnight;
        MJnight(:,Iavg-i) = Jnight;
    end
        
     
     notdone = (i>0);
end
toc