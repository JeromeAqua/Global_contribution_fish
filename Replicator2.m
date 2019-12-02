%Replicator code
%global C F J M A B %the global variables are the proportions of each population using each strategy
Niter = 3*10^5; % [-] number of iterations of the replicator equation
Iavg = Niter; %100000; % [-] How many of the last time steps do we save?
dtfact = 0.05; %Max percentage of change per time step
sig = 0.2; % Parameter used for the Gaussian filter default = 0.3
sigf = 0.7; % Parameter used for the Gaussian filter on the fitness matrices
reinit = 1; %Do we start from the last simulation or do we initialize strategy matrices?

minimort = 0.01; % [day^-1] Background mortality
minimortC = 0.1; % [day^-1] Background mortality for small copepods

for k=4
    P = Parameters3(k);

%Coefficient to prevent extinction of strategies - and bugs in the OMZ
coeff = 10^-12;%10^-6; % [-] 10^-5
dC0 = coeff*P.C; % [gC m^-3] Minimum concentration of organisms in a strategy
dF0 = coeff*P.F;
dJ0 = coeff*P.J;
dA0 = coeff*P.A;
dM0 = coeff*P.M;
dB0 = coeff*P.B;
dP0 = coeff*P.P;

%Initialization of the different strategy matrices
if reinit==1
    C = ones(P.n).*P.MaskC + dC0; C = C/sum(sum(C)); % [-] random distributions at first - we can also use ones(n) to have equal distribution
    PC = ones(P.n).*P.MaskP + dP0; PC = PC/sum(sum(PC));
    F = ones(P.n).*P.MaskF + dF0; F = F/sum(sum(F)); 
    M = ones(P.n).*P.MaskM + dM0; M = M/sum(sum(M));
    A = ones(P.n).*P.MaskA + dA0; A = A/sum(sum(A));
    B = ones(P.n).*P.MaskB + dB0; B = B/sum(sum(B));
    J = ones(P.n).*P.MaskJ + dJ0; J = J/sum(sum(J));
    
    A(P.MaskA==0) = 0;
    PC(P.MaskP==0) = 0;
    C(P.MaskC==0) = 0;
    F(P.MaskF==0) = 0;
    M(P.MaskM==0) = 0;
    B(P.MaskB==0) = 0;
    J(P.MaskJ==0) = 0;
end

% Saved fitnesses for the last time steps
FitC = zeros(1,Iavg); % [day^-1] 
FitP = zeros(1,Iavg);
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
MPday = MAday; MPnight = MAday;


Cday = P.n*P.C*sum(C,2)'; % [gC m^-3] Average concentration in each layer during day for copepod
Cnight = P.n*P.C*sum(C,1); % [gC m^-3] Average concentration in each layer during night
Pday = P.n*P.P*sum(PC,2)'; % [gC m^-3] Average concentration in each layer during day for copepod
Pnight = P.n*P.P*sum(PC,1); % [gC m^-3] Average concentration in each layer during night
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

mtry = P.sigma*repmat(P.LD'/P.LD(1),1,P.n)+(1-P.sigma)*repmat(P.LN/P.LN(1),P.n,1);
tic
while notdone
    i = i - 1;

%Denominators for ingestion rates calculations
    NF1 = P.IDF + P.EDFd.*repmat(pref('forage','detritus'),1,P.n).*repmat(P.D',1,P.n)     +P.EDFC.*pref('forage','copepod').*repmat(Cday',1,P.n)+ P.EDFP.*pref('forage','predcop').*repmat(Pday',1,P.n)+...
                  P.EDFb.*pref('forage','benthos').*repmat(P.Benthos',1,P.n)+P.EDFM.*pref('forage','meso').*repmat(Mday',1,P.n); % [gC day^-1] Denominator for ingestion function of forage fish during day
    NF0 = P.INF + P.ENFd.*repmat(pref('forage','detritus')',P.n,1).*repmat(P.D,P.n,1)      +P.ENFC.*pref('forage','copepod').*repmat(Cnight,P.n,1)+P.ENFP.*pref('forage','predcop').*repmat(Pnight,P.n,1)+...
                  P.ENFb.* pref('forage','benthos').*repmat(P.Benthos,P.n,1)+P.ENFM.*pref('forage','meso').*repmat(Mnight,P.n,1); % [gC day^-1] Denominator for the ingestion function
    NA1 = P.IDA + P.EDAF.*pref('top','forage').*repmat(Fday',1,P.n)+P.EDAJ.*pref('top','tactile').*repmat(Jday',1,P.n)+...
                  P.EDAB.*pref('top','bathy').*repmat(Bday',1,P.n) +P.EDAM.*pref('top','meso').*repmat(Mday',1,P.n); % [gC day^-1] Denominator for ingestion function of top predator during day
    NA0 = P.INA + P.ENAF.*pref('top','forage').*repmat(Fnight,P.n,1)+P.ENAJ.*pref('top','tactile').*repmat(Jnight,P.n,1)+...
                  P.ENAB.*pref('top','bathy').*repmat(Bnight,P.n,1)+P.ENAM.*pref('top','meso').*repmat(Mnight,P.n,1); % [gC day^-1] Denominator for the ingestion function           
    NC1 = P.IDC + P.EDCp.*pref('copepod','phyto').*repmat(P.R',1,P.n)+P.EDCd.*repmat(pref('copepod','detritus'),1,P.n).*repmat(P.D',1,P.n); % [gC day^-1] Denominator for ingestion function of copepods during day
    NC0 = P.INC + P.ENCp.*pref('copepod','phyto').*repmat(P.R,P.n,1)+P.ENCd.*repmat(pref('copepod','detritus')',P.n,1).*repmat(P.D,P.n,1); % [gC day^-1] Denominator for the ingestion function                   
    NP1 = P.IDP + P.EDPp.*pref('predcop','phyto').*repmat(P.R',1,P.n)+P.EDPd.*repmat(pref('predcop','detritus'),1,P.n).*repmat(P.D',1,P.n); % [gC day^-1] Denominator for ingestion function of copepods during day
    NP0 = P.INP + P.ENPp.*pref('predcop','phyto').*repmat(P.R,P.n,1)+P.ENPd.*repmat(pref('predcop','detritus')',P.n,1).*repmat(P.D,P.n,1); % [gC day^-1] Denominator for the ingestion function                      
    %J: No denominator because functional response type I
    NM1 = P.IDM + P.EDMd.*repmat(pref('meso','detritus'),1,P.n).*repmat(P.D',1,P.n)+P.EDMC.*pref('meso','copepod').*repmat(Cday',1,P.n)+P.EDMP.*pref('meso','predcop').*repmat(Pday',1,P.n); % [gC day^-1] Denominator for ingestion function of mesopelagic during day
    NM0 = P.INM + P.ENMd.*repmat(pref('meso','detritus')',P.n,1).*repmat(P.D,P.n,1)+P.ENMC.*pref('meso','copepod').*repmat(Cnight,P.n,1)+P.ENMP.*pref('meso','predcop').*repmat(Pnight,P.n,1); % [gC day^-1] Denominator for the ingestion function                   
    NB1 = P.IDB + P.EDBd.*repmat(pref('bathy','detritus'),1,P.n).*repmat(P.D',1,P.n)+P.EDBb.*pref('bathy','benthos').*repmat(P.Benthos',1,P.n)+...
                  P.EDBC.*pref('bathy','copepod').*repmat(Cday',1,P.n)+P.EDBP.*pref('bathy','predcop').*repmat(Pday',1,P.n)+P.EDBM.*pref('bathy','meso').*repmat(Mday',1,P.n); % [gC day^-1] Denominator for ingestion function of top predator during day
    NB0 = P.INB + P.ENBd.*repmat(pref('bathy','detritus')',P.n,1).*repmat(P.D,P.n,1)+P.ENBb.*pref('bathy','benthos').*repmat(P.Benthos,P.n,1)+...
                  P.ENBC.*pref('bathy','copepod').*repmat(Cnight,P.n,1)+P.ENBP.*pref('bathy','predcop').*repmat(Pnight,P.n,1)+P.ENBM.*pref('bathy','meso').*repmat(Mnight,P.n,1); % [gC day^-1] Denominator for the ingestion function  


%Ingestion rates
    IFC1 = P.IDF.*P.EDFC*pref('forage','copepod').* repmat(Cday' ,1,P.n)./NF1; % [gC day^-1] Ingestion rate of copepods during daytime by forage fish
    IFC0 = P.INF.*P.ENFC*pref('forage','copepod').* repmat(Cnight,P.n,1)./NF0;
    IFP1 = P.IDP.*P.EDFP*pref('forage','predcop').* repmat(Pday' ,1,P.n)./NF1; % [gC day^-1] Ingestion rate of copepods during daytime by forage fish
    IFP0 = P.INP.*P.ENFP*pref('forage','predcop').* repmat(Pnight,P.n,1)./NF0;
    IFD1 = P.IDF.*P.EDFd*repmat(pref('forage','detritus'),1,P.n).*repmat(P.D',1,P.n)  ./NF1;
    IFD0 = P.INF.*P.ENFd*repmat(pref('forage','detritus')',P.n,1).*repmat(P.D ,P.n,1)  ./NF0;
    IFb1 = P.IDF.*P.EDFb*pref('forage','benthos') .*repmat(P.Benthos',1,P.n)./NF1; % [gC day^-1] small b to show that it is for benthos and not bathypelagic fish
    IFb0 = P.INF.*P.ENFb*pref('forage','benthos') .*repmat(P.Benthos, P.n,1)./NF0;
    IFM1 = P.IDF.*P.EDFM*pref('forage','meso')    .*repmat(Mday' ,1,P.n)./NF1;
    IFM0 = P.INF.*P.ENFM*pref('forage','meso')    .*repmat(Mnight,P.n,1)./NF0;

    IAF1 = P.IDA.*P.EDAF*pref('top','forage').* repmat(Fday' ,1,P.n)./NA1; % [gC day^-1] Ingestion rate of forage fish during daytime by top predators
    IAF0 = P.INA.*P.ENAF*pref('top','forage').* repmat(Fnight,P.n,1)./NA0;
    IAJ1 = P.IDA.*P.EDAJ*pref('top','tactile').*repmat(Jday',1,P.n)  ./NA1;
    IAJ0 = P.INA.*P.ENAJ*pref('top','tactile').*repmat(Jnight ,P.n,1)  ./NA0;
    IAM1 = P.IDA.*P.EDAM*pref('top','meso') .*repmat(Mday',1,P.n)./NA1; % [gC day^-1]
    IAM0 = P.INA.*P.ENAM*pref('top','meso') .*repmat(Mnight, P.n,1)./NA0;
    IAB1 = P.IDA.*P.EDAB*pref('top','bathy')    .*repmat(Bday' ,1,P.n)./NA1;
    IAB0 = P.INA.*P.ENAB*pref('top','bathy')    .*repmat(Bnight,P.n,1)./NA0;

    ICR1 = P.IDC.*P.EDCp*pref('copepod','phyto') .*repmat(P.R',1,P.n)./NC1; % [gC day^-1] Ingestion rate of phytoplankton during daytime by copepods
    ICR0 = P.INC.*P.ENCp*pref('copepod','phyto') .*repmat(P.R, P.n,1)./NC0;
    ICD1 = P.IDC.*P.EDCd*repmat(pref('copepod','detritus'),1,P.n).*repmat(P.D',1,P.n)  ./NC1;
    ICD0 = P.INC.*P.ENCd*repmat(pref('copepod','detritus')',P.n,1).*repmat(P.D ,P.n,1)  ./NC0;
    
    IPR1 = P.IDP.*P.EDPp*pref('predcop','phyto') .*repmat(P.R',1,P.n)./NP1; % [gC day^-1] Ingestion rate of phytoplankton during daytime by copepods
    IPR0 = P.INP.*P.ENPp*pref('predcop','phyto') .*repmat(P.R, P.n,1)./NP0;
    IPD1 = P.IDP.*P.EDPd*repmat(pref('predcop','detritus'),1,P.n).*repmat(P.D',1,P.n)  ./NP1;
    IPD0 = P.INP.*P.ENPd*repmat(pref('predcop','detritus')',P.n,1).*repmat(P.D ,P.n,1)  ./NP0;

    IJC1 = P.EDJC*pref('tactile','copepod').*repmat(Cday',1,P.n); % [gC day^-1] Ingestion rate of copepods during daytime by tactile predators - note the Type I functional response
    IJC0 = P.ENJC*pref('tactile','copepod').*repmat(Cnight,P.n,1);
    IJP1 = P.EDJP*pref('tactile','predcop').*repmat(Pday',1,P.n); % [gC day^-1] Ingestion rate of copepods during daytime by tactile predators - note the Type I functional response
    IJP0 = P.ENJP*pref('tactile','predcop').*repmat(Pnight,P.n,1);
    IJM1 = P.EDJM*pref('tactile','meso').*repmat(Mday',1,P.n);
    IJM0 = P.ENJM*pref('tactile','meso').*repmat(Mnight,P.n,1);

    IMC1 = P.IDM.*P.EDMC*pref('meso','copepod') .*repmat(Cday',1,P.n)./NM1; % [gC day^-1] Ingestion rate of copepods during daytime by mesopelagic fish
    IMC0 = P.INM.*P.ENMC*pref('meso','copepod') .*repmat(Cnight, P.n,1)./NM0;
    IMP1 = P.IDM.*P.EDMP*pref('meso','predcop') .*repmat(Pday',1,P.n)./NM1; % [gC day^-1] Ingestion rate of copepods during daytime by mesopelagic fish
    IMP0 = P.INM.*P.ENMP*pref('meso','predcop') .*repmat(Pnight, P.n,1)./NM0;
    IMD1 = P.IDM.*P.EDMd*repmat(pref('meso','detritus'),1,P.n).*repmat(P.D',1,P.n)  ./NM1; 
    IMD0 = P.INM.*P.ENMd*repmat(pref('meso','detritus')',P.n,1).*repmat(P.D ,P.n,1)  ./NM0; 

    IBD1 = P.IDB.*P.EDBd*repmat(pref('bathy','detritus'),1,P.n).* repmat(P.D' ,1,P.n)./NB1; % [gC day^-1] Ingestion rate of forage fish during daytime by top predators
    IBD0 = P.INB.*P.ENBd*repmat(pref('bathy','detritus')',P.n,1).* repmat(P.D,P.n,1)./NB0;
    IBb1 = P.IDB.*P.EDBb*pref('bathy','benthos').*repmat(P.Benthos',1,P.n)  ./NB1;
    IBb0 = P.INB.*P.ENBb*pref('bathy','benthos').*repmat(P.Benthos ,P.n,1)  ./NB0;
    IBC1 = P.IDB.*P.EDBC*pref('bathy','copepod') .*repmat(Cday',1,P.n)./NB1; % [gC day^-1]
    IBC0 = P.INB.*P.ENBC*pref('bathy','copepod') .*repmat(Cnight, P.n,1)./NB0;
    IBP1 = P.IDB.*P.EDBP*pref('bathy','predcop') .*repmat(Pday',1,P.n)./NB1; % [gC day^-1]
    IBP0 = P.INB.*P.ENBP*pref('bathy','predcop') .*repmat(Pnight, P.n,1)./NB0;
    IBM1 = P.IDB.*P.EDBM*pref('bathy','meso')    .*repmat(Mday' ,1,P.n)./NB1;
    IBM0 = P.INB.*P.ENBM*pref('bathy','meso')    .*repmat(Mnight,P.n,1)./NB0;
    
    
%Remove all the NaN of ingestion rates when they do not feed at all
    IFC1(isnan(IFC1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by forage fish
    IFC0(isnan(IFC0)) = 0;
    IFP1(isnan(IFP1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by forage fish
    IFP0(isnan(IFP0)) = 0;
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
    IPR1(isnan(IPR1)) = 0; % [gC day^-1] Ingestion rate of phytoplankton during daytime by copepods
    IPR0(isnan(IPR0)) = 0;
    IPD1(isnan(IPD1)) = 0;
    IPD0(isnan(IPD0)) = 0;
    IJC1(isnan(IJC1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by tactile predators - note the Type I functional response
    IJC0(isnan(IJC0)) = 0;
    IJP1(isnan(IJP1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by tactile predators - note the Type I functional response
    IJP0(isnan(IJP0)) = 0;
    IJM1(isnan(IJM1)) = 0;
    IJM0(isnan(IJM0)) = 0;
    IMC1(isnan(IMC1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by mesopelagic fish
    IMC0(isnan(IMC0)) = 0; 
    IMP1(isnan(IMP1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by mesopelagic fish
    IMP0(isnan(IMP0)) = 0; 
    IMD1(isnan(IMD1)) = 0; 
    IMD0(isnan(IMD0)) = 0;
    IBD1(isnan(IBD1)) = 0; % [gC day^-1] Ingestion rate of forage fish during daytime by top predators
    IBD0(isnan(IBD0)) = 0;
    IBb1(isnan(IBb1)) = 0;
    IBb0(isnan(IBb0)) = 0;
    IBC1(isnan(IBC1)) = 0; % [gC day^-1]
    IBC0(isnan(IBC0)) = 0;
    IBP1(isnan(IBP1)) = 0; % [gC day^-1]
    IBP0(isnan(IBP0)) = 0;
    IBM1(isnan(IBM1)) = 0;
    IBM0(isnan(IBM0)) = 0 ;
    

%Assimilation rates
    IC = (P.sigma*(P.fCR*ICR1+P.fCd*ICD1)+(1-P.sigma)*(P.fCR*ICR0+P.fCd*ICD0))/P.wC; % [day^-1] Total assimilation rate per individual per strategy for copepods
    IP = (P.sigma*(P.fPR*IPR1+P.fPd*IPD1)+(1-P.sigma)*(P.fPR*IPR0+P.fPd*IPD0))/P.wP; % [day^-1] Total assimilation rate per individual per strategy for copepods
    IF = P.fF*(P.sigma*(IFD1+IFb1+IFC1+IFD1+IFM1)+(1-P.sigma)*(IFD0+IFb0+IFC0+IFP0+IFM0))/P.wF; % [day^-1] Total assimilation rate per individual per strategy for forage fish
    IA = P.fA*(P.sigma*(IAF1+IAJ1+IAM1+IAB1)+(1-P.sigma)*(IAF0+IAJ0+IAM0+IAB0))/P.wA; % [day^-1] Total assimilation rate per individual per strategy for top predator
    IJ = P.fJ*(P.sigma*(IJC1+IJP1+IJM1)+(1-P.sigma)*(IJC0+IJP0+IJM0))/P.wJ; % [day^-1] Total assimilation rate per individual per strategy for tactile predator
    IM = (P.sigma*(P.fMd*IMD1+P.fMC*IMC1+P.fMC*IMP1)+(1-P.sigma)*(P.fMd*IMD0+P.fMC*IMC0+P.fMC*IMP0))/P.wM; % [day^-1] Total assimilation rate per individual per strategy for mesopelagic fish
    IB = P.fB*(P.sigma*(IBD1+IBb1+IBC1+IBP1+IBM1)+(1-P.sigma)*(IBD0+IBb0+IBC0+IBP0+IBM0))/P.wB; % [day^-1] Total assimilation rate per individual per strategy for bathypelagic fish

%Mortality rates due to predation - we calculate the mortality rate imposed by all strategies (predators) at each depth before redistributing it equally among prey
%Copepod
    mCday = (P.IDF.*P.EDFC*pref('forage','copepod')./NF1*P.n^2*P.F.*F/P.wF +...
             P.EDJC*pref('tactile','copepod')*P.n^2*P.J.*J/P.wJ + ...
             P.IDM.*P.EDMC*pref('meso','copepod')  ./NM1*P.n^2*P.M.*M/P.wM +...
             P.IDB.*P.EDBC*pref('bathy','copepod') ./NB1*P.n^2*P.B.*B/P.wB); % [day^-1] size n*n How much each strategy eats copepods during daytime
   
    mCday(isnan(mCday)) = 0;
         
    mCD = sum(mCday,2)'; % [day^-1] size 1*n What is the mortality rate experienced at each depth during day
    MortDa = repmat(mCD',1,P.n); % [day^-1] Mortality rate experienced by the different copepod strategies during daytime

    mCnight = (P.INF.*P.ENFC*pref('forage','copepod')./NF0*P.n^2*P.F.*F/P.wF +...
               P.ENJC*pref('tactile','copepod')*P.n^2*P.J.*J/P.wJ + ...
               P.INM.*P.ENMC*pref('meso','copepod') ./NM0*P.n^2*P.M.*M/P.wM + ...
               P.INB.*P.ENBC*pref('bathy','copepod') ./NB0*P.n^2*P.B.*B/P.wB); % [day^-1]
    
    mCnight(isnan(mCnight)) = 0;
    
    mCN = sum(mCnight,1); % [day^-1]
    MortNi = repmat(mCN,P.n,1); % [day^-1] Mortality rate experienced by the different bathypelagic fish strategies during nighttime

    MortC = minimortC+ P.sigma*MortDa + (1-P.sigma)*MortNi; % [day^-1] Total mortality rate experienced by the different copepod strategies

%Predatory Copepod
    mPday = (P.IDF.*P.EDFP*pref('forage','predcop')./NF1*P.n^2*P.F.*F/P.wF +...
             P.EDJP*pref('tactile','predcop')*P.n^2*P.J.*J/P.wJ + ...
             P.IDM.*P.EDMP*pref('meso','predcop')  ./NM1*P.n^2*P.M.*M/P.wM +...
             P.IDB.*P.EDBP*pref('bathy','predcop') ./NB1*P.n^2*P.B.*B/P.wB); % [day^-1] size n*n How much each strategy eats copepods during daytime
   
    mPday(isnan(mPday)) = 0;
         
    mPD = sum(mPday,2)'; % [day^-1] size 1*n What is the mortality rate experienced at each depth during day
    MortDa = repmat(mPD',1,P.n); % [day^-1] Mortality rate experienced by the different copepod strategies during daytime

    mPnight = (P.INF.*P.ENFP*pref('forage','predcop')./NF0*P.n^2*P.F.*F/P.wF +...
               P.ENJP*pref('tactile','predcop')*P.n^2*P.J.*J/P.wJ + ...
               P.INM.*P.ENMP*pref('meso','predcop') ./NM0*P.n^2*P.M.*M/P.wM + ...
               P.INB.*P.ENBP*pref('bathy','predcop') ./NB0*P.n^2*P.B.*B/P.wB); % [day^-1]
    
    mPnight(isnan(mPnight)) = 0;
    
    mPN = sum(mPnight,1); % [day^-1]
    MortNi = repmat(mPN,P.n,1); % [day^-1] Mortality rate experienced by the different bathypelagic fish strategies during nighttime

    MortP = minimort+ P.sigma*MortDa + (1-P.sigma)*MortNi; % [day^-1] Total mortality rate experienced by the different copepod strategies

%Forage fish
    mFday = (P.IDA.*P.EDAF*pref('top','forage')./NA1*P.n^2*P.A.*A/P.wA); % [day^-1] size n*n How much each strategy eats forage fish during daytime
    mFday(isnan(mFday)) = 0;
    mFD = sum(mFday,2)'; % [day^-1] size 1*n What is the mortality rate experienced at each depth during day
    MortDa = repmat(mFD',1,P.n); % [day^-1] Mortality rate experienced by the different forage fish strategies during daytime

    mFnight = (P.INA.*P.ENAF*pref('top','forage')./NA0*P.n^2*P.A.*A/P.wA); % [day^-1]
    mFnight(isnan(mFnight)) = 0;
    mFN = sum(mFnight,1); % [day^-1]
    MortNi = repmat(mFN,P.n,1); % [day^-1] Mortality rate experienced by the different forage fish strategies during nighttime

    MortF = minimort + P.sigma*MortDa + (1-P.sigma)*MortNi; % [day^-1] Total mortality rate experienced by the different forage fish strategies

%Tactile predator
    mJday = (P.IDA.*P.EDAJ*pref('top','tactile')./NA1*P.n^2*P.A.*A/P.wA); % [day^-1] size n*n How much each strategy eats tactile pred during daytime
    mJday(isnan(mJday)) = 0;
    mJD = sum(mJday,2)'; % [day^-1] size 1*n What is the mortality rate experienced at each depth during day
    MortDa = repmat(mJD',1,P.n); % [day^-1] Mortality rate experienced by the different tactile pred strategies during daytime

    mJnight = (P.INA.*P.ENAJ*pref('top','tactile')./NA0*P.n^2*P.A.*A/P.wA); % [day^-1]
    mJnight(isnan(mJnight)) = 0;
    mJN = sum(mJnight,1); % [day^-1]
    MortNi = repmat(mJN,P.n,1); % [day^-1] Mortality rate experienced by the different tactile pred strategies during nighttime

    MortJ = minimort + P.sigma*MortDa + (1-P.sigma)*MortNi; % [day^-1] Total mortality rate experienced by the different tactile predator strategies

%Mesopelagic fish
    mMday = (P.IDF.*P.EDFM*pref('forage','meso')./NF1*P.n^2*P.F.*F/P.wF +...
             P.EDJM*pref('tactile','meso')*P.n^2*P.J.*J/P.wJ +...
             P.IDA.*P.EDAM*pref('top','meso')./NA1*P.n^2*P.A.*A/P.wA +...
             P.IDB.*P.EDBM*pref('bathy','meso')./NB1*P.n^2*P.B.*B/P.wB); %./repmat(Mday' ,1,P.n); % [day^-1] size n*n How much each strategy eats mesopelagic during daytime
    mMday(isnan(mMday)) = 0;
    mMD = sum(mMday,2)'; % [day^-1] size 1*n What is the mortality rate experienced at each depth during day
    MortDa = repmat(mMD',1,P.n); % [day^-1] Mortality rate experienced by the different mesopelagic fish strategies during daytime

    mMnight = (P.INF.*P.ENFM*pref('forage','meso')./NF0*P.n^2*P.F.*F/P.wF +...
               P.ENJM*pref('tactile','meso')*P.n^2*P.J.*J/P.wJ + ...
               P.INA.*P.ENAM*pref('top','meso')./NA0*P.n^2*P.A.*A/P.wA +...
               P.INB.*P.ENBM*pref('bathy','meso')./NB0*P.n^2*P.B.*B/P.wB); % [day^-1]
    mMnight(isnan(mMnight)) = 0;
    mMN = sum(mMnight,1); % [day^-1]
    MortNi = repmat(mMN,P.n,1); % [day^-1] Mortality rate experienced by the different bathypelagic fish strategies during nighttime

    MortM = minimort+ P.sigma*MortDa + (1-P.sigma)*MortNi; % [day^-1] Total mortality rate experienced by the different copepod strategies

%Bathypelagic fish
    mBday = (P.IDA.*P.EDAM*pref('top','bathy')./NA1*P.n^2*P.A.*A/P.wA); % [day^-1] size n*n How much each strategy eats bathypelagic fish during daytime
    mBday(isnan(mBday)) = 0;
    mBD = sum(mBday,2)'; % [day^-1] size 1*n What is the mortality rate experienced at each depth during day
    MortDa = repmat(mBD',1,P.n); % [day^-1] Mortality rate experienced by the different bathypelagic fish strategies during daytime

    mBnight = (P.INA.*P.ENAM*pref('top','bathy')./NA0*P.n^2*P.A.*A/P.wA); % [day^-1]
    mBnight(isnan(mBnight)) = 0;
    mBN = sum(mBnight,1); % [day^-1]
    MortNi = repmat(mBN,P.n,1); % [day^-1] Mortality rate experienced by the different bathypelagic fish strategies during nighttime

    MortB = minimort + P.sigma*MortDa + (1-P.sigma)*MortNi; % [day^-1] Total mortality rate experienced by the different bathypelagic fish strategies


%Fitnesses
    fitA = (IA  - P.CA/P.wA - P.metA) ./(0.01*((P.sigma*Aday'+(1-P.sigma)*Anight)/(P.n*P.A)).^2); % [day^-1] Fitness of top predator - Frequency-dependent mortality rate
    fitC = (IC  - P.CC/P.wC - P.metC)./MortC; %-(0.001*((P.sigma*Cday'+(1-P.sigma)*Cnight)/(P.n*P.C)).^2);%- 0.2 ; % [day^-1]
    fitP = (IP  - P.CP/P.wP - P.metP)./MortP; %-(0.001*((P.sigma*Pday'+(1-P.sigma)*Pnight)/(P.n*P.P)).^2);%- 0.2 ; % [day^-1]
    fitJ = (IJ  - P.CJ/P.wJ - P.metJ)./MortJ; %-(0.001*((P.sigma*Jday'+(1-P.sigma)*Jnight)/(P.n*P.J)).^2);%-0.1 ; % [day^-1]
    fitF = (IF  - P.CF/P.wF - P.metF)./MortF; %-(0.001*((P.sigma*Fday'+(1-P.sigma)*Fnight)/(P.n*P.F)).^2);%-0.05; % [day^-1]
    fitM = (IM  - P.CM/P.wM - P.metM)./MortM; %-(0.001*((P.sigma*Mday'+(1-P.sigma)*Mnight)/(P.n*P.M)).^2);%-0.05; % [day^-1]
    fitB = (IB  - P.CB/P.wB - P.metB)./MortB; %-(0.001*((P.sigma*Bday'+(1-P.sigma)*Bnight)/(P.n*P.B)).^2);%-0.03; % [day^-1]
    
%     GA = (IA  - P.CA/P.wA - P.metA);
%     MA = (0.01*((P.sigma*Aday'+(1-P.sigma)*Anight)/(P.n*P.A)).^2);
%     GC = (IC  - P.CC/P.wC - P.metC);
%     GP = (IP  - P.CP/P.wP - P.metP);
%     GJ = (IJ  - P.CJ/P.wJ - P.metJ);
%     GF = (IF  - P.CF/P.wF - P.metF);
%     GM = (IM  - P.CM/P.wM - P.metM);
%     GB = (IB  - P.CB/P.wB - P.metB);
%     
%     fitA(fitA<0) = GA(fitA<0)./MA(fitA<0);
%     fitC(fitC<0) = GC(fitC<0)./MortC(fitC<0);
%     fitP(fitP<0) = GP(fitP<0)./MortP(fitP<0);
%     fitJ(fitJ<0) = GJ(fitJ<0)./MortJ(fitJ<0);
%     fitF(fitF<0) = GF(fitF<0)./MortF(fitF<0);
%     fitM(fitM<0) = GM(fitM<0)./MortM(fitM<0);
%     fitB(fitB<0) = GB(fitB<0)./MortB(fitB<0);
    
    
    
    if i<Niter-1
    x = 0.8;
    fitA = (1-x)*fitA + x*fita2;
    fitB = (1-x)*fitB + x*fitb2;
    fitC = (1-x)*fitC + x*fitc2;
    fitF = (1-x)*fitF + x*fitf2;
    fitM = (1-x)*fitM + x*fitm2;  
    fitJ = (1-x)*fitJ + x*fitj2;
    fitP = (1-x)*fitP + x*fitp2;
    end
    
% %Applying a Gaussian filter - tries
%     fitA = imgaussfilt(fitA,sigf);
%     fitC = imgaussfilt(fitC,sigf);
%     fitP = imgaussfilt(fitP,sigf);
%     fitB = imgaussfilt(fitB,sigf);
%     fitF = imgaussfilt(fitF,sigf);
%     fitJ = imgaussfilt(fitJ,sigf);
%     fitM = imgaussfilt(fitM,sigf);    
    
%     
%     fitF = sign(fitF).*log10(1+abs(fitF.*P.MaskF)); % Transformation to make it flatter - just a try for now
%     fitA = sign(fitA).*log10(1+abs(fitA.*P.MaskA));
%     fitB = sign(fitB).*log10(1+abs(fitB.*P.MaskB));
%     fitC = sign(fitC).*log10(1+abs(fitC.*P.MaskC));
%     fitP = sign(fitP).*log10(1+abs(fitP.*P.MaskP));
%     fitM = sign(fitM).*log10(1+abs(fitM.*P.MaskM));
%     fitJ = sign(fitJ).*log10(1+abs(fitJ.*P.MaskJ));

% sigmo = @(x) 1./(1+exp(-50*x))-1/2;
%     fitF = sigmo(fitF);
%     fitA = sigmo(fitA);
%     fitM = sigmo(fitM);
%     fitB = sigmo(fitB);
%     fitJ = sigmo(fitJ);
%     fitC = sigmo(fitC);
%     fitP = sigmo(fitP);

%%% NOW IS THE REPLICATOR PART
    FAmax = max(max(fitA)); FAmin = min(min(fitA));
    FBmax = max(max(fitB)); FBmin = min(min(fitB));
    FCmax = max(max(fitC)); FCmin = min(min(fitC));
    FPmax = max(max(fitP)); FPmin = min(min(fitP));
    FJmax = max(max(fitJ)); FJmin = min(min(fitJ));
    FMmax = max(max(fitM)); FMmin = min(min(fitM));
    FFmax = max(max(fitF)); FFmin = min(min(fitF));


%the proportionality factor is dynamic so that maximum increase is at most 2% per time step
     factA = abs(dtfact/max([FAmax]));%, -FAmin]); % [day]
     factB = abs(dtfact/max([FBmax]));%, -FBmin]);
     factC = abs(dtfact/max([FCmax]));%, -FCmin]);
     factP = abs(dtfact/max([FPmax]));%, -FPmin]);
     factJ = abs(dtfact/max([FJmax]));%, -FJmin]);
     factM = abs(dtfact/max([FMmax]));%, -FMmin]);
     factF = abs(dtfact/max([FFmax]));%, -FFmin]);

%increment, the core of the replicator equation
    A = A.*(1 + factA*fitA.*P.MaskA); % [-] Proportion of all strategies, before renormalization
    B = B.*(1 + factB*fitB.*P.MaskB);
    C = C.*(1 + factC*fitC.*P.MaskC);
    PC = PC.*(1 + factP*fitP.*P.MaskP);
    J = J.*(1 + factJ*fitJ.*P.MaskJ);
    M = M.*(1 + factM*fitM.*P.MaskM);
    F = F.*(1 + factF*fitF.*P.MaskF);
    
%no extinction, so that every strategy can still emerge
    A(A<dA0) = dA0; % [-]
    B(B<dB0) = dB0;
    C(C<dC0) = dC0;
    PC(PC<dP0) = dP0;
    J(J<dJ0) = dJ0;
    M(M<dM0) = dM0;
    F(F<dF0) = dF0;
%     %Renormalization
%     A = A/sum(sum(A)); % [-] Good matrices of strategies after each replicator time step
%     B = B/sum(sum(B));
%     C = C/sum(sum(C));
%     M = M/sum(sum(M));
%     J = J/sum(sum(J));
%     F = F/sum(sum(F));  
%     PC = PC/sum(sum(PC)); 
%     
% %Applying a Gaussian filter - tries
%     A = imgaussfilt(A,sig);
%     C = imgaussfilt(C,sig);
%     PC = imgaussfilt(PC,sig);
%     B = imgaussfilt(B,sig);
%     F = imgaussfilt(F,sig);
%     J = imgaussfilt(J,sig);
%     M = imgaussfilt(M,sig);
%     
    
    %Removing the impossible places -- just a try for now
    F(P.MaskF==0) = 0;
    A(P.MaskA==0) = 0;
    C(P.MaskC==0) = 0;
    PC(P.MaskP==0) = 0;
    J(P.MaskJ==0) = 0;
    B(P.MaskB==0) = 0;
    M(P.MaskM==0) = 0;

    %Renormalization
    A = A/sum(sum(A)); % [-] Good matrices of strategies after each replicator time step
    B = B/sum(sum(B));
    C = C/sum(sum(C));
    M = M/sum(sum(M));
    J = J/sum(sum(J));
    F = F/sum(sum(F));
    PC = PC/sum(sum(PC));
%Calculation of the day and night concentrations
     Cday = P.n*P.C*sum(C,2)'; % [gC m^-3] Average concentration in each layer during day for copepod
     Cnight = P.n*P.C*sum(C,1); % [gC m^-3] Average concentration in each layer during night
     Pday = P.n*P.P*sum(PC,2)'; % [gC m^-3] Average concentration in each layer during day for copepod
     Pnight = P.n*P.P*sum(PC,1); % [gC m^-3] Average concentration in each layer during night
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
     
     %Try for now fitnesses are also from the previous time step
     fita2 = fitA;
     fitb2 = fitB;
     fitc2 = fitC;
     fitf2 = fitF;
     fitj2 = fitJ;
     fitm2 = fitM;
     fitp2 = fitP;
  
     
     
%Save historic of convergence
    if i<Iavg
        MAday(:,Iavg-i) = Aday;
        MBday(:,Iavg-i) = Bday;
        MCday(:,Iavg-i) = Cday;
        MPday(:,Iavg-i) = Pday;
        MMday(:,Iavg-i) = Mday;
        MFday(:,Iavg-i) = Fday;
        MJday(:,Iavg-i) = Jday;
        MAnight(:,Iavg-i) = Anight;
        MBnight(:,Iavg-i) = Bnight;
        MCnight(:,Iavg-i) = Cnight;
        MPnight(:,Iavg-i) = Pnight;
        MMnight(:,Iavg-i) = Mnight;
        MFnight(:,Iavg-i) = Fnight;
        MJnight(:,Iavg-i) = Jnight;
        
        FitA(Iavg-i) = max(max(fitA));
        FitB(Iavg-i) = max(max(fitB));
        FitC(Iavg-i) = max(max(fitC));
        FitP(Iavg-i) = max(max(fitP));
        FitJ(Iavg-i) = max(max(fitJ));
        FitF(Iavg-i) = max(max(fitF));
        FitM(Iavg-i) = max(max(fitM));
    end
    
     notdone = (i>0);
     
%      if i==Niter-1
%          A = max(fitA,dA0); A = A/sum(sum(A));
%          B = max(fitB,dB0); B = B/sum(sum(B));
%          C = max(fitC,dC0); C = C/sum(sum(C)); 
%          F = max(fitF,dF0); F = F/sum(sum(F));
%          J = max(fitJ,dJ0); J = J/sum(sum(J));
%          M = max(fitM,dM0); M = M/sum(sum(M));
%          P = max(fitP,dP0); P = P/sum(sum(P));
%      end
end
toc
filename = strcat('Run_',num2str(k),'.mat');
save(filename)

Plot_DVM;
drawnow
end