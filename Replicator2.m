%Replicator code
%global C F J M A B %the global variables are the proportions of each population using each strategy
Niter = 10^5; % [-] number of iterations of the replicator equation
Iavg = Niter; %100000; % [-] How many of the last time steps do we save?
dtfact = 0.05; %Max percentage of change per time step
sig = 0.2; % Parameter used for the Gaussian filter default = 0.3
sigf = 0.7; % Parameter used for the Gaussian filter on the fitness matrices
reinit = 1; %Do we start from the last simulation or do we initialize strategy matrices?

minimort = 0.01; % [day^-1] Background mortality
minimortC = 0.1; % [day^-1] Background mortality for small copepods

for k=[2:4]
    P = Parameters3(k);

%Coefficient to prevent extinction of strategies - and bugs in the OMZ
coeff = 10^-12;%10^-6; % [-] 10^-5
dC0 = coeff*P.C; % [gC m^-3] Minimum concentration of organisms in a strategy
dF0 = coeff*P.F;
dJ0 = coeff*P.J;
dA0 = coeff*P.A;
dM0 = coeff*P.M;
dP0 = coeff*P.P;

%Initialization of the different strategy matrices
if reinit==1
    C = ones(P.n).*P.MaskC + dC0; C = C/sum(sum(C)); % [-] random distributions at first - we can also use ones(n) to have equal distribution
    PC = ones(P.n).*P.MaskP + dP0; PC = PC/sum(sum(PC));
    F = ones(P.n).*P.MaskF + dF0; F = F/sum(sum(F)); 
    M = ones(P.n).*P.MaskM + dM0; M = M/sum(sum(M));
    A = ones(P.n).*P.MaskA + dA0; A = A/sum(sum(A));
    J = ones(P.n).*P.MaskJ + dJ0; J = J/sum(sum(J));
    D = zeros(P.n,7); % first time step will be without detritus, but that's ok as they come later
    
    A(P.MaskA==0) = 0;
    PC(P.MaskP==0) = 0;
    C(P.MaskC==0) = 0;
    F(P.MaskF==0) = 0;
    M(P.MaskM==0) = 0;
    J(P.MaskJ==0) = 0;
end

% Saved fitnesses for the last time steps
FitC = zeros(1,Iavg); % [day^-1] 
FitP = zeros(1,Iavg);
FitF = zeros(1,Iavg);
FitA = zeros(1,Iavg);
FitJ = zeros(1,Iavg);
FitM = zeros(1,Iavg);
%Saved distributions for the last time steps
MAday = zeros(P.n,Iavg); % [-]
MAnight = MAday;
MCday = MAday; MCnight = MAday;
MJday = MAday; MJnight = MAday;
MMday = MAday; MMnight = MAday;
MFday = MAday; MFnight = MAday;
MPday = MAday; MPnight = MAday;

MA = zeros(P.n);
MC = MA; MJ = MA; MM = MA; MF = MA; MP = MA;

MD = zeros(P.n,7,Iavg);


Cday = P.n*P.C*sum(C,2)'; % [gC m^-3] Average concentration in each layer during day for copepod
Cnight = P.n*P.C*sum(C,1); % [gC m^-3] Average concentration in each layer during night
Pday = P.n*P.P*sum(PC,2)'; % [gC m^-3] Average concentration in each layer during day for copepod
Pnight = P.n*P.P*sum(PC,1); % [gC m^-3] Average concentration in each layer during night
Fday = P.n*P.F*sum(F,2)'; % [gC m^-3] Average concentration in each layer during day for forage fish
Fnight = P.n*P.F*sum(F,1); % [gC m^-3] Average concentration in each layer during night
Aday = P.n*P.A*sum(A,2)'; % [gC m^-3] Average concentration in each layer during day for top predator
Anight = P.n*P.A*sum(A,1); % [gC m^-3] Average concentration in each layer during night
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
    NF1 = P.IDF + sum(P.EDFd.*repmat(pref('forage','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n),3)  +P.EDFC.*pref('forage','copepod').*repmat(Cday',1,P.n)+ P.EDFP.*pref('forage','predcop').*repmat(Pday',1,P.n)+...
                  P.EDFM.*pref('forage','meso').*repmat(Mday',1,P.n); % [gC day^-1] Denominator for ingestion function of forage fish during day
    NF0 = P.INF + sum(P.ENFd.*repmat(pref('forage','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n),[2,1,3]),3)      +P.ENFC.*pref('forage','copepod').*repmat(Cnight,P.n,1)+P.ENFP.*pref('forage','predcop').*repmat(Pnight,P.n,1)+...
                  P.ENFM.*pref('forage','meso').*repmat(Mnight,P.n,1); % [gC day^-1] Denominator for the ingestion function
    NA1 = P.IDA + P.EDAF.*pref('top','forage').*repmat(Fday',1,P.n)+P.EDAJ.*pref('top','tactile').*repmat(Jday',1,P.n)+...
                  P.EDAM.*pref('top','meso').*repmat(Mday',1,P.n); % [gC day^-1] Denominator for ingestion function of top predator during day
    NA0 = P.INA + P.ENAF.*pref('top','forage').*repmat(Fnight,P.n,1)+P.ENAJ.*pref('top','tactile').*repmat(Jnight,P.n,1)+...
                  P.ENAM.*pref('top','meso').*repmat(Mnight,P.n,1); % [gC day^-1] Denominator for the ingestion function           
    NC1 = P.IDC + P.EDCp.*pref('copepod','phyto').*repmat(P.R',1,P.n)+sum(P.EDCd.*repmat(pref('copepod','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n),3); % [gC day^-1] Denominator for ingestion function of copepods during day
    NC0 = P.INC + P.ENCp.*pref('copepod','phyto').*repmat(P.R,P.n,1)+sum(P.ENCd.*repmat(pref('copepod','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n),[2,1,3]),3); % [gC day^-1] Denominator for the ingestion function                   
    NP1 = P.IDP + P.EDPp.*pref('predcop','phyto').*repmat(P.R',1,P.n)+sum(P.EDPd.*repmat(pref('predcop','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n),3); % [gC day^-1] Denominator for ingestion function of copepods during day
    NP0 = P.INP + P.ENPp.*pref('predcop','phyto').*repmat(P.R,P.n,1)+sum(P.ENPd.*repmat(pref('predcop','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n),[2,1,3]),3); % [gC day^-1] Denominator for the ingestion function                      
    %J: No denominator because functional response type I
    NM1 = P.IDM + sum(P.EDMd.*repmat(pref('meso','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n),3)+P.EDMC.*pref('meso','copepod').*repmat(Cday',1,P.n)+P.EDMP.*pref('meso','predcop').*repmat(Pday',1,P.n); % [gC day^-1] Denominator for ingestion function of mesopelagic during day
    NM0 = P.INM + sum(P.ENMd.*repmat(pref('meso','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n),[2,1,3]),3)+P.ENMC.*pref('meso','copepod').*repmat(Cnight,P.n,1)+P.ENMP.*pref('meso','predcop').*repmat(Pnight,P.n,1); % [gC day^-1] Denominator for the ingestion function                   
    
%Ingestion rates

    % First of detritus in case a rescaling is needed
     IFD1 = P.IDF.*P.EDFd.*repmat(pref('forage','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n)  ./NF1;
     IFD0 = P.INF.*P.ENFd.*repmat(pref('forage','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n),[2,1,3])  ./NF0;
     IMD1 = P.IDM.*P.EDMd.*repmat(pref('meso','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n)  ./NM1; 
     IMD0 = P.INM.*P.ENMd.*repmat(pref('meso','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n),[2,1,3])  ./NM0; 
     ICD1 = P.IDC.*P.EDCd.*repmat(pref('copepod','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n)  ./NC1;
     ICD0 = P.INC.*P.ENCd.*repmat(pref('copepod','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n),[2,1,3])  ./NC0;
     IPD1 = P.IDP.*P.EDPd.*repmat(pref('predcop','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n)  ./NP1;
     IPD0 = P.INP.*P.ENPd.*repmat(pref('predcop','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n),[2,1,3])  ./NP0;
     
        IFD1(isnan(IFD1)) = 0;
        IFD0(isnan(IFD0)) = 0;
        IMD1(isnan(IMD1)) = 0;
        IMD0(isnan(IMD0)) = 0;
        IPD1(isnan(IPD1)) = 0;
        IPD0(isnan(IPD0)) = 0;
        ICD1(isnan(ICD1)) = 0;
        ICD0(isnan(ICD0)) = 0;
     
        %Make sure that they do not eat more than there actually is
        ConsdayD = squeeze(sum(IFD1.*repmat(Fday',1,1,7)/P.wF+IPD1.*repmat(Pday',1,1,7)/P.wP+ICD1.*repmat(Cday',1,1,7)/P.wC+IMD1.*repmat(Mday',1,1,7)/P.wM,2)); 
        ConsnigD = squeeze(sum(permute(IFD0,[2 1 3]).*repmat(Fnight',1,1,7)/P.wF+permute(IPD0,[2 1 3]).*repmat(Pnight',1,1,7)/P.wP+permute(ICD0,[2 1 3]).*repmat(Cnight',1,1,7)/P.wC+...
                                permute(IMD0,[2 1 3]).*repmat(Mnight',1,1,7)/P.wM,2));
        ConsD = P.sigma*ConsdayD + (1-P.sigma)*ConsnigD;
        
        resc = ones(P.n,7);
        max_ing = 0.8; %maximum percentage of detritus that we allow to be eaten in one day
        for depth=1:P.n
            for detr = 1:7
                if ConsD(depth,detr) > max_ing*D(depth,detr)                  
                    resc(depth,detr) = max_ing*D(depth,detr)/ConsD(depth,detr);
                    IFD1(depth,:,detr) = IFD1(depth,:,detr)*resc(depth,detr);
                    IFD0(:,depth,detr) = IFD0(:,depth,detr)*resc(depth,detr);
                    IMD1(depth,:,detr) = IMD1(depth,:,detr)*resc(depth,detr);
                    IMD0(:,depth,detr) = IMD0(:,depth,detr)*resc(depth,detr);
                    ICD1(depth,:,detr) = ICD1(depth,:,detr)*resc(depth,detr);
                    ICD0(:,depth,detr) = ICD0(:,depth,detr)*resc(depth,detr);
                    IPD1(depth,:,detr) = IPD1(depth,:,detr)*resc(depth,detr);
                    IPD0(:,depth,detr) = IPD0(:,depth,detr)*resc(depth,detr);
                end
            end
        end
        
    
        ConsdayD = ConsdayD.*resc;
        ConsnigD = ConsnigD.*resc;
        ConsD = ConsD.*resc;
        
        RESC = repmat(reshape(resc,P.n,1,7),1,P.n,1);
        
        % Denominators again but with the good rescaling
    NF1 = P.IDF + sum(P.EDFd.*repmat(pref('forage','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n).*RESC,3)  +P.EDFC.*pref('forage','copepod').*repmat(Cday',1,P.n)+ P.EDFP.*pref('forage','predcop').*repmat(Pday',1,P.n)+...
                  P.EDFM.*pref('forage','meso').*repmat(Mday',1,P.n); % [gC day^-1] Denominator for ingestion function of forage fish during day
    NF0 = P.INF + sum(P.ENFd.*repmat(pref('forage','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n),[2,1,3]).*RESC,3)      +P.ENFC.*pref('forage','copepod').*repmat(Cnight,P.n,1)+P.ENFP.*pref('forage','predcop').*repmat(Pnight,P.n,1)+...
                  P.ENFM.*pref('forage','meso').*repmat(Mnight,P.n,1); % [gC day^-1] Denominator for the ingestion function
    NC1 = P.IDC + P.EDCp.*pref('copepod','phyto').*repmat(P.R',1,P.n)+sum(P.EDCd.*repmat(pref('copepod','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n).*RESC,3); % [gC day^-1] Denominator for ingestion function of copepods during day
    NC0 = P.INC + P.ENCp.*pref('copepod','phyto').*repmat(P.R,P.n,1)+sum(P.ENCd.*repmat(pref('copepod','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n),[2,1,3]).*RESC,3); % [gC day^-1] Denominator for the ingestion function                   
    NP1 = P.IDP + P.EDPp.*pref('predcop','phyto').*repmat(P.R',1,P.n)+sum(P.EDPd.*repmat(pref('predcop','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n).*RESC,3); % [gC day^-1] Denominator for ingestion function of copepods during day
    NP0 = P.INP + P.ENPp.*pref('predcop','phyto').*repmat(P.R,P.n,1)+sum(P.ENPd.*repmat(pref('predcop','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n),[2,1,3]).*RESC,3); % [gC day^-1] Denominator for the ingestion function                      
    NM1 = P.IDM + sum(P.EDMd.*repmat(pref('meso','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n).*RESC,3)+P.EDMC.*pref('meso','copepod').*repmat(Cday',1,P.n)+P.EDMP.*pref('meso','predcop').*repmat(Pday',1,P.n); % [gC day^-1] Denominator for ingestion function of mesopelagic during day
    NM0 = P.INM + sum(P.ENMd.*repmat(pref('meso','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n).*RESC,[2,1,3]),3)+P.ENMC.*pref('meso','copepod').*repmat(Cnight,P.n,1)+P.ENMP.*pref('meso','predcop').*repmat(Pnight,P.n,1); % [gC day^-1] Denominator for the ingestion function                   



    IFC1 = P.IDF.*P.EDFC*pref('forage','copepod').* repmat(Cday' ,1,P.n)./NF1; % [gC day^-1] Ingestion rate of copepods during daytime by forage fish
    IFC0 = P.INF.*P.ENFC*pref('forage','copepod').* repmat(Cnight,P.n,1)./NF0;
    IFP1 = P.IDP.*P.EDFP*pref('forage','predcop').* repmat(Pday' ,1,P.n)./NF1; % [gC day^-1] Ingestion rate of copepods during daytime by forage fish
    IFP0 = P.INP.*P.ENFP*pref('forage','predcop').* repmat(Pnight,P.n,1)./NF0;
    IFM1 = P.IDF.*P.EDFM*pref('forage','meso')    .*repmat(Mday' ,1,P.n)./NF1;
    IFM0 = P.INF.*P.ENFM*pref('forage','meso')    .*repmat(Mnight,P.n,1)./NF0;

    IAF1 = P.IDA.*P.EDAF*pref('top','forage').* repmat(Fday' ,1,P.n)./NA1; % [gC day^-1] Ingestion rate of forage fish during daytime by top predators
    IAF0 = P.INA.*P.ENAF*pref('top','forage').* repmat(Fnight,P.n,1)./NA0;
    IAJ1 = P.IDA.*P.EDAJ*pref('top','tactile').*repmat(Jday',1,P.n)  ./NA1;
    IAJ0 = P.INA.*P.ENAJ*pref('top','tactile').*repmat(Jnight ,P.n,1)  ./NA0;
    IAM1 = P.IDA.*P.EDAM*pref('top','meso') .*repmat(Mday',1,P.n)./NA1; % [gC day^-1]
    IAM0 = P.INA.*P.ENAM*pref('top','meso') .*repmat(Mnight, P.n,1)./NA0;
    
    ICR1 = P.IDC.*P.EDCp*pref('copepod','phyto') .*repmat(P.R',1,P.n)./NC1; % [gC day^-1] Ingestion rate of phytoplankton during daytime by copepods
    ICR0 = P.INC.*P.ENCp*pref('copepod','phyto') .*repmat(P.R, P.n,1)./NC0;
        
    IPR1 = P.IDP.*P.EDPp*pref('predcop','phyto') .*repmat(P.R',1,P.n)./NP1; % [gC day^-1] Ingestion rate of phytoplankton during daytime by copepods
    IPR0 = P.INP.*P.ENPp*pref('predcop','phyto') .*repmat(P.R, P.n,1)./NP0;
   
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

        
%Remove all the NaN of ingestion rates when they do not feed at all
    IFC1(isnan(IFC1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by forage fish
    IFC0(isnan(IFC0)) = 0;
    IFP1(isnan(IFP1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by forage fish
    IFP0(isnan(IFP0)) = 0;
    IFD1(isnan(IFD1)) = 0;
    IFD0(isnan(IFD0)) = 0;
    IFM1(isnan(IFM1)) = 0;
    IFM0(isnan(IFM0)) = 0;
    IAF1(isnan(IAF1)) = 0; % [gC day^-1] Ingestion rate of forage fish during daytime by top predators
    IAF0(isnan(IAF0)) = 0;
    IAJ1(isnan(IAJ1)) = 0;
    IAJ0(isnan(IAJ0)) = 0;
    IAM1(isnan(IAM1)) = 0 ; % [gC day^-1]
    IAM0(isnan(IAM0)) = 0;
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
    
    

%Assimilation rates
    IC = (P.sigma*(P.fCR*ICR1+P.fCd*sum(ICD1,3))+(1-P.sigma)*(P.fCR*ICR0+P.fCd*sum(ICD0,3)))/P.wC; % [day^-1] Total assimilation rate per individual per strategy for copepods
    IP = (P.sigma*(P.fPR*IPR1+P.fPd*sum(IPD1,3))+(1-P.sigma)*(P.fPR*IPR0+P.fPd*sum(IPD0,3)))/P.wP; % [day^-1] Total assimilation rate per individual per strategy for copepods
    IF = P.fF*(P.sigma*(IFP1+IFC1+sum(IFD1,3)+IFM1)+(1-P.sigma)*(sum(IFD0,3)+IFC0+IFP0+IFM0))/P.wF; % [day^-1] Total assimilation rate per individual per strategy for forage fish
    IA = P.fA*(P.sigma*(IAF1+IAJ1+IAM1)+(1-P.sigma)*(IAF0+IAJ0+IAM0))/P.wA; % [day^-1] Total assimilation rate per individual per strategy for top predator
    IJ = P.fJ*(P.sigma*(IJC1+IJP1+IJM1)+(1-P.sigma)*(IJC0+IJP0+IJM0))/P.wJ; % [day^-1] Total assimilation rate per individual per strategy for tactile predator
    IM = (P.sigma*(P.fMd*sum(IMD1,3)+P.fMC*IMC1+P.fMC*IMP1)+(1-P.sigma)*(P.fMd*sum(IMD0,3)+P.fMC*IMC0+P.fMC*IMP0))/P.wM; % [day^-1] Total assimilation rate per individual per strategy for mesopelagic fish
    
%Mortality rates due to predation - we calculate the mortality rate imposed by all strategies (predators) at each depth before redistributing it equally among prey
%Copepod
    mCday = (P.IDF.*P.EDFC*pref('forage','copepod')./NF1*P.n^2*P.F.*F/P.wF +...
             P.EDJC*pref('tactile','copepod')*P.n^2*P.J.*J/P.wJ + ...
             P.IDM.*P.EDMC*pref('meso','copepod')  ./NM1*P.n^2*P.M.*M/P.wM); % [day^-1] size n*n How much each strategy eats copepods during daytime
   
    mCday(isnan(mCday)) = 0;
         
    mCD = sum(mCday,2)'; % [day^-1] size 1*n What is the mortality rate experienced at each depth during day
    MortDa = repmat(mCD',1,P.n); % [day^-1] Mortality rate experienced by the different copepod strategies during daytime

    mCnight = (P.INF.*P.ENFC*pref('forage','copepod')./NF0*P.n^2*P.F.*F/P.wF +...
               P.ENJC*pref('tactile','copepod')*P.n^2*P.J.*J/P.wJ + ...
               P.INM.*P.ENMC*pref('meso','copepod') ./NM0*P.n^2*P.M.*M/P.wM); % [day^-1]
    
    mCnight(isnan(mCnight)) = 0;
    
    mCN = sum(mCnight,1); % [day^-1]
    MortNi = repmat(mCN,P.n,1); % [day^-1] Mortality rate experienced by the different bathypelagic fish strategies during nighttime

    MortC = minimortC+ P.sigma*MortDa + (1-P.sigma)*MortNi; % [day^-1] Total mortality rate experienced by the different copepod strategies

%Predatory Copepod
    mPday = (P.IDF.*P.EDFP*pref('forage','predcop')./NF1*P.n^2*P.F.*F/P.wF +...
             P.EDJP*pref('tactile','predcop')*P.n^2*P.J.*J/P.wJ + ...
             P.IDM.*P.EDMP*pref('meso','predcop')  ./NM1*P.n^2*P.M.*M/P.wM); % [day^-1] size n*n How much each strategy eats copepods during daytime
   
    mPday(isnan(mPday)) = 0;
         
    mPD = sum(mPday,2)'; % [day^-1] size 1*n What is the mortality rate experienced at each depth during day
    MortDa = repmat(mPD',1,P.n); % [day^-1] Mortality rate experienced by the different copepod strategies during daytime

    mPnight = (P.INF.*P.ENFP*pref('forage','predcop')./NF0*P.n^2*P.F.*F/P.wF +...
               P.ENJP*pref('tactile','predcop')*P.n^2*P.J.*J/P.wJ + ...
               P.INM.*P.ENMP*pref('meso','predcop') ./NM0*P.n^2*P.M.*M/P.wM); % [day^-1]
    
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
             P.IDA.*P.EDAM*pref('top','meso')./NA1*P.n^2*P.A.*A/P.wA); %./repmat(Mday' ,1,P.n); % [day^-1] size n*n How much each strategy eats mesopelagic during daytime
    mMday(isnan(mMday)) = 0;
    mMD = sum(mMday,2)'; % [day^-1] size 1*n What is the mortality rate experienced at each depth during day
    MortDa = repmat(mMD',1,P.n); % [day^-1] Mortality rate experienced by the different mesopelagic fish strategies during daytime

    mMnight = (P.INF.*P.ENFM*pref('forage','meso')./NF0*P.n^2*P.F.*F/P.wF +...
               P.ENJM*pref('tactile','meso')*P.n^2*P.J.*J/P.wJ + ...
               P.INA.*P.ENAM*pref('top','meso')./NA0*P.n^2*P.A.*A/P.wA); % [day^-1]
    mMnight(isnan(mMnight)) = 0;
    mMN = sum(mMnight,1); % [day^-1]
    MortNi = repmat(mMN,P.n,1); % [day^-1] Mortality rate experienced by the different bathypelagic fish strategies during nighttime

    MortM = minimort+ P.sigma*MortDa + (1-P.sigma)*MortNi; % [day^-1] Total mortality rate experienced by the different copepod strategies


%Fitnesses
    fitA = (IA  - P.CA/P.wA - P.metA) ./(0.01*((P.sigma*Aday'+(1-P.sigma)*Anight)/(P.n*P.A)).^2+10^-4); % [day^-1] Fitness of top predator - Frequency-dependent mortality rate
    fitC = (IC  - P.CC/P.wC - P.metC)./MortC; %-(0.001*((P.sigma*Cday'+(1-P.sigma)*Cnight)/(P.n*P.C)).^2);%- 0.2 ; % [day^-1]
    fitP = (IP  - P.CP/P.wP - P.metP)./MortP; %-(0.001*((P.sigma*Pday'+(1-P.sigma)*Pnight)/(P.n*P.P)).^2);%- 0.2 ; % [day^-1]
    fitJ = (IJ  - P.CJ/P.wJ - P.metJ)./MortJ; %-(0.001*((P.sigma*Jday'+(1-P.sigma)*Jnight)/(P.n*P.J)).^2);%-0.1 ; % [day^-1]
    fitF = (IF  - P.CF/P.wF - P.metF)./MortF; %-(0.001*((P.sigma*Fday'+(1-P.sigma)*Fnight)/(P.n*P.F)).^2);%-0.05; % [day^-1]
    fitM = (IM  - P.CM/P.wM - P.metM)./MortM; %-(0.001*((P.sigma*Mday'+(1-P.sigma)*Mnight)/(P.n*P.M)).^2);%-0.05; % [day^-1]
    
%     GA = (IA  - P.CA/P.wA - P.metA);
%     MA = (0.01*((P.sigma*Aday'+(1-P.sigma)*Anight)/(P.n*P.A)).^2);
%     GC = (IC  - P.CC/P.wC - P.metC);
%     GP = (IP  - P.CP/P.wP - P.metP);
%     GJ = (IJ  - P.CJ/P.wJ - P.metJ);
%     GF = (IF  - P.CF/P.wF - P.metF);
%     GM = (IM  - P.CM/P.wM - P.metM);
%     
%     fitA(fitA<0) = GA(fitA<0)./MA(fitA<0);
%     fitC(fitC<0) = GC(fitC<0)./MortC(fitC<0);
%     fitP(fitP<0) = GP(fitP<0)./MortP(fitP<0);
%     fitJ(fitJ<0) = GJ(fitJ<0)./MortJ(fitJ<0);
%     fitF(fitF<0) = GF(fitF<0)./MortF(fitF<0);
%     fitM(fitM<0) = GM(fitM<0)./MortM(fitM<0);
    
    
    
    if i<Niter-1
    x = 0.8;
    fitA = (1-x)*fitA + x*fita2;
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
%     fitF = imgaussfilt(fitF,sigf);
%     fitJ = imgaussfilt(fitJ,sigf);
%     fitM = imgaussfilt(fitM,sigf);    
    
%     
%     fitF = sign(fitF).*log10(1+abs(fitF.*P.MaskF)); % Transformation to make it flatter - just a try for now
%     fitA = sign(fitA).*log10(1+abs(fitA.*P.MaskA));
%     fitC = sign(fitC).*log10(1+abs(fitC.*P.MaskC));
%     fitP = sign(fitP).*log10(1+abs(fitP.*P.MaskP));
%     fitM = sign(fitM).*log10(1+abs(fitM.*P.MaskM));
%     fitJ = sign(fitJ).*log10(1+abs(fitJ.*P.MaskJ));

% sigmo = @(x) 1./(1+exp(-50*x))-1/2;
%     fitF = sigmo(fitF);
%     fitA = sigmo(fitA);
%     fitM = sigmo(fitM);
%     fitJ = sigmo(fitJ);
%     fitC = sigmo(fitC);
%     fitP = sigmo(fitP);


%%%%%%%%%%%%%%%%%% NOW FECAL PELLET FLUX CALCULATION %%%%%%%%%%%%%%%%%%%%%%

%first calculation of the production at each depth - for each size (k+1 sizes): k for copepods, 1 for fish

C0F = (1-P.fF)*(P.sigma*sum(sum(IFD1,3)+IFC1+IFP1+IFM1,2).*Fday'+(1-P.sigma)*sum(sum(IFD0,3)+IFC0+IFP0+IFM0)'.*Fnight')/P.wF ./ (P.alpha(:,5) + P.SR(5)/P.dZ); % [gC/m3] Detritus concentration where each fish is creating it - from formaula 7 Stamieszkin et al 2015
C0C = (P.sigma*((1-P.fCR)*sum(ICR1,2)+(1-P.fCd)*sum(sum(ICD1,3),2)).*Cday'+(1-P.sigma)*((1-P.fCR)*sum(ICR0)'+(1-P.fCd)*sum(sum(ICD0,3))').*Cnight')/P.wC ./ (P.alpha(:,2) + P.SR(2)/P.dZ);
C0P = (P.sigma*((1-P.fPR)*sum(IPR1,2)+(1-P.fPd)*sum(sum(IPD1,3),2)).*Pday'+(1-P.sigma)*((1-P.fPR)*sum(IPR0)'+(1-P.fPd)*sum(sum(IPD0,3))').*Pnight')/P.wP ./ (P.alpha(:,3) + P.SR(3)/P.dZ);
C0A = (1-P.fA)*(P.sigma*sum(IAF1+IAJ1+IAM1,2).*Aday'+(1-P.sigma)*sum(IAF0+IAJ0+IAM0)'.*Anight')/P.wA ./ (P.alpha(:,6) + P.SR(6)/P.dZ);
C0M = (P.sigma*(sum((1-P.fMd)*sum(IMD1,3)+(1-P.fMC)*IMC1+(1-P.fMC)*IMP1,2).*Mday')+(1-P.sigma)*sum((1-P.fMd)*sum(IMD0,3)+(1-P.fMC)*IMC0+(1-P.fMC)*IMP0)'.*Mnight')/P.wM ./ (P.alpha(:,4) + P.SR(4)/P.dZ);
C0J = (1-P.fJ)*(P.sigma*sum(IJC1+IJP1+IJM1,2).*Jday'+(1-P.sigma)*sum(IJC0+IJP0+IJM0)'.*Jnight')/P.wJ ./ (P.alpha(:,7) + P.SR(7)/P.dZ);

Dnew = [P.BD', zeros(P.n,6)];
        
D0 = [C0C C0P C0M C0F C0A C0J]; % [gC /m3] Concentration of detritus where it is produced, i.e. "beginning of Martin curves" - in the order C P M F A J

        %Calculation of the curves sinking from the sources     
        for j=2:7  %for each detritus sizes (i.e. each producing population)
                for ii=P.n:-1:1 % go backbward to not count detritus twice 
                    Dnew(ii:P.n,j) = Dnew(ii:P.n,j) + D0(ii,j-1).*exp(P.alpha(ii:P.n,j)/P.SR(j).*(-P.zi(ii:P.n)'+P.zi(ii)));
                end
        end
        
        %Removal of the detritus eaten previously
        for j=1:7  %for each detritus sizes (backgr+zpk + fish)
                for ii=P.n:-1:1 % go backbward to not count detritus twice 
                    Dnew(ii:P.n,j) = Dnew(ii:P.n,j) - ConsD(ii,j)*P.dZ/P.SR(j).*exp(P.alpha(ii:P.n,j)/P.SR(j).*(-P.zi(ii:P.n)'+P.zi(ii)));  
                end
        end

        Dnew(Dnew<0) = 0;
        
        pD = 0.9; %0.9   
        D = pD*D + (1-pD)*Dnew;  
        

%%% NOW IS THE REPLICATOR PART
    FAmax = max(max(fitA)); FAmin = min(min(fitA));
    FCmax = max(max(fitC)); FCmin = min(min(fitC));
    FPmax = max(max(fitP)); FPmin = min(min(fitP));
    FJmax = max(max(fitJ)); FJmin = min(min(fitJ));
    FMmax = max(max(fitM)); FMmin = min(min(fitM));
    FFmax = max(max(fitF)); FFmin = min(min(fitF));


%the proportionality factor is dynamic so that maximum increase is at most 2% per time step
     factA = abs(dtfact/max([FAmax]));%, -FAmin]); % [day]
     factC = abs(dtfact/max([FCmax]));%, -FCmin]);
     factP = abs(dtfact/max([FPmax]));%, -FPmin]);
     factJ = abs(dtfact/max([FJmax]));%, -FJmin]);
     factM = abs(dtfact/max([FMmax]));%, -FMmin]);
     factF = abs(dtfact/max([FFmax]));%, -FFmin]);

%increment, the core of the replicator equation
 %  A = A.*(1 + factA*fitA.*P.MaskA); % [-] Proportion of all strategies, before renormalization
    C = C.*(1 + factC*fitC.*P.MaskC);
    PC = PC.*(1 + factP*fitP.*P.MaskP);
    J = J.*(1 + factJ*fitJ.*P.MaskJ);
    M = M.*(1 + factM*fitM.*P.MaskM);
    F = F.*(1 + factF*fitF.*P.MaskF);
    
%no extinction, so that every strategy can still emerge
    A(A<dA0) = dA0; % [-]
    C(C<dC0) = dC0;
    PC(PC<dP0) = dP0;
    J(J<dJ0) = dJ0;
    M(M<dM0) = dM0;
    F(F<dF0) = dF0;
    
    
%A =    diag([linspace(1,0,25), linspace(0,0,25)]);   %%% just a try to see the results with top predator distributed like that
    
    
    
%     %Renormalization
%     A = A/sum(sum(A)); % [-] Good matrices of strategies after each replicator time step
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
    M(P.MaskM==0) = 0;

    %Renormalization
    A = A/sum(sum(A)); % [-] Good matrices of strategies after each replicator time step
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
     Mday = P.n*P.M*sum(M,2)'; % [gC m^-3] Average concentration in each layer during day for mesopelagic fish
     Mnight = P.n*P.M*sum(M,1); % [gC m^-3] Average concentration in each layer during night
     Jday = P.n*P.J*sum(J,2)'; % [gC m^-3] Average concentration in each layer during day for tactile predator
     Jnight = P.n*P.J*sum(J,1); % [gC m^-3] Average concentration in each layer during night
     
     %Try for now fitnesses are also from the previous time step
     fita2 = fitA;
     fitc2 = fitC;
     fitf2 = fitF;
     fitj2 = fitJ;
     fitm2 = fitM;
     fitp2 = fitP;
  
     
     
%Save historic of convergence
    if i<Iavg
        MAday(:,Iavg-i) = Aday;
        MCday(:,Iavg-i) = Cday;
        MPday(:,Iavg-i) = Pday;
        MMday(:,Iavg-i) = Mday;
        MFday(:,Iavg-i) = Fday;
        MJday(:,Iavg-i) = Jday;
        MAnight(:,Iavg-i) = Anight;
        MCnight(:,Iavg-i) = Cnight;
        MPnight(:,Iavg-i) = Pnight;
        MMnight(:,Iavg-i) = Mnight;
        MFnight(:,Iavg-i) = Fnight;
        MJnight(:,Iavg-i) = Jnight;
        
        FitA(Iavg-i) = max(max(fitA));
        FitC(Iavg-i) = max(max(fitC));
        FitP(Iavg-i) = max(max(fitP));
        FitJ(Iavg-i) = max(max(fitJ));
        FitF(Iavg-i) = max(max(fitF));
        FitM(Iavg-i) = max(max(fitM));
        
        MA = MA + A/Iavg;
        MC = MC + C/Iavg;
        MF = MF + F/Iavg;
        MJ = MJ + J/Iavg;
        MM = MM + M/Iavg;
        MP = MP + PC/Iavg;
        
        MD(:,:,Iavg-i) = D;
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
filename = strcat('Run_F_',num2str(k),'.mat');
save(filename)

Plot_DVM;
drawnow
end