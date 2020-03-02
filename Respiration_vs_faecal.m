% load('TEST_RUNS_lat_0_long_-140.mat')
% 
%  Niter= a; Iavg= b;  P= c; MAday= d; MAnight= e; MCday=f; MCnight=g; MPday=h; MPnight=i; MFday=j; MFnight=k; MJday=l; MJnight=m;
%  MMday=n; MMnight=o; DegPOC_depth=p; DIC_dept=q; Dmean=r; MA=s; MC=t; MF=u; MJ=v; MM=w; MP=x; MD=y; FitA=z; FitC=aa; FitF=bb; FitJ=cc; FitM=dd; FitP=ee;
 
 imean = Niter/2;
MAday = mean(MAday(:,end-imean:end),2);
MAnight = mean(MAnight(:,end-imean:end),2);
MMday = mean(MMday(:,end-imean:end),2);
MMnight = mean(MMnight(:,end-imean:end),2);
MCday = mean(MCday(:,end-imean:end),2);
MCnight = mean(MCnight(:,end-imean:end),2);
MPday = mean(MPday(:,end-imean:end),2);
MPnight = mean(MPnight(:,end-imean:end),2);
MJday = mean(MJday(:,end-imean:end),2);
MJnight = mean(MJnight(:,end-imean:end),2);
MFday = mean(MFday(:,end-imean:end),2);
MFnight = mean(MFnight(:,end-imean:end),2);
 
 
 %% Respiration  [gC m^-2 day^-1]
 
 Resp = sum(DIC*P.dZ); 
 Resp_C = Resp(1);
 Resp_P = Resp(2);
 Resp_M = Resp(3);
 Resp_F = Resp(4);
 Resp_A = Resp(5);
 Resp_J = Resp(6);
 
%% Faecal pellet production

%Denominators for ingestion rates calculations
    NF1 = P.IDF + sum(P.EDFd.*repmat(pref('forage','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n),3)  +P.EDFC.*pref('forage','copepod').*repmat(MCday,1,P.n)+ P.EDFP.*pref('forage','predcop').*repmat(MPday,1,P.n)+...
                  P.EDFM.*pref('forage','meso').*repmat(MMday,1,P.n); % [gC day^-1] Denominator for ingestion function of forage fish during day
    NF0 = P.INF + sum(P.ENFd.*repmat(pref('forage','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n),[2,1,3]),3)      +P.ENFC.*pref('forage','copepod').*repmat(MCnight',P.n,1)+P.ENFP.*pref('forage','predcop').*repmat(MPnight',P.n,1)+...
                  P.ENFM.*pref('forage','meso').*repmat(MMnight',P.n,1); % [gC day^-1] Denominator for the ingestion function
    NA1 = P.IDA + P.EDAF.*pref('top','forage').*repmat(MFday,1,P.n)+P.EDAJ.*pref('top','tactile').*repmat(MJday,1,P.n)+...
                  P.EDAM.*pref('top','meso').*repmat(MMday,1,P.n); % [gC day^-1] Denominator for ingestion function of top predator during day
    NA0 = P.INA + P.ENAF.*pref('top','forage').*repmat(MFnight',P.n,1)+P.ENAJ.*pref('top','tactile').*repmat(MJnight',P.n,1)+...
                  P.ENAM.*pref('top','meso').*repmat(MMnight',P.n,1); % [gC day^-1] Denominator for the ingestion function           
    NC1 = P.IDC + P.EDCp.*pref('copepod','phyto').*repmat(P.R',1,P.n)+sum(P.EDCd.*repmat(pref('copepod','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n),3); % [gC day^-1] Denominator for ingestion function of copepods during day
    NC0 = P.INC + P.ENCp.*pref('copepod','phyto').*repmat(P.R,P.n,1)+sum(P.ENCd.*repmat(pref('copepod','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n),[2,1,3]),3); % [gC day^-1] Denominator for the ingestion function                   
    NP1 = P.IDP + P.EDPp.*pref('predcop','phyto').*repmat(P.R',1,P.n)+sum(P.EDPd.*repmat(pref('predcop','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n),3); % [gC day^-1] Denominator for ingestion function of copepods during day
    NP0 = P.INP + P.ENPp.*pref('predcop','phyto').*repmat(P.R,P.n,1)+sum(P.ENPd.*repmat(pref('predcop','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n),[2,1,3]),3); % [gC day^-1] Denominator for the ingestion function                      
    %J: No denominator because functional response type I
    NM1 = P.IDM + sum(P.EDMd.*repmat(pref('meso','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n),3)+P.EDMC.*pref('meso','copepod').*repmat(MCday,1,P.n)+P.EDMP.*pref('meso','predcop').*repmat(MPday,1,P.n); % [gC day^-1] Denominator for ingestion function of mesopelagic during day
    NM0 = P.INM + sum(P.ENMd.*repmat(pref('meso','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n),[2,1,3]),3)+P.ENMC.*pref('meso','copepod').*repmat(MCnight',P.n,1)+P.ENMP.*pref('meso','predcop').*repmat(MPnight',P.n,1); % [gC day^-1] Denominator for the ingestion function                   
    
%Ingestion rates

    % First of detritus in case a rescaling is needed
     IFD1 = P.IDF.*P.EDFd.*repmat(pref('forage','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n)  ./NF1;
     IFD0 = P.INF.*P.ENFd.*repmat(pref('forage','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n),[2,1,3])  ./NF0;
     IMD1 = P.IDM.*P.EDMd.*repmat(pref('meso','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n)  ./NM1; 
     IMD0 = P.INM.*P.ENMd.*repmat(pref('meso','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n),[2,1,3])  ./NM0; 
     ICD1 = P.IDC.*P.EDCd.*repmat(pref('copepod','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n)  ./NC1;
     ICD0 = P.INC.*P.ENCd.*repmat(pref('copepod','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n),[2,1,3])  ./NC0;
     IPD1 = P.IDP.*P.EDPd.*repmat(pref('predcop','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n)  ./NP1;
     IPD0 = P.INP.*P.ENPd.*repmat(pref('predcop','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n),[2,1,3])  ./NP0;
     
        IFD1(isnan(IFD1)) = 0; % [gC / day]
        IFD0(isnan(IFD0)) = 0;
        IMD1(isnan(IMD1)) = 0;
        IMD0(isnan(IMD0)) = 0;
        IPD1(isnan(IPD1)) = 0;
        IPD0(isnan(IPD0)) = 0;
        ICD1(isnan(ICD1)) = 0;
        ICD0(isnan(ICD0)) = 0;
     
        %Make sure that they do not eat more than there actually is
        ConsdayD = squeeze(sum(IFD1.*repmat(MFday,1,1,7)/P.wF+IPD1.*repmat(MPday,1,1,7)/P.wP+ICD1.*repmat(MCday,1,1,7)/P.wC+IMD1.*repmat(MMday,1,1,7)/P.wM,2)); 
        ConsnigD = squeeze(sum(permute(IFD0,[2 1 3]).*repmat(MFnight,1,1,7)/P.wF+permute(IPD0,[2 1 3]).*repmat(MPnight,1,1,7)/P.wP+permute(ICD0,[2 1 3]).*repmat(MCnight,1,1,7)/P.wC+...
                                permute(IMD0,[2 1 3]).*repmat(MMnight,1,1,7)/P.wM,2));
        ConsD = P.sigma*ConsdayD + (1-P.sigma)*ConsnigD;
        
        resc = ones(P.n,7);
        max_ing = 0.8; %maximum percentage of detritus that we allow to be eaten in one day
        for depth=1:P.n
            for detr = 1:7
                if ConsD(depth,detr) > max_ing*Dmean(depth,detr)                  
                    resc(depth,detr) = max_ing*Dmean(depth,detr)/ConsD(depth,detr);
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
    NF1 = P.IDF + sum(P.EDFd.*repmat(pref('forage','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n).*RESC,3)  +P.EDFC.*pref('forage','copepod').*repmat(MCday,1,P.n)+ P.EDFP.*pref('forage','predcop').*repmat(MPday,1,P.n)+...
                  P.EDFM.*pref('forage','meso').*repmat(MMday,1,P.n); % [gC day^-1] Denominator for ingestion function of forage fish during day
    NF0 = P.INF + sum(P.ENFd.*repmat(pref('forage','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n),[2,1,3]).*RESC,3)      +P.ENFC.*pref('forage','copepod').*repmat(MCnight',P.n,1)+P.ENFP.*pref('forage','predcop').*repmat(MPnight',P.n,1)+...
                  P.ENFM.*pref('forage','meso').*repmat(MMnight',P.n,1); % [gC day^-1] Denominator for the ingestion function
    NC1 = P.IDC + P.EDCp.*pref('copepod','phyto').*repmat(P.R',1,P.n)+sum(P.EDCd.*repmat(pref('copepod','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n).*RESC,3); % [gC day^-1] Denominator for ingestion function of copepods during day
    NC0 = P.INC + P.ENCp.*pref('copepod','phyto').*repmat(P.R,P.n,1)+sum(P.ENCd.*repmat(pref('copepod','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n),[2,1,3]).*RESC,3); % [gC day^-1] Denominator for the ingestion function                   
    NP1 = P.IDP + P.EDPp.*pref('predcop','phyto').*repmat(P.R',1,P.n)+sum(P.EDPd.*repmat(pref('predcop','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n).*RESC,3); % [gC day^-1] Denominator for ingestion function of copepods during day
    NP0 = P.INP + P.ENPp.*pref('predcop','phyto').*repmat(P.R,P.n,1)+sum(P.ENPd.*repmat(pref('predcop','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n),[2,1,3]).*RESC,3); % [gC day^-1] Denominator for the ingestion function                      
    NM1 = P.IDM + sum(P.EDMd.*repmat(pref('meso','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n).*RESC,3)+P.EDMC.*pref('meso','copepod').*repmat(MCday,1,P.n)+P.EDMP.*pref('meso','predcop').*repmat(MPday,1,P.n); % [gC day^-1] Denominator for ingestion function of mesopelagic during day
    NM0 = P.INM + sum(P.ENMd.*repmat(pref('meso','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n).*RESC,[2,1,3]),3)+P.ENMC.*pref('meso','copepod').*repmat(MCnight',P.n,1)+P.ENMP.*pref('meso','predcop').*repmat(MPnight',P.n,1); % [gC day^-1] Denominator for the ingestion function                   


    IFC1 = P.IDF.*P.EDFC*pref('forage','copepod').* repmat(MCday ,1,P.n)./NF1; % [gC day^-1] Ingestion rate of copepods during daytime by forage fish
    IFC0 = P.INF.*P.ENFC*pref('forage','copepod').* repmat(MCnight',P.n,1)./NF0;
    IFP1 = P.IDP.*P.EDFP*pref('forage','predcop').* repmat(MPday ,1,P.n)./NF1; % [gC day^-1] Ingestion rate of copepods during daytime by forage fish
    IFP0 = P.INP.*P.ENFP*pref('forage','predcop').* repmat(MPnight',P.n,1)./NF0;
    IFM1 = P.IDF.*P.EDFM*pref('forage','meso')    .*repmat(MMday ,1,P.n)./NF1;
    IFM0 = P.INF.*P.ENFM*pref('forage','meso')    .*repmat(MMnight',P.n,1)./NF0;

    IAF1 = P.IDA.*P.EDAF*pref('top','forage').* repmat(MFday ,1,P.n)./NA1; % [gC day^-1] Ingestion rate of forage fish during daytime by top predators
    IAF0 = P.INA.*P.ENAF*pref('top','forage').* repmat(MFnight',P.n,1)./NA0;
    IAJ1 = P.IDA.*P.EDAJ*pref('top','tactile').*repmat(MJday,1,P.n)  ./NA1;
    IAJ0 = P.INA.*P.ENAJ*pref('top','tactile').*repmat(MJnight' ,P.n,1)  ./NA0;
    IAM1 = P.IDA.*P.EDAM*pref('top','meso') .*repmat(MMday,1,P.n)./NA1; % [gC day^-1]
    IAM0 = P.INA.*P.ENAM*pref('top','meso') .*repmat(MMnight', P.n,1)./NA0;
    
    ICR1 = P.IDC.*P.EDCp*pref('copepod','phyto') .*repmat(P.R',1,P.n)./NC1; % [gC day^-1] Ingestion rate of phytoplankton during daytime by copepods
    ICR0 = P.INC.*P.ENCp*pref('copepod','phyto') .*repmat(P.R, P.n,1)./NC0;
        
    IPR1 = P.IDP.*P.EDPp*pref('predcop','phyto') .*repmat(P.R',1,P.n)./NP1; % [gC day^-1] Ingestion rate of phytoplankton during daytime by copepods
    IPR0 = P.INP.*P.ENPp*pref('predcop','phyto') .*repmat(P.R, P.n,1)./NP0;
   
    IJC1 = P.EDJC*pref('tactile','copepod').*repmat(MCday,1,P.n); % [gC day^-1] Ingestion rate of copepods during daytime by tactile predators - note the Type I functional response
    IJC0 = P.ENJC*pref('tactile','copepod').*repmat(MCnight',P.n,1);
    IJP1 = P.EDJP*pref('tactile','predcop').*repmat(MPday,1,P.n); % [gC day^-1] Ingestion rate of copepods during daytime by tactile predators - note the Type I functional response
    IJP0 = P.ENJP*pref('tactile','predcop').*repmat(MPnight',P.n,1);
    IJM1 = P.EDJM*pref('tactile','meso').*repmat(MMday,1,P.n);
    IJM0 = P.ENJM*pref('tactile','meso').*repmat(MMnight',P.n,1);

    IMC1 = P.IDM.*P.EDMC*pref('meso','copepod') .*repmat(MCday,1,P.n)./NM1; % [gC day^-1] Ingestion rate of copepods during daytime by mesopelagic fish
    IMC0 = P.INM.*P.ENMC*pref('meso','copepod') .*repmat(MCnight', P.n,1)./NM0;
    IMP1 = P.IDM.*P.EDMP*pref('meso','predcop') .*repmat(MPday,1,P.n)./NM1; % [gC day^-1] Ingestion rate of copepods during daytime by mesopelagic fish
    IMP0 = P.INM.*P.ENMP*pref('meso','predcop') .*repmat(MPnight', P.n,1)./NM0;

        
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
    
    

%Faecal pellet production rates
    FecC = (P.sigma*((1-P.fCR)*ICR1+(1-P.fCd)*sum(ICD1,3))+(1-P.sigma)*((1-P.fCR)*ICR0+(1-P.fCd)*sum(ICD0,3)))/P.wC; % [day^-1]
    FecP = (P.sigma*((1-P.fPR)*IPR1+(1-P.fPd)*sum(IPD1,3))+(1-P.sigma)*((1-P.fPR)*IPR0+(1-P.fPd)*sum(IPD0,3)))/P.wP; % [day^-1]
    FecF = (1-P.fF)*(P.sigma*(IFP1+IFC1+sum(IFD1,3)+IFM1)+(1-P.sigma)*(sum(IFD0,3)+IFC0+IFP0+IFM0))/P.wF; % [day^-1]
    FecA = (1-P.fA)*(P.sigma*(IAF1+IAJ1+IAM1)+(1-P.sigma)*(IAF0+IAJ0+IAM0))/P.wA; % [day^-1]
    FecJ = (1-P.fJ)*(P.sigma*(IJC1+IJP1+IJM1)+(1-P.sigma)*(IJC0+IJP0+IJM0))/P.wJ; % [day^-1]
    FecM = (P.sigma*((1-P.fMd)*sum(IMD1,3)+(1-P.fMC)*IMC1+(1-P.fMC)*IMP1)+(1-P.sigma)*((1-P.fMd)*sum(IMD0,3)+(1-P.fMC)*IMC0+(1-P.fMC)*IMP0))/P.wM; % [day^-1]

    
    sourceC = P.C*P.n*(sum(FecC.*C,2)*P.sigma+(1-P.sigma)*sum(FecC.*C,1)'); % [gC / m^3 / day]
    sourceP = P.P*P.n*(sum(FecP.*PC,2)*P.sigma+(1-P.sigma)*sum(FecC.*PC,1)');
    sourceM = P.M*P.n*(sum(FecM.*M,2)*P.sigma+(1-P.sigma)*sum(FecC.*M,1)');
    sourceF = P.F*P.n*(sum(FecF.*F,2)*P.sigma+(1-P.sigma)*sum(FecC.*F,1)');
    sourceA = P.A*P.n*(sum(FecA.*A,2)*P.sigma+(1-P.sigma)*sum(FecC.*A,1)');
    sourceJ = P.J*P.n*(sum(FecJ.*J,2)*P.sigma+(1-P.sigma)*sum(FecC.*J,1)');
    
 
 Fec_C = sum(sourceC*P.dZ); % [gC / m^2 / day]
 Fec_P = sum(sourceP*P.dZ);
 Fec_M = sum(sourceM*P.dZ);
 Fec_F = sum(sourceF*P.dZ);
 Fec_A = sum(sourceA*P.dZ);
 Fec_J = sum(sourceJ*P.dZ);
 
 Fec_C_Net = sum((sourceC-ConsD(:,2))*P.dZ); % [gC / m^2 / day] fecal pellets that will be degraded and not re eaten
 Fec_P_Net = sum((sourceP-ConsD(:,3))*P.dZ);
 Fec_M_Net = sum((sourceM-ConsD(:,4))*P.dZ);
 Fec_F_Net = sum((sourceF-ConsD(:,5))*P.dZ);
 Fec_A_Net = sum((sourceA-ConsD(:,6))*P.dZ);
 Fec_J_Net = sum((sourceJ-ConsD(:,7))*P.dZ);
  
 disp(['Fecal: ', num2str([Fec_C Fec_P Fec_M Fec_F Fec_A Fec_J])])
 
 disp(['Resp: ', num2str([Resp_C Resp_P Resp_M Resp_F Resp_A Resp_J])])
 
 disp(['Ratio: ', num2str([Fec_C_Net Fec_P_Net Fec_M_Net Fec_F_Net Fec_A_Net Fec_J_Net]./[Resp_C Resp_P Resp_M Resp_F Resp_A Resp_J])])
 
 