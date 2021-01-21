%% First find the limit of the euphotic zone
zeupho = -log(0.1) ./ P.klight; % [m] Depth at which we receive 1% of the surface light. Solve 0.01Is = Is exp(-l*z)

%% Then compute the fluxes
P = Pbasic;
P.scarc = [300 500 2000 2000 3000 800]; % [m/day] sinking rate of carcasses for C P M F A J

FluxEZ = zeros(1,nrun);
Respibelow = zeros(1,nrun);
FPbelow = zeros(1,nrun);
Carcbelow = zeros(1,nrun);
OLbelow = zeros(1,nrun);

OL = zeros(nrun,P.n,6);

for RUN = 1: nrun
    
    P.R = factorbio(1,RUN)*Pbasic.R;
    P.A = factorbio(6,RUN)*Pbasic.A;
    P.C = factorbio(2,RUN)*Pbasic.C;
    P.P = factorbio(3,RUN)*Pbasic.P;
    P.M = factorbio(4,RUN)*Pbasic.M;
    P.F = factorbio(5,RUN)*Pbasic.F;
    P.J = factorbio(7,RUN)*Pbasic.J;
    P.fCR = factor_assimCR(RUN)*Pbasic.fCR;
    P.fCd = factor_assimCD(RUN)*Pbasic.fCd;
    P.fPR = factor_assimPR(RUN)*Pbasic.fPR;
    P.fPd = factor_assimPD(RUN)*Pbasic.fPd;
    P.fMC = factor_assimMC(RUN)*Pbasic.fMC;
    P.fMd = factor_assimMD(RUN)*Pbasic.fMd;
    P.fF =  factor_assimF(RUN)*Pbasic.fF;
    P.fA =  factor_assimA(RUN)*Pbasic.fA;
    P.fJ =  factor_assimJ(RUN)*Pbasic.fJ;
    P.SR =  factorsink(:,RUN)'.*Pbasic.SR;
    P.alpha = factordegrad(RUN)*Pbasic.alpha;

            %%%%%%%%%%%%%%%%%%%%%%%%% CARCASSES CREATION RATE %%%%%%%%%%%%%
            
            Dead_C = reshape(P.n*P.C*(P.sigma*sum(P.minimortC.*squeeze(Glob_CT(RUN,:,:)),2)+(1-P.sigma)*sum(P.minimortC.*squeeze(Glob_CT(RUN,:,:)),1)'),1,P.n); % [gC / m3 / day] How much carcasses are created per day
                     
            Dead_P = reshape(P.n*P.P*(P.sigma*sum((P.minimortP+0.1*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_PT(RUN,:,:)),2)+(1-P.sigma)*sum((P.minimortP+0.1*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_PT(RUN,:,:)),1)'),1,P.n);
              
            Dead_M = reshape(P.n*P.M*(P.sigma*sum((P.minimortM+0.01*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_MT(RUN,:,:)),2)+(1-P.sigma)*sum((P.minimortM+0.01*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_MT(RUN,:,:)),1)'),1,P.n);
            
            Dead_F = reshape(P.n*P.F*(P.sigma*sum(P.minimortF.*squeeze(Glob_FT(RUN,:,:)),2)+(1-P.sigma)*sum(P.minimortF.*squeeze(Glob_FT(RUN,:,:)),1)'),1,P.n); 
                     
            Dead_A = reshape(P.n*P.A*(P.sigma*sum(P.minimortA.*squeeze(Glob_AT(RUN,:,:)),2)+(1-P.sigma)*sum(P.minimortA.*squeeze(Glob_AT(RUN,:,:)),1)'),1,P.n);
                      
            Dead_J = reshape(P.n*P.J*(P.sigma*sum(P.minimortJ.*squeeze(Glob_JT(RUN,:,:)),2)+(1-P.sigma)*sum(P.minimortJ.*squeeze(Glob_JT(RUN,:,:)),1)'),1,P.n); 
            
            
            % Now using it to compute the creation rate of DIC thanks to
            % carcasse degradation
            Dead_z = zeros(P.n,6);
            Dead_z(1,1) = (Dead_C(1)) /(P.scarc(1)/P.dZ+P.alpha(1, 2)); %For now unit is [gC / m3]
            Dead_z(1,2) = (Dead_P(1)) /(P.scarc(2)/P.dZ+P.alpha(1, 3));
            Dead_z(1,3) = (Dead_M(1)) /(P.scarc(3)/P.dZ+P.alpha(1, 4));
            Dead_z(1,4) = (Dead_F(1)) /(P.scarc(4)/P.dZ+P.alpha(1, 5));
            Dead_z(1,5) = (Dead_A(1)) /(P.scarc(5)/P.dZ+P.alpha(1, 6));
            Dead_z(1,6) = (Dead_J(1)) /(P.scarc(6)/P.dZ+P.alpha(1, 7));
            
                for depthindex = 2:P.n
                    Dead_z(depthindex,1) = (Dead_C(depthindex) + P.scarc(1)*Dead_z(depthindex-1,1)/P.dZ )/(P.scarc(1)/P.dZ+P.alpha(depthindex, 2));
                    Dead_z(depthindex,2) = (Dead_P(depthindex) + P.scarc(2)*Dead_z(depthindex-1,2)/P.dZ )/(P.scarc(2)/P.dZ+P.alpha(depthindex, 3));
                    Dead_z(depthindex,3) = (Dead_M(depthindex) + P.scarc(3)*Dead_z(depthindex-1,3)/P.dZ )/(P.scarc(3)/P.dZ+P.alpha(depthindex, 4));
                    Dead_z(depthindex,4) = (Dead_F(depthindex) + P.scarc(4)*Dead_z(depthindex-1,4)/P.dZ )/(P.scarc(4)/P.dZ+P.alpha(depthindex, 5));
                    Dead_z(depthindex,5) = (Dead_A(depthindex) + P.scarc(5)*Dead_z(depthindex-1,5)/P.dZ )/(P.scarc(5)/P.dZ+P.alpha(depthindex, 6));
                    Dead_z(depthindex,6) = (Dead_J(depthindex) + P.scarc(6)*Dead_z(depthindex-1,6)/P.dZ )/(P.scarc(6)/P.dZ+P.alpha(depthindex, 7));                
                end 
            
            Deg_carcasse = Dead_z.*reshape(P.alpha(:,2:end),1,1,P.n,6); % [gC / m3 / day] Now is the rate of degradation of carcasses at each depth
            
            Source_carc = cat(3, Dead_C, Dead_P, Dead_M, Dead_F, Dead_A, Dead_J); % [gC / m3 / day] How much carcasse is created per day   
    
            %%%%%%%%%%%%%%%%%%%%%%%OTHER LOSSES CALCULATION %%%%%%%%%%%%%
            
            %%% For C %%%
            sdac = McdT(RUN)*P.fCd + McrT(RUN)*P.fCR - sum(squeeze(DIC_globT(RUN,:,1))) -...
                   MpcT(RUN) - MfcT(RUN) - MmcT(RUN) - MjcT(RUN) -...
                   sum(P.dZ*P.n*P.C*(P.sigma*sum(P.minimortC.*squeeze(Glob_CT(RUN,:,:)),2)+(1-P.sigma)*sum(P.minimortC.*squeeze(Glob_CT(RUN,:,:)),1)')); % [gC / m2 / day]
                     
            CCday = squeeze(Glob_CdayT(RUN,:))/max(max(squeeze(Glob_CnightT(RUN,:))),max(squeeze(Glob_CdayT(RUN,:)))); % [-] Fraction of C at z during day
            CCnight = squeeze(Glob_CnightT(RUN,:))/max(max(squeeze(Glob_CnightT(RUN,:))),max(squeeze(Glob_CdayT(RUN,:)))); % [-]
            
            CCmean = CCday*P.sigma + (1-P.sigma)*CCnight; %mean distribution of C in the WC
            
            Q10C  = P.QC.^((P.T-P.T0C)/10);
            
            SDA_Cz = sdac .*(reshape(Q10C.*CCmean,1,P.n)) / sum(Q10C'.*CCmean.*P.dZ); % other losses at each depth [gC / m3 /day]
                     
                     
                     
            %%% For P %%%        
            sdap = MpdT(RUN)*P.fPd + MprT(RUN)*P.fPR + MpcT(RUN)*P.fPR - sum(squeeze(DIC_globT(RUN,:,2))) -...
                   MmpT(RUN) - MfpT(RUN) - MjpT(RUN) -...
                   sum(P.dZ*P.n*P.P*(P.sigma*sum((P.minimortP+0.1*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_PT(RUN,:,:)),2)+(1-P.sigma)*sum((P.minimortP+0.1*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_PT(RUN,:,:)),1)'));
              
            PPday = squeeze(Glob_PdayT(RUN,:))/max(max(squeeze(Glob_PnightT(RUN,:))),max(squeeze(Glob_PdayT(RUN,:)))); % [-] Fraction of P at z during day
            PPnight = squeeze(Glob_PnightT(RUN,:))/max(max(squeeze(Glob_PnightT(RUN,:))),max(squeeze(Glob_PdayT(RUN,:)))); % [-]
            
            PPmean = PPday*P.sigma + (1-P.sigma)*PPnight; %mean distribution of P in the WC
            
            Q10P  = P.QP.^((P.T-P.T0P)/10);
            
            SDA_Pz = sdap .*(reshape(Q10P.*PPmean,1,P.n)) / sum(Q10P'.*PPmean.*P.dZ); % other losses at each depth [gC / m3 /day]
            
            
            %%% For M %%%      
            sdam = MmcT(RUN)*P.fMC + MmpT(RUN)*P.fMC - sum(squeeze(DIC_globT(RUN,:,3))) -...
                   MfmT(RUN) - MamT(RUN) -...
                   sum(P.dZ*P.n*P.M*(P.sigma*sum((P.minimortM+0.5*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_MT(RUN,:,:)),2)+(1-P.sigma)*sum((P.minimortM+0.5*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_MT(RUN,:,:)),1)'));

            MMday = squeeze(Glob_MdayT(RUN,:))/max(max(squeeze(Glob_MnightT(RUN,:))),max(squeeze(Glob_MdayT(RUN,:)))); % [-] Fraction of M at z during day
            MMnight = squeeze(Glob_MnightT(RUN,:))/max(max(squeeze(Glob_MnightT(RUN,:))),max(squeeze(Glob_MdayT(RUN,:)))); % [-]
            
            MMmean = MMday*P.sigma + (1-P.sigma)*MMnight; %mean distribution of M in the WC
            
            Q10M  = P.QM.^((P.T-P.T0M)/10);
            
            SDA_Mz = sdam .*(reshape(Q10M.*MMmean,1,P.n)) / sum(Q10M'.*MMmean.*P.dZ); % other losses at each depth [gC / m3 /day]
                     
                     
            %%% For F %%%         
            sdaf = MfcT(RUN)*P.fF + MfpT(RUN)*P.fF + MfmT(RUN)*P.fF - sum(squeeze(DIC_globT(RUN,:,4))) -...
                   MafT(RUN) -...
                   sum(P.dZ*P.n*P.F*(P.sigma*sum(P.minimortF.*squeeze(Glob_FT(RUN,:,:)),2)+(1-P.sigma)*sum(P.minimortF.*squeeze(Glob_FT(RUN,:,:)),1)')); 
                     
            FFday = squeeze(Glob_FdayT(RUN,:))/max(max(squeeze(Glob_FnightT(RUN,:))),max(squeeze(Glob_FdayT(RUN,:)))); % [-] Fraction of F at z during day
            FFnight = squeeze(Glob_FnightT(RUN,:))/max(max(squeeze(Glob_FnightT(RUN,:))),max(squeeze(Glob_FdayT(RUN,:)))); % [-]
            
            FFmean = FFday*P.sigma + (1-P.sigma)*FFnight; %mean distribution of F in the WC
            
            Q10F  = P.QF.^((P.T-P.T0F)/10);
            
            SDA_Fz = sdaf .*(reshape(Q10F.*FFmean,1,P.n)) / sum(Q10F'.*FFmean.*P.dZ); % other losses at each depth [gC / m3 /day]
             
            
            %%% For A %%%
            sdaa = MafT(RUN)*P.fA + MamT(RUN)*P.fA + MajT(RUN)*P.fA - sum(squeeze(DIC_globT(RUN,:,5))) -...
                   sum(P.dZ*P.n*P.A*(P.sigma*sum(P.minimortA.*squeeze(Glob_AT(RUN,:,:)),2)+(1-P.sigma)*sum(P.minimortA.*squeeze(Glob_AT(RUN,:,:)),1)'));
                     
                     
            AAday = squeeze(Glob_AdayT(RUN,:))/max(max(squeeze(Glob_AnightT(RUN,:))),max(squeeze(Glob_AdayT(RUN,:)))); % [-] Fraction of A at z during day
            AAnight = squeeze(Glob_AnightT(RUN,:))/max(max(squeeze(Glob_AnightT(RUN,:))),max(squeeze(Glob_AdayT(RUN,:)))); % [-]
            
            AAmean = AAday*P.sigma + (1-P.sigma)*AAnight; %mean distribution of A in the WC
            
            Q10A  = P.QA.^((P.T-P.T0A)/10);
            
            SDA_Az = sdaa .*(reshape(Q10A.*AAmean,1,P.n)) / sum(Q10A'.*AAmean.*P.dZ); % other losses at each depth [gC / m3 /day]
             
            
            %%% For J %%%
            sdaj = MjcT(RUN)*P.fJ + MjpT(RUN)*P.fJ - sum(squeeze(DIC_globT(RUN,:,6))) -...
                   MajT(RUN) -...
                   sum(P.dZ*P.n*P.J*(P.sigma*sum(P.minimortJ.*squeeze(Glob_JT(RUN,:,:)),2)+(1-P.sigma)*sum(P.minimortJ.*squeeze(Glob_JT(RUN,:,:)),1)')); 
            
            JJday = squeeze(Glob_JdayT(RUN,:))/max(max(squeeze(Glob_JnightT(RUN,:))),max(squeeze(Glob_JdayT(RUN,:)))); % [-] Fraction of J at z during day
            JJnight = squeeze(Glob_JnightT(RUN,:))/max(max(squeeze(Glob_JnightT(RUN,:))),max(squeeze(Glob_JdayT(RUN,:)))); % [-]
            
            JJmean = JJday*P.sigma + (1-P.sigma)*JJnight; %mean distribution of J in the WC
            
            Q10J  = P.QJ.^((P.T-P.T0J)/10);
            
            SDA_Jz = sdaj .*(reshape(Q10J.*JJmean,1,P.n)) / sum(Q10J'.*JJmean.*P.dZ); % other losses at each depth [gC / m3 /day]
            
            
            OL(RUN,:,1) = SDA_Cz;
            OL(RUN,:,2) = SDA_Pz;
            OL(RUN,:,3) = SDA_Mz;
            OL(RUN,:,4) = SDA_Fz;
            OL(RUN,:,5) = SDA_Az;
            OL(RUN,:,6) = SDA_Jz;
            
             %%%%%%%%%%%%%%%%%%%%%%%%% SINKING FLUX %%%%%%%%%%%%%%%%%%%%%%%%%%%
             
             s = P.SR.*squeeze(D_globT(RUN,:,:)); % [gC / m^2 / day] sinking flux of fecal pellets
             s2 = P.scarc.*Dead_z; % [gC / m2 / day] % sinking flux of carcasses P.dZ*squeeze(SOURCE(j,i,:,carc_considered)-Deg_carcasse(j,i,:,carc_considered)); %
               
             s = sum(s(:,2:end),2);%2:end),2); % these two lines are just to look at sinking rates by functional groups
             s2 = sum(s2(:,:),2);%1:end),2);
             
             sinking_flux = interp1(P.zi, s, zeupho);%interp1(P.zi, sum(s(:,2:end),2), zeupho);
             sinking_flux2 = interp1(P.zi, s2, zeupho);%interp1(P.zi, sum(s2,2), zeupho);
             FluxEZ(RUN) = sinking_flux2 +sinking_flux; % [gC / m2 / day]


             %%%%%%%%%%%%%%%%%%%%%%%%% RESPIRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

             R = squeeze(DIC_globT(RUN,:,:)); %Here last number to choose if we want to compute how much was created by each group
             
             Respibelow(RUN) = sum(sum(R(P.zi'>zeupho,:))); % [gC m^-2 day^-1]
             
             
             %%%%%%%%%%%%%%%%%%%%%%%% FECAL PELLETS %%%%%%%%%%%%%%%%%%%%%%%%%%%
           
             S = squeeze(Glob_sourceT(RUN,:,:)); %Here last number to choose if we want to compute how much was created by each group
             DDD = squeeze(Glob_DconsT(RUN,:,:));
             
             int = sum(S-DDD,2); % [gC / m^3 / day] total source term at each depth -DDD
             FPbelow(RUN) = sum(int(P.zi'>zeupho)*P.dZ); % [gC / m^2 / day]
             
             
             %%%%%%%%%%%%%%%%%%%%%%%  CARCASSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             S2 = squeeze(Source_carc); %Here last number to choose if we want to compute how much was created by each group
             
             int = sum(S2,2); % [gC / m^3 / day] total source term at each depth -DDD
             Carcbelow(RUN) = sum(int(P.zi'>zeupho)*P.dZ); % [gC / m^2 / day]
             
             
             %%%%%%%%%%%%%%%%%%%%%%%  OTHER LOSSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             R2 = squeeze(OL(RUN,:,:)); %Here last number to choose if we want to compute how much was created by each group
             
             OLbelow(RUN) = sum(sum(R2(P.zi'>zeupho,:))); % [gC m^-2 day^-1]

end

%% Now plot the fluxes

BP = 10^3*[FluxEZ' FPbelow' Respibelow' Carcbelow' OLbelow'];
fluxes = {'Sinking flux', 'Fecal pellets', 'Respiration', 'Carcasses', 'Other losses'};
boxplot(BP,fluxes)
hold on
plot(BP(1,:),'*r')
ylabel({'Flux or net production below the','euphotic zone [mgC / m^2 / day]'})