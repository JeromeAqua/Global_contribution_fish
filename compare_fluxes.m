load C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\ocean_data_jerome.mat
load C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\data_jerome2.mat
longitude = [0:2:178, -180:2:-2]; %What we will use for our runs
latitude = -90:2:90;

% idxlattest = [68 46 47 49 40];
% idxlongtest = [78 111 111 111 113];

idxlattest = [36    53    48    55    55    48    64    65    68    47    66    65    69    63    46    47    49    52    40,...
            49    57    42    44    46    52    55    57    44    46    47    43    45];
idxlongtest = [ 5    43    44    45    59    68    74    78    78    81    83    88    88    89   111   111   111   111   113,...
            138   165   167   168   169   170   170   170   175   175   175   176   176];

add_on = P.ZMAX:200:5000; % [m]
Zdeep = [P.zi, add_on]; % [m] a deeper water column - to prevent extrapolations as NaN when the traps are deeper than P.ZMAX

POC_observed = [];
POC_computed = [];
Z_traps = [];

for ii=1:size(idxlattest,2)
    try
    filename = strcat('TEST_RUNS_lat_',num2str(latitude(idxlattest(ii))),'_long_',num2str(longitude(idxlongtest(ii))),'.mat');
    load(filename,'-regexp', '^(?!ii|POC_observed|POC_computed|Z_traps)\w')
    
    Niter= a; Iavg= b;  P= c; MAday= d; MAnight= e; MCday=f; MCnight=g; MPday=h; MPnight=i; MFday=j; MFnight=k; MJday=l; MJnight=m;
    MMday=n; MMnight=o; DegPOC_depth=p; DIC_dept=q; Dmean=r; MA=s; MC=t; MF=u; MJ=v; MM=w; MP=x; MD=y; FitA=z; FitC=aa; FitF=bb; FitJ=cc; FitM=dd; FitP=ee;
    
    add_on_D =  repmat(Dmean(end,:),size(add_on,2),1).*exp(-repmat(P.alpha(end,:)./P.SR,size(add_on,2),1).*repmat(add_on'-P.zi(end),1,7));
    
    [idxpoc] = find( ~isnan(squeeze(FPOC_obs(idxlattest(ii),idxlongtest(ii),1:15))));
    poc = squeeze(12.01*exp(FPOC_obs(idxlattest(ii),idxlongtest(ii),idxpoc)) * 10^-3 /365.25); % [gC / m^2 / day]
    
    POC_observed = [POC_observed; poc];
    Z_traps = [Z_traps, zt(idxpoc)];
    
    flux = [DegPOC_depth;add_on_D.*P.SR];
    
    for k=1:size(idxpoc,1)
        f = 0;
        for j=1:7
            ftemp = interp1(Zdeep, flux(:,j)', zt(idxpoc(k))); % [gC /m^2 /day]
            f = f + ftemp;
        end
        POC_computed = [POC_computed; f];
    end
    catch
        ii
    end
end

POC_observed = 10^3*POC_observed; % [mgC /m^2 /day]
POC_computed = 10^3*POC_computed; % [mgC /m^2 /day]

%%
a = colormap(flipud(jet));
figure
for kkk=1:length(POC_observed)
    plot(POC_observed(kkk), POC_computed(kkk),'o','MarkerFaceColor', a(floor(64*Z_traps(kkk)/2000),:)) %k');%
    hold on
end
mm = min([POC_observed; POC_computed]);
MM = max([POC_observed; POC_computed]);
plot([mm MM], [mm MM], 'k')
xlabel('Observed POC flux [gC/m^2/day]')
ylabel('Modeled POC flux [gC / m^2 /day]')

colorbar
caxis([50 2000])
error = 1/size(POC_observed,1)*sum((POC_observed-POC_computed).^2)

% hh = [];% 5to know the respective contributions of the different
% %populations to POC
%  for j=1:7
%             hh =[hh, interp1(Zdeep, flux(:,j)', zt(idxpoc(k)))]; % [gC /m^2 /day]
%       
%         end
% hh