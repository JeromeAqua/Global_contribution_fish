load C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\ocean_data_jerome.mat
load C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\data_jerome2.mat
longitude = [0:2:178, -180:2:-2]; %What we will use for our runs
latitude = -90:2:90;

idxlattest = [68 46 47 49 40];
idxlongtest = [78 111 111 111 113];

add_on = P.ZMAX:200:5000; % [m]
Zdeep = [P.zi, add_on]; % [m] a deeper water column - to prevent extrapolations as NaN when the traps are deeper than P.ZMAX

POC_observed = [];
POC_computed = [];
Z_traps = [];

for ii=2%:5
    filename = strcat('TEST_2_lat_',num2str(latitude(idxlattest(ii))),'_long_',num2str(longitude(idxlongtest(ii))),'.mat');
    load(filename,'-regexp', '^(?!ii|POC_observed|POC_computed|Z_traps)\w')
    
    Dmean = [Dmean; repmat(Dmean(end,:),size(add_on,2),1).*exp(-repmat(P.alpha(end,:)./P.SR,size(add_on,2),1).*repmat(add_on'-P.zi(end),1,7))];
    
    [idxpoc] = find( ~isnan(squeeze(FPOC_obs(idxlattest(ii),idxlongtest(ii),1:15))));
    poc = squeeze(12.01*exp(FPOC_obs(idxlattest(ii),idxlongtest(ii),idxpoc)) * 10^-3 /365); % [gC / m^2 / day]
    
    POC_observed = [POC_observed; poc];
    Z_traps = [Z_traps; zt(idxpoc)];
    
    flux = Dmean.*P.SR;
    
    for k=1:size(idxpoc,1)
        f = 0;
        for j=1:7
            ftemp = interp1(Zdeep, flux(:,j)', zt(idxpoc(k))); % [gC /m^2 /day]
            f = f + ftemp;
        end
        POC_computed = [POC_computed; f];
    end
end