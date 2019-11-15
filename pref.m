function out = pref(predator,prey)
% [-] Preference values for each prey and predator
    
if strcmp(predator,'copepod')
    if strcmp(prey,'phyto')
        out = 1;
    elseif strcmp(prey,'detritus')
        out = 0.1;
    else
        out = NaN;
    end   
    
elseif strcmp(predator,'predcop')
    if strcmp(prey,'phyto')
        out = 0.0;
    elseif strcmp(prey,'detritus')
        out = 1;
    else
        out = NaN;
    end 
    
elseif strcmp(predator,'forage')
    if strcmp(prey,'detritus')
        out = 0;
    elseif strcmp(prey,'benthos')
        out = 0.0;
    elseif strcmp(prey,'copepod')
        out = 1;
    elseif strcmp(prey,'predcop')
        out = 1;
    elseif strcmp(prey,'meso')
        out = 1;
    else
        out = NaN;
    end  
    
elseif strcmp(predator,'top')
    if strcmp(prey,'forage')
        out = 1;
    elseif strcmp(prey,'tactile')
        out = 0.1;
    elseif strcmp(prey,'meso')
        out = 0.2;
    elseif strcmp(prey,'bathy')
        out = 0.5;
    else
        out = NaN;
    end        
              
elseif strcmp(predator,'tactile')
    if strcmp(prey,'copepod')
        out = 1;
    elseif strcmp(prey,'predcop')
        out = 1;
    elseif strcmp(prey,'meso')
        out = 0;
    else
        out = NaN;
    end                
             
elseif strcmp(predator,'meso')
    if strcmp(prey,'detritus')
        out = 0;
    elseif strcmp(prey,'copepod')
        out = 0.1;
    elseif strcmp(prey,'predcop')
        out = 1;
    else
        out = NaN;
    end                             
                
elseif strcmp(predator,'bathy')
    if strcmp(prey,'detritus')
        out = 0;
    elseif strcmp(prey,'benthos')
        out = 1;
    elseif strcmp(prey,'copepod')
        out = 0;
    elseif strcmp(prey,'predcop')
        out = 0;
    elseif strcmp(prey,'meso')
        out = 0.1;
    else
        out = NaN;
    end                 
                                        
end

if strcmp(prey,'detritus')
    n = 50;
    ZMAX = 1500;
    zext = linspace(0,ZMAX,n+1);
    zi = (zext(2:end)+zext(1:end-1))/2;
    fz = 20*ones(n,1);
    
    tresh = 400; % [m]
    fz(zi<tresh) = 1+19*zi(zi<tresh)/tresh;
    out = out*fz;
end

end



