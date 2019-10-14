function out = pref(predator,prey)
% [-] Preference values for each prey and predator
    
if strcmp(predator,'copepod')
    if strcmp(prey,'phyto')
        out = 1;
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
    elseif strcmp(prey,'meso')
        out = 0.1;
    else
        out = NaN;
    end  
    
elseif strcmp(predator,'top')
    if strcmp(prey,'forage')
        out = 1;
    elseif strcmp(prey,'tactile')
        out = 0.2;
    elseif strcmp(prey,'meso')
        out = 1;
    elseif strcmp(prey,'bathy')
        out = 0.5;
    else
        out = NaN;
    end        
              
elseif strcmp(predator,'tactile')
    if strcmp(prey,'copepod')
        out = 1;
    elseif strcmp(prey,'meso')
        out = 0;
    else
        out = NaN;
    end                
             
elseif strcmp(predator,'meso')
    if strcmp(prey,'detritus')
        out = 0.5;
    elseif strcmp(prey,'copepod')
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
        out = 0.0;
    elseif strcmp(prey,'meso')
        out = 0.8;
    else
        out = NaN;
    end                 
                                        
end

end



