function fit = fit_expcurve_2(parms, data)

nSubs   = size(data,1); % rows = different subject data
nTrials = length(data); % colns = trial data for each subject

fit = 0;

for i = 1:nSubs
    
    y = expFun(parms,nTrials);
    fit = fit + rmsd(y,data);
    
end % for i ...


end % of function... 
