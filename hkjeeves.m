function [optParms, fit, nFuncEvals] = hkjeeves(objectiveFunction, parms, varargin)
% HKJEEVES Unconstrained optimization using Hooke & Jeeves algorithm
%   [OPTPARMS, FIT, NFUNCEVALS] = HKJEEVES(OBJECTIVEFUNCTION, PARMS) 
%   [OPTPARMS, FIT, NFUNCEVALS] = HKJEEVES(OBJECTIVEFUNCTION, PARMS,...
%                                          IP, LOWERBOUND, UPPERBOUND, PROBLEM, TOL, MAXITER,...
%                                          PARMSTEP, AMPLIFICATION, REDUCTION, OBJFUN_ARGS)
%
%   OBJECTIVEFUNCTION	: objective function
%   PARMS               : initial point
%
%   Optional Inputs
%   IP              : (0): no plot (default), (>0) plot figure ip with pause (<0) plot figure ip - only works if nparms = 2
%   LOWERBOUND, UPPERBOUND	    : lower and upper bounds on parameters (default = parms * (1 +/- 2))
%   PROBLEM         : (-1): minimum (default), (1): maximum
%   TOL             : tolerance (default = 1e-4)
%   MAXITER         : maximum number of stages (default = 50*(1+4*~(ip>0)))
%   PARMSTEP        : stepsize vector for the independent variables (default = max(0.01*abs(guess+~guess),0.1))
%   AMPLIFICATION	: stepsize enlargement factor (1,oo) (default = 1.5)
%   REDUCTION		: stepsize reduction factor (0,1) (default = 0.5)
%
%   OBJFUN_ARGS : parameters needed to function

%   Outputs
%   OPTPARMS		: optimal parameters
%   FIT         	: optimal fit of objective function
%   NFUNCEVALS		: number of objective function evaluations

%% Parse Optional Arguments
optargs = {0,...                                % IP
           -parms - ~parms,...                  % LOWERBOUND
           2*parms + ~parms,...                 % UPPERBOUND  
           -1,...                               % PROBLEM [MINIMUM]
           1e-4,...                             % TOL
           50 * (1 + 4 * ~(varargin{1} > 0)),...% MAXITER
           max(0.01*abs(parms+~parms),0.1),...  % PARMSTEP
           1.5,...                              % AMPLIFICATION
           0.5,...                              % REDUCTION
           varargin(10:end)};                   % OBJFUN_ARGS
       
newVals = cellfun(@(x) ~isempty(x), varargin); % skip any new inputs if they are empty
optargs(newVals) = varargin(newVals); % now put these defaults into the valuesToUse cell array, and overwrite the ones specified in varargin.

% Place optional args in memorable variable names
[ip, lowerParmBound, upperParmBound, problem, tol, maxIter, parmstep, parmStepAmplification, parmStepReduction] = optargs{1:9}; 
if size(optargs, 2) > 9
    funargs = optargs(10:end);
end

%% Set up variables
parms                 = parms(:);
parmstep              = abs(parmstep(:));
parmStepAmplification = parmStepAmplification(:);
parmStepReduction     = parmStepReduction(:);

%% Run Algorithm
initialFit = feval(objectiveFunction, parms, funargs{:}) * problem; % Get initial function estimate using guess
nParms     = size(parms,1);

% Initalise
startParms      = parms; % Starting to-be-optimized parameters
currentOptParms = parms;      % Current optimal parameters
fit             = initialFit;    % Starting function value
startingFit     = initialFit;

iteration     = 0;
nFuncEvals    = 1;
top    = 0; % end of idx list
bottom = 0; % begin of idx list
idx    = zeros(nParms+1,1); %  variables index vector with enlarged stepsize
idx(bottom+1) = top;

while iteration < maxIter
    next = bottom; % index of enlarged variable
    norma = 0;

    % exploration
    for i=1:nParms
        parmstep_i = parmstep(i);
        
        for direction=[1 -1] % On initial loop take a step up, then take a step down
            startParms(i) = currentOptParms(i) + parmstep_i * direction;
            
            % Reset any out of bounds parms
            startParms(startParms < lowerParmBound') = lowerParmBound(startParms < lowerParmBound');
            startParms(startParms > upperParmBound') = upperParmBound(startParms > upperParmBound');
            
            currentFit = feval(objectiveFunction,startParms, funargs{:}) * problem;
            nFuncEvals = nFuncEvals+1;
            
            if currentFit > initialFit     % success
                initialFit = currentFit;
                currentOptParms(i) = startParms(i);
                idx(next+1) = i;
                next = i;
                break;
            end
        end
        startParms(i) = currentOptParms(i);
        
        if parmStepReduction(i) < 1
            norma = norma + parmstep_i * parmstep_i/(startParms(i) * startParms(i) + (abs(startParms(i))<tol));
        end
    end
    
    iteration = iteration + 1;
    if ip % If ip == 1 then show iteration and fit
        disp([iteration fit])
    end
    
    if ~isreal(fit)
        warning('Fit is imaginary. Check parameters')
        return
    elseif isnan(fit)
        warning('Fit is NaN. Check parameters')
        return
    end
    

%     Why is this break criterion so complex? I think it's iterating based
%     on parameter reduction as well as fit, but if certain parameters
%     don't change their step size then the algorithm can't break based on
%     tolerance
    if sqrt(norma) < tol && abs(initialFit-fit) < tol*(0.1+abs(fit))
        break;
    end
    % progression
    if next == bottom
        parmstep = parmstep .* parmStepReduction;
    else
        good = 1;
        idx(next+1) = top;
        
        while good,
            fit = initialFit;
            
            next = idx(bottom+1);
            while next ~= top,
                if numel(parmStepAmplification) > 1
                    startParms(next) = startParms(next) + parmStepAmplification(next) .* (startParms(next) - parms(next));
                else
                    startParms(next) = startParms(next) + parmStepAmplification * (startParms(next) - parms(next));
                end
                parms(next) = currentOptParms(next);
                next=idx(next+1);
            end
            
            if idx(bottom+1) ~= top,
                startParms(startParms < lowerParmBound') = lowerParmBound(startParms < lowerParmBound');
                startParms(startParms > upperParmBound') = lowerParmBound(startParms > upperParmBound');
                currentFit = feval(objectiveFunction, startParms, funargs{:} ) * problem;
                
                nFuncEvals=nFuncEvals+1;
                
                if currentFit > initialFit,
                    initialFit=currentFit;
                    good=1;
                else
                    good = 0;
                end
                
                next = idx(bottom+1);
                while next ~= top,
                    if good
                        currentOptParms(next)=startParms(next);
                    else
                        startParms(next)=currentOptParms(next);
                    end
                    next=idx(next+1);
                end
            end
        end
    end
end

fprintf('%5.3f\t%5.3f\n',  startingFit, initialFit)
fit=initialFit*problem;
optParms = currentOptParms;

if iteration == maxIter,
    disp('Warning Hkjeeves: reached maximum number of iterations!');
end