function lpp = compLogTrapPrior(params,LB,PLB,PUB,UB)
%Compute the log prior probability for the combination of parameter values
%of a trapezoid-shaped prior. 
%We construct the prior PDFs as trapezoids (piecewise continuous): They are
%maximal between the plausible bounds PLB and PUB, and linearly going to 
%zero on either side between plausible bounds and hard bounds: from PLB to 
%LB, and from PUB to UB. 
%The combined log prior probability is the sum of the log probabilities for
%each parameter.   

%This function is vectorized and works with implicit expansion, so long as 
%the dimensions of params and bounds match. E.g. params may be NxP (PxN), 
%where all bounds must be 1xP (Px1). The vector output size is Nx1 (1xN).

%% Compute the height such that the area under the prior pdf integrates to 1.
maxHeightPriors = 1./(0.5*(PLB-LB)+(PUB-PLB)+0.5*(UB-PUB));

%% The method uses 'min' and 'max' rather than indexing, because it seems faster.   

%Compute all probabilities as if they were on the upwards slope ((params > LB) & (params < PLB))
pp_upward_slope = (params-LB)./(PLB-LB).*maxHeightPriors;

%Compute all probabilities as if they were on the downward slope ((params < UB) & (params > PUB))  
pp_downward_slope = (1-((params-PUB)./(UB-PUB))).*maxHeightPriors;

%Combine both and ... 
% - set lower bound to zero (for values outside/on hard bounds ((params <= LB) | (params >= UB)))
% - set upper bound to maxHeightPriors (for values between the plausible bounds ((params >= PLB) & (params <= PUB)))
pp = min(maxHeightPriors,max(0,min(pp_upward_slope,pp_downward_slope)));

%% Compute the log prior probabilities

%Find the parameters dimension (first non-sigleton dim of bounds)
sum_dim = find(size(LB)~=1,1,'first');

%Compute the overall log prior probability by summing the log probabilities of the parameters   
lpp = sum(log(pp),sum_dim);

end %[EOF]
