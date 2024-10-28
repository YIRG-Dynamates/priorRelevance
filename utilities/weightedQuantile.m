function q = weightedQuantile(x,w,qr)
% Compute requested quantiles (qr) for weighted (w) values in vector x
%
% David Meijer, 21-12-2023

assert(numel(x) == numel(w),'vectors x and w must have same length');
assert(all(qr(:) >= 0) && all(qr(:) <= 1),'requested quantiles (qr) must be between 0 and 1');

%For each value in x, get the index that refers to its location in a list of unique values in x (N.B.  x = unique_x(idx_unique_x))
[unique_x,~,idx_unique_x] = unique(x(:));                     

%Sum the weights for each unique x
w_unique_x = accumarray(idx_unique_x,w(:));

%Find the middle cumulative probability of each set of unique values
cdf_upper_step_weights = cumsum(w_unique_x);
cdf_lower_step_weights = [0; cdf_upper_step_weights(1:(end-1))];      
cdf_middle_probs = .5*(cdf_lower_step_weights+cdf_upper_step_weights)/sum(w_unique_x);

%Ensure that the requested quantiles do not exceed the cumulative probabilities   
%This is another way to set the extrapolation method to the most extreme values 
qr_cropped = min(max(qr(:),cdf_middle_probs(1)),cdf_middle_probs(end));

%Linearly interpolate at the requested quantiles (qr)   
[~,idx_last_lower] = min(diff([cdf_middle_probs', 2] <= qr_cropped,1,2),[],2);      %[1 1 1 0];        
[~,idx_first_higher] = max(diff([-1, cdf_middle_probs'] >= qr_cropped,1,2),[],2);   %[0 1 1 1];
normalized_distance = (qr_cropped - cdf_middle_probs(idx_last_lower)) ./ (cdf_middle_probs(idx_first_higher) - cdf_middle_probs(idx_last_lower));
q = (1-normalized_distance).*unique_x(idx_last_lower) + normalized_distance.*unique_x(idx_first_higher);
q(idx_last_lower == idx_first_higher) = unique_x(idx_last_lower(idx_last_lower == idx_first_higher));                   %Overwrite entries where normalized_distance == NaN

%Resize to match input
q = reshape(q,size(qr));

end %[EoF

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simple Tester Code %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% x = 1:10;
% w = 1:10;
% qr = [0.25 0.5 0.75];
% 
% q1 = weightedQuantile(x,w,qr)
% q2 = weightedQuantile(x,fliplr(w),qr)
% 
% q3 = weightedQuantile(x,ones(1,10),[0 qr 1])
% q4 = quantile(x,[0 qr 1])
