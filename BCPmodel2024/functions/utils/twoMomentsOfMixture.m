function [expectation, variance] = twoMomentsOfMixture(expectations, variances, weights, dim)
% Compute the expectation and variance of 1D mixture distributions over the 
% specified dimension 'dim'. By default dim=1, i.e. each column of the input
% matrices represents a mixture distribution. 
%
% https://en.wikipedia.org/wiki/Mixture_distribution#Moments
%
% David Meijer
% 31-08-2023

if nargin < 4
    dim = 1;
end

% Check that all input variables have the same size
if ~isempty(variances)
    assert(isequal(size(expectations), size(variances)), 'All input variables must have equal size');
    assert(all(variances >= 0,'all'), 'All variances must be >= 0');
end
assert(isequal(size(expectations), size(weights)), 'All input variables must have equal size');
assert(all(weights >= 0,'all'), 'All weights must be >= 0');

% Catch annoying special case
if isempty(expectations)
    expectation = [];
    variance = [];
    return;
end

% Normalize the weights
sum_of_weights = sum(weights,dim);
assert(all(sum_of_weights >= 0,'all'), 'At least one of the weights in the mixture must be > 0');
weights = weights ./ sum_of_weights;

% Avoid 0*inf = NaN issues etc
expectations(weights == 0) = 0;

% Compute the weighted expectation
expectation = sum(weights.*expectations, dim);

% Also compute the weighted variance if requested
if nargout > 1
    assert(~isempty(variances),'Cannot compute the variance without the variances of the components.')
    
    % Avoid 0*inf = NaN issues etc
    variances(weights == 0) = 0;
    
    % Compute the weighted variance
    variance = sum(weights.*(variances+expectations.^2), dim)-expectation.^2;
    
    % Correct for small numerical errors
    variance = max(variance,0);
end

end %[EoF]
